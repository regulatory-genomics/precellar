mod segment;
use crate::region::Region;
use crate::Modality;
pub use segment::{Segment, SegmentInfo, SegmentInfoElem, SplitError};

use anyhow::{bail, Result};
use file_download::download::Downloader;
use indexmap::IndexMap;
use noodles_fastq as fastq;
use seqarc::{AlignedRecordGroup, DecompressionOptions, SequenceRecordGroup};
use serde::{Deserialize, Serialize, Serializer};
use serde_yaml::Value;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::ops::{Deref, DerefMut};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// A seekable byte source used by archive readers.
pub trait ReadSeek: std::io::Read + std::io::Seek {}
impl<T: std::io::Read + std::io::Seek> ReadSeek for T {}

/// Specification of a sequencing library.
#[derive(Debug, Clone, PartialEq)]
pub struct SeqSpec(IndexMap<String, Read>);

/// A grouped sequencing source represented by one multi-stream archive.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadGroup {
    /// Read definitions keyed by archive stream name (`R1`, `R2`, `I1`, `I2`).
    pub reads: IndexMap<String, Read>,
    /// The files backing the group.
    pub files: Vec<File>,
}

impl Deref for SeqSpec {
    type Target = IndexMap<String, Read>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SeqSpec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Serialize for SeqSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.values().collect::<Vec<_>>().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for SeqSpec {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let values = Vec::<Value>::deserialize(deserializer)?;
        let mut sequences = IndexMap::new();
        for value in values {
            let is_group = value.get("read_id").and_then(Value::as_mapping).is_some();
            if is_group {
                let group: ReadGroupValue =
                    serde_yaml::from_value(value).map_err(serde::de::Error::custom)?;
                for (stream, read_id) in group.read_id {
                    let read = Read {
                        read_id: read_id.clone(),
                        name: group.name.clone(),
                        modality: group.modality,
                        primer_id: group.primer_id.get(&stream).cloned().ok_or_else(|| {
                            serde::de::Error::custom(format!(
                                "missing primer_id for stream {stream}"
                            ))
                        })?,
                        min_len: *group.min_len.get(&stream).unwrap_or(&0),
                        max_len: *group.max_len.get(&stream).unwrap_or(&0),
                        strand: *group.strand.get(&stream).unwrap_or(&Strand::Pos),
                        files: Some(
                            group
                                .files
                                .iter()
                                .cloned()
                                .map(|mut file| {
                                    file.url.push('#');
                                    file.url.push_str(&stream);
                                    file
                                })
                                .collect(),
                        ),
                    };
                    if sequences.insert(read_id.clone(), read).is_some() {
                        return Err(serde::de::Error::custom(format!(
                            "Duplicate read id: {}",
                            &read_id
                        )));
                    }
                }
                continue;
            }
            let read: Read = serde_yaml::from_value(value).map_err(serde::de::Error::custom)?;
            let read_id = read.read_id.clone();
            if sequences.insert(read_id.clone(), read).is_some() {
                return Err(serde::de::Error::custom(format!(
                    "Duplicate read id: {}",
                    &read_id
                )));
            }
        }
        Ok(Self(sequences))
    }
}

#[derive(Deserialize)]
struct ReadGroupValue {
    read_id: HashMap<String, String>,
    name: Option<String>,
    modality: Modality,
    primer_id: HashMap<String, String>,
    #[serde(default)]
    min_len: HashMap<String, u32>,
    #[serde(default)]
    max_len: HashMap<String, u32>,
    #[serde(default)]
    strand: HashMap<String, Strand>,
    #[serde(default)]
    files: Vec<File>,
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Read {
    pub read_id: String,
    pub name: Option<String>,
    pub modality: Modality,
    pub primer_id: String,
    pub min_len: u32,
    pub max_len: u32,
    pub strand: Strand, // whether the read orientation is the same as the library
    pub files: Option<Vec<File>>,
}

impl Default for Read {
    fn default() -> Self {
        Self {
            read_id: "".to_string(),
            name: None,
            modality: Modality::DNA,
            primer_id: "".to_string(),
            min_len: 0,
            max_len: 0,
            strand: Strand::Pos,
            files: None,
        }
    }
}

pub struct FastqReader {
    reader: Box<dyn RecordSource + Send>,
    min_len: u32,
    max_len: u32,
}

trait RecordSource {
    fn read_record(&mut self, record: &mut fastq::Record) -> Result<usize>;
}

struct TextRecordSource {
    reader: fastq::io::Reader<Box<dyn BufRead + Send>>,
}

impl RecordSource for TextRecordSource {
    fn read_record(&mut self, record: &mut fastq::Record) -> Result<usize> {
        Ok(self.reader.read_record(record)?)
    }
}

struct SeqArcRecordSource {
    reader: Arc<Mutex<SeqArcSharedReader>>,
    stream: String,
}

impl RecordSource for SeqArcRecordSource {
    fn read_record(&mut self, record: &mut fastq::Record) -> Result<usize> {
        let mut reader = self
            .reader
            .lock()
            .map_err(|_| anyhow::anyhow!("SeqArc reader lock poisoned"))?;
        let Some(next) = reader.next_stream_record(&self.stream)? else {
            return Ok(0);
        };
        let converted = next;
        let length = converted.sequence().len();
        *record = converted;
        Ok(length)
    }
}

struct SeqArcSharedReader {
    readers: Vec<SeqArcReader>,
    current: usize,
    requested: HashSet<String>,
    served: HashSet<String>,
    group: Option<AlignedRecordGroup>,
}

impl SeqArcSharedReader {
    fn next_stream_record(&mut self, stream: &str) -> Result<Option<fastq::Record>> {
        if self.group.is_none() {
            loop {
                let Some(reader) = self.readers.get_mut(self.current) else {
                    return Ok(None);
                };
                match reader.next_group() {
                    Some(group) => {
                        self.group = Some(group?);
                        self.served.clear();
                        break;
                    }
                    None => self.current += 1,
                }
            }
        }
        let Some(group) = self.group.as_ref() else {
            return Ok(None);
        };
        let record = fastq_record_from_seqarc_group_stream(group, stream)?;
        self.served.insert(stream.to_string());
        if self.served.len() == self.requested.len() {
            self.group = None;
        }
        Ok(Some(record))
    }
}

/// Shared SeqArc archive readers for one assay execution.
#[derive(Clone, Default)]
pub struct SeqArcContext {
    archives: Arc<Mutex<HashMap<String, Arc<Mutex<SeqArcSharedReader>>>>>,
}

impl SeqArcContext {
    fn open_stream(&self, files: &[File], stream: String) -> Result<FastqReader> {
        let key = files
            .iter()
            .map(|file| file.url.split('#').next().unwrap_or(&file.url))
            .collect::<Vec<_>>()
            .join("\0");
        let shared = {
            let mut archives = self
                .archives
                .lock()
                .map_err(|_| anyhow::anyhow!("SeqArc context lock poisoned"))?;
            if let Some(shared) = archives.get(&key) {
                let shared = shared.clone();
                shared
                    .lock()
                    .map_err(|_| anyhow::anyhow!("SeqArc reader lock poisoned"))?
                    .requested
                    .insert(stream.clone());
                shared
            } else {
                let mut readers = Vec::with_capacity(files.len());
                for file in files {
                    let reader = SeqArcReader::open(file)?
                        .ok_or_else(|| anyhow::anyhow!("failed to open SeqArc input"))?;
                    anyhow::ensure!(
                        reader.preserve_ids(),
                        "SeqArc archive does not preserve identifiers"
                    );
                    anyhow::ensure!(
                        reader.preserve_quality(),
                        "SeqArc archive does not preserve quality scores"
                    );
                    readers.push(reader);
                }
                let shared = Arc::new(Mutex::new(SeqArcSharedReader {
                    readers,
                    current: 0,
                    requested: std::iter::once(stream.clone()).collect(),
                    served: HashSet::new(),
                    group: None,
                }));
                archives.insert(key, shared.clone());
                shared
            }
        };
        Ok(FastqReader {
            reader: Box::new(SeqArcRecordSource {
                reader: shared,
                stream,
            }),
            min_len: 0,
            max_len: 0,
        })
    }
}

/// A streaming reader for one SeqArc archive.
///
/// SeqArc records are exposed as complete logical groups so paired biological
/// reads and separate index reads remain synchronized and the archive is only
/// decompressed once.
pub struct SeqArcReader {
    inner: seqarc::AlignedRecordIter<BufReader<Box<dyn ReadSeek + Send + Sync>>>,
}

/// Build a selective-I/O mask for the requested stream names.
pub fn seqarc_selection_for_streams<'a, I>(streams: I) -> Result<seqarc::StreamSelection>
where
    I: IntoIterator<Item = &'a str>,
{
    let mut selection = seqarc::StreamSelection {
        biological: [false, false],
        index: [false, false],
    };
    for stream in streams {
        match stream {
            "R1" | "1" => selection.biological[0] = true,
            "R2" | "2" => selection.biological[1] = true,
            "I1" => selection.index[0] = true,
            "I2" => selection.index[1] = true,
            _ => bail!("unknown SeqArc stream '{stream}'"),
        }
    }
    anyhow::ensure!(
        selection.biological.iter().any(|selected| *selected)
            || selection.index.iter().any(|selected| *selected),
        "at least one SeqArc stream must be selected"
    );
    Ok(selection)
}

impl SeqArcReader {
    /// Open a SeqArc archive with the default decoding options.
    pub fn open(file: &File) -> Result<Option<Self>> {
        Self::open_with_selection(file, None)
    }

    /// Open a SeqArc archive while materializing only the selected streams.
    pub fn open_with_selection(
        file: &File,
        selection: Option<seqarc::StreamSelection>,
    ) -> Result<Option<Self>> {
        if !file.is_seqarc() {
            return Ok(None);
        }
        let input = file.open_seekable()?;
        let inner = seqarc::decompress_aligned_iter(
            BufReader::new(input),
            &DecompressionOptions {
                threads: 1,
                stream_selection: selection,
                ..Default::default()
            },
        )?;
        Ok(Some(Self { inner }))
    }

    /// Read the next complete logical record group.
    pub fn next_group(&mut self) -> Option<Result<AlignedRecordGroup>> {
        self.next().map(|result| result.map(|(_, group)| group))
    }

    /// Read a selected stream from the next logical record.
    pub fn next_stream_record(&mut self, stream: &str) -> Option<Result<fastq::Record>> {
        self.next_group().map(|result| {
            let group = result?;
            let record = match stream {
                "R1" | "1" => match group.sequence {
                    SequenceRecordGroup::Single(record) => record,
                    SequenceRecordGroup::Paired(record, _) => record,
                },
                "R2" | "2" => match group.sequence {
                    SequenceRecordGroup::Paired(_, record) => record,
                    SequenceRecordGroup::Single(_) => {
                        bail!("SeqArc archive is single-end; R2 is unavailable")
                    }
                },
                stream if stream.starts_with('I') => {
                    let index = stream[1..]
                        .parse::<usize>()
                        .map_err(|_| anyhow::anyhow!("invalid SeqArc index stream: {stream}"))?;
                    group
                        .index_reads
                        .into_iter()
                        .nth(
                            index
                                .checked_sub(1)
                                .ok_or_else(|| anyhow::anyhow!("index streams are 1-based"))?,
                        )
                        .ok_or_else(|| {
                            anyhow::anyhow!("SeqArc index stream {stream} is unavailable")
                        })?
                }
                stream => bail!("unknown SeqArc stream '{stream}' (expected R1, R2, I1, or I2)"),
            };
            fastq_record_from_seqarc(record)
        })
    }

    /// Layout of the biological records in the archive.
    pub fn layout(&self) -> seqarc::AlignmentLayout {
        self.inner.layout()
    }

    /// Number of separate index-read streams.
    pub fn index_streams(&self) -> usize {
        self.inner.index_streams()
    }

    /// Number of embedded index ranges already restored into mate sequences.
    pub fn embedded_index_count(&self) -> usize {
        self.inner.embedded_index_count()
    }

    /// Whether identifiers are available.
    pub fn preserve_ids(&self) -> bool {
        self.inner.preserve_ids()
    }

    /// Whether qualities are available.
    pub fn preserve_quality(&self) -> bool {
        self.inner.preserve_quality()
    }

    /// Validate that a named stream exists in this archive.
    pub fn validate_stream(&self, stream: &str) -> Result<()> {
        match stream {
            "R1" | "1" => Ok(()),
            "R2" | "2" if self.layout() == seqarc::AlignmentLayout::PairedEnd => Ok(()),
            "R2" | "2" => bail!("SeqArc archive is single-end; stream R2 is unavailable"),
            "I1" | "I2" => {
                let index = stream[1..].parse::<usize>()?;
                anyhow::ensure!(
                    index > 0 && index <= self.index_streams(),
                    "SeqArc archive provides {} separate index streams; {} was requested",
                    self.index_streams(),
                    stream
                );
                Ok(())
            }
            _ => bail!("unknown SeqArc stream '{stream}'"),
        }
    }
}

impl Iterator for SeqArcReader {
    type Item = Result<(u64, AlignedRecordGroup)>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|result| result.map_err(anyhow::Error::msg))
    }
}

/// Convert a SeqArc read record into the FASTQ record type used by seqspec.
pub fn fastq_record_from_seqarc(record: seqarc::ReadRecord) -> Result<fastq::Record> {
    let quality = record
        .quality
        .ok_or_else(|| anyhow::anyhow!("SeqArc archive does not preserve quality scores"))?;
    if quality.len() != record.sequence.len() {
        bail!("SeqArc sequence and quality lengths differ");
    }
    let definition = record.id;
    let name_end = definition
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .map_or(0, |name| name.len());
    let name = definition[..name_end].to_vec();
    let description = definition[name_end..]
        .iter()
        .skip_while(|byte| byte.is_ascii_whitespace())
        .copied()
        .collect::<Vec<_>>();
    Ok(fastq::Record::new(
        fastq::record::Definition::new(name, description),
        record.sequence,
        quality,
    ))
}

/// Convert a decoded SeqArc logical group to FASTQ records in R1, R2, I1, I2
/// order. Missing biological or index streams are omitted.
pub fn fastq_records_from_seqarc_group(group: AlignedRecordGroup) -> Result<Vec<fastq::Record>> {
    let mut records = Vec::new();
    match group.sequence {
        SequenceRecordGroup::Single(record) => records.push(fastq_record_from_seqarc(record)?),
        SequenceRecordGroup::Paired(first, second) => {
            records.push(fastq_record_from_seqarc(first)?);
            records.push(fastq_record_from_seqarc(second)?);
        }
    }
    for record in group.index_reads {
        records.push(fastq_record_from_seqarc(record)?);
    }
    Ok(records)
}

/// Convert one stream from a decoded SeqArc logical group to a FASTQ record.
pub fn fastq_record_from_seqarc_group_stream(
    group: &AlignedRecordGroup,
    stream: &str,
) -> Result<fastq::Record> {
    let record = match stream {
        "R1" | "1" => match &group.sequence {
            SequenceRecordGroup::Single(record) => record,
            SequenceRecordGroup::Paired(record, _) => record,
        },
        "R2" | "2" => match &group.sequence {
            SequenceRecordGroup::Paired(_, record) => record,
            SequenceRecordGroup::Single(_) => {
                bail!("SeqArc archive is single-end; R2 is unavailable")
            }
        },
        stream if stream.starts_with('I') => {
            let index = stream[1..]
                .parse::<usize>()
                .map_err(|_| anyhow::anyhow!("invalid SeqArc index stream: {stream}"))?;
            group
                .index_reads
                .get(
                    index
                        .checked_sub(1)
                        .ok_or_else(|| anyhow::anyhow!("index streams are 1-based"))?,
                )
                .ok_or_else(|| anyhow::anyhow!("SeqArc index stream {stream} is unavailable"))?
        }
        stream => bail!("unknown SeqArc stream '{stream}' (expected R1, R2, I1, or I2)"),
    };
    fastq_record_from_seqarc(record.clone())
}

impl FastqReader {
    pub fn new(reader: Box<dyn BufRead + Send>, min_len: u32, max_len: u32) -> Self {
        Self {
            reader: Box::new(TextRecordSource {
                reader: fastq::io::Reader::new(reader),
            }),
            min_len,
            max_len,
        }
    }

    pub fn read_record(&mut self, record: &mut fastq::Record) -> Result<usize> {
        let n = self.reader.read_record(record)?;
        if self.max_len > 0 {
            record.quality_scores_mut().truncate(self.max_len as usize);
            record.sequence_mut().truncate(self.max_len as usize);
        }
        if n != 0 && record.sequence().len() < self.min_len as usize {
            bail!(
                "Record sequence is too short: {} < {}",
                record.sequence().len(),
                self.min_len
            );
        }
        Ok(n)
    }

    fn read_record_unchecked(&mut self, record: &mut fastq::Record) -> Result<usize> {
        Ok(self.reader.read_record(record)?)
    }

    pub fn records(&mut self) -> FastqRecords<'_> {
        FastqRecords { inner: self }
    }
}

pub struct FastqRecords<'a> {
    inner: &'a mut FastqReader,
}

impl<'a> Iterator for FastqRecords<'a> {
    type Item = Result<fastq::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = fastq::Record::default();

        match self.inner.read_record(&mut buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(buf)),
            Err(e) => Some(Err(e)),
        }
    }
}

impl Read {
    /// Open the FASTQ files for reading, returning I/O failures to the caller.
    /// If the read has multiple FASTQ files, they are concatenated. If it has
    /// no FASTQ files, return `Ok(None)`.
    pub fn try_open(&self) -> Result<Option<FastqReader>> {
        self.try_open_with_context(&SeqArcContext::default())
    }

    /// Open this read using a shared SeqArc context.
    pub fn try_open_with_context(&self, context: &SeqArcContext) -> Result<Option<FastqReader>> {
        let files = self
            .files
            .clone()
            .unwrap_or_default()
            .into_iter()
            .filter(|file| file.filetype == "fastq" || file.is_seqarc())
            .collect::<Vec<_>>();
        if files.is_empty() {
            return Ok(None);
        }
        if files.iter().any(File::is_seqarc) {
            anyhow::ensure!(
                files.iter().all(File::is_seqarc),
                "a SeqArc input cannot be mixed with FASTQ files in one Read"
            );
            let stream = files[0]
                .url
                .split_once('#')
                .and_then(|(_, fragment)| fragment.split('&').next())
                .filter(|stream| !stream.is_empty())
                .unwrap_or("R1")
                .to_string();
            let mut reader = context.open_stream(&files, stream)?;
            reader.min_len = self.min_len;
            reader.max_len = self.max_len;
            return Ok(Some(reader));
        }
        let readers = files
            .into_iter()
            .map(|file| file.open())
            .collect::<Result<Vec<_>>>()?;
        let reader = multi_reader::MultiReader::new(readers.into_iter());
        Ok(Some(FastqReader::new(
            Box::new(BufReader::new(reader)),
            self.min_len,
            self.max_len,
        )))
    }

    /// Open the fastq files for reading, and return a fastq reader.
    /// If the read has multiple fastq files, they will be concatenated.
    /// If the read has no fastq files, return None.
    pub fn open(&self) -> Option<FastqReader> {
        self.try_open().expect("failed to open FASTQ files")
    }

    /// Get the actual length of the read by reading the first N record from the fastq file.
    pub fn infer_length(&self, n: usize) -> Result<std::ops::Range<usize>> {
        let mut reader = self
            .open()
            .ok_or_else(|| anyhow::anyhow!("No fastq files found."))?;
        let mut record = fastq::Record::default();
        let mut min_len = usize::MAX;
        let mut max_len = 0;
        for _ in 0..n {
            if reader.read_record_unchecked(&mut record)? == 0 {
                break;
            }
            let len = record.sequence().len();
            min_len = min_len.min(len);
            max_len = max_len.max(len);
        }
        anyhow::ensure!(
            min_len != usize::MAX,
            "No reads found while inferring read length"
        );
        Ok(min_len..max_len)
    }

    /// Check if the read is reverse.
    ///
    /// For `Strand::Unstranded` (used in long-read assays), this returns `false` (forward),
    /// which means the segment layout from the YAML is interpreted in forward orientation.
    /// The actual per-read orientation is detected dynamically during barcode extraction.
    pub fn is_reverse(&self) -> bool {
        match self.strand {
            Strand::Neg => true,
            Strand::Pos | Strand::Unstranded => false,
        }
    }

    pub(crate) fn get_segments<'a>(
        &'a self,
        region: &'a Region,
        truncate_by_length: bool,
    ) -> Option<SegmentInfo> {
        if !region.sequence_type.is_joined() {
            return None;
        }

        let subregions = region.subregions.iter();
        let subregions: Box<dyn Iterator<Item = _>> = if self.is_reverse() {
            Box::new(subregions.rev())
        } else {
            Box::new(subregions)
        };

        let segment_info = subregions
            .skip_while(|region| {
                let region = region.read().unwrap();
                let found =
                    region.region_type.is_sequencing_primer() && region.region_id == self.primer_id;
                !found
            }) // Skip until we find the primer region
            .skip(1)
            .map(|x| x.read().unwrap().deref().clone());
        let mut segment_info = SegmentInfo::new(segment_info, self.is_reverse());
        if truncate_by_length {
            segment_info = segment_info.truncate_max(self.max_len as usize);
        }
        if segment_info.len() > 0 {
            Some(segment_info)
        } else {
            None
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Strand {
    Pos,
    Neg,
    Unstranded,
}

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct File {
    pub file_id: String,
    pub filename: String,
    pub filetype: String,
    pub filesize: u64,
    pub url: String,
    pub urltype: UrlType,
    pub md5: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grouped_seqspec_expands_stream_definitions() {
        let yaml = r#"
- read_id:
    R1: atac-R1
    R2: atac-R2
    I1: atac-I1
  name: null
  modality: atac
  primer_id:
    R1: p1
    R2: p2
    I1: p3
  min_len:
    R1: 151
    R2: 151
    I1: 15
  max_len:
    R1: 151
    R2: 151
    I1: 15
  strand:
    R1: pos
    R2: neg
    I1: neg
  files:
  - file_id: reads.seqarc
    filename: reads.seqarc
    filetype: seqarc
    filesize: 0
    url: reads.seqarc
    urltype: local
    md5: ""
"#;
        let spec: SeqSpec = serde_yaml::from_str(yaml).unwrap();
        assert_eq!(spec.len(), 3);
        assert_eq!(
            spec.get("atac-R1").unwrap().files.as_ref().unwrap()[0].url,
            "reads.seqarc#R1"
        );
        assert_eq!(spec.get("atac-I1").unwrap().min_len, 15);
        assert!(spec.get("atac-R2").unwrap().files.as_ref().unwrap()[0].is_seqarc());
    }

    #[test]
    fn seqarc_file_extension_is_detected_and_fragment_is_ignored_for_opening() {
        let file = File {
            file_id: "reads.seqarc".into(),
            filename: "reads.seqarc".into(),
            filetype: "fastq".into(),
            filesize: 0,
            url: "reads.seqarc#R2".into(),
            urltype: UrlType::Local,
            md5: String::new(),
        };
        assert!(file.is_seqarc());
    }

    #[test]
    fn seqarc_file_entry_uses_stream_fragment() {
        let path = std::env::temp_dir().join(format!(
            "precellar-seqarc-test-{}-reads.seqarc",
            std::process::id()
        ));
        std::fs::write(&path, b"placeholder").unwrap();
        let file = File::from_seqarc(&path, "R1", false).unwrap();
        assert_eq!(file.filetype, "seqarc");
        assert!(file.url.ends_with("reads.seqarc#R1"));
        assert!(file.is_seqarc());
        assert!(File::from_seqarc(&path, "R3", false).is_err());
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn infer_length_ignores_stale_length_constraints() {
        let path = std::env::temp_dir().join(format!(
            "precellar-infer-length-{}-reads.fastq",
            std::process::id()
        ));
        std::fs::write(&path, b"@read\nACGT\n+\nIIII\n").unwrap();
        let file = File::from_fastq(&path, false).unwrap();
        let read = Read {
            read_id: "R1".into(),
            files: Some(vec![file]),
            min_len: 150,
            max_len: 150,
            ..Default::default()
        };
        assert_eq!(read.infer_length(10).unwrap(), 4..4);
        let _ = std::fs::remove_file(path);
    }
}

impl File {
    pub fn from_fastq<P: AsRef<Path>>(path: P, compute_md5: bool) -> Result<Self> {
        let file = std::fs::File::open(&path)?;
        let filename = path.as_ref().file_name().unwrap().to_str().unwrap();
        let md5 = if compute_md5 {
            crate::utils::md5sum(&path)?
        } else {
            "0".to_string()
        };
        Ok(File {
            file_id: filename.to_string(),
            filename: filename.to_string(),
            filetype: "fastq".to_string(),
            filesize: file.metadata()?.len(),
            url: path.as_ref().to_str().unwrap().to_string(),
            urltype: UrlType::Local,
            md5,
        })
    }

    /// Create a SeqArc file entry for one archive stream.
    pub fn from_seqarc<P: AsRef<Path>>(path: P, stream: &str, compute_md5: bool) -> Result<Self> {
        let path = path.as_ref();
        let file = std::fs::File::open(path)?;
        let filename = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| anyhow::anyhow!("invalid SeqArc filename: {path:?}"))?;
        let md5 = if compute_md5 {
            crate::utils::md5sum(path)?
        } else {
            "0".to_string()
        };
        anyhow::ensure!(
            matches!(stream, "R1" | "R2" | "I1" | "I2"),
            "invalid SeqArc stream '{stream}'; expected R1, R2, I1, or I2"
        );
        Ok(Self {
            file_id: filename.to_string(),
            filename: filename.to_string(),
            filetype: "seqarc".to_string(),
            filesize: file.metadata()?.len(),
            url: format!("{}#{stream}", path.display()),
            urltype: UrlType::Local,
            md5,
        })
    }

    /// Returns true for an explicitly declared or extension-detected SeqArc file.
    pub fn is_seqarc(&self) -> bool {
        self.filetype.eq_ignore_ascii_case("seqarc")
            || self
                .url
                .split('#')
                .next()
                .and_then(|url| Path::new(url).extension())
                .is_some_and(|extension| extension.eq_ignore_ascii_case("seqarc"))
    }

    pub(crate) fn normalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            let (path, fragment) = self.url.split_once('#').unwrap_or((&self.url, ""));
            let mut path = path.to_string();
            let normalized = crate::utils::normalize_path(work_dir, &mut path)?
                .to_str()
                .unwrap()
                .to_owned();
            self.url = if fragment.is_empty() {
                normalized
            } else {
                format!("{normalized}#{fragment}")
            };
        }
        Ok(())
    }

    pub(crate) fn unnormalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            let (path, fragment) = self.url.split_once('#').unwrap_or((&self.url, ""));
            let mut path = path.to_string();
            let normalized = crate::utils::unnormalize_path(work_dir, &mut path)?
                .to_str()
                .unwrap()
                .to_owned();
            self.url = if fragment.is_empty() {
                normalized
            } else {
                format!("{normalized}#{fragment}")
            };
        }
        Ok(())
    }

    /// Open the file for reading.
    /// If the file is remote, it will be downloaded to the cache directory.
    /// If the file is local, it will be opened directly.
    /// The base_dir is used to resolve relative paths.
    pub fn open(&self) -> Result<Box<dyn std::io::Read + Send + Sync>> {
        let source_url = self.url.split('#').next().unwrap_or(&self.url);
        match self.urltype {
            UrlType::Local => Ok(Box::new(crate::utils::open_file(source_url)?)),
            _ => {
                let cache_dir = home::home_dir().unwrap().join(".cache/seqspec");
                let downloader = Downloader::new(Some(cache_dir))?;
                let file_path =
                    downloader.retrieve(source_url, Some(&self.filename), None, false)?;
                Ok(Box::new(crate::utils::open_file(file_path)?))
            }
        }
    }

    /// Open a local or cached file with seek support.
    pub fn open_seekable(&self) -> Result<Box<dyn ReadSeek + Send + Sync>> {
        let source_url = self.url.split('#').next().unwrap_or(&self.url);
        match self.urltype {
            UrlType::Local => Ok(Box::new(std::fs::File::open(source_url)?)),
            _ => {
                let cache_dir = home::home_dir().unwrap().join(".cache/seqspec");
                let downloader = Downloader::new(Some(cache_dir))?;
                let file_path =
                    downloader.retrieve(source_url, Some(&self.filename), None, false)?;
                Ok(Box::new(std::fs::File::open(file_path)?))
            }
        }
    }
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum UrlType {
    Local,
    Ftp,
    Http,
    Https,
}
