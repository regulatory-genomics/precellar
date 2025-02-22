use crate::region::{Region, SequenceType};
use crate::{Modality, RegionType};

use anyhow::Result;
use cached_path::Cache;
use indexmap::{IndexMap, IndexSet};
use noodles::fastq;
use serde::{Deserialize, Serialize, Serializer};
use std::ops::{Deref, DerefMut};
use std::path::Path;
use std::sync::{Arc, RwLock};
use std::{
    io::{BufRead, BufReader},
    ops::Range,
};

/// Specification of a sequencing library.
#[derive(Debug, Clone, PartialEq)]
pub struct SeqSpec(IndexMap<String, Read>);

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
        let reads = Vec::<Read>::deserialize(deserializer)?;
        let mut sequences = IndexMap::new();
        for read in reads {
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

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Read {
    pub read_id: String,
    pub name: Option<String>,
    pub modality: Modality,
    pub primer_id: String,
    pub min_len: u32,
    pub max_len: u32,
    pub strand: Strand,
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
    reader: fastq::Reader<Box<dyn BufRead>>,
    min_len: u32,
    max_len: u32,
}

impl FastqReader {
    pub fn new(reader: Box<dyn BufRead>, min_len: u32, max_len: u32) -> Self {
        Self {
            reader: fastq::Reader::new(reader),
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
        Ok(n)
    }

    pub fn records(&mut self) -> FastqRecords {
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
    /// Open the fastq files for reading, and return a fastq reader.
    /// If the read has multiple fastq files, they will be concatenated.
    /// If the read has no fastq files, return None.
    pub fn open(&self) -> Option<FastqReader> {
        let files = self
            .files
            .clone()
            .unwrap_or_default()
            .into_iter()
            .filter(|file| file.filetype == "fastq")
            .collect::<Vec<_>>();
        if files.is_empty() {
            return None;
        }
        let reader =
            multi_reader::MultiReader::new(files.into_iter().map(move |file| file.open().unwrap()));
        Some(FastqReader::new(
            Box::new(BufReader::new(reader)),
            self.min_len,
            self.max_len,
        ))
    }

    /// Get the actual length of the read by reading the first record from the fastq file.
    pub fn actual_len(&self) -> Result<usize> {
        let mut reader = self.open().expect("No fastq files found.").reader;
        let mut record = fastq::Record::default();
        reader.read_record(&mut record)?;
        Ok(record.sequence().len())
    }

    /// Check if the read is reverse.
    pub fn is_reverse(&self) -> bool {
        match self.strand {
            Strand::Neg => true,
            Strand::Pos => false,
        }
    }

    pub(crate) fn get_segments<'a>(&'a self, region: &'a Region) -> Option<SegmentInfo> {
        if region.sequence_type != SequenceType::Joined {
            return None;
        }

        let mut found_primer = false;
        let subregions = region.subregions.iter();
        let subregions: Box<dyn Iterator<Item = _>> = if self.is_reverse() {
            Box::new(subregions.rev())
        } else {
            Box::new(subregions)
        };

        let result = self.get_read_span(
            subregions
                .skip_while(|region| {
                    let region = region.read().unwrap();
                    let found = region.region_type.is_sequencing_primer()
                        && region.region_id == self.primer_id;
                    if found {
                        found_primer = true;
                    }
                    !found
                })
                .skip(1),
        );

        if found_primer {
            Some(result)
        } else {
            None
        }
    }

    /// Helper function to get the region index for a read.
    fn get_read_span<'a, I>(&self, mut regions: I) -> SegmentInfo
    where
        I: Iterator<Item = &'a Arc<RwLock<Region>>>,
    {
        let mut segments = Vec::new();
        let read_len = self.max_len;
        let mut cur_pos = 0;
        let mut readlen_info = ReadSpan::Covered;
        while let Some(region) = regions.next() {
            let mut region = region.read().unwrap();
            let mut region_id = region.region_id.clone();
            let mut region_type = SegmentType::R(region.region_type);
            if region.is_fixed_length() {
                // Fixed-length region
                let end = cur_pos + region.min_len;
                if end >= read_len {
                    segments.push(SegmentInfoElem::new(
                        region_id,
                        region_type,
                        cur_pos..read_len,
                    ));
                    if end > read_len {
                        readlen_info = ReadSpan::NotEnough;
                    }
                    break;
                } else {
                    segments.push(SegmentInfoElem::new(region_id, region_type, cur_pos..end));
                    cur_pos = end;
                }
            } else {
                // Variable-length region
                if let Some(nucl) = region.region_type.is_poly_nucl() {
                    if let Some(next_region) = regions.next() {
                        region = next_region.read().unwrap();
                        region_id = region.region_id.clone();
                        region_type = SegmentType::PolyNucl((nucl, region.region_type));
                    }
                }
                if cur_pos + region.min_len >= read_len {
                    // Variable-length region and read is shorter
                    segments.push(SegmentInfoElem::new(
                        region_id,
                        region_type,
                        cur_pos..read_len,
                    ));
                    readlen_info = ReadSpan::Span((read_len - cur_pos) as usize);
                    break;
                } else if cur_pos + region.max_len < read_len {
                    // Variable-length region and read is longer than max length
                    segments.push(SegmentInfoElem::new(
                        region_id,
                        region_type,
                        cur_pos..cur_pos + region.max_len,
                    ));
                    if let Some(next_region) = regions.next() {
                        let next_region = next_region.read().unwrap();
                        readlen_info = ReadSpan::ReadThrough(next_region.region_id.clone());
                    }
                    break;
                } else {
                    // Variable-length region and read is within the length range
                    segments.push(SegmentInfoElem::new(
                        region_id,
                        region_type,
                        cur_pos..read_len,
                    ));
                    if let Some(next_region) = regions.next() {
                        let next_region = next_region.read().unwrap();
                        readlen_info = ReadSpan::MayReadThrough(next_region.region_id.clone());
                    }
                    break;
                }
            }
        }
        SegmentInfo {
            segments,
            readlen_info,
        }
    }
}

/// A read may be divided into multiple segments, each corresponding to a region
/// with biological meaning. SegmentInfo stores the information about the regions
/// that the read is divided into.
#[derive(Debug, Clone)]
pub struct SegmentInfo {
    pub segments: Vec<SegmentInfoElem>,
    pub readlen_info: ReadSpan,
}

#[derive(Clone)]
pub enum SegmentType {
    R(RegionType),
    PolyNucl((u8, RegionType)), // (nucleotide, region_type)
}

impl std::fmt::Debug for SegmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SegmentType::R(region_type) => write!(f, "{:?}", region_type),
            SegmentType::PolyNucl((nucl, region_type)) => {
                write!(f, "poly{}+{:?}", *nucl as char, region_type)
            }
        }
    }
}

impl SegmentType {
    pub fn is_barcode(&self) -> bool {
        match self {
            SegmentType::R(region_type) => region_type.is_barcode(),
            _ => false,
        }
    }

    pub fn is_umi(&self) -> bool {
        match self {
            SegmentType::R(region_type) => region_type.is_umi(),
            _ => false,
        }
    }

    pub fn is_target(&self) -> bool {
        match self {
            SegmentType::R(region_type) => region_type.is_target(),
            SegmentType::PolyNucl((_, region_type)) => region_type.is_target(),
        }
    }

    pub fn poly_nucl(&self) -> Option<u8> {
        match self {
            SegmentType::PolyNucl((nucl, _)) => Some(*nucl),
            _ => None,
        }
    }
}

/// Information about a segment in a read.
#[derive(Debug, Clone)]
pub struct SegmentInfoElem {
    pub region_id: String,
    pub region_type: SegmentType,
    pub range: Range<u32>,
}

impl SegmentInfoElem {
    pub fn new(region_id: String, region_type: SegmentType, range: Range<u32>) -> Self {
        Self {
            region_id,
            region_type: region_type,
            range,
        }
    }
}

/// Information about the region index for a read.
#[derive(Debug, Clone)]
pub enum ReadSpan {
    Covered,                // The read is fully contained within the target region
    Span(usize),            // The read spans the target region
    NotEnough,              // Read is too short to reach the target region
    ReadThrough(String),    // Read is longer than the target region
    MayReadThrough(String), // Read may be longer than the target region
}

#[derive(Deserialize, Serialize, Debug, Copy, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Strand {
    Pos,
    Neg,
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

    pub(crate) fn normalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::normalize_path(work_dir, &mut self.url)?
                .to_str()
                .unwrap()
                .to_owned();
        }
        Ok(())
    }

    pub(crate) fn unnormalize_path<P: AsRef<Path>>(&mut self, work_dir: P) -> Result<()> {
        if self.urltype == UrlType::Local {
            self.url = crate::utils::unnormalize_path(work_dir, &mut self.url)?
                .to_str()
                .unwrap()
                .to_owned();
        }
        Ok(())
    }

    /// Open the file for reading.
    /// If the file is remote, it will be downloaded to the cache directory.
    /// If the file is local, it will be opened directly.
    /// The base_dir is used to resolve relative paths.
    pub fn open(&self) -> Result<Box<dyn std::io::Read>> {
        match self.urltype {
            UrlType::Local => Ok(Box::new(crate::utils::open_file(&self.url)?)),
            _ => {
                let mut cache = Cache::new().unwrap();
                cache.dir = home::home_dir().unwrap().join(".cache/seqspec");
                let file = cache.cached_path(&self.url).unwrap();
                Ok(Box::new(crate::utils::open_file(file)?))
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

pub(crate) struct ReadValidator<'a> {
    region: &'a Region,
    range: Option<Range<usize>>,
    n_total: usize,
    n_matched: usize,
    onlist: Option<IndexSet<Vec<u8>>>,
    strand: Strand,
    tolerance: f64,
}

impl<'a> ReadValidator<'a> {
    pub fn id(&self) -> &str {
        &self.region.region_id
    }

    pub fn sequence_type(&self) -> SequenceType {
        self.region.sequence_type
    }

    pub fn new(region: &'a Region) -> Result<Self> {
        let onlist = if let Some(onlist) = &region.onlist {
            Some(onlist.read()?)
        } else {
            None
        };
        Ok(Self {
            region,
            range: None,
            n_total: 0,
            n_matched: 0,
            onlist,
            strand: Strand::Pos,
            tolerance: 0.2,
        })
    }

    pub fn with_range(mut self, range: Range<usize>) -> Self {
        self.range = Some(range);
        self
    }

    pub fn with_strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        if !(0.0..=1.0).contains(&tolerance) {
            panic!("Tolerance must be between 0 and 1");
        }
        self.tolerance = tolerance;
        self
    }

    pub fn frac_matched(&self) -> f64 {
        self.n_matched as f64 / self.n_total as f64
    }

    /// Validate a sequence. Return true if the sequence is valid.
    /// For fixed sequences, the sequence is checked against the fixed sequence with
    /// certain tolerance (default is 1 mismatches).
    /// Note that for onlist, we check for the exact match.
    pub fn validate(&mut self, seq: &[u8]) -> ValidateResult {
        self.n_total += 1;
        let seq = if let Some(range) = &self.range {
            if range.end > seq.len() {
                // If range end is out of bounds, use the full sequence
                eprintln!("Warning: Range end {} exceeds sequence length {}", range.end, seq.len());
                seq
            } else {
                &seq[range.clone()]
            }
        } else {
            seq
        };
        let seq = if self.strand == Strand::Neg {
            crate::utils::rev_compl(seq)
        } else {
            seq.to_vec()
        };
        let len = seq.len();
        if len < self.region.min_len as usize {
            ValidateResult::TooShort((seq.len(), self.region.min_len as usize))
        } else if seq.len() > self.region.max_len as usize {
            ValidateResult::TooLong((seq.len(), self.region.max_len as usize))
        } else if self.sequence_type() == SequenceType::Fixed {
            let d = hamming::distance(&seq, self.region.sequence.as_bytes()) as f64;
            if d <= self.tolerance * len as f64 {
                self.n_matched += 1;
                ValidateResult::Valid
            } else {
                ValidateResult::MisMatch(d as usize)
            }
        } else if let Some(onlist) = &self.onlist {
            if onlist.contains(&seq) {
                self.n_matched += 1;
                ValidateResult::Valid
            } else {
                ValidateResult::OnlistFail
            }
        } else {
            self.n_matched += 1;
            ValidateResult::Valid
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) enum ValidateResult {
    TooShort((usize, usize)),
    TooLong((usize, usize)),
    OnlistFail,
    MisMatch(usize),
    Valid,
}

impl std::fmt::Display for ValidateResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ValidateResult::TooShort((n, min_len)) => {
                write!(f, "Sequence too short: {} < {}", n, min_len)
            }
            ValidateResult::TooLong((n, max_len)) => {
                write!(f, "Sequence too long: {} > {}", n, max_len)
            }
            ValidateResult::OnlistFail => {
                write!(f, "Not on the onlist")
            }
            ValidateResult::MisMatch(n) => {
                write!(f, "Mismatch: {}", n)
            }
            ValidateResult::Valid => {
                write!(f, "Valid")
            }
        }
    }
}
