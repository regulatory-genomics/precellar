mod read;
mod region;
pub mod utils;

use log::{debug, warn};
use noodles::fastq;
use read::ReadValidator;
pub use read::{FastqReader, File, Read, SegmentInfo, SegmentInfoElem, Strand, UrlType};
use read::{ReadSpan, ValidateResult};
use region::LibSpec;
pub use region::{Onlist, Region, RegionId, RegionType, SequenceType};

use anyhow::{anyhow, bail, Context, Result};
use serde::{Deserialize, Deserializer, Serialize};
use serde_yaml::{self, Value};
use std::path::Path;
use std::{
    fs,
    path::PathBuf,
    str::FromStr,
    sync::{Arc, RwLock},
};
use utils::{create_file, Compression};

#[derive(Deserialize, Serialize, Debug, Clone, PartialEq)]
pub struct Assay {
    #[serde(skip, default)]
    pub file: Option<PathBuf>,
    pub seqspec_version: String,
    pub assay_id: String,
    pub name: String,
    pub doi: String,
    pub date: String,
    pub description: String,
    pub modalities: Vec<Modality>,
    pub lib_struct: String,
    pub library_protocol: LibraryProtocol,
    pub library_kit: LibraryKit,
    pub sequence_protocol: SequenceProtocol,
    pub sequence_kit: SequenceKit,
    pub sequence_spec: read::SeqSpec,
    pub library_spec: LibSpec,
}

impl Assay {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let yaml_str = fs::read_to_string(&path)
            .with_context(|| format!("Failed to read file: {:?}", path.as_ref()))?;
        let mut assay: Assay = serde_yaml::from_str(&yaml_str)?;
        assay.file = Some(path.as_ref().to_path_buf());
        assay.normalize_all_paths();
        Ok(assay)
    }

    pub fn from_url(url: &str) -> Result<Self> {
        let yaml_str = reqwest::blocking::get(url)?.text()?;
        Ok(serde_yaml::from_str(&yaml_str)?)
    }

    /// Return the modality if there is only one modality.
    pub fn modality(&self) -> Result<Modality> {
        if self.modalities.len() == 1 {
            Ok(self.modalities[0].clone())
        } else {
            bail!(
                "Multiple modalities found: {}",
                self.modalities
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }

    /// Normalize all paths in the sequence spec.
    /// This is used to bring relative paths to absolute paths.
    pub fn normalize_all_paths(&mut self) {
        if let Some(file) = &self.file {
            let base_dir = file.parent().unwrap();
            self.sequence_spec.values_mut().for_each(|read| {
                if let Some(files) = &mut read.files {
                    files.iter_mut().for_each(|file| {
                        if let Err(e) = file.normalize_path(base_dir) {
                            warn!("{}", e);
                        }
                    });
                }
            });
            self.library_spec.regions().for_each(|region| {
                if let Some(onlist) = &mut region.write().unwrap().onlist {
                    if let Err(e) = onlist.normalize_path(base_dir) {
                        warn!("{}", e);
                    }
                }
            });
        }
    }

    /// Unnormalize all paths in the sequence spec.
    /// This is used to bring absolute paths to relative paths.
    pub fn unnormalize_all_paths<P: AsRef<Path>>(&mut self, base_dir: P) {
        self.sequence_spec.values_mut().for_each(|read| {
            if let Some(files) = &mut read.files {
                files.iter_mut().for_each(|file| {
                    if let Err(e) = file.unnormalize_path(base_dir.as_ref()) {
                        warn!("Failed to unnormalize path: {}", e);
                    }
                });
            }
        });
        self.library_spec.regions().for_each(|region| {
            if let Some(onlist) = &mut region.write().unwrap().onlist {
                if let Err(e) = onlist.unnormalize_path(base_dir.as_ref()) {
                    warn!("Failed to unnormalize path: {}", e);
                }
            }
        });
    }

    /// Add default Illumina reads to the sequence spec.
    pub fn add_illumina_reads(
        &mut self,
        modality: Modality,
        read_len: usize,
        forward_strand_workflow: bool,
    ) -> Result<()> {
        fn advance_until(
            iterator: &mut std::slice::Iter<'_, Arc<RwLock<Region>>>,
            f: fn(&Region) -> bool,
        ) -> Option<(Arc<RwLock<Region>>, Vec<Arc<RwLock<Region>>>)> {
            let mut regions = Vec::new();
            for next_region in iterator.by_ref() {
                let r = next_region.read().unwrap();
                if f(&r) {
                    return Some((next_region.clone(), regions));
                } else {
                    regions.push(next_region.clone());
                }
            }
            None
        }

        fn get_length(regions: &[Arc<RwLock<Region>>], reverse: bool) -> usize {
            if reverse {
                regions
                    .iter()
                    .skip_while(|region| {
                        region.read().unwrap().sequence_type == SequenceType::Fixed
                    })
                    .map(|region| region.read().unwrap().len().unwrap() as usize)
                    .sum()
            } else {
                regions
                    .iter()
                    .rev()
                    .skip_while(|region| {
                        region.read().unwrap().sequence_type == SequenceType::Fixed
                    })
                    .map(|region| region.read().unwrap().len().unwrap() as usize)
                    .sum()
            }
        }

        fn is_read1(region: &Region) -> bool {
            region.region_type == RegionType::NexteraRead1
                || region.region_type == RegionType::TruseqRead1
        }

        fn is_read2(region: &Region) -> bool {
            region.region_type == RegionType::NexteraRead2
                || region.region_type == RegionType::TruseqRead2
        }

        fn is_p5(region: &Region) -> bool {
            region.region_type == RegionType::IlluminaP5
        }

        fn is_p7(region: &Region) -> bool {
            region.region_type == RegionType::IlluminaP7
        }

        self.delete_all_reads(modality);
        let regions = self
            .library_spec
            .get_modality(&modality)
            .ok_or_else(|| anyhow!("Modality not found: {:?}", modality))?
            .clone();
        let regions = regions.read().unwrap();
        let mut sub_regions = regions.subregions.iter();
        while let Some(current_region) = sub_regions.next() {
            let current_region = current_region.read().unwrap();
            if is_p5(&current_region) {
                if let Some((next_region, acc)) = advance_until(&mut sub_regions, is_read1) {
                    let next_region = next_region.read().unwrap();
                    self.update_read::<PathBuf>(
                        &format!("{}-R1", modality),
                        Some(modality),
                        Some(&next_region.region_id),
                        Some(false),
                        None,
                        Some(read_len),
                        Some(read_len),
                        false,
                    )?;
                    if forward_strand_workflow {
                        let acc_len = get_length(acc.as_slice(), false);
                        if acc_len > 0 {
                            self.update_read::<PathBuf>(
                                &format!("{}-I2", modality),
                                Some(modality),
                                Some(&current_region.region_id),
                                Some(false),
                                None,
                                Some(acc_len),
                                Some(acc_len),
                                false,
                            )?;
                        }
                    } else {
                        let acc_len = get_length(acc.as_slice(), true);
                        if acc_len > 0 {
                            self.update_read::<PathBuf>(
                                &format!("{}-I2", modality),
                                Some(modality),
                                Some(&next_region.region_id),
                                Some(true),
                                None,
                                Some(acc_len),
                                Some(acc_len),
                                false,
                            )?;
                        }
                    }
                }
            } else if is_read2(&current_region) {
                if let Some((_, acc)) = advance_until(&mut sub_regions, is_p7) {
                    let acc_len = get_length(acc.as_slice(), false);
                    self.update_read::<PathBuf>(
                        &format!("{}-R2", modality),
                        Some(modality),
                        Some(&current_region.region_id),
                        Some(true),
                        None,
                        Some(read_len),
                        Some(read_len),
                        false,
                    )?;
                    if acc_len > 0 {
                        self.update_read::<PathBuf>(
                            &format!("{}-I1", modality),
                            Some(modality),
                            Some(&current_region.region_id),
                            Some(false),
                            None,
                            Some(acc_len),
                            Some(acc_len),
                            false,
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Update read information. If the read does not exist, it will be created.
    pub fn update_read<P: AsRef<Path>>(
        &mut self,
        read_id: &str,
        modality: Option<Modality>,
        primer_id: Option<&str>,
        is_reverse: Option<bool>,
        fastqs: Option<&[P]>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        compute_md5: bool,
    ) -> Result<()> {
        let mut read_buffer = if let Some(r) = self.sequence_spec.get(read_id) {
            r.clone()
        } else {
            assert!(
                modality.is_some(),
                "modality must be provided for a new read"
            );
            assert!(
                primer_id.is_some(),
                "primer_id must be provided for a new read"
            );
            assert!(
                is_reverse.is_some(),
                "is_reverse must be provided for a new read"
            );
            Read {
                read_id: read_id.to_string(),
                ..Default::default()
            }
        };

        if let Some(rev) = is_reverse {
            read_buffer.strand = if rev { Strand::Neg } else { Strand::Pos };
        }
        if let Some(modality) = modality {
            read_buffer.modality = modality;
        }
        if let Some(primer_id) = primer_id {
            read_buffer.primer_id = primer_id.to_string();
        }
        if let Some(fastq) = fastqs {
            read_buffer.files = Some(
                fastq
                    .iter()
                    .map(|f| File::from_fastq(f, compute_md5))
                    .collect::<Result<Vec<File>>>()?,
            );
        }

        if (min_len.is_none() || max_len.is_none()) && read_buffer.files.is_some() {
            let len = read_buffer.actual_len()?;
            read_buffer.min_len = min_len.unwrap_or(len) as u32;
            read_buffer.max_len = max_len.unwrap_or(len) as u32;
        } else {
            read_buffer.min_len = min_len.unwrap() as u32;
            read_buffer.max_len = max_len.unwrap() as u32;
        }

        self.verify(&read_buffer)?;
        if let Some(region) = self.library_spec.get(&read_buffer.primer_id) {
            if !region.read().unwrap().region_type.is_sequencing_primer() {
                warn!(
                    "primer_id '{}' is not a sequencing primer (type: {:?})",
                    read_buffer.primer_id,
                    region.read().unwrap().region_type
                );
            }
        }
        self.sequence_spec.insert(read_id.to_string(), read_buffer);

        Ok(())
    }

    pub fn delete_read(&mut self, read_id: &str) -> Option<Read> {
        self.sequence_spec.shift_remove(read_id)
    }

    pub fn delete_all_reads(&mut self, modality: Modality) {
        let ids: Vec<_> = self
            .iter_reads(modality)
            .map(|read| read.read_id.to_string())
            .collect();
        ids.into_iter().for_each(|id| {
            self.delete_read(&id);
        });
    }

    /// Get the index of atomic regions of each read in the sequence spec.
    pub fn get_segments_by_modality(
        &self,
        modality: Modality,
    ) -> impl Iterator<Item = (&Read, SegmentInfo)> {
        self.sequence_spec.values().filter_map(move |read| {
            if read.modality == modality {
                let parent_region_index = self
                    .get_segments(&read.read_id)
                    .unwrap_or_else(|| panic!("Cannot find index for Read: {}", read.read_id));
                Some((read, parent_region_index))
            } else {
                None
            }
        })
    }

    /// Get atomic regions of a read in the sequence spec.
    pub fn get_segments(&self, read_id: &str) -> Option<SegmentInfo> {
        let read = self.sequence_spec.get(read_id)?;
        let library_parent_region = self.library_spec.get_parent(&read.primer_id)?;
        read.get_segments(&library_parent_region.read().unwrap())
    }

    pub fn iter_reads(&self, modality: Modality) -> impl Iterator<Item = &Read> {
        self.sequence_spec
            .values()
            .filter(move |read| read.modality == modality)
    }

    /// Validate reads in the sequence spec based on the fixed sequences and
    /// save the valid and invalid reads to different files.
    pub fn validate<P: AsRef<Path>>(
        &self,
        modality: Modality,
        dir: P,
        tolerance: f64,
    ) -> Result<()> {
        fs::create_dir_all(&dir)?;
        let (mut readers, reads): (Vec<_>, Vec<_>) = self
            .iter_reads(modality)
            .flat_map(|read| {
                let reader = read.open()?;
                let region = self
                    .library_spec
                    .get_parent(&read.primer_id)
                    .ok_or_else(|| anyhow!("Primer not found: {}", read.primer_id))
                    .unwrap();
                let index = read.get_segments(&region.read().unwrap())?;
                let regions: Vec<_> = index
                    .segments
                    .iter()
                    .map(|elem_info| {
                        let region = self.library_spec.get(&elem_info.region_id).unwrap();
                        (region.read().unwrap(), elem_info.range.clone())
                    })
                    .collect();
                Some((reader, (read, regions)))
            })
            .collect();

        let mut validators = reads
            .iter()
            .map(|(read, regions)| {
                regions
                    .iter()
                    .map(|(region, range)| {
                        Ok(ReadValidator::new(region)?
                            .with_range(range.start as usize..range.end as usize)
                            .with_strand(read.strand)
                            .with_tolerance(tolerance))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;

        let mut outputs: Vec<_> = reads
            .iter()
            .map(|(read, _)| {
                let output_valid = dir.as_ref().join(format!("{}.fq.zst", read.read_id));
                let output_valid = create_file(output_valid, Some(Compression::Zstd), Some(9), 8)?;
                let output_valid = fastq::io::Writer::new(output_valid);
                let output_other = dir
                    .as_ref()
                    .join(format!("Invalid_{}.fq.zst", read.read_id));
                let output_other = create_file(output_other, Some(Compression::Zstd), Some(9), 8)?;
                let output_other = fastq::io::Writer::new(output_other);

                anyhow::Ok((output_valid, output_other))
            })
            .collect::<Result<_>>()?;

        loop {
            let records: Vec<_> = readers
                .iter_mut()
                .flat_map(|reader| {
                    let mut record = fastq::Record::default();
                    if reader.read_record(&mut record).unwrap() == 0 {
                        None
                    } else {
                        Some(record)
                    }
                })
                .collect();

            if records.len() != readers.len() {
                break;
            }

            let valid = records
                .iter()
                .zip(validators.iter_mut())
                .all(|(record, validators)| {
                    validators.iter_mut().all(|validator| {
                        matches!(
                            validator.validate(record.sequence()),
                            ValidateResult::OnlistFail | ValidateResult::Valid
                        )
                    })
                });

            records.iter().zip(outputs.iter_mut()).try_for_each(
                |(record, (output_valid, output_other))| {
                    if valid {
                        output_valid.write_record(record)?;
                    } else {
                        output_other.write_record(record)?;
                    }
                    anyhow::Ok(())
                },
            )?;
        }
        Ok(())
    }

    /// Verify reads in the sequence spec.
    fn verify(&self, read: &Read) -> Result<()> {
        let region = self
            .library_spec
            .get_parent(&read.primer_id)
            .ok_or_else(|| anyhow!("Primer not found: {}", read.primer_id))?;
        // Check if the primer exists
        if let Some(index) = read.get_segments(&region.read().unwrap()) {
            match index.readlen_info {
                ReadSpan::Covered | ReadSpan::Span(_) => {}
                ReadSpan::NotEnough => {
                    warn!("'{}' does not cover the region", read.read_id);
                }
                ReadSpan::MayReadThrough(id) => {
                    debug!(
                        "'{}' may read through and contain sequences from: '{}'",
                        read.read_id, id
                    );
                }
                ReadSpan::ReadThrough(id) => {
                    warn!("'{}' length exceeds maximum length of the variable-length region (insertion), \
                    truncating the reads to the maximum length of the region. \
                    Read reads through and contains sequences from: '{}'.
                    If this is not the desired behavior, please adjust the region lengths.", read.read_id, id);
                }
            }
        
            if let Some(mut reader) = read.open() {
                let regions = index
                    .segments
                    .iter()
                    .map(|info| {
                        let region = self.library_spec.get(&info.region_id).unwrap();
                        (region.read().unwrap(), &info.range)
                    })
                    .collect::<Vec<_>>();
                /* */
                let mut validators = regions
                    .iter()
                    .map(|(region, range)| {
                        Ok(ReadValidator::new(region)?
                            .with_range(range.start as usize..range.end as usize)
                            .with_strand(read.strand))
                    })
                    .collect::<Result<Vec<_>>>()?;

                reader.records().take(1000).try_for_each(|record| {
                    let record = record?;
                    for validator in &mut validators {
                        let result = validator.validate(record.sequence());
                        match result {
                            ValidateResult::TooLong(_) | ValidateResult::TooShort(_) => {
                                bail!("{}: {}", read.read_id, result);
                            }
                            _ => {}
                        }
                    }
                    anyhow::Ok(())
                })?;

                for validator in validators {
                    let percent_matched = validator.frac_matched() * 100.0;
                    if percent_matched < 50.0 {
                        warn!(
                            "Read '{}' has low percentage of matched records for region '{}'. \
                        Percentage of matched records: {:.2}%",
                            read.read_id,
                            validator.id(),
                            percent_matched
                        );
                    }
                }
            }
            
        }
        Ok(())
    }
}

#[derive(Serialize, Debug, Copy, Clone, PartialEq, Eq, Hash)]
#[serde(rename_all = "lowercase")]
pub enum Modality {
    DNA,
    RNA,
    Tag,
    Protein,
    ATAC,
    Crispr,
}

impl Modality {
    pub fn from_str_case_insensitive(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "RNA" => Some(Modality::RNA),
            "ATAC" => Some(Modality::ATAC),
            "PROTEIN" => Some(Modality::Protein),
            "DNA" => Some(Modality::DNA),
            "TAG" => Some(Modality::Tag),
            "CRISPR" => Some(Modality::Crispr),
            _ => None,
        }
    }
}
impl<'de> Deserialize<'de> for Modality {
    fn deserialize<D>(deserializer: D) -> Result<Modality, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Modality::from_str(&s).map_err(serde::de::Error::custom),
            _ => Err(serde::de::Error::custom(format!(
                "invalid value: {:?}",
                value
            ))),
        }
    }
}

impl FromStr for Modality {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "dna" => Ok(Modality::DNA),
            "rna" => Ok(Modality::RNA),
            "tag" => Ok(Modality::Tag),
            "protein" => Ok(Modality::Protein),
            "atac" => Ok(Modality::ATAC),
            "crispr" => Ok(Modality::Crispr),
            _ => bail!("Invalid modality: {}", s),
        }
    }
}

impl std::fmt::Display for Modality {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Modality::DNA => "dna",
            Modality::RNA => "rna",
            Modality::Tag => "tag",
            Modality::Protein => "protein",
            Modality::ATAC => "atac",
            Modality::Crispr => "crispr",
        };
        write!(f, "{}", s)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum LibraryProtocol {
    Standard(String),
    Custom(Vec<ProtocolItem>),
}

#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
pub struct ProtocolItem {
    pub protocol_id: String,
    pub name: Option<String>,
    pub modality: Modality,
}

impl<'de> Deserialize<'de> for LibraryProtocol {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(LibraryProtocol::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq
                    .into_iter()
                    .map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(LibraryProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!(
                "invalid value: {:?}",
                value
            ))),
        }
    }
}

impl Serialize for LibraryProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        match self {
            LibraryProtocol::Standard(s) => s.serialize(serializer),
            LibraryProtocol::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum LibraryKit {
    Standard(String),
    Custom(Vec<KitItem>),
}

#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
pub struct KitItem {
    pub kit_id: String,
    pub name: Option<String>,
    pub modality: Modality,
}

impl<'de> Deserialize<'de> for LibraryKit {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(LibraryKit::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq
                    .into_iter()
                    .map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(LibraryKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!(
                "invalid value: {:?}",
                value
            ))),
        }
    }
}

impl Serialize for LibraryKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        match self {
            LibraryKit::Standard(s) => s.serialize(serializer),
            LibraryKit::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SequenceProtocol {
    Custom(Vec<ProtocolItem>),
    Standard(String),
}

impl<'de> Deserialize<'de> for SequenceProtocol {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(SequenceProtocol::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq
                    .into_iter()
                    .map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(SequenceProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!(
                "invalid value: {:?}",
                value
            ))),
        }
    }
}

impl Serialize for SequenceProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        match self {
            SequenceProtocol::Standard(s) => s.serialize(serializer),
            SequenceProtocol::Custom(items) => items.serialize(serializer),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum SequenceKit {
    Standard(String),
    Custom(Vec<KitItem>),
}

impl<'de> Deserialize<'de> for SequenceKit {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(SequenceKit::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq
                    .into_iter()
                    .map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(SequenceKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!(
                "invalid value: {:?}",
                value
            ))),
        }
    }
}

impl Serialize for SequenceKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        match self {
            SequenceKit::Standard(s) => s.serialize(serializer),
            SequenceKit::Custom(items) => items.serialize(serializer),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const YAML_FILE: &str = "../seqspec_templates/10x_rna_atac.yaml";
    const YAML_FILE_2: &str = "../seqspec_templates/smartseq2.yaml";

    #[test]
    fn test_parse() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");
        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        println!("{:?}", assay);
    }

    #[test]
    fn test_serialize() {
        fn se_de(yaml_str: &str) {
            let assay: Assay = serde_yaml::from_str(yaml_str).expect("Failed to parse YAML");
            let yaml_str_ = serde_yaml::to_string(&assay).unwrap();
            let assay_ = serde_yaml::from_str(&yaml_str_).expect("Failed to parse YAML");
            assert_eq!(assay, assay_);
        }

        se_de(&fs::read_to_string(YAML_FILE).expect("Failed to read file"));
    }

    #[test]
    fn test_index() {
        let yaml_str = fs::read_to_string(YAML_FILE_2).expect("Failed to read file");
        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
        for (read, index) in assay.get_segments_by_modality(Modality::RNA) {
            eprintln!(
                "{}: {:?}",
                read.read_id,
                index
                    .segments
                    .into_iter()
                    .map(|x| (x.region_type, x.range))
                    .collect::<Vec<_>>()
            );
        }
        for (read, index) in assay.get_segments_by_modality(Modality::ATAC) {
            println!(
                "{}: {:?}",
                read.read_id,
                index
                    .segments
                    .into_iter()
                    .map(|x| (x.region_type, x.range))
                    .collect::<Vec<_>>()
            );
        }
        for (read, index) in assay.get_segments_by_modality(Modality::Protein) {
            println!(
                "{}: {:?}",
                read.read_id,
                index
                    .segments
                    .into_iter()
                    .map(|x| (x.region_type, x.range))
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn test_update_read_primer_warning() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");
        let mut assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        // Capture log messages
        let (log_sender, log_receiver) = std::sync::mpsc::channel();
        let _guard = env_logger::builder()
            .format(move |_, record| {
                let _ = log_sender.send(record.args().to_string());
                Ok(())
            })
            .is_test(true)
            .try_init();

        // Test with non-sequencing primer (barcode) - should warn
        let result = assay.update_read::<PathBuf>(
            "test_read1",
            Some(Modality::RNA),
            Some("rna-cell_barcode"),  // cell barcode is not a sequencing primer
            Some(false),
            None,
            Some(16),
            Some(16),
            false,
        );
        assert!(result.is_ok());
        // Verify warning was logged for barcode primer
        let warning = log_receiver.try_recv().expect("Should have received warning");
        assert!(warning.contains("primer_id 'rna-cell_barcode' is not a sequencing primer"));
        assert!(warning.contains("type: Barcode"));

        // Test with sequencing primer - should not warn
        let result = assay.update_read::<PathBuf>(
            "test_read2", 
            Some(Modality::RNA),
            Some("rna-illumina_p5"),  // this is a sequencing primer
            Some(false),
            None,
            Some(29),
            Some(29),
            false,
        );
        assert!(result.is_ok());

        // Verify no warning was logged for sequencing primer
        assert!(log_receiver.try_recv().is_err(), "Should not have received warning for sequencing primer");

        // Verify reads were added with correct primers
        let read1 = assay.sequence_spec.get("test_read1").unwrap();
        let read2 = assay.sequence_spec.get("test_read2").unwrap();
        assert_eq!(read1.primer_id, "rna-cell_barcode");
        assert_eq!(read2.primer_id, "rna-illumina_p5");
    }

    #[test]
    // This test is to test variable length reads. Still in progress.
    fn test_update_read_with_fastq() {
        // Initialize logger
        let _ = env_logger::builder().is_test(true).try_init();

        // Use existing FASTQ files
        let fastq_path1 = PathBuf::from("../data/test_1.fastq");
        let fastq_path2 = PathBuf::from("../data/test_2.fastq");

        // Load test YAML
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");
        let mut assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        // Update read with first FASTQ file
        assay.update_read::<PathBuf>(
            "test_read1",
            Some(Modality::RNA),
            Some("rna-truseq_read1"),
            Some(false),
            Some(&[fastq_path1.clone()]),
            None,  // Let it determine length from FASTQ
            Some(16),
            false,
        ).expect("Failed to update read with test_1.fastq");

        // Update read with second FASTQ file
        assay.update_read::<PathBuf>(
            "test_read2",
            Some(Modality::RNA),
            Some("rna-truseq_read2"),
            Some(true),
            Some(&[fastq_path2.clone()]),
            None,
            None,
            false,
        ).expect("Failed to update read with test_2.fastq");

        // Verify reads
        let read1 = assay.sequence_spec.get("test_read1").expect("Read1 not found");
        let read2 = assay.sequence_spec.get("test_read2").expect("Read2 not found");
        
        // Print read information
        println!("Read1:");
        println!("  Length: {}", read1.min_len);
        //println!("  File: {:?}", read1.files.as_ref().unwrap()[0].path);
        println!("Read2:");
        println!("  Length: {}", read2.min_len);
        //println!("  File: {:?}", read2.files.as_ref().unwrap()[0].path);

        // Assert lengths (32 based on the FASTQ content you showed)
        assert_eq!(read1.min_len, 32, "Incorrect length for read1");
        assert_eq!(read2.min_len, 70, "Incorrect length for read2");
    }
    #[test]
    fn test_case_insensitive_yaml() {
        // Test YAML with different cases
        let yaml = r#"
region_id: "test_region"
region_type: "BARCODE"  # uppercase
name: "Test Region"
sequence_type: "fixed"
sequence: "ACGT"
min_len: 4
max_len: 4
regions: []
"#;
        let region: Region = serde_yaml::from_str(yaml).expect("Failed to parse YAML");
        assert_eq!(region.region_type, RegionType::Barcode);

        // Test mixed case
        let yaml = r#"
region_id: "test_region"
region_type: "TruSeq_Read1"  # mixed case
name: "Test Region"
sequence_type: "fixed"
sequence: "ACGT"
min_len: 4
max_len: 4
regions: []
"#;
        let region: Region = serde_yaml::from_str(yaml).expect("Failed to parse YAML");
        assert_eq!(region.region_type, RegionType::TruseqRead1);
    }

    #[test]
    fn test_get_segments_by_modality() {
        let yaml_str = fs::read_to_string(YAML_FILE_2).expect("Failed to read file");
        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        // Test RNA modality segments
        let rna_segments: Vec<_> = assay.get_segments_by_modality(Modality::RNA)
            .map(|(read, info)| {
                (
                    read.read_id.clone(),
                    info.segments.into_iter().map(|seg| {
                        (
                            seg.region_id.clone(),
                            seg.region_type,
                            seg.range.start..seg.range.end
                        )
                    }).collect::<Vec<_>>()
                )
            })
            .collect();

        // Print detailed segment information for debugging
        println!("\nRNA Segments:");
        for (read_id, segments) in &rna_segments {
            println!("Read {}", read_id);
            for (region_id, region_type, range) in segments {
                println!("  - Region: {}, Type: {:?}, Range: {:?}", 
                    region_id, region_type, range);
            }
        }

        // Test ATAC modality segments
        let atac_segments: Vec<_> = assay.get_segments_by_modality(Modality::ATAC)
            .map(|(read, info)| {
                (
                    read.read_id.clone(),
                    info.segments.into_iter().map(|seg| {
                        (
                            seg.region_id.clone(),
                            seg.region_type,
                            seg.range.start..seg.range.end
                        )
                    }).collect::<Vec<_>>()
                )
            })
            .collect();

        println!("\nATAC Segments:");
        for (read_id, segments) in &atac_segments {
            println!("Read {}", read_id);
            for (region_id, region_type, range) in segments {
                println!("  - Region: {}, Type: {:?}, Range: {:?}", 
                    region_id, region_type, range);
            }
        }

        // Add assertions to verify expected segments
        // Example assertions (adjust according to your expected values):
        assert!(!rna_segments.is_empty(), "Should have RNA segments");
        //assert!(!atac_segments.is_empty(), "Should have ATAC segments");

        // Verify specific RNA read segments
        if let Some(rna_read) = rna_segments.iter().find(|(id, _)| id == "RNA R1") {
            let segments = &rna_read.1;
            assert!(!segments.is_empty(), "RNA-R1 should have segments");
            // Add more specific assertions about the segments
        }

    }
}
