mod read;
mod region;
pub mod utils;

use log::warn;
use read::{ReadSpan, RegionIndex, ValidateResult};
pub use read::{Read, File, UrlType, Strand};
use read::ReadValidator;
use region::LibSpec;
pub use region::{Region, RegionType, SequenceType, Onlist};
use noodles::fastq;

use serde::{Deserialize, Deserializer, Serialize};
use serde_yaml::{self, Value};
use utils::open_file_for_write;
use std::{fs, path::PathBuf, str::FromStr, sync::Arc};
use anyhow::{bail, anyhow, Result};
use std::path::Path;

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
        let yaml_str = fs::read_to_string(&path)?;
        let mut assay: Assay = serde_yaml::from_str(&yaml_str)?;
        assay.file = Some(path.as_ref().to_path_buf());
        assay.normalize_all_paths();
        Ok(assay)
    }

    pub fn from_url(url: &str) -> Result<Self> {
        let yaml_str = reqwest::blocking::get(url)?.text()?;
        Ok(serde_yaml::from_str(&yaml_str)?)
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
            self.library_spec.regions_mut().for_each(|region| {
                if let Some(onlist) = &mut Arc::<Region>::get_mut(region).unwrap().onlist {
                    if let Err(e) = onlist.normalize_path(base_dir) {
                        warn!("{}", e);
                    }
                }
            });
        }
    }

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
        /*
        self.library_spec.regions_mut().for_each(|region| {
            if let Some(onlist) = &mut Arc::<Region>::get_mut(region).unwrap().onlist {
                if let Err(e) = onlist.unnormalize_path(base_dir.as_ref()) {
                    warn!("Failed to unnormalize path: {}", e);
                }
            }
        });
        */
    }

    /// Add default Illumina reads to the sequence spec.
    pub fn add_illumina_reads(&mut self, modality: Modality, read_len: usize, forward_strand_workflow: bool) -> Result<()> {
        fn advance_until(iterator: &mut std::slice::Iter<'_, Arc<Region>>, f: fn(&Region) -> bool) -> Option<(Arc<Region>, Vec<Arc<Region>>)> {
            let mut regions = Vec::new();
            while let Some(next_region) = iterator.next() {
                if f(next_region) {
                    return Some((next_region.clone(), regions))
                } else {
                    regions.push(next_region.clone());
                }
            }
            None
        }

        fn get_length(regions: &[Arc<Region>], reverse: bool) -> usize {
            if reverse {
                regions.iter()
                    .skip_while(|region| region.sequence_type == SequenceType::Fixed)
                    .map(|region| region.len().unwrap() as usize).sum()
            } else {
                regions.iter().rev()
                    .skip_while(|region| region.sequence_type == SequenceType::Fixed)
                    .map(|region| region.len().unwrap() as usize).sum()
            }
        }

        fn is_read1(region: &Region) -> bool {
            region.region_type == RegionType::NexteraRead1 || region.region_type == RegionType::TruseqRead1
        }

        fn is_read2(region: &Region) -> bool {
            region.region_type == RegionType::NexteraRead2 || region.region_type == RegionType::TruseqRead2
        }

        fn is_p5(region: &Region) -> bool {
            region.region_type == RegionType::IlluminaP5
        }

        fn is_p7(region: &Region) -> bool {
            region.region_type == RegionType::IlluminaP7
        }

        self.delete_all_reads(modality);
        let regions = self.library_spec.get_modality(&modality).ok_or_else(|| anyhow!("Modality not found: {:?}", modality))?.clone();
        let mut regions = regions.subregions.iter();
        while let Some(current_region) = regions.next() {
            if is_p5(&current_region) {
                if let Some((next_region, acc)) = advance_until(&mut regions, is_read1) {
                    self.update_read::<PathBuf>(
                        &format!("{}-R1", modality.to_string()),
                        Some(modality),
                        Some(&next_region.region_id),
                        Some(false),
                        None, Some(read_len), Some(read_len)
                    )?;
                    if forward_strand_workflow {
                        let acc_len = get_length(acc.as_slice(), false);
                        if acc_len > 0 {
                            self.update_read::<PathBuf>(
                                &format!("{}-I1", modality.to_string()),
                                Some(modality),
                                Some(&current_region.region_id),
                                Some(false),
                                None, Some(acc_len), Some(acc_len)
                            )?;
                        }
                    } else {
                        let acc_len = get_length(acc.as_slice(), true);
                        if acc_len > 0 {
                            self.update_read::<PathBuf>(
                                &format!("{}-I1", modality.to_string()),
                                Some(modality),
                                Some(&next_region.region_id),
                                Some(true),
                                None, Some(acc_len), Some(acc_len)
                            )?;
                        }
                    }
                }
            } else if is_read2(&current_region) {
                if let Some((_, acc)) = advance_until(&mut regions, is_p7) {
                    let acc_len = get_length(acc.as_slice(), false);
                    self.update_read::<PathBuf>(
                        &format!("{}-R2", modality.to_string()),
                        Some(modality),
                        Some(&current_region.region_id),
                        Some(true),
                            None, Some(read_len), Some(read_len)
                    )?;
                    if acc_len > 0 {
                        self.update_read::<PathBuf>(
                            &format!("{}-I2", modality.to_string()),
                            Some(modality),
                            Some(&current_region.region_id),
                            Some(false),
                            None, Some(acc_len), Some(acc_len)
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
    ) -> Result<()> {
        let mut read_buffer = if let Some(r) = self.sequence_spec.get(read_id) {
            r.clone()
        } else {
            assert!(modality.is_some(), "modality must be provided for a new read");
            assert!(primer_id.is_some(), "primer_id must be provided for a new read");
            assert!(is_reverse.is_some(), "is_reverse must be provided for a new read");
            Read { read_id: read_id.to_string(), ..Default::default() }
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
            read_buffer.files = Some(fastq.iter().map(|path| File::from_fastq(path)).collect::<Result<Vec<File>>>()?);
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
        self.sequence_spec.insert(read_id.to_string(), read_buffer);

        Ok(())
    }

    pub fn delete_read(&mut self, read_id: &str) -> Option<Read> {
        self.sequence_spec.shift_remove(read_id)
    }

    pub fn delete_all_reads(&mut self, modality: Modality) {
        let ids: Vec<_> = self.iter_reads(modality).map(|read| {
            read.read_id.to_string()
        }).collect();
        ids.into_iter().for_each(|id| {
            self.delete_read(&id);
        });
    }

    /// Get the index of atomic regions of each read in the sequence spec.
    pub fn get_index_by_modality(&self, modality: Modality) -> impl Iterator<Item = (&Read, RegionIndex)> {
        self.sequence_spec.values().filter_map(move |read| if read.modality == modality {
            let index = self.get_index(&read.read_id)
                .expect(&format!("Cannot find index for Read: {}", read.read_id));
            Some((read, index))
        } else {
            None
        })
    }

    /// Get the index of atomic regions of a read in the sequence spec.
    pub fn get_index(&self, read_id: &str) -> Option<RegionIndex> {
        let read = self.sequence_spec.get(read_id)?;
        let region = self.library_spec.get_parent(&read.primer_id)?;
        read.get_index(region)
    }

    pub fn iter_reads(&self, modality: Modality) -> impl Iterator<Item = &Read> {
        self.sequence_spec.values().filter(move |read| read.modality == modality)
    }

    /// Demultiplex reads into separate files based on the primer ID.
    pub fn validate<P: AsRef<Path>>(&self, read: &Read, dir: P) -> Result<()> {
        let region = self.library_spec.get_parent(&read.primer_id)
            .ok_or_else(|| anyhow!("Primer not found: {}", read.primer_id))?;
        if let Some(index) = read.get_index(region) {
            let output_valid = dir.as_ref().join(format!("{}.fq.zst", read.read_id));
            let output_valid = open_file_for_write(output_valid, None, None, 8)?;
            let mut output_valid = fastq::io::Writer::new(output_valid);
            let output_other = dir.as_ref().join(format!("Invalid_{}.fq.zst", read.read_id));
            let output_other = open_file_for_write(output_other, None, None, 8)?;
            let mut output_other = fastq::io::Writer::new(output_other);
            if let Some(mut reader) = read.open() {
                let mut validators: Vec<_> = index.index.iter().map(|(region_id, _, range)| {
                    let region = self.library_spec.get(region_id).unwrap();
                    ReadValidator::new(region)
                        .with_range(range.start as usize ..range.end as usize)
                        .with_strand(read.strand)
                }).collect();

                reader.records().try_for_each(|record| {
                    let record = record?;
                    let valid = validators.iter_mut().all(|validator| {
                        match validator.validate(record.sequence()) {
                            ValidateResult::OnlistFail | ValidateResult::Valid => true,
                            _ => false,
                        }
                    });
                    if valid {
                        output_valid.write_record(&record)?;
                    } else {
                        output_other.write_record(&record)?;
                    }
                    anyhow::Ok(())
                })?;
            }
        }
        Ok(())
    }

    fn verify(&self, read: &Read) -> Result<()> {
        let region = self.library_spec.get_parent(&read.primer_id)
            .ok_or_else(|| anyhow!("Primer not found: {}", read.primer_id))?;
        // Check if the primer exists
        if let Some(index) = read.get_index(region) {
            match index.readlen_info {
                ReadSpan::Covered | ReadSpan::Span(_) => {},
                ReadSpan::NotEnough => {
                    warn!("'{}' does not cover the region", read.read_id);
                },
                ReadSpan::MayReadThrough(id) => {
                    warn!("'{}' may read through and contain sequences from: '{}'", read.read_id, id);
                },
                ReadSpan::ReadThrough(id) => {
                    warn!("'{}' length exceeds maximum length of the variable-length region (insertion), \
                    truncating the reads to the maximum length of the region. \
                    Read reads through and contains sequences from: '{}'.
                    If this is not the desired behavior, please adjust the region lengths.", read.read_id, id);
                }
            }

            if let Some(mut reader) = read.open() {
                let mut validators: Vec<_> = index.index.iter().map(|(region_id, _, range)| {
                    let region = self.library_spec.get(region_id).unwrap();
                    ReadValidator::new(region)
                        .with_range(range.start as usize ..range.end as usize)
                        .with_strand(read.strand)
                }).collect();

                reader.records().take(500).try_for_each(|record| {
                    let record = record?;
                    for validator in &mut validators {
                        let result = validator.validate(record.sequence());
                        match result {
                            ValidateResult::TooLong(_) | ValidateResult::TooShort(_) => {
                                bail!("{}", result);
                            }
                            _ => {},
                        }
                    }
                    anyhow::Ok(())
                })?;

                for validator in validators {
                    let percent_matched = validator.frac_matched() * 100.0;
                    if percent_matched < 50.0 {
                        warn!("Read '{}' has low percentage of matched bases for region '{}'. \
                        Percentage of matched bases: {:.2}%", read.read_id, validator.id(), percent_matched);
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

impl<'de> Deserialize<'de> for Modality {
    fn deserialize<D>(deserializer: D) -> Result<Modality, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Modality::from_str(&s).map_err(serde::de::Error::custom),
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
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

impl ToString for Modality {
    fn to_string(&self) -> String {
        match self {
            Modality::DNA => "dna".to_string(),
            Modality::RNA => "rna".to_string(),
            Modality::Tag => "tag".to_string(),
            Modality::Protein => "protein".to_string(),
            Modality::ATAC => "atac".to_string(),
            Modality::Crispr => "crispr".to_string(),
        }
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
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(LibraryProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for LibraryProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
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
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(LibraryKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for LibraryKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
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

impl <'de> Deserialize<'de> for SequenceProtocol {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = Value::deserialize(deserializer)?;
        match value {
            Value::String(s) => Ok(SequenceProtocol::Standard(s)),
            Value::Sequence(seq) => {
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<ProtocolItem>>();
                Ok(SequenceProtocol::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for SequenceProtocol {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
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
                let items = seq.into_iter().map(|item| serde_yaml::from_value(item).unwrap())
                    .collect::<Vec<KitItem>>();
                Ok(SequenceKit::Custom(items))
            }
            _ => Err(serde::de::Error::custom(format!("invalid value: {:?}", value))),
        }
    }
}

impl Serialize for SequenceKit {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer
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

    #[test]
    fn test_parse() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");
        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");

        println!("{:?}", assay);
    }

    #[test]
    fn test_serialize() {
        fn se_de(yaml_str: &str) {
            let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
            let yaml_str_ = serde_yaml::to_string(&assay).unwrap();
            let assay_ = serde_yaml::from_str(&yaml_str_).expect("Failed to parse YAML");
            assert_eq!(assay, assay_);
        }

        se_de(&fs::read_to_string(YAML_FILE).expect("Failed to read file"));
    }

    #[test]
    fn test_index() {
        let yaml_str = fs::read_to_string(YAML_FILE).expect("Failed to read file");

        let assay: Assay = serde_yaml::from_str(&yaml_str).expect("Failed to parse YAML");
        for (read, index) in assay.get_index_by_modality(Modality::RNA) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
        for (read, index) in assay.get_index_by_modality(Modality::ATAC) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
        for (read, index) in assay.get_index_by_modality(Modality::Protein) {
            println!("{}: {:?}", read.read_id, index.index.into_iter().map(|x| (x.1, x.2)).collect::<Vec<_>>());
        }
    }
}