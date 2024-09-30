use std::{collections::HashMap, path::PathBuf, str::FromStr};

use pyo3::prelude::*;
use seqspec::{Assay, File, Modality, Read, Region, Strand, UrlType};
use anyhow::Result;
use termtree::Tree;
use cached_path::Cache;

/** A SeqSpec object.
    
    A SeqSpec object is used to annotate sequencing libraries produced by genomics assays.
    Genomic library structure depends on both the assay and sequencer (and kits) used to
    generate and bind the assay-specific construct to the sequencing adapters to generate
    a sequencing library. SeqSpec is specific to both a genomics assay and sequencer
    and provides a standardized format for describing the structure of sequencing
    libraries and the resulting sequencing reads. See https://github.com/pachterlab/seqspec for more details.

    Parameters
    ----------
    path
        The local path or url to the seqspec file.

    See Also
    --------
    align
*/
#[pyclass]
#[repr(transparent)]
pub struct SeqSpec(pub(crate) Assay);

#[pymethods]
impl SeqSpec {
    #[new]
    #[pyo3(signature = (path))] 
    pub fn new(path: &str) -> Result<Self> {
        let cache = Cache::new()?;
        let file = cache.cached_path(path)?;
        let assay = Assay::from_path(file)?;
        Ok(SeqSpec(assay))
    }

    /// Add a fastq file containing reads to the AnnData object.
    ///
    /// Parameters
    /// ----------
    /// read_id: str
    ///     The id of the read.
    /// modality: str
    ///    The modality of the read.
    /// primer_id: str
    ///   The id of the primer.
    /// is_reverse: bool
    ///   Whether the read is reverse.
    /// fastq: Path | list[Path]
    ///    The path to the fastq file containing the reads.
    /// min_len: int
    ///   The minimum length of the read. If not provided, the minimum length is inferred from the fastq file.
    /// max_len: int
    ///   The maximum length of the read. If not provided, the maximum length is inferred from the fastq file.
    #[pyo3(
        signature = (read_id, *, modality, primer_id, is_reverse, fastq, min_len=None, max_len=None),
        text_signature = "($self, read_id, *, modality, primer_id, is_reverse, fastq, min_len=None, max_len=None)",
    )]
    pub fn add_read(
        &mut self,
        read_id: &str,
        modality: &str,
        primer_id: &str,
        is_reverse: bool,
        fastq: Bound<'_, PyAny>,
        min_len: Option<usize>,
        max_len: Option<usize>,
    ) -> Result<()> {
        let fastq = if fastq.is_instance_of::<pyo3::types::PyList>() {
            fastq.extract::<Vec<PathBuf>>()?
        } else {
            vec![fastq.extract::<PathBuf>()?]
        };

        let assay = &mut self.0;

        let mut reads = assay.sequence_spec.take().unwrap_or(Vec::new());
        reads = reads.into_iter().filter(|r| r.read_id != read_id).collect();

        let mut read = Read::default();
        if is_reverse {
            read.strand = Strand::Neg;
        } else {
            read.strand = Strand::Pos;
        }
        read.read_id = read_id.to_string();
        read.modality = Modality::from_str(modality)?;
        read.primer_id = primer_id.to_string();
        read.files = Some(fastq.into_iter().map(|path| make_file_path(path)).collect::<Result<Vec<File>>>()?);

        if min_len.is_none() || max_len.is_none() {
            let len = precellar::io::get_read_length(&read, "./")?;
            read.min_len = min_len.unwrap_or(len) as u32;
            read.max_len = max_len.unwrap_or(len) as u32;
        } else {
            read.min_len = min_len.unwrap() as u32;
            read.max_len = max_len.unwrap() as u32;
        }

        reads.push(read);
        assay.sequence_spec = Some(reads);
        Ok(())
    }

    /*
    /// Identify the position of elements in a spec.
    ///
    /// Parameters
    /// ----------
    /// read_id: str
    ///     The id of the read.
    #[pyo3(
        signature = (modality=None),
        text_signature = "($self, modality=None)",
    )]
    pub fn index(&mut self, modality: &str) -> Result<()> {
    */

    #[pyo3(text_signature = "($self)")]
    pub fn to_yaml(&self) -> String {
        serde_yaml::to_string(&self.0).unwrap()
    }

    fn __repr__(&self) -> String {
        let assay = &self.0;
        let mut read_list = HashMap::new();
        if let Some(reads) = assay.sequence_spec.as_ref() {
            for read in reads {
                read_list.entry(read.primer_id.clone()).or_insert(Vec::new()).push(read);
            }
        }

        let tree = Tree::new("".to_string()).with_leaves(
            assay.library_spec.as_ref()
                .unwrap_or(&Vec::new()).iter()
                .map(|region| build_tree(region, &read_list))
            );
        format!("{}", tree)
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

fn build_tree(region: &Region, read_list: &HashMap<String, Vec<&Read>>) -> Tree<String> {
    let id = &region.region_id;
    let len = if region.min_len == region.max_len {
        region.min_len.to_string()
    } else {
        format!("{}-{}", region.min_len, region.max_len)
    };
    let label = if let Some(reads) = read_list.get(id) {
        let s = reads.iter().map(|r| format_read(r)).collect::<Vec<String>>().join(", ");
        format!("{}({}) [{}]", id, len, s)
    } else {
        format!("{}({})", id, len)
    };
    Tree::new(label)
        .with_leaves(region.regions.as_ref().unwrap_or(&Vec::new())
        .iter().map(|child| build_tree(child, read_list)))
}

fn format_read(read: &Read) -> String {
    let len = if read.min_len == read.max_len {
        read.min_len.to_string()
    } else {
        format!("{}-{}", read.min_len, read.max_len)
    };
    if read.is_reverse() {
        format!("↑{}({})", read.read_id, len)
    } else {
        format!("↓{}({})", read.read_id, len)
    }
}

fn make_file_path(path: PathBuf) -> Result<File> {
    let file = std::fs::File::open(&path)?;
    Ok(File {
        file_id: path.file_name().unwrap().to_str().unwrap().to_string(),
        filename: path.file_name().unwrap().to_str().unwrap().to_string(),
        filetype: "fastq".to_string(),
        filesize: file.metadata()?.len(),
        url: path.to_str().unwrap().to_string(),
        urltype: UrlType::Local,
        md5: "0".to_string(),
    })
}