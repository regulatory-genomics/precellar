use std::collections::HashMap;

use pyo3::prelude::*;
use seqspec::{Assay, Read, Region};
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

    /// Update read information in the SeqSpec object.
    /// 
    /// This method updates the read information in the SeqSpec object.
    /// If the read does not exist, it will be created.
    ///
    /// Parameters
    /// ----------
    /// read_id: str
    ///     The id of the read.
    /// modality: str | None
    ///    The modality of the read.
    /// primer_id: str | None
    ///   The id of the primer.
    /// is_reverse: bool | None
    ///   Whether the read is reverse.
    /// fastq: Path | list[Path] | None
    ///    The path to the fastq file containing the reads.
    /// min_len: int | None
    ///   The minimum length of the read. If not provided, the minimum length is inferred from the fastq file.
    /// max_len: int | None
    ///   The maximum length of the read. If not provided, the maximum length is inferred from the fastq file.
    #[pyo3(
        signature = (read_id, *, modality=None, primer_id=None, is_reverse=None, fastq=None, min_len=None, max_len=None),
        text_signature = "($self, read_id, *, modality=None, primer_id=None, is_reverse=None, fastq=None, min_len=None, max_len=None)",
    )]
    pub fn update_read(
        &mut self,
        read_id: &str,
        modality: Option<&str>,
        primer_id: Option<&str>,
        is_reverse: Option<bool>,
        fastq: Option<Bound<'_, PyAny>>,
        min_len: Option<usize>,
        max_len: Option<usize>,
    ) -> Result<()> {
        let fastqs = fastq.map(|f| if f.is_instance_of::<pyo3::types::PyList>() {
            f.extract::<Vec<String>>().unwrap()
        } else {
            vec![f.extract::<String>().unwrap()]
        });

        self.0.update_read(
            read_id, modality, primer_id, is_reverse,
            fastqs.as_ref().map(|x| x.as_slice()), min_len, max_len,
        )
    }

    /// Delete a read from the SeqSpec object.
    #[pyo3(signature = (read_id), text_signature = "($self, read_id)")]
    pub fn delete_read(&mut self, read_id: &str) -> Result<()> {
        let reads = self.0.sequence_spec.take();
        if let Some(r) = reads {
            self.0.sequence_spec = Some(
                r.into_iter().filter(|r| r.read_id != read_id).collect::<Vec<_>>()
            );
        }
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
    let orientation = if read.is_reverse() { "↑" } else { "↓" };
    let has_files = if read.files.as_ref().map(|x| !x.is_empty()).unwrap_or(false) { "✓" } else { "✗" };
    format!("{}{}({}){}", orientation, read.read_id, len, has_files)
}