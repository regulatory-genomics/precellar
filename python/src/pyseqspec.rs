use anyhow::Result;
use glob::glob;
use pyo3::prelude::*;
use seqspec::{Modality, Read, Region, SequenceKit, SequenceProtocol};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    str::FromStr,
};
use termtree::Tree;

pub(crate) fn extract_assays<'py>(assays: Bound<'py, PyAny>) -> PyResult<Vec<seqspec::Assay>> {
    fn extract<'py>(assay: Bound<'py, PyAny>) -> PyResult<seqspec::Assay> {
        let assay = match assay.extract::<PathBuf>() {
            Ok(p) => seqspec::Assay::from_path(&p)?,
            _ => assay.extract::<PyRef<Assay>>()?.0.clone(),
        };
        Ok(assay)
    }
    if assays.is_instance_of::<pyo3::types::PyList>() {
        let assays: Vec<Bound<'py, PyAny>> = assays.extract()?;
        assays.into_iter()
            .map(extract)
            .collect()
    } else {
        let assay = extract(assays)?;
        Ok(vec![assay])
    }
}

/** A Assay object.

    A Assay object is used to annotate sequencing libraries produced by genomics assays.
    Genomic library structure depends on both the assay and sequencer (and kits) used to
    generate and bind the assay-specific construct to the sequencing adapters to generate
    a sequencing library. Assay is specific to both a genomics assay and sequencer
    and provides a standardized format for describing the structure of sequencing
    libraries and the resulting sequencing reads. See https://github.com/pachterlab/seqspec for more details.

    Parameters
    ----------
    path: Path | str
        The local path or url to the seqspec file.

    See Also
    --------
    align
*/
#[pyclass]
#[repr(transparent)]
pub struct Assay(pub(crate) seqspec::Assay);

#[pymethods]
impl Assay {
    #[new]
    #[pyo3(signature = (path))]
    pub fn new(path: Bound<'_, PyAny>) -> Result<Self> {
        let assay = if let Ok(path) = path.extract::<&str>() {
            if url::Url::parse(path).is_ok() {
                seqspec::Assay::from_url(path)?
            } else {
                seqspec::Assay::from_path(path)?
            }
        } else {
            let path = path.extract::<PathBuf>()?;
            seqspec::Assay::from_path(path)?
        };
        Ok(Assay(assay))
    }

    /// The assay_id of the Assay object.
    #[getter]
    fn get_id(&self) -> &str {
        &self.0.assay_id
    }

    #[setter]
    fn set_id(&mut self, assay_id: &str) {
        self.0.assay_id = assay_id.to_string();
    }

    /// Name of the Assay object.
    #[getter]
    fn get_name(&self) -> &str {
        &self.0.name
    }

    #[setter]
    fn set_name(&mut self, name: &str) {
        self.0.name = name.to_string();
    }

    /// Description of the Assay object.
    #[getter]
    fn get_description(&self) -> &str {
        &self.0.description
    }

    #[setter]
    fn set_description(&mut self, description: &str) {
        self.0.description = description.to_string();
    }

    /// The doi of the Assay object.
    #[getter]
    fn get_doi(&self) -> &str {
        &self.0.doi
    }

    #[setter]
    fn set_doi(&mut self, doi: &str) {
        self.0.doi = doi.to_string();
    }

    /// The sequence protocol of the Assay object.
    #[getter]
    fn get_sequence_protocol(&self) -> String {
        match &self.0.sequence_protocol {
            SequenceProtocol::Standard(x) => x.clone(),
            SequenceProtocol::Custom(x) => format!("custom: {:?}", x),
        }
    }

    #[setter]
    fn set_sequence_protocol(&mut self, sequence_protocol: &str) {
        self.0.sequence_protocol = SequenceProtocol::Standard(sequence_protocol.to_string());
    }

    /// The sequence kit of the Assay object.
    #[getter]
    fn get_sequence_kit(&self) -> String {
        match &self.0.sequence_kit {
            SequenceKit::Standard(x) => x.clone(),
            SequenceKit::Custom(x) => format!("custom: {:?}", x),
        }
    }

    #[setter]
    fn set_sequence_kit(&mut self, sequence_kit: &str) {
        self.0.sequence_kit = SequenceKit::Standard(sequence_kit.to_string());
    }

    /// Filename of the seqspec file.
    ///
    /// Returns
    /// -------
    /// Path | None
    #[getter]
    pub fn filename(&self) -> Option<PathBuf> {
        self.0.file.clone()
    }

    /// Return all modality names in the Assay object.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///  A list of all modality names in the Assay object.
    pub fn modalities(&self) -> Vec<String> {
        self.0.modalities.iter().map(|x| x.to_string()).collect()
    }

    /// Add default Illumina reads to the Assay object.
    ///
    /// This method adds default Illumina reads to the Assay object.
    ///
    /// Parameters
    /// ----------
    /// modality: str
    ///    The modality of the read.
    /// length: int
    ///   The length of the read.
    /// forward_strand_workflow: bool
    ///   Whether the reads are sequenced in the forward strand workflow.
    #[pyo3(
        signature = (modality, *, length=150, forward_strand_workflow=false),
        text_signature = "($self, modality, *, length=150, forward_strand_workflow=False)",
    )]
    pub fn add_illumina_reads(
        &mut self,
        modality: &str,
        length: usize,
        forward_strand_workflow: bool,
    ) -> Result<()> {
        let modality = Modality::from_str(modality)?;
        self.0
            .add_illumina_reads(modality, length, forward_strand_workflow)
    }

    /// Update read information in the Assay object.
    ///
    /// This method updates the read information in the Assay object.
    /// If the read does not exist, it will be created.
    ///
    /// Parameters
    /// ----------
    /// read_id: str
    ///     The id of the read.
    /// modality: str | None
    ///     The modality of the read.
    /// primer_id: str | None
    ///     The id of the primer.
    /// is_reverse: bool | None
    ///     Whether the read is reverse.
    /// fastq: Path | list[Path] | None
    ///     The path to the fastq file containing the reads. It can be:
    ///     (1) a single path, or (2) a list of paths, or (3) a glob pattern, e.g., `*.fastq`.
    /// min_len: int | None
    ///     The minimum length of the read. If not provided, the minimum length is inferred from the fastq file.
    /// max_len: int | None
    ///     The maximum length of the read. If not provided, the maximum length is inferred from the fastq file.
    /// compute_md5: bool
    ///     Whether to compute the md5 hash of the fastq file.
    /// infer_read_length: bool
    ///     Whether to infer the read length from the fastq file.
    /// infer_read_length_sample: int | None
    ///     The number of records to read from the fastq file to infer the read length.
    #[pyo3(
        signature = (read_id, *, modality=None, primer_id=None, is_reverse=None, fastq=None, min_len=None, max_len=None, compute_md5=false, infer_read_length=true, infer_read_length_sample=10000),
        text_signature = "($self, read_id, *, modality=None, primer_id=None, is_reverse=None, fastq=None, min_len=None, max_len=None, compute_md5=False, infer_read_length=True, infer_read_length_sample=10000)",
    )]
    fn update_read(
        &mut self,
        read_id: &str,
        modality: Option<&str>,
        primer_id: Option<&str>,
        is_reverse: Option<bool>,
        fastq: Option<Bound<'_, PyAny>>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        compute_md5: bool,
        infer_read_length: bool,
        infer_read_length_sample: usize,
    ) -> Result<()> {
        let fastqs = fastq.map(|f| {
            if f.is_instance_of::<pyo3::types::PyList>() {
                f.extract::<Vec<PathBuf>>().unwrap()
            } else {
                if let Ok(path) = f.extract::<&str>() {
                    glob(path)
                        .expect("Failed to read glob pattern")
                        .map(Result::unwrap)
                        .collect::<Vec<_>>()
                } else {
                    let path = f.extract::<PathBuf>().unwrap();
                    vec![path]
                }
            }
        });
        let modality = modality.map(|x| Modality::from_str(x)).transpose()?;

        self.0.update_read(
            read_id,
            modality,
            primer_id,
            is_reverse,
            fastqs.as_ref().map(|x| x.as_slice()),
            min_len,
            max_len,
            compute_md5,
            infer_read_length,
            infer_read_length_sample,
        )
    }

    /// Delete a read from the Assay object.
    #[pyo3(signature = (read_id), text_signature = "($self, read_id)")]
    fn delete_read(&mut self, read_id: &str) {
        self.0.delete_read(read_id);
    }

    /// Return the whitelist of cell barcodes.
    #[pyo3(
        signature = (modality),
        text_signature = "($self, modality)",
    )]
    fn whitelist(&self, modality: &str) -> Option<Vec<String>> {
        let modality = Modality::from_str(modality).unwrap();
        let barcodes = self
            .0
            .library_spec
            .cat_barcodes(&modality)?
            .into_iter()
            .map(|bc| std::str::from_utf8(&bc).unwrap().to_string())
            .collect();
        Some(barcodes)
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

    /// Convert the Assay object to a yaml string.
    ///
    /// This method converts the Assay object to a yaml string. If you want to save the
    /// Assay object to a yaml file, use the `.save()` method as it will convert all paths
    /// to relative paths so that the files can be moved without breaking the Assay object.
    ///
    /// Returns
    /// -------
    /// str
    ///   The yaml string representation of the Assay object.
    ///
    /// See Also
    /// --------
    /// save
    #[pyo3(text_signature = "($self)")]
    fn to_yaml(&mut self) -> Result<String> {
        Ok(serde_yaml::to_string(&self.0)?)
    }

    /// Save the Assay object to a yaml file.
    ///
    /// This method saves the Assay object to a yaml file. All paths in the Assay object
    /// will be converted to relative paths with respect to the location of the output file.
    ///
    /// Parameters
    /// ----------
    /// out: Path
    ///    The path to save the yaml file.
    #[pyo3(text_signature = "($self, out)")]
    fn save(&self, out: PathBuf) -> Result<()> {
        let mut assay = self.0.clone();
        let parent = out.parent().unwrap();
        if parent == Path::new("") {
            assay.unnormalize_all_paths("./");
        } else {
            assay.unnormalize_all_paths(parent);
        }
        std::fs::write(&out, serde_yaml::to_string(&assay)?)?;
        Ok(())
    }

    fn __repr__(&self) -> String {
        let assay = &self.0;
        let mut read_list = HashMap::new();
        for read in assay.sequence_spec.values() {
            read_list
                .entry(read.primer_id.clone())
                .or_insert(Vec::new())
                .push(read);
        }

        let root_label = self.0.assay_id.clone();
        let tree = Tree::new(root_label).with_leaves(
            assay
                .library_spec
                .modalities()
                .map(|region| build_tree(&region.read().unwrap(), &read_list)),
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
        let s = reads
            .iter()
            .map(|r| format_read(r))
            .collect::<Vec<String>>()
            .join(", ");
        format!("{}({}) [{}]", id, len, s)
    } else {
        format!("{}({})", id, len)
    };
    Tree::new(label).with_leaves(
        region
            .subregions
            .iter()
            .map(|child| build_tree(&child.read().unwrap(), read_list)),
    )
}

fn format_read(read: &Read) -> String {
    let len = if read.min_len == read.max_len {
        read.min_len.to_string()
    } else {
        format!("{}-{}", read.min_len, read.max_len)
    };
    let orientation = if read.is_reverse() { "↑" } else { "↓" };
    let has_files = if read.files.as_ref().map(|x| !x.is_empty()).unwrap_or(false) {
        "✓"
    } else {
        "✗"
    };
    format!("{}{}({}){}", orientation, read.read_id, len, has_files)
}
