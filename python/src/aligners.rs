use std::{ops::{Deref, DerefMut}, path::PathBuf};
use anyhow::Result;
use noodles::sam::Header;
use precellar::{align::Aligner, transcript::{AlignmentAnnotator, Transcript}};
use pyo3::prelude::*;
use star_aligner::{StarAligner, StarOpts};
use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex, PairedEndStats};

pub enum AlignerRef<'py> {
    STAR(PyRefMut<'py, STAR>),
    BWA(PyRefMut<'py, BWAMEM2>),
}

impl AlignerRef<'_> {
    pub fn header(&self) -> Header {
        match self {
            AlignerRef::STAR(aligner) => aligner.header(),
            AlignerRef::BWA(aligner) => aligner.header(),
        }
    }

    pub fn transcript_annotator(&self) -> Option<AlignmentAnnotator> {
        match self {
            AlignerRef::STAR(aligner) => {
                let transcriptome: Vec<_> = aligner.get_transcriptome().unwrap().iter()
                    .map(|t| Transcript::try_from(t.clone()).unwrap()).collect();
                Some(AlignmentAnnotator::new(transcriptome))
            }
            AlignerRef::BWA(_) => None,
        }
    }
}

impl<'py> TryFrom<Bound<'py, PyAny>> for AlignerRef<'py> {
    type Error = PyErr;

    fn try_from(value: Bound<'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(aligner) = value.extract::<PyRefMut<'_, STAR>>() {
            Ok(AlignerRef::STAR(aligner))
        } else if let Ok(aligner) = value.extract::<PyRefMut<'_, BWAMEM2>>() {
            Ok(AlignerRef::BWA(aligner))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected a STAR or BWA aligner",
            ))
        }
    }
}

/** The STAR aligner.

    STAR aligner is a fast and accurate RNA-seq aligner. It is used to align RNA-seq reads to a reference genome.

    Parameters
    ----------
    index_path : str
        The path to the STAR index directory.
*/
#[pyclass]
#[repr(transparent)]
pub struct STAR(StarAligner);

impl Deref for STAR {
    type Target = StarAligner;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for STAR {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[pymethods]
impl STAR {
    #[new]
    #[pyo3(
        signature = (index_path),
        text_signature = "(index_path)",
    )]
    pub fn new(index_path: PathBuf) -> Result<Self> {
        let opts = StarOpts::new(index_path);
        Ok(STAR(StarAligner::new(opts)?))
    }
}

/** The BWA-MEM2 aligner.

    BWA-MEM2 is a fast and accurate genome aligner. It is used to align reads to a reference genome.

    Parameters
    ----------
    index_path : str
        The path to the BWA-MEM2 index directory.
*/
#[pyclass]
#[repr(transparent)]
pub struct BWAMEM2(BurrowsWheelerAligner);

impl Deref for BWAMEM2 {
    type Target = BurrowsWheelerAligner;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for BWAMEM2 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[pymethods]
impl BWAMEM2 {
    #[new]
    #[pyo3(
        signature = (index_path),
        text_signature = "(index_path)",
    )]
    pub fn new(index_path: PathBuf) -> Result<Self> {
        let aligner = BurrowsWheelerAligner::new(
            FMIndex::read(index_path).unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default(),
        );
        Ok(BWAMEM2(aligner))
    }

    /// The minimum seed length of the aligner.
    /// 
    /// Returns
    /// -------
    /// int
    ///    The minimum seed length.
    #[getter]
    pub fn get_min_seed_length(&self) -> u16 {
        self.0.opts.min_seed_len()
    }

    #[setter]
    pub fn set_min_seed_length(&mut self, min_seed_length: u16) {
        self.0.opts.set_min_seed_len(min_seed_length);
    }

    /// Whether to output log messages.
    pub fn logging(&mut self, enable: bool) {
        if enable {
            self.0.opts.enable_log();
        } else {
            self.0.opts.disable_log();
        }
    }
}

#[pymodule]
pub(crate) fn register_aligners(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "aligners")?;

    m.add_class::<STAR>()?;
    m.add_class::<BWAMEM2>()?;

    parent_module.add_submodule(&m)
}