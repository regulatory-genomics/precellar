use std::{ops::{Deref, DerefMut}, path::PathBuf};
use anyhow::Result;
use noodles::sam::Header;
use precellar::{align::{Aligner, AnnotatedFastq}, transcriptome::{Transcript, TxAligner}};
use pyo3::prelude::*;
use seqspec::ChemistryStrandedness;
use star_aligner::{StarAligner, StarOpts};
use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex};
use precellar::align::{Minimap2Aligner, Minimap2Opts};

pub enum AlignerRef<'py> {
    STAR(PyRefMut<'py, STAR>),
    BWA(PyRefMut<'py, BWAMEM2>),
    Minimap2(PyRefMut<'py, MINIMAP2>),
}

impl AlignerRef<'_> {
    pub fn header(&self) -> Header {
        match self {
            AlignerRef::STAR(aligner) => aligner.header(),
            AlignerRef::BWA(aligner) => aligner.header(),
            AlignerRef::Minimap2(aligner) => aligner.header(),
        }
    }

    pub fn transcript_annotator(&self, strandness: Option<ChemistryStrandedness>) -> Option<TxAligner> {
        match self {
            AlignerRef::STAR(aligner) => {
                let transcriptome: Vec<_> = aligner.get_transcriptome().unwrap().iter()
                    .map(|t| Transcript::try_from(t.clone()).unwrap()).collect();
                Some(TxAligner::new(transcriptome, self.header(), strandness))
            }
            AlignerRef::BWA(_) => None,
            AlignerRef::Minimap2(_) => None,
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
        } else if let Ok(aligner) = value.extract::<PyRefMut<'_, MINIMAP2>>() {
            Ok(AlignerRef::Minimap2(aligner))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected a STAR, BWA, or MINIMAP2 aligner",
            ))
        }
    }
}

impl Aligner for AlignerRef<'_> {
    fn header(&self) -> noodles::sam::Header {
        self.header()
    }

    fn align_reads(
            &mut self,
            num_threads: u16,
            records: Vec<AnnotatedFastq>,
        ) -> Vec<(Option<precellar::align::MultiMapR>, Option<precellar::align::MultiMapR>)> {
        match self {
            AlignerRef::STAR(aligner) => aligner.align_reads(num_threads, records),
            AlignerRef::BWA(aligner) => Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records),
            AlignerRef::Minimap2(aligner) => Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records),
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
        let opts = StarOpts::new(index_path).with_sam_attributes("NH HI AS nM");
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
        );
        Ok(BWAMEM2(aligner))
    }

    /// The maximum number of occurrences of a seed in the reference.
    /// Skip a seed if its occurrence is larger than this value. The default is 500.
    #[getter]
    pub fn get_max_occurrence(&self) -> u16 {
        self.0.opts.max_occurrence()
    }

    #[setter]
    pub fn set_max_occurrence(&mut self, max_occurence: u16) {
        self.0.opts.set_max_occurrence(max_occurence);
    }

    /// The minimum seed length of the aligner. The shorter the seed more
    /// sensitive the search will be. The default value is 19.
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

/** The Minimap2 aligner.

    Minimap2 is a versatile aligner for long reads (Oxford Nanopore, PacBio),
    splice alignment, assembly-to-assembly alignment, and more.

    Parameters
    ----------
    index_path : str
        The path to the Minimap2 index file (.mmi).
    preset : str | None
        The minimap2 preset to use. Available presets:
        - 'map-ont': Oxford Nanopore genomic reads (default)
        - 'map-pb': PacBio CLR genomic reads
        - 'map-hifi': PacBio HiFi/CCS genomic reads
        - 'splice': Long-read spliced alignment (RNA-seq)
        - 'splice:hq': High-quality long-read spliced alignment
        - 'asm5': Assembly-to-assembly alignment (divergence ~5%)
        - 'asm10': Assembly-to-assembly alignment (divergence ~10%)
        - 'asm20': Assembly-to-assembly alignment (divergence ~20%)
        - 'short': Short single-end reads
        - 'sr': Short paired-end reads
        If None, defaults to 'map-ont'.
*/
#[pyclass]
#[repr(transparent)]
pub struct MINIMAP2(Minimap2Aligner);

impl Deref for MINIMAP2 {
    type Target = Minimap2Aligner;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for MINIMAP2 {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[pymethods]
impl MINIMAP2 {
    #[new]
    #[pyo3(
        signature = (index_path, *, preset=None),
        text_signature = "(index_path, *, preset=None)",
    )]
    pub fn new(index_path: PathBuf, preset: Option<&str>) -> Result<Self> {
        let mut opts = Minimap2Opts::new(index_path);

        // Parse preset string to minimap2::Preset enum
        if let Some(preset_str) = preset {
            let preset_enum = match preset_str.to_lowercase().as_str() {
                "map-ont" => minimap2::Preset::MapOnt,
                "map-pb" => minimap2::Preset::MapPb,
                "map-hifi" => minimap2::Preset::MapHifi,
                "splice" => minimap2::Preset::Splice,
                "splice:hq" => minimap2::Preset::SpliceHq,
                "asm5" => minimap2::Preset::Asm5,
                "asm10" => minimap2::Preset::Asm10,
                "asm20" => minimap2::Preset::Asm20,
                "short" => minimap2::Preset::Short,
                "sr" => minimap2::Preset::Sr,
                _ => return Err(anyhow::anyhow!(
                    "Invalid preset '{}'. Valid presets: map-ont, map-pb, map-hifi, splice, splice:hq, asm5, asm10, asm20, short, sr",
                    preset_str
                )),
            };
            opts = opts.with_preset(preset_enum);
        }

        Ok(MINIMAP2(Minimap2Aligner::new(opts)?))
    }

    /// Get the currently configured preset name.
    ///
    /// Returns
    /// -------
    /// str | None
    ///     The preset name, or None if using default (map-ont).
    #[getter]
    pub fn get_preset(&self) -> Option<String> {
        self.0.get_opts().preset().map(|p| format!("{:?}", p).to_lowercase())
    }
}

#[pymodule]
pub(crate) fn register_aligners(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "aligners")?;

    m.add_class::<STAR>()?;
    m.add_class::<BWAMEM2>()?;
    m.add_class::<MINIMAP2>()?;

    parent_module.add_submodule(&m)
}