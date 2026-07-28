use crate::align::parse_minimap2_preset;
use anyhow::{bail, Context, Result};
use bwa_mem2::{AlignerOpts, BurrowsWheelerAligner, FMIndex};
use minibwa::{Index as MiniBwaIndex, MiniBwaSR, Options as MiniBwaOptions};
use noodles_sam::Header;
use precellar::align::{Minimap2Aligner, Minimap2Opts};
use precellar::{
    align::{Aligner, AlignmentInput},
    transcriptome::{Transcript, TxAligner},
};
use pyo3::prelude::*;
use seqspec::ChemistryStrandedness;
use star_aligner::{StarAligner, StarOpts};
use std::{
    ops::{Deref, DerefMut},
    path::PathBuf,
};

pub enum AlignerRef<'py> {
    STAR(PyRefMut<'py, STAR>),
    BWA(PyRefMut<'py, BWAMEM2>),
    MiniBwa(PyRefMut<'py, MINIBWA>),
    Minimap2(PyRefMut<'py, MINIMAP2>),
}

impl AlignerRef<'_> {
    pub fn header(&self) -> Header {
        match self {
            AlignerRef::STAR(aligner) => aligner.header(),
            AlignerRef::BWA(aligner) => aligner.header(),
            AlignerRef::MiniBwa(aligner) => aligner.header(),
            AlignerRef::Minimap2(aligner) => aligner.header(),
        }
    }

    pub fn transcript_annotator(
        &self,
        strandness: Option<ChemistryStrandedness>,
    ) -> Option<TxAligner> {
        match self {
            AlignerRef::STAR(aligner) => {
                let transcriptome: Vec<_> = aligner
                    .get_transcriptome()
                    .unwrap()
                    .iter()
                    .map(|t| Transcript::try_from(t.clone()).unwrap())
                    .collect();
                Some(TxAligner::new(transcriptome, self.header(), strandness))
            }
            AlignerRef::BWA(_) => None,
            AlignerRef::MiniBwa(_) => None,
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
        } else if let Ok(aligner) = value.extract::<PyRefMut<'_, MINIBWA>>() {
            Ok(AlignerRef::MiniBwa(aligner))
        } else if let Ok(aligner) = value.extract::<PyRefMut<'_, MINIMAP2>>() {
            Ok(AlignerRef::Minimap2(aligner))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected a Star, BwaMem2, MiniBwa, or Minimap2 aligner",
            ))
        }
    }
}

impl Aligner for AlignerRef<'_> {
    fn header(&self) -> noodles_sam::Header {
        self.header()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(
        Option<precellar::align::MultiMapR>,
        Option<precellar::align::MultiMapR>,
    )> {
        match self {
            AlignerRef::STAR(aligner) => {
                Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records)
            }
            AlignerRef::BWA(aligner) => {
                Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records)
            }
            AlignerRef::MiniBwa(aligner) => {
                Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records)
            }
            AlignerRef::Minimap2(aligner) => {
                Aligner::align_reads(aligner.deref_mut().deref_mut(), num_threads, records)
            }
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
#[pyclass(name = "Star")]
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
        The path prefix for the BWA-MEM2 index files.
    fasta : str | None
        Reference FASTA used to build the index when `<index_path>.0123` does not exist.
    build_if_missing : bool
        Build a persistent index when `<index_path>.0123` does not exist. Requires `fasta`.
        Defaults to `True`.
*/
#[pyclass(name = "BwaMem2")]
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
        signature = (index_path, *, fasta=None, build_if_missing=true),
        text_signature = "(index_path, *, fasta=None, build_if_missing=True)",
    )]
    pub fn new(
        index_path: PathBuf,
        fasta: Option<PathBuf>,
        build_if_missing: bool,
    ) -> Result<Self> {
        let sentinel = path_with_suffix(&index_path, ".0123");
        let index = if sentinel.exists() {
            FMIndex::read(&index_path)?
        } else {
            if !build_if_missing {
                bail!(
                    "BWA-MEM2 index does not exist at '{}'. Provide a prebuilt index or set build_if_missing=True with fasta=...",
                    index_path.display()
                );
            }
            let fasta = fasta.context("fasta is required when build_if_missing=True")?;
            log::info!(
                "Creating BWA-MEM2 index for fasta: {:?} with prefix: {:?}",
                fasta,
                index_path
            );
            FMIndex::new(fasta, &index_path)?
        };
        Ok(BWAMEM2(BurrowsWheelerAligner::new(
            index,
            AlignerOpts::default(),
        )))
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
    pub fn set_logging_enabled(&mut self, enable: bool) {
        if enable {
            self.0.opts.enable_log();
        } else {
            self.0.opts.disable_log();
        }
    }
}

fn path_with_suffix(path: &std::path::Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    value.into()
}

/** The MiniBWA aligner.

    MiniBWA is a fast short-read aligner and successor to BWA-MEM. It is used to
    align reads to a minibwa index built from a reference genome.

    Parameters
    ----------
    index_prefix : str
        The path prefix for the MiniBWA index files, without the .l2b or .mbw extension.
    fasta : str | None
        Reference FASTA used to build the index when it cannot be loaded.
    build_if_missing : bool
        Build a persistent index when the index cannot be loaded. Requires `fasta`.
        Defaults to `True`.
    num_threads : int
        Number of threads used for index construction. Defaults to 8.
    preset : str | None
        Optional minibwa preset. Available presets are 'sr', 'lr', and 'adap'.
    methylation : bool
        Whether to load an index built for methylation-aware mapping.
*/
#[pyclass(name = "MiniBwa")]
#[repr(transparent)]
pub struct MINIBWA(MiniBwaSR);

impl Deref for MINIBWA {
    type Target = MiniBwaSR;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for MINIBWA {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[pymethods]
impl MINIBWA {
    #[new]
    #[pyo3(
        signature = (index_prefix, *, fasta=None, build_if_missing=true, num_threads=8, preset=None, methylation=false),
        text_signature = "(index_prefix, *, fasta=None, build_if_missing=True, num_threads=8, preset=None, methylation=False)",
    )]
    pub fn new(
        index_prefix: PathBuf,
        fasta: Option<PathBuf>,
        build_if_missing: bool,
        num_threads: i32,
        preset: Option<&str>,
        methylation: bool,
    ) -> Result<Self> {
        let options = match preset {
            Some(preset) => MiniBwaOptions::preset(preset)?,
            None => MiniBwaOptions::default(),
        };
        let index = match MiniBwaIndex::load(&index_prefix, methylation) {
            Ok(index) => index,
            Err(load_error) if !build_if_missing => bail!(
                "Failed to load MiniBWA index at '{}': {}. Provide a loadable index or set build_if_missing=True with fasta=...",
                index_prefix.display(),
                load_error
            ),
            Err(load_error) => {
                let fasta = fasta.with_context(|| {
                    format!(
                        "Failed to load MiniBWA index at '{}': {}. fasta is required when build_if_missing=True",
                        index_prefix.display(),
                        load_error
                    )
                })?;
                log::info!(
                    "Creating MiniBWA index for fasta: {:?} with prefix: {:?}",
                    fasta,
                    index_prefix
                );
                MiniBwaIndex::build(fasta, &index_prefix, num_threads, methylation).with_context(
                    || {
                        format!(
                            "Failed to build MiniBWA index at '{}' after loading failed: {}",
                            index_prefix.display(),
                            load_error
                        )
                    },
                )?
            }
        };
        Ok(MINIBWA(MiniBwaSR::new(index, options)?))
    }
}

/** The Minimap2 aligner.

    Minimap2 is a versatile aligner for long reads (Oxford Nanopore, PacBio),
    splice alignment, assembly-to-assembly alignment, and more.

    Parameters
    ----------
    index_path : str
        The path to the Minimap2 index file (.mmi).
    fasta : str | None
        Reference FASTA used to build `index_path` when it does not exist.
    build_if_missing : bool
        Build a persistent index when `index_path` does not exist. Requires `fasta`.
        Defaults to `True`.
    preset : str
        The minimap2 preset to use. Available presets:
        - Long Reads DNA Mapping:
          - 'map-ont': Oxford Nanopore genomic reads (default)
          - 'map-pb': PacBio CLR genomic reads
          - 'map-hifi': PacBio HiFi/CCS genomic reads
          - 'lr:hq': Long reads, high quality
        - Spliced / RNA-seq Alignment:
          - 'splice': Long-read spliced alignment (RNA-seq)
          - 'splice:hq': High-quality long-read spliced alignment
          - 'splice:sr': Short-read RNA-seq
        - Long Assembly to Reference Mapping:
          - 'asm5': Assembly-to-assembly alignment (divergence ~5%)
          - 'asm10': Assembly-to-assembly alignment (divergence ~10%)
          - 'asm20': Assembly-to-assembly alignment (divergence ~20%)
        - Short Reads Mapping:
          - 'short': Short single-end reads
          - 'sr': Short paired-end reads
        - All-vs-All Overlap Mapping:
          - 'ava-pb': PacBio all-vs-all overlap
          - 'ava-ont': ONT all-vs-all overlap
*/
#[pyclass(name = "Minimap2")]
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
        signature = (index_path, *, fasta=None, build_if_missing=true, preset="map-ont"),
        text_signature = "(index_path, *, fasta=None, build_if_missing=True, preset='map-ont')",
    )]
    pub fn new(
        index_path: PathBuf,
        fasta: Option<PathBuf>,
        build_if_missing: bool,
        preset: &str,
    ) -> Result<Self> {
        let preset = parse_minimap2_preset(preset)?;

        if is_fasta_file(&index_path) {
            bail!(
                "Minimap2 index_path must point to an index file, not a FASTA. Pass the output index path as index_path and use fasta=... with build_if_missing=True"
            );
        }

        if !index_path.exists() {
            if !build_if_missing {
                bail!(
                    "Minimap2 index does not exist at '{}'. Provide a prebuilt index or set build_if_missing=True with fasta=...",
                    index_path.display()
                );
            }
            let fasta = fasta.context("fasta is required when build_if_missing=True")?;
            let output_index = index_path
                .to_str()
                .context("Minimap2 index path must be valid UTF-8")?;
            log::info!(
                "Creating minimap2 index for fasta: {:?} with preset: {:?}",
                fasta,
                preset
            );
            minimap2::Aligner::builder()
                .preset(preset.clone())
                .with_index(&fasta, Some(output_index))
                .map_err(|error| anyhow::anyhow!("Failed to create minimap2 index: {}", error))?;
        }

        let opts = Minimap2Opts::new(index_path).with_preset(preset);
        Ok(Self(Minimap2Aligner::new(opts)?))
    }

    /// Get the currently configured preset name.
    ///
    /// Returns
    /// -------
    /// str | None
    ///     The preset name, or None if using default (map-ont).
    #[getter]
    pub fn get_preset(&self) -> Option<String> {
        self.0
            .get_opts()
            .preset()
            .map(|p| format!("{:?}", p).to_lowercase())
    }
}

fn is_fasta_file(path: &std::path::Path) -> bool {
    let path = if path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    {
        path.with_extension("")
    } else {
        path.to_path_buf()
    };

    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn"
            )
        })
}

#[pymodule]
pub(crate) fn register_aligners(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "aligners")?;

    m.add_class::<STAR>()?;
    m.add_class::<BWAMEM2>()?;
    m.add_class::<MINIBWA>()?;
    m.add_class::<MINIMAP2>()?;

    parent_module.add_submodule(&m)
}
