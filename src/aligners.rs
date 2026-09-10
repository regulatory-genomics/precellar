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
    collections::BTreeMap,
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
};

fn reference_sequences_from_star_index(index_path: &Path) -> Result<BTreeMap<String, usize>> {
    let path = index_path.join("chrNameLength.txt");
    let contents = std::fs::read_to_string(&path).with_context(|| {
        format!(
            "failed to read STAR reference sequences from '{}'",
            path.display()
        )
    })?;
    let mut references = BTreeMap::new();

    for (line_index, line) in contents.lines().enumerate() {
        let mut fields = line.split_whitespace();
        let Some(name) = fields.next() else {
            continue;
        };
        let length = fields
            .next()
            .with_context(|| format!("missing length at {}:{}", path.display(), line_index + 1))?
            .parse::<usize>()
            .with_context(|| {
                format!(
                    "invalid reference length at {}:{}",
                    path.display(),
                    line_index + 1
                )
            })?;
        if references.insert(name.to_owned(), length).is_some() {
            bail!(
                "duplicate reference sequence '{name}' in '{}'",
                path.display()
            );
        }
    }

    if references.is_empty() {
        bail!(
            "STAR reference sequence table is empty: '{}'",
            path.display()
        );
    }
    Ok(references)
}

fn summarize_references(references: &[String]) -> String {
    if references.is_empty() {
        "none".to_owned()
    } else {
        references
            .iter()
            .take(5)
            .cloned()
            .collect::<Vec<_>>()
            .join(", ")
    }
}

fn validate_reference_compatibility(index_path: &Path, header: &Header) -> Result<()> {
    let aligner_references = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| (name.to_string(), reference.length().get()))
        .collect::<BTreeMap<_, _>>();
    let star_references = reference_sequences_from_star_index(index_path)?;

    let only_aligner = aligner_references
        .iter()
        .filter(|(name, _)| !star_references.contains_key(*name))
        .map(|(name, length)| format!("{name}:{length}"))
        .collect::<Vec<_>>();
    let only_star = star_references
        .iter()
        .filter(|(name, _)| !aligner_references.contains_key(*name))
        .map(|(name, length)| format!("{name}:{length}"))
        .collect::<Vec<_>>();
    let length_mismatches = aligner_references
        .iter()
        .filter_map(|(name, aligner_length)| {
            star_references
                .get(name)
                .filter(|star_length| *star_length != aligner_length)
                .map(|star_length| format!("{name}:aligner={aligner_length},STAR={star_length}"))
        })
        .collect::<Vec<_>>();

    if !only_aligner.is_empty() || !only_star.is_empty() || !length_mismatches.is_empty() {
        bail!(
            "aligner and STAR transcriptome references are incompatible: aligner={}, STAR={}; only in aligner: [{}]; only in STAR: [{}]; length mismatches: [{}]. Rebuild both indexes from the same reference FASTA, and ensure the GTF uses the same chromosome names (for example, 'chr1' versus '1'). STAR index: '{}'",
            aligner_references.len(),
            star_references.len(),
            summarize_references(&only_aligner),
            summarize_references(&only_star),
            summarize_references(&length_mismatches),
            index_path.display(),
        );
    }
    log::info!(
        "Validated {} reference sequences between the aligner and STAR index '{}'",
        aligner_references.len(),
        index_path.display()
    );
    Ok(())
}

fn make_transcript_annotator(
    transcriptome: impl IntoIterator<Item = star_aligner::transcript::Transcript>,
    header: Header,
    strandness: Option<ChemistryStrandedness>,
) -> Result<TxAligner> {
    let transcripts = transcriptome
        .into_iter()
        .map(Transcript::try_from)
        .collect::<Result<Vec<_>>>()?;
    Ok(TxAligner::new(transcripts, header, strandness))
}

pub(crate) fn transcript_annotator_from_star_index(
    index_path: &Path,
    header: Header,
    strandness: Option<ChemistryStrandedness>,
) -> Result<TxAligner> {
    validate_reference_compatibility(index_path, &header)?;
    let transcriptome = star_aligner::transcript::Transcriptome::from_path(index_path)
        .with_context(|| {
            format!(
                "failed to load transcript annotation from STAR index '{}'",
                index_path.display()
            )
        })?;
    make_transcript_annotator(transcriptome.iter().cloned(), header, strandness)
}

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
    ) -> Result<Option<TxAligner>> {
        match self {
            AlignerRef::STAR(aligner) => {
                let transcriptome = aligner
                    .get_transcriptome()
                    .context("failed to load transcript annotation from STAR index")?;
                Ok(Some(make_transcript_annotator(
                    transcriptome.iter().cloned(),
                    self.header(),
                    strandness,
                )?))
            }
            AlignerRef::BWA(_) => Ok(None),
            AlignerRef::MiniBwa(_) => Ok(None),
            AlignerRef::Minimap2(_) => Ok(None),
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

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::{fs, num::NonZeroUsize};

    #[test]
    fn loads_transcript_annotation_independently_of_aligner() {
        let directory = tempfile::tempdir().unwrap();
        fs::write(
            directory.path().join("geneInfo.tab"),
            "1\ngene1\tGene 1\tprotein_coding\n",
        )
        .unwrap();
        fs::write(
            directory.path().join("transcriptInfo.tab"),
            "1\ntx1\t0\t9\t9\t1\t1\t0\t0\n",
        )
        .unwrap();
        fs::write(directory.path().join("exonInfo.tab"), "1\n0\t9\t0\n").unwrap();
        fs::write(directory.path().join("chrNameLength.txt"), "chr1\t100\n").unwrap();
        fs::write(directory.path().join("chrStart.txt"), "0\n100\n").unwrap();

        let header = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            )
            .build();
        let annotator = transcript_annotator_from_star_index(
            directory.path(),
            header,
            Some(ChemistryStrandedness::Reverse),
        )
        .unwrap();

        let transcripts: Vec<_> = annotator.transcripts().collect();
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].id, "tx1");
        assert_eq!(transcripts[0].gene_id, "gene1");
        assert_eq!(transcripts[0].gene_name, "Gene 1");
    }

    #[test]
    fn rejects_incompatible_aligner_and_star_references() {
        let directory = tempfile::tempdir().unwrap();
        fs::write(
            directory.path().join("chrNameLength.txt"),
            "chr1\t101\nchr2\t200\n",
        )
        .unwrap();

        let header = Header::builder()
            .add_reference_sequence(
                "1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(201).unwrap()),
            )
            .build();
        let error = validate_reference_compatibility(directory.path(), &header).unwrap_err();
        let message = error.to_string();

        assert!(message.contains("only in aligner: [1:100]"));
        assert!(message.contains("only in STAR: [chr1:101]"));
        assert!(message.contains("chr2:aligner=201,STAR=200"));
        assert!(message.contains("same reference FASTA"));
    }
}
