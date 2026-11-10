//! Terminal UCF sinks: aligned records go straight into UCF, with no
//! intermediate fragment TSV or BAM round trip.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use anyhow::{anyhow, bail, Result};
use noodles_sam as sam;
use precellar::fragment::{IntoFragOpts, IntoFragments};
use precellar::qc::QcFragment;
use precellar::ucf::sc_atac::{
    sort_snv_rows, strand_flags, Dictionary, ScAtacRow, ScAtacSnvRow,
};
use precellar::ucf::sc_rna::{
    classify_record, sort_rows as sort_sc_rna_rows, AlignmentPolicy, MissingBarcodePolicy,
    ScRnaDisposition, ScRnaRetention, ScRnaRow,
};
use pyo3::prelude::*;
use seqspec::utils::{create_file, Compression};
use ucf_core::BuildInfo;
use ucf_rna::{
    aligned_blocks_from_cigar, pack_rna_coordinate, RnaCigarOpV1, RnaError, ScRnaCellV1,
    ScRnaMetadataV1, ScRnaRowInput, ScRnaWriteOptions, ScRnaWriteRequest, ScRnaWriter,
    UNMAPPED_MATE_LOCUS_ID,
};
use ucf_sc_atac::{
    CellInfoV2, ChromInfoV1, ScAtacError, ScAtacMetadataV2, ScAtacRecordV2, ScAtacWriteOptions,
    ScAtacWriter, DEFAULT_SC_ATAC_SCHEMA_V2,
};

use crate::align::AlignProgressBar;
use crate::aligners::AlignerRef;
use crate::sinks::{fragment_file_header, AlignmentContext, AlignmentSink, OutputReport};

type Shared<T> = Rc<RefCell<T>>;

fn shared<T: Default>(value: T) -> Shared<T> {
    Rc::new(RefCell::new(value))
}

fn rna_chromosome_table(header: &sam::Header) -> Result<Vec<ucf_rna::RnaChromosomeV1>> {
    header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(index, (name, reference))| {
            Ok(ucf_rna::RnaChromosomeV1 {
                chrom_id: u64::try_from(index)?,
                name: name.to_string(),
                length_bp: reference.length().get() as u64,
            })
        })
        .collect()
}

fn chromosome_table(header: &sam::Header) -> Result<(HashMap<String, u32>, Vec<ChromInfoV1>)> {
    let mut ids = HashMap::new();
    let mut chroms = Vec::new();
    for (index, (name, reference)) in header.reference_sequences().iter().enumerate() {
        let chrom_id = u32::try_from(index)?;
        ids.insert(name.to_string(), chrom_id);
        chroms.push(ChromInfoV1 {
            chrom_id,
            name: name.to_string(),
            length_bp: reference.length().get() as u64,
        });
    }
    if chroms.is_empty() {
        bail!("the alignment header declares no reference sequences");
    }
    Ok((ids, chroms))
}

fn build_info(context: &AlignmentContext) -> BuildInfo {
    BuildInfo {
        created_by: format!("precellar {}", env!("CARGO_PKG_VERSION")),
        created_at_utc: epoch_seconds(),
        source_paths: context.source_paths.clone(),
    }
}

fn epoch_nanos() -> u128 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or_default()
}

fn epoch_seconds() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_secs().to_string())
        .unwrap_or_default()
}

fn default_snv_path(output: &Path) -> PathBuf {
    PathBuf::from(format!("{}.snv.tsv.zst", output.display()))
}

fn record_error(slot: &Shared<Option<anyhow::Error>>, error: impl Into<anyhow::Error>) {
    let mut slot = slot.borrow_mut();
    if slot.is_none() {
        *slot = Some(error.into());
    }
}

fn single_assay_id(assays: &[seqspec::Assay]) -> Result<Option<String>> {
    match assays {
        [assay] => Ok(Some(assay.assay_id.clone())),
        _ => bail!(
            "UCF sinks write one assay per output; found {} assays",
            assays.len()
        ),
    }
}

/// One output in the pair of files published by the scATAC sink.
///
/// The staging file is removed on drop unless it was published. The caller
/// publishes the main UCF and optional sidecar together with rollback for
/// ordinary rename errors.
struct AtomicOutput {
    output: PathBuf,
    staging: PathBuf,
    committed: bool,
}

impl AtomicOutput {
    fn new(output: &Path) -> Result<Self> {
        let parent = output
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let name = output
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("sidecar");
        let staging = parent.join(format!(
            ".{name}.{}-{}.tmp",
            std::process::id(),
            epoch_nanos()
        ));
        Ok(Self {
            output: output.to_path_buf(),
            staging,
            committed: false,
        })
    }

    fn path(&self) -> &Path {
        &self.staging
    }

}

impl Drop for AtomicOutput {
    fn drop(&mut self) {
        if !self.committed {
            let _ = std::fs::remove_file(&self.staging);
        }
    }
}

fn comparable_output_path(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                normalized.pop();
            }
            _ => normalized.push(component.as_os_str()),
        }
    }
    if normalized.exists() {
        return Ok(std::fs::canonicalize(normalized)?);
    }
    if let (Some(parent), Some(name)) = (normalized.parent(), normalized.file_name()) {
        if let Ok(parent) = std::fs::canonicalize(parent) {
            return Ok(parent.join(name));
        }
    }
    Ok(normalized)
}

fn validate_output_paths(main: &Path, sidecar: &Path) -> Result<()> {
    if comparable_output_path(main)? == comparable_output_path(sidecar)? {
        bail!("UCF and SNV sidecar outputs must be different paths");
    }
    Ok(())
}

fn write_sidecar_line(writer: &mut dyn Write, line: &str) -> std::io::Result<()> {
    writeln!(writer, "{line}")
}

fn backup_existing_output(output: &Path) -> Result<Option<PathBuf>> {
    let metadata = match std::fs::symlink_metadata(output) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    if !metadata.file_type().is_file() && !metadata.file_type().is_symlink() {
        bail!("output path is not a regular file: {}", output.display());
    }
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = output
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("output");
    let backup = parent.join(format!(
        ".{name}.backup-{}-{}",
        std::process::id(),
        epoch_nanos()
    ));
    std::fs::rename(output, &backup)
        .map_err(|error| anyhow!("cannot stage existing {}: {error}", output.display()))?;
    Ok(Some(backup))
}

fn restore_output(
    output: &Path,
    backup: Option<PathBuf>,
    published: bool,
) -> Result<()> {
    if published && output.exists() {
        std::fs::remove_file(output)
            .map_err(|error| anyhow!("cannot remove published {}: {error}", output.display()))?;
    }
    if let Some(backup) = backup {
        std::fs::rename(&backup, output).map_err(|error| {
            anyhow!("cannot restore previous {}: {error}", output.display())
        })?;
    }
    Ok(())
}

fn publish_output_pair(
    main: &mut AtomicOutput,
    mut sidecar: Option<&mut AtomicOutput>,
) -> Result<()> {
    if let Some(sidecar) = sidecar.as_deref() {
        validate_output_paths(&main.output, &sidecar.output)?;
    }

    let mut main_backup = None;
    let mut sidecar_backup = None;
    let mut main_published = false;
    let mut sidecar_published = false;
    let result = (|| {
        main_backup = backup_existing_output(&main.output)?;
        if let Some(sidecar) = sidecar.as_deref_mut() {
            sidecar_backup = backup_existing_output(&sidecar.output)?;
            std::fs::rename(&sidecar.staging, &sidecar.output).map_err(|error| {
                anyhow!("cannot publish {}: {error}", sidecar.output.display())
            })?;
            sidecar.committed = true;
            sidecar_published = true;
        }
        std::fs::rename(&main.staging, &main.output)
            .map_err(|error| anyhow!("cannot publish {}: {error}", main.output.display()))?;
        main.committed = true;
        main_published = true;
        Ok(())
    })();

    if let Err(error) = result {
        let rollback = restore_output(&main.output, main_backup, main_published)
            .and_then(|_| match sidecar.as_deref() {
                Some(sidecar) => restore_output(
                    &sidecar.output,
                    sidecar_backup,
                    sidecar_published,
                ),
                None => Ok(()),
            });
        return match rollback {
            Ok(()) => Err(error),
            Err(rollback_error) => Err(anyhow!("{error}; rollback failed: {rollback_error}")),
        };
    }

    if let Some(backup) = main_backup {
        let _ = std::fs::remove_file(backup);
    }
    if let Some(backup) = sidecar_backup {
        let _ = std::fs::remove_file(backup);
    }
    Ok(())
}

fn parse_alignment_policy(value: &str) -> Result<AlignmentPolicy> {
    match value {
        "primary" => Ok(AlignmentPolicy::Primary),
        "all" => Ok(AlignmentPolicy::All),
        _ => bail!("alignments must be 'primary' or 'all'"),
    }
}

fn parse_missing_barcode_policy(value: &str) -> Result<MissingBarcodePolicy> {
    match value {
        "drop" => Ok(MissingBarcodePolicy::Drop),
        "error" => Ok(MissingBarcodePolicy::Error),
        _ => bail!("on_missing_barcode must be 'drop' or 'error'"),
    }
}

fn add_retention(metrics: &mut serde_json::Map<String, serde_json::Value>, retention: &ScRnaRetention) {
    metrics.insert("records_stored".to_owned(), retention.stored.into());
    metrics.insert("records_unmapped".to_owned(), retention.unmapped.into());
    metrics.insert("records_secondary".to_owned(), retention.secondary.into());
    metrics.insert(
        "records_supplementary".to_owned(),
        retention.supplementary.into(),
    );
    metrics.insert(
        "records_missing_barcode".to_owned(),
        retention.missing_barcode.into(),
    );
    metrics.insert(
        "records_missing_reference".to_owned(),
        retention.missing_reference.into(),
    );
    metrics.insert("records_empty_cigar".to_owned(), retention.empty_cigar.into());
}

#[pyclass(module = "precellar.sinks", name = "ScAtacUcfSink")]
pub(crate) struct ScAtacUcfSink {
    output: PathBuf,
    shift_left: i64,
    shift_right: i64,
    compute_snv: bool,
    snv_output: Option<PathBuf>,
    block_size: u32,
    chunk_rows: usize,
    sort_chunk_size: usize,
    zstd_level: i32,
}

#[pymethods]
impl ScAtacUcfSink {
    /// Write deduplicated fragments directly as an scATAC UCF v2 file.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination UCF path. The file appears only after a successful write.
    /// shift_left : int, keyword-only, default=4
    ///     Shift applied to the left end of paired-end fragments.
    /// shift_right : int, keyword-only, default=-5
    ///     Shift applied to the right end of paired-end fragments.
    /// compute_snv : bool, keyword-only, default=False
    ///     Also write the fragment SNV sidecar.
    /// snv_output : pathlib.Path | str, optional
    ///     Sidecar path. Defaults to `<output>.snv.tsv.zst`.
    /// block_size : int, keyword-only, default=65536
    ///     Genomic block size stored in UCF chunk summaries.
    /// chunk_rows : int, keyword-only, default=65536
    ///     Target rows per independently decodable UCF chunk.
    /// sort_chunk_size : int, keyword-only, default=1000000
    ///     Records buffered in memory before an external-sort run is spilled.
    /// zstd_level : int, keyword-only, default=2
    ///     UCF payload compression level.
    ///
    /// Notes
    /// -----
    /// One assay per output. The assay id is used as the UCF sample id and the
    /// cell dictionary contains observed barcodes only. The SNV sidecar is a
    /// separate file so existing scATAC UCF readers keep working unchanged.
    #[new]
    #[pyo3(signature = (
        output, *, shift_left=4, shift_right=-5, compute_snv=false,
        snv_output=None, block_size=65536, chunk_rows=65536,
        sort_chunk_size=1_000_000, zstd_level=2
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        output: PathBuf,
        shift_left: i64,
        shift_right: i64,
        compute_snv: bool,
        snv_output: Option<PathBuf>,
        block_size: u32,
        chunk_rows: usize,
        sort_chunk_size: usize,
        zstd_level: i32,
    ) -> Result<Self> {
        if block_size == 0 {
            bail!("block_size must be greater than zero");
        }
        if chunk_rows == 0 {
            bail!("chunk_rows must be greater than zero");
        }
        if sort_chunk_size == 0 {
            bail!("sort_chunk_size must be greater than zero");
        }
        if !(1..=22).contains(&zstd_level) {
            bail!("zstd_level must be in 1..=22");
        }
        if snv_output.is_some() && !compute_snv {
            bail!("snv_output requires compute_snv=True");
        }
        if compute_snv {
            let sidecar = snv_output
                .clone()
                .unwrap_or_else(|| default_snv_path(&output));
            validate_output_paths(&output, &sidecar)?;
        }
        Ok(Self {
            output,
            shift_left,
            shift_right,
            compute_snv,
            snv_output,
            block_size,
            chunk_rows,
            sort_chunk_size,
            zstd_level,
        })
    }
}

impl AlignmentSink for ScAtacUcfSink {
    fn sample_id(&self, assays: &[seqspec::Assay]) -> Result<Option<String>> {
        single_assay_id(assays)
    }

    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        let (chrom_ids, chromosomes) = chromosome_table(header)?;
        let cells: Shared<Dictionary> = shared(Dictionary::default());
        let qc: Shared<QcFragment> = shared(QcFragment::default());
        let error: Shared<Option<anyhow::Error>> = shared(None);
        {
            let mut qc = qc.borrow_mut();
            for mito in &context.mito_dna {
                qc.add_mito_dna(mito);
            }
        }
        let opts = IntoFragOpts {
            shift_left: self.shift_left,
            shift_right: self.shift_right,
            temp_dir: context.temp_dir.clone(),
            compute_snv: self.compute_snv,
            chunk_size: self.sort_chunk_size,
            ..Default::default()
        };

        let unique = (&mut *alignments).into_fragments(header, opts);
        let records = {
            let cells = Rc::clone(&cells);
            let qc = Rc::clone(&qc);
            let error_slot = Rc::clone(&error);
            let mut groups = (&unique).into_iter();
            let mut buffer: std::vec::IntoIter<ScAtacSnvRow> = Vec::new().into_iter();
            let rows = std::iter::from_fn(move || loop {
                if let Some(row) = buffer.next() {
                    return Some(row);
                }
                let Some(group) = groups.next() else {
                    return None;
                };
                let mut out = Vec::with_capacity(group.len());
                for fragment in &group {
                    if let Err(error) = py.check_signals() {
                        record_error(&error_slot, error);
                        return None;
                    }
                    let Some(&chrom_id) = chrom_ids.get(&fragment.chrom) else {
                        record_error(
                            &error_slot,
                            anyhow!("fragment references unknown chromosome {}", fragment.chrom),
                        );
                        return None;
                    };
                    let (Ok(start), Ok(end)) =
                        (u32::try_from(fragment.start), u32::try_from(fragment.end))
                    else {
                        record_error(
                            &error_slot,
                            anyhow!("fragment coordinate exceeds u32: {fragment}"),
                        );
                        return None;
                    };
                    qc.borrow_mut().update(fragment);
                    let cell_id = cells.borrow_mut().id(fragment.barcode.as_deref().unwrap_or("."));
                    out.push(ScAtacSnvRow {
                        row: ScAtacRow {
                            chrom_id,
                            start,
                            end,
                            cell_id,
                            multiplicity: fragment.count,
                            flags: strand_flags(fragment.strand),
                        },
                        line: if self.compute_snv {
                            fragment.to_string()
                        } else {
                            String::new()
                        },
                    });
                }
                buffer = out.into_iter();
            });
            sort_snv_rows(rows, context.temp_dir.clone(), self.sort_chunk_size)?
        };
        if let Some(error) = error.borrow_mut().take() {
            return Err(error);
        }

        let cell_entries = cells
            .borrow()
            .entries()
            .iter()
            .enumerate()
            .map(|(index, barcode)| {
                Ok(CellInfoV2 {
                    cell_id: u32::try_from(index)?,
                    sample_id: context.sample_id.clone(),
                    barcode: barcode.clone(),
                    cell_key: format!("{}:{barcode}", context.sample_id),
                    order_rank: None,
                    tile_id: None,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        let sidecar_path = if self.compute_snv {
            Some(
                self.snv_output
                    .clone()
                    .unwrap_or_else(|| default_snv_path(&self.output)),
            )
        } else {
            None
        };
        let mut main_output = AtomicOutput::new(&self.output)?;
        let mut sidecar = match &sidecar_path {
            Some(path) => Some(AtomicOutput::new(path)?),
            None => None,
        };
        let mut snv_writer = match &sidecar {
            Some(sidecar) => Some(create_file(
                sidecar.path().to_path_buf(),
                Some(Compression::Zstd),
                None,
                context.num_threads as u32,
            )?),
            None => None,
        };
        if let Some(writer) = snv_writer.as_mut() {
            writeln!(
                writer,
                "{}",
                fragment_file_header(self.compute_snv, self.shift_left, self.shift_right)
            )?;
        }

        let metadata = ScAtacMetadataV2::new(
            DEFAULT_SC_ATAC_SCHEMA_V2,
            2,
            self.block_size,
            0,
            0,
            chromosomes,
        );
        let mut records = records;
        let stream = std::iter::from_fn(|| {
            let next = match records.next()? {
                Ok(next) => next,
                Err(error) => return Some(Err(ScAtacError::Invalid(error.to_string()))),
            };
            if let Some(sidecar) = snv_writer.as_mut() {
                if let Err(error) = write_sidecar_line(sidecar.as_mut(), &next.line) {
                    return Some(Err(ScAtacError::Invalid(format!(
                        "write SNV sidecar: {error}"
                    ))));
                }
            }
            Some(Ok(ScAtacRecordV2::new(
                next.row.chrom_id,
                next.row.start,
                next.row.end,
                next.row.cell_id,
                next.row.multiplicity,
                next.row.flags,
            )))
        });
        ScAtacWriter::write_file_v2_iter_fallible_with_options_auto_counts_zstd(
            main_output.path(),
            &metadata,
            &cell_entries,
            stream,
            build_info(context),
            ScAtacWriteOptions {
                chunk_target_rows: self.chunk_rows,
                span_blocks: false,
            },
            self.zstd_level,
        )?;
        if let Some(error) = error.borrow_mut().take() {
            return Err(error);
        }
        if let Some(mut writer) = snv_writer.take() {
            writer.flush()?;
            drop(writer);
            if let Some(sidecar) = sidecar.as_ref() {
                std::fs::File::open(sidecar.path())?.sync_all()?;
            }
        }
        std::fs::File::open(main_output.path())?.sync_all()?;
        publish_output_pair(&mut main_output, sidecar.as_mut())?;

        let serde_json::Value::Object(mut metrics) =
            std::mem::take(&mut *qc.borrow_mut()).into()
        else {
            bail!("fragment QC metrics are not a JSON object");
        };
        metrics.insert("cells".to_owned(), serde_json::Value::from(cell_entries.len()));
        metrics.insert(
            "sort_chunk_size".to_owned(),
            serde_json::Value::from(self.sort_chunk_size),
        );
        if let Some(path) = &sidecar_path {
            metrics.insert(
                "snv_sidecar".to_owned(),
                serde_json::Value::String(path.display().to_string()),
            );
        }
        Ok(Some(OutputReport {
            name: "sc_atac_ucf".to_owned(),
            metrics: serde_json::Value::Object(metrics),
        }))
    }
}

#[pyclass(module = "precellar.sinks", name = "ScRnaUcfSink")]
pub(crate) struct ScRnaUcfSink {
    output: PathBuf,
    alignments: AlignmentPolicy,
    on_missing_barcode: MissingBarcodePolicy,
    chunk_rows: usize,
    sort_chunk_size: usize,
    zstd_level: i32,
}

#[pymethods]
impl ScRnaUcfSink {
    /// Write mapped alignments directly as a cell-aware scRNA UCF.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination UCF path. The file appears only after a successful write.
    /// alignments : {"primary", "all"}, keyword-only, default="primary"
    ///     Store primary alignments only, or every mapped alignment. Excluded
    ///     secondary/supplementary alignments are counted in the report.
    /// on_missing_barcode : {"drop", "error"}, keyword-only, default="drop"
    ///     Skip mapped records without a CB tag and count them, or fail.
    /// chunk_rows : int, keyword-only, default=65536
    ///     Target rows per independently decodable UCF chunk.
    /// sort_chunk_size : int, keyword-only, default=1000000
    ///     Records buffered in memory before an external-sort run is spilled.
    /// zstd_level : int, keyword-only, default=2
    ///     UCF payload compression level.
    ///
    /// Notes
    /// -----
    /// One assay per output. Stores coordinates, CIGAR, aligned blocks, flags,
    /// NH, MAPQ, cell barcode, and exact UMI. Read sequence, quality, QNAME, and
    /// arbitrary SAM tags stay in the source BAM. The family is
    /// `sc_rna_alignment`, which is distinct from bulk `rna_alignment`.
    #[new]
    #[pyo3(signature = (
        output, *, alignments="primary", on_missing_barcode="drop",
        chunk_rows=65536, sort_chunk_size=1_000_000, zstd_level=2
    ))]
    fn new(
        output: PathBuf,
        alignments: &str,
        on_missing_barcode: &str,
        chunk_rows: usize,
        sort_chunk_size: usize,
        zstd_level: i32,
    ) -> Result<Self> {
        if chunk_rows == 0 {
            bail!("chunk_rows must be greater than zero");
        }
        if sort_chunk_size == 0 {
            bail!("sort_chunk_size must be greater than zero");
        }
        if !(1..=22).contains(&zstd_level) {
            bail!("zstd_level must be in 1..=22");
        }
        Ok(Self {
            output,
            alignments: parse_alignment_policy(alignments)?,
            on_missing_barcode: parse_missing_barcode_policy(on_missing_barcode)?,
            chunk_rows,
            sort_chunk_size,
            zstd_level,
        })
    }
}

#[derive(Default)]
struct ScRnaTotals {
    blocks: u64,
    cigars: u64,
    max_span: u32,
}

pub(crate) fn register_ucf_sinks(sinks: &Bound<'_, PyModule>) -> PyResult<()> {
    sinks.add_class::<ScAtacUcfSink>()?;
    sinks.add_class::<ScRnaUcfSink>()?;
    Ok(())
}

impl AlignmentSink for ScRnaUcfSink {
    fn sample_id(&self, assays: &[seqspec::Assay]) -> Result<Option<String>> {
        single_assay_id(assays)
    }

    fn consume<'py, 'stream>(
        &self,
        py: Python<'py>,
        header: &sam::Header,
        alignments: &mut AlignProgressBar<'stream, AlignerRef<'py>>,
        context: &mut AlignmentContext,
    ) -> Result<Option<OutputReport>> {
        let chromosomes = rna_chromosome_table(header)?;
        let cells: Shared<Dictionary> = shared(Dictionary::default());
        let umis: Shared<Dictionary> = shared(Dictionary::default());
        let totals: Shared<ScRnaTotals> = shared(ScRnaTotals::default());
        let retention: Shared<ScRnaRetention> = shared(ScRnaRetention::default());
        let error: Shared<Option<anyhow::Error>> = shared(None);

        let rows = {
            let cells = Rc::clone(&cells);
            let umis = Rc::clone(&umis);
            let totals = Rc::clone(&totals);
            let retention = Rc::clone(&retention);
            let error_slot = Rc::clone(&error);
            let mut records: std::vec::IntoIter<ScRnaRow> = Vec::new().into_iter();
            let rows = std::iter::from_fn(move || loop {
                if let Some(row) = records.next() {
                    return Some(row);
                }
                let Some(chunk) = alignments.next() else {
                    return None;
                };
                if let Err(err) = py.check_signals() {
                    record_error(&error_slot, err);
                    return None;
                }
                let mut out = Vec::new();
                for (read1, read2) in chunk {
                    for record in read1.into_iter().chain(read2) {
                        for alignment in record.iter() {
                            match classify_record(
                                alignment,
                                |barcode| cells.borrow_mut().id(barcode),
                                |umi| umis.borrow_mut().id(umi),
                                self.alignments,
                                self.on_missing_barcode,
                            ) {
                                Ok(disposition) => {
                                    retention.borrow_mut().record(&disposition);
                                    if let ScRnaDisposition::Stored(row) = disposition {
                                        out.push(row);
                                    }
                                }
                                Err(err) => {
                                    record_error(&error_slot, err);
                                    return None;
                                }
                            }
                        }
                    }
                }
                {
                    let mut totals = totals.borrow_mut();
                    for row in &out {
                        totals.blocks = totals.blocks.saturating_add(u64::from(row.block_count));
                        totals.cigars = totals.cigars.saturating_add(row.cigar.len() as u64);
                        totals.max_span = totals.max_span.max(row.reference_span);
                    }
                }
                records = out.into_iter();
            });
            sort_sc_rna_rows(rows, context.temp_dir.clone(), self.sort_chunk_size)?
        };
        if let Some(error) = error.borrow_mut().take() {
            return Err(error);
        }

        let cell_entries = cells
            .borrow()
            .entries()
            .iter()
            .enumerate()
            .map(|(index, barcode)| {
                Ok(ScRnaCellV1::new(
                    u32::try_from(index)?,
                    context.sample_id.clone(),
                    barcode.clone(),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let umi_entries: Vec<String> = umis.borrow().entries().to_vec();
        let totals = totals.borrow();
        let retention = retention.borrow();
        let metadata = ScRnaMetadataV1 {
            schema: "sc_rna_alignment_v1".to_owned(),
            schema_version: 1,
            source_format: "precellar".to_owned(),
            total_record_count: retention.stored + retention.excluded(),
            stored_alignment_count: retention.stored,
            unmapped_record_count: retention.unmapped,
            total_block_count: totals.blocks,
            total_cigar_op_count: totals.cigars,
            max_alignment_span_bp: totals.max_span,
            cell_count: cell_entries.len() as u64,
            umi_count: umi_entries.len() as u64,
            chromosomes,
        };

        let mut records = rows;
        let stream = std::iter::from_fn(|| {
            let row = match records.next()? {
                Ok(row) => row,
                Err(err) => return Some(Err(RnaError::Invalid(err.to_string()))),
            };
            let mut cigar_ops = Vec::with_capacity(row.cigar.len());
            for &(len, code) in &row.cigar {
                match RnaCigarOpV1::new(len, code) {
                    Ok(op) => cigar_ops.push(op),
                    Err(err) => return Some(Err(err)),
                }
            }
            if let Err(err) = aligned_blocks_from_cigar(&cigar_ops) {
                return Some(Err(err));
            }
            Some(Ok(ScRnaRowInput {
                locus_id: pack_rna_coordinate(u64::from(row.chrom_id), u64::from(row.start)),
                mate_locus_id: match (row.mate_chrom_id, row.mate_start) {
                    (Some(chrom), Some(start)) => {
                        pack_rna_coordinate(u64::from(chrom), u64::from(start))
                    }
                    _ => UNMAPPED_MATE_LOCUS_ID,
                },
                template_length: row.template_length,
                sam_flags: row.flags,
                nh: row.nh,
                mapq: row.mapq,
                cell_id: row.cell_id,
                umi_id: row.umi_id,
                cigar: row.cigar,
            }))
        });
        ScRnaWriter::write_file(
            &self.output,
            ScRnaWriteRequest {
                metadata: &metadata,
                cells: &cell_entries,
                umis: &umi_entries,
                build: build_info(context),
                options: ScRnaWriteOptions {
                    chunk_rows: self.chunk_rows,
                    zstd_level: self.zstd_level,
                },
            },
            stream,
        )?;

        let mut metrics = serde_json::Map::new();
        metrics.insert("cells".to_owned(), cell_entries.len().into());
        metrics.insert("umis".to_owned(), umi_entries.len().into());
        metrics.insert(
            "alignments".to_owned(),
            serde_json::Value::String(match self.alignments {
                AlignmentPolicy::Primary => "primary",
                AlignmentPolicy::All => "all",
            }
            .to_owned()),
        );
        metrics.insert(
            "on_missing_barcode".to_owned(),
            serde_json::Value::String(match self.on_missing_barcode {
                MissingBarcodePolicy::Drop => "drop",
                MissingBarcodePolicy::Error => "error",
            }
            .to_owned()),
        );
        metrics.insert(
            "sort_chunk_size".to_owned(),
            serde_json::Value::from(self.sort_chunk_size),
        );
        add_retention(&mut metrics, &retention);
        Ok(Some(OutputReport {
            name: "sc_rna_ucf".to_owned(),
            metrics: serde_json::Value::Object(metrics),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::Other, "injected sidecar failure"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn rejects_identical_ucf_and_snv_paths() {
        assert!(validate_output_paths(
            Path::new("results/sample.ucf"),
            Path::new("./results/./sample.ucf"),
        )
        .is_err());
    }

    #[test]
    fn sidecar_write_failure_is_returned() {
        let mut writer = FailingWriter;
        let error = write_sidecar_line(&mut writer, "fragment").unwrap_err();
        assert!(error.to_string().contains("injected sidecar failure"));
    }

    #[test]
    fn publishing_pair_rolls_back_if_main_publish_fails() {
        let directory = tempfile::tempdir().unwrap();
        let main_output = directory.path().join("main.ucf");
        let sidecar_output = directory.path().join("main.snv.tsv.zst");
        std::fs::create_dir(&main_output).unwrap();
        std::fs::write(&sidecar_output, b"old sidecar").unwrap();

        let (main_stage, sidecar_stage, publish_error) = {
            let mut main = AtomicOutput::new(&main_output).unwrap();
            let mut sidecar = AtomicOutput::new(&sidecar_output).unwrap();
            let main_stage = main.path().to_owned();
            let sidecar_stage = sidecar.path().to_owned();
            std::fs::write(sidecar.path(), b"new sidecar").unwrap();
            let publish_error = publish_output_pair(&mut main, Some(&mut sidecar)).unwrap_err();
            (main_stage, sidecar_stage, publish_error)
        };
        assert!(publish_error.to_string().contains("not a regular file"));
        assert!(main_output.is_dir());
        assert_eq!(std::fs::read(&sidecar_output).unwrap(), b"old sidecar");
        assert!(!main_stage.exists());
        assert!(!sidecar_stage.exists());
    }
}
