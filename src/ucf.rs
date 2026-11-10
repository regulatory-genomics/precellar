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
use precellar::ucf::sc_rna::{row_from_record, sort_rows as sort_sc_rna_rows, ScRnaRow};
use pyo3::prelude::*;
use seqspec::utils::{create_file, Compression};
use ucf_core::BuildInfo;
use ucf_rna::{
    aligned_blocks_from_cigar, pack_rna_coordinate, RnaCigarOpV1, ScRnaCellV1, ScRnaMetadataV1,
    ScRnaRowInput, ScRnaWriteOptions, ScRnaWriteRequest, ScRnaWriter, UNMAPPED_MATE_LOCUS_ID,
};
use ucf_sc_atac::{
    CellInfoV2, ChromInfoV1, ScAtacMetadataV2, ScAtacRecordV2, ScAtacWriteOptions, ScAtacWriter,
    DEFAULT_SC_ATAC_SCHEMA_V2,
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

fn build_info() -> BuildInfo {
    BuildInfo {
        created_by: format!("precellar {}", env!("CARGO_PKG_VERSION")),
        created_at_utc: "0".to_owned(),
        source_paths: Vec::new(),
    }
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

#[pyclass(module = "precellar.sinks", name = "ScAtacUcfSink")]
pub(crate) struct ScAtacUcfSink {
    output: PathBuf,
    shift_left: i64,
    shift_right: i64,
    compute_snv: bool,
    snv_output: Option<PathBuf>,
    block_size: u32,
    chunk_rows: usize,
    zstd_level: i32,
}

#[pymethods]
impl ScAtacUcfSink {
    /// Write deduplicated fragments directly as an scATAC UCF v2 file.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination UCF path.
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
        snv_output=None, block_size=65536, chunk_rows=65536, zstd_level=2
    ))]
    fn new(
        output: PathBuf,
        shift_left: i64,
        shift_right: i64,
        compute_snv: bool,
        snv_output: Option<PathBuf>,
        block_size: u32,
        chunk_rows: usize,
        zstd_level: i32,
    ) -> Result<Self> {
        if block_size == 0 {
            bail!("block_size must be greater than zero");
        }
        if chunk_rows == 0 {
            bail!("chunk_rows must be greater than zero");
        }
        if !(1..=22).contains(&zstd_level) {
            bail!("zstd_level must be in 1..=22");
        }
        if snv_output.is_some() && !compute_snv {
            bail!("snv_output requires compute_snv=True");
        }
        Ok(Self {
            output,
            shift_left,
            shift_right,
            compute_snv,
            snv_output,
            block_size,
            chunk_rows,
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
                        record_error(&error_slot,
                            anyhow!("fragment references unknown chromosome {}", fragment.chrom),
                        );
                        return None;
                    };
                    let (Ok(start), Ok(end)) =
                        (u32::try_from(fragment.start), u32::try_from(fragment.end))
                    else {
                        record_error(&error_slot,
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
            sort_snv_rows(rows, context.temp_dir.clone(), 10_000_000)
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
        let mut snv_writer = match &sidecar_path {
            Some(path) => {
                let mut writer = create_file(
                    path.clone(),
                    Some(Compression::Zstd),
                    None,
                    context.num_threads as u32,
                )?;
                writeln!(
                    writer,
                    "{}",
                    fragment_file_header(self.compute_snv, self.shift_left, self.shift_right)
                )?;
                Some(writer)
            }
            None => None,
        };

        let metadata = ScAtacMetadataV2::new(
            DEFAULT_SC_ATAC_SCHEMA_V2,
            2,
            self.block_size,
            0,
            0,
            chromosomes,
        );
        let mut records = records;
        let error_slot = Rc::clone(&error);
        let stream = std::iter::from_fn(|| {
            let next = records.next()?;
            if let Some(writer) = snv_writer.as_mut() {
                if let Err(error) = writeln!(writer, "{}", next.line) {
                    record_error(&error_slot, error);
                    return None;
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
            &self.output,
            &metadata,
            &cell_entries,
            stream,
            build_info(),
            ScAtacWriteOptions {
                chunk_target_rows: self.chunk_rows,
                span_blocks: false,
            },
            self.zstd_level,
        )?;
        if let Some(error) = error.borrow_mut().take() {
            return Err(error);
        }
        if let Some(mut writer) = snv_writer {
            writer.flush()?;
        }

        let mut metrics: serde_json::Value =
            std::mem::take(&mut *qc.borrow_mut()).into();
        metrics["cells"] = serde_json::Value::from(cell_entries.len());
        if let Some(path) = &sidecar_path {
            metrics["snv_sidecar"] = serde_json::Value::String(path.display().to_string());
        }
        Ok(Some(OutputReport {
            name: "sc_atac_ucf".to_owned(),
            metrics,
        }))
    }
}

#[pyclass(module = "precellar.sinks", name = "ScRnaUcfSink")]
pub(crate) struct ScRnaUcfSink {
    output: PathBuf,
    chunk_rows: usize,
    zstd_level: i32,
}

#[pymethods]
impl ScRnaUcfSink {
    /// Write mapped alignments directly as a cell-aware scRNA UCF.
    ///
    /// Parameters
    /// ----------
    /// output : pathlib.Path | str
    ///     Destination UCF path.
    /// chunk_rows : int, keyword-only, default=65536
    ///     Target rows per independently decodable UCF chunk.
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
    #[pyo3(signature = (output, *, chunk_rows=65536, zstd_level=2))]
    fn new(output: PathBuf, chunk_rows: usize, zstd_level: i32) -> Result<Self> {
        if chunk_rows == 0 {
            bail!("chunk_rows must be greater than zero");
        }
        if !(1..=22).contains(&zstd_level) {
            bail!("zstd_level must be in 1..=22");
        }
        Ok(Self {
            output,
            chunk_rows,
            zstd_level,
        })
    }
}

#[derive(Default)]
struct ScRnaTotals {
    rows: u64,
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
        let error: Shared<Option<anyhow::Error>> = shared(None);

        let rows = {
            let cells = Rc::clone(&cells);
            let umis = Rc::clone(&umis);
            let totals = Rc::clone(&totals);
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
                            match row_from_record(
                                alignment,
                                |barcode| cells.borrow_mut().id(barcode),
                                |umi| umis.borrow_mut().id(umi),
                            ) {
                                Ok(Some(row)) => out.push(row),
                                Ok(None) => {}
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
                        totals.rows = totals.rows.saturating_add(1);
                        totals.blocks = totals.blocks.saturating_add(u64::from(row.block_count));
                        totals.cigars = totals.cigars.saturating_add(row.cigar.len() as u64);
                        totals.max_span = totals.max_span.max(row.reference_span);
                    }
                }
                records = out.into_iter();
            });
            sort_sc_rna_rows(rows, context.temp_dir.clone(), 10_000_000)
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
        let metadata = ScRnaMetadataV1 {
            schema: "sc_rna_alignment_v1".to_owned(),
            schema_version: 1,
            source_format: "precellar".to_owned(),
            total_record_count: totals.rows,
            stored_alignment_count: totals.rows,
            unmapped_record_count: 0,
            total_block_count: totals.blocks,
            total_cigar_op_count: totals.cigars,
            max_alignment_span_bp: totals.max_span,
            cell_count: cell_entries.len() as u64,
            umi_count: umi_entries.len() as u64,
            chromosomes,
        };

        let mut records = rows;
        let error_slot = Rc::clone(&error);
        let stream = std::iter::from_fn(|| {
            let row = records.next()?;
            let mut cigar_ops = Vec::with_capacity(row.cigar.len());
            for &(len, code) in &row.cigar {
                match RnaCigarOpV1::new(len, code) {
                    Ok(op) => cigar_ops.push(op),
                    Err(err) => {
                        record_error(&error_slot, err);
                        return None;
                    }
                }
            }
            if let Err(err) = aligned_blocks_from_cigar(&cigar_ops) {
                record_error(&error_slot, err);
                return None;
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
                build: build_info(),
                options: ScRnaWriteOptions {
                    chunk_rows: self.chunk_rows,
                    zstd_level: self.zstd_level,
                },
            },
            stream,
        )?;
        if let Some(error) = error.borrow_mut().take() {
            return Err(error);
        }
        let mut metrics = serde_json::Map::new();
        metrics.insert("rows".to_owned(), totals.rows.into());
        metrics.insert("cells".to_owned(), cell_entries.len().into());
        metrics.insert("umis".to_owned(), umi_entries.len().into());
        Ok(Some(OutputReport {
            name: "sc_rna_ucf".to_owned(),
            metrics: serde_json::Value::Object(metrics),
        }))
    }
}
