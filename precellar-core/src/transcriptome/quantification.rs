use std::{collections::HashMap, ops::Range, path::PathBuf, sync::mpsc};

use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp, ArrayData, AxisArraysOp, ElemCollectionOp,
};
use anndata_hdf5::H5;
use anyhow::{ensure, Result};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use bitcode::{Decode, Encode};
use indexmap::IndexMap;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::info;
use nalgebra_sparse::CsrMatrix;
use noodles_sam::Header;
use polars::df;
use rayon::prelude::*;

use super::TxAligner;
use crate::{
    align::MultiMapR,
    genome::GenomeBaseIndex,
    qc::QcGeneQuant,
    transcriptome::{GeneCounter, TxAlignment},
};

pub struct Quantifier {
    annotator: TxAligner,
    genes: IndexMap<String, String>,
    chunk_size: usize,
    pub temp_dir: Option<PathBuf>,
    pub num_threads: usize,
}

impl Quantifier {
    pub fn new(annotator: TxAligner) -> Result<Self> {
        let genes = annotator
            .transcripts()
            .map(|t| (t.gene_id.clone(), t.gene_name.clone()))
            .sorted_by(|a, b| a.1.cmp(&b.1)) // sort by gene names
            .collect();
        Ok(Self {
            annotator,
            genes,
            chunk_size: 10000000,
            num_threads: 16,
            temp_dir: None,
        })
    }

    pub fn quantify<'a, I, P>(
        &'a self,
        header: &'a Header,
        records: I,
        output: P,
    ) -> Result<QcGeneQuant>
    where
        I: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
        P: AsRef<std::path::Path>,
    {
        // The stream is consumed lazily inside the external sorter. Its end
        // marks the start of the `quantify` stage (sort collection, counting
        // and writing), which `end_quantify` closes after the output is written.
        let records = crate::timing::trace_tail(records);
        let mut qc = QcGeneQuant::default();
        let genome_index = GenomeBaseIndex::new(header);

        let mut umi_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut spliced_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut unspliced_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut ambiguous_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut coverage_cache: ExternalChunkBuilder<Vec<Vec<(usize, i32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();

        {
            // Align reads to transcripts
            let tx_alignments = self.annotator.align(records, &mut qc).flat_map(|aln| {
                if aln
                    .as_ref()
                    .and_then(|x| x.uniquely_mapped_gene())
                    .is_none()
                {
                    None
                } else {
                    aln
                }
            });

            let mut barcodes_seen = BarcodeTable::default();
            let sortable_aln = tx_alignments.map(|aln| barcodes_seen.record(aln));
            let mut runs = sort_runs(sortable_aln, self.chunk_size, self.temp_dir.as_deref())?;

            let sty = ProgressStyle::with_template(
                "{percent}%|{wide_bar:.cyan/blue}| {human_pos:>}/{human_len:} [{elapsed}<{eta}, {per_sec}]",
            )
            .unwrap();
            let bar = ProgressBar::new(barcodes_seen.entries.iter().map(|x| x.count as u64).sum())
                .with_style(sty);

            info!("Generating gene counts ...");
            let entries = barcodes_seen.entries;
            let mut start = 0;
            while start < entries.len() {
                let alignments = {
                    let _span = crate::timing::stage(crate::timing::Stage::CountCollect);
                    let (end, count) = batch_end(&entries, start, self.num_threads * 1_000_000);
                    let range = start..end;
                    let alignments = collect_batch(&mut runs, &entries, range.clone())?;
                    barcodes.extend(entries[range].iter().map(|entry| entry.name.clone()));
                    bar.inc(count as u64);
                    start = end;
                    alignments
                };
                let counts = {
                    let _span = crate::timing::stage(crate::timing::Stage::CellCount);
                    self.count_cells(&genome_index, alignments)
                };
                let _span = crate::timing::stage(crate::timing::Stage::WriteCache);
                umi_cache.add(counts.umi).unwrap();
                spliced_cache.add(counts.spliced).unwrap();
                unspliced_cache.add(counts.unspliced).unwrap();
                ambiguous_cache.add(counts.ambiguous).unwrap();
                coverage_cache.add(counts.coverage).unwrap();
            }

            ensure!(
                runs.iter().all(|run| run.pending.is_none()),
                "unconsumed sort records"
            );
            bar.abandon();
        }

        let _span = crate::timing::stage(crate::timing::Stage::WriteDisk);

        let umi_count = umi_cache.finish()?.map(|x| {
            let x = x.unwrap();
            qc.num_unique_umi += x
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            make_csr(x, num_cols)
        });
        adata.set_x_from_iter(umi_count)?;

        let spliced_count = spliced_cache.finish()?.map(|x| {
            let x = x.unwrap();
            qc.num_spliced += x
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            make_csr(x, num_cols)
        });
        adata.layers().add_iter("spliced", spliced_count)?;

        let unspliced_count = unspliced_cache.finish()?.map(|x| {
            let x = x.unwrap();
            qc.num_unspliced += x
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            make_csr(x, num_cols)
        });
        adata.layers().add_iter("unspliced", unspliced_count)?;

        let ambiguous_count = ambiguous_cache.finish()?.map(|x| {
            let x = x.unwrap();
            qc.num_ambiguous += x
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            make_csr(x, num_cols)
        });
        adata.layers().add_iter("ambiguous", ambiguous_count)?;

        let coverage = coverage_cache.finish()?.map(|x| {
            let x = x.unwrap();
            make_csr(x, genome_index.len())
        });
        adata.obsm().add_iter("fragment_single", coverage)?;

        adata
            .uns()
            .add("reference_sequences", genome_index.get_chrom_sizes())?;
        adata.set_obs_names(barcodes.into())?;
        adata.set_var_names(self.genes.keys().map(|g| g.clone()).collect())?;
        adata.set_var(df!(
            "gene_name" => self.genes.values().map(|g| g.clone()).collect::<Vec<_>>()
        )?)?;

        adata.close()?;
        crate::timing::end_quantify();
        Ok(qc)
    }

    /// Counts a batch of cells, one entry per cell in each returned vector.
    fn count_cells(&self, index: &GenomeBaseIndex, alignments: Vec<CellAlignments>) -> CellCounts {
        let per_cell: Vec<_> = crate::threads::run(|| {
            alignments
                .into_par_iter()
                .map(|parts| {
                    let umi_group =
                        GeneCounter::from_alignments(parts.iter().flatten(), &self.genes);
                    let coverage = index.count_fragments(umi_group.to_fragments());
                    (
                        umi_group.to_counts(),
                        umi_group.to_spliced_counts(),
                        coverage,
                    )
                })
                .collect()
        });

        let mut counts = CellCounts::default();
        for (umi, spliced, coverage) in per_cell {
            counts.umi.push(umi);
            counts.spliced.push(spliced.spliced);
            counts.unspliced.push(spliced.unspliced);
            counts.ambiguous.push(spliced.ambiguous);
            counts.coverage.push(coverage);
        }
        counts
    }
}

/// Counts for a batch of cells, each vector holding one entry per cell.
#[derive(Default)]
struct CellCounts {
    umi: Vec<Vec<(usize, u32)>>,
    spliced: Vec<Vec<(usize, u32)>>,
    unspliced: Vec<Vec<(usize, u32)>>,
    ambiguous: Vec<Vec<(usize, u32)>>,
    coverage: Vec<Vec<(usize, i32)>>,
}

#[derive(Default)]
struct BarcodeTable {
    ids: HashMap<String, u32>,
    entries: Vec<BarcodeEntry>,
}

struct BarcodeEntry {
    name: String,
    count: usize,
}

impl BarcodeTable {
    fn record(&mut self, alignment: TxAlignment) -> SortableAlignment {
        let barcode = alignment
            .barcode()
            .expect("quantification input must have a barcode");
        let barcode_id = match self.ids.get(barcode) {
            Some(&id) => id,
            None => {
                let id = u32::try_from(self.entries.len()).expect("too many barcode IDs");
                self.ids.insert(barcode.to_owned(), id);
                self.entries.push(BarcodeEntry {
                    name: barcode.to_owned(),
                    count: 0,
                });
                id
            }
        };
        self.entries[barcode_id as usize].count += 1;
        SortableAlignment {
            barcode_id,
            alignment,
        }
    }
}

// Both the table and the runs use the same ID order. A cell larger than the
// target occupies its own batch instead of being split across two batches.
fn batch_end(entries: &[BarcodeEntry], start: usize, target: usize) -> (usize, usize) {
    let mut end = start + 1;
    let mut count = entries[start].count;
    while end < entries.len() && entries[end].count <= target.saturating_sub(count) {
        count += entries[end].count;
        end += 1;
    }
    (end, count)
}

#[derive(Encode, Decode)]
struct SortableAlignment {
    barcode_id: u32,
    alignment: TxAlignment,
}

struct SortedRun {
    chunk: ExternalChunk<SortableAlignment>,
    pending: Option<SortableAlignment>,
}

impl SortedRun {
    fn read_batch(&mut self, range: &Range<usize>) -> Result<Vec<(usize, Vec<TxAlignment>)>> {
        let mut groups: Vec<(usize, Vec<TxAlignment>)> = Vec::new();
        while let Some(record) = self
            .pending
            .take()
            .map(Ok)
            .or_else(|| self.chunk.next())
            .transpose()?
        {
            let id = record.barcode_id as usize;
            if id >= range.end {
                self.pending = Some(record);
                break;
            }
            ensure!(id >= range.start, "sort record precedes the current batch");
            match groups.last_mut() {
                Some((last_id, records)) if *last_id == id => records.push(record.alignment),
                _ => groups.push((id, vec![record.alignment])),
            }
        }
        Ok(groups)
    }
}

// One cell may have a group from each run. Move only these group descriptors;
// cell counting traverses their records directly instead of copying them.
type CellAlignments = Vec<Vec<TxAlignment>>;

fn collect_batch(
    runs: &mut [SortedRun],
    entries: &[BarcodeEntry],
    range: Range<usize>,
) -> Result<Vec<CellAlignments>> {
    let groups: Vec<_> = crate::threads::run(|| {
        runs.par_iter_mut()
            .map(|run| run.read_batch(&range))
            .collect::<Result<_>>()
    })?;
    let mut cells: Vec<CellAlignments> = (0..range.len()).map(|_| Vec::new()).collect();
    for (id, records) in groups.into_iter().flatten() {
        cells[id - range.start].push(records);
    }
    for (cell, entry) in cells.iter().zip(&entries[range]) {
        ensure!(
            cell.iter().map(Vec::len).sum::<usize>() == entry.count,
            "unexpected number of alignments for {}",
            entry.name
        );
    }
    Ok(cells)
}

fn sort_runs(
    alignments: impl Iterator<Item = SortableAlignment>,
    chunk_size: usize,
    temp_dir: Option<&std::path::Path>,
) -> Result<Vec<SortedRun>> {
    // Spill in the background while the caller feeds the pipeline. Bound the
    // number of in-flight buffers, each of which holds up to one full sort run.
    let max_in_flight = crate::threads::current_num_threads().min(4);
    let pool = crate::threads::build_pool(max_in_flight)?;
    let temp_dir = temp_dir.map_or_else(tempfile::tempdir, tempfile::tempdir_in)?;
    let (tx, rx) = mpsc::channel();
    let mut runs = Vec::new();
    let mut in_flight = 0;
    let mut buffer = Vec::new();

    let submit = |mut records: Vec<SortableAlignment>| {
        let tx = tx.clone();
        let path = temp_dir.path().to_owned();
        pool.spawn(move || {
            let result = (|| -> Result<_> {
                records.sort_unstable_by_key(|record| record.barcode_id);
                let mut builder = ExternalChunkBuilder::new(tempfile::tempfile_in(path)?, 2)?;
                for record in records {
                    builder.add(record)?;
                }
                Ok(SortedRun {
                    chunk: builder.finish()?,
                    pending: None,
                })
            })();
            let _ = tx.send(result);
        });
    };

    for record in alignments {
        buffer.push(record);
        if buffer.len() >= chunk_size {
            submit(std::mem::take(&mut buffer));
            in_flight += 1;
            if in_flight == max_in_flight {
                runs.push(rx.recv()??);
                in_flight -= 1;
            }
        }
    }
    if !buffer.is_empty() {
        submit(buffer);
    }
    drop(tx);
    for result in rx {
        runs.push(result?);
    }
    Ok(runs)
}

fn make_csr<T>(data: Vec<Vec<(usize, T)>>, num_cols: usize) -> ArrayData
where
    ArrayData: From<CsrMatrix<T>>,
    ArrayData: From<anndata::data::CsrNonCanonical<T>>,
{
    let (nrows, ncols, indptr, indices, data) = to_csr_data(data, num_cols);
    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcriptome::align::{AlignType, SplicedRecord};

    fn read(barcode: &str, umi: &str) -> TxAlignment {
        TxAlignment::SeAligned {
            read: SplicedRecord::for_test("chr1", vec![(0, 10)]),
            alignment: vec![("tx1".into(), AlignType::Exonic)],
            gene_id: "gene1".into(),
            barcode: barcode.into(),
            umi: Some(umi.into()),
        }
    }

    #[test]
    fn collect_across_runs_and_batch_boundaries() -> Result<()> {
        let mut table = BarcodeTable::default();
        // The first cell spans three runs, and the second spans two. The
        // first run has to retain its next cell when batch one stops.
        let input = [
            ("TTTT", "AAAA"),
            ("AAAA", "GGGG"),
            ("TTTT", "AAAA"),
            ("CCCC", "TTTT"),
            ("TTTT", "AAAT"),
            ("AAAA", "GGGG"),
        ];
        let mut runs = sort_runs(
            input
                .into_iter()
                .map(|(bc, umi)| table.record(read(bc, umi))),
            2,
            None,
        )?;
        assert_eq!(runs.len(), 3);
        let entries = table.entries;
        assert_eq!(
            entries.iter().map(|x| x.name.as_str()).collect::<Vec<_>>(),
            ["TTTT", "AAAA", "CCCC"]
        );
        assert_eq!(batch_end(&entries, 0, 3), (1, 3));
        assert_eq!(batch_end(&entries, 0, 2), (1, 3));
        assert_eq!(batch_end(&entries, 1, 3), (3, 3));

        let first = collect_batch(&mut runs, &entries, 0..1)?;
        assert_eq!(first[0].len(), 3);
        let genes = IndexMap::from([("gene1".into(), "gene1".into())]);
        let count = GeneCounter::from_alignments(first[0].iter().flatten(), &genes);
        assert_eq!(count.to_counts(), [(0, 1)]);
        let rest = collect_batch(&mut runs, &entries, 1..3)?;
        assert_eq!(rest.iter().map(|x| x.len()).collect::<Vec<_>>(), [2, 1]);
        assert!(runs.iter().all(|run| run.pending.is_none()));
        Ok(())
    }

    #[test]
    fn no_runs_for_empty_input() -> Result<()> {
        assert!(sort_runs(std::iter::empty(), 2, None)?.is_empty());
        Ok(())
    }
}
