use std::path::{Path, PathBuf};

use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp, ArrayData, AxisArraysOp, ElemCollectionOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::{ExternalChunkBuilder, ExternalSorterBuilder};
use indexmap::IndexMap;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::info;
use nalgebra_sparse::CsrMatrix;
use noodles::sam::Header;
use polars::df;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};

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
        let mut qc = QcGeneQuant::default();
        let genome_index = GenomeBaseIndex::new(header);

        let mut umi_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut spliced_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut unspliced_cache: ExternalChunkBuilder<Vec<Vec<(usize, u32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;
        let mut coverage_cache: ExternalChunkBuilder<Vec<Vec<(usize, i32)>>> =
            ExternalChunkBuilder::new(tempfile::tempfile()?, 3)?;

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();

        {
            // Align reads to transcripts
            let tx_alignments = self.annotator.align(records).flat_map(|aln| {
                qc.update(aln.as_ref());
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

            // Sort alignments by cell barcodes
            let sorted_aln = sort_alignments(tx_alignments, self.chunk_size, self.temp_dir.as_ref());

            let sty = ProgressStyle::with_template(
                "{percent}%|{wide_bar:.cyan/blue}| {human_pos:>}/{human_len:} [{elapsed}<{eta}, {per_sec}]",
            )
            .unwrap();
            let bar = ProgressBar::new(sorted_aln.len() as u64).with_style(sty);

            info!("Generating gene counts ...");
            let aln_cell_group = sorted_aln.chunk_by(|x| x.barcode().unwrap().to_string());
            let chunk_size = self.num_threads * 1_000_000;
            let aln_chunks = aln_cell_group
                .into_iter()
                .scan(0, |count, (barcode, group)| {
                    barcodes.push(barcode);
                    let aln = group.collect::<Vec<_>>();
                    let n = aln.len();
                    *count += n;
                    bar.inc(n as u64);
                    Some((*count / chunk_size, aln))
                })
                .chunk_by(|(chunk_id, _)| *chunk_id);

            aln_chunks.into_iter().for_each(|(_, chunk)| {
                let ((t, s), (uns, cov)) = self.count_one_cell(&genome_index, chunk.map(|x| x.1));
                umi_cache.add(t).unwrap();
                spliced_cache.add(s).unwrap();
                unspliced_cache.add(uns).unwrap();
                coverage_cache.add(cov).unwrap();
            });

            bar.abandon();
        }

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
        Ok(qc)
    }

    fn count_one_cell(
        &self,
        index: &GenomeBaseIndex,
        alignments: impl Iterator<Item = Vec<TxAlignment>>,
    ) -> (
        (Vec<Vec<(usize, u32)>>, Vec<Vec<(usize, u32)>>),
        (Vec<Vec<(usize, u32)>>, Vec<Vec<(usize, i32)>>),
    ) {
        let alignments = alignments.collect::<Vec<_>>();
        let chunk_size = (alignments.len() / self.num_threads).max(1);
        alignments
            .par_chunks(chunk_size)
            .flat_map_iter(|chunk| {
                chunk.into_iter().map(|aln| {
                    let umi_group = GeneCounter::new(aln, &self.genes);
                    let coverage = index.count_fragments(umi_group.to_fragments());
                    let (spliced, unspliced) = umi_group.to_spliced_counts();
                    ((umi_group.to_counts(), spliced), (unspliced, coverage))
                })
            })
            .unzip()
    }
}

fn sort_alignments<I>(
    alignments: I,
    chunk_size: usize,
    temp_dir: Option<impl AsRef<Path>>,
) -> impl ExactSizeIterator<Item = TxAlignment>
where
    I: Iterator<Item = TxAlignment>,
{
    let mut sorter = ExternalSorterBuilder::new()
        .with_chunk_size(chunk_size)
        .with_compression(2);
    if let Some(tmp) = temp_dir {
        sorter = sorter.with_tmp_dir(tmp);
    }
    sorter
        .build()
        .unwrap()
        .sort_by_async(alignments, |a, b| a.barcode().cmp(&b.barcode()))
        .unwrap()
        .map(|x| x.unwrap())
}

fn make_csr<T>(data: Vec<Vec<(usize, T)>>, num_cols: usize) -> ArrayData
where
    ArrayData: From<CsrMatrix<T>>,
    ArrayData: From<anndata::data::CsrNonCanonical<T>>,
{
    let (nrows, ncols, indptr, indices, data) = to_csr_data(data, num_cols);
    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
}
