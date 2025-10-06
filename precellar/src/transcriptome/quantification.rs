use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp, ArrayData, AxisArraysOp, ElemCollectionOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::ExternalSorterBuilder;
use indexmap::IndexMap;
use itertools::Itertools;
use nalgebra_sparse::CsrMatrix;
use noodles::sam::Header;
use polars::df;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::path::PathBuf;

use super::TxAligner;
use crate::{
    align::MultiMapR,
    genome::GenomeBaseIndex,
    qc::QcGeneQuant,
    transcriptome::{CorrectUMI, TxAlignment},
};

#[derive(Debug)]
pub struct Quantifier {
    annotator: TxAligner,
    genes: IndexMap<String, String>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
}

impl Quantifier {
    pub fn new(annotator: TxAligner) -> Self {
        let genes = annotator
            .transcripts()
            .map(|t| (t.gene_id.clone(), t.gene_name.clone()))
            .sorted_by(|a, b| a.1.cmp(&b.1)) // sort by gene names
            .collect();
        Self {
            annotator,
            genes,
            temp_dir: None,
            chunk_size: 50000000,
        }
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
        let mut num_unique_umi = 0;
        let mut num_spliced = 0;
        let mut num_unspliced = 0;
        let genome_index = GenomeBaseIndex::new(header);

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();

        {
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

            let alignments_cell_group =
                sort_alignments(tx_alignments, self.temp_dir.as_ref(), self.chunk_size)
                    .chunk_by(|x| x.barcode().unwrap().to_string());
            let alignments_chunks = alignments_cell_group
                .into_iter()
                .map(|(barcode, group)| {
                    barcodes.push(barcode);
                    group.collect::<Vec<_>>()
                })
                .chunks(512);
            let (total, spliced, unspliced, coverage): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
                alignments_chunks
                    .into_iter()
                    .flat_map(|chunk| {
                        let results: Vec<_> = chunk
                            .collect::<Vec<_>>()
                            .into_par_iter()
                            .map(|aln| self.count_one_cell(&genome_index, aln))
                            .collect();
                        results
                    })
                    .multiunzip();

            num_unique_umi += total
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            num_spliced += spliced
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            num_unspliced += unspliced
                .iter()
                .flat_map(|x| x.iter().map(|(_, c)| *c as u64))
                .sum::<u64>();
            let coverage = make_csr(coverage, genome_index.len());
            adata.set_x(make_csr(total, num_cols))?;
            adata.layers().add("spliced", make_csr(spliced, num_cols))?;
            adata
                .layers()
                .add("unspliced", make_csr(unspliced, num_cols))?;
            adata.obsm().add("fragment_single", coverage)?;
            adata
                .uns()
                .add("reference_sequences", genome_index.get_chrom_sizes())?;
        }

        adata.set_obs_names(barcodes.into())?;
        adata.set_var_names(self.genes.keys().map(|g| g.clone()).collect())?;
        adata.set_var(df!(
            "gene_name" => self.genes.values().map(|g| g.clone()).collect::<Vec<_>>()
        )?)?;

        adata.close()?;
        qc.num_unique_umi = num_unique_umi;
        qc.num_spliced = num_spliced;
        qc.num_unspliced = num_unspliced;
        Ok(qc)
    }

    fn count_one_cell(
        &self,
        index: &GenomeBaseIndex,
        alignments: impl IntoIterator<Item = TxAlignment>,
    ) -> (
        Vec<(usize, u32)>,
        Vec<(usize, u32)>,
        Vec<(usize, u32)>,
        Vec<(usize, i32)>,
    ) {
        let umi_group = alignments.into_iter().correct_umi(&self.genes);

        let coverage = index.count_fragments(umi_group.to_fragments());
        let (spliced, unspliced) = umi_group.to_spliced_counts();
        (umi_group.to_counts(), spliced, unspliced, coverage)
    }
}

fn sort_alignments<I, P>(
    alignments: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = TxAlignment>
where
    I: Iterator<Item = TxAlignment>,
    P: AsRef<std::path::Path>,
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
        .sort_by(alignments, |a, b| a.barcode().cmp(&b.barcode()))
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
