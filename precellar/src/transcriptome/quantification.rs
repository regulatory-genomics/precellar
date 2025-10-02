use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::ExternalSorterBuilder;
use bincode::{Decode, Encode};
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::sam::Header;
use polars::df;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{collections::HashSet, path::PathBuf};

use crate::{align::MultiMapR, qc::QcGeneQuant};

use super::{
    annotate::AnnotationRegion, de_dups::count_unique_umi, AlignmentAnnotator, AnnotatedAlignment,
};

#[derive(Debug, Encode, Decode)]
pub struct GeneAlignment {
    pub idx: usize,
    pub umi: Option<String>,
    pub align_type: AnnotationRegion,
}

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator,
    genes: IndexMap<String, String>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
    mito_genes: HashSet<usize>,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        let genes = annotator
            .transcripts
            .iter()
            .map(|(_, t)| (t.gene_id.clone(), t.gene_name.clone()))
            .collect();
        Self {
            annotator,
            genes,
            temp_dir: None,
            chunk_size: 50000000,
            mito_genes: HashSet::new(),
        }
    }

    pub fn add_mito_dna(&mut self, mito_chr: &str) {
        let iter = self.annotator.transcripts.iter().flat_map(|(_, t)| {
            if t.chrom == mito_chr {
                Some(self.genes.get_full(&t.gene_id).unwrap().0)
            } else {
                None
            }
        });
        self.mito_genes.extend(iter);
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

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();
        let mut exon_count: Vec<u64> = Vec::new();
        let mut intron_count: Vec<u64> = Vec::new();
        let mut mito_count: Vec<u64> = Vec::new();

        {
            let tx_alignments = records.flat_map(|recs| {
                qc.total_reads += recs.len() as u64;
                recs.into_par_iter()
                    .filter_map(|(r1, r2)| self.make_gene_alignment(header, r1, r2))
                    .collect::<Vec<_>>()
            });
            let tx_alignments_chunks =
                sort_alignments(tx_alignments, self.temp_dir.as_ref(), self.chunk_size)
                    .chunk_by(|x| x.0.clone());
            let tx_alignments_chunks = tx_alignments_chunks
                .into_iter()
                .map(|(barcode, group)| {
                    barcodes.push(barcode);
                    group.map(|x| x.1).collect::<Vec<_>>()
                })
                .chunks(500);
            let counts = tx_alignments_chunks.into_iter().map(|chunk| {
                let results: Vec<_> = chunk
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .map(|alignments| count_unique_umi(alignments, &self.mito_genes))
                    .collect();
                results.iter().for_each(|r| {
                    qc.total_umi += r.total_umi;
                    qc.unique_umi += r.unique_umi;
                    exon_count.push(r.uniq_exon);
                    intron_count.push(r.uniq_intron);
                    mito_count.push(r.uniq_mito);
                });
                let (nrows, ncols, indptr, indices, data) = to_csr_data(
                    results
                        .into_iter()
                        .map(|r| r.into_counts().collect()),
                    num_cols,
                );
                from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
            });

            adata.set_x_from_iter(counts)?;
        }

        adata.set_obs_names(barcodes.into())?;
        adata.set_obs(df!(
            "exon_count" => exon_count,
            "intron_count" => intron_count,
            "mitochondrial_count" => mito_count
        )?)?;

        adata.set_var_names(self.genes.keys().map(|g| g.clone()).collect())?;
        adata.set_var(df!(
            "gene_name" => self.genes.values().map(|g| g.clone()).collect::<Vec<_>>()
        )?)?;

        adata.close()?;
        Ok(qc)
    }

    fn make_gene_alignment(
        &self,
        header: &Header,
        rec1: Option<MultiMapR>,
        rec2: Option<MultiMapR>,
    ) -> Option<(String, GeneAlignment)> {
        let barcode;
        let umi;
        let anno = if rec1.is_some() && rec2.is_some() {
            let rec1 = rec1.unwrap();
            let rec2 = rec2.unwrap();
            barcode = rec1.barcode().unwrap()?;
            umi = rec1.umi().unwrap();
            self.annotator.annotate_alignments_pe(header, rec1, rec2)
        } else {
            let rec = rec1.or(rec2).unwrap();
            barcode = rec.barcode().unwrap()?;
            umi = rec.umi().unwrap();
            self.annotator.annotate_alignments_se(header, rec)
        }?;

        let gene_id;
        let align_type;

        match anno {
            AnnotatedAlignment::PeMapped(a1, a2, anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = match (a1.region, a2.region) {
                    (AnnotationRegion::Intronic, _) => AnnotationRegion::Intronic,
                    (_, AnnotationRegion::Intronic) => AnnotationRegion::Intronic,
                    _ => AnnotationRegion::Exonic,
                };
                gene_id = self.genes.get_full(gene).unwrap().0;
            }
            AnnotatedAlignment::SeMapped(anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = anno.region;
                gene_id = self.genes.get_full(gene).unwrap().0;
            }
        }

        let alignment = GeneAlignment {
            idx: gene_id,
            umi,
            align_type,
        };
        Some((barcode, alignment))
    }
}

fn sort_alignments<I, P>(
    alignments: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = (String, GeneAlignment)>
where
    I: Iterator<Item = (String, GeneAlignment)>,
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
        .sort_by(alignments, |a, b| a.0.cmp(&b.0))
        .unwrap()
        .map(|x| x.unwrap())
}
