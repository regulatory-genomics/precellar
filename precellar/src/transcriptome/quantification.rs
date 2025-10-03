use anndata::{
    data::utils::{from_csr_data, to_csr_data},
    AnnData, AnnDataOp, ArrayData, AxisArraysOp,
};
use anndata_hdf5::H5;
use anyhow::Result;
use bed_utils::extsort::ExternalSorterBuilder;
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::sam::Header;
use polars::df;
use rayon::{
    iter::{IntoParallelIterator, ParallelIterator},
    slice::ParallelSlice,
};
use std::{
    collections::{BTreeMap, HashMap},
    path::PathBuf,
};

use super::AlignmentAnnotator;
use crate::{align::MultiMapR, qc::QcGeneQuant, transcriptome::AnnotatedAlignment};

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator,
    genes: IndexMap<String, String>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        let genes = annotator
            .transcripts
            .iter()
            .map(|(_, t)| (t.gene_id.clone(), t.gene_name.clone()))
            .sorted_by(|a, b| a.1.cmp(&b.1))  // sort by gene names
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

        let adata: AnnData<H5> = AnnData::new(output)?;
        let num_cols = self.genes.len();

        let mut barcodes = Vec::new();

        {
            let alignments = annotate_alignments(&self.annotator, header, records, &mut qc)
                .filter(|x| x.1.uniquely_mapped_gene().is_some());
            let alignments_cell_group =
                sort_alignments(alignments, self.temp_dir.as_ref(), self.chunk_size)
                    .chunk_by(|x| x.0.clone());
            let alignments_chunks = alignments_cell_group
                .into_iter()
                .map(|(barcode, group)| {
                    barcodes.push(barcode);
                    group.map(|x| x.1).collect::<Vec<_>>()
                })
                .chunks(512);
            let (total, exonic, intronic, reads): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
                alignments_chunks
                    .into_iter()
                    .flat_map(|chunk| {
                        let results: Vec<_> = chunk
                            .collect::<Vec<_>>()
                            .into_par_iter()
                            .map(|aln| self.count_one_cell(aln))
                            .collect();
                        results
                    })
                    .multiunzip();

            let total = make_csr(total, num_cols);
            let exonic = make_csr(exonic, num_cols);
            let intronic = make_csr(intronic, num_cols);
            adata.set_x(total)?;
            adata.layers().add("intronic", intronic)?;
            adata.layers().add("exonic", exonic)?;
        }

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
        alignments: impl IntoIterator<Item = AnnotatedAlignment>,
    ) -> (
        Vec<(usize, u32)>,
        Vec<(usize, u32)>,
        Vec<(usize, u32)>,
        Vec<AnnotatedAlignment>,
    ) {
        let mut exonic_count: BTreeMap<usize, u32> = BTreeMap::new();
        let mut intronic_count: BTreeMap<usize, u32> = BTreeMap::new();
        let mut uniq_alignments = Vec::new();

        // Count of Gene -> {UMI: (count, alignment)}
        let mut exonic_gene_umi_map = HashMap::new();
        let mut intronic_gene_umi_map = HashMap::new();
        alignments.into_iter().for_each(|alignment| {
            let gene_idx = self
                .genes
                .get_full(alignment.uniquely_mapped_gene().unwrap())
                .unwrap()
                .0;
            if let Some(umi) = alignment.umi() {
                let umi = umi.as_bytes().to_vec();
                if alignment.is_intronic() {
                    intronic_gene_umi_map
                        .entry(gene_idx)
                        .or_insert_with(HashMap::new)
                        .entry(umi)
                        .and_modify(|(count, _)| *count += 1)
                        .or_insert((1, alignment));
                } else if alignment.is_exonic() {
                    exonic_gene_umi_map
                        .entry(gene_idx)
                        .or_insert_with(HashMap::new)
                        .entry(umi)
                        .and_modify(|(count, _)| *count += 1)
                        .or_insert((1, alignment));
                }
            } else {
                // No UMI, just count the alignment
                if alignment.is_intronic() {
                    let entry = intronic_count.entry(gene_idx).or_insert(0);
                    *entry = entry.checked_add(1).unwrap();
                } else if alignment.is_exonic() {
                    let entry = exonic_count.entry(gene_idx).or_insert(0);
                    *entry = entry.checked_add(1).unwrap();
                }
                uniq_alignments.push(alignment);
            }
        });

        // Update counts for alignments with UMIs
        exonic_gene_umi_map.into_iter().for_each(|(gene, umi_map)| {
            let aln = correct_umi(umi_map);
            let entry = exonic_count.entry(gene).or_insert(0);
            *entry = entry.checked_add(aln.len() as u32).unwrap();
            uniq_alignments.extend(aln);
        });
        intronic_gene_umi_map
            .into_iter()
            .for_each(|(gene, umi_map)| {
                let aln = correct_umi(umi_map);
                let entry = intronic_count.entry(gene).or_insert(0);
                *entry = entry.checked_add(aln.len() as u32).unwrap();
                uniq_alignments.extend(aln);
            });

        let mut total = exonic_count.clone();
        for (gene, count) in intronic_count.iter() {
            let entry = total.entry(*gene).or_insert(0);
            *entry = entry.checked_add(*count).unwrap();
        }
        (
            total.into_iter().collect(),
            exonic_count.into_iter().collect(),
            intronic_count.into_iter().collect(),
            uniq_alignments,
        )
    }
}

fn annotate_alignments<'a, I>(
    annotator: &'a AlignmentAnnotator,
    header: &'a Header,
    records: I,
    qc: &'a mut QcGeneQuant,
) -> impl Iterator<Item = (String, AnnotatedAlignment)> + 'a
where
    I: Iterator<Item = Vec<(Option<MultiMapR>, Option<MultiMapR>)>> + 'a,
{
    records
        .flat_map(|recs| {
            let alignments = recs
                .par_chunks(4096)
                .flat_map_iter(|chunk| {
                    chunk.into_iter().map(|(rec1, rec2)| {
                        let barcode;
                        let anno = if rec1.is_some() && rec2.is_some() {
                            let rec1 = rec1.as_ref().unwrap();
                            let rec2 = rec2.as_ref().unwrap();
                            barcode = rec1.barcode().unwrap()?;
                            annotator.annotate_alignments_pe(header, rec1, rec2)
                        } else {
                            let rec = if rec1.is_some() {
                                rec1.as_ref().unwrap()
                            } else {
                                rec2.as_ref().unwrap()
                            };
                            barcode = rec.barcode().unwrap()?;
                            annotator.annotate_alignments_se(header, rec)
                        }?;
                        if anno.is_confidently_mapped() {
                            Some((barcode, anno))
                        } else {
                            None
                        }
                    })
                })
                .collect::<Vec<_>>();
            alignments.iter().for_each(|x| {
                let x = x.as_ref().map(|(_, a)| a);
                qc.update(x);
            });
            alignments
        })
        .flatten()
}

fn sort_alignments<I, P>(
    alignments: I,
    temp_dir: Option<P>,
    chunk_size: usize,
) -> impl ExactSizeIterator<Item = (String, AnnotatedAlignment)>
where
    I: Iterator<Item = (String, AnnotatedAlignment)>,
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

fn correct_umi<T>(umi_count: HashMap<Vec<u8>, (u64, T)>) -> Vec<T> {
    let umi_map = make_umi_map(&umi_count);
    let mut uniq_umi = HashMap::new();
    umi_count.into_iter().for_each(|(umi, (_, label))| {
        let corrected_umi = umi_map.get(umi.as_slice()).unwrap_or(&umi);
        uniq_umi.insert(corrected_umi.clone(), label);
    });
    uniq_umi.drain().map(|(_, x)| x).collect()
}

/// Returns a map from each UMI to its corrected UMI by correcting Hamming-distance-one UMIs.
fn make_umi_map<T>(umi_count: &HashMap<Vec<u8>, (u64, T)>) -> HashMap<Vec<u8>, Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for (umi, orig_count) in umi_count {
        let mut test_umi = umi.clone();

        let mut best_dest_count = orig_count.0;
        let mut best_dest_umi = umi.to_vec();

        for pos in 0..umi.len() {
            // Try each nucleotide at this position
            for test_char in nucs {
                if *test_char == umi[pos] {
                    // Skip the identitical nucleotide
                    continue;
                }
                test_umi[pos] = *test_char;

                // Test for the existence of this mutated UMI
                let test_count = *umi_count.get(&test_umi).map_or(&0, |x| &x.0);

                // If there's a 1-HD UMI w/ greater count, move to that UMI.
                // If there's a 1-HD UMI w/ equal count, move to the lexicographically larger UMI.
                if test_count > best_dest_count
                    || (test_count == best_dest_count && test_umi > best_dest_umi)
                {
                    best_dest_umi = test_umi.clone();
                    best_dest_count = test_count;
                }
            }
            // Reset this position to the unmutated sequence
            test_umi[pos] = umi[pos];
        }
        if *umi != best_dest_umi {
            corrections.insert(umi.to_vec(), best_dest_umi);
        }
    }
    corrections
}

fn make_csr(data: Vec<Vec<(usize, u32)>>, num_cols: usize) -> ArrayData {
    let (nrows, ncols, indptr, indices, data) = to_csr_data(data, num_cols);
    from_csr_data(nrows, ncols, indptr, indices, data).unwrap()
}
