use bed_utils::extsort::ExternalSorterBuilder;
use indexmap::IndexMap;
use itertools::Itertools;
use noodles::sam::Header;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::{align::MultiMapR, transcript::Gene};

use super::{
    annotate::AnnotationRegion, de_dups::count_unique_umi, AlignmentAnnotator, AnnotatedAlignment,
};

#[derive(Debug, Serialize, Deserialize)]
pub struct GeneAlignment {
    pub idx: usize,
    pub umi: Option<String>,
    pub align_type: AnnotationRegion,
}

#[derive(Debug)]
pub struct Quantifier {
    annotator: AlignmentAnnotator,
    genes: IndexMap<String, Gene>,
    temp_dir: Option<PathBuf>,
    chunk_size: usize,
}

impl Quantifier {
    pub fn new(annotator: AlignmentAnnotator) -> Self {
        let genes = annotator
            .transcripts
            .iter()
            .map(|(_, t)| {
                let g = t.gene.clone();
                (g.id.clone(), g)
            })
            .collect();
        Self {
            annotator,
            genes,
            temp_dir: None,
            chunk_size: 50000000,
        }
    }

    pub fn quantify<'a, I, P>(&'a self, header: &'a Header, records: I, out_dir: P)
    where
        I: Iterator<Item = Vec<(MultiMapR, Option<MultiMapR>)>> + 'a,
        P: AsRef<std::path::Path>,
    {
        // create output files
        std::fs::create_dir_all(&out_dir).unwrap();

        let mut output_feature = seqspec::utils::create_file(
            out_dir.as_ref().join("features.tsv.gz"),
            Some(seqspec::utils::Compression::Gzip),
            Some(7),
            8,
        )
        .unwrap();
        self.genes.values().for_each(|g| {
            writeln!(output_feature, "{}\t{}", g.id, g.name).unwrap();
        });

        let mut output_barcode = seqspec::utils::create_file(
            out_dir.as_ref().join("barcodes.tsv.gz"),
            Some(seqspec::utils::Compression::Gzip),
            Some(7),
            8,
        )
        .unwrap();

        let mut output_mat = seqspec::utils::create_file(
            out_dir.as_ref().join("matrix.mtx.gz"),
            Some(seqspec::utils::Compression::Gzip),
            Some(7),
            8,
        )
        .unwrap();

        let mut n_barcodes  = 0usize; 
        let alignments = records.flat_map(|recs| {
            recs.into_iter()
                .filter_map(|(r1, r2)| self.make_gene_alignment(header, r1, r2))
        });
        let alignments_barcode = 
            sort_alignments(alignments, self.temp_dir.as_ref(), self.chunk_size)
            .chunk_by(|x| x.0.clone());
        let counts = alignments_barcode
            .into_iter()
            .enumerate()
            .flat_map(|(i, (barcode, alignments))| {
                n_barcodes += 1;
                let counts = count_unique_umi(alignments.map(|(_, a)| a));
                writeln!(output_barcode, "{}", barcode).unwrap();
                counts.into_iter().map(move |(gene, count)|
                    [gene + 1, i + 1, count]
                )
            })
            .collect::<Vec<_>>();

        writeln!(output_mat, "%%MatrixMarket matrix coordinate integer general\n%\n{} {} {}", self.genes.len(), n_barcodes, counts.len()).unwrap();
        for count in counts {
            writeln!(output_mat, "{} {} {}", count[0], count[1], count[2]).unwrap();
        }

    }

    fn make_gene_alignment(
        &self,
        header: &Header,
        rec1: MultiMapR,
        rec2: Option<MultiMapR>,
    ) -> Option<(String, GeneAlignment)> {
        let barcode = rec1.barcode().unwrap()?;
        let umi = rec1.umi().unwrap();
        let anno = if let Some(rec2) = rec2 {
            self.annotator.annotate_alignments_pe(header, rec1, rec2)
        } else {
            self.annotator.annotate_alignments_se(header, rec1)
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
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
            }
            AnnotatedAlignment::SeMapped(anno) => {
                let genes = anno.genes;
                if genes.len() != 1 {
                    return None;
                }
                let gene = genes.iter().next().unwrap();
                align_type = anno.region;
                gene_id = self.genes.get_full(&gene.id).unwrap().0;
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
