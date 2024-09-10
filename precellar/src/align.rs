use crate::barcode::{BarcodeCorrector, Whitelist};
use crate::seqspec::{open_file_for_read, RegionType, SeqSpec};

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use anyhow::{bail, Result};
use multi_reader::MultiReader;
use noodles::{sam, fastq};
use noodles::sam::alignment::{
    Record, record_buf::RecordBuf, record::data::field::tag::Tag,
};
use noodles::sam::alignment::record_buf::data::field::value::Value;
use bwa::BurrowsWheelerAligner;
use log::info;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

pub trait Alinger {
    fn chunk_size(&self) -> usize;

    fn header(&self) -> sam::Header;

    fn align_reads(&self, records: &mut [fastq::Record]) -> impl ExactSizeIterator<Item = sam::Record>;

    fn align_read_pairs(&self, records: &mut [(fastq::Record, fastq::Record)]) ->
        impl ExactSizeIterator<Item = (sam::Record, sam::Record)>;
}

impl Alinger for BurrowsWheelerAligner {
    fn chunk_size(&self) -> usize {
        self.chunk_size()
    }

    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(&self, records: &mut [fastq::Record]) -> impl ExactSizeIterator<Item = sam::Record> {
        self.align_reads(records)
    }

    fn align_read_pairs(&self, records: &mut [(fastq::Record, fastq::Record)]) ->
        impl ExactSizeIterator<Item = (sam::Record, sam::Record)> {
        self.align_read_pairs(records)
    }
}

pub struct FastqProcessor<A> {
    seqspec: SeqSpec,
    aligner: A,
    current_modality: Option<String>,
}

impl<A: Alinger> FastqProcessor<A> {
    pub fn new(seqspec: SeqSpec, aligner: A) -> Self {
        Self { seqspec, aligner, current_modality: None }
    }

    pub fn modality(&self) -> &str {
        self.current_modality.as_ref().expect("modality not set, please call set_modality first")
    }

    pub fn set_modality(mut self, modality: &str) -> Self {
        self.current_modality = Some(modality.into());
        self
    }

    pub fn gen_barcoded_alignments(&self) -> Box<dyn Iterator<Item = RecordBuf> + '_> {
        info!("Counting barcodes...");
        let whitelist = self.count_barcodes().unwrap();

        let fq_records = self.gen_raw_fastq_records();
        let is_paired = fq_records.is_paired_end();
        let fq_chunks = fq_records.chunk(self.aligner.chunk_size());
        let corrector = BarcodeCorrector::default();
        let header = self.aligner.header();

        info!("Aligning reads...");
        if is_paired {
            Box::new(fq_chunks.flat_map(move |data| {
                let (barcodes, mut reads): (Vec<_>, Vec<_>) = data.into_iter()
                    .map(|(barcode, (read1, read2))| (barcode, (read1.unwrap(), read2.unwrap()))).unzip();
                let alignments: Vec<_> = self.aligner.align_read_pairs(&mut reads).collect();
                barcodes.into_par_iter().zip(alignments).flat_map(|(barcode, (ali1, ali2))| {
                    let corrected_barcode = corrector.correct(
                        whitelist.get_barcode_counts(),
                        std::str::from_utf8(barcode.sequence()).unwrap(),
                        barcode.quality_scores()
                    ).ok();
                    let ali1_ = add_cell_barcode(
                        &header,
                        &ali1,
                        std::str::from_utf8(barcode.sequence()).unwrap(),
                        barcode.quality_scores(),
                        corrected_barcode.as_ref().map(|x| x.as_str())
                    ).unwrap();
                    let ali2_ = add_cell_barcode(
                        &header,
                        &ali2,
                        std::str::from_utf8(barcode.sequence()).unwrap(),
                        barcode.quality_scores(),
                        corrected_barcode.as_ref().map(|x| x.as_str())
                    ).unwrap();
                    [ali1_, ali2_]
                }).collect::<Vec<_>>()
            }))
        } else {
            Box::new(fq_chunks.flat_map(move |data| {
                let (barcodes, mut reads): (Vec<_>, Vec<_>) = data.into_iter()
                    .map(|(barcode, (read1, _))| (barcode, read1.unwrap())).unzip();
                let alignments: Vec<_> = self.aligner.align_reads(&mut reads).collect();
                barcodes.into_par_iter().zip(alignments).map(|(barcode, alignment)| {
                    let corrected_barcode = corrector.correct(
                        whitelist.get_barcode_counts(),
                        std::str::from_utf8(barcode.sequence()).unwrap(),
                        barcode.quality_scores()
                    ).ok();
                    add_cell_barcode(
                        &header,
                        &alignment,
                        std::str::from_utf8(barcode.sequence()).unwrap(),
                        barcode.quality_scores(),
                        corrected_barcode.as_ref().map(|x| x.as_str())
                    ).unwrap()
                }).collect::<Vec<_>>()
            }))
        }
    }

    pub fn gen_raw_alignments(&self) -> Box<dyn Iterator<Item = sam::Record> + '_> {
        let fq_records = self.gen_raw_fastq_records();
        let is_paired = fq_records.is_paired_end();
        let fq_chunks = fq_records.chunk(self.aligner.chunk_size());
        if is_paired {
            Box::new(fq_chunks.flat_map(|data| {
                let mut reads: Vec<_> = data.into_iter().map(|(_, (read1, read2))| (read1.unwrap(), read2.unwrap())).collect();
                self.aligner.align_read_pairs(&mut reads).flat_map(|(r1, r2)| [r1, r2]).collect::<Vec<_>>()
            }))
        } else {
            Box::new(fq_chunks.flat_map(|data| {
                let mut reads: Vec<_> = data.into_iter().map(|(_, (read1, _))| read1.unwrap()).collect();
                self.aligner.align_reads(reads.as_mut()).collect::<Vec<_>>()
            }))
        }
    }

    pub fn gen_raw_fastq_records(&self) -> FastqRecords<BufReader<MultiReader<Box<(dyn std::io::Read + 'static)>, std::vec::IntoIter<Box<(dyn std::io::Read + 'static)>>>>> {
        let modality = self.modality();
        let fq_list: HashMap<_, _> = self.seqspec.modality(modality).unwrap()
            .iter_regions()
            .filter(|region| region.region_type == RegionType::Fastq &&
                region.iter_regions().any(|r|
                    r.region_type == RegionType::GDNA ||
                    r.region_type == RegionType::CDNA ||
                    r.region_type == RegionType::Barcode ||
                    r.region_type == RegionType::UMI
                )
            ).map(|fq| (&fq.id, fq)).collect();
        let mut read_list = HashMap::new();
        self.seqspec.sequence_spec.get(modality).unwrap().iter()
            .for_each(|read| if fq_list.contains_key(&read.primer_id) {
                read_list.entry(&read.primer_id).or_insert_with(Vec::new).push(read);
            });
        let data = read_list.into_iter().map(|(id, reads)| {
            let strand = reads[0].strand;
            let fq = fq_list.get(id).unwrap();
            let regions = fq.subregion_range();
            let readers: Vec<_> = reads.iter().map(|read| open_file_for_read(&read.id)).collect();
            let readers = MultiReader::new(readers.into_iter());
            let readers = BufReader::new(readers);
            (fq.id.clone(), strand, regions, readers)
        });
        FastqRecords::new(data)
    }

    fn count_barcodes(&self) -> Result<Whitelist> {
        let modality = self.modality();
        let mut whitelist = self.get_whitelist()?.unwrap();
        let region = self.seqspec.modality(self.modality()).unwrap()
            .iter_regions().find(|r|
                r.region_type == RegionType::Fastq && r.iter_regions().any(|x| x.region_type == RegionType::Barcode)
            ).unwrap();

        let readers = self.seqspec.get_read_by_primer_id(modality, &region.id)
            .into_iter().map(|read| open_file_for_read(&read.id));
        let range = region.subregion_range().find(|x| x.0 == RegionType::Barcode).unwrap().1;
        fastq::Reader::new(BufReader::new(MultiReader::new(readers)))
            .records().for_each(|record| {
                let record = record.unwrap();
                let n = record.sequence().len();
                let slice = range.start..range.end.unwrap_or(n);
                let barcode = record.sequence().get(slice.clone()).unwrap();
                let barcode_qual = record.quality_scores().get(slice).unwrap();
                whitelist.count_barcode(std::str::from_utf8(barcode).unwrap(), barcode_qual);
            });

        Ok(whitelist)
    }

    fn get_whitelist(&self) -> Result<Option<Whitelist>> {
        let regions: Vec<_> = self.seqspec.modality(self.modality()).unwrap()
            .iter_regions().filter(|r| r.region_type == RegionType::Barcode).collect();
        if regions.len() != 1 {
            bail!("Expecting exactly one barcode region, found {}", regions.len());
        }
        let region = regions[0];
        if region.sequence_type.as_str() == "onlist" {
            Ok(Some(Whitelist::new(region.sequence_type.fetch_onlist()?)))
        } else {
            Ok(None)
        }
    }
}

pub struct FastqRecords<R> {
    ids: Vec<String>,
    strandness: Vec<bool>,
    subregions: Vec<Vec<(RegionType, crate::seqspec::Range)>>,
    records: Vec<fastq::Reader<R>>,
}

pub type Barcode = fastq::Record;
pub type UMI = fastq::Record;

impl<R: BufRead> FastqRecords<R> {
    fn new<I, It>(iter: I) -> Self
    where
        I: Iterator<Item = (String, bool, It, R)>,
        It: Iterator<Item = (RegionType, crate::seqspec::Range)>,
    {
        let mut ids = Vec::new();
        let mut strandness = Vec::new();
        let mut subregions = Vec::new();
        let mut records = Vec::new();
        iter.for_each(|(f, s, sr, r)| {
            ids.push(f);
            strandness.push(s);
            subregions.push(sr.into_iter().filter(|x|
                x.0 == RegionType::Barcode || x.0 == RegionType::CDNA || x.0 == RegionType::GDNA
            ).collect());
            records.push(fastq::Reader::new(r));
        });
        Self { ids, strandness, subregions, records }
    }

    fn chunk(self, chunk_size: usize) -> FastqRecordChunk<R> {
        FastqRecordChunk { fq: self, chunk_size }
    }

    fn is_paired_end(&self) -> bool {
        let mut read1 = false;
        let mut read2 = false;
        self.subregions.iter().enumerate().for_each(|(i, sr)| {
            sr.iter().for_each(|(region_type, _)| {
                match region_type {
                    RegionType::CDNA | RegionType::GDNA => {
                        if self.strandness[i] {
                            read1 = true;
                        } else {
                            read2 = true;
                        }
                    },
                    _ => (),
                }
            });
        });
        read1 && read2
    }
}

impl<R: BufRead> Iterator for FastqRecords<R> {
    type Item = (Barcode, (Option<fastq::Record>, Option<fastq::Record>));

    fn next(&mut self) -> Option<Self::Item> {
        let mut id_without_record = Vec::new();
        let records = self.records.iter_mut().enumerate().map(|(i, reader)| {
            let mut record = fastq::Record::default();
            let s = reader.read_record(&mut record).expect("error reading fastq record");
            if s == 0 {
                id_without_record.push(self.ids[i].as_str());
                None
            } else {
                Some(record)
            }
        }).collect::<Vec<_>>();
        if id_without_record.len() == records.len() {
            return None;
        } else if id_without_record.len() > 0 {
            panic!("Missing records in these files: {}", id_without_record.join(","));
        }

        let mut barcode = None;
        let mut read1 = None;
        let mut read2 = None;
        records.iter().enumerate().for_each(|(i, r)| {
            let record = r.as_ref().unwrap();
            self.subregions[i].iter().for_each(|(region_type, range)| {
                let fq = slice_fastq_record(record, range.start, range.end.unwrap_or(record.sequence().len()));
                match region_type {
                    RegionType::Barcode => barcode = Some(fq),
                    RegionType::CDNA | RegionType::GDNA => {
                        if self.strandness[i] {
                            read1 = Some(fq);
                        } else {
                            read2 = Some(fq);
                        }
                    },
                    _ => (),
                }
            });
        });
        let barcode = barcode.expect("barcode should be present");
        Some((barcode, (read1, read2)))
    }
}

pub struct FastqRecordChunk<R> {
    fq: FastqRecords<R>,
    chunk_size: usize,
}

impl<R: BufRead> Iterator for FastqRecordChunk<R> {
    type Item = Vec<(Barcode, (Option<fastq::Record>, Option<fastq::Record>))>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk = Vec::new();
        let mut accumulated_length = 0;

        for record in self.fq.by_ref() {
            accumulated_length += record.1.0.as_ref().map_or(0, |x| x.sequence().len()) +
                record.1.1.as_ref().map_or(0, |x| x.sequence().len());
            chunk.push(record);
            if accumulated_length >= self.chunk_size {
                break;
            }
        }

        if chunk.is_empty() {
            None
        } else {
            Some(chunk)
        }
    }
}


pub fn add_cell_barcode<R: Record>(
    header: &sam::Header,
    record: &R,
    ori_barcode: &str,
    ori_qual: &[u8],
    correct_barcode: Option<&str>,
) -> std::io::Result<RecordBuf> {
    let mut record_buf = RecordBuf::try_from_alignment_record(header, record)?;
    let data = record_buf.data_mut();
    data.insert(Tag::CELL_BARCODE_SEQUENCE, Value::String(ori_barcode.into()));
    data.insert(Tag::CELL_BARCODE_QUALITY_SCORES, Value::String(ori_qual.into()));
    if let Some(barcode) = correct_barcode {
        data.insert(Tag::CELL_BARCODE_ID, Value::String(barcode.into()));
    }
    Ok(record_buf)
}

fn slice_fastq_record(record: &fastq::Record, start: usize, end: usize) -> fastq::Record {
    fastq::Record::new(
        record.definition().clone(),
        &record.sequence()[start..end],
        record.quality_scores().get(start..end).unwrap(),
    )
}

#[cfg(test)]
mod tests {
    use bwa::{AlignerOpts, FMIndex, PairedEndStats};

    use super::*;

    #[test]
    fn test_seqspec_io() {
        let spec = SeqSpec::from_path("tests/data/spec.yaml").unwrap();
        let aligner = BurrowsWheelerAligner::new(
            FMIndex::read("tests/data/hg38").unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default()
        );
        let processor = FastqProcessor::new(spec, aligner)
            .set_modality("atac");

        processor.gen_barcoded_alignments().take(6).for_each(|x| {
            println!("{:?}", x);
        });

        /*
        let whitelist = processor
            .set_modality("atac")
            .count_barcodes().unwrap();
        println!("{:?}", whitelist);
        */
    }
}