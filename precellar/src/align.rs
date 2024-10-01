use crate::barcode::{BarcodeCorrector, Whitelist};
use seqspec::{Assay, Modality, Read, Region, RegionType, SequenceType};
use crate::qc::{AlignQC, Metrics};

use bstr::BString;
use kdam::{tqdm, BarExt};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use smallvec::SmallVec;
use anyhow::{bail, Result};
use noodles::{bam, sam, fastq};
use noodles::sam::alignment::{
    Record, record_buf::RecordBuf, record::data::field::tag::Tag,
};
use noodles::sam::alignment::record_buf::data::field::value::Value;
use bwa_mem2::BurrowsWheelerAligner;
use log::info;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use either::Either;

pub trait Alinger {
    fn chunk_size(&self) -> usize;

    fn header(&self) -> sam::Header;

    fn align_reads(&mut self, records: &mut [fastq::Record]) -> impl ExactSizeIterator<Item = sam::Record>;

    fn align_read_pairs(&mut self, records: &mut [(fastq::Record, fastq::Record)]) ->
        impl ExactSizeIterator<Item = (sam::Record, sam::Record)>;
}

impl Alinger for BurrowsWheelerAligner {
    fn chunk_size(&self) -> usize {
        self.chunk_size()
    }

    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(&mut self, records: &mut [fastq::Record]) -> impl ExactSizeIterator<Item = sam::Record> {
        self.align_reads(records)
    }

    fn align_read_pairs(&mut self, records: &mut [(fastq::Record, fastq::Record)]) ->
        impl ExactSizeIterator<Item = (sam::Record, sam::Record)> {
        self.align_read_pairs(records)
    }
}

pub struct FastqProcessor<A> {
    base_dir: PathBuf,
    assay: Assay,
    aligner: A,
    current_modality: Option<Modality>,
    mito_dna: HashSet<usize>,
    metrics: HashMap<Modality, Metrics>,
    align_qc: HashMap<Modality, Arc<Mutex<AlignQC>>>,
    barcode_correct_prob: f64,  // if the posterior probability of a correction
                                // exceeds this threshold, the barcode will be corrected.
                                // cellrange uses 0.975 for ATAC and 0.9 for multiome.
}

impl<A: Alinger> FastqProcessor<A> {
    pub fn new(assay: Assay, aligner: A) -> Self {
        Self {
            assay, aligner, current_modality: None, metrics: HashMap::new(),
            align_qc: HashMap::new(), mito_dna: HashSet::new(), base_dir: PathBuf::from("./"),
            barcode_correct_prob: 0.975,
        }
    }

    pub fn with_base_dir<P: AsRef<Path>>(mut self, base_dir: P) -> Self {
        self.base_dir = base_dir.as_ref().to_path_buf();
        self
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        self.barcode_correct_prob = prob;
        self
    }

    pub fn modality(&self) -> Modality {
        self.current_modality.expect("modality not set, please call set_modality first")
    }

    pub fn add_mito_dna(&mut self, mito_dna: &str) {
        self.aligner.header().reference_sequences().get_index_of(&BString::from(mito_dna))
            .map(|x| self.mito_dna.insert(x));
    }

    pub fn with_modality(mut self, modality: Modality) -> Self {
        self.current_modality = Some(modality.into());
        self
    }

    pub fn get_report(&self) -> Metrics {
        let mut metrics = self.metrics.get(&self.modality()).map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(&self.modality()) {
            align_qc.lock().unwrap().report(&mut metrics);
        }
        metrics
    }

    pub fn gen_barcoded_alignments(&mut self) ->
        impl Iterator<Item = Either<Vec<RecordBuf>, Vec<(RecordBuf, RecordBuf)>>> + '_
    {
        info!("Counting barcodes...");
        let whitelist = self.count_barcodes().unwrap();
        info!("Found {} barcodes. {:.2}% of them have an exact match in whitelist", whitelist.total_count, whitelist.frac_exact_match() * 100.0);

        info!("Aligning reads...");
        let fq_records = self.gen_raw_fastq_records();
        let is_paired = fq_records.is_paired_end();
        let corrector = BarcodeCorrector::default().with_bc_confidence_threshold(self.barcode_correct_prob);
        let header = self.aligner.header();
        self.align_qc.insert(
            self.modality(),
            Arc::new(Mutex::new(AlignQC {
                mito_dna: self.mito_dna.clone(),
                ..AlignQC::default()
            }))
        );
        let align_qc = self.align_qc.get(&self.modality()).unwrap().clone();

        let mut progress_bar = tqdm!(total = whitelist.total_count);
        fq_records.chunk(self.aligner.chunk_size()).map(move |data| if is_paired {
            let (barcodes, mut reads): (Vec<_>, Vec<_>) = data.into_iter()
                .map(|(barcode, (read1, read2))| (barcode, (read1.unwrap(), read2.unwrap()))).unzip();
            let alignments: Vec<_> = self.aligner.align_read_pairs(&mut reads).collect();
            let results = barcodes.into_par_iter().zip(alignments).map(|(barcode, (ali1, ali2))| {
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
                {
                    let mut align_qc_lock = align_qc.lock().unwrap();
                    align_qc_lock.update(&ali1_, &header);
                    align_qc_lock.update(&ali2_, &header);
                }
                (ali1_, ali2_)
            }).collect::<Vec<_>>();
            progress_bar.update(results.len()).unwrap();
            Either::Right(results)
        } else {
            let (barcodes, mut reads): (Vec<_>, Vec<_>) = data.into_iter()
                .map(|(barcode, (read1, _))| (barcode, read1.unwrap())).unzip();
            let alignments: Vec<_> = self.aligner.align_reads(&mut reads).collect();
            let results = barcodes.into_par_iter().zip(alignments).map(|(barcode, alignment)| {
                let corrected_barcode = corrector.correct(
                    whitelist.get_barcode_counts(),
                    std::str::from_utf8(barcode.sequence()).unwrap(),
                    barcode.quality_scores()
                ).ok();
                let ali = add_cell_barcode(
                    &header,
                    &alignment,
                    std::str::from_utf8(barcode.sequence()).unwrap(),
                    
                    barcode.quality_scores(),
                    corrected_barcode.as_ref().map(|x| x.as_str())
                ).unwrap();
                { align_qc.lock().unwrap().update(&ali, &header); }
                ali
            }).collect::<Vec<_>>();
            progress_bar.update(results.len()).unwrap();
            Either::Left(results)
        })
    }

    pub fn gen_raw_alignments(&mut self) -> 
        impl Iterator<Item = Either<Vec<sam::Record>, Vec<(sam::Record, sam::Record)>>> + '_
    {
        let fq_records = self.gen_raw_fastq_records();
        let is_paired = fq_records.is_paired_end();
        fq_records.chunk(self.aligner.chunk_size()).map(move |data| if is_paired {
            let mut reads: Vec<_> = data.into_iter().map(|(_, (read1, read2))|
                (read1.unwrap(), read2.unwrap())).collect();
            let alignments = self.aligner.align_read_pairs(&mut reads).collect();
            Either::Right(alignments)
        } else {
            let mut reads: Vec<_> = data.into_iter().map(|(_, (read1, _))|
                read1.unwrap()).collect();
            let alignments = self.aligner.align_reads(reads.as_mut()).collect();
            Either::Left(alignments)
        })
    }

    pub fn gen_raw_fastq_records(&self) -> FastqRecords<impl BufRead> {
        let modality = self.modality();
        let data = self.assay.get_index_of(modality)
            .map(|(read, regions)| (read, regions, crate::io::read_fastq(read, self.base_dir.clone())));
        FastqRecords::new(data)
    }

    fn count_barcodes(&mut self) -> Result<Whitelist> {
        let modality = self.modality();
        let mut whitelist = self.get_whitelist()?;

        let (read, index) = self.assay.get_index_of(modality).into_iter()
            .find(|(_, index)| index.into_iter().any(|x| x.0.region_type.is_barcode()))
            .expect("No barcode region found");
        let range = index.into_iter().find(|x| x.0.region_type.is_barcode()).unwrap().1;

        crate::io::read_fastq(read, &self.base_dir).records().for_each(|record| {
            let mut record = record.unwrap();
            record = slice_fastq_record(&record, range.start as usize, range.end as usize);
            if read.is_reverse() {
                record = rev_compl_fastq_record(record);
            }
            whitelist.count_barcode(std::str::from_utf8(record.sequence()).unwrap(), record.quality_scores());
        });

        self.metrics.entry(modality).or_insert_with(Metrics::default)
            .insert("frac_q30_bases_barcode".to_string(), whitelist.frac_q30_bases());

        Ok(whitelist)
    }

    fn get_whitelist(&self) -> Result<Whitelist> {
        let regions: Vec<_> = self.assay.get_region_by_modality(self.modality()).unwrap()
            .iter_regions().filter(|r| r.region_type.is_barcode()).collect();
        if regions.len() != 1 {
            bail!("Expecting exactly one barcode region, found {}", regions.len());
        }
        let region = regions[0];
        if region.sequence_type == SequenceType::Onlist {
            Ok(Whitelist::new(crate::io::read_onlist(region.onlist.as_ref().unwrap())?))
        } else {
            Ok(Whitelist::empty())
        }
    }
}

pub struct FastqRecords<R>(Vec<FastqRecord<R>>);

struct FastqRecord<R> {
    id: String,
    is_reverse: bool,
    subregion: Vec<(RegionType, Range<u32>)>,
    reader: fastq::Reader<R>,
    min_len: usize,
    max_len: usize,
}

pub type Barcode = fastq::Record;
pub type UMI = fastq::Record;

impl<R: BufRead> FastqRecords<R> {
    fn new<'a, I>(iter: I) -> Self
    where
        I: Iterator<Item = (&'a Read, Vec<(&'a Region, Range<u32>)>, fastq::Reader<R>)>,
    {
        let records = iter.map(|(read, regions, reader)|
            FastqRecord {
                id: read.read_id.to_string(),
                is_reverse: read.is_reverse(),
                subregion: regions.into_iter().filter_map(|x| {
                    let region_type = x.0.region_type;
                    if region_type.is_barcode() || region_type.is_dna() {
                        Some((region_type, x.1))
                    } else {
                        None
                    }
                }).collect(),
                reader,
                min_len: read.min_len as usize,
                max_len: read.max_len as usize,
            }
        ).collect();
        Self(records)
    }

    fn chunk(self, chunk_size: usize) -> FastqRecordChunk<R> {
        FastqRecordChunk { fq: self, chunk_size }
    }

    fn is_paired_end(&self) -> bool {
        let mut read1 = false;
        let mut read2 = false;
        self.0.iter().for_each(|x| {
            x.subregion.iter().for_each(|(region_type, _)| if region_type.is_dna() {
                if x.is_reverse {
                    read1 = true;
                } else {
                    read2 = true;
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
        let records: SmallVec<[_; 4]> = self.0.iter_mut().map(|x| {
            let mut record = fastq::Record::default();
            if x.reader.read_record(&mut record).expect("error reading fastq record") == 0 {
                id_without_record.push(x.id.as_str());
                None
            } else {
                let n = record.sequence().len();
                if n < x.min_len || n > x.max_len {
                    panic!("Read length ({}) out of range: {}-{}", n, x.min_len, x.max_len);
                }
                Some(record)
            }
        }).collect();
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
            self.0[i].subregion.iter().for_each(|(region_type, range)| {
                let fq = slice_fastq_record(record, range.start as usize, range.end as usize);
                let is_reverse = self.0[i].is_reverse;
                if region_type.is_barcode() {
                    if is_reverse {
                        barcode = Some(rev_compl_fastq_record(fq));
                    } else {
                        barcode = Some(fq);
                    }
                } else if region_type.is_dna() {
                    if is_reverse {
                        read1 = Some(fq);
                    } else {
                        read2 = Some(fq);
                    }
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


fn add_cell_barcode<R: Record>(
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

fn rev_compl_fastq_record(mut record: fastq::Record) -> fastq::Record {
    *record.quality_scores_mut() = record.quality_scores().iter().rev().copied().collect();
    *record.sequence_mut() = record.sequence().iter().rev().map(|&x| match x {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => x,
    }).collect();
    record
}

pub struct NameCollatedRecords<'a, R> {
    records: bam::io::reader::Records<'a, R>,
    prev_record: Option<(BString, bam::Record)>,
    checker: HashSet<BString>,
}

impl<'a, R: std::io::Read> NameCollatedRecords<'a, R> {
    pub fn new(records: bam::io::reader::Records<'a, R>) -> Self {
        Self {
            records,
            prev_record: None,
            checker: HashSet::new(),
        }
    }

    fn check(&mut self, name: &BString) {
        assert!(!self.checker.contains(name), "bam file must be name collated or name sorted");
        self.checker.insert(name.to_owned());
    }
}

impl<'a, R: std::io::Read> Iterator for NameCollatedRecords<'a, R> {
    type Item = (bam::Record, bam::Record);

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.records.next()?.unwrap();
        let name = record.name().unwrap().to_owned();
        if let Some((prev_name, prev_record)) = self.prev_record.take() {
            if name == prev_name {
                Some((prev_record, record))
            } else {
                panic!("Expecting paired end reads with the same name, found {} and {}", prev_name, name);
            }
        } else {
            self.check(&name);
            self.prev_record = Some((name, record));
            self.next()
        }
    }
}


#[cfg(test)]
mod tests {
    use bwa_mem2::{AlignerOpts, FMIndex, PairedEndStats};

    use super::*;

    #[test]
    fn test_seqspec_io() {
        let spec = Assay::from_path("tests/data/spec.yaml").unwrap();
        let aligner = BurrowsWheelerAligner::new(
            FMIndex::read("tests/data/hg38").unwrap(),
            AlignerOpts::default(),
            PairedEndStats::default()
        );
        let mut processor = FastqProcessor::new(spec, aligner)
            .with_modality(Modality::ATAC);

        processor.gen_barcoded_alignments().take(6).for_each(|x| {
            println!("{:?}", x);
        });

        println!("{}", processor.get_report());
    }
}