use crate::barcode::{BarcodeCorrector, Whitelist};
use crate::seqspec::{Modality, RegionType, SeqSpec};
use crate::qc::{AlignQC, Metrics};

use bstr::BString;
use indicatif::ProgressStyle;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use anyhow::{bail, Result};
use noodles::{bam, sam, fastq};
use noodles::sam::alignment::{
    Record, record_buf::RecordBuf, record::data::field::tag::Tag,
};
use noodles::sam::alignment::record_buf::data::field::value::Value;
use bwa::BurrowsWheelerAligner;
use log::info;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use either::Either;

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
    base_dir: PathBuf,
    seqspec: SeqSpec,
    aligner: A,
    current_modality: Option<String>,
    mito_dna: HashSet<usize>,
    metrics: HashMap<Modality, Metrics>,
    align_qc: HashMap<Modality, Arc<Mutex<AlignQC>>>,
    barcode_correct_prob: f64,  // if the posterior probability of a correction
                                // exceeds this threshold, the barcode will be corrected.
                                // cellrange uses 0.975 for ATAC and 0.9 for multiome.
}

impl<A: Alinger> FastqProcessor<A> {
    pub fn new(seqspec: SeqSpec, aligner: A) -> Self {
        Self {
            seqspec, aligner, current_modality: None, metrics: HashMap::new(),
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

    pub fn modality(&self) -> &str {
        self.current_modality.as_ref().expect("modality not set, please call set_modality first")
    }

    pub fn add_mito_dna(&mut self, mito_dna: &str) {
        self.aligner.header().reference_sequences().get_index_of(&BString::from(mito_dna))
            .map(|x| self.mito_dna.insert(x));
    }

    pub fn with_modality(mut self, modality: &str) -> Self {
        self.current_modality = Some(modality.into());
        self
    }

    pub fn get_report(&self) -> Metrics {
        let mut metrics = self.metrics.get(self.modality()).map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(self.modality()) {
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
            self.modality().into(),
            Arc::new(Mutex::new(AlignQC {
                mito_dna: self.mito_dna.clone(),
                ..AlignQC::default()
            }))
        );
        let align_qc = self.align_qc.get(self.modality()).unwrap().clone();

        let style = ProgressStyle::with_template(
            "[{elapsed}] {bar:40.cyan/blue} {human_pos:>7}/{human_len:7} records (eta: {eta})"
        ).unwrap();
        let progress_bar = indicatif::ProgressBar::new(whitelist.total_count as u64)
            .with_style(style);
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
            progress_bar.inc(results.len() as u64);
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
            progress_bar.inc(results.len() as u64);
            Either::Left(results)
        })
    }

    pub fn gen_raw_alignments(&self) -> 
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

    pub fn gen_raw_fastq_records(&self) -> FastqRecords<impl BufRead>  {
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
                read_list.insert(&read.primer_id, read);
            });
        let data = read_list.into_iter().map(|(id, read)| {
            let is_reverse = read.is_reverse();
            let fq = fq_list.get(id).unwrap();
            let regions = if is_reverse {
                fq.subregion_range_rev()
            } else {
                fq.subregion_range()
            };
            (fq.id.clone(), is_reverse, regions, read.read_fastq(self.base_dir.clone()))
        });
        FastqRecords::new(data)
    }

    fn count_barcodes(&mut self) -> Result<Whitelist> {
        let modality = self.modality();
        let mut whitelist = self.get_whitelist()?;
        let region_with_barcode = self.seqspec.modality(self.modality()).unwrap()
            .iter_regions().find(|r|
                r.region_type == RegionType::Fastq && r.iter_regions().any(|x| x.region_type == RegionType::Barcode)
            ).unwrap();

        let sequence_read = self.seqspec.get_read_by_primer_id(modality, &region_with_barcode.id).unwrap();
        let subregions = if sequence_read.is_reverse() {
            region_with_barcode.subregion_range_rev()
        } else {
            region_with_barcode.subregion_range()
        };
        let range = subregions.into_iter().find(|x| x.0 == RegionType::Barcode).unwrap().1;
        
        sequence_read.read_fastq(&self.base_dir).records().for_each(|record| {
            let mut record = record.unwrap();
            let n = record.sequence().len();
            record = slice_fastq_record(&record, range.start, range.end.unwrap_or(n));
            if sequence_read.is_reverse() {
                record = rev_compl_fastq_record(record);
            }
            whitelist.count_barcode(std::str::from_utf8(record.sequence()).unwrap(), record.quality_scores());
        });

        self.metrics.entry(modality.to_string()).or_insert_with(Metrics::default)
            .insert("frac_q30_bases_barcode".to_string(), whitelist.frac_q30_bases());

        Ok(whitelist)
    }

    fn get_whitelist(&self) -> Result<Whitelist> {
        let regions: Vec<_> = self.seqspec.modality(self.modality()).unwrap()
            .iter_regions().filter(|r| r.region_type == RegionType::Barcode).collect();
        if regions.len() != 1 {
            bail!("Expecting exactly one barcode region, found {}", regions.len());
        }
        let region = regions[0];
        if region.sequence_type.as_str() == "onlist" {
            Ok(Whitelist::new(region.sequence_type.fetch_onlist()?))
        } else {
            Ok(Whitelist::empty())
        }
    }
}

pub struct FastqRecords<R> {
    ids: Vec<String>,
    is_reverse: Vec<bool>,
    subregions: Vec<Vec<(RegionType, crate::seqspec::Range)>>,
    readers: Vec<fastq::Reader<R>>,
}

pub type Barcode = fastq::Record;
pub type UMI = fastq::Record;

impl<R: BufRead> FastqRecords<R> {
    fn new<I>(iter: I) -> Self
    where
        I: Iterator<Item = (String, bool, Vec<(RegionType, crate::seqspec::Range)>, fastq::Reader<R>)>,
    {
        let mut ids = Vec::new();
        let mut is_reverse = Vec::new();
        let mut subregions = Vec::new();
        let mut readers = Vec::new();
        iter.for_each(|(f, s, sr, r)| {
            ids.push(f);
            is_reverse.push(s);
            subregions.push(sr.into_iter().filter(|x|
                x.0 == RegionType::Barcode || x.0 == RegionType::CDNA || x.0 == RegionType::GDNA
            ).collect());
            readers.push(r);
        });
        Self { ids, is_reverse, subregions, readers }
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
                        if self.is_reverse[i] {
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
        let records = self.readers.iter_mut().enumerate().map(|(i, reader)| {
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
                let is_reverse = self.is_reverse[i];
                match region_type {
                    RegionType::Barcode => if is_reverse {
                        barcode = Some(rev_compl_fastq_record(fq));
                    } else {
                        barcode = Some(fq);
                    },
                    RegionType::CDNA | RegionType::GDNA => if is_reverse {
                        read1 = Some(fq);
                    } else {
                        read2 = Some(fq);
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
        let mut processor = FastqProcessor::new(spec, aligner)
            .with_modality("atac");

        processor.gen_barcoded_alignments().take(6).for_each(|x| {
            ()
        });

        println!("{}", processor.get_report());
    }
}