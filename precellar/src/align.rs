use crate::barcode::{BarcodeCorrector, OligoFrequncy, Whitelist};
use crate::qc::{AlignQC, Metrics};
use anyhow::{bail, Result};
use bstr::BString;
use bwa_mem2::BurrowsWheelerAligner;
use either::Either;
use indexmap::IndexMap;
use kdam::{tqdm, BarExt};
use log::info;
use noodles::sam::alignment::record_buf::data::field::value::Value;
use noodles::sam::alignment::{record::data::field::tag::Tag, record_buf::RecordBuf, Record};
use noodles::{bam, fastq, sam};
use rayon::iter::ParallelIterator;
use seqspec::{Assay, Modality, Read, RegionId, RegionIndex, RegionType};
use smallvec::SmallVec;
use star_aligner::StarAligner;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::ops::Range;
use std::sync::{Arc, Mutex};

pub trait AsIterator {
    type Item;
    type AsIter<'a>: Iterator<Item = &'a Self::Item> where Self: 'a;

    fn as_iter(&self) -> Self::AsIter<'_>;
}

impl AsIterator for RecordBuf {
    type Item = RecordBuf;
    type AsIter<'a> = std::iter::Once<&'a RecordBuf>;

    fn as_iter(&self) -> Self::AsIter<'_> {
        std::iter::once(&self)
    }
}

impl AsIterator for Vec<RecordBuf> {
    type Item = RecordBuf;
    type AsIter<'a> = std::slice::Iter<'a, RecordBuf>;

    fn as_iter(&self) -> Self::AsIter<'_> {
        self.iter()
    }
}

pub trait Aligner {
    type AlignOutput: AsIterator<Item = RecordBuf>;

    fn chunk_size(&self) -> usize;

    fn header(&self) -> sam::Header;

    fn align_reads(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<Self::AlignOutput>;

    fn align_read_pairs(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<(Self::AlignOutput, Self::AlignOutput)>;
}

pub struct DummyAligner;

impl Aligner for DummyAligner {
    type AlignOutput = RecordBuf;

    fn chunk_size(&self) -> usize {
        0
    }

    fn header(&self) -> sam::Header {
        sam::Header::default()
    }

    fn align_reads(&mut self, _: &sam::Header, _: Vec<AnnotatedRecord>) -> Vec<Self::AlignOutput> {
        Vec::new()
    }

    fn align_read_pairs(&mut self, _: &sam::Header, _: Vec<AnnotatedRecord>) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
        Vec::new()
    }
}

impl Aligner for BurrowsWheelerAligner {
    type AlignOutput = RecordBuf;

    fn chunk_size(&self) -> usize {
        self.chunk_size()
    }

    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<Self::AlignOutput> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records
            .into_iter()
            .map(|rec| ((rec.barcode.unwrap(), rec.umi), rec.read1.unwrap()))
            .unzip();

        // TODO: add UMI
        self.align_reads(reads.as_mut_slice()).enumerate().map(|(i, alignment)| {
            let (bc, umi) = info.get(i).unwrap();
            add_cell_barcode(
                header,
                &alignment,
                bc.raw.sequence(),
                bc.raw.quality_scores(),
                bc.corrected.as_deref(),
            )
            .unwrap()
        }).collect()
    }

    fn align_read_pairs(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records 
            .into_iter()
            .map(|rec| {
                (
                    (rec.barcode.unwrap(), rec.umi),
                    (rec.read1.unwrap(), rec.read2.unwrap()),
                )
            })
            .unzip();
        self.align_read_pairs(&mut reads).enumerate().map(|(i, (ali1, ali2))| {
            let (bc, umi) = info.get(i).unwrap();
            let ali1_ = add_cell_barcode(
                &header,
                &ali1,
                bc.raw.sequence(),
                bc.raw.quality_scores(),
                bc.corrected.as_deref(),
            )
            .unwrap();
            let ali2_ = add_cell_barcode(
                &header,
                &ali2,
                bc.raw.sequence(),
                bc.raw.quality_scores(),
                bc.corrected.as_deref(),
            )
            .unwrap();
            (ali1_, ali2_)
        }).collect()
    }
}

impl Aligner for StarAligner {
    type AlignOutput = Vec<RecordBuf>;

    fn chunk_size(&self) -> usize {
        0
    }

    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn align_reads(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<Self::AlignOutput> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records
            .into_iter()
            .map(|rec| ((rec.barcode.unwrap(), rec.umi), rec.read1.unwrap()))
            .unzip();

        // TODO: StarAligner can expose a method to align a single read instead of a batch,
        // so that barcode and UMI processing can be done in parallel.
        StarAligner::align_reads(self, reads.as_mut_slice())
            .collect::<Vec<_>>().into_iter().enumerate().map(|(i, alignment)| {
                let (bc, umi) = info.get(i).unwrap();
                alignment.into_iter().map(|x| 
                    add_cell_barcode(
                        header,
                        &x,
                        bc.raw.sequence(),
                        bc.raw.quality_scores(),
                        bc.corrected.as_deref(),
                    )
                    .unwrap()
                ).collect()
            }).collect()
    }

    fn align_read_pairs(&mut self, header: &sam::Header, records: Vec<AnnotatedRecord>) -> Vec<(Self::AlignOutput, Self::AlignOutput)> {
        let (info, mut reads): (Vec<_>, Vec<_>) = records 
            .into_iter()
            .map(|rec| {
                (
                    (rec.barcode.unwrap(), rec.umi),
                    (rec.read1.unwrap(), rec.read2.unwrap()),
                )
            })
            .unzip();
        StarAligner::align_read_pairs(self, &mut reads)
            .collect::<Vec<_>>().into_iter().enumerate().map(|(i, (ali1, ali2))| {
                let (bc, umi) = info.get(i).unwrap();
                let ali1_ = ali1.into_iter().map(|x| add_cell_barcode(
                    &header,
                    &x,
                    bc.raw.sequence(),
                    bc.raw.quality_scores(),
                    bc.corrected.as_deref(),
                )
                .unwrap()).collect();
                let ali2_ = ali2.into_iter().map(|x| add_cell_barcode(
                    &header,
                    &x,
                    bc.raw.sequence(),
                    bc.raw.quality_scores(),
                    bc.corrected.as_deref(),
                )
                .unwrap()).collect();
                (ali1_, ali2_)
            }).collect()
    }
}

pub struct FastqProcessor<A> {
    assay: Assay,
    aligner: A,
    current_modality: Option<Modality>,
    mito_dna: HashSet<usize>,
    metrics: HashMap<Modality, Metrics>,
    align_qc: HashMap<Modality, Arc<Mutex<AlignQC>>>,
    barcode_correct_prob: f64, // if the posterior probability of a correction
    // exceeds this threshold, the barcode will be corrected.
    // cellrange uses 0.975 for ATAC and 0.9 for multiome.
    mismatch_in_barcode: usize, // The number of mismatches allowed in barcode
}

impl<A: Aligner> FastqProcessor<A> {
    pub fn new(assay: Assay, aligner: A) -> Self {
        Self {
            assay,
            aligner,
            current_modality: None,
            metrics: HashMap::new(),
            align_qc: HashMap::new(),
            mito_dna: HashSet::new(),
            barcode_correct_prob: 0.975,
            mismatch_in_barcode: 1,
        }
    }

    pub fn with_barcode_correct_prob(mut self, prob: f64) -> Self {
        self.barcode_correct_prob = prob;
        self
    }

    pub fn modality(&self) -> Modality {
        self.current_modality
            .expect("modality not set, please call set_modality first")
    }

    pub fn add_mito_dna(&mut self, mito_dna: &str) {
        self.aligner
            .header()
            .reference_sequences()
            .get_index_of(&BString::from(mito_dna))
            .map(|x| self.mito_dna.insert(x));
    }

    pub fn with_modality(mut self, modality: Modality) -> Self {
        self.current_modality = Some(modality);
        self
    }

    pub fn get_report(&self) -> Metrics {
        let mut metrics = self
            .metrics
            .get(&self.modality())
            .map_or(Metrics::default(), |x| x.clone());
        if let Some(align_qc) = self.align_qc.get(&self.modality()) {
            align_qc.lock().unwrap().report(&mut metrics);
        }
        metrics
    }

    pub fn gen_barcoded_alignments(
        &mut self,
    ) -> impl Iterator<Item = Either<Vec<A::AlignOutput>, Vec<(A::AlignOutput, A::AlignOutput)>>> + '_ {
        let fq_reader = self.gen_barcoded_fastq(true);
        let is_paired = fq_reader.is_paired_end();

        info!("Aligning reads...");
        let header = self.aligner.header();
        self.align_qc.insert(
            self.modality(),
            Arc::new(Mutex::new(AlignQC {
                mito_dna: self.mito_dna.clone(),
                ..AlignQC::default()
            })),
        );
        let align_qc = self.align_qc.get(&self.modality()).unwrap().clone();

        let mut progress_bar = tqdm!(total = fq_reader.total_reads.unwrap_or(0));
        let fq_reader = VectorChunk::new(fq_reader, self.aligner.chunk_size());
        fq_reader.map(move |data| {
            let mut align_qc_lock = align_qc.lock().unwrap();
            if is_paired {
                let results: Vec<_> = self.aligner.align_read_pairs(&header, data);
                results.iter().for_each(|(ali1, ali2)| {
                    ali1.as_iter().for_each(|x| align_qc_lock.update(x, &header));
                    ali2.as_iter().for_each(|x| align_qc_lock.update(x, &header));
                });
                progress_bar.update(results.len()).unwrap();
                Either::Right(results)
            } else {
                let results: Vec<_> = self.aligner.align_reads(&header, data);
                results.iter().for_each(|ali| {
                    ali.as_iter().for_each(|x| align_qc_lock.update(x, &header));
                });
                progress_bar.update(results.len()).unwrap();
                Either::Left(results)
            }
        })
    }

    pub fn gen_barcoded_fastq(&mut self, correct_barcode: bool) -> AnnotatedFastqReader {
        let modality = self.modality();

        let whitelists = if correct_barcode {
            info!("Counting barcodes...");
            // let whitelist = self.count_barcodes().unwrap();
            let whitelists = self.count_barcodes().unwrap();
            for (id, whitelist) in whitelists.iter() {
                info!(
                    "Found {} barcodes. {:.2}% of them have an exact match in whitelist {}",
                    whitelist.total_count,
                    whitelist.frac_exact_match() * 100.0,
                    id
                );
            }
            whitelists
        } else {
            IndexMap::new()
        };

        let corrector = BarcodeCorrector::default()
            .with_max_missmatch(self.mismatch_in_barcode)
            .with_bc_confidence_threshold(self.barcode_correct_prob);

        let mut fq_reader: AnnotatedFastqReader = self
            .assay
            .get_index_by_modality(modality)
            .filter_map(|(read, index)| {
                let annotator = FastqAnnotator::new(read, index, &whitelists, corrector.clone())?;
                Some((annotator, read.open().unwrap()))
            })
            .collect();
        if !whitelists.is_empty() {
            fq_reader.total_reads = Some(whitelists[0].total_count);
        }
        fq_reader
    }

    fn count_barcodes(&mut self) -> Result<IndexMap<RegionId, Whitelist>> {
        let modality = self.modality();
        let mut whitelists = self.get_whitelists()?;

        fn count(
            read: &Read,
            barcode_region_index: RegionIndex,
            whitelist: &mut Whitelist,
        ) -> Result<()> {
            let range = &barcode_region_index
                .index
                .iter()
                .find(|x| x.1.is_barcode())
                .unwrap()
                .2;
            read.open().unwrap().records().for_each(|record| {
                let mut record = record.unwrap();
                record = slice_fastq_record(&record, range.start as usize, range.end as usize);
                if read.is_reverse() {
                    record = rev_compl_fastq_record(record);
                }
                whitelist.count_barcode(record.sequence(), record.quality_scores());
            });
            Ok(())
        }

        for (i, (read, barcode_region_index)) in self
            .assay
            .get_index_by_modality(modality)
            .filter(|(_, region_index)| region_index.index.iter().any(|x| x.1.is_barcode()))
            .enumerate()
        {
            count(read, barcode_region_index, &mut whitelists[i])?;
        }

        self.metrics.entry(modality).or_default().insert(
            "frac_q30_bases_barcode".to_string(),
            whitelists.values().map(|x| x.frac_q30_bases()).sum::<f64>() / whitelists.len() as f64,
        );
        Ok(whitelists)
    }

    fn get_whitelists(&self) -> Result<IndexMap<RegionId, Whitelist>> {
        let regions = self
            .assay
            .library_spec
            .get_modality(&self.modality())
            .unwrap()
            .read()
            .map_err(|_| anyhow::anyhow!("Cannot obtain lock"))?;
        regions
            .subregions
            .iter()
            .filter_map(|r| {
                let r = r.read().unwrap();
                if r.region_type.is_barcode() {
                    if let Some(onlist) = r.onlist.as_ref() {
                        let list = onlist
                            .read()
                            .map(|list| (r.region_id.to_string(), Whitelist::new(list)));
                        Some(list)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect()
    }
}

pub struct AnnotatedFastqReader {
    buffer: fastq::Record,
    total_reads: Option<usize>,
    inner: Vec<(FastqAnnotator, fastq::Reader<Box<dyn BufRead>>)>,
}

impl AnnotatedFastqReader {
    pub fn get_all_barcodes(&self) -> Vec<(&str, usize)> {
        self.inner
            .iter()
            .flat_map(|(annotator, _)| {
                annotator
                    .subregions
                    .iter()
                    .filter(|(_, region_type, _)| region_type.is_barcode())
                    .map(|(id, _, r)| (id.as_str(), r.len()))
            })
            .collect()
    }

    pub fn get_all_umi(&self) -> Vec<(&str, usize)> {
        self.inner
            .iter()
            .flat_map(|(annotator, _)| {
                annotator
                    .subregions
                    .iter()
                    .filter(|(_, region_type, _)| region_type.is_umi())
                    .map(|(id, _, r)| (id.as_str(), r.len()))
            })
            .collect()
    }

    pub fn is_paired_end(&self) -> bool {
        let mut has_read1 = false;
        let mut has_read2 = false;
        self.inner.iter().for_each(|x| {
            x.0.subregions.iter().for_each(|(_, region_type, _)| {
                if region_type.is_target() {
                    if x.0.is_reverse {
                        has_read1 = true;
                    } else {
                        has_read2 = true;
                    }
                }
            });
        });
        has_read1 && has_read2
    }
}

impl FromIterator<(FastqAnnotator, fastq::Reader<Box<dyn BufRead>>)> for AnnotatedFastqReader {
    fn from_iter<T: IntoIterator<Item = (FastqAnnotator, fastq::Reader<Box<dyn BufRead>>)>>(
        iter: T,
    ) -> Self {
        Self {
            buffer: fastq::Record::default(),
            total_reads: None,
            inner: iter.into_iter().collect(),
        }
    }
}

impl Iterator for AnnotatedFastqReader {
    type Item = AnnotatedRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut missing = None;
        let records: SmallVec<[_; 4]> = self
            .inner
            .iter_mut()
            .flat_map(|(annotator, reader)| {
                if reader
                    .read_record(&mut self.buffer)
                    .expect("error reading fastq record")
                    == 0
                {
                    missing = Some(annotator.id.as_str());
                    None
                } else {
                    Some(annotator.annotate(&self.buffer).unwrap())
                }
            })
            .collect();

        if records.is_empty() {
            None
        } else if missing.is_some() {
            panic!("Missing records in this file: {}", missing.unwrap());
        } else {
            Some(
                records
                    .into_iter()
                    .reduce(|mut this, other| {
                        this.join(other);
                        this
                    })
                    .unwrap(),
            )
        }
    }
}

/// A FastqAnnotator that splits the reads into subregions, e.g., barcode, UMI, and
/// return annotated reads.
struct FastqAnnotator {
    whitelists: IndexMap<String, OligoFrequncy>,
    corrector: BarcodeCorrector,
    id: String,
    is_reverse: bool,
    subregions: Vec<(String, RegionType, Range<u32>)>,
    min_len: usize,
    max_len: usize,
}

impl FastqAnnotator {
    pub fn new(
        read: &Read,
        index: RegionIndex,
        whitelists: &IndexMap<String, Whitelist>,
        corrector: BarcodeCorrector,
    ) -> Option<Self> {
        let subregions: Vec<_> = index
            .index
            .into_iter()
            .filter(|x| x.1.is_barcode() || x.1.is_umi() || x.1.is_target()) // only barcode and target regions
            .collect();
        if subregions.is_empty() {
            None
        } else {
            let whitelists = subregions
                .iter()
                .flat_map(|(id, _, _)| {
                    let v = whitelists.get(id)?;
                    Some((id.clone(), v.get_barcode_counts().clone()))
                })
                .collect();
            let anno = Self {
                whitelists,
                corrector,
                id: read.read_id.clone(),
                is_reverse: read.is_reverse(),
                subregions,
                min_len: read.min_len as usize,
                max_len: read.max_len as usize,
            };
            Some(anno)
        }
    }

    fn annotate(&self, record: &fastq::Record) -> Result<AnnotatedRecord> {
        let n = record.sequence().len();
        if n < self.min_len || n > self.max_len {
            bail!(
                "Read length ({}) out of range: {}-{}",
                n,
                self.min_len,
                self.max_len
            );
        }

        let mut barcode: Option<Barcode> = None;
        let mut umi = None;
        let mut read1 = None;
        let mut read2 = None;
        self.subregions.iter().for_each(|(id, region_type, range)| {
            let mut fq = slice_fastq_record(record, range.start as usize, range.end as usize);
            if self.is_reverse && (region_type.is_barcode() || region_type.is_umi()) {
                fq = rev_compl_fastq_record(fq);
            }
            if region_type.is_umi() {
                umi = Some(fq);
            } else if region_type.is_barcode() {
                let corrected =
                    self.whitelists
                        .get(id)
                        .map_or(Some(fq.sequence().to_vec()), |counts| {
                            self.corrector
                                .correct(counts, fq.sequence(), fq.quality_scores())
                                .ok()
                                .map(|x| x.to_vec())
                        });
                if let Some(bc) = &mut barcode {
                    bc.extend(&Barcode { raw: fq, corrected });
                } else {
                    barcode = Some(Barcode { raw: fq, corrected });
                }
            } else if region_type.is_target() {
                if self.is_reverse {
                    if let Some(s) = &mut read2 {
                        extend_fastq_record(s, &fq);
                    } else {
                        read2 = Some(fq);
                    }
                } else if let Some(s) = &mut read1 {
                    extend_fastq_record(s, &fq);
                } else {
                    read1 = Some(fq);
                }
            }
        });
        Ok(AnnotatedRecord {
            barcode,
            umi,
            read1,
            read2,
        })
    }
}

pub struct Barcode {
    pub raw: fastq::Record,
    pub corrected: Option<Vec<u8>>,
}

impl Barcode {
    pub fn extend(&mut self, other: &Self) {
        extend_fastq_record(&mut self.raw, &other.raw);
        if let Some(c2) = &other.corrected {
            if let Some(c1) = &mut self.corrected {
                c1.extend_from_slice(c2);
            }
        } else {
            self.corrected = None;
        }
    }
}

pub type UMI = fastq::Record;

/// An annotated fastq record with barcode, UMI, and sequence.
pub struct AnnotatedRecord {
    pub barcode: Option<Barcode>,
    pub umi: Option<UMI>,
    pub read1: Option<fastq::Record>,
    pub read2: Option<fastq::Record>,
}

impl AnnotatedRecord {
    pub fn len(&self) -> usize {
        self.read1.as_ref().map_or(0, |x| x.sequence().len())
            + self.read2.as_ref().map_or(0, |x| x.sequence().len())
    }
    pub fn is_empty(&self) -> bool {
        self.read1.is_none() && self.read2.is_none()
    }
}

impl AnnotatedRecord {
    pub fn join(&mut self, other: Self) {
        if let Some(bc) = &mut self.barcode {
            if let Some(x) = other.barcode.as_ref() {
                bc.extend(x)
            }
        } else {
            self.barcode = other.barcode;
        }

        if let Some(umi) = &mut self.umi {
            if let Some(x) = other.umi.as_ref() {
                extend_fastq_record(umi, x)
            }
        } else {
            self.umi = other.umi;
        }

        if self.read1.is_some() {
            if other.read1.is_some() {
                panic!("Read1 already exists");
            }
        } else {
            self.read1 = other.read1;
        }

        if self.read2.is_some() {
            if other.read2.is_some() {
                panic!("Read2 already exists");
            }
        } else {
            self.read2 = other.read2;
        }
    }
}

pub struct VectorChunk<I> {
    inner: I,
    chunk_size: usize,
}

impl<I> VectorChunk<I> {
    pub fn new(inner: I, chunk_size: usize) -> Self {
        Self { inner, chunk_size }
    }
}

impl<I: Iterator<Item = AnnotatedRecord>> Iterator for VectorChunk<I> {
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk = Vec::new();
        let mut accumulated_length = 0;

        for record in self.inner.by_ref() {
            accumulated_length += record.len();
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
    ori_barcode: &[u8],
    ori_qual: &[u8],
    correct_barcode: Option<&[u8]>,
) -> std::io::Result<RecordBuf> {
    let mut record_buf = RecordBuf::try_from_alignment_record(header, record)?;
    let data = record_buf.data_mut();
    data.insert(
        Tag::CELL_BARCODE_SEQUENCE,
        Value::String(ori_barcode.into()),
    );
    data.insert(
        Tag::CELL_BARCODE_QUALITY_SCORES,
        Value::String(ori_qual.into()),
    );
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

pub fn extend_fastq_record(this: &mut fastq::Record, other: &fastq::Record) {
    this.sequence_mut().extend_from_slice(other.sequence());
    this.quality_scores_mut()
        .extend_from_slice(other.quality_scores());
}

fn rev_compl_fastq_record(mut record: fastq::Record) -> fastq::Record {
    *record.quality_scores_mut() = record.quality_scores().iter().rev().copied().collect();
    *record.sequence_mut() = record
        .sequence()
        .iter()
        .rev()
        .map(|&x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => x,
        })
        .collect();
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
        assert!(
            !self.checker.contains(name),
            "bam file must be name collated or name sorted"
        );
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
                panic!(
                    "Expecting paired end reads with the same name, found {} and {}",
                    prev_name, name
                );
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
            PairedEndStats::default(),
        );
        let mut processor = FastqProcessor::new(spec, aligner).with_modality(Modality::ATAC);

        processor.gen_barcoded_alignments().take(6).for_each(|x| {
            println!("{:?}", x);
        });

        println!("{}", processor.get_report());
    }
}
