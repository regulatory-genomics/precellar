/// This module provides an abstraction for aligning sequencing reads using different alignment tools like BWA and STAR.
mod minimap2;
pub use minimap2::{Minimap2Aligner, Minimap2Opts};

use super::fastq::{AlignmentInput, Barcode};
use crate::barcode::{get_barcode, get_umi};

use anyhow::{bail, ensure, Result};
pub use bwa_mem2::BurrowsWheelerAligner;
pub use minibwa::{Index as MiniBwaIndex, MiniBwaSR, Options as MiniBwaOptions};
use noodles_sam::alignment::Record;
pub use star_aligner::StarAligner;

use log;
use noodles_fastq as fastq;
use noodles_sam as sam;
use noodles_sam::alignment::record::data::field::tag::Tag;
use noodles_sam::alignment::record_buf::{data::field::value::Value, RecordBuf};
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSlice;
use std::sync::{Arc, Mutex};

pub type MultiMapR = MultiMap<RecordBuf>;

/// Represents a set of alignments (primary and optional secondary alignments) for a single sequencing read.
#[derive(Debug, Clone)]
pub struct MultiMap<R> {
    /// The primary alignment for the read.
    pub primary: R,
    /// Optional secondary alignments for the read.
    pub others: Option<Vec<R>>,
}

impl<R: Record> MultiMap<R> {
    /// Constructs a new `MultiMap`.
    ///
    /// # Arguments
    /// * `primary` - The primary alignment for the read.
    /// * `others` - Optional secondary alignments for the read.
    pub fn new(primary: R, others: Option<Vec<R>>) -> Self {
        Self { primary, others }
    }

    /// Return the number of records.
    pub fn len(&self) -> usize {
        self.others.as_ref().map_or(0, |x| x.len()) + 1
    }

    /// Return the cell barcode if it exists.
    pub fn barcode(&self) -> Result<Option<String>> {
        get_barcode(&self.primary)
    }

    /// Return the UMI if it exists.
    pub fn umi(&self) -> Result<Option<String>> {
        get_umi(&self.primary)
    }

    /// Whether the read is confidently mapped. A read is confidently mapped if it
    /// is mapped to a single location.
    pub fn is_confidently_mapped(&self) -> bool {
        self.others.is_none() && !self.primary.flags().unwrap().is_unmapped()
    }

    /// Returns an iterator over all alignments (primary and secondary).
    pub fn iter(&self) -> impl Iterator<Item = &R> {
        std::iter::once(&self.primary).chain(self.others.iter().flatten())
    }
}

impl<R> From<R> for MultiMap<R> {
    fn from(record: R) -> Self {
        Self {
            primary: record,
            others: None,
        }
    }
}

impl<R: Record> TryFrom<Vec<R>> for MultiMap<R> {
    type Error = anyhow::Error;

    fn try_from(mut vec: Vec<R>) -> Result<Self, Self::Error> {
        let n = vec.len();
        if n == 0 {
            Err(anyhow::anyhow!("No alignments"))
        } else if n == 1 {
            Ok(MultiMap::from(vec.into_iter().next().unwrap()))
        } else {
            let mut primary = None;
            vec.iter().enumerate().try_for_each(|(i, rec)| {
                if !rec.flags()?.is_secondary() {
                    if primary.is_some() {
                        bail!("Multiple primary alignments");
                    } else {
                        primary = Some(i);
                    }
                }
                Ok(())
            })?;
            ensure!(primary.is_some(), "No primary alignment");

            Ok(MultiMap::new(vec.swap_remove(primary.unwrap()), Some(vec)))
        }
    }
}

/// Records per alignment task. Small enough that rayon can rebalance work
/// between threads while a batch runs, large enough that per-task overhead
/// stays negligible.
const ALIGN_CHUNK_SIZE: usize = 1024;

/// Thread pool for chunked alignment, sized by the `num_threads` of the latest
/// call. It is kept across calls so that every batch reuses the same threads.
static ALIGN_POOL: Mutex<Option<Arc<rayon::ThreadPool>>> = Mutex::new(None);

fn align_pool(num_threads: usize) -> Arc<rayon::ThreadPool> {
    let mut slot = ALIGN_POOL.lock().unwrap_or_else(|e| e.into_inner());
    match slot.as_ref() {
        Some(pool) if pool.current_num_threads() == num_threads => pool.clone(),
        _ => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .thread_name(|i| format!("precellar-align-{i}"))
                .build()
                .expect("failed to build the alignment thread pool");
            let pool = Arc::new(pool);
            *slot = Some(pool.clone());
            pool
        }
    }
}

/// Aligns `records` on exactly `num_threads` threads, in tasks of
/// [`ALIGN_CHUNK_SIZE`] records. Each thread clones `aligner` at most once per
/// call and reuses the clone for all its tasks. The output keeps input order.
fn align_chunked<A, F>(
    aligner: &A,
    num_threads: u16,
    records: &[AlignmentInput],
    align_one: F,
) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)>
where
    A: Clone + Send + Sync,
    F: Fn(&mut A, &AlignmentInput) -> (Option<MultiMapR>, Option<MultiMapR>) + Sync,
{
    let pool = align_pool(usize::from(num_threads.max(1)));
    // One slot per pool thread, indexed by `rayon::current_thread_index`. A
    // thread runs one task at a time, so its lock is never contended.
    let workers: Vec<Mutex<Option<A>>> = (0..pool.current_num_threads())
        .map(|_| Mutex::new(None))
        .collect();
    pool.install(|| {
        records
            .par_chunks(ALIGN_CHUNK_SIZE)
            .map(|chunk| {
                let index = rayon::current_thread_index()
                    .expect("alignment tasks run on the alignment pool");
                let mut worker = workers[index].lock().unwrap();
                let worker = worker.get_or_insert_with(|| aligner.clone());
                chunk
                    .iter()
                    .map(|rec| align_one(worker, rec))
                    .collect::<Vec<_>>()
            })
            .flatten()
            .collect()
    })
}

/// Trait defining the behavior of aligners like BWA and STAR.
pub trait Aligner {
    /// Retrieves the SAM header associated with the aligner.
    fn header(&self) -> sam::Header;

    /// Aligns a batch of sequencing reads.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads to use for alignment.
    /// * `records` - Vector of annotated FASTQ records to align.
    ///
    /// # Returns
    /// A vector of tuples where each tuple contains the alignments.
    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)>;
}

/// Select one single-end read and retain whether it came from read 1.
fn select_single_read<T>(read1: Option<T>, read2: Option<T>) -> Option<(T, bool)> {
    match (read1, read2) {
        (Some(read), None) => Some((read, true)),
        (None, Some(read)) => Some((read, false)),
        (Some(_), Some(_)) | (None, None) => None,
    }
}

fn restore_single_slot<T>(value: T, is_read1: bool) -> (Option<T>, Option<T>) {
    if is_read1 {
        (Some(value), None)
    } else {
        (None, Some(value))
    }
}

impl Aligner for BurrowsWheelerAligner {
    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        if records.is_empty() {
            return Vec::new();
        }
        let all_paired = records
            .iter()
            .all(|record| record.read1.is_some() && record.read2.is_some());
        let all_single = records
            .iter()
            .all(|record| record.read1.is_some() ^ record.read2.is_some());
        assert!(
            all_paired || all_single,
            "an alignment batch cannot mix paired and single-end records"
        );

        if all_paired {
            let (info, mut reads): (Vec<_>, Vec<_>) = records
                .into_iter()
                .map(|rec| {
                    let AlignmentInput {
                        read1,
                        read2,
                        metadata,
                    } = rec;
                    (
                        (metadata.barcode, metadata.umi),
                        (read1.unwrap(), read2.unwrap()),
                    )
                })
                .unzip();
            self.align_read_pairs(num_threads, &mut reads)
                .enumerate()
                .map(|(i, (mut ali1, mut ali2))| {
                    let (bc, umi) = info.get(i).unwrap();
                    attach_read_metadata(&mut ali1, bc, umi.as_ref());
                    attach_read_metadata(&mut ali2, bc, umi.as_ref());
                    (Some(ali1.into()), Some(ali2.into()))
                })
                .collect()
        } else {
            let (info, mut reads): (Vec<_>, Vec<_>) = records
                .into_iter()
                .map(|rec| {
                    let AlignmentInput {
                        read1,
                        read2,
                        metadata,
                    } = rec;
                    let (read, slot) = select_single_read(read1, read2)
                        .expect("single-end batch contains an invalid read layout");
                    ((metadata.barcode, metadata.umi, slot), read)
                })
                .unzip();

            self.align_reads(num_threads, reads.as_mut_slice())
                .enumerate()
                .map(|(i, mut alignment)| {
                    let (bc, umi, slot) = info.get(i).unwrap();
                    attach_read_metadata(&mut alignment, bc, umi.as_ref());
                    restore_single_slot(alignment.into(), *slot)
                })
                .collect()
        }
    }
}

impl Aligner for MiniBwaSR {
    fn header(&self) -> sam::Header {
        self.get_sam_header()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        if records.is_empty() {
            return Vec::new();
        }
        let all_paired = records
            .iter()
            .all(|record| record.read1.is_some() && record.read2.is_some());
        let all_single = records
            .iter()
            .all(|record| record.read1.is_some() ^ record.read2.is_some());
        assert!(
            all_paired || all_single,
            "an alignment batch cannot mix paired and single-end records"
        );

        if all_paired {
            let (info, mut reads): (Vec<_>, Vec<_>) = records
                .into_iter()
                .map(|rec| {
                    let AlignmentInput {
                        read1,
                        read2,
                        metadata,
                    } = rec;
                    (
                        (metadata.barcode, metadata.umi),
                        (read1.unwrap(), read2.unwrap()),
                    )
                })
                .unzip();

            self.align_read_pairs(num_threads, &mut reads)
                .unwrap()
                .enumerate()
                .map(|(i, (ali1, ali2))| {
                    // Extract the primary alignment and discard the rest
                    let mut ali1 = ali1.into_iter().next().unwrap();
                    let mut ali2 = ali2.into_iter().next().unwrap();
                    let (bc, umi) = info.get(i).unwrap();
                    attach_read_metadata(&mut ali1, bc, umi.as_ref());
                    attach_read_metadata(&mut ali2, bc, umi.as_ref());
                    (Some(ali1.into()), Some(ali2.into()))
                })
                .collect()
        } else {
            let (info, mut reads): (Vec<_>, Vec<_>) = records
                .into_iter()
                .map(|rec| {
                    let AlignmentInput {
                        read1,
                        read2,
                        metadata,
                    } = rec;
                    let (read, slot) = select_single_read(read1, read2)
                        .expect("single-end batch contains an invalid read layout");
                    ((metadata.barcode, metadata.umi, slot), read)
                })
                .unzip();

            self.align_reads(num_threads, reads.as_mut_slice())
                .unwrap()
                .enumerate()
                .map(|(i, alignment)| {
                    // Extract the primary alignment and discard the rest
                    let mut alignment = alignment.into_iter().next().unwrap();
                    let (bc, umi, slot) = info.get(i).unwrap();
                    attach_read_metadata(&mut alignment, bc, umi.as_ref());
                    restore_single_slot(alignment.into(), *slot)
                })
                .collect()
        }
    }
}

impl Aligner for StarAligner {
    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        align_chunked(self, num_threads, &records, |aligner, rec| {
            let bc = &rec.metadata.barcode;
            let read1 = rec.read1.as_ref();
            let read2 = rec.read2.as_ref();

            if read1.is_some() && read2.is_some() {
                let (mut ali1, mut ali2) = aligner
                    .align_read_pair(&read1.unwrap(), &read2.unwrap())
                    .unwrap();
                ali1.iter_mut()
                    .chain(ali2.iter_mut())
                    .for_each(|alignment| {
                        attach_read_metadata(alignment, bc, rec.metadata.umi.as_ref());
                    });
                (
                    Some(ali1.try_into().unwrap()),
                    Some(ali2.try_into().unwrap()),
                )
            } else if let Some((read, slot)) = select_single_read(read1, read2) {
                let mut ali = aligner.align_read(read).unwrap();
                ali.iter_mut().for_each(|alignment| {
                    attach_read_metadata(alignment, bc, rec.metadata.umi.as_ref());
                });
                restore_single_slot(ali.try_into().unwrap(), slot)
            } else {
                log::warn!(
                    "Found record with no reads (read1 and read2 are both None). Barcode: {:?}",
                    String::from_utf8_lossy(bc.raw.sequence())
                );
                (None, None)
            }
        })
    }
}

impl Aligner for Minimap2Aligner {
    fn header(&self) -> sam::Header {
        self.get_header().clone()
    }

    fn align_reads(
        &mut self,
        num_threads: u16,
        records: Vec<AlignmentInput>,
    ) -> Vec<(Option<MultiMapR>, Option<MultiMapR>)> {
        // Cloning is efficient: only clones Arc pointers to the shared index.
        align_chunked(self, num_threads, &records, |thread_aligner, rec| {
            let bc = &rec.metadata.barcode;
            let read1 = rec.read1.as_ref();
            let read2 = rec.read2.as_ref();

            if read1.is_some() && read2.is_some() {
                let (mut ali1, mut ali2) = thread_aligner
                    .align_read_pair(&read1.unwrap(), &read2.unwrap())
                    .unwrap();
                ali1.iter_mut()
                    .chain(ali2.iter_mut())
                    .for_each(|alignment| {
                        attach_read_metadata(alignment, bc, rec.metadata.umi.as_ref());
                    });
                (
                    Some(ali1.try_into().unwrap()),
                    Some(ali2.try_into().unwrap()),
                )
            } else if let Some((read, slot)) = select_single_read(read1, read2) {
                let mut ali = thread_aligner.align_read(read).unwrap();
                ali.iter_mut().for_each(|alignment| {
                    attach_read_metadata(alignment, bc, rec.metadata.umi.as_ref());
                });
                restore_single_slot(ali.try_into().unwrap(), slot)
            } else {
                log::warn!(
                    "Found record with no reads (read1 and read2 are both None). Barcode: {:?}",
                    String::from_utf8_lossy(bc.raw.sequence())
                );
                (None, None)
            }
        })
    }
}

// Centralize FASTQ metadata policy so every alignment backend emits the same SAM tags.
fn attach_read_metadata(
    record_buf: &mut RecordBuf,
    barcode: &Barcode,
    umi: Option<&fastq::Record>,
) {
    let data = record_buf.data_mut();
    data.insert(
        Tag::CELL_BARCODE_SEQUENCE,
        Value::String(barcode.raw.sequence().into()),
    );
    data.insert(
        Tag::CELL_BARCODE_QUALITY_SCORES,
        Value::String(barcode.raw.quality_scores().into()),
    );

    if let Some(corrected) = barcode.corrected.as_deref() {
        data.insert(Tag::CELL_BARCODE_ID, Value::String(corrected.into()));
    }
    if let Some(umi) = umi {
        data.insert(Tag::UMI_SEQUENCE, Value::String(umi.sequence().into()));
        data.insert(
            Tag::UMI_QUALITY_SCORES,
            Value::String(umi.quality_scores().into()),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::fastq::ReadMetadata;
    use bwa_mem2::{AlignerOpts, FMIndex};
    use noodles_fastq::record::Definition;

    fn reference_sequence(len: usize) -> Vec<u8> {
        let mut state = 0x1234_5678_u32;
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[((state >> 16) & 3) as usize]
            })
            .collect()
    }

    fn write_reference(directory: &std::path::Path) -> (std::path::PathBuf, Vec<u8>) {
        let fasta = directory.join("reference.fa");
        let sequence = reference_sequence(2_000);
        std::fs::write(
            &fasta,
            format!(">chr1\n{}\n", String::from_utf8_lossy(&sequence)),
        )
        .unwrap();
        (fasta, sequence)
    }

    fn fastq_record(name: &str, sequence: &[u8]) -> fastq::Record {
        fastq::Record::new(
            Definition::new(name.as_bytes(), ""),
            sequence.to_vec(),
            vec![b'I'; sequence.len()],
        )
    }

    fn alignment_input(name: &str, sequence: &[u8], is_read1: bool) -> AlignmentInput {
        let read = fastq_record(name, sequence);
        let (read1, read2) = restore_single_slot(read, is_read1);
        AlignmentInput {
            read1,
            read2,
            metadata: ReadMetadata {
                barcode: Barcode {
                    raw: fastq_record("barcode", b"ACGT"),
                    corrected: Some(b"ACGT".to_vec()),
                },
                umi: None,
            },
        }
    }

    fn mixed_single_end_inputs(sequence: &[u8]) -> Vec<AlignmentInput> {
        vec![
            alignment_input("forward", sequence, true),
            alignment_input("reverse", sequence, false),
        ]
    }

    fn assert_mixed_single_end_slots(output: &[(Option<MultiMapR>, Option<MultiMapR>)]) {
        assert_eq!(output.len(), 2);
        assert!(output[0].0.is_some());
        assert!(output[0].1.is_none());
        assert!(output[1].0.is_none());
        assert!(output[1].1.is_some());
        assert_eq!(
            output[0].0.as_ref().unwrap().barcode().unwrap().as_deref(),
            Some("ACGT")
        );
        assert_eq!(
            output[1].1.as_ref().unwrap().barcode().unwrap().as_deref(),
            Some("ACGT")
        );
    }

    /// Mock aligner that counts its clones.
    struct CountingAligner {
        clones: Arc<std::sync::atomic::AtomicUsize>,
    }

    impl Clone for CountingAligner {
        fn clone(&self) -> Self {
            self.clones
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            Self {
                clones: self.clones.clone(),
            }
        }
    }

    #[test]
    fn align_chunked_caps_threads_and_keeps_order() {
        let inputs: Vec<_> = (0..10 * ALIGN_CHUNK_SIZE)
            .map(|i| alignment_input(&i.to_string(), b"ACGT", true))
            .collect();
        for num_threads in [3u16, 2] {
            let aligner = CountingAligner {
                clones: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            };
            let threads = Mutex::new(std::collections::HashSet::new());
            let output = align_chunked(&aligner, num_threads, &inputs, |_, rec| {
                threads.lock().unwrap().insert(std::thread::current().id());
                let mut record = RecordBuf::default();
                *record.name_mut() = Some(rec.read1.as_ref().unwrap().name().to_vec().into());
                (Some(MultiMap::new(record, None)), None)
            });

            let names: Vec<_> = output
                .iter()
                .map(|(ali, _)| ali.as_ref().unwrap().primary.name().unwrap().to_string())
                .collect();
            let expected: Vec<_> = (0..inputs.len()).map(|i| i.to_string()).collect();
            assert_eq!(names, expected);
            let used = threads.lock().unwrap().len();
            assert!(
                used <= num_threads as usize,
                "{used} threads for {num_threads}"
            );
            let clones = aligner.clones.load(std::sync::atomic::Ordering::Relaxed);
            assert!(
                clones <= num_threads as usize,
                "{clones} clones for {num_threads}"
            );
        }
    }

    #[test]
    fn single_read_helpers_preserve_read1_and_read2_slots() {
        assert_eq!(select_single_read(Some("r1"), None), Some(("r1", true)));
        assert_eq!(select_single_read(None, Some("r2")), Some(("r2", false)));
        assert_eq!(select_single_read(Some("r1"), Some("r2")), None);
        assert_eq!(select_single_read::<&str>(None, None), None);

        assert_eq!(
            restore_single_slot("alignment", true),
            (Some("alignment"), None)
        );
        assert_eq!(
            restore_single_slot("alignment", false),
            (None, Some("alignment"))
        );
    }

    #[test]
    fn minimap2_smoke_preserves_mixed_single_end_slots() {
        let directory = tempfile::tempdir().unwrap();
        let (fasta, reference) = write_reference(directory.path());
        let query = &reference[400..1_000];
        let opts = Minimap2Opts::new(fasta).with_preset(::minimap2::Preset::MapOnt);
        let mut aligner = Minimap2Aligner::new(opts).unwrap();

        let output = Aligner::align_reads(&mut aligner, 1, mixed_single_end_inputs(query));
        assert_mixed_single_end_slots(&output);
    }

    #[test]
    fn minibwa_smoke_preserves_mixed_single_end_slots() {
        let directory = tempfile::tempdir().unwrap();
        let (fasta, reference) = write_reference(directory.path());
        let prefix = directory.path().join("minibwa-index");
        let index = MiniBwaIndex::build(&fasta, &prefix, 1, false).unwrap();
        let mut options = MiniBwaOptions::default();
        options.set_min_seed_len(8);
        options.set_min_chain_score(1);
        options.set_min_dp_score(1);
        let mut aligner = MiniBwaSR::new(index, options).unwrap();

        let output = Aligner::align_reads(
            &mut aligner,
            1,
            mixed_single_end_inputs(&reference[400..500]),
        );
        assert_mixed_single_end_slots(&output);
    }

    #[test]
    fn bwa_mem2_smoke_preserves_mixed_single_end_slots() {
        let directory = tempfile::tempdir().unwrap();
        let (fasta, reference) = write_reference(directory.path());
        let prefix = directory.path().join("bwa-index");
        let index = FMIndex::new(&fasta, &prefix).unwrap();
        let mut options = AlignerOpts::default();
        options.set_min_seed_len(8);
        let mut aligner = BurrowsWheelerAligner::new(index, options);

        let output = Aligner::align_reads(
            &mut aligner,
            1,
            mixed_single_end_inputs(&reference[400..500]),
        );
        assert_mixed_single_end_slots(&output);
    }

    #[test]
    #[ignore = "requires a prebuilt STAR index in PRECELLAR_TEST_STAR_INDEX"]
    fn star_smoke_preserves_mixed_single_end_slots() {
        let index = std::env::var_os("PRECELLAR_TEST_STAR_INDEX")
            .expect("PRECELLAR_TEST_STAR_INDEX must point to a prebuilt STAR index");
        let mut aligner = StarAligner::new(star_aligner::StarOpts::new(index)).unwrap();
        let sequence = reference_sequence(100);

        let output = Aligner::align_reads(&mut aligner, 1, mixed_single_end_inputs(&sequence));
        assert_mixed_single_end_slots(&output);
    }
}
