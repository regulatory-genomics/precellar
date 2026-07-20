use crate::align::AnnotatedFastq;
use anyhow::{ensure, Context, Result};
use piscem_rs::index::contig_table::ContigTable;
use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::chain_state::SketchHitInfoChained;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy as PiscemSkippingStrategy};
use piscem_rs::mapping::hits::{
    MappingType, SimpleHit, SketchHitInfo, INVALID_FRAG_LEN, INVALID_MATE_POS,
};
use piscem_rs::mapping::map_fragment::{map_pe_fragment, map_se_fragment};
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use rayon::prelude::*;
use sshash_lib::{dispatch_on_k, Dictionary, Kmer, KmerBits};
use std::path::Path;

/// How piscem advances between informative k-mers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkippingStrategy {
    Strict,
    Permissive,
}

impl From<SkippingStrategy> for PiscemSkippingStrategy {
    fn from(value: SkippingStrategy) -> Self {
        match value {
            SkippingStrategy::Strict => Self::Strict,
            SkippingStrategy::Permissive => Self::Permissive,
        }
    }
}

/// Configuration for the in-memory piscem mapper.
#[derive(Debug, Clone)]
pub struct PiscemOptions {
    pub num_threads: usize,
    pub skipping_strategy: SkippingStrategy,
    pub max_hit_occ: usize,
    pub max_hit_occ_recover: usize,
    pub max_read_occ: usize,
    pub load_ec: bool,
    pub max_ec_card: u32,
    pub load_poison: bool,
    pub structural_constraints: bool,
}

impl Default for PiscemOptions {
    fn default() -> Self {
        Self {
            num_threads: 16,
            skipping_strategy: SkippingStrategy::Permissive,
            max_hit_occ: 256,
            max_hit_occ_recover: 1024,
            max_read_occ: 2500,
            load_ec: true,
            max_ec_card: 4096,
            load_poison: true,
            structural_constraints: false,
        }
    }
}

/// Piscem's classification of a biological read or read pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PseudoalignmentType {
    Unmapped,
    Single,
    FirstOrphan,
    SecondOrphan,
    Pair,
}

impl From<MappingType> for PseudoalignmentType {
    fn from(value: MappingType) -> Self {
        match value {
            MappingType::Unmapped => Self::Unmapped,
            MappingType::SingleMapped => Self::Single,
            MappingType::MappedFirstOrphan => Self::FirstOrphan,
            MappingType::MappedSecondOrphan => Self::SecondOrphan,
            MappingType::MappedPair => Self::Pair,
        }
    }
}

/// Outcome for one annotated input record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PseudoalignmentStatus {
    Mapped,
    Unmapped,
    Poisoned,
    MissingBarcode,
    MissingCorrectedBarcode,
    InvalidBarcode,
    MissingUmi,
    InvalidUmi,
    MissingBiologicalRead,
}

/// One accepted mapping target returned by piscem.
#[derive(Debug, Clone, PartialEq)]
pub struct PseudoalignmentHit {
    pub target_id: u32,
    pub forward: bool,
    pub position: i32,
    pub mate_forward: Option<bool>,
    pub mate_position: Option<i32>,
    pub fragment_length: Option<i32>,
    /// Number of supporting k-mer hits. Proper pairs currently report zero in piscem.
    pub num_hits: u32,
    /// Piscem support score. Proper pairs currently report zero.
    pub score: f32,
}

impl From<SimpleHit> for PseudoalignmentHit {
    fn from(hit: SimpleHit) -> Self {
        let has_mate = hit.mate_pos != INVALID_MATE_POS;
        Self {
            target_id: hit.tid,
            forward: hit.is_fw,
            position: hit.pos,
            mate_forward: has_mate.then_some(hit.mate_is_fw),
            mate_position: has_mate.then_some(hit.mate_pos),
            fragment_length: (hit.fragment_length != INVALID_FRAG_LEN)
                .then_some(hit.fragment_length),
            num_hits: hit.num_hits,
            score: hit.score,
        }
    }
}

/// In-memory pseudoalignment for one input record.
#[derive(Debug, Clone, PartialEq)]
pub struct Pseudoalignment {
    pub barcode: Option<Vec<u8>>,
    pub umi: Option<Vec<u8>>,
    pub status: PseudoalignmentStatus,
    pub mapping_type: PseudoalignmentType,
    pub hits: Vec<PseudoalignmentHit>,
}

impl Pseudoalignment {
    fn rejected(
        barcode: Option<Vec<u8>>,
        umi: Option<Vec<u8>>,
        status: PseudoalignmentStatus,
    ) -> Self {
        Self {
            barcode,
            umi,
            status,
            mapping_type: PseudoalignmentType::Unmapped,
            hits: Vec::new(),
        }
    }
}

/// Aggregate counters for an in-memory mapping batch.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PiscemStats {
    pub input_records: u64,
    pub submitted_records: u64,
    pub mapped_records: u64,
    pub unmapped_records: u64,
    pub poisoned_records: u64,
    pub missing_barcode: u64,
    pub missing_corrected_barcode: u64,
    pub invalid_barcode: u64,
    pub missing_umi: u64,
    pub invalid_umi: u64,
    pub missing_biological_read: u64,
}

impl PiscemStats {
    fn from_records(records: &[Pseudoalignment]) -> Self {
        let mut stats = Self {
            input_records: records.len() as u64,
            ..Self::default()
        };
        for record in records {
            match record.status {
                PseudoalignmentStatus::Mapped => {
                    stats.submitted_records += 1;
                    stats.mapped_records += 1;
                }
                PseudoalignmentStatus::Unmapped => {
                    stats.submitted_records += 1;
                    stats.unmapped_records += 1;
                }
                PseudoalignmentStatus::Poisoned => {
                    stats.submitted_records += 1;
                    stats.poisoned_records += 1;
                }
                PseudoalignmentStatus::MissingBarcode => stats.missing_barcode += 1,
                PseudoalignmentStatus::MissingCorrectedBarcode => {
                    stats.missing_corrected_barcode += 1
                }
                PseudoalignmentStatus::InvalidBarcode => stats.invalid_barcode += 1,
                PseudoalignmentStatus::MissingUmi => stats.missing_umi += 1,
                PseudoalignmentStatus::InvalidUmi => stats.invalid_umi += 1,
                PseudoalignmentStatus::MissingBiologicalRead => stats.missing_biological_read += 1,
            }
        }
        stats
    }
}

/// Ordered pseudoalignments and their aggregate counters.
#[derive(Debug, Clone, PartialEq)]
pub struct PiscemBatch {
    pub records: Vec<Pseudoalignment>,
    pub stats: PiscemStats,
}

/// Loaded piscem index and a dedicated mapping thread pool.
pub struct PiscemMapper {
    index: ReferenceIndex<Dictionary, ContigTable>,
    options: PiscemOptions,
    thread_pool: rayon::ThreadPool,
}

impl PiscemMapper {
    /// Load a standard SSHash-backed piscem index.
    pub fn open(index_prefix: impl AsRef<Path>, options: PiscemOptions) -> Result<Self> {
        validate_options(&options)?;
        let index_prefix = index_prefix.as_ref();
        let index = ReferenceIndex::load(index_prefix, options.load_ec, options.load_poison)
            .with_context(|| {
                format!(
                    "failed to load piscem index from {}",
                    index_prefix.display()
                )
            })?;
        ensure!(
            is_supported_k(index.k()),
            "unsupported piscem k-mer size {}; expected an odd value from 3 through 63",
            index.k()
        );
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(options.num_threads)
            .build()
            .context("failed to create piscem mapping thread pool")?;
        Ok(Self {
            index,
            options,
            thread_pool,
        })
    }

    pub fn k(&self) -> usize {
        self.index.k()
    }

    pub fn m(&self) -> usize {
        self.index.m()
    }

    pub fn num_targets(&self) -> usize {
        self.index.num_refs()
    }

    pub fn has_ec_table(&self) -> bool {
        self.index.has_ec_table()
    }

    pub fn target_name(&self, target_id: usize) -> Option<&str> {
        (target_id < self.index.num_refs()).then(|| self.index.ref_name(target_id))
    }

    pub fn target_length(&self, target_id: usize) -> Option<u64> {
        (target_id < self.index.num_refs()).then(|| self.index.ref_len(target_id))
    }

    /// Pseudoalign an annotated batch without serializing it as FASTQ.
    pub fn map_reads(&self, records: Vec<AnnotatedFastq>) -> PiscemBatch {
        let k = self.index.k();
        let mapped = dispatch_on_k!(k, K => {
            if self.options.structural_constraints {
                self.map_reads_with::<K, SketchHitInfoChained>(records)
            } else {
                self.map_reads_with::<K, SketchHitInfoSimple>(records)
            }
        });
        let stats = PiscemStats::from_records(&mapped);
        PiscemBatch {
            records: mapped,
            stats,
        }
    }

    fn map_reads_with<const K: usize, S>(
        &self,
        records: Vec<AnnotatedFastq>,
    ) -> Vec<Pseudoalignment>
    where
        Kmer<K>: KmerBits,
        S: SketchHitInfo + Send,
    {
        self.thread_pool.install(|| {
            records
                .into_par_iter()
                .map_init(
                    || PiscemWorker::<K, S>::new(&self.index, &self.options),
                    PiscemWorker::map_record,
                )
                .collect()
        })
    }
}

fn validate_options(options: &PiscemOptions) -> Result<()> {
    ensure!(
        options.num_threads > 0,
        "num_threads must be greater than zero"
    );
    ensure!(
        options.max_hit_occ > 0,
        "max_hit_occ must be greater than zero"
    );
    Ok(())
}

fn is_supported_k(k: usize) -> bool {
    (3..=63).contains(&k) && k % 2 == 1
}

struct PiscemWorker<'a, const K: usize, S>
where
    Kmer<K>: KmerBits,
    S: SketchHitInfo,
{
    index: &'a ReferenceIndex<Dictionary, ContigTable>,
    hit_searcher: HitSearcher<'a, Dictionary, ContigTable>,
    query: PiscemStreamingQuery<'a, K, Dictionary>,
    cache_out: MappingCache<S>,
    cache_left: MappingCache<S>,
    cache_right: MappingCache<S>,
    poison_state: PoisonState<'a>,
    skipping_strategy: PiscemSkippingStrategy,
}

impl<'a, const K: usize, S> PiscemWorker<'a, K, S>
where
    Kmer<K>: KmerBits,
    S: SketchHitInfo,
{
    fn new(index: &'a ReferenceIndex<Dictionary, ContigTable>, options: &PiscemOptions) -> Self {
        let mut cache_out = MappingCache::new(K);
        let mut cache_left = MappingCache::new(K);
        let mut cache_right = MappingCache::new(K);
        for cache in [&mut cache_out, &mut cache_left, &mut cache_right] {
            apply_options(cache, options);
        }
        Self {
            index,
            hit_searcher: HitSearcher::new(index),
            query: PiscemStreamingQuery::new(index.dict()),
            cache_out,
            cache_left,
            cache_right,
            poison_state: PoisonState::new(index.poison_table()),
            skipping_strategy: options.skipping_strategy.into(),
        }
    }

    fn map_record(&mut self, record: AnnotatedFastq) -> Pseudoalignment {
        let PreparedRecord {
            barcode,
            umi,
            read1,
            read2,
        } = match prepare_record(record) {
            Ok(record) => record,
            Err(rejected) => return rejected,
        };

        let mapping_type = match (read1.as_ref(), read2.as_ref()) {
            (Some(read1), Some(read2)) => {
                self.poison_state.paired_for_mapping = true;
                map_pe_fragment(
                    read1.sequence(),
                    read2.sequence(),
                    &mut self.hit_searcher,
                    &mut self.query,
                    &mut self.cache_left,
                    &mut self.cache_right,
                    &mut self.cache_out,
                    self.index,
                    &mut self.poison_state,
                    self.skipping_strategy,
                );
                self.cache_out.map_type
            }
            (Some(read), None) | (None, Some(read)) => {
                self.poison_state.paired_for_mapping = false;
                map_se_fragment(
                    read.sequence(),
                    &mut self.hit_searcher,
                    &mut self.query,
                    &mut self.cache_out,
                    self.index,
                    &mut self.poison_state,
                    self.skipping_strategy,
                );
                self.cache_out.map_type
            }
            (None, None) => unreachable!("validated records contain a biological read"),
        };

        if self.poison_state.is_poisoned() {
            return Pseudoalignment::rejected(
                Some(barcode),
                Some(umi),
                PseudoalignmentStatus::Poisoned,
            );
        }

        let status = if mapping_type == MappingType::Unmapped {
            PseudoalignmentStatus::Unmapped
        } else {
            PseudoalignmentStatus::Mapped
        };
        let mut hits: Vec<_> = self
            .cache_out
            .accepted_hits
            .iter()
            .copied()
            .map(PseudoalignmentHit::from)
            .collect();
        hits.sort_by(compare_hits);
        Pseudoalignment {
            barcode: Some(barcode),
            umi: Some(umi),
            status,
            mapping_type: mapping_type.into(),
            hits,
        }
    }
}

#[derive(Debug)]
struct PreparedRecord {
    barcode: Vec<u8>,
    umi: Vec<u8>,
    read1: Option<noodles_fastq::Record>,
    read2: Option<noodles_fastq::Record>,
}

fn prepare_record(record: AnnotatedFastq) -> std::result::Result<PreparedRecord, Pseudoalignment> {
    let AnnotatedFastq {
        barcode,
        umi,
        read1,
        read2,
    } = record;

    let Some(barcode) = barcode else {
        return Err(Pseudoalignment::rejected(
            None,
            None,
            PseudoalignmentStatus::MissingBarcode,
        ));
    };
    let Some(barcode) = barcode.corrected else {
        return Err(Pseudoalignment::rejected(
            None,
            None,
            PseudoalignmentStatus::MissingCorrectedBarcode,
        ));
    };
    if !is_acgt(&barcode) {
        return Err(Pseudoalignment::rejected(
            Some(barcode),
            None,
            PseudoalignmentStatus::InvalidBarcode,
        ));
    }

    let Some(umi) = umi else {
        return Err(Pseudoalignment::rejected(
            Some(barcode),
            None,
            PseudoalignmentStatus::MissingUmi,
        ));
    };
    let umi = umi.sequence().to_vec();
    if !is_acgt(&umi) {
        return Err(Pseudoalignment::rejected(
            Some(barcode),
            Some(umi),
            PseudoalignmentStatus::InvalidUmi,
        ));
    }

    // A record with one nonempty biological read is intentionally mapped as single-end.
    let read1 = read1.filter(|read| !read.sequence().is_empty());
    let read2 = read2.filter(|read| !read.sequence().is_empty());
    if read1.is_none() && read2.is_none() {
        return Err(Pseudoalignment::rejected(
            Some(barcode),
            Some(umi),
            PseudoalignmentStatus::MissingBiologicalRead,
        ));
    }

    Ok(PreparedRecord {
        barcode,
        umi,
        read1,
        read2,
    })
}

fn apply_options<S: SketchHitInfo>(cache: &mut MappingCache<S>, options: &PiscemOptions) {
    cache.max_hit_occ = options.max_hit_occ;
    cache.max_hit_occ_recover = options.max_hit_occ_recover;
    cache.attempt_occ_recover = options.max_hit_occ_recover > options.max_hit_occ;
    cache.max_read_occ = options.max_read_occ;
    cache.max_ec_card = options.max_ec_card;
}

fn compare_hits(left: &PseudoalignmentHit, right: &PseudoalignmentHit) -> std::cmp::Ordering {
    left.target_id
        .cmp(&right.target_id)
        .then_with(|| right.forward.cmp(&left.forward))
        .then_with(|| left.position.cmp(&right.position))
        .then_with(|| left.mate_position.cmp(&right.mate_position))
        .then_with(|| left.mate_forward.cmp(&right.mate_forward))
        .then_with(|| left.fragment_length.cmp(&right.fragment_length))
        .then_with(|| left.num_hits.cmp(&right.num_hits))
        .then_with(|| left.score.total_cmp(&right.score))
}

fn is_acgt(sequence: &[u8]) -> bool {
    !sequence.is_empty()
        && sequence
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::Barcode;
    use noodles_fastq as fastq;
    use piscem_rs::index::build::{build_index, BuildConfig};
    use piscem_rs::index::poison_table::{LabeledPoisonOcc, PoisonTable};
    use piscem_rs::mapping::kmer_value::CanonicalKmer;
    use std::io::BufWriter;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn fastq_record(sequence: &[u8]) -> fastq::Record {
        fastq::Record::new(
            fastq::record::Definition::default(),
            sequence,
            vec![b'I'; sequence.len()],
        )
    }

    fn annotated_record(
        barcode: Option<Option<&[u8]>>,
        umi: Option<&[u8]>,
        read1: Option<&[u8]>,
        read2: Option<&[u8]>,
    ) -> AnnotatedFastq {
        AnnotatedFastq {
            barcode: barcode.map(|corrected| Barcode {
                raw: fastq_record(corrected.unwrap_or(b"AAAA")),
                corrected: corrected.map(Vec::from),
            }),
            umi: umi.map(fastq_record),
            read1: read1.map(fastq_record),
            read2: read2.map(fastq_record),
        }
    }

    fn build_test_index() -> (TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let input_prefix = dir.path().join("tiny_input");
        let output_prefix = dir.path().join("tiny_index");
        std::fs::write(
            input_prefix.with_extension("cf_seg"),
            "0\tACGTTGCACTGATCGTACGA\n",
        )
        .unwrap();
        std::fs::write(
            input_prefix.with_extension("cf_seq"),
            "Reference:N_Sequence:tx0\t0+\n",
        )
        .unwrap();
        std::fs::write(input_prefix.with_extension("json"), "{}\n").unwrap();
        build_index(&BuildConfig {
            input_prefix,
            output_prefix: output_prefix.clone(),
            k: 5,
            m: 3,
            build_ec_table: true,
            num_threads: 1,
            canonical: false,
            seed: 1,
            single_mphf: true,
            emit_tiny: Some(false),
            tmp_dir: Some(dir.path().join("sshash_tmp")),
            ram_limit_gib: Some(1),
        })
        .unwrap();
        (dir, output_prefix)
    }

    #[test]
    fn converts_mapping_types() {
        assert_eq!(
            PseudoalignmentType::from(MappingType::MappedPair),
            PseudoalignmentType::Pair
        );
        assert_eq!(
            PseudoalignmentType::from(MappingType::MappedSecondOrphan),
            PseudoalignmentType::SecondOrphan
        );
    }

    #[test]
    fn converts_hits_without_losing_mapping_fields() {
        let hit = SimpleHit {
            tid: 7,
            is_fw: true,
            pos: 42,
            mate_is_fw: false,
            mate_pos: 99,
            fragment_length: 120,
            score: 3.5,
            num_hits: 8,
            ..SimpleHit::default()
        };
        let converted = PseudoalignmentHit::from(hit);
        assert_eq!(converted.target_id, 7);
        assert!(converted.forward);
        assert_eq!(converted.position, 42);
        assert_eq!(converted.mate_position, Some(99));
        assert_eq!(converted.mate_forward, Some(false));
        assert_eq!(converted.fragment_length, Some(120));
        assert_eq!(converted.num_hits, 8);
        assert_eq!(converted.score, 3.5);
    }

    #[test]
    fn converts_missing_mate_sentinels_to_none() {
        let converted = PseudoalignmentHit::from(SimpleHit::default());
        assert_eq!(converted.mate_position, None);
        assert_eq!(converted.mate_forward, None);
        assert_eq!(converted.fragment_length, None);
    }

    #[test]
    fn sorts_hits_deterministically() {
        let mut hits = vec![
            PseudoalignmentHit::from(SimpleHit {
                tid: 2,
                pos: 10,
                ..SimpleHit::default()
            }),
            PseudoalignmentHit::from(SimpleHit {
                tid: 1,
                pos: 20,
                is_fw: false,
                ..SimpleHit::default()
            }),
            PseudoalignmentHit::from(SimpleHit {
                tid: 1,
                pos: 10,
                is_fw: true,
                ..SimpleHit::default()
            }),
        ];
        hits.sort_by(compare_hits);
        assert_eq!(
            hits.iter()
                .map(|hit| (hit.target_id, hit.forward, hit.position))
                .collect::<Vec<_>>(),
            vec![(1, true, 10), (1, false, 20), (2, false, 10)]
        );
    }

    #[test]
    fn validates_technical_sequences() {
        assert!(is_acgt(b"ACGTacgt"));
        assert!(!is_acgt(b""));
        assert!(!is_acgt(b"ACNT"));
    }

    #[test]
    fn validates_mapper_options_and_supported_k_values() {
        let mut options = PiscemOptions::default();
        options.max_hit_occ = 0;
        assert!(validate_options(&options)
            .unwrap_err()
            .to_string()
            .contains("max_hit_occ"));

        assert!(is_supported_k(3));
        assert!(is_supported_k(63));
        assert!(!is_supported_k(2));
        assert!(!is_supported_k(32));
        assert!(!is_supported_k(65));
    }

    #[test]
    fn maps_records_through_a_real_index_in_input_order() {
        let (_dir, index_prefix) = build_test_index();
        let options = PiscemOptions {
            num_threads: 2,
            load_ec: false,
            ..PiscemOptions::default()
        };
        let mapper = PiscemMapper::open(&index_prefix, options).unwrap();
        assert_eq!(mapper.k(), 5);
        assert_eq!(mapper.m(), 3);
        assert_eq!(mapper.num_targets(), 1);
        assert_eq!(mapper.target_name(0), Some("tx0"));
        assert_eq!(mapper.target_length(0), Some(20));
        assert_eq!(mapper.target_name(1), None);

        let batch = mapper.map_reads(vec![
            annotated_record(
                Some(Some(b"AAAA")),
                Some(b"CCCC"),
                Some(b"ACGTTGCACTGATCGTACGA"),
                None,
            ),
            annotated_record(
                Some(Some(b"AAAC")),
                Some(b"CCCG"),
                Some(b"TTTTTTTTTTTTTTTTTTTT"),
                None,
            ),
            annotated_record(
                Some(Some(b"AAAG")),
                Some(b"CCCT"),
                Some(b""),
                Some(b"ACGTTGCACTGATCGTACGA"),
            ),
        ]);

        assert_eq!(batch.records.len(), 3);
        assert_eq!(
            batch.records[0].barcode.as_deref(),
            Some(b"AAAA".as_slice())
        );
        assert_eq!(batch.records[0].status, PseudoalignmentStatus::Mapped);
        assert_eq!(batch.records[0].mapping_type, PseudoalignmentType::Single);
        assert!(batch.records[0].hits.iter().all(|hit| hit.target_id == 0));
        assert_eq!(
            batch.records[1].barcode.as_deref(),
            Some(b"AAAC".as_slice())
        );
        assert_eq!(batch.records[1].status, PseudoalignmentStatus::Unmapped);
        assert_eq!(
            batch.records[2].barcode.as_deref(),
            Some(b"AAAG".as_slice())
        );
        assert_eq!(batch.records[2].status, PseudoalignmentStatus::Mapped);
        assert_eq!(batch.records[2].mapping_type, PseudoalignmentType::Single);
        assert_eq!(batch.stats.input_records, 3);
        assert_eq!(batch.stats.submitted_records, 3);
        assert_eq!(batch.stats.mapped_records, 2);
        assert_eq!(batch.stats.unmapped_records, 1);
    }

    #[test]
    fn structural_mapping_uses_the_same_real_index() {
        let (_dir, index_prefix) = build_test_index();
        let options = PiscemOptions {
            num_threads: 1,
            load_ec: false,
            structural_constraints: true,
            ..PiscemOptions::default()
        };
        let mapper = PiscemMapper::open(index_prefix, options).unwrap();
        let batch = mapper.map_reads(vec![annotated_record(
            Some(Some(b"AAAA")),
            Some(b"CCCC"),
            Some(b"ACGTTGCACTGATCGTACGA"),
            None,
        )]);
        assert_eq!(batch.records[0].status, PseudoalignmentStatus::Mapped);
    }

    #[test]
    fn ec_loading_is_independent_of_the_cardinality_threshold() {
        let (_dir, index_prefix) = build_test_index();
        let options = PiscemOptions {
            num_threads: 1,
            load_ec: true,
            max_ec_card: 0,
            ..PiscemOptions::default()
        };
        let mapper = PiscemMapper::open(index_prefix, options).unwrap();
        assert!(mapper.has_ec_table());
    }

    #[test]
    fn reports_reads_rejected_by_a_loaded_poison_table() {
        let (_dir, index_prefix) = build_test_index();
        let table = PoisonTable::build_from_occs(vec![LabeledPoisonOcc {
            // Canonical 2-bit encoding of ACGTT (reverse complement AACGT).
            canonical_kmer: CanonicalKmer::new(27),
            unitig_id: 0,
            unitig_pos: 0,
        }])
        .unwrap();
        let file = std::fs::File::create(index_prefix.with_extension("poison")).unwrap();
        table.save(&mut BufWriter::new(file)).unwrap();

        let mapper = PiscemMapper::open(
            index_prefix,
            PiscemOptions {
                num_threads: 1,
                load_ec: false,
                load_poison: true,
                ..PiscemOptions::default()
            },
        )
        .unwrap();
        let batch = mapper.map_reads(vec![annotated_record(
            Some(Some(b"AAAA")),
            Some(b"CCCC"),
            Some(b"ACGTTGCACTGATCGTACGA"),
            None,
        )]);
        assert_eq!(batch.records[0].status, PseudoalignmentStatus::Poisoned);
        assert_eq!(batch.stats.poisoned_records, 1);
        assert_eq!(batch.stats.submitted_records, 1);
    }

    #[test]
    fn rejects_missing_and_invalid_metadata_before_mapping() {
        let missing_barcode =
            prepare_record(annotated_record(None, Some(b"CCCC"), Some(b"ACGT"), None)).unwrap_err();
        assert_eq!(
            missing_barcode.status,
            PseudoalignmentStatus::MissingBarcode
        );

        let uncorrected = prepare_record(annotated_record(
            Some(None),
            Some(b"CCCC"),
            Some(b"ACGT"),
            None,
        ))
        .unwrap_err();
        assert_eq!(
            uncorrected.status,
            PseudoalignmentStatus::MissingCorrectedBarcode
        );

        let invalid_umi = prepare_record(annotated_record(
            Some(Some(b"AAAA")),
            Some(b"CCNC"),
            Some(b"ACGT"),
            None,
        ))
        .unwrap_err();
        assert_eq!(invalid_umi.status, PseudoalignmentStatus::InvalidUmi);
    }

    #[test]
    fn retains_a_nonempty_second_read_when_the_first_is_empty() {
        let prepared = prepare_record(annotated_record(
            Some(Some(b"AAAA")),
            Some(b"CCCC"),
            Some(b""),
            Some(b"ACGT"),
        ))
        .unwrap();
        assert!(prepared.read1.is_none());
        assert_eq!(prepared.read2.unwrap().sequence(), b"ACGT");
    }

    #[test]
    fn aggregates_batch_statuses() {
        let records = vec![
            Pseudoalignment::rejected(None, None, PseudoalignmentStatus::MissingBarcode),
            Pseudoalignment::rejected(
                Some(b"AAAA".to_vec()),
                Some(b"CCCC".to_vec()),
                PseudoalignmentStatus::Unmapped,
            ),
            Pseudoalignment {
                barcode: Some(b"AAAA".to_vec()),
                umi: Some(b"CCCC".to_vec()),
                status: PseudoalignmentStatus::Mapped,
                mapping_type: PseudoalignmentType::Single,
                hits: Vec::new(),
            },
        ];
        let stats = PiscemStats::from_records(&records);
        assert_eq!(stats.input_records, 3);
        assert_eq!(stats.submitted_records, 2);
        assert_eq!(stats.mapped_records, 1);
        assert_eq!(stats.unmapped_records, 1);
        assert_eq!(stats.missing_barcode, 1);
    }
}
