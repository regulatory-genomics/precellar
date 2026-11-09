//! Coarse stage timing for profiling a run.
//!
//! Timing is off unless the `PRECELLAR_TRACE` environment variable names an
//! output path. While enabled, each stage accumulates four values: the time its
//! first span started, the time its last span ended, the total active time and
//! the number of spans. Nothing is written during the run; the caller flushes
//! once at the end with [`flush`].
//!
//! Times are seconds since the trace origin, which is the moment [`init`] runs
//! (Python module import). The output header also records the origin as a Unix
//! timestamp so an external CPU sampler can be aligned on it.
//!
//! Only batch-level work is timed. Per-record work such as `qc.update` or
//! pushing into the external sorter is deliberately left out: its cost per call
//! is far below timer resolution, and summing millions of spans would add
//! serial time to the very phases being profiled. Such work shows up as a gap
//! between the surrounding stages.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::OnceLock;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Number of stages; keep in sync with [`Stage`].
const NUM_STAGES: usize = 4;

/// Stages reported by the trace, in output order.
///
/// Stages are deliberately coarse. `prepare` covers everything before the
/// first batch: assay parsing, read-length inference, reference index loading,
/// whitelist counting and FASTQ reader setup. `quantify` covers everything
/// after the last batch is annotated: external-sort collection, per-cell
/// counting and writing the output.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    /// Up to the first FASTQ batch.
    Prepare = 0,
    /// Reading, decompression and barcode annotation of one batch.
    FastqAnnotation = 1,
    /// Genome alignment of one batch, including alignment QC.
    Alignment = 2,
    /// Transcriptome annotation of one batch.
    TxAlign = 3,
}

impl Stage {
    /// Name used in the trace output.
    pub fn name(self) -> &'static str {
        match self {
            Stage::Prepare => "prepare",
            Stage::FastqAnnotation => "fastq_annotation",
            Stage::Alignment => "alignment",
            Stage::TxAlign => "tx_align",
        }
    }
}

/// The tail of a run, which is not a [`Stage`] because it spans the sink rather
/// than one iterator. Reported after `tx_align`.
const QUANTIFY: usize = 4;

/// The stages written to the trace, in output order.
fn stage_names() -> [&'static str; NUM_STAGES + 1] {
    [
        Stage::Prepare.name(),
        Stage::FastqAnnotation.name(),
        Stage::Alignment.name(),
        Stage::TxAlign.name(),
        "quantify",
    ]
}

struct Origin {
    instant: Instant,
    unix: f64,
}

fn origin() -> &'static Origin {
    static ORIGIN: OnceLock<Origin> = OnceLock::new();
    ORIGIN.get_or_init(|| Origin {
        unix: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs_f64())
            .unwrap_or(0.0),
        instant: Instant::now(),
    })
}

fn enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

static ENABLED: AtomicBool = AtomicBool::new(false);
static PREPARE_FINISHED: AtomicBool = AtomicBool::new(false);
static TRACE_PATH: OnceLock<String> = OnceLock::new();

/// `u64::MAX` when the stage has not run yet.
static FIRST_START_NS: [AtomicU64; NUM_STAGES + 1] =
    [const { AtomicU64::new(u64::MAX) }; NUM_STAGES + 1];
/// Zero until the stage has finished at least one span.
static LAST_END_NS: [AtomicU64; NUM_STAGES + 1] = [const { AtomicU64::new(0) }; NUM_STAGES + 1];
static ACTIVE_NS: [AtomicU64; NUM_STAGES + 1] = [const { AtomicU64::new(0) }; NUM_STAGES + 1];
static COUNT: [AtomicU64; NUM_STAGES + 1] = [const { AtomicU64::new(0) }; NUM_STAGES + 1];

/// Starts the trace clock and arms the timer if `PRECELLAR_TRACE` is set.
///
/// Safe to call more than once; only the first call has an effect.
pub fn init() {
    let _ = origin();
    if ENABLED.load(Ordering::Relaxed) {
        return;
    }
    let Ok(path) = std::env::var("PRECELLAR_TRACE") else {
        return;
    };
    if path.is_empty() {
        return;
    }
    let _ = TRACE_PATH.set(path);
    ENABLED.store(true, Ordering::Relaxed);
    // `prepare` starts at the origin by definition.
    FIRST_START_NS[Stage::Prepare as usize].store(0, Ordering::Relaxed);
}

fn elapsed_ns() -> u64 {
    origin().instant.elapsed().as_nanos() as u64
}

/// One timed span; the stage is charged when the guard is dropped.
pub struct StageSpan {
    index: usize,
    start: Option<Instant>,
}

impl Drop for StageSpan {
    fn drop(&mut self) {
        let Some(start) = self.start else { return };
        let end_ns = elapsed_ns();
        let active_ns = start.elapsed().as_nanos() as u64;
        let start_ns = end_ns.saturating_sub(active_ns);
        let count = COUNT[self.index].fetch_add(1, Ordering::Relaxed);
        if count == 0 {
            FIRST_START_NS[self.index].store(start_ns, Ordering::Relaxed);
        }
        LAST_END_NS[self.index].fetch_max(end_ns, Ordering::Relaxed);
        ACTIVE_NS[self.index].fetch_add(active_ns, Ordering::Relaxed);
    }
}

/// Times one span of `stage` until the returned guard is dropped.
pub fn stage(stage: Stage) -> StageSpan {
    StageSpan {
        index: stage as usize,
        start: enabled().then(Instant::now),
    }
}

/// Ends the `prepare` stage at the first FASTQ batch.
///
/// Called from the batch iterator; only the first call has any effect.
pub fn finish_prepare() {
    if !enabled() || PREPARE_FINISHED.swap(true, Ordering::Relaxed) {
        return;
    }
    let index = Stage::Prepare as usize;
    let now = elapsed_ns();
    COUNT[index].store(1, Ordering::Relaxed);
    LAST_END_NS[index].store(now, Ordering::Relaxed);
    ACTIVE_NS[index].store(now, Ordering::Relaxed);
    FIRST_START_NS[index].store(0, Ordering::Relaxed);
}

/// `u64::MAX` until the alignment stream has ended.
static QUANTIFY_START_NS: AtomicU64 = AtomicU64::new(u64::MAX);

/// Marks the start of the sink tail. Only the first call has an effect.
pub fn begin_quantify() {
    if !enabled() {
        return;
    }
    let now = elapsed_ns();
    let _ = QUANTIFY_START_NS.compare_exchange(u64::MAX, now, Ordering::Relaxed, Ordering::Relaxed);
}

/// Marks the end of the sink tail and records it as a single span.
pub fn end_quantify() {
    if !enabled() {
        return;
    }
    let start = QUANTIFY_START_NS.load(Ordering::Relaxed);
    let end = elapsed_ns();
    if start == u64::MAX || start > end {
        return;
    }
    FIRST_START_NS[QUANTIFY].store(start, Ordering::Relaxed);
    LAST_END_NS[QUANTIFY].store(end, Ordering::Relaxed);
    ACTIVE_NS[QUANTIFY].store(end - start, Ordering::Relaxed);
    COUNT[QUANTIFY].store(1, Ordering::Relaxed);
}

/// Iterator adapter that starts the `quantify` stage when its input runs out.
///
/// Quantification pulls the alignment stream lazily from inside the external
/// sorter, so no call site sees the moment the last batch has been annotated.
/// The first `None` from the stream marks that moment.
pub struct TraceTail<I> {
    inner: I,
}

impl<I: Iterator> Iterator for TraceTail<I> {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        let item = self.inner.next();
        if item.is_none() {
            begin_quantify();
        }
        item
    }
}

/// Wraps the alignment stream so the `quantify` stage starts when it ends.
pub fn trace_tail<I: Iterator>(inner: I) -> TraceTail<I> {
    TraceTail { inner }
}

fn seconds(ns: u64) -> f64 {
    ns as f64 / 1e9
}

fn summary() -> String {
    let names = stage_names();
    let mut parts = Vec::new();
    for index in 0..=NUM_STAGES {
        let count = COUNT[index].load(Ordering::Relaxed);
        if count == 0 {
            continue;
        }
        parts.push(format!(
            "{}: {:.2}-{:.2}s active {:.2}s in {} span(s)",
            names[index],
            seconds(FIRST_START_NS[index].load(Ordering::Relaxed)),
            seconds(LAST_END_NS[index].load(Ordering::Relaxed)),
            seconds(ACTIVE_NS[index].load(Ordering::Relaxed)),
            count
        ));
    }
    parts.join("; ")
}

/// Writes the trace and logs the summary. Does nothing when timing is off.
pub fn flush() {
    if !enabled() {
        return;
    }
    let summary = summary();
    log::info!("Stage timing: {summary}");

    let Some(path) = TRACE_PATH.get() else {
        return;
    };
    let names = stage_names();
    let mut text = format!(
        "# precellar stage timing\n# origin_unix {:.6}\nstage\tfirst_start_s\tlast_end_s\tactive_s\tcount\n",
        origin().unix
    );
    for index in 0..=NUM_STAGES {
        let count = COUNT[index].load(Ordering::Relaxed);
        if count == 0 {
            continue;
        }
        text.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{:.6}\t{}\n",
            names[index],
            seconds(FIRST_START_NS[index].load(Ordering::Relaxed)),
            seconds(LAST_END_NS[index].load(Ordering::Relaxed)),
            seconds(ACTIVE_NS[index].load(Ordering::Relaxed)),
            count
        ));
    }
    if let Err(error) = std::fs::write(path, text) {
        log::warn!("Could not write the stage trace to {path}: {error}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_by_default() {
        // Tests run without PRECELLAR_TRACE, so no stage is recorded.
        let _span = stage(Stage::Alignment);
        assert!(!enabled());
    }

    #[test]
    fn stage_indices_match_output_order() {
        let names = stage_names();
        assert_eq!(names[Stage::Prepare as usize], "prepare");
        assert_eq!(names[Stage::FastqAnnotation as usize], "fastq_annotation");
        assert_eq!(names[Stage::Alignment as usize], "alignment");
        assert_eq!(names[Stage::TxAlign as usize], "tx_align");
        assert_eq!(names[QUANTIFY], "quantify");
    }
}
