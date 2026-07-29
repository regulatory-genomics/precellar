//! In-memory pseudoalignment backends.

mod piscem;

pub use piscem::{
    build_splici_index, build_transcriptome_index, PiscemBatch, PiscemIndexBuildOptions,
    PiscemMapper, PiscemOptions, PiscemStats, Pseudoalignment, PseudoalignmentHit,
    PseudoalignmentStatus, PseudoalignmentType, SkippingStrategy, SpliciIndex,
    SpliciIndexBuildOptions,
};
