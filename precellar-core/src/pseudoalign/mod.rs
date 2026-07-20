//! In-memory pseudoalignment backends.

mod piscem;

pub use piscem::{
    PiscemBatch, PiscemMapper, PiscemOptions, PiscemStats, Pseudoalignment, PseudoalignmentHit,
    PseudoalignmentStatus, PseudoalignmentType, SkippingStrategy,
};
