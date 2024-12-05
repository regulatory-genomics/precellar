use std::path::Path;

use noodles::sam;
use precellar::align::{Aligner, BurrowsWheelerAligner, StarAligner};
use seqspec::Modality;

pub enum AlignerType {
    STAR(StarAligner),
    BWA(BurrowsWheelerAligner),
}

impl AlignerType {
    pub fn from_name<P: AsRef<Path>>(name: &str, path: P) -> Self {
        match name.to_lowercase().as_str() {
            "star" => AlignerType::STAR(StarAligner::from_path(path)),
            "bwa" => AlignerType::BWA(BurrowsWheelerAligner::from_path(path)),
            _ => unimplemented!(),
        }
    }

    pub fn from_modality<P: AsRef<Path>>(modality: Modality, path: P) -> Self {
        match modality {
            Modality::RNA => AlignerType::STAR(StarAligner::from_path(path)),
            Modality::ATAC => AlignerType::BWA(BurrowsWheelerAligner::from_path(path)),
            _ => unimplemented!(),
        }
    }

    pub fn header(&self) -> sam::Header {
        match self {
            AlignerType::STAR(aligner) => aligner.header(),
            AlignerType::BWA(aligner) => aligner.header(),
        }
    }
}
