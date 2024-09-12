// In Illumina sequencing, there are typically two types of duplicates:
// PCR duplicates and optical duplicates. Optical duplicates are sequences from
// one cluster in fact but identified by software from multiple adjacent clusters.
// They can be identified without alignment. You just check sequence and
// the coordinates on the image.
// PCR duplicates are usually identified after alignment.
// Dedup typically works by identifying read pairs having identical 5'-end
// coordinates (3'-end coordinates are not considered). The chance of
// two 5'-end coordinates being identical depends on the coverage, but usually small.
// You can calculate the theoretical false dedup rate using my formula
// (0.28*m/s/L, where m is the number of read pairs, s is the standard deviation
// of insert size distribution and L is the length of the genome;
// sharp insert size distribution (small s) increases false dedup rate unfortunately).
// If the empirical duplicate rate is much higher than this 0.28*m/s/L, they are true PCR duplicates.
// Dedup for single-end reads has high false positive rate given deep data.
// PCR duplicates do not necessarily have the same sequence OR the same length.
// The 5' end of the read is the most reliable position because
// if it's a PCR duplicate, it will be guaranteed to be the same at that end
// but not at the 3' end.

use indexmap::IndexMap;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::Record;
use noodles::sam::{alignment::record::data::field::{Tag, Value}, Header};
use bed_utils::bed::Strand;
use itertools::Itertools;
use std::hash::Hash;
use anyhow::{Result, Context};
use serde::{Serialize, Deserialize};

use crate::fragment::Fragment;

#[derive(Serialize, Deserialize, Debug)]
pub struct AlignmentInfo {
    pub(crate) this: SingleAlignment,
    other: Option<SingleAlignment>,
}

impl AlignmentInfo {
    pub fn from_read<R: Record>(rec: &R, header: &Header) -> Result<Self> {
        Ok(Self {
            this: SingleAlignment::new(rec, header)?,
            other: None,
        })
    }

    pub fn from_read_pair<R: Record>(records: &(R, R), header: &Header) -> Result<Self> {
        let this = SingleAlignment::new(&records.0, header)?;
        let other = SingleAlignment::new(&records.1, header)?;
        Ok(Self { this, other: Some(other) })
    }

    pub fn to_fragment(&self, header: &Header, count: usize) -> Option<Fragment> {
        let ref_id1: usize = self.this.reference_sequence_id;
        let chrom = header.reference_sequences().get_index(ref_id1).unwrap().0.to_string();
        let barcode = Some(self.this.barcode.clone());
        if self.other.is_none() {
            Some(Fragment {
                chrom,
                start: self.this.alignment_start as u64 - 1,
                end: self.this.alignment_end as u64,
                barcode,
                count: count.try_into().unwrap(),
                strand: Some(if self.this.is_reverse_complemented() {
                    Strand::Reverse
                } else {
                    Strand::Forward
                }),
            })
        } else {
            let other = self.other.as_ref().unwrap();
            let ref_id2: usize = other.reference_sequence_id.try_into().unwrap();
            if ref_id1 != ref_id2 { return None; }

            let rec1_5p = self.this.alignment_5p();
            let rec2_5p = other.alignment_5p();
            let (start, end) = if rec1_5p < rec2_5p {
                (rec1_5p, rec2_5p)
            } else {
                (rec2_5p, rec1_5p)
            };

            Some(Fragment {
                chrom,
                start: start as u64 - 1,
                end: end as u64,
                barcode,
                count: count.try_into().unwrap(),
                strand: None,
            })
        }
    }
}

pub fn remove_duplicates<I>(data: I, header: &Header) -> Vec<Fragment>
where
    I: IntoIterator<Item = AlignmentInfo>,
{
    let mut fingerprints = IndexMap::new();
    let collected = data.into_iter().filter(|ali| {
        let fingerprint = FingerPrint::from(ali);
        let is_duplicate = fingerprints.contains_key(&fingerprint);
        if is_duplicate {
            fingerprints.entry(fingerprint).and_modify(|x| *x += 1);
        } else {
            fingerprints.insert(fingerprint, 1);
        }
        !is_duplicate
    }).collect::<Vec<_>>();
    collected.into_iter().zip_eq(fingerprints.values()).flat_map(|(ali, count)|
        ali.to_fragment(header, *count)
    ).collect()
}

// Library type    orientation   Vizualization according to first strand
// FF_firststrand  matching      3' <==2==----<==1== 5'
//                               5' ---------------- 3'
//
// FF_secondstrand matching      3' ---------------- 5'
//                               5' ==1==>----==2==> 3'
//
// RR_firststrand  matching      3' <==1==----<==2== 5'
//                               5' ---------------- 3'
//
// RR_secondstrand matching      3' ---------------- 5'
//                               5' ==2==>----==1==> 3'
//
// FR_firststrand  inward        3' ----------<==1== 5'
//                               5' ==2==>---------- 3'
//
// FR_secondstrand inward        3' ----------<==2== 5'
//                               5' ==1==>---------- 3'
//
// RF_firststrand  outward       3' <==2==---------- 5'
//                               5' ----------==1==> 3'
//
// RF_secondstrand outward       3' <==1==---------- 5'
//                               5' ----------==2==> 3'
#[derive(Eq, PartialEq, Debug, Hash, Serialize, Deserialize)]
pub enum Orientation { FR, FF, RR, RF }

/// Minimal information about an alignment extracted from the BAM record.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct SingleAlignment {
    pub name: String,
    pub reference_sequence_id: usize,
    pub alignment_start: u32,
    pub alignment_end: u32,
    pub unclipped_start: u32,
    pub unclipped_end: u32,
    pub barcode: String,
    pub umi: Option<String>,
    pub flags: u16,
}

fn get_barcode<R: Record>(rec: &R) -> Result<Option<String>> {
    Ok(rec.data().get(&Tag::CELL_BARCODE_ID).transpose()?.and_then(|x| match x {
        Value::String(barcode) => Some(barcode.to_string()),
        _ => None,
    }))
}

fn get_umi<R: Record>(rec: &R) -> Result<Option<String>> {
    Ok(rec.data().get(&Tag::UMI_ID).transpose()?.and_then(|x| match x {
        Value::String(umi) => Some(umi.to_string()),
        _ => None,
    }))
}

impl SingleAlignment {
    pub fn new<R: Record>(
        rec: &R,
        header: &Header,
    ) -> Result<Self> {
        let cigar = rec.cigar();
        let start: usize = rec.alignment_start().unwrap().unwrap().try_into()?;
        let alignment_start: u32 = start.try_into()?;
        let alignment_span: u32 = cigar.alignment_span()?.try_into()?;
        let alignment_end = alignment_start + alignment_span - 1;
        let clip_groups = cigar.iter().map(Result::unwrap).chunk_by(|op| {
            let kind = op.kind();
            kind == Kind::HardClip || kind == Kind::SoftClip
        });
        let mut clips = clip_groups.into_iter();
        let clipped_start: u32 = clips.next().map_or(0, |(is_clip, x)| if is_clip {
            x.map(|x| x.len() as u32).sum()
        } else {
            0
        });
        let clipped_end: u32 = clips.last().map_or(0, |(is_clip, x)| if is_clip {
            x.map(|x| x.len() as u32).sum()
        } else {
            0
        });
        Ok(Self {
            name: std::str::from_utf8(rec.name().context("no read name")?)?.to_string(),
            reference_sequence_id: rec.reference_sequence_id(header).context("no reference sequence id")??,
            alignment_start,
            alignment_end,
            unclipped_start: alignment_start - clipped_start,
            unclipped_end: alignment_end + clipped_end,
            barcode: get_barcode(rec)?.with_context(|| "no barcode")?,
            umi: get_umi(rec)?,
            flags: rec.flags()?.bits(),
        })
    }

    pub fn is_reverse_complemented(&self) -> bool {
        Flags::from_bits_retain(self.flags).is_reverse_complemented()
    }

    pub fn is_first_segment(&self) -> bool {
        Flags::from_bits_retain(self.flags).is_first_segment()
    }

    pub fn coord_5p(&self) -> u32 {
        if self.is_reverse_complemented() {
            self.unclipped_end
        } else {
            self.unclipped_start
        }
    }

    pub fn alignment_5p(&self) -> u32 {
        if self.is_reverse_complemented() {
            self.alignment_end
        } else {
            self.alignment_start
        }
    }
}

/// Reads are considered duplicates if and only if they have the same fingerprint.
#[derive(Eq, PartialEq, Debug, Hash)]
enum FingerPrint {
    Single {
        reference_id: usize,
        coord_5p: u32,
        orientation: Orientation,
        barcode: Option<String>,
    },
    Paired {
        left_reference_id: usize,
        right_reference_id: usize,
        left_coord_5p: u32,
        right_coord_5p: u32,
        orientation: Orientation,
        barcode: Option<String>,
    },
}

impl From<&AlignmentInfo> for FingerPrint {
    fn from(ali: &AlignmentInfo) -> FingerPrint {
        if ali.other.is_none() {
            let orientation = if ali.this.is_reverse_complemented() {
                Orientation::RR
            } else {
                Orientation::FF
            };
            FingerPrint::Single {
                reference_id: ali.this.reference_sequence_id,
                coord_5p: ali.this.coord_5p(),
                orientation,
                barcode: ali.this.umi.clone(),
            }
        } else {
            let mut this = &ali.this;
            let mut other = ali.other.as_ref().unwrap();
            if this.umi != other.umi { panic!("UMI mismatch"); }

            let this_is_leftmost = if this.reference_sequence_id == other.reference_sequence_id {
                if this.coord_5p() < other.coord_5p() { true } else { false }
            } else {
                this.reference_sequence_id < other.reference_sequence_id
            };
            if !this_is_leftmost {
                this = ali.other.as_ref().unwrap();
                other = &ali.this;
            }
            let orientation = if this.is_reverse_complemented() == other.is_reverse_complemented() {
                if this.is_reverse_complemented() {
                    if this.is_first_segment() { Orientation::RR } else { Orientation::FF }
                } else {
                    if this.is_first_segment() { Orientation::FF } else { Orientation::RR }
                }
            } else {
                if this.is_reverse_complemented() { Orientation::RF } else { Orientation::FR }
            };
            FingerPrint::Paired {
                left_reference_id: this.reference_sequence_id,
                right_reference_id: other.reference_sequence_id,
                left_coord_5p: this.coord_5p(),
                right_coord_5p: other.coord_5p(),
                orientation,
                barcode: this.umi.clone(),
            }
        }
    }
}