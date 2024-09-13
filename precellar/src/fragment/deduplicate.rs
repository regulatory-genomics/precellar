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

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::Record;
use noodles::sam::{alignment::record::data::field::{Tag, Value}, Header};
use bed_utils::bed::Strand;
use itertools::Itertools;
use std::collections::HashMap;
use std::hash::Hash;
use anyhow::{Result, Context};
use serde::{Serialize, Deserialize};

use crate::fragment::Fragment;

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
pub enum Orientation { F, R, FR, FF, RR, RF }

/// Reads are considered duplicates if and only if they have the same fingerprint.
#[derive(Eq, PartialEq, Debug, Hash)]
pub(crate) enum FingerPrint {
    Single {
        reference_id: u16,
        coord_5p: u32,
        orientation: Orientation,
        umi: Option<String>,
    },
    Paired {
        reference_id: u16,
        left_coord_5p: u32,
        right_coord_5p: u32,
        orientation: Orientation,
        umi: Option<String>,
    },
}

/// Minimal information about an alignment extracted from the BAM record.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct AlignmentMini {
    alignment_start: u32,
    alignment_end: u32,
    unclipped_start: u32,
    unclipped_end: u32,
    flags: u16,
}

impl AlignmentMini {
    fn new<R: Record>(rec: &R) -> Result<Self> {
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
            alignment_start,
            alignment_end,
            unclipped_start: alignment_start - clipped_start,
            unclipped_end: alignment_end + clipped_end,
            flags: rec.flags()?.bits(),
        })
    }

    pub fn is_reverse_complemented(&self) -> bool {
        Flags::from_bits_retain(self.flags).is_reverse_complemented()
    }

    pub fn is_first_segment(&self) -> bool {
        Flags::from_bits_retain(self.flags).is_first_segment()
    }

    pub fn unclipped_5p(&self) -> u32 {
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
 
#[derive(Serialize, Deserialize, Debug)]
pub struct AlignmentInfo {
    name: String,
    reference_sequence_id: u16,
    reference_sequence: String,
    read1: AlignmentMini,
    read2: Option<AlignmentMini>,
    pub(crate) barcode: String,
    umi: Option<String>,
}

impl AlignmentInfo {
    pub fn from_read<R: Record>(rec: &R, header: &Header) -> Result<Option<Self>> {
        let barcode = get_barcode(rec)?;
        if barcode.is_none() { return Ok(None); }

        let name = rec.name().unwrap().to_string();
        let reference_sequence_id: u16 = rec.reference_sequence_id(header).context("no reference sequence id")??.try_into()?;
        let reference_sequence = header.reference_sequences().get_index(reference_sequence_id as usize).unwrap().0.to_string();
        let umi = get_umi(rec)?;
        Ok(Some(Self {
            name,
            reference_sequence_id,
            reference_sequence,
            read1: AlignmentMini::new(rec)?,
            read2: None,
            barcode: barcode.unwrap(),
            umi,
        }))
    }

    pub fn from_read_pair<R: Record>(records: (&R, &R), header: &Header) -> Result<Option<Self>> {
        let rec1 = records.0;
        let rec2 = records.1;
        let barcode1 = get_barcode(rec1)?;
        let barcode2 = get_barcode(rec2)?;
        if barcode1 != barcode2 || barcode1.is_none() { return Ok(None); }

        let name1 = rec1.name().unwrap();
        let name2 = rec2.name().unwrap();
        assert!(name1 == name2, "Read names do not match");

        let reference_sequence_id1: u16 = rec1.reference_sequence_id(header).context("no reference sequence id")??.try_into()?;
        let reference_sequence_id2: u16 = rec2.reference_sequence_id(header).context("no reference sequence id")??.try_into()?;
        if reference_sequence_id1 != reference_sequence_id2 { return Ok(None); }
        Ok(Some(Self {
            name: name1.to_string(),
            reference_sequence_id: reference_sequence_id1,
            reference_sequence: header.reference_sequences()
                .get_index(reference_sequence_id1 as usize).unwrap().0.to_string(),
            read1: AlignmentMini::new(rec1)?,
            read2: Some(AlignmentMini::new(rec2)?),
            barcode: barcode1.unwrap(),
            umi: get_umi(rec1)?,
        }))
    }
}

impl From<(AlignmentInfo, usize)> for Fragment {
    fn from(value: (AlignmentInfo, usize)) -> Fragment {
        if value.0.read2.is_none() {
            Fragment {
                chrom: value.0.reference_sequence.clone(),
                start: value.0.read1.alignment_start as u64 - 1,
                end: value.0.read1.alignment_end as u64,
                barcode: Some(value.0.barcode.clone()),
                count: value.1.try_into().unwrap(),
                strand: Some(if value.0.read1.is_reverse_complemented() {
                    Strand::Reverse
                } else {
                    Strand::Forward
                }),
            }
        } else {
            let rec1_5p = value.0.read1.alignment_5p();
            let rec2_5p = value.0.read2.unwrap().alignment_5p();
            let (start, end) = if rec1_5p < rec2_5p {
                (rec1_5p, rec2_5p)
            } else {
                (rec2_5p, rec1_5p)
            };

            Fragment {
                chrom: value.0.reference_sequence.clone(),
                start: start as u64 - 1,
                end: end as u64,
                barcode: Some(value.0.barcode.clone()),
                count: value.1.try_into().unwrap(),
                strand: None,
            }
        }
    }
}

impl From<&AlignmentInfo> for FingerPrint {
    fn from(ali: &AlignmentInfo) -> FingerPrint {
        if ali.read2.is_none() {
            let orientation = if ali.read1.is_reverse_complemented() {
                Orientation::R
            } else {
                Orientation::F
            };
            FingerPrint::Single {
                reference_id: ali.reference_sequence_id,
                coord_5p: ali.read1.unclipped_5p(),
                orientation,
                umi: ali.umi.clone(),
            }
        } else {
            let this;
            let other;
            if ali.read1.unclipped_5p() < ali.read2.as_ref().unwrap().unclipped_5p() {
                this = &ali.read1;
                other = ali.read2.as_ref().unwrap();
            } else {
                other = &ali.read1;
                this = ali.read2.as_ref().unwrap();
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
                reference_id: ali.reference_sequence_id,
                left_coord_5p: this.unclipped_5p(),
                right_coord_5p: other.unclipped_5p(),
                orientation,
                umi: ali.umi.clone(),
            }
        }
    }
}


pub(crate) fn remove_duplicates<I>(data: I) -> HashMap<FingerPrint, Fragment>
where
    I: IntoIterator<Item = AlignmentInfo>,
{
    let mut result: HashMap<FingerPrint, Fragment> = HashMap::new();
    data.into_iter().for_each(|ali| {
        let fingerprint = FingerPrint::from(&ali);
        result.entry(fingerprint).and_modify(|x| x.count +=1)
            .or_insert(Fragment::from((ali, 1)));
    });
    result
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

