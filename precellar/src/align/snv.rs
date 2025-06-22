use anyhow::{bail, Result};
use bincode::{Decode, Encode};
use bstr::ByteSlice;
use lexical::parse_partial;
use noodles::sam::alignment::{
    record::{cigar::op::Kind, data::field::Tag, Cigar},
    record_buf::data::field::Value,
    RecordBuf,
};
use std::fmt::Write;

#[derive(Encode, Decode, Debug, Clone, PartialEq, Eq)]
pub struct SNV {
    pub relative_position: u32, // 0-based position of the mutation in the read
    pub ty: Mutation,
}

impl PartialOrd for SNV {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.relative_position.cmp(&other.relative_position))
    }
}

impl Ord for SNV {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.relative_position.cmp(&other.relative_position)
    }
}

/* This function formats a vector of position sorted SNVs into a string representation.
  The format is as follows:
  - The first number is the offset from the start of the read.
  - Each SNV is represented by its relative position and type:
    - Substitutions are represented by the base character.
    - Deletions are represented by '^' followed by '-' for each deleted base.
    - Insertions are represented by '+' followed by the inserted bases.
*/
pub fn format_snvs(snv: &[SNV]) -> Option<String> {
    if snv.is_empty() {
        return None;
    }

    let mut output = String::new();
    let mut offset = 0;
    snv.iter().for_each(|s| {
        if s.relative_position < offset {
            panic!(
                "SNV positions are not sorted or have gaps: {} < {}. {:?}",
                s.relative_position, offset, snv
            );
        }
        let n = s.relative_position - offset;
        write!(output, "{}", n).unwrap();
        offset = s.relative_position;

        match &s.ty {
            Mutation::Substitution(base) => {
                write!(output, "{}", *base as char).unwrap();
                offset += 1; // Increment offset for substitution
            }
            Mutation::Deletion(len) => {
                write!(output, "^").unwrap();
                for _ in 0..*len {
                    write!(output, "-").unwrap(); // Use '-' to represent deleted bases
                }
            }
            Mutation::Insertion(bases) => {
                let bases_str = String::from_utf8(bases.clone()).unwrap();
                write!(output, "+{}", bases_str).unwrap();
                offset += bases.len() as u32; // Increment offset for insertion
            }
        }
    });
    Some(output)
}

#[derive(Encode, Decode, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Mutation {
    Substitution(u8),   // The alternative base (A, T, C, G)
    Deletion(u16),      // The length of the deletion
    Insertion(Vec<u8>), // the inserted bases
}

impl Mutation {
    pub fn len(&self) -> usize {
        match self {
            Mutation::Substitution(_) => 1,
            Mutation::Deletion(_) => 0,
            Mutation::Insertion(bases) => bases.len(),
        }
    }

    pub fn is_deletion(&self) -> bool {
        matches!(self, Mutation::Deletion(_))
    }

    pub fn is_insertion(&self) -> bool {
        matches!(self, Mutation::Insertion(_))
    }

    pub fn is_substitution(&self) -> bool {
        matches!(self, Mutation::Substitution(_))
    }
}

impl From<&Mutation> for String {
    fn from(m: &Mutation) -> Self {
        match m {
            Mutation::Substitution(base) => String::from_utf8(vec![*base]).unwrap(),
            Mutation::Deletion(len) => format!("D{}", len),
            Mutation::Insertion(bases) => format!("I{}", String::from_utf8(bases.clone()).unwrap()),
        }
    }
}

pub fn get_snv(rec: &RecordBuf) -> Result<Vec<SNV>> {
    let mut snv = Vec::new();
    let md = match rec.data().get(&Tag::MISMATCHED_POSITIONS).unwrap() {
        Value::String(s) => s,
        _ => {
            bail!(
                "MD tag not found or not a string in record: {}",
                rec.name().unwrap()
            );
        }
    }
    .as_bytes();
    let mut md: Vec<_> = parse_md_tag_str(md)
        .flat_map(|x| match x.unwrap() {
            MDKind::Match(m) => {
                if m == 0 {
                    None
                } else {
                    Some(MDKind::Match(m))
                }
            }
            k => Some(k),
        })
        .collect();
    md.reverse();
    let mut n_soft_clip = 0;
    rec.cigar()
        .iter()
        .map(Result::unwrap)
        .take_while(|op| {
            let kind = op.kind();
            kind == Kind::HardClip || kind == Kind::SoftClip
        })
        .for_each(|op| {
            if op.kind() == Kind::SoftClip {
                n_soft_clip += op.len();
            }
        });

    let mut idx = 0;
    rec.cigar()
        .iter()
        .map(Result::unwrap)
        .skip_while(|op| {
            let kind = op.kind();
            kind == Kind::HardClip || kind == Kind::SoftClip
        })
        .try_for_each(|op| {
            let l = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMismatch => {
                    let mut acc = l;
                    while acc > 0 {
                        match md.pop() {
                            Some(MDKind::Match(m)) => {
                                if m as usize > acc {
                                    md.push(MDKind::Match(m - acc as u16));
                                    idx += acc;
                                    acc = 0;
                                } else {
                                    idx += m as usize;
                                    acc -= m as usize;
                                }
                            }
                            Some(MDKind::Substitution(_)) => {
                                if rec.quality_scores().as_ref()[idx + n_soft_clip] >= 30 {
                                    let alt_base = rec.sequence().as_ref()[idx + n_soft_clip];
                                    snv.push(SNV {
                                        relative_position: u32::try_from(idx)?,
                                        ty: Mutation::Substitution(alt_base),
                                    });
                                }
                                idx += 1;
                                acc -= 1;
                            }
                            _ => {
                                bail!(err_msg("Expecting MATCH", rec));
                            }
                        }
                    }
                    if acc != 0 {
                        bail!(err_msg("Length mismatch", rec));
                    }
                }
                Kind::Deletion => {
                    match md.pop() {
                        Some(MDKind::Deletion(deletion_bases)) => {
                            if deletion_bases.len() != l {
                                bail!(err_msg("Deletion length mismatch", rec));
                            }
                        }
                        _ => {
                            bail!(err_msg("Expecting DELETION", rec));
                        }
                    }
                    snv.push(SNV {
                        relative_position: u32::try_from(idx)?,
                        ty: Mutation::Deletion(l as u16),
                    });
                }
                Kind::Insertion => {
                    if rec.quality_scores().as_ref()[idx + n_soft_clip] >= 30 {
                        let alt_base = rec.sequence().as_ref()
                            [idx + n_soft_clip..idx + n_soft_clip + l]
                            .to_vec();
                        snv.push(SNV {
                            relative_position: u32::try_from(idx)?,
                            ty: Mutation::Insertion(alt_base),
                        });
                    }
                    idx += l;
                }
                Kind::Skip | Kind::Pad | Kind::HardClip => {}
                _ => {
                    // For other kinds of CIGAR operations, we do not record mutations.
                    // This includes soft clips, hard clips, etc.
                    idx += l;
                }
            };
            Ok(())
        })?;
    Ok(snv)
}

fn err_msg(msg: &str, rec: &RecordBuf) -> String {
    format!(
        "{}. CIGAR: {:?}, MD tag: {:?}",
        msg,
        rec.cigar(),
        rec.data().get(&Tag::MISMATCHED_POSITIONS).unwrap(),
    )
}

#[derive(Debug, PartialEq, Eq)]
enum MDKind {
    Match(u16),
    Substitution(u8),
    Deletion(Vec<u8>),
}

impl MDKind {
    fn from_str(s: &[u8]) -> Result<Option<(Self, &[u8])>> {
        if s.is_empty() {
            return Ok(None);
        }
        if s[0] == b'^' {
            // Deletion
            let mut l = 0;
            while 1 + l < s.len() && s[1 + l].is_ascii_alphabetic() {
                l += 1;
            }
            if l == 0 {
                bail!("Invalid MD tag: Deletion without bases");
            } else {
                let md = MDKind::Deletion(s[1..1 + l].to_vec());
                Ok(Some((md, &s[1 + l..])))
            }
        } else if s[0].is_ascii_digit() {
            let (m, n) = parse_partial::<u16, _>(s)?;
            let md = MDKind::Match(m);
            Ok(Some((md, &s[n..])))
        } else if s[0].is_ascii_alphabetic() {
            let md = MDKind::Substitution(s[0]);
            Ok(Some((md, &s[1..])))
        } else {
            bail!("Invalid MD tag: {}", String::from_utf8_lossy(s));
        }
    }
}

fn parse_md_tag_str(s: &[u8]) -> impl Iterator<Item = Result<MDKind>> + '_ {
    let mut s = s;
    std::iter::from_fn(move || {
        if s.is_empty() {
            return None;
        }
        match MDKind::from_str(s) {
            Ok(Some((md, rest))) => {
                s = rest;
                Some(Ok(md))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    })
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use bstr::BString;
    use noodles::sam::{
        self as sam,
        header::record::value::{map::ReferenceSequence, Map},
    };

    #[test]
    fn test_parse_md_tag_str() {
        let md = b"10A2^CTC3T";
        let mut iter = parse_md_tag_str(md);
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(10));
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Substitution(b'A'));
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(2));
        assert_eq!(
            iter.next().unwrap().unwrap(),
            MDKind::Deletion(vec![b'C', b'T', b'C'])
        );
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(3));
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Substitution(b'T'));
        assert!(iter.next().is_none());

        let md = b"61G0^GAGGGGGAGGGTGG52";
        let mut iter = parse_md_tag_str(md);
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(61));
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Substitution(b'G'));
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(0));
        assert_eq!(
            iter.next().unwrap().unwrap(),
            MDKind::Deletion(b"GAGGGGGAGGGTGG".to_vec()),
        );
        assert_eq!(iter.next().unwrap().unwrap(), MDKind::Match(52));
    }

    #[test]
    fn test_call_snv() -> Result<()> {
        let reference_sequences = [(
            BString::from("chr1"),
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(248956422)?),
        )]
        .into_iter()
        .collect();
        let header = sam::Header::builder()
            .set_reference_sequences(reference_sequences)
            .build();

        let sam = [
            "1\t99\tchr1\t10003\t0\t50M\t=\t10083\t129\tCCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC\tFFFF:FFFFF:FF:FF:FFFFFFFFFFFFFFFFF,FF:FF,FFFFFFFFF\tMD:Z:0A49",
            "2\t163\tchr1\t10004\t0\t49M\t=\t10064\t110\tCCCTAACCCAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC\tFFFFFFFFF,FFFFF:FFFFFFF:FFFFFFFFF:FFFFFFFFFFFFFFF\tMD:Z:9T39",
            "3\t163\tchr1\t10009\t0\t49M\t=\t10176\t213\tACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCAAACCCTAA\tFFFFFFFFFFFFFFFFFFFFFF,F,,FF,F:,,FF:F,F:F,F::F:::\tMD:Z:34T5T8",
            "4\t99\tchr1\t10010\t0\t50M\t=\t10296\t335\tCCCTAACCCTAACCCTAACCCTAACCCTAACCCTTACCCTTACCCTTACC\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tMD:Z:34A5A5A3",
            "5\t99\tchr1\t10010\t0\t5M2I4M\t=\t10296\t335\tCCCTGAACCCT\tFFFFFFFFFFF\tMD:Z:4A4",
            "6\t99\tchr1\t10010\t0\t8M2I4M\t=\t10296\t335\tCCCTGTTTAACCCT\tFFFFFFFFFFFFFF\tMD:Z:4A7",
            "7\t99\tchr1\t10010\t0\t5S8M2I4M\t=\t10296\t335\tAAAAACCCTGTTTAACCCT\tFFFFFFFFFFFFFFFFFFF\tMD:Z:4A7",
            "8\t99\tchr1\t10010\t0\t5H8M2I4M\t=\t10296\t335\tCCCTGTTTAACCCT\tFFFFFFFFFFFFFFFFFFF\tMD:Z:4A7",
            "9\t99\tchr1\t10010\t0\t23S61M14I1M14D52M\t=\t10296\t335\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tMD:Z:61G0^GAGGGGGAGGGTGG52",
        ];
        let ground_truth = vec![
            "0C", ".", "34A5A", "34T5T5T", "4G0+AA", "4G3+AA", "4G3+AA", "4G3+AA", "61+AAAAAAAAAAAAAA0A0^--------------",
        ];
        let rec: Vec<_> = sam
            .into_iter()
            .map(|x| {
                let r = noodles::sam::Record::try_from(x.as_bytes())?;
                let r = noodles::sam::alignment::RecordBuf::try_from_alignment_record(&header, &r)?;
                Ok(r)
            })
            .collect::<Result<_>>()?;

        rec.into_iter()
            .zip(ground_truth.into_iter())
            .enumerate()
            .for_each(|(i, (r, gt))| {
                let snv = get_snv(&r).unwrap();
                let snv = format_snvs(&snv).unwrap_or(".".to_string());
                assert_eq!(snv, gt, "SNV mismatch at record {}", i + 1);
            });

        Ok(())
    }
}
