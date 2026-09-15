//! GTF annotation parser.
//!
//! Produces the same internal [`Transcript`] representation as the STAR index
//! loader, so that both annotation sources feed a single `TxAligner`:
//!
//! ```text
//! STAR index annotation ─┐
//!                        ├─> precellar Transcript ─> TxAligner
//! GTF parser            ─┘
//! ```
//!
//! Only `exon` lines are read. Transcript bounds are derived from the exon
//! envelope, which matches what STAR itself does when building its index (it
//! asserts the transcript end equals the last exon end). This also means GTFs
//! that omit `transcript` lines parse without special handling.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{anyhow, bail, ensure, Context, Result};
use bed_utils::bed::Strand;
use log::{info, warn};
use noodles_gff::feature::record::Strand as GtfStrand;
use noodles_gtf as gtf;
use noodles_sam::Header;

use crate::transcriptome::{Exons, Transcript};

const EXON: &[u8] = b"exon";
const TRANSCRIPT_ID: &[u8] = b"transcript_id";
const GENE_ID: &[u8] = b"gene_id";
const GENE_NAME: &[u8] = b"gene_name";

/// Exons accumulated for a single `transcript_id` while scanning the file.
struct PartialTranscript {
    gene_id: String,
    gene_name: String,
    chrom: String,
    strand: Strand,
    /// 0-based, half-open, in the order encountered in the file.
    exons: Vec<(u64, u64)>,
    /// Line number of the first exon, used to make errors locatable.
    first_line: usize,
}

/// Reads transcript annotation from a GTF file.
///
/// The file may be plain, gzip-, or zstd-compressed. Exon records are grouped
/// by `transcript_id`.
pub fn read_transcripts<P: AsRef<Path>>(path: P) -> Result<Vec<Transcript>> {
    let reader = seqspec::utils::open_file(path.as_ref())
        .with_context(|| format!("failed to open GTF '{}'", path.as_ref().display()))?;
    parse_transcripts(std::io::BufReader::new(reader)).with_context(|| {
        format!(
            "failed to parse transcript annotation from GTF '{}'",
            path.as_ref().display()
        )
    })
}

/// Parses transcript annotation from any GTF byte stream.
pub fn parse_transcripts<R: BufRead>(reader: R) -> Result<Vec<Transcript>> {
    let mut reader = gtf::io::Reader::new(reader);
    let mut line = gtf::Line::default();
    let mut partials: HashMap<String, PartialTranscript> = HashMap::new();
    let mut line_no = 0;
    let mut num_exons = 0;

    loop {
        line_no += 1;
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let Some(record) = line.as_record() else {
            continue; // comment or track line
        };
        let record = record.with_context(|| format!("malformed GTF record at line {line_no}"))?;
        if record.ty() != EXON {
            continue;
        }
        num_exons += 1;

        let attributes = record
            .attributes()
            .with_context(|| format!("malformed attributes at line {line_no}"))?;
        let attribute = |key: &[u8]| -> Result<Option<String>> {
            attributes
                .get(key)
                .transpose()
                .with_context(|| format!("malformed attributes at line {line_no}"))?
                .map(|value| {
                    let mut values = value.iter();
                    let first = values.next().ok_or_else(|| {
                        anyhow!(
                            "empty '{}' at line {}",
                            String::from_utf8_lossy(key),
                            line_no
                        )
                    })?;
                    ensure!(
                        values.next().is_none(),
                        "multiple '{}' values at line {}",
                        String::from_utf8_lossy(key),
                        line_no
                    );
                    Ok(String::from_utf8_lossy(first).into_owned())
                })
                .transpose()
        };

        let transcript_id = attribute(TRANSCRIPT_ID)?
            .ok_or_else(|| anyhow!("missing 'transcript_id' at line {line_no}"))?;
        let gene_id =
            attribute(GENE_ID)?.ok_or_else(|| anyhow!("missing 'gene_id' at line {line_no}"))?;
        // STAR's geneInfo.tab falls back to the gene ID when no symbol is
        // available, so mirror that here.
        let gene_name = attribute(GENE_NAME)?.unwrap_or_else(|| gene_id.clone());

        let strand = match record
            .strand()
            .with_context(|| format!("malformed strand at line {line_no}"))?
        {
            GtfStrand::Forward => Strand::Forward,
            GtfStrand::Reverse => Strand::Reverse,
            // STAR's loader also rejects transcripts without a definite strand.
            other => bail!(
                "exon at line {} has strand {:?}, but transcript annotation requires '+' or '-'",
                line_no,
                other
            ),
        };

        // GTF is 1-based and closed; the internal representation is 0-based and
        // half-open.
        let start = record.start()?.get() as u64 - 1;
        let end = record.end()?.get() as u64;
        ensure!(
            start < end,
            "exon at line {line_no} has start greater than end"
        );
        let chrom = String::from_utf8_lossy(record.reference_sequence_name()).into_owned();

        match partials.get_mut(&transcript_id) {
            Some(partial) => {
                ensure!(
                    partial.chrom == chrom,
                    "transcript '{}' spans multiple chromosomes ('{}' at line {}, '{}' at line {})",
                    transcript_id,
                    partial.chrom,
                    partial.first_line,
                    chrom,
                    line_no
                );
                ensure!(
                    partial.strand == strand,
                    "transcript '{}' spans both strands (line {} and line {})",
                    transcript_id,
                    partial.first_line,
                    line_no
                );
                ensure!(
                    partial.gene_id == gene_id,
                    "transcript '{}' is assigned to multiple genes ('{}' at line {}, '{}' at line {})",
                    transcript_id,
                    partial.gene_id,
                    partial.first_line,
                    gene_id,
                    line_no
                );
                partial.exons.push((start, end));
            }
            None => {
                partials.insert(
                    transcript_id,
                    PartialTranscript {
                        gene_id,
                        gene_name,
                        chrom,
                        strand,
                        exons: vec![(start, end)],
                        first_line: line_no,
                    },
                );
            }
        }
    }

    ensure!(
        !partials.is_empty(),
        "found no 'exon' records in GTF ({num_exons} exon lines seen)"
    );

    // No sort here: `TxAligner::new` collects into a `GIntervalMap`, whose
    // `Lapper` sorts by coordinate on construction, so any order imposed here
    // would be discarded.
    let transcripts = partials
        .into_iter()
        .map(|(id, partial)| partial.into_transcript(id))
        .collect::<Result<Vec<_>>>()?;

    info!(
        "Loaded {} transcripts ({} exons) from GTF",
        transcripts.len(),
        num_exons
    );
    Ok(transcripts)
}

impl PartialTranscript {
    fn into_transcript(mut self, id: String) -> Result<Transcript> {
        // Exons on the reverse strand are listed in descending order in GTF, so
        // always re-sort. Fully duplicated exons are dropped; genuine overlaps
        // are an error, since `Exons` requires non-overlapping intervals.
        self.exons.sort_unstable();
        self.exons.dedup();
        for pair in self.exons.windows(2) {
            ensure!(
                pair[1].0 >= pair[0].1,
                "transcript '{}' (line {}) has overlapping exons: [{}, {}) and [{}, {})",
                id,
                self.first_line,
                pair[0].0,
                pair[0].1,
                pair[1].0,
                pair[1].1
            );
        }

        // The transcript interval is the exon envelope.
        let start = self.exons.first().unwrap().0;
        let end = self.exons.last().unwrap().1;
        let exons = Exons::new(self.exons)
            .with_context(|| format!("invalid exons for transcript '{id}' (line {})", self.first_line))?;

        Ok(Transcript {
            id,
            chrom: self.chrom,
            start,
            end,
            strand: self.strand,
            gene_id: self.gene_id,
            gene_name: self.gene_name,
            exons,
        })
    }
}

/// Checks GTF transcripts against the alignment header.
///
/// Unlike the STAR index check, this is one-directional: a GTF only annotates
/// the assembled chromosomes, so it is normally a strict subset of the header
/// (scaffolds carry no annotation). Only two things are hard errors: nothing
/// matching at all, which usually means a `chr1` versus `1` naming mismatch,
/// and a transcript running past the end of its reference, which indicates a
/// coordinate-system mismatch.
pub fn validate_against_header(transcripts: &[Transcript], header: &Header) -> Result<()> {
    let references: HashMap<&str, u64> = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            (
                std::str::from_utf8(name.as_ref()).unwrap_or_default(),
                reference.length().get() as u64,
            )
        })
        .collect();

    let mut matched = 0;
    let mut unmatched: Vec<&str> = Vec::new();
    let mut gtf_chroms: Vec<&str> = transcripts.iter().map(|t| t.chrom.as_str()).collect();
    gtf_chroms.sort_unstable();
    gtf_chroms.dedup();

    for chrom in &gtf_chroms {
        if references.contains_key(chrom) {
            matched += 1;
        } else {
            unmatched.push(chrom);
        }
    }

    if matched == 0 {
        bail!(
            "no GTF chromosome matches the alignment header: GTF has {} ({}), header has {} ({}). \
             Check that both use the same naming convention, for example 'chr1' versus '1'.",
            gtf_chroms.len(),
            summarize(&gtf_chroms),
            references.len(),
            summarize(&references.keys().copied().collect::<Vec<_>>()),
        );
    }

    for transcript in transcripts {
        if let Some(length) = references.get(transcript.chrom.as_str()) {
            ensure!(
                transcript.end <= *length,
                "transcript '{}' ends at {} on '{}', past the reference length {}. \
                 The GTF and the aligner index are built from different assemblies.",
                transcript.id,
                transcript.end,
                transcript.chrom,
                length
            );
        }
    }

    if !unmatched.is_empty() {
        warn!(
            "{} of {} GTF chromosomes are absent from the alignment header and will be ignored: {}",
            unmatched.len(),
            gtf_chroms.len(),
            summarize(&unmatched),
        );
    }
    info!(
        "Matched {}/{} GTF chromosomes against the alignment header",
        matched,
        gtf_chroms.len()
    );
    Ok(())
}

fn summarize(names: &[&str]) -> String {
    if names.is_empty() {
        return "none".to_owned();
    }
    let shown = names.iter().take(5).copied().collect::<Vec<_>>().join(", ");
    if names.len() > 5 {
        format!("{shown}, ... ")
    } else {
        shown
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    use std::num::NonZeroUsize;

    fn exon(chrom: &str, start: u64, end: u64, strand: &str, attrs: &str) -> String {
        format!("{chrom}\tsrc\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
    }

    /// The conversion that must match STAR exactly: 1-based closed GTF
    /// coordinates become 0-based half-open, transcript bounds come from the
    /// exon envelope rather than any `transcript` line, and reverse-strand
    /// exons (listed descending in GTF) are stored ascending.
    #[test]
    fn builds_transcripts_from_exon_lines() {
        let fwd = r#"gene_id "g1"; transcript_id "t1"; gene_name "G1";"#;
        let rev = r#"gene_id "g2"; transcript_id "t2";"#;
        let src = format!("#!comment\nchr1\tsrc\tgene\t1\t9999\t.\t+\t.\t{fwd}\n")
            + &exon("chr1", 101, 150, "+", fwd)
            + &exon("chr1", 301, 400, "+", fwd)
            + &exon("chr2", 601, 700, "-", rev)
            + &exon("chr2", 501, 550, "-", rev);

        let transcripts = parse_transcripts(src.as_bytes()).unwrap();
        assert_eq!(transcripts.len(), 2);
        let by_id: std::collections::HashMap<_, _> =
            transcripts.iter().map(|t| (t.id.as_str(), t)).collect();

        let t1 = by_id["t1"];
        assert_eq!((t1.start, t1.end), (100, 400)); // envelope, not the gene line
        assert_eq!(
            t1.exons()
                .iter()
                .map(|e| (e.start(), e.end()))
                .collect::<Vec<_>>(),
            vec![(100, 150), (300, 400)]
        );
        assert_eq!(t1.gene_name, "G1");

        let t2 = by_id["t2"];
        assert_eq!(t2.strand, Strand::Reverse);
        // Descending in the file, ascending in the structure.
        assert_eq!(
            t2.exons()
                .iter()
                .map(|e| (e.start(), e.end()))
                .collect::<Vec<_>>(),
            vec![(500, 550), (600, 700)]
        );
        // No gene_name attribute, so it falls back to gene_id, as STAR does for
        // a single-column geneInfo.tab.
        assert_eq!(t2.gene_name, "g2");
    }

    /// Malformed annotation must fail loudly rather than silently produce a
    /// transcript set that quietly miscounts. Each of these would otherwise
    /// corrupt downstream counts.
    #[test]
    fn rejects_malformed_annotation() {
        let case = |src: String, expected: &str| {
            let err = parse_transcripts(src.as_bytes()).unwrap_err().to_string();
            assert!(err.contains(expected), "expected {expected:?}, got {err:?}");
        };
        let attrs = r#"gene_id "g1"; transcript_id "t1";"#;

        case(
            exon("chr1", 1, 10, "+", r#"gene_id "g1";"#),
            "missing 'transcript_id' at line 1",
        );
        case(
            exon("chr1", 1, 10, ".", attrs),
            "requires '+' or '-'",
        );
        // Overlapping exons would break the binary search in `find_exons`.
        case(
            exon("chr1", 1, 20, "+", attrs) + &exon("chr1", 10, 30, "+", attrs),
            "overlapping exons",
        );
        case(
            exon("chr1", 1, 10, "+", attrs) + &exon("chr2", 1, 10, "+", attrs),
            "spans multiple chromosomes",
        );
        case(
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n".to_owned(),
            "no 'exon' records",
        );
    }

    /// A GTF annotates assembled chromosomes only, so being a subset of the
    /// header is normal; a disjoint naming scheme (`chr1` vs `1`) and
    /// coordinates past the reference end are not.
    #[test]
    fn validates_chromosomes_against_header() {
        fn header(entries: &[(&str, usize)]) -> Header {
            let mut builder = Header::builder();
            for (name, length) in entries {
                builder = builder.add_reference_sequence(
                    *name,
                    Map::<ReferenceSequence>::new(NonZeroUsize::new(*length).unwrap()),
                );
            }
            builder.build()
        }

        let attrs = r#"gene_id "g1"; transcript_id "t1";"#;
        let transcripts = parse_transcripts(exon("chr1", 1, 10, "+", attrs).as_bytes()).unwrap();

        // Unannotated scaffolds in the header are fine.
        validate_against_header(&transcripts, &header(&[("chr1", 1000), ("GL456210.1", 500)]))
            .unwrap();

        let err = validate_against_header(&transcripts, &header(&[("1", 1000)]))
            .unwrap_err()
            .to_string();
        assert!(err.contains("no GTF chromosome matches"), "{err}");
        assert!(err.contains("'chr1' versus '1'"), "{err}");

        let long = parse_transcripts(exon("chr1", 1, 2000, "+", attrs).as_bytes()).unwrap();
        let err = validate_against_header(&long, &header(&[("chr1", 1000)]))
            .unwrap_err()
            .to_string();
        assert!(err.contains("past the reference length"), "{err}");
    }
}
