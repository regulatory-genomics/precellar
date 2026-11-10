use noodles_sam::{
    self as sam,
    header::record::value::{map::ReferenceSequence, Map},
};
use precellar_core::ucf::sc_atac::{sort_rows as sort_sc_atac_rows, Dictionary, ScAtacRow};
use precellar_core::ucf::sc_rna::{row_from_record, sort_rows as sort_sc_rna_rows};
use std::num::NonZeroUsize;

fn header() -> sam::Header {
    let reference_sequences = [(
        bstr::BString::from("chr1"),
        Map::<ReferenceSequence>::new(NonZeroUsize::try_from(248_956_422).unwrap()),
    )]
    .into_iter()
    .collect();
    sam::Header::builder()
        .set_reference_sequences(reference_sequences)
        .build()
}

fn record(sam_line: &str) -> noodles_sam::alignment::RecordBuf {
    let header = header();
    let record = noodles_sam::Record::try_from(sam_line.as_bytes()).unwrap();
    noodles_sam::alignment::RecordBuf::try_from_alignment_record(&header, &record).unwrap()
}

const PAIRED: &str = "1\t99\tchr1\t10003\t60\t4H50M2H\t=\t10083\t129\tCCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC\tFFFF:FFFFF:FF:FF:FFFFFFFFFFFFFFFFF,FF:FF,FFFFFFFF\tCB:Z:ACGTACGT\tRX:Z:TTTTGGGG\tNH:i:1\tMD:Z:50";

#[test]
fn sc_rna_row_keeps_cell_umi_and_cigar_semantics() {
    let mut cells = Dictionary::default();
    let mut umis = Dictionary::default();
    let row = row_from_record(
        &record(PAIRED),
        |barcode| cells.id(barcode),
        |umi| umis.id(umi),
    )
    .unwrap()
    .expect("mapped record with barcode");

    assert_eq!(row.chrom_id, 0);
    assert_eq!(row.start, 10_002);
    assert_eq!(row.mate_chrom_id, Some(0));
    assert_eq!(row.mate_start, Some(10_082));
    assert_eq!(row.template_length, 129);
    assert_eq!(row.flags, 99);
    assert_eq!(row.nh, 1);
    assert_eq!(row.mapq, 60);
    assert_eq!(row.reference_span, 50);
    assert_eq!(row.block_count, 1);
    assert_eq!(row.cigar, vec![(4, 5), (50, 0), (2, 5)]);
    assert_eq!(cells.entries(), &["ACGTACGT".to_owned()]);
    assert_eq!(umis.entries(), &["TTTTGGGG".to_owned()]);
    assert_eq!(row.cell_id, 0);
    assert_eq!(row.umi_id, 0);
}

#[test]
fn sc_rna_keeps_spliced_blocks_and_sentinel_without_umi() {
    let line = "1\t0\tchr1\t100\t60\t10M1000N10M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tFFFFFFFFFFFFFFFFFFFF\tCB:Z:AAACCG";
    let row = row_from_record(&record(line), |_| 0, |_| 0)
        .unwrap()
        .expect("mapped record");
    assert_eq!(row.reference_span, 1_020);
    assert_eq!(row.block_count, 2);
    assert_eq!(row.umi_id, u32::MAX);
    assert_eq!(row.cigar, vec![(10, 0), (1000, 3), (10, 0)]);
}

#[test]
fn unmapped_records_are_skipped() {
    let line = "1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF";
    assert!(row_from_record(&record(line), |_| 0, |_| 0)
        .unwrap()
        .is_none());
}

#[test]
fn sc_atac_rows_sort_in_ucf_order() {
    let rows = vec![
        ScAtacRow {
            chrom_id: 1,
            start: 10,
            end: 20,
            cell_id: 0,
            multiplicity: 1,
            flags: 1,
        },
        ScAtacRow {
            chrom_id: 0,
            start: 30,
            end: 40,
            cell_id: 1,
            multiplicity: 1,
            flags: 2,
        },
        ScAtacRow {
            chrom_id: 0,
            start: 10,
            end: 20,
            cell_id: 2,
            multiplicity: 3,
            flags: 1,
        },
        ScAtacRow {
            chrom_id: 0,
            start: 10,
            end: 25,
            cell_id: 0,
            multiplicity: 1,
            flags: 0,
        },
    ];
    let sorted = sort_sc_atac_rows(rows.into_iter(), None, 2).collect::<Vec<_>>();
    let keys = sorted
        .iter()
        .map(|row| (row.chrom_id, row.start, row.cell_id))
        .collect::<Vec<_>>();
    assert_eq!(keys, vec![(0, 10, 0), (0, 10, 2), (0, 30, 1), (1, 10, 0)]);
}

#[test]
fn sc_rna_rows_sort_by_locus_then_cell() {
    let mut rows = Vec::new();
    for (start, cell) in [(50u32, 1u32), (10, 2), (10, 0), (90, 0)] {
        rows.push(precellar_core::ucf::sc_rna::ScRnaRow {
            chrom_id: 0,
            start,
            mate_chrom_id: None,
            mate_start: None,
            template_length: 0,
            flags: 0,
            nh: 1,
            mapq: 60,
            cell_id: cell,
            umi_id: 0,
            reference_span: 10,
            block_count: 1,
            cigar: vec![(10, 0)],
        });
    }
    let sorted = sort_sc_rna_rows(rows.into_iter(), None, 2).collect::<Vec<_>>();
    let keys = sorted.iter().map(|row| (row.start, row.cell_id)).collect::<Vec<_>>();
    assert_eq!(keys, vec![(10, 0), (10, 2), (50, 1), (90, 0)]);
}
