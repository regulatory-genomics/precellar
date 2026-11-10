use noodles_sam::{
    self as sam,
    header::record::value::{map::ReferenceSequence, Map},
};
use precellar_core::ucf::sc_atac::{sort_rows as sort_sc_atac_rows, Dictionary, ScAtacRow};
use precellar_core::ucf::sc_rna::{
    classify_record, sort_rows as sort_sc_rna_rows, AlignmentPolicy, MissingBarcodePolicy,
    ScRnaDisposition, ScRnaRow,
};
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

fn stored(disposition: ScRnaDisposition) -> ScRnaRow {
    match disposition {
        ScRnaDisposition::Stored(row) => row,
        other => panic!("expected a stored row, got {other:?}"),
    }
}

fn classify(sam_line: &str) -> ScRnaDisposition {
    classify_record(
        &record(sam_line),
        |_| 0,
        |_| 0,
        AlignmentPolicy::Primary,
        MissingBarcodePolicy::Drop,
    )
    .unwrap()
}

const PAIRED: &str = "1\t99\tchr1\t10003\t60\t4H50M2H\t=\t10083\t129\tCCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC\tFFFF:FFFFF:FF:FF:FFFFFFFFFFFFFFFFF,FF:FF,FFFFFFFF\tCB:Z:ACGTACGT\tRX:Z:TTTTGGGG\tNH:i:1\tMD:Z:50";

#[test]
fn sc_rna_row_keeps_cell_umi_and_cigar_semantics() {
    let mut cells = Dictionary::default();
    let mut umis = Dictionary::default();
    let row = stored(
        classify_record(
            &record(PAIRED),
            |barcode| cells.id(barcode),
            |umi| umis.id(umi),
            AlignmentPolicy::Primary,
            MissingBarcodePolicy::Drop,
        )
        .unwrap(),
    );

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
    let row = stored(classify(line));
    assert_eq!(row.reference_span, 1_020);
    assert_eq!(row.block_count, 2);
    assert_eq!(row.umi_id, u32::MAX);
    assert_eq!(row.cigar, vec![(10, 0), (1000, 3), (10, 0)]);
}

#[test]
fn unmapped_records_are_classified_not_stored() {
    let line = "1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF";
    assert!(matches!(classify(line), ScRnaDisposition::Unmapped));
}

#[test]
fn primary_policy_excludes_secondary_and_supplementary() {
    let secondary = "1\t256\tchr1\t100\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AAACCG";
    let supplementary = "1\t2048\tchr1\t100\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tCB:Z:AAACCG";
    assert!(matches!(
        classify(secondary),
        ScRnaDisposition::Secondary
    ));
    assert!(matches!(
        classify(supplementary),
        ScRnaDisposition::Supplementary
    ));

    for line in [secondary, supplementary] {
        let disposition = classify_record(
            &record(line),
            |_| 0,
            |_| 0,
            AlignmentPolicy::All,
            MissingBarcodePolicy::Drop,
        )
        .unwrap();
        assert!(matches!(disposition, ScRnaDisposition::Stored(_)));
    }
}

#[test]
fn missing_barcode_is_dropped_or_reported() {
    let line = "1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    assert!(matches!(
        classify(line),
        ScRnaDisposition::MissingBarcode
    ));
    let error = classify_record(
        &record(line),
        |_| 0,
        |_| 0,
        AlignmentPolicy::Primary,
        MissingBarcodePolicy::Error,
    )
    .unwrap_err();
    assert!(error.to_string().contains("cell barcode"), "{error}");
}

#[test]
fn retention_counts_every_class() {
    use precellar_core::ucf::sc_rna::ScRnaRetention;

    let mut retention = ScRnaRetention::default();
    retention.record(&classify(PAIRED));
    retention.record(&classify("1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tFFFF"));
    assert_eq!(retention.stored, 1);
    assert_eq!(retention.unmapped, 1);
    assert_eq!(retention.excluded(), 1);
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
    let sorted = sort_sc_atac_rows(rows.into_iter(), None, 2)
        .unwrap()
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
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
        rows.push(ScRnaRow {
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
    let sorted = sort_sc_rna_rows(rows.into_iter(), None, 2)
        .unwrap()
        .collect::<anyhow::Result<Vec<_>>>()
        .unwrap();
    let keys = sorted
        .iter()
        .map(|row| (row.start, row.cell_id))
        .collect::<Vec<_>>();
    assert_eq!(keys, vec![(10, 0), (10, 2), (50, 1), (90, 0)]);
}

#[test]
fn prepare_temp_dir_creates_the_directory_and_rejects_files() {
    use precellar_core::ucf::prepare_temp_dir;

    let directory = tempfile::tempdir().unwrap();
    let missing = directory.path().join("nested/tmp");
    prepare_temp_dir(Some(&missing)).unwrap();
    assert!(missing.is_dir());

    let file = directory.path().join("not-a-dir");
    std::fs::write(&file, b"x").unwrap();
    assert!(prepare_temp_dir(Some(&file)).is_err());
}
