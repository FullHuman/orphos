mod common;
use crate::common::{
    normalize_output, prodigal_available, run_original_prodigal, run_orphos, similarity,
};
use tempfile::NamedTempFile;

const SAMPLE_FASTA: &str =
    ">seq\nATGAAAAAACTATTAACCTCTCTGCTGCTGTTGGCCGCAGCCGCCACCCTGACCACCGAGGCCATCAAGAACCTGGGCTGA";

#[test]
#[ignore]
fn gbk_format_equivalence_meta_mode() {
    if !prodigal_available() {
        eprintln!("Skipping: original prodigal not in PATH");
        return;
    }
    let input = NamedTempFile::new().unwrap();
    std::fs::write(input.path(), SAMPLE_FASTA).unwrap();
    let orig_out = NamedTempFile::new().unwrap();
    let rust_out = NamedTempFile::new().unwrap();

    run_original_prodigal(
        input.path().to_str().unwrap(),
        orig_out.path().to_str().unwrap(),
        "gbk",
        "meta",
    )
    .unwrap();
    run_orphos(
        input.path().to_str().unwrap(),
        rust_out.path().to_str().unwrap(),
        "gbk",
        "meta",
    )
    .unwrap();

    let a = normalize_output(&std::fs::read_to_string(orig_out.path()).unwrap());
    let b = normalize_output(&std::fs::read_to_string(rust_out.path()).unwrap());
    let sim = similarity(&a, &b);
    println!("GBK meta similarity: {:.2}%", sim * 100.0);
    assert!(sim > 0.95, "Similarity too low: {:.2}%", sim * 100.0);
}
