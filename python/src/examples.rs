use anyhow::Result;
use pyo3::prelude::*;
use std::{collections::HashMap, path::PathBuf};

#[pyfunction]
fn txg_multiome(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/869db",
        "md5:8fd4a6e69dfd4f442b19e8f34d62e481",
        "10x_multiome_rna_R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/n9e87",
        "md5:8d9eda7dfe65753aa9ea4219b0bfeb3e",
        "10x_multiome_rna_R2.fq.zst",
    )?;
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/ebk5v",
        "md5:f9b17e828835ab4e18f68cadbb72d6d3",
        "10x_multiome_atac_R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/qaxu6/",
        "md5:af1817ebc06ff0a33f2e07ff5fcfb901",
        "10x_multiome_atac_R2.fq.zst",
    )?;
    let atac_i2 = retrieve_file(
        py,
        "https://osf.io/download/xgyt6/",
        "md5:ecb5dc088bcd88c5b6e319e81dc4ad5e",
        "10x_multiome_atac_I2.fq.zst",
    )?;
    Ok(HashMap::from([
        ("rna-R1", rna_r1),
        ("rna-R2", rna_r2),
        ("atac-R1", atac_r1),
        ("atac-R2", atac_r2),
        ("atac-I2", atac_i2),
    ]))
}

#[pyfunction]
fn sci_rna_seq3(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/wmhj5",
        "md5:4fe3a8f5410f498138ae48498431577a",
        "sci_rna_seq3_R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/nc3py",
        "md5:7905d8bbcaf5834cebabfb262e896f83",
        "sci_rna_seq3_R2.fq.zst",
    )?;
    Ok(HashMap::from([("R1", rna_r1), ("R2", rna_r2)]))
}

#[pyfunction]
fn dsc_atac(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/wtn72",
        "md5:d9fb3a624be9d39a55fbbc0de6cf7607",
        "dsc_atac_R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/khcxq/",
        "md5:8b6c0d58eac183657ac9e1d3b69682d9",
        "dsc_atac_R2.fq.zst",
    )?;
    Ok(HashMap::from([("R1", atac_r1), ("R2", atac_r2)]))
}

fn retrieve_file(py: Python<'_>, url: &str, known_hash: &str, fname: &str) -> Result<PathBuf> {
    let pooch = PyModule::import(py, "pooch")?;
    let kwargs = pyo3::types::PyDict::new(py);
    kwargs.set_item("progressbar", true)?;
    Ok(pooch
        .getattr("retrieve")?
        .call((url, known_hash, fname), Some(&kwargs))?
        .extract()?)
}

#[pymodule]
pub(crate) fn register_examples(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "examples")?;

    m.add_function(wrap_pyfunction!(txg_multiome, &m)?)?;
    m.add_function(wrap_pyfunction!(sci_rna_seq3, &m)?)?;
    m.add_function(wrap_pyfunction!(dsc_atac, &m)?)?;

    parent_module.add_submodule(&m)
}
