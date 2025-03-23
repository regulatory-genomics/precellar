use anyhow::Result;
use pyo3::prelude::*;
use std::{collections::HashMap, path::PathBuf};

#[pyfunction]
fn txg_rna_v3(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/xpbr3",
        "md5:64e23447a85383944dcf241bd15b5bb6",
        "10x_rna_v3_R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/gya7p",
        "md5:8abfb22a2f3621aca0e56d494144895c",
        "10x_rna_v3_R2.fq.zst",
    )?;
    Ok(HashMap::from([("R1", rna_r1), ("R2", rna_r2)]))
}

#[pyfunction]
fn mars_seq(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/s45kg",
        "md5:3daa417056e31aa8fe2f3de504c50ec4",
        "mars_seq_rna_R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/r3vdu",
        "md5:933a6136018f59d44bacc31570d56c26",
        "mars_seq_rna_R2.fq.zst",
    )?;
    Ok(HashMap::from([("R1", rna_r1), ("R2", rna_r2)]))
}

#[pyfunction]
fn txg_atac(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let atac_i2 = retrieve_file(
        py,
        "https://osf.io/download/nv34x",
        "md5:6b2a6e536b8c351ae08f574f74c43e83",
        "10x_atac_I2.fq.zst",
    )?;
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/7dm48",
        "md5:0ee922b6f5a4774eb3cc425e93f7c532",
        "10x_atac_R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/vf5ya",
        "md5:c9257abf04f5f9af7d0e2bd21de761c0",
        "10x_atac_R2.fq.zst",
    )?;
    Ok(HashMap::from([("I2", atac_i2), ("R1", atac_r1), ("R2", atac_r2)]))
}

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

#[pyfunction]
fn share_seq(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_i1 = retrieve_file(
        py,
        "https://osf.io/download/p974h",
        "md5:1f354afb600fc1dda4738572d8c1b633",
        "share_seq_rna-I1.fq.zst",
    )?;
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/ysvdu",
        "md5:1dd61b3dfdc6ca8553f234fb2cbce742",
        "share_seq_rna-R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/8h4mz",
        "md5:5712bd18c13b4bb70ea38e1b3222d7a9",
        "share_seq_rna-R2.fq.zst",
    )?;
    let atac_i1 = retrieve_file(
        py,
        "https://osf.io/download/wcvpe",
        "md5:88d1a7220538963fffc30b56c6914753",
        "share_seq_atac-I1.fq.zst",
    )?;
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/egvt3",
        "md5:ab85a5648ba58ecccd2c67fe181cdc73",
        "share_seq_atac-R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/cnwpy/",
        "md5:639fcd024dfc4775c2bbef435eb4483b",
        "share_seq_atac-R2.fq.zst",
    )?;
    Ok(HashMap::from([
        ("rna-I1", rna_i1),
        ("rna-R1", rna_r1),
        ("rna-R2", rna_r2),
        ("atac-I1", atac_i1),
        ("atac-R1", atac_r1),
        ("atac-R2", atac_r2),
    ]))
}

#[pyfunction]
fn snare_seq(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/f4jat",
        "md5:816d88d618fa51be87e149074932346c",
        "snare_seq_rna-R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/pa4dm",
        "md5:6261446f297b44cc196816a09911750e",
        "snare_seq_rna-R2.fq.zst",
    )?;
    let atac_i1 = retrieve_file(
        py,
        "https://osf.io/download/b2stp",
        "md5:1264525b00e6598ba72389b90b114b57",
        "snare_seq_atac-I1.fq.zst",
    )?;
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/r7hjx",
        "md5:c35ef90ddc110d32d7bcb6f14db3fad4",
        "snare_seq_atac-R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/rqbkw/",
        "md5:73c6b2ccef40a0109c762104f12068dd",
        "snare_seq_atac-R2.fq.zst",
    )?;
    Ok(HashMap::from([
        ("rna-R1", rna_r1),
        ("rna-R2", rna_r2),
        ("atac-I1", atac_i1),
        ("atac-R1", atac_r1),
        ("atac-R2", atac_r2),
    ]))
}

#[pyfunction]
fn droplet_paired_tag(py: Python<'_>) -> Result<HashMap<&str, PathBuf>> {
    let rna_r1 = retrieve_file(
        py,
        "https://osf.io/download/zat7v",
        "md5:d69f4d0e04f8257868f63d904cf3ca83",
        "droplet_paired_tag_rna_R1.fq.zst",
    )?;
    let rna_r2 = retrieve_file(
        py,
        "https://osf.io/download/3r9tg",
        "md5:08ff0fe213600c9ce73e3d507125a964",
        "droplet_paired_tag_rna_R2.fq.zst",
    )?;
    let atac_r1 = retrieve_file(
        py,
        "https://osf.io/download/mxquz",
        "md5:8d02974f25f0ab28bfc92f6cc8e8bcbf",
        "droplet_paired_tag_atac_R1.fq.zst",
    )?;
    let atac_r2 = retrieve_file(
        py,
        "https://osf.io/download/e6pyq/",
        "md5:c0859bceb4ddc9b12f6c25874edd6903",
        "droplet_paired_tag_atac_R2.fq.zst",
    )?;
    let atac_i2 = retrieve_file(
        py,
        "https://osf.io/download/m63pw/",
        "md5:44efed45d305c60e0b92631b0e44cc73",
        "droplet_paired_tag_atac_I2.fq.zst",
    )?;
    Ok(HashMap::from([
        ("rna-R1", rna_r1),
        ("rna-R2", rna_r2),
        ("atac-R1", atac_r1),
        ("atac-R2", atac_r2),
        ("atac-I2", atac_i2),
    ]))
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

    m.add_function(wrap_pyfunction!(txg_rna_v3, &m)?)?;
    m.add_function(wrap_pyfunction!(txg_atac, &m)?)?;
    m.add_function(wrap_pyfunction!(sci_rna_seq3, &m)?)?;
    m.add_function(wrap_pyfunction!(dsc_atac, &m)?)?;
    m.add_function(wrap_pyfunction!(txg_multiome, &m)?)?;
    m.add_function(wrap_pyfunction!(share_seq, &m)?)?;
    m.add_function(wrap_pyfunction!(snare_seq, &m)?)?;
    m.add_function(wrap_pyfunction!(droplet_paired_tag, &m)?)?;
    m.add_function(wrap_pyfunction!(mars_seq, &m)?)?;

    parent_module.add_submodule(&m)
}
