mod fastq;
mod bam;

use pyo3::prelude::*;

#[pymodule]
pub(crate) fn register_utils(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let utils = PyModule::new(parent_module.py(), "utils")?;

    utils.add_function(wrap_pyfunction!(fastq::strip_barcode_from_fastq, &utils)?)?;
    utils.add_function(wrap_pyfunction!(fastq::multiplex_fastq, &utils)?)?;

    utils.add_function(wrap_pyfunction!(bam::bam_to_fastq, &utils)?)?;

    parent_module.add_submodule(&utils)
}