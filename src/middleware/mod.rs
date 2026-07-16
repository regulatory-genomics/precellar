pub(crate) mod tfseq;

pub(crate) use tfseq::FloatingBarcodeFinder;

use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(crate) fn register_middleware(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let middleware = PyModule::new(parent_module.py(), "middleware")?;
    let tfseq = PyModule::new(parent_module.py(), "tfseq")?;
    tfseq::register(&tfseq)?;
    middleware.add_submodule(&tfseq)?;
    parent_module.add_submodule(&middleware)?;

    let sys = PyModule::import(parent_module.py(), "sys")?;
    let modules = sys.getattr("modules")?.cast_into::<PyDict>()?;
    modules.set_item("precellar.middleware", &middleware)?;
    modules.set_item("precellar.middleware.tfseq", &tfseq)
}
