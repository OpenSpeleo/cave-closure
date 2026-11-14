//! Independent linear survey loop closure in a caller-supplied local ENU frame.
//! Gaussian covariance models and numerical backends are separate choices.
mod graph;
pub mod json;
pub mod loop_closure;
mod model;
pub mod normalize;
mod prepare;
#[cfg(feature = "python")]
mod python;
mod solve;
pub use model::*;
pub use prepare::{PreparedProblem, prepare};
#[cfg(feature = "generate-stubs")]
pub use python::stub_info;
pub use solve::{SolveControl, Workspace};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("validation: {0}")]
    Validation(String),
    #[error("numerical failure: {0}")]
    Numerical(String),
    #[error("incomplete solve: {0:?}")]
    Incomplete(Termination),
    #[error("JSON: {0}")]
    Json(#[from] serde_json::Error),
}
pub(crate) fn invalid(message: impl Into<String>) -> Error {
    Error::Validation(message.into())
}
pub(crate) fn finite(values: impl IntoIterator<Item = f64>) -> bool {
    values.into_iter().all(f64::is_finite)
}
