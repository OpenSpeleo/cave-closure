use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type Point = [f64; 3];

/// Externally opaque identifiers. Integer indexes are never station identities.
macro_rules! id {
    ($name:ident) => {
        #[derive(
            Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, JsonSchema,
        )]
        #[serde(transparent)]
        pub struct $name(pub Box<str>);
        impl From<&str> for $name {
            fn from(value: &str) -> Self {
                Self(value.into())
            }
        }
    };
}
id!(StationId);
id!(ObservationId);
id!(ReadingId);
id!(RelationshipId);

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    pub kind: FrameKind,
    pub length_unit: LengthUnit,
    #[serde(default)]
    pub provenance: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum FrameKind {
    LocalEnu,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum LengthUnit {
    #[serde(rename = "m")]
    Metres,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum LoopClosureMode {
    Spatial,
    Horizontal {
        preserved_z: BTreeMap<StationId, f64>,
    },
}
impl LoopClosureMode {
    pub fn dimension(&self) -> usize {
        if matches!(self, Self::Spatial) { 3 } else { 2 }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Station {
    pub id: StationId,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Uncertainty {
    Isotropic { sigma_m: f64 },
    Covariance { matrix: [[f64; 3]; 3] },
    LengthProportional { length_m: f64, variance_per_m: f64 },
    Relative { weight: f64 },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Observation {
    Displacement {
        id: ObservationId,
        #[serde(rename = "from")]
        from: StationId,
        to: StationId,
        value: Point,
        uncertainty: Uncertainty,
        #[serde(default)]
        reading_ids: Vec<ReadingId>,
    },
}
impl Observation {
    pub fn parts(
        &self,
    ) -> (
        &ObservationId,
        &StationId,
        &StationId,
        &Point,
        &Uncertainty,
        &[ReadingId],
    ) {
        let Self::Displacement {
            id,
            from,
            to,
            value,
            uncertainty,
            reading_ids,
        } = self;
        (id, from, to, value, uncertainty, reading_ids)
    }
}

#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    X,
    Y,
    Z,
}
impl Axis {
    pub fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Control {
    Hard {
        id: String,
        station: StationId,
        axes: Vec<Axis>,
        value: Vec<f64>,
    },
    Soft {
        id: String,
        station: StationId,
        axes: Vec<Axis>,
        value: Vec<f64>,
        covariance: Vec<Vec<f64>>,
    },
}
impl Control {
    pub fn parts(&self) -> (&str, &StationId, &[Axis], &[f64]) {
        match self {
            Self::Hard {
                id,
                station,
                axes,
                value,
            }
            | Self::Soft {
                id,
                station,
                axes,
                value,
                ..
            } => (id, station, axes, value),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Relationship {
    Equivalent {
        id: RelationshipId,
        #[serde(rename = "from")]
        from: StationId,
        to: StationId,
    },
    HardOffset {
        id: RelationshipId,
        #[serde(rename = "from")]
        from: StationId,
        to: StationId,
        value: Point,
    },
    Recorded {
        id: RelationshipId,
        #[serde(rename = "from")]
        from: StationId,
        to: StationId,
        source: String,
    },
}
impl Relationship {
    pub fn parts(&self) -> (&RelationshipId, &StationId, &StationId, Option<Point>) {
        match self {
            Self::Equivalent { id, from, to } => (id, from, to, Some([0.; 3])),
            Self::HardOffset {
                id,
                from,
                to,
                value,
            } => (id, from, to, Some(*value)),
            Self::Recorded { id, from, to, .. } => (id, from, to, None),
        }
    }
}

/// Numerical algorithm used for covariance-weighted loop closure.
#[cfg_attr(feature = "python", pyo3_stub_gen::derive::gen_stub_pyclass_enum)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(
        eq,
        hash,
        frozen,
        from_py_object,
        module = "cave_closure._native",
        rename_all = "SCREAMING_SNAKE_CASE"
    )
)]
#[derive(
    ts_rs::TS, Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema,
)]
#[allow(clippy::enum_variant_names)]
#[derive(Hash)]
pub enum LoopClosureAlgorithm {
    #[default]
    #[serde(rename = "LSE_REDUCED_PCG")]
    LseReducedPcg,
    #[serde(rename = "LSE_UNREDUCED_PCG")]
    LseUnreducedPcg,
    #[serde(rename = "LSE_SPARSE_QR")]
    LseSparseQr,
    #[serde(rename = "LSE_SPARSE_CHOLESKY")]
    LseSparseCholesky,
}

impl LoopClosureAlgorithm {
    pub const ALL: [Self; 4] = [
        Self::LseReducedPcg,
        Self::LseUnreducedPcg,
        Self::LseSparseQr,
        Self::LseSparseCholesky,
    ];

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::LseReducedPcg => "LSE_REDUCED_PCG",
            Self::LseUnreducedPcg => "LSE_UNREDUCED_PCG",
            Self::LseSparseQr => "LSE_SPARSE_QR",
            Self::LseSparseCholesky => "LSE_SPARSE_CHOLESKY",
        }
    }

    pub(crate) fn backend(self) -> Result<LinearBackend, crate::Error> {
        match self {
            Self::LseReducedPcg | Self::LseUnreducedPcg => Ok(LinearBackend::Pcg),
            #[cfg(feature = "sparse-direct")]
            Self::LseSparseQr => Ok(LinearBackend::SparseQr),
            #[cfg(feature = "sparse-direct")]
            Self::LseSparseCholesky => Ok(LinearBackend::SparseCholesky),
            #[cfg(not(feature = "sparse-direct"))]
            Self::LseSparseQr | Self::LseSparseCholesky => Err(crate::invalid(
                "algorithm requires the sparse-direct feature",
            )),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum LinearBackend {
    Pcg,
    #[cfg(feature = "sparse-direct")]
    SparseQr,
    #[cfg(feature = "sparse-direct")]
    SparseCholesky,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NormalizedProblem {
    #[schemars(range(min = 1, max = 1))]
    pub schema_version: u32,
    pub frame: Frame,
    pub mode: LoopClosureMode,
    pub stations: Vec<Station>,
    pub observations: Vec<Observation>,
    pub controls: Vec<Control>,
    pub relationships: Vec<Relationship>,
    pub algorithm: LoopClosureAlgorithm,
}

#[derive(Debug, Clone)]
pub struct PrepareOptions {
    pub require_control: bool,
    pub allow_shared_reading_approximation: bool,
    pub max_stations: usize,
    pub max_observations: usize,
}
impl Default for PrepareOptions {
    fn default() -> Self {
        Self {
            require_control: false,
            allow_shared_reading_approximation: false,
            max_stations: 1_000_000,
            max_observations: 4_000_000,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SolveConfig {
    pub max_iterations: usize,
    pub max_refinement_steps: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
    pub backward_error_target: f64,
}
impl Default for SolveConfig {
    fn default() -> Self {
        Self {
            max_iterations: 10_000,
            max_refinement_steps: 3,
            relative_tolerance: 1e-13,
            absolute_tolerance: 1e-12,
            backward_error_target: 1e-10,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Termination {
    Converged,
    IterationLimit,
    Stagnation,
    Cancelled,
    NumericalBreakdown,
    ResourceLimit,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Coordinate {
    pub station: StationId,
    pub value: Point,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ObservationResult {
    pub id: ObservationId,
    pub correction: Point,
    pub residual: Point,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Gauge {
    pub station: StationId,
    pub axis: Axis,
    pub value: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Diagnostics {
    pub objective: f64,
    pub cycle_rank: usize,
    pub stochastic_redundancy: usize,
    pub artificial_gauges: Vec<Gauge>,
    pub warnings: Vec<String>,
    pub calibrated_statistics_eligible: bool,
    pub variance_factor: Option<f64>,
    pub components: usize,
    pub original_unknowns: usize,
    pub kernel_unknowns: usize,
    pub iterations: usize,
    pub refinement_steps: usize,
    pub true_residual_infinity: f64,
    pub backward_error: f64,
    pub hard_constraint_max_error: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SolveReport {
    #[schemars(range(min = 1, max = 1))]
    pub schema_version: u32,
    pub status: Termination,
    pub coordinates: Vec<Coordinate>,
    pub observations: Vec<ObservationResult>,
    pub diagnostics: Diagnostics,
    pub frame: Frame,
}
impl SolveReport {
    pub fn require_converged(&self) -> Result<&Self, crate::Error> {
        if self.status == Termination::Converged {
            Ok(self)
        } else {
            Err(crate::Error::Incomplete(self.status))
        }
    }
    pub fn to_json(&self) -> Result<String, crate::Error> {
        if self
            .diagnostics
            .artificial_gauges
            .iter()
            .any(|g| !g.value.is_finite())
            || !self.diagnostics.objective.is_finite()
            || !self.diagnostics.backward_error.is_finite()
            || !self.diagnostics.true_residual_infinity.is_finite()
            || !self.diagnostics.hard_constraint_max_error.is_finite()
            || self
                .diagnostics
                .variance_factor
                .is_some_and(|v| !v.is_finite())
            || self
                .coordinates
                .iter()
                .any(|c| c.value.iter().any(|v| !v.is_finite()))
            || self.observations.iter().any(|o| {
                o.residual
                    .iter()
                    .chain(o.correction.iter())
                    .any(|v| !v.is_finite())
            })
        {
            return Err(crate::Error::Numerical(
                "non-finite result cannot be serialized".into(),
            ));
        }
        Ok(serde_json::to_string(self)?)
    }
}
