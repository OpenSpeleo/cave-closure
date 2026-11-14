use crate::*;
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
use std::{
    sync::{Arc, mpsc},
    time::Duration,
};

pyo3_stub_gen::create_exception!(cave_closure._native, ValidationError, PyValueError);
pyo3_stub_gen::create_exception!(cave_closure._native, SolveError, PyRuntimeError);
pyo3_stub_gen::create_exception!(cave_closure._native, IncompleteResultError, SolveError);

fn error(e: Error) -> PyErr {
    match e {
        Error::Validation(_) | Error::Json(_) => ValidationError::new_err(e.to_string()),
        Error::Incomplete(_) => IncompleteResultError::new_err(e.to_string()),
        _ => SolveError::new_err(e.to_string()),
    }
}
#[gen_stub_pymethods]
#[pymethods]
impl LoopClosureAlgorithm {
    #[staticmethod]
    fn variants() -> Vec<Self> {
        Self::ALL.to_vec()
    }

    #[getter]
    fn name(&self) -> &'static str {
        self.as_str()
    }
}

#[gen_stub_pyclass]
#[pyclass(frozen, module = "cave_closure._native")]
struct CancellationToken {
    inner: SolveControl,
}
#[gen_stub_pymethods]
#[pymethods]
impl CancellationToken {
    #[new]
    fn new() -> Self {
        Self {
            inner: SolveControl::default(),
        }
    }
    fn cancel(&self) {
        self.inner.cancel();
    }
    #[getter]
    fn cancelled(&self) -> bool {
        self.inner.is_cancelled()
    }
}
enum PreparedInner {
    Standard(Box<PreparedProblem>),
    L1(Box<l1::L1Prepared>),
    LoopClosurePreparation(loop_closure::LoopClosurePreparation),
}
#[gen_stub_pyclass]
#[pyclass(frozen, module = "cave_closure._native")]
struct Prepared {
    inner: Arc<PreparedInner>,
}
#[gen_stub_pymethods]
#[pymethods]
impl Prepared {
    /// Inspect the immutable numerical inputs (one entry per independent axis).
    fn normalized_json(&self) -> PyResult<Vec<String>> {
        let sources: Vec<&NormalizedProblem> = match self.inner.as_ref() {
            PreparedInner::Standard(p) => vec![p.normalized_problem()],
            PreparedInner::L1(p) => return p.input_json().map(|s| vec![s]).map_err(error),
            PreparedInner::LoopClosurePreparation(
                loop_closure::LoopClosurePreparation::Standard(p),
            ) => {
                vec![p.normalized_problem()]
            }
            PreparedInner::LoopClosurePreparation(loop_closure::LoopClosurePreparation::Axes {
                axes,
                ..
            }) => axes.iter().map(|p| p.normalized_problem()).collect(),
            PreparedInner::LoopClosurePreparation(loop_closure::LoopClosurePreparation::Graph(
                graph,
            )) => {
                let source = graph.normalized_problem().map_err(error)?;
                return serde_json::to_string(&source)
                    .map(|s| vec![s])
                    .map_err(|e| error(e.into()));
            }
        };
        sources
            .into_iter()
            .map(|p| serde_json::to_string(p).map_err(|e| error(e.into())))
            .collect()
    }
    #[pyo3(signature=(*, max_iterations=None, cancellation=None))]
    fn solve(
        &self,
        py: Python<'_>,
        max_iterations: Option<usize>,
        cancellation: Option<&CancellationToken>,
    ) -> PyResult<Report> {
        let prepared = Arc::clone(&self.inner);
        run_worker(py, cancellation, move |control| match prepared.as_ref() {
            PreparedInner::Standard(p) => p.solve(
                &SolveConfig {
                    max_iterations: max_iterations.unwrap_or(10000),
                    ..Default::default()
                },
                &mut Workspace::default(),
                &control,
            ),
            PreparedInner::LoopClosurePreparation(p) => p.solve(max_iterations, &control),
            PreparedInner::L1(p) => p.solve(max_iterations, &control),
        })
    }
}
fn run_worker(
    py: Python<'_>,
    cancellation: Option<&CancellationToken>,
    operation: impl FnOnce(SolveControl) -> Result<SolveReport, Error> + Send + 'static,
) -> PyResult<Report> {
    let control = cancellation.map_or_else(SolveControl::default, |c| c.inner.clone());
    let worker_control = control.clone();
    let (sender, mut receiver) = mpsc::sync_channel(1);
    let worker = std::thread::spawn(move || {
        let _ = sender.send(operation(worker_control));
    });
    let result = loop {
        let (returned, event) = py.detach(move || {
            let event = receiver.recv_timeout(Duration::from_millis(20));
            (receiver, event)
        });
        receiver = returned;
        match event {
            Ok(result) => break result,
            Err(mpsc::RecvTimeoutError::Timeout) => {
                if let Err(signal) = py.check_signals() {
                    control.cancel();
                    py.detach(|| {
                        let _ = worker.join();
                    });
                    return Err(signal);
                }
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                return Err(PyRuntimeError::new_err("solver worker failed"));
            }
        }
    };
    py.detach(|| worker.join())
        .map_err(|_| PyRuntimeError::new_err("solver worker panicked"))?;
    py.check_signals()?;
    result.map(|inner| Report { inner }).map_err(error)
}

fn loop_closure_input(
    input: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<loop_closure::LoopClosureInput> {
    use loop_closure::{LoopClosureInput, LoopClosureOptions};
    use pyo3::types::PyDict;
    macro_rules! field {
        ($map:expr,$name:literal) => {
            $map.get_item($name)?
                .ok_or_else(|| ValidationError::new_err(concat!("missing ", $name)))?
                .extract()?
        };
    }
    let options = input
        .get_item("options")?
        .ok_or_else(|| ValidationError::new_err("missing options"))?
        .cast_into::<PyDict>()?;
    if input.len() != 7 || options.len() != 15 {
        return Err(ValidationError::new_err(
            "unexpected loop-closure request fields",
        ));
    }
    let algorithm: String = field!(options, "algorithm");
    let algorithm = serde_json::from_value(serde_json::Value::String(algorithm))
        .map_err(|e| ValidationError::new_err(e.to_string()))?;
    Ok(LoopClosureInput {
        station_ids: field!(input, "station_ids"),
        measurements: field!(input, "measurements"),
        initial_positions: field!(input, "initial_positions"),
        controls: field!(input, "controls"),
        relationships: field!(input, "relationships"),
        source_order: field!(input, "source_order"),
        options: LoopClosureOptions {
            mode: field!(options, "mode"),
            algorithm,
            weighting: field!(options, "weighting"),
            method: field!(options, "method"),
            sigma_m: field!(options, "sigma_m"),
            weight_cap: field!(options, "weight_cap"),
            sigma_length_m: field!(options, "sigma_length_m"),
            sigma_azimuth_deg: field!(options, "sigma_azimuth_deg"),
            average_azimuths: field!(options, "average_azimuths"),
            anchor_floating: field!(options, "anchor_floating"),
            preserve_centroids: field!(options, "preserve_centroids"),
            fix_control_depths: field!(options, "fix_control_depths"),
            max_iterations: field!(options, "max_iterations"),
            absolute_tolerance: field!(options, "absolute_tolerance"),
            relative_tolerance: field!(options, "relative_tolerance"),
        },
    })
}
fn l1_input(input: &Bound<'_, pyo3::types::PyDict>) -> PyResult<l1::L1Input> {
    use pyo3::types::PyDict;
    macro_rules! field {
        ($map:expr,$name:literal) => {
            $map.get_item($name)?
                .ok_or_else(|| ValidationError::new_err(concat!("missing ", $name)))?
                .extract()?
        };
    }
    let options = input
        .get_item("options")?
        .ok_or_else(|| ValidationError::new_err("missing options"))?
        .cast_into::<PyDict>()?;
    if input.len() != 5 || options.len() != 5 {
        return Err(ValidationError::new_err("unexpected L1 request fields"));
    }
    Ok(l1::L1Input {
        station_ids: field!(input, "station_ids"),
        observations: field!(input, "observations"),
        initial_positions: field!(input, "initial_positions"),
        anchor_indices: field!(input, "anchor_indices"),
        options: l1::L1Options {
            mode: field!(options, "mode"),
            length_floor_m: field!(options, "length_floor_m"),
            minimum_anchor_count: field!(options, "minimum_anchor_count"),
            deduplicate_edges: field!(options, "deduplicate_edges"),
            max_iterations: field!(options, "max_iterations"),
        },
    })
}

/// Copy and validate an owned inverse-length L1 graph.
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
fn prepare_l1_loop_closure(
    py: Python<'_>,
    input: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<Prepared> {
    let input = l1_input(input)?;
    let prepared = py
        .detach(move || l1::prepare_l1_loop_closure(input))
        .map_err(error)?;
    Ok(Prepared {
        inner: Arc::new(PreparedInner::L1(Box::new(prepared))),
    })
}
/// One crossing for validation, LP assembly, HiGHS solves and complete output.
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
#[pyo3(signature=(input,*,max_iterations=None,cancellation=None))]
fn solve_l1_loop_closure_json(
    py: Python<'_>,
    input: &Bound<'_, pyo3::types::PyDict>,
    max_iterations: Option<usize>,
    cancellation: Option<&CancellationToken>,
) -> PyResult<String> {
    let input = l1_input(input)?;
    let report = run_worker(py, cancellation, move |control| {
        l1::prepare_l1_loop_closure(input)?.solve(max_iterations, &control)
    })?;
    py.detach(move || report.inner.to_json()).map_err(error)
}

/// Copy a loop-closure request and prepare all numerical factors in Rust.
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
fn prepare_loop_closure(
    py: Python<'_>,
    input: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<Prepared> {
    let input = loop_closure_input(input)?;
    let prepared = py
        .detach(move || loop_closure::prepare_loop_closure(input, &SolveControl::default()))
        .map_err(error)?;
    Ok(Prepared {
        inner: Arc::new(PreparedInner::LoopClosurePreparation(prepared)),
    })
}
/// One crossing: copy, prepare, solve and export the complete result in Rust.
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
#[pyo3(signature=(input, *, max_iterations=None, cancellation=None))]
fn solve_loop_closure_json(
    py: Python<'_>,
    input: &Bound<'_, pyo3::types::PyDict>,
    max_iterations: Option<usize>,
    cancellation: Option<&CancellationToken>,
) -> PyResult<String> {
    let input = loop_closure_input(input)?;
    let report = run_worker(py, cancellation, move |control| {
        loop_closure::prepare_loop_closure(input, &control)?.solve(max_iterations, &control)
    })?;
    py.detach(move || report.inner.to_json()).map_err(error)
}

#[gen_stub_pyclass]
#[pyclass(frozen, module = "cave_closure._native")]
struct Report {
    inner: SolveReport,
}
#[gen_stub_pymethods]
#[pymethods]
impl Report {
    #[getter]
    fn status(&self) -> String {
        format!("{:?}", self.inner.status)
    }
    #[getter]
    fn coordinates(&self) -> Vec<(String, [f64; 3])> {
        self.inner
            .coordinates
            .iter()
            .map(|c| (c.station.0.to_string(), c.value))
            .collect()
    }
    #[getter]
    fn observations(&self) -> Vec<(String, [f64; 3])> {
        self.inner
            .observations
            .iter()
            .map(|o| (o.id.0.to_string(), o.residual))
            .collect()
    }
    fn require_converged(&self) -> PyResult<()> {
        self.inner.require_converged().map(|_| ()).map_err(error)
    }
    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(error)
    }
}
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
#[pyo3(signature=(input, *, algorithm=None, require_control=false, allow_shared_reading_approximation=false, max_stations=1000000, max_observations=4000000))]
fn prepare_json(
    input: &str,
    algorithm: Option<LoopClosureAlgorithm>,
    require_control: bool,
    allow_shared_reading_approximation: bool,
    max_stations: usize,
    max_observations: usize,
) -> PyResult<Prepared> {
    let mut p = json::from_json(input).map_err(error)?;
    if let Some(algorithm) = algorithm {
        p.algorithm = algorithm;
    }
    prepare(
        &p,
        &PrepareOptions {
            require_control,
            allow_shared_reading_approximation,
            max_stations,
            max_observations,
        },
    )
    .map(|inner| Prepared {
        inner: Arc::new(PreparedInner::Standard(Box::new(inner))),
    })
    .map_err(error)
}
type ControlRow = (String, usize, Vec<usize>, Vec<f64>, Option<Vec<Vec<f64>>>);
type RelationshipRow = (String, String, usize, usize, [f64; 3]);
#[gen_stub_pyfunction(module = "cave_closure._native")]
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature=(station_ids, observation_ids, endpoints, displacements, covariances, *, mode, algorithm=LoopClosureAlgorithm::LseReducedPcg, preserved_z=None, controls=vec![], relationships=vec![], reading_ids=None, provenance=vec![], require_control=false, allow_shared_reading_approximation=false))]
fn prepare_arrays(
    station_ids: Vec<String>,
    observation_ids: Vec<String>,
    endpoints: Vec<[usize; 2]>,
    displacements: Vec<[f64; 3]>,
    covariances: Vec<[[f64; 3]; 3]>,
    mode: &str,
    algorithm: LoopClosureAlgorithm,
    preserved_z: Option<Vec<f64>>,
    controls: Vec<ControlRow>,
    relationships: Vec<RelationshipRow>,
    reading_ids: Option<Vec<Vec<String>>>,
    provenance: Vec<String>,
    require_control: bool,
    allow_shared_reading_approximation: bool,
) -> PyResult<Prepared> {
    let n = station_ids.len();
    let m = observation_ids.len();
    if n > 1_000_000
        || m > 4_000_000
        || endpoints.len() != m
        || displacements.len() != m
        || covariances.len() != m
        || reading_ids.as_ref().is_some_and(|r| r.len() != m)
    {
        return Err(ValidationError::new_err(
            "inconsistent array lengths or configured input resource limit exceeded",
        ));
    }
    let ids: Vec<StationId> = station_ids
        .iter()
        .map(|s| StationId(s.clone().into()))
        .collect();
    let station = |i: usize| {
        ids.get(i)
            .cloned()
            .ok_or_else(|| ValidationError::new_err("station index out of bounds"))
    };
    let mode = match (mode, preserved_z) {
        ("spatial", None) => LoopClosureMode::Spatial,
        ("horizontal", Some(z)) if z.len() == n => LoopClosureMode::Horizontal {
            preserved_z: ids.iter().cloned().zip(z).collect(),
        },
        _ => {
            return Err(ValidationError::new_err(
                "mode must be spatial or horizontal with a complete preserved_z array",
            ));
        }
    };
    let observations = (0..m)
        .map(|i| {
            Ok(Observation::Displacement {
                id: ObservationId(observation_ids[i].clone().into()),
                from: station(endpoints[i][0])?,
                to: station(endpoints[i][1])?,
                value: displacements[i],
                uncertainty: Uncertainty::Covariance {
                    matrix: covariances[i],
                },
                reading_ids: reading_ids.as_ref().map_or_else(Vec::new, |r| {
                    r[i].iter().map(|s| ReadingId(s.clone().into())).collect()
                }),
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let controls = controls
        .into_iter()
        .map(|(id, i, axes, value, covariance)| {
            let axes = axes
                .into_iter()
                .map(|a| {
                    [Axis::X, Axis::Y, Axis::Z]
                        .get(a)
                        .copied()
                        .ok_or_else(|| ValidationError::new_err("invalid axis"))
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(match covariance {
                Some(covariance) => Control::Soft {
                    id,
                    station: station(i)?,
                    axes,
                    value,
                    covariance,
                },
                None => Control::Hard {
                    id,
                    station: station(i)?,
                    axes,
                    value,
                },
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let relationships = relationships
        .into_iter()
        .map(|(id, kind, u, v, value)| {
            let (id, from, to) = (RelationshipId(id.into()), station(u)?, station(v)?);
            Ok(match kind.as_str() {
                "equivalent" => Relationship::Equivalent { id, from, to },
                "hard_offset" => Relationship::HardOffset {
                    id,
                    from,
                    to,
                    value,
                },
                "recorded" => Relationship::Recorded {
                    id,
                    from,
                    to,
                    source: "adapter".into(),
                },
                _ => return Err(ValidationError::new_err("unknown relationship kind")),
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let p = NormalizedProblem {
        schema_version: 1,
        frame: Frame {
            kind: FrameKind::LocalEnu,
            length_unit: LengthUnit::Metres,
            provenance,
        },
        mode,
        stations: ids.into_iter().map(|id| Station { id }).collect(),
        observations,
        controls,
        relationships,
        algorithm,
    };
    prepare(
        &p,
        &PrepareOptions {
            require_control,
            allow_shared_reading_approximation,
            ..Default::default()
        },
    )
    .map(|inner| Prepared {
        inner: Arc::new(PreparedInner::Standard(Box::new(inner))),
    })
    .map_err(error)
}
#[pymodule(gil_used = false)]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("ValidationError", m.py().get_type::<ValidationError>())?;
    m.add("SolveError", m.py().get_type::<SolveError>())?;
    m.add(
        "IncompleteResultError",
        m.py().get_type::<IncompleteResultError>(),
    )?;
    m.add_class::<LoopClosureAlgorithm>()?;
    m.add_class::<Prepared>()?;
    m.add_class::<Report>()?;
    m.add_class::<CancellationToken>()?;
    m.add_function(wrap_pyfunction!(prepare_json, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_arrays, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_loop_closure, m)?)?;
    m.add_function(wrap_pyfunction!(solve_loop_closure_json, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_l1_loop_closure, m)?)?;
    m.add_function(wrap_pyfunction!(solve_l1_loop_closure_json, m)?)?;
    Ok(())
}

/// Gather native API metadata using the package manifest.
#[cfg(feature = "generate-stubs")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    let mut info = pyo3_stub_gen::StubInfo::from_pyproject_toml(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("pyproject.toml"),
    )?;
    // stub-gen treats custom exception bases as builtins. This base belongs to
    // our native module, so retain the actual runtime exception hierarchy.
    let module = info.modules.get_mut("cave_closure._native").unwrap();
    let exception = module
        .class
        .get_mut(&std::any::TypeId::of::<IncompleteResultError>())
        .unwrap();
    exception.bases = vec![pyo3_stub_gen::TypeInfo::locally_defined(
        "SolveError",
        "cave_closure._native".into(),
    )];
    Ok(info)
}
