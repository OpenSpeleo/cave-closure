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
#[gen_stub_pyclass]
#[pyclass(frozen, module = "cave_closure._native")]
struct Prepared {
    inner: Arc<PreparedProblem>,
}
#[gen_stub_pymethods]
#[pymethods]
impl Prepared {
    #[pyo3(signature=(*, max_iterations=10000, cancellation=None))]
    fn solve(
        &self,
        py: Python<'_>,
        max_iterations: usize,
        cancellation: Option<&CancellationToken>,
    ) -> PyResult<Report> {
        let prepared = Arc::clone(&self.inner);
        let control = cancellation.map_or_else(SolveControl::default, |c| c.inner.clone());
        let worker_control = control.clone();
        let (sender, mut receiver) = mpsc::sync_channel(1);
        let worker = std::thread::spawn(move || {
            let result = prepared.solve(
                &SolveConfig {
                    max_iterations,
                    ..Default::default()
                },
                &mut Workspace::default(),
                &worker_control,
            );
            let _ = sender.send(result);
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
        inner: Arc::new(inner),
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
        inner: Arc::new(inner),
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
