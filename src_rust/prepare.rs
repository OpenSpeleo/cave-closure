use crate::*;
use nalgebra::{DMatrix, DVector};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Clone)]
pub(crate) struct Factor {
    pub id: ObservationId,
    pub u: usize,
    pub v: usize,
    pub measured: DVector<f64>,
    pub precision: DMatrix<f64>,
}
#[derive(Clone)]
pub(crate) struct Soft {
    pub root: usize,
    pub axes: Vec<usize>,
    pub value: DVector<f64>,
    pub precision: DMatrix<f64>,
}
/// Immutable, validated owned input; share across concurrent solves with distinct workspaces.
#[derive(Clone)]
pub struct PreparedProblem {
    pub(crate) source: NormalizedProblem,
    pub(crate) ids: Vec<StationId>,
    pub(crate) roots: Vec<usize>,
    pub(crate) offsets: Vec<DVector<f64>>,
    pub(crate) fixed: Vec<Vec<Option<f64>>>,
    pub(crate) factors: Vec<Factor>,
    pub(crate) soft: Vec<Soft>,
    pub(crate) gauges: Vec<Gauge>,
    pub(crate) warnings: Vec<String>,
    pub(crate) components: usize,
    pub(crate) cycle_rank: usize,
    pub(crate) rows: usize,
    pub(crate) calibrated: bool,
}

impl PreparedProblem {
    /// Immutable normalized source evidence owned by this preparation.
    pub fn normalized_problem(&self) -> &NormalizedProblem {
        &self.source
    }
}

pub(crate) fn precision(cov: DMatrix<f64>, context: &str) -> Result<DMatrix<f64>, Error> {
    if cov.nrows() == 0 || cov.nrows() != cov.ncols() || !finite(cov.iter().copied()) {
        return Err(invalid(format!(
            "{context}: expected finite square covariance"
        )));
    }
    for i in 0..cov.nrows() {
        for j in 0..i {
            if cov[(i, j)] != cov[(j, i)] {
                return Err(invalid(format!("{context}: covariance is not symmetric")));
            }
        }
    }
    let chol = cov
        .cholesky()
        .ok_or_else(|| invalid(format!("{context}: covariance is not positive definite")))?;
    let result = chol.inverse();
    if !finite(result.iter().copied()) {
        return Err(invalid(format!("{context}: precision overflow")));
    }
    Ok(result)
}
fn covariance(uncertainty: &Uncertainty) -> Result<DMatrix<f64>, Error> {
    let variance = match uncertainty {
        Uncertainty::Isotropic { sigma_m } if sigma_m.is_finite() && *sigma_m > 0. => {
            sigma_m * sigma_m
        }
        Uncertainty::LengthProportional {
            length_m,
            variance_per_m,
        } if length_m.is_finite()
            && *length_m > 0.
            && variance_per_m.is_finite()
            && *variance_per_m > 0. =>
        {
            length_m * variance_per_m
        }
        Uncertainty::Relative { weight } if weight.is_finite() && *weight > 0. => 1. / weight,
        Uncertainty::Covariance { matrix } => {
            let cov = DMatrix::from_fn(3, 3, |i, j| matrix[i][j]);
            precision(cov.clone(), "spatial covariance")?;
            return Ok(cov);
        }
        _ => {
            return Err(invalid(
                "uncertainty requires positive finite parameters; zero-length proportional variance is undefined",
            ));
        }
    };
    if !variance.is_finite() || variance <= 0. {
        return Err(invalid("variance overflow or underflow"));
    }
    Ok(DMatrix::identity(3, 3) * variance)
}
fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 128. * f64::EPSILON * (a.abs() + b.abs() + scale).max(1.)
}
/// Validate and own a problem, preserving its independent observations.
///
/// ```
/// use cave_closure::{json, prepare, PrepareOptions, SolveConfig, SolveControl, Workspace};
/// let input = include_str!("../tests/fixtures/controlled_tree.json");
/// let problem = json::from_json(input)?;
/// let prepared = prepare(&problem, &PrepareOptions::default())?;
/// let report = prepared.solve(&SolveConfig::default(), &mut Workspace::default(), &SolveControl::default())?;
/// report.require_converged()?;
/// # Ok::<(), cave_closure::Error>(())
/// ```
pub fn prepare(
    problem: &NormalizedProblem,
    options: &PrepareOptions,
) -> Result<PreparedProblem, Error> {
    if problem.schema_version != 1 {
        return Err(invalid("unsupported schema version"));
    }
    if problem.stations.len() > options.max_stations
        || problem.observations.len() > options.max_observations
    {
        return Err(invalid("configured input resource limit exceeded"));
    }
    problem.algorithm.backend()?;
    let mut source = problem.clone();
    source.stations.sort_by(|a, b| a.id.cmp(&b.id));
    source
        .observations
        .sort_by(|a, b| a.parts().0.cmp(b.parts().0));
    source
        .relationships
        .sort_by(|a, b| a.parts().0.cmp(b.parts().0));
    source.controls.sort_by(|a, b| a.parts().0.cmp(b.parts().0));
    let ids: Vec<_> = source.stations.iter().map(|s| s.id.clone()).collect();
    let index: BTreeMap<_, _> = ids
        .iter()
        .enumerate()
        .map(|(i, id)| (id.clone(), i))
        .collect();
    if ids.len() != index.len() || ids.iter().any(|id| id.0.is_empty()) {
        return Err(invalid("empty or duplicate station identity"));
    }
    let locate = |id: &StationId| {
        index
            .get(id)
            .copied()
            .ok_or_else(|| invalid(format!("unknown station {}", id.0)))
    };
    let n = ids.len();
    let d = source.mode.dimension();
    if let LoopClosureMode::Horizontal { preserved_z } = &source.mode
        && (preserved_z.len() != n
            || ids.iter().any(|id| !preserved_z.contains_key(id))
            || !finite(preserved_z.values().copied()))
    {
        return Err(invalid(
            "horizontal mode requires exactly one finite preserved Z per station",
        ));
    }
    let mut exact = vec![Vec::new(); n];
    let mut relationship_ids = BTreeSet::new();
    for rel in &source.relationships {
        let (id, u, v, value) = rel.parts();
        if id.0.is_empty() || !relationship_ids.insert(id) {
            return Err(invalid("empty or duplicate relationship identity"));
        }
        let (u, v) = (locate(u)?, locate(v)?);
        if let Some(value) = value {
            if !finite(value) {
                return Err(invalid(format!("relationship {}: non-finite offset", id.0)));
            }
            let delta = DVector::from_column_slice(&value[..d]);
            exact[u].push((v, delta.clone()));
            exact[v].push((u, -delta));
        }
    }
    let mut roots = vec![usize::MAX; n];
    let mut offsets = vec![DVector::zeros(d); n];
    let mut path_scales = vec![0.; n];
    let mut representatives = Vec::new();
    for start in 0..n {
        if roots[start] != usize::MAX {
            continue;
        }
        let root = representatives.len();
        representatives.push(start);
        roots[start] = root;
        let mut queue = VecDeque::from([start]);
        while let Some(u) = queue.pop_front() {
            for (v, delta) in &exact[u] {
                let value = &offsets[u] + delta;
                let scale = path_scales[u] + delta.amax();
                if !finite(value.iter().copied()) || !scale.is_finite() {
                    return Err(invalid("hard offset accumulation overflow"));
                }
                if roots[*v] == usize::MAX {
                    roots[*v] = root;
                    offsets[*v] = value;
                    path_scales[*v] = scale;
                    queue.push_back(*v);
                } else if (0..d).any(|a| !close(offsets[*v][a], value[a], scale + path_scales[*v]))
                {
                    return Err(invalid("inconsistent cycle of exact offsets"));
                }
            }
        }
    }
    let nr = representatives.len();
    let mut fixed = vec![vec![None; d]; nr];
    let mut identified = vec![vec![false; d]; nr];
    let mut control_ids = BTreeSet::new();
    let mut soft = Vec::new();
    for control in &source.controls {
        let (id, station, axes, value) = control.parts();
        let i = locate(station)?;
        if id.is_empty()
            || !control_ids.insert(id)
            || axes.is_empty()
            || axes.len() != value.len()
            || !finite(value.iter().copied())
        {
            return Err(invalid(format!(
                "control {id}: invalid identity, axes or values"
            )));
        }
        let axes: Vec<_> = axes.iter().map(|a| a.index()).collect();
        if axes.iter().any(|a| *a >= d) || axes.iter().collect::<BTreeSet<_>>().len() != axes.len()
        {
            return Err(invalid(format!("control {id}: repeated or inactive axis")));
        }
        let transformed = DVector::from_iterator(
            axes.len(),
            axes.iter().zip(value).map(|(a, v)| v - offsets[i][*a]),
        );
        match control {
            Control::Hard { .. } => {
                for (k, a) in axes.iter().enumerate() {
                    let slot = &mut fixed[roots[i]][*a];
                    if slot.is_some_and(|old| !close(old, transformed[k], path_scales[i])) {
                        return Err(invalid(format!("control {id}: incompatible hard controls")));
                    }
                    *slot = Some(transformed[k]);
                }
            }
            Control::Soft { covariance, .. } => {
                let m = axes.len();
                if covariance.len() != m || covariance.iter().any(|r| r.len() != m) {
                    return Err(invalid(format!("control {id}: covariance shape")));
                }
                let precision = precision(DMatrix::from_fn(m, m, |a, b| covariance[a][b]), id)?;
                soft.push(Soft {
                    root: roots[i],
                    axes: axes.clone(),
                    value: transformed,
                    precision,
                });
            }
        }
        for a in axes {
            identified[roots[i]][a] = true;
        }
    }
    let mut adj = vec![BTreeSet::new(); nr];
    let mut factors = Vec::new();
    let mut observation_ids = BTreeSet::new();
    let mut reading_ids = BTreeSet::new();
    let mut reused = false;
    let mut calibrated = true;
    for observation in &source.observations {
        let (id, u, v, value, uncertainty, readings) = observation.parts();
        let (u, v) = (locate(u)?, locate(v)?);
        if id.0.is_empty() || !observation_ids.insert(id) || !finite(*value) {
            return Err(invalid(
                "empty/duplicate observation identity or non-finite displacement",
            ));
        }
        let mut local = BTreeSet::new();
        for reading in readings {
            if reading.0.is_empty() || !local.insert(reading) {
                return Err(invalid(format!(
                    "observation {}: duplicate/empty reading",
                    id.0
                )));
            }
            if !reading_ids.insert(reading) {
                reused = true;
            }
        }
        let cov = covariance(uncertainty)?;
        // Marginal XY covariance, never a submatrix of spatial precision.
        let precision = precision(cov.view((0, 0), (d, d)).into_owned(), &id.0)?;
        calibrated &= !matches!(uncertainty, Uncertainty::Relative { .. });
        if roots[u] != roots[v] {
            adj[roots[u]].insert(roots[v]);
            adj[roots[v]].insert(roots[u]);
        }
        factors.push(Factor {
            id: id.clone(),
            u,
            v,
            measured: DVector::from_column_slice(&value[..d]),
            precision,
        });
    }
    if reused && !options.allow_shared_reading_approximation {
        return Err(invalid(
            "shared raw readings require correlation support or explicit independent-factor approximation",
        ));
    }
    let mut warnings = Vec::new();
    if reused {
        warnings
            .push("Shared readings treated as independent by explicit approximation policy".into());
        calibrated = false;
    }
    let mut gauges = Vec::new();
    let mut seen = vec![false; nr];
    let mut components = 0;
    for start in 0..nr {
        if seen[start] {
            continue;
        }
        components += 1;
        seen[start] = true;
        let mut members = vec![start];
        let mut cursor = 0;
        while cursor < members.len() {
            let u = members[cursor];
            cursor += 1;
            for v in &adj[u] {
                if !seen[*v] {
                    seen[*v] = true;
                    members.push(*v);
                }
            }
        }
        for a in 0..d {
            if members.iter().any(|u| identified[*u][a]) {
                continue;
            }
            if options.require_control {
                return Err(invalid(format!("uncontrolled translation direction {a}")));
            }
            fixed[start][a] = Some(0.);
            gauges.push(Gauge {
                station: ids[representatives[start]].clone(),
                axis: [Axis::X, Axis::Y, Axis::Z][a],
                value: 0.,
            });
        }
    }
    let rows = factors
        .len()
        .checked_mul(d)
        .and_then(|m| m.checked_add(soft.iter().map(|s| s.axes.len()).sum()))
        .ok_or_else(|| invalid("stochastic row count overflow"))?;
    let edges = source.observations.len()
        + source
            .relationships
            .iter()
            .filter(|r| r.parts().3.is_some())
            .count();
    let cycle_rank = edges.saturating_add(components).saturating_sub(n);
    Ok(PreparedProblem {
        source,
        ids,
        roots,
        offsets,
        fixed,
        factors,
        soft,
        gauges,
        warnings,
        components,
        cycle_rank,
        rows,
        calibrated,
    })
}
