//! Inverse-length weighted least-absolute-deviations loop closure.
//!
//! The model is assembled and solved here, with HiGHS through its C interface.
//! No Python optimizer or consumer survey types participate in a computation.
use crate::*;
use std::collections::BTreeSet;

pub type L1Observation = (String, usize, usize, Option<Point>);

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct L1Options {
    pub mode: String,
    pub length_floor_m: f64,
    pub minimum_anchor_count: usize,
    pub deduplicate_edges: bool,
    pub max_iterations: Option<usize>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct L1Input {
    pub station_ids: Vec<String>,
    pub observations: Vec<L1Observation>,
    pub initial_positions: Vec<Point>,
    pub anchor_indices: Vec<usize>,
    pub options: L1Options,
}

/// Owned, immutable evidence; each solve creates independent LP workspaces.
#[derive(Clone)]
pub struct L1Prepared {
    input: L1Input,
}

fn pair(u: usize, v: usize) -> (usize, usize) {
    (u.min(v), u.max(v))
}
fn source_length(value: Point) -> f64 {
    // Vector3D.length in Compass uses left-associated libm powers, unlike the
    // compensated Python sum used by the historical CG model.
    let square = |x: f64| x.powf(std::hint::black_box(2.0));
    (square(value[0]) + square(value[1]) + square(value[2])).powf(std::hint::black_box(0.5))
}
fn source_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut high, mut low) = (0.0_f64, 0.0_f64);
    for x in values {
        let next = high + x;
        low += if high.abs() >= x.abs() {
            (high - next) + x
        } else {
            (x - next) + high
        };
        high = next;
    }
    high + low
}

pub fn prepare_l1_loop_closure(mut input: L1Input) -> Result<L1Prepared, Error> {
    let n = input.station_ids.len();
    if n == 0 || n > 1_000_000 || input.observations.len() > 4_000_000 {
        return Err(invalid("configured L1 input resource limit exceeded"));
    }
    if input.initial_positions.len() != n
        || !finite(input.initial_positions.iter().flatten().copied())
    {
        return Err(invalid("invalid L1 initial positions"));
    }
    if input.station_ids.iter().any(|id| id.is_empty())
        || input.station_ids.iter().collect::<BTreeSet<_>>().len() != n
    {
        return Err(invalid("invalid L1 station identities"));
    }
    if input.anchor_indices.iter().any(|&i| i >= n)
        || input.anchor_indices.iter().collect::<BTreeSet<_>>().len() != input.anchor_indices.len()
    {
        return Err(invalid("invalid L1 anchor indexes"));
    }
    if !["horizontal", "spatial"].contains(&input.options.mode.as_str())
        || !input.options.length_floor_m.is_finite()
        || input.options.length_floor_m <= 0.0
        || input.options.minimum_anchor_count == 0
    {
        return Err(invalid("invalid L1 options"));
    }
    if input
        .options
        .max_iterations
        .is_some_and(|limit| limit > i32::MAX as usize)
    {
        return Err(invalid("L1 iteration limit exceeds backend integer range"));
    }
    let mut identities = BTreeSet::new();
    let mut pairs = BTreeSet::new();
    let mut rows = Vec::new();
    for (id, u, v, delta) in &input.observations {
        if id.is_empty() || !identities.insert(id.clone()) || *u >= n || *v >= n {
            return Err(invalid("invalid L1 observation identity or endpoint"));
        }
        let delta = delta.unwrap_or_else(|| {
            std::array::from_fn(|a| input.initial_positions[*v][a] - input.initial_positions[*u][a])
        });
        if !finite(delta) || !source_length(delta).is_finite() {
            return Err(invalid("nonfinite L1 measurement"));
        }
        if !input.options.deduplicate_edges || pairs.insert(pair(*u, *v)) {
            rows.push((id.clone(), *u, *v, Some(delta)));
        }
    }
    input.observations = rows;
    Ok(L1Prepared { input })
}

impl L1Prepared {
    pub fn input_json(&self) -> Result<String, Error> {
        Ok(serde_json::to_string(&self.input)?)
    }
    pub fn solve(
        &self,
        max_iterations: Option<usize>,
        control: &SolveControl,
    ) -> Result<SolveReport, Error> {
        if control.is_cancelled() {
            return Err(Error::Incomplete(Termination::Cancelled));
        }
        let input = &self.input;
        let anchors: BTreeSet<_> = input.anchor_indices.iter().copied().collect();
        let mut free: Vec<_> = (0..input.station_ids.len())
            .filter(|i| !anchors.contains(i))
            .collect();
        free.sort_by_key(|&i| &input.station_ids[i]);
        let dimension = if input.options.mode == "spatial" {
            3
        } else {
            2
        };
        let mut positions = input.initial_positions.clone();
        let mut warnings = Vec::new();
        let mut iterations = 0;
        if anchors.len() >= input.options.minimum_anchor_count && !free.is_empty() {
            let mut ordered: Vec<_> = anchors.iter().copied().collect();
            ordered.sort_by_key(|&i| &input.station_ids[i]);
            let origin: Point = std::array::from_fn(|a| {
                source_sum(ordered.iter().map(|&i| positions[i][a])) / ordered.len() as f64
            });
            let mut index = vec![None; positions.len()];
            for (i, &station) in free.iter().enumerate() {
                index[station] = Some(i);
            }
            let rows: Vec<_> = input
                .observations
                .iter()
                .filter(|(_, u, v, _)| !(anchors.contains(u) && anchors.contains(v)))
                .collect();
            if !rows.is_empty() {
                let n = free.len();
                let m = rows.len();
                let mut columns = vec![Vec::<(usize, f64)>::new(); n + 2 * m];
                let mut costs = vec![0.0; n + 2 * m];
                let mut rhs = [vec![0.0; m], vec![0.0; m], vec![0.0; m]];
                for (r, (_, u, v, delta)) in rows.iter().enumerate() {
                    let delta = delta.expect("prepared measurements are populated");
                    let weight = 1.0 / source_length(delta).max(input.options.length_floor_m);
                    match (index[*u], index[*v]) {
                        (Some(left), Some(right)) => {
                            // Retain the source dense assignment's last-write
                            // behavior for a self-observation.
                            if left != right {
                                columns[left].push((r, -1.0));
                            }
                            columns[right].push((r, 1.0));
                            for a in 0..3 {
                                rhs[a][r] = delta[a];
                            }
                        }
                        (None, Some(right)) => {
                            columns[right].push((r, 1.0));
                            for a in 0..3 {
                                rhs[a][r] = positions[*u][a] - origin[a] + delta[a];
                            }
                        }
                        (Some(left), None) => {
                            columns[left].push((r, -1.0));
                            for a in 0..3 {
                                rhs[a][r] = delta[a] - (positions[*v][a] - origin[a]);
                            }
                        }
                        _ => unreachable!(),
                    }
                    columns[n + r].push((r, -1.0));
                    columns[n + m + r].push((r, 1.0));
                    costs[n + r] = weight;
                    costs[n + m + r] = weight;
                }
                let limit = max_iterations.or(input.options.max_iterations);
                if limit.is_some_and(|value| value > i32::MAX as usize) {
                    return Err(invalid("L1 iteration limit exceeds backend integer range"));
                }
                for a in 0..dimension {
                    let (solution, used) =
                        solve_axis(&columns, &costs, &rhs[a], n, limit, control)?;
                    iterations += used;
                    for (i, &station) in free.iter().enumerate() {
                        positions[station][a] = solution[i] + origin[a];
                    }
                }
            }
        } else if anchors.len() < input.options.minimum_anchor_count {
            warnings.push("Fewer than the required anchors; initial coordinates retained".into());
        }
        let coordinates = positions
            .iter()
            .enumerate()
            .map(|(i, &value)| Coordinate {
                station: StationId(input.station_ids[i].clone().into()),
                value,
            })
            .collect();
        let mut objective = 0.0;
        let observations = input
            .observations
            .iter()
            .map(|(id, u, v, delta)| {
                let delta = delta.expect("prepared measurements are populated");
                let residual: Point =
                    std::array::from_fn(|a| positions[*v][a] - positions[*u][a] - delta[a]);
                objective += source_sum(residual[..dimension].iter().map(|v| v.abs()))
                    / source_length(delta).max(input.options.length_floor_m);
                ObservationResult {
                    id: ObservationId(id.clone().into()),
                    correction: residual,
                    residual,
                }
            })
            .collect();
        warnings.push(
            "L1 relative objective; Gaussian variance and redundancy statistics are not reported"
                .into(),
        );
        let report = SolveReport {
            schema_version: 1,
            status: Termination::Converged,
            coordinates,
            observations,
            frame: Frame {
                kind: FrameKind::LocalEnu,
                length_unit: LengthUnit::Metres,
                provenance: vec![
                    "Inverse-length weighted L1 loop closure, HiGHS dual simplex".into(),
                ],
            },
            diagnostics: Diagnostics {
                objective,
                cycle_rank: 0,
                stochastic_redundancy: 0,
                artificial_gauges: vec![],
                warnings,
                calibrated_statistics_eligible: false,
                variance_factor: None,
                components: 0,
                original_unknowns: free.len() * dimension,
                kernel_unknowns: free.len() * dimension,
                iterations,
                refinement_steps: 0,
                true_residual_infinity: 0.0,
                backward_error: 0.0,
                hard_constraint_max_error: 0.0,
            },
        };
        report.to_json()?;
        Ok(report)
    }
}

#[cfg(not(feature = "l1"))]
fn solve_axis(
    _: &[Vec<(usize, f64)>],
    _: &[f64],
    _: &[f64],
    _: usize,
    _: Option<usize>,
    _: &SolveControl,
) -> Result<(Vec<f64>, usize), Error> {
    Err(invalid("L1 loop closure requires the l1 feature"))
}

#[cfg(feature = "l1")]
fn solve_axis(
    columns: &[Vec<(usize, f64)>],
    costs: &[f64],
    rhs: &[f64],
    free: usize,
    limit: Option<usize>,
    control: &SolveControl,
) -> Result<(Vec<f64>, usize), Error> {
    use highs_sys::*;
    use std::{
        ffi::{c_char, c_void},
        ptr::{self, NonNull},
    };
    struct Model(NonNull<c_void>);
    impl Drop for Model {
        fn drop(&mut self) {
            unsafe {
                Highs_destroy(self.0.as_ptr());
            }
        }
    }
    fn check(status: HighsInt) -> Result<(), Error> {
        if status == STATUS_ERROR {
            Err(Error::Numerical(
                "HiGHS rejected the L1 model or options".into(),
            ))
        } else {
            Ok(())
        }
    }
    unsafe extern "C" fn interrupt(
        _: i32,
        _: *const c_char,
        _: *const HighsCallbackDataOut,
        data: *mut HighsCallbackDataIn,
        user: *mut c_void,
    ) {
        if !data.is_null() && !user.is_null() {
            // The solve owns the token for the entire synchronous HiGHS run.
            let token = unsafe { &*user.cast::<SolveControl>() };
            if token.is_cancelled() {
                unsafe {
                    (*data).user_interrupt = 1;
                }
            }
        }
    }
    if control.is_cancelled() {
        return Err(Error::Incomplete(Termination::Cancelled));
    }
    let mut starts = Vec::with_capacity(columns.len() + 1);
    let mut indexes = Vec::new();
    let mut values = Vec::new();
    for column in columns {
        starts.push(values.len() as HighsInt);
        for &(row, value) in column {
            indexes.push(row as HighsInt);
            values.push(value);
        }
    }
    starts.push(values.len() as HighsInt);
    let mut lower = vec![0.0; columns.len()];
    lower[..free].fill(-1e30);
    let upper = vec![1e30; columns.len()];
    let mut solution = vec![0.0; columns.len()];
    // All buffers are owned, immutable inputs (or uniquely owned output), and
    // remain alive until the model is destroyed. Counts are resource-bounded.
    unsafe {
        let model = Model(
            NonNull::new(Highs_create())
                .ok_or_else(|| Error::Numerical("HiGHS allocation failed".into()))?,
        );
        let h = model.0.as_ptr();
        check(Highs_setBoolOptionValue(h, c"output_flag".as_ptr(), 0))?;
        check(Highs_setBoolOptionValue(h, c"log_to_console".as_ptr(), 0))?;
        check(Highs_setStringOptionValue(
            h,
            c"presolve".as_ptr(),
            c"on".as_ptr(),
        ))?;
        check(Highs_setStringOptionValue(
            h,
            c"parallel".as_ptr(),
            c"off".as_ptr(),
        ))?;
        check(Highs_setIntOptionValue(h, c"threads".as_ptr(), 1))?;
        check(Highs_setIntOptionValue(h, c"simplex_strategy".as_ptr(), 1))?;
        if let Some(limit) = limit {
            check(Highs_setIntOptionValue(
                h,
                c"simplex_iteration_limit".as_ptr(),
                limit as HighsInt,
            ))?;
        }
        check(Highs_setCallback(
            h,
            Some(interrupt),
            (control as *const SolveControl).cast_mut().cast(),
        ))?;
        check(Highs_startCallback(h, kHighsCallbackSimplexInterrupt))?;
        check(Highs_startCallback(h, kHighsCallbackIpmInterrupt))?;
        check(Highs_passLp(
            h,
            columns.len() as HighsInt,
            rhs.len() as HighsInt,
            values.len() as HighsInt,
            MATRIX_FORMAT_COLUMN_WISE,
            OBJECTIVE_SENSE_MINIMIZE,
            0.0,
            costs.as_ptr(),
            lower.as_ptr(),
            upper.as_ptr(),
            rhs.as_ptr(),
            rhs.as_ptr(),
            starts.as_ptr(),
            indexes.as_ptr(),
            values.as_ptr(),
        ))?;
        let completion = Highs_run(h);
        if control.is_cancelled() {
            return Err(Error::Incomplete(Termination::Cancelled));
        }
        check(completion)?;
        match Highs_getModelStatus(h) {
            MODEL_STATUS_OPTIMAL => (),
            MODEL_STATUS_REACHED_ITERATION_LIMIT => {
                return Err(Error::Incomplete(Termination::IterationLimit));
            }
            MODEL_STATUS_REACHED_INTERRUPT => {
                return Err(Error::Incomplete(Termination::Cancelled));
            }
            _ => {
                return Err(Error::Numerical(
                    "L1 model did not reach an optimal solution".into(),
                ));
            }
        }
        check(Highs_getSolution(
            h,
            solution.as_mut_ptr(),
            ptr::null_mut(),
            ptr::null_mut(),
            ptr::null_mut(),
        ))?;
        let mut iterations = 0;
        check(Highs_getIntInfoValue(
            h,
            c"simplex_iteration_count".as_ptr(),
            &mut iterations,
        ))?;
        if !finite(solution.iter().copied()) {
            return Err(Error::Numerical("nonfinite L1 solution".into()));
        }
        Ok((solution, iterations.max(0) as usize))
    }
}
