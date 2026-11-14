//! Historical graph CG and Compass traverse-quality adjustment, entirely native.
use crate::*;
use nalgebra::DVector;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

type Shot = (ObservationId, usize, usize, Point);

#[derive(Clone)]
pub struct GraphConfig {
    pub reference: bool,
    pub horizontal: bool,
    pub max_iterations: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub weight_cap: f64,
}

#[derive(Clone)]
pub struct GraphProblem {
    pub ids: Vec<StationId>,
    pub shots: Vec<Shot>,
    pub positions: Vec<Point>,
    pub roots: Vec<usize>,
    pub anchors: BTreeSet<usize>,
    pub config: GraphConfig,
}

fn path(
    start: usize,
    stop: usize,
    adj: &[Vec<usize>],
    anchors: &BTreeSet<usize>,
) -> Option<Vec<usize>> {
    let mut parents = vec![usize::MAX; adj.len()];
    parents[start] = start;
    let mut queue = VecDeque::from([start]);
    while let Some(current) = queue.pop_front() {
        for &next in &adj[current] {
            if parents[next] != usize::MAX || (anchors.contains(&next) && next != stop) {
                continue;
            }
            parents[next] = current;
            if next == stop {
                let mut result = vec![stop];
                while *result.last().unwrap() != start {
                    result.push(parents[*result.last().unwrap()]);
                }
                result.reverse();
                return Some(result);
            }
            queue.push_back(next);
        }
    }
    None
}
fn key(u: usize, v: usize) -> (usize, usize) {
    (u.min(v), u.max(v))
}
fn length(p: Point) -> f64 {
    // Preserve the compensated summation and libm powers of the source model.
    let (mut high, mut low) = (0.0_f64, 0.0_f64);
    for x in p {
        let x = x.powf(std::hint::black_box(2.0));
        let next = high + x;
        low += if high.abs() >= x.abs() {
            (high - next) + x
        } else {
            (x - next) + high
        };
        high = next;
    }
    (high + low).powf(std::hint::black_box(0.5))
}

impl GraphProblem {
    /// Inspect the exact relative-weight model used by the graph solver.
    pub fn normalized_problem(&self) -> Result<NormalizedProblem, Error> {
        let quality = self.quality(&SolveControl::default())?;
        Ok(NormalizedProblem {
            schema_version: 1,
            algorithm: LoopClosureAlgorithm::LseUnreducedPcg,
            frame: Frame {
                kind: FrameKind::LocalEnu,
                length_unit: LengthUnit::Metres,
                provenance: vec!["Native graph relative-precision model".into()],
            },
            mode: if self.config.horizontal {
                LoopClosureMode::Horizontal {
                    preserved_z: self
                        .ids
                        .iter()
                        .cloned()
                        .zip(self.positions.iter().map(|p| p[2]))
                        .collect(),
                }
            } else {
                LoopClosureMode::Spatial
            },
            stations: self.ids.iter().cloned().map(|id| Station { id }).collect(),
            observations: self
                .shots
                .iter()
                .map(|(id, u, v, value)| Observation::Displacement {
                    id: id.clone(),
                    from: self.ids[self.roots[*u]].clone(),
                    to: self.ids[self.roots[*v]].clone(),
                    value: *value,
                    uncertainty: Uncertainty::Relative {
                        weight: self.weight(self.roots[*u], self.roots[*v], *value, &quality),
                    },
                    reading_ids: vec![],
                })
                .collect(),
            controls: self
                .anchors
                .iter()
                .map(|&i| Control::Hard {
                    id: format!("anchor:{i}"),
                    station: self.ids[i].clone(),
                    axes: if self.config.horizontal {
                        vec![Axis::X, Axis::Y]
                    } else {
                        vec![Axis::X, Axis::Y, Axis::Z]
                    },
                    value: self.positions[i][..if self.config.horizontal { 2 } else { 3 }].to_vec(),
                })
                .collect(),
            relationships: self
                .roots
                .iter()
                .enumerate()
                .filter(|(i, r)| *i != **r)
                .map(|(i, &r)| Relationship::Equivalent {
                    id: RelationshipId(format!("alias:{i}").into()),
                    from: self.ids[i].clone(),
                    to: self.ids[r].clone(),
                })
                .collect(),
        })
    }
    fn quality(&self, control: &SolveControl) -> Result<BTreeMap<(usize, usize), f64>, Error> {
        let mut result = BTreeMap::<_, f64>::new();
        if self.config.reference {
            return Ok(result);
        }
        let mut adjacency = vec![vec![]; self.ids.len()];
        let mut delta = BTreeMap::new();
        for (_, from, to, value) in &self.shots {
            let (u, v) = (self.roots[*from], self.roots[*to]);
            adjacency[u].push(v);
            adjacency[v].push(u);
            delta.insert((u, v), *value);
            delta.insert((v, u), value.map(|x| -x));
        }
        let mut anchors: Vec<_> = self.anchors.iter().copied().collect();
        anchors.sort_by_key(|&i| &self.ids[i]);
        for (index, &start) in anchors.iter().enumerate() {
            if control.is_cancelled() {
                return Err(Error::Incomplete(Termination::Cancelled));
            }
            for &stop in &anchors[index + 1..] {
                let Some(path) = path(start, stop, &adjacency, &self.anchors) else {
                    continue;
                };
                if path.len() < 2 {
                    continue;
                }
                let mut position = self.positions[start];
                for pair in path.windows(2) {
                    for (axis, value) in position.iter_mut().enumerate() {
                        *value += delta[&(pair[0], pair[1])][axis];
                    }
                }
                let misclosure = std::array::from_fn(|a| position[a] - self.positions[stop][a]);
                let quality = length(misclosure) / (path.len() - 1) as f64;
                for pair in path.windows(2) {
                    result
                        .entry(key(pair[0], pair[1]))
                        .and_modify(|v| *v = v.min(quality))
                        .or_insert(quality);
                }
            }
        }
        Ok(result)
    }
    fn weight(
        &self,
        u: usize,
        v: usize,
        displacement: Point,
        quality: &BTreeMap<(usize, usize), f64>,
    ) -> f64 {
        let distance = length(displacement);
        if self.config.reference {
            return if distance == 0.0 {
                self.config.weight_cap
            } else {
                self.config.weight_cap.min(1.0 / distance)
            };
        }
        let distance = distance.max(0.1);
        let base = 1.0 / (distance * distance);
        match quality.get(&key(u, v)) {
            Some(q) if *q > 1e-6 => base * (1.0 / (q * q)),
            _ => base,
        }
    }
    pub fn solve(
        &self,
        max_iterations: Option<usize>,
        control: &SolveControl,
    ) -> Result<SolveReport, Error> {
        if control.is_cancelled() {
            return Err(Error::Incomplete(Termination::Cancelled));
        }
        let quality = self.quality(control)?;
        let mut seen = BTreeSet::new();
        let mut free: Vec<_> = self
            .roots
            .iter()
            .copied()
            .filter(|i| seen.insert(*i) && !self.anchors.contains(i))
            .collect();
        if !self.config.reference {
            free.sort_by_key(|&i| &self.ids[i]);
        }
        let mut positions = self.positions.clone();
        let mut iterations = 0;
        let mut warnings = Vec::new();
        if (self.config.reference || self.anchors.len() >= 2) && !free.is_empty() {
            let mut origin = [0.0; 3];
            if !self.config.reference {
                let mut anchors: Vec<_> = self.anchors.iter().copied().collect();
                anchors.sort_by_key(|&i| &self.ids[i]);
                for (a, value) in origin.iter_mut().enumerate() {
                    *value = anchors.iter().map(|&i| positions[i][a]).sum::<f64>()
                        / anchors.len() as f64;
                }
            }
            let mut index = vec![None; self.ids.len()];
            for (i, &station) in free.iter().enumerate() {
                index[station] = Some(i);
            }
            let n = free.len();
            let mut matrix = vec![BTreeMap::<usize, f64>::new(); n];
            let mut rhs = [DVector::zeros(n), DVector::zeros(n), DVector::zeros(n)];
            for (_, from, to, delta) in &self.shots {
                let (u, v) = (self.roots[*from], self.roots[*to]);
                let weight = self.weight(u, v, *delta, &quality);
                match (index[u], index[v]) {
                    (Some(i), Some(j)) => {
                        *matrix[i].entry(i).or_default() += weight;
                        *matrix[j].entry(j).or_default() += weight;
                        *matrix[i].entry(j).or_default() -= weight;
                        *matrix[j].entry(i).or_default() -= weight;
                        for a in 0..3 {
                            rhs[a][i] -= weight * delta[a];
                            rhs[a][j] += weight * delta[a];
                        }
                    }
                    (Some(i), None) => {
                        *matrix[i].entry(i).or_default() += weight;
                        for a in 0..3 {
                            let fixed = positions[v][a] - origin[a];
                            if self.config.reference {
                                rhs[a][i] -= weight * delta[a];
                                rhs[a][i] += weight * fixed;
                            } else {
                                rhs[a][i] += weight * (fixed - delta[a]);
                            }
                        }
                    }
                    (None, Some(j)) => {
                        *matrix[j].entry(j).or_default() += weight;
                        for a in 0..3 {
                            let fixed = positions[u][a] - origin[a];
                            if self.config.reference {
                                rhs[a][j] += weight * delta[a];
                                rhs[a][j] += weight * fixed;
                            } else {
                                rhs[a][j] += weight * (fixed + delta[a]);
                            }
                        }
                    }
                    _ => (),
                }
            }
            if matrix.iter().any(|r| !r.is_empty()) {
                for a in 0..if self.config.horizontal { 2 } else { 3 } {
                    let initial = DVector::from_iterator(
                        n,
                        free.iter().map(|&i| positions[i][a] - origin[a]),
                    );
                    let (mut values, used, converged) = conjugate_gradient(
                        &matrix,
                        &rhs[a],
                        initial,
                        &self.config,
                        max_iterations,
                        control,
                    )?;
                    iterations = iterations.max(used);
                    if !converged {
                        if self.config.reference {
                            warnings.push(format!(
                                "Reference CG axis {a} returned its current iterate"
                            ));
                        } else {
                            values = direct(&matrix, &rhs[a], control)?;
                            warnings.push(format!("CG axis {a} used sparse direct fallback"));
                        }
                    }
                    if !finite(values.iter().copied()) {
                        return Err(Error::Numerical("nonfinite graph solution".into()));
                    }
                    for (i, &station) in free.iter().enumerate() {
                        positions[station][a] = values[i] + origin[a];
                    }
                }
            }
        }
        let coordinates: Vec<_> = self
            .ids
            .iter()
            .enumerate()
            .map(|(i, id)| Coordinate {
                station: id.clone(),
                value: std::array::from_fn(|a| {
                    if a == 2 && self.config.horizontal {
                        self.positions[i][a]
                    } else {
                        positions[self.roots[i]][a]
                    }
                }),
            })
            .collect();
        let mut objective = 0.0;
        let observations = self
            .shots
            .iter()
            .map(|(id, u, v, delta)| {
                let residual: Point = std::array::from_fn(|a| {
                    coordinates[*v].value[a] - coordinates[*u].value[a] - delta[a]
                });
                let dimensions = if self.config.horizontal { 2 } else { 3 };
                objective += self.weight(self.roots[*u], self.roots[*v], *delta, &quality)
                    * residual[..dimensions].iter().map(|x| x * x).sum::<f64>();
                ObservationResult {
                    id: id.clone(),
                    correction: residual,
                    residual,
                }
            })
            .collect();
        let report = SolveReport {
            schema_version: 1,
            status: Termination::Converged,
            coordinates,
            observations,
            frame: Frame {
                kind: FrameKind::LocalEnu,
                length_unit: LengthUnit::Metres,
                provenance: vec!["Native source-compatible graph loop closure".into()],
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
                original_unknowns: free.len(),
                kernel_unknowns: free.len(),
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
fn multiply(matrix: &[BTreeMap<usize, f64>], value: &DVector<f64>) -> DVector<f64> {
    DVector::from_iterator(
        matrix.len(),
        matrix.iter().map(|row| {
            row.iter()
                .fold(0.0, |sum, (&j, &w)| w.mul_add(value[j], sum))
        }),
    )
}
// Match the native BLAS reduction used by the reference implementation on
// Apple platforms. No Python or array-library runtime participates in solving.
#[cfg(target_os = "macos")]
#[link(name = "Accelerate", kind = "framework")]
unsafe extern "C" {
    fn cblas_ddot(n: i32, x: *const f64, inc_x: i32, y: *const f64, inc_y: i32) -> f64;
}
fn graph_dot(x: &DVector<f64>, y: &DVector<f64>) -> f64 {
    #[cfg(target_os = "macos")]
    {
        // Both vectors are owned contiguous allocations with equal lengths;
        // input resource limits keep the count within the BLAS integer range.
        unsafe { cblas_ddot(x.len() as i32, x.as_ptr(), 1, y.as_ptr(), 1) }
    }
    #[cfg(not(target_os = "macos"))]
    {
        x.dot(y)
    }
}

fn conjugate_gradient(
    matrix: &[BTreeMap<usize, f64>],
    rhs: &DVector<f64>,
    mut x: DVector<f64>,
    config: &GraphConfig,
    limit: Option<usize>,
    control: &SolveControl,
) -> Result<(DVector<f64>, usize, bool), Error> {
    let mut r = rhs - multiply(matrix, &x);
    let mut p = r.clone();
    let mut rho = graph_dot(&r, &r);
    let tolerance = if config.reference {
        config.absolute_tolerance
    } else {
        config
            .absolute_tolerance
            .max(config.relative_tolerance * graph_dot(rhs, rhs).sqrt())
    };
    let limit = limit.unwrap_or(config.max_iterations);
    for iteration in 0..limit {
        if control.is_cancelled() {
            return Err(Error::Incomplete(Termination::Cancelled));
        }
        if rho.sqrt() < tolerance {
            return Ok((x, iteration, true));
        }
        let ap = multiply(matrix, &p);
        let denominator = graph_dot(&p, &ap);
        if !denominator.is_finite()
            || (config.reference && denominator.abs() < 1e-15)
            || denominator == 0.0
        {
            return Ok((x, iteration, false));
        }
        let alpha = rho / denominator;
        x.axpy(alpha, &p, 1.0);
        r.axpy(-alpha, &ap, 1.0);
        let updated = graph_dot(&r, &r);
        p.scale_mut(updated / rho);
        p += &r;
        rho = updated;
    }
    Ok((x, limit, false))
}
#[cfg(feature = "sparse-direct")]
fn direct(
    matrix: &[BTreeMap<usize, f64>],
    rhs: &DVector<f64>,
    control: &SolveControl,
) -> Result<DVector<f64>, Error> {
    use faer::{
        Mat,
        prelude::*,
        sparse::{SparseColMat, Triplet},
    };
    if control.is_cancelled() {
        return Err(Error::Incomplete(Termination::Cancelled));
    }
    let entries: Vec<_> = matrix
        .iter()
        .enumerate()
        .flat_map(|(i, row)| row.iter().map(move |(&j, &v)| Triplet::new(i, j, v)))
        .collect();
    let failure = |e| Error::Numerical(format!("graph sparse factorization: {e:?}"));
    let matrix =
        SparseColMat::<usize, f64>::try_new_from_triplets(matrix.len(), matrix.len(), &entries)
            .map_err(failure)?;
    let value = Mat::from_fn(rhs.len(), 1, |i, _| rhs[i]);
    let solution = matrix
        .sp_lu()
        .map_err(|e| Error::Numerical(format!("graph LU: {e:?}")))?
        .solve(&value);
    if control.is_cancelled() {
        return Err(Error::Incomplete(Termination::Cancelled));
    }
    Ok(DVector::from_iterator(
        rhs.len(),
        (0..rhs.len()).map(|i| solution[(i, 0)]),
    ))
}
#[cfg(not(feature = "sparse-direct"))]
fn direct(
    _: &[BTreeMap<usize, f64>],
    _: &DVector<f64>,
    _: &SolveControl,
) -> Result<DVector<f64>, Error> {
    Err(crate::invalid(
        "graph direct fallback requires sparse-direct",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn network() -> GraphProblem {
        GraphProblem {
            ids: ["A", "B", "C", "D", "E"]
                .into_iter()
                .map(StationId::from)
                .collect(),
            shots: [("AB", 0, 1), ("BC", 1, 2), ("CD", 2, 3), ("BE", 1, 4)]
                .into_iter()
                .map(|(id, u, v)| (ObservationId::from(id), u, v, [10.0, 0.0, 0.0]))
                .collect(),
            positions: vec![
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
                [35.0, 0.0, 0.0],
                [20.0, 0.0, 0.0],
            ],
            roots: (0..5).collect(),
            anchors: BTreeSet::from([0, 3, 4]),
            config: GraphConfig {
                reference: false,
                horizontal: false,
                max_iterations: 60_000,
                absolute_tolerance: 1e-8,
                relative_tolerance: 1e-5,
                weight_cap: 10.0,
            },
        }
    }
    #[test]
    fn best_traverse_protects_shared_good_shots() {
        let p = network();
        let q = p.quality(&SolveControl::default()).unwrap();
        assert_eq!(q[&key(0, 1)], 0.0);
        assert_eq!(q[&key(1, 4)], 0.0);
        assert!((q[&key(1, 2)] - 5.0 / 3.0).abs() < 1e-12);
        assert!(p.weight(0, 1, [10.0, 0.0, 0.0], &q) > p.weight(1, 2, [10.0, 0.0, 0.0], &q));
    }
    #[test]
    fn traverses_cannot_cross_an_intermediate_anchor() {
        let mut p = network();
        p.shots.push(("ED".into(), 4, 3, [10.0, 0.0, 0.0]));
        let q = p.quality(&SolveControl::default()).unwrap();
        assert!((q[&key(1, 2)] - 5.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn shortest_path_is_selected_by_edge_count() {
        let adjacency = vec![vec![1, 2], vec![0, 3], vec![0, 4], vec![1], vec![2, 3]];
        assert_eq!(
            path(0, 3, &adjacency, &BTreeSet::from([0, 3])),
            Some(vec![0, 1, 3])
        );
    }
    #[test]
    fn blocked_anchor_forces_an_alternative_path() {
        let adjacency = vec![vec![1, 2], vec![0, 3], vec![0, 4], vec![1], vec![2, 3]];
        assert_eq!(
            path(0, 3, &adjacency, &BTreeSet::from([0, 1, 3])),
            Some(vec![0, 2, 4, 3])
        );
    }
    #[test]
    fn disconnected_anchor_pair_has_no_traverse() {
        assert!(path(0, 2, &[vec![1], vec![0], vec![]], &BTreeSet::from([0, 2])).is_none());
    }
    #[test]
    fn parallel_edge_quality_uses_the_last_measurement() {
        let mut p = network();
        p.shots.push(("AB-repeat".into(), 0, 1, [30.0, 0.0, 0.0]));
        let q = p.quality(&SolveControl::default()).unwrap();
        assert!((q[&key(0, 1)] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn zero_quality_uses_unit_quality_factor() {
        let p = network();
        let q = BTreeMap::from([(key(0, 1), 0.0)]);
        assert_eq!(p.weight(0, 1, [10.0, 0.0, 0.0], &q), 0.01);
    }
    #[test]
    fn quality_threshold_does_not_invent_infinite_precision() {
        let p = network();
        let q = BTreeMap::from([(key(0, 1), 1e-6)]);
        assert_eq!(p.weight(0, 1, [10.0, 0.0, 0.0], &q), 0.01);
    }
    #[test]
    fn unvisited_edges_keep_base_precision() {
        let p = network();
        assert_eq!(p.weight(0, 1, [10.0, 0.0, 0.0], &BTreeMap::new()), 0.01);
    }
    #[test]
    fn graph_length_floor_is_applied_before_squaring() {
        let p = network();
        let w = p.weight(0, 1, [0.001, 0.0, 0.0], &BTreeMap::new());
        assert!((w - 100.0).abs() < 1e-12);
    }
    #[test]
    fn reference_method_ignores_traverse_quality() {
        let mut p = network();
        p.config.reference = true;
        assert!(p.quality(&SolveControl::default()).unwrap().is_empty());
        assert_eq!(
            p.weight(0, 1, [10.0, 0.0, 0.0], &BTreeMap::from([(key(0, 1), 5.0)])),
            0.1
        );
    }
    #[test]
    fn native_fallback_preserves_zero_rhs() {
        let p = network();
        let r = p.solve(Some(0), &SolveControl::default());
        if cfg!(feature = "sparse-direct") {
            assert!(r.is_ok());
        } else {
            assert!(r.is_err());
        }
    }
}

#[cfg(test)]
mod rounding_contract_tests {
    use super::*;
    #[test]
    fn compensated_length_retains_small_squared_terms() {
        // Naive summation loses both unit squares and returns the original
        // length. Compensation retains them, giving the next float above it.
        let large = 100_000_000.0_f64;
        let expected = f64::from_bits(large.to_bits() + 1);
        assert_eq!(length([large, 1.0, 1.0]), expected);
    }
    #[test]
    fn length_retains_source_libm_rounding_at_power_of_two_boundaries() {
        // The source model uses libm pow, which returns 32 on macOS and the
        // adjacent float below it on Linux for this boundary case.
        let upper = 32.0_f64;
        let lower = f64::from_bits(upper.to_bits() - 1);
        let actual = length([31.587966364419366, 5.0199981036162455, 1.0]);
        assert!(
            (lower..=upper).contains(&actual),
            "unexpected length: {actual}"
        );
    }
    #[test]
    fn sparse_product_retains_fused_rounding() {
        let factor = 2.0_f64.powi(-27);
        let rows = vec![
            BTreeMap::from([(0, -1.0), (1, 1.0 + factor)]),
            BTreeMap::new(),
        ];
        let value = DVector::from_vec(vec![1.0, 1.0 - factor]);
        assert_eq!(multiply(&rows, &value)[0], -2.0_f64.powi(-54));
    }
}
