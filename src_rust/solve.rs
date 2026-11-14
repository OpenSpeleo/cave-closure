use crate::*;
use nalgebra::{DMatrix, DVector};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};

#[derive(Clone, Default)]
pub struct SolveControl {
    cancelled: Arc<AtomicBool>,
}
impl SolveControl {
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
}
/// A workspace belongs to one solve at a time. No zero-allocation guarantee.
#[derive(Default)]
pub struct Workspace {
    coordinates: Vec<DVector<f64>>,
}

#[derive(Clone)]
struct System {
    h: Vec<BTreeMap<usize, DMatrix<f64>>>,
    g: Vec<DVector<f64>>,
}
impl System {
    fn new(sizes: &[usize]) -> Self {
        Self {
            h: sizes.iter().map(|_| BTreeMap::new()).collect(),
            g: sizes.iter().map(|n| DVector::zeros(*n)).collect(),
        }
    }
    fn add(&mut self, u: usize, v: usize, value: DMatrix<f64>) {
        if value.nrows() == 0 || value.ncols() == 0 {
            return;
        }
        *self.h[u]
            .entry(v)
            .or_insert_with(|| DMatrix::zeros(value.nrows(), value.ncols())) += value;
    }
    fn apply(&self, x: &[DVector<f64>]) -> Vec<DVector<f64>> {
        self.h
            .iter()
            .zip(&self.g)
            .map(|(row, g)| {
                let mut y = DVector::zeros(g.len());
                for (j, block) in row {
                    y += block * &x[*j];
                }
                y
            })
            .collect()
    }
    fn residual(&self, x: &[DVector<f64>]) -> Vec<DVector<f64>> {
        self.apply(x)
            .iter()
            .zip(&self.g)
            .map(|(hx, g)| g - hx)
            .collect()
    }
    fn norm(&self) -> f64 {
        self.h
            .iter()
            .zip(&self.g)
            .map(|(row, g)| {
                (0..g.len())
                    .map(|i| {
                        row.values()
                            .map(|b| (0..b.ncols()).map(|j| b[(i, j)].abs()).sum::<f64>())
                            .sum::<f64>()
                    })
                    .fold(0., f64::max)
            })
            .fold(0., f64::max)
    }
    fn is_finite(&self) -> bool {
        self.g.iter().all(|v| finite(v.iter().copied()))
            && self
                .h
                .iter()
                .all(|r| r.values().all(|b| finite(b.iter().copied())))
    }
}
fn norm(v: &[DVector<f64>]) -> f64 {
    if v.iter().any(|v| !finite(v.iter().copied())) {
        return f64::INFINITY;
    }
    v.iter().map(|v| v.amax()).fold(0., f64::max)
}
fn dot(a: &[DVector<f64>], b: &[DVector<f64>]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a.dot(b)).sum()
}
fn backward(h: &System, x: &[DVector<f64>]) -> f64 {
    let r = norm(&h.residual(x));
    let denominator = h.norm() * norm(x) + norm(&h.g);
    if !denominator.is_finite() {
        return f64::INFINITY;
    }
    if denominator == 0. {
        if r == 0. { 0. } else { f64::MAX }
    } else {
        r / denominator
    }
}
struct Design {
    jac: Vec<(usize, DMatrix<f64>)>,
    constant: DVector<f64>,
    precision: DMatrix<f64>,
}
struct Elimination {
    u: usize,
    diagonal: nalgebra::linalg::Cholesky<f64, nalgebra::Dyn>,
    rhs: DVector<f64>,
    neighbors: Vec<(usize, DMatrix<f64>)>,
}
fn reduce(system: &mut System, control: &SolveControl) -> Vec<Elimination> {
    let mut pending: BTreeSet<_> = (0..system.g.len())
        .filter(|u| {
            !system.g[*u].is_empty() && system.h[*u].keys().filter(|v| **v != *u).count() <= 2
        })
        .collect();
    let mut records = Vec::new();
    while let Some(u) = pending.pop_first() {
        if control.is_cancelled() {
            break;
        }
        let neighbors: Vec<_> = system.h[u]
            .iter()
            .filter(|(v, _)| **v != u)
            .map(|(v, b)| (*v, b.clone()))
            .collect();
        if neighbors.len() > 2 || system.g[u].is_empty() {
            continue;
        }
        let Some(diag) = system.h[u].get(&u) else {
            continue;
        };
        let Some(chol) = diag.clone().cholesky() else {
            continue;
        };
        // Retain unsafe pivots; never regularize or alter the objective.
        let scale = system.h[u].values().map(DMatrix::amax).fold(0., f64::max);
        if chol
            .l()
            .diagonal()
            .iter()
            .any(|v| v * v <= 64. * f64::EPSILON * scale)
        {
            continue;
        }
        let rhs = system.g[u].clone();
        let solved_rhs = chol.solve(&rhs);
        let updates: Vec<_> = neighbors
            .iter()
            .map(|(i, hi)| {
                let gi = hi.transpose() * &solved_rhs;
                let blocks: Vec<_> = neighbors
                    .iter()
                    .map(|(j, hj)| (*j, -hi.transpose() * chol.solve(hj)))
                    .collect();
                (*i, gi, blocks)
            })
            .collect();
        if updates.iter().any(|(_, g, blocks)| {
            !finite(g.iter().copied()) || blocks.iter().any(|(_, b)| !finite(b.iter().copied()))
        }) {
            continue;
        }
        for (i, g, blocks) in updates {
            system.g[i] -= g;
            for (j, b) in blocks {
                system.add(i, j, b);
            }
            system.h[i].remove(&u);
            pending.insert(i);
        }
        system.h[u].clear();
        system.g[u] = DVector::zeros(0);
        records.push(Elimination {
            u,
            diagonal: chol,
            rhs,
            neighbors,
        });
    }
    records
}
fn pcg(
    h: &System,
    config: &SolveConfig,
    control: &SolveControl,
) -> (Vec<DVector<f64>>, Termination, usize) {
    let mut x: Vec<_> = h.g.iter().map(|g| DVector::zeros(g.len())).collect();
    let mut inverse = Vec::new();
    for (u, g) in h.g.iter().enumerate() {
        if g.is_empty() {
            inverse.push(None);
            continue;
        }
        let Some(diag) = h.h[u].get(&u) else {
            return (x, Termination::NumericalBreakdown, 0);
        };
        let Some(chol) = diag.clone().cholesky() else {
            return (x, Termination::NumericalBreakdown, 0);
        };
        inverse.push(Some(chol));
    }
    let precondition = |r: &[DVector<f64>]| {
        r.iter()
            .zip(&inverse)
            .map(|(r, c)| c.as_ref().map_or_else(|| DVector::zeros(0), |c| c.solve(r)))
            .collect::<Vec<_>>()
    };
    let mut r = h.g.clone();
    let matrix_norm = h.norm();
    let rhs_norm = norm(&h.g);
    // True residuals incur cancellation at the scale of H*x, even when g is
    // small. Use the same scale as the original-system qualification below.
    let tolerance = |x: &[DVector<f64>]| {
        config.absolute_tolerance + config.relative_tolerance * (matrix_norm * norm(x) + rhs_norm)
    };
    if norm(&r) <= tolerance(&x) && backward(h, &x) <= config.backward_error_target {
        return (x, Termination::Converged, 0);
    }
    let mut z = precondition(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut previous = f64::INFINITY;
    let mut stagnant = 0;
    for k in 0..config.max_iterations {
        if control.is_cancelled() {
            return (x, Termination::Cancelled, k);
        }
        let hp = h.apply(&p);
        let curvature = dot(&p, &hp);
        if !curvature.is_finite() || curvature <= 0. || !rz.is_finite() || rz <= 0. {
            return (x, Termination::NumericalBreakdown, k);
        }
        let alpha = rz / curvature;
        for u in 0..x.len() {
            x[u] += alpha * &p[u];
        }
        // Recompute the true residual, rather than trusting recurrence convergence.
        r = h.residual(&x);
        if !finite(x.iter().flat_map(|v| v.iter().copied()))
            || !finite(r.iter().flat_map(|v| v.iter().copied()))
        {
            return (x, Termination::NumericalBreakdown, k + 1);
        }
        let current = norm(&r);
        if current <= tolerance(&x) && backward(h, &x) <= config.backward_error_target {
            return (x, Termination::Converged, k + 1);
        }
        if current >= previous * (1. - 8. * f64::EPSILON) {
            stagnant += 1;
        } else {
            stagnant = 0;
        }
        if stagnant >= 64 {
            return (x, Termination::Stagnation, k + 1);
        }
        previous = current;
        z = precondition(&r);
        let next = dot(&r, &z);
        let beta = next / rz;
        for u in 0..p.len() {
            p[u] = &z[u] + beta * &p[u];
        }
        rz = next;
    }
    (x, Termination::IterationLimit, config.max_iterations)
}

impl PreparedProblem {
    pub fn solve(
        &self,
        config: &SolveConfig,
        workspace: &mut Workspace,
        control: &SolveControl,
    ) -> Result<SolveReport, Error> {
        if !finite([
            config.relative_tolerance,
            config.absolute_tolerance,
            config.backward_error_target,
        ]) || config.relative_tolerance < 0.
            || config.absolute_tolerance < 0.
            || config.backward_error_target <= 0.
        {
            return Err(invalid("invalid stopping tolerances"));
        }
        let d = self.source.mode.dimension();
        let nr = self.fixed.len();
        let axes: Vec<Vec<_>> = self
            .fixed
            .iter()
            .map(|f| (0..d).filter(|a| f[*a].is_none()).collect())
            .collect();
        let sizes: Vec<_> = axes.iter().map(Vec::len).collect();
        let unknowns = sizes.iter().sum::<usize>();
        // Center each connected component before assembly. Translation is exact
        // algebraically and avoids injecting geographic offsets into the RHS.
        let mut adjacency = vec![BTreeSet::new(); nr];
        for f in &self.factors {
            let (u, v) = (self.roots[f.u], self.roots[f.v]);
            adjacency[u].insert(v);
            adjacency[v].insert(u);
        }
        let mut suggestions = vec![vec![None; d]; nr];
        for soft in &self.soft {
            for (i, a) in soft.axes.iter().enumerate() {
                suggestions[soft.root][*a].get_or_insert(soft.value[i]);
            }
        }
        let mut centers = vec![DVector::zeros(d); nr];
        let mut components = Vec::new();
        let mut visited = vec![false; nr];
        for start in 0..nr {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            let mut members = vec![start];
            let mut cursor = 0;
            while cursor < members.len() {
                let u = members[cursor];
                cursor += 1;
                for v in &adjacency[u] {
                    if !visited[*v] {
                        visited[*v] = true;
                        members.push(*v);
                    }
                }
            }
            let mut center = DVector::zeros(d);
            for a in 0..d {
                center[a] = members
                    .iter()
                    .find_map(|u| self.fixed[*u][a])
                    .or_else(|| members.iter().find_map(|u| suggestions[*u][a]))
                    .unwrap_or(0.);
            }
            for u in &members {
                centers[*u] = center.clone();
            }
            components.push(members);
        }
        let fixed: Vec<_> = self
            .fixed
            .iter()
            .enumerate()
            .map(|(u, f)| {
                DVector::from_iterator(
                    d,
                    f.iter()
                        .enumerate()
                        .map(|(a, v)| v.map_or(0., |v| v - centers[u][a])),
                )
            })
            .collect();
        let selector = |root: usize| {
            DMatrix::from_fn(
                d,
                sizes[root],
                |a, k| if axes[root][k] == a { 1. } else { 0. },
            )
        };
        let mut designs = Vec::new();
        for f in &self.factors {
            let (u, v) = (self.roots[f.u], self.roots[f.v]);
            let constant =
                &fixed[v] - &fixed[u] + &self.offsets[f.v] - &self.offsets[f.u] - &f.measured;
            let jac = if u == v {
                Vec::new()
            } else {
                vec![(u, -selector(u)), (v, selector(v))]
            };
            designs.push(Design {
                jac,
                constant,
                precision: f.precision.clone(),
            });
        }
        for s in &self.soft {
            let k = s.axes.len();
            let jac = DMatrix::from_fn(k, sizes[s.root], |i, j| {
                if s.axes[i] == axes[s.root][j] { 1. } else { 0. }
            });
            let constant = DVector::from_iterator(
                k,
                s.axes
                    .iter()
                    .enumerate()
                    .map(|(i, a)| fixed[s.root][*a] + centers[s.root][*a] - s.value[i]),
            );
            designs.push(Design {
                jac: vec![(s.root, jac)],
                constant,
                precision: s.precision.clone(),
            });
        }
        let mut original = System::new(&sizes);
        for design in &designs {
            if control.is_cancelled() {
                return Err(Error::Incomplete(Termination::Cancelled));
            }
            for (u, ju) in &design.jac {
                original.g[*u] -= ju.transpose() * &design.precision * &design.constant;
                for (v, jv) in &design.jac {
                    original.add(*u, *v, ju.transpose() * &design.precision * jv);
                }
            }
        }
        if !original.is_finite() || !original.norm().is_finite() {
            return Err(Error::Numerical("assembly overflow".into()));
        }
        // Qualify grounding in the assembled floating-point system. A weak
        // unary control can disappear when added to a strong network diagonal.
        let validate_normal_grounding = match self.source.algorithm.backend()? {
            #[cfg(feature = "sparse-direct")]
            LinearBackend::SparseQr => false,
            _ => true,
        };
        for members in &components {
            let controlled_axes =
                (0..d).map(|axis| members.iter().any(|u| self.fixed[*u][axis].is_some()));
            for (a, is_controlled) in controlled_axes.enumerate() {
                if !validate_normal_grounding || is_controlled {
                    continue;
                }
                let mut energy = 0.;
                for u in members {
                    let i = axes[*u]
                        .iter()
                        .position(|axis| *axis == a)
                        .expect("free translation axis");
                    for (v, block) in &original.h[*u] {
                        if let Some(j) = axes[*v].iter().position(|axis| *axis == a) {
                            energy += block[(i, j)];
                        }
                    }
                }
                if !energy.is_finite() || energy <= 0. {
                    return Err(Error::Numerical(
                        "grounding was lost during floating-point assembly".into(),
                    ));
                }
            }
        }
        let (mut q, mut status, mut iterations, kernel_unknowns) =
            match self.source.algorithm.backend()? {
                LinearBackend::Pcg => {
                    let mut h = original.clone();
                    let records = if self.source.algorithm == LoopClosureAlgorithm::LseReducedPcg {
                        reduce(&mut h, control)
                    } else {
                        Vec::new()
                    };
                    let core = h.g.iter().map(DVector::len).sum();
                    // Reconstruction can amplify kernel residuals. Solve the kernel
                    // more tightly, then qualify the original-system residual below.
                    let inner_config = SolveConfig {
                        relative_tolerance: config.relative_tolerance * 0.01,
                        absolute_tolerance: config.absolute_tolerance * 0.01,
                        ..config.clone()
                    };
                    let (mut x, status, k) = pcg(&h, &inner_config, control);
                    for record in records.iter().rev() {
                        let mut rhs = record.rhs.clone();
                        for (v, b) in &record.neighbors {
                            rhs -= b * &x[*v];
                        }
                        x[record.u] = record.diagonal.solve(&rhs);
                    }
                    (x, status, k, core)
                }
                #[cfg(feature = "sparse-direct")]
                backend => {
                    let x = direct(&original, &designs, &sizes, backend, control)?;
                    (x, Termination::Converged, 0, unknowns)
                }
            };
        let mut refinement_steps = 0;
        if matches!(self.source.algorithm.backend()?, LinearBackend::Pcg)
            && (self.source.algorithm == LoopClosureAlgorithm::LseReducedPcg)
        {
            for _ in 0..config.max_refinement_steps {
                if status != Termination::Converged || control.is_cancelled() {
                    break;
                }
                let residual = original.residual(&q);
                let before = norm(&residual);
                if before
                    <= config.absolute_tolerance
                        + config.relative_tolerance
                            * (original.norm() * norm(&q) + norm(&original.g))
                {
                    break;
                }
                let mut correction = original.clone();
                correction.g = residual;
                let records = reduce(&mut correction, control);
                let inner = SolveConfig {
                    max_iterations: config.max_iterations.saturating_sub(iterations),
                    absolute_tolerance: config.absolute_tolerance * 0.01,
                    relative_tolerance: config.relative_tolerance * 0.01,
                    ..config.clone()
                };
                let (mut delta, completion, steps) = pcg(&correction, &inner, control);
                for record in records.iter().rev() {
                    let mut rhs = record.rhs.clone();
                    for (v, block) in &record.neighbors {
                        rhs -= block * &delta[*v];
                    }
                    delta[record.u] = record.diagonal.solve(&rhs);
                }
                for (value, change) in q.iter_mut().zip(delta) {
                    *value += change;
                }
                iterations += steps;
                refinement_steps += 1;
                status = completion;
                if norm(&original.residual(&q)) >= before {
                    break;
                }
            }
        }
        if control.is_cancelled() {
            status = Termination::Cancelled;
        }
        if !finite(q.iter().flat_map(|v| v.iter().copied())) {
            return Err(Error::Numerical("non-finite reconstruction".into()));
        }
        let error = backward(&original, &q);
        let true_residual = norm(&original.residual(&q));
        if status == Termination::Converged
            && (!error.is_finite()
                || error > config.backward_error_target
                || true_residual
                    > config.absolute_tolerance
                        + config.relative_tolerance
                            * (original.norm() * norm(&q) + norm(&original.g)))
        {
            status = Termination::NumericalBreakdown;
        }
        workspace.coordinates.clear();
        for root in 0..nr {
            let mut value = &fixed[root] + &centers[root];
            for (k, a) in axes[root].iter().enumerate() {
                value[*a] = q[root][k] + centers[root][*a];
            }
            workspace.coordinates.push(value);
        }
        let coordinates: Vec<_> = self
            .ids
            .iter()
            .enumerate()
            .map(|(i, id)| {
                let v = &workspace.coordinates[self.roots[i]] + &self.offsets[i];
                Coordinate {
                    station: id.clone(),
                    value: [
                        v[0],
                        v[1],
                        if d == 3 {
                            v[2]
                        } else {
                            let LoopClosureMode::Horizontal { preserved_z } = &self.source.mode
                            else {
                                unreachable!()
                            };
                            preserved_z[id]
                        },
                    ],
                }
            })
            .collect();
        let mut objective = 0.;
        let mut results = Vec::new();
        for f in &self.factors {
            let u = &coordinates[f.u].value;
            let v = &coordinates[f.v].value;
            let r = DVector::from_iterator(d, (0..d).map(|a| v[a] - u[a] - f.measured[a]));
            objective += r.dot(&(&f.precision * &r));
            let residual = [r[0], r[1], if d == 3 { r[2] } else { 0. }];
            results.push(ObservationResult {
                id: f.id.clone(),
                correction: residual,
                residual,
            });
        }
        for s in &self.soft {
            let r = DVector::from_iterator(
                s.axes.len(),
                s.axes
                    .iter()
                    .enumerate()
                    .map(|(i, a)| workspace.coordinates[s.root][*a] - s.value[i]),
            );
            objective += r.dot(&(&s.precision * &r));
        }
        let mut hard_error: f64 = 0.;
        let index: BTreeMap<_, _> = self.ids.iter().enumerate().map(|(i, id)| (id, i)).collect();
        for rel in &self.source.relationships {
            let (_, u, v, offset) = rel.parts();
            if let Some(offset) = offset {
                for (a, expected) in offset.iter().enumerate().take(d) {
                    hard_error = hard_error.max(
                        (coordinates[index[v]].value[a]
                            - coordinates[index[u]].value[a]
                            - expected)
                            .abs(),
                    );
                }
            }
        }
        for c in &self.source.controls {
            if let Control::Hard {
                station,
                axes,
                value,
                ..
            } = c
            {
                for (axis, expected) in axes.iter().zip(value) {
                    hard_error = hard_error
                        .max((coordinates[index[station]].value[axis.index()] - expected).abs());
                }
            }
        }
        let redundancy = self
            .rows
            .checked_sub(unknowns)
            .ok_or_else(|| Error::Numerical("rank accounting failed".into()))?;
        let report = SolveReport {
            schema_version: 1,
            status,
            coordinates,
            observations: results,
            frame: self.source.frame.clone(),
            diagnostics: Diagnostics {
                objective,
                cycle_rank: self.cycle_rank,
                stochastic_redundancy: redundancy,
                artificial_gauges: self.gauges.clone(),
                warnings: self.warnings.clone(),
                calibrated_statistics_eligible: self.calibrated && status == Termination::Converged,
                variance_factor: if self.calibrated
                    && redundancy > 0
                    && status == Termination::Converged
                {
                    Some(objective / redundancy as f64)
                } else {
                    None
                },
                components: self.components,
                original_unknowns: unknowns,
                kernel_unknowns,
                iterations,
                refinement_steps,
                true_residual_infinity: true_residual,
                backward_error: error,
                hard_constraint_max_error: hard_error,
            },
        };
        report.to_json()?; // Prevent non-finite success and non-finite serialization.
        q.clear();
        Ok(report)
    }
}

#[cfg(feature = "sparse-direct")]
fn direct(
    h: &System,
    designs: &[Design],
    sizes: &[usize],
    backend: LinearBackend,
    control: &SolveControl,
) -> Result<Vec<DVector<f64>>, Error> {
    use faer::{
        Mat, Side,
        prelude::*,
        sparse::{SparseColMat, Triplet},
    };
    let n: usize = sizes.iter().sum();
    let mut offsets = vec![0; sizes.len()];
    for u in 1..sizes.len() {
        offsets[u] = offsets[u - 1] + sizes[u - 1];
    }
    if control.is_cancelled() {
        return Err(Error::Incomplete(Termination::Cancelled));
    }
    if n == 0 {
        return Ok(sizes.iter().map(|s| DVector::zeros(*s)).collect());
    }
    let failure = |e| Error::Numerical(format!("sparse factorization: {e:?}"));
    let solution = match backend {
        LinearBackend::SparseCholesky => {
            let mut triplets = Vec::new();
            for (u, row) in h.h.iter().enumerate() {
                for (v, b) in row {
                    for i in 0..b.nrows() {
                        for j in 0..b.ncols() {
                            triplets.push(Triplet::new(offsets[u] + i, offsets[*v] + j, b[(i, j)]));
                        }
                    }
                }
            }
            let matrix = SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
                .map_err(failure)?;
            let rhs = Mat::from_fn(n, 1, |i, _| {
                let u = offsets.partition_point(|o| *o <= i) - 1;
                h.g[u][i - offsets[u]]
            });
            matrix
                .sp_cholesky(Side::Lower)
                .map_err(|e| Error::Numerical(format!("Cholesky failed: {e:?}")))?
                .solve(&rhs)
        }
        LinearBackend::SparseQr => {
            let m: usize = designs.iter().map(|d| d.constant.len()).sum();
            let mut triplets = Vec::new();
            let mut rhs = Mat::zeros(m, 1);
            let mut row = 0;
            for d in designs {
                let l = d
                    .precision
                    .clone()
                    .cholesky()
                    .ok_or_else(|| Error::Numerical("whitening failed".into()))?
                    .l()
                    .transpose();
                let target = -&l * &d.constant;
                for i in 0..target.len() {
                    rhs[(row + i, 0)] = target[i];
                }
                for (u, jac) in &d.jac {
                    let a = &l * jac;
                    for i in 0..a.nrows() {
                        for j in 0..a.ncols() {
                            triplets.push(Triplet::new(row + i, offsets[*u] + j, a[(i, j)]));
                        }
                    }
                }
                row += target.len();
            }
            let matrix = SparseColMat::<usize, f64>::try_new_from_triplets(m, n, &triplets)
                .map_err(failure)?;
            matrix
                .sp_qr()
                .map_err(|e| Error::Numerical(format!("QR failed: {e:?}")))?
                .solve_lstsq(&rhs)
        }
        _ => unreachable!(),
    };
    if control.is_cancelled() {
        return Err(Error::Incomplete(Termination::Cancelled));
    }
    // Exact cancellation can make the normal right-hand side zero while QR
    // returns roundoff-sized coordinates. The exact solution is then zero;
    // retain it so a relative backward-error check does not report error 1.
    // Factorization above still qualifies the selected backend and grounding.
    if norm(&h.g) == 0. {
        return Ok(sizes.iter().map(|s| DVector::zeros(*s)).collect());
    }
    Ok(sizes
        .iter()
        .enumerate()
        .map(|(u, s)| DVector::from_iterator(*s, (0..*s).map(|i| solution[(offsets[u] + i, 0)])))
        .collect())
}
