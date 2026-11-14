use super::{AxisDiagnostics, System};
use crate::{Error, SolveControl, Termination};
use std::collections::{BTreeMap, HashMap, VecDeque};

pub(super) struct KernelSolution {
    pub coordinates: [Vec<f64>; 2],
    pub axes: [AxisDiagnostics; 2],
}

pub(super) fn solve_reduced_laplacian(
    system: &System,
    x0_x: &[f64],
    x0_y: &[f64],
    control: &SolveControl,
) -> Result<KernelSolution, Error> {
    let System {
        diag_x,
        diag_y,
        atb_x,
        atb_y,
        adj,
    } = system;
    let n = diag_x.len();
    let mut axes = [AxisDiagnostics::default(); 2];
    // A node is only eliminated when its current pivot is *safely positive*
    // relative to the node's original diagonal. A pivot that is non-positive or
    // has collapsed to roundoff (a tiny fraction of its starting scale) means
    // the node sits in the system's null space - a floating axis. Such a pivot
    // is mathematically 0 but, with varied edge weights, lands at ~±1e-12 from
    // cancellation, so an absolute 1e-30 floor would wrongly eliminate it and
    // divide by ~0, shifting the whole component off its warm-started (measured)
    // translation. These stay in the core for the warm-started PCG instead.
    const PIVOT_REL: f64 = 1e-9;

    // Mutable working copies of the system (diagonals, rhs, adjacency).
    let mut dx = diag_x.to_vec();
    let mut dy = diag_y.to_vec();
    let mut cx = atb_x.to_vec();
    let mut cy = atb_y.to_vec();
    let mut nbr: Vec<BTreeMap<usize, (f64, f64)>> = vec![BTreeMap::new(); n];
    for (i, row) in adj.iter().enumerate() {
        for &(j, ox, oy) in row {
            let e = nbr[i].entry(j).or_insert((0.0, 0.0));
            e.0 += ox;
            e.1 += oy;
        }
    }

    // Back-substitution record captured at the moment a node is eliminated.
    struct Elim {
        u: usize,
        nbrs: Vec<(usize, f64, f64)>,
        dxu: f64,
        dyu: f64,
        cxu: f64,
        cyu: f64,
    }
    let mut eliminated = vec![false; n];
    let mut order: Vec<Elim> = Vec::new();

    let mut queue: VecDeque<usize> = (0..n).filter(|&i| nbr[i].len() <= 2).collect();
    while let Some(u) = queue.pop_front() {
        super::check_cancelled(control)?;
        // Stale entry: already gone, or degree grew above 2 since enqueued.
        if eliminated[u] || nbr[u].len() > 2 {
            continue;
        }
        let (pdx, pdy) = (dx[u], dy[u]);
        // Scale-relative, sign-aware floor: skip when either axis pivot is
        // non-positive or tiny vs. its original diagonal (null-space / floating).
        if pdx <= PIVOT_REL * diag_x[u] || pdy <= PIVOT_REL * diag_y[u] {
            continue; // null-space / ill-conditioned pivot → leave for the core solve
        }
        let nbrs: Vec<(usize, f64, f64)> =
            nbr[u].iter().map(|(&v, &(ox, oy))| (v, ox, oy)).collect();
        let (pcx, pcy) = (cx[u], cy[u]);
        eliminated[u] = true;
        order.push(Elim {
            u,
            nbrs: nbrs.clone(),
            dxu: pdx,
            dyu: pdy,
            cxu: pcx,
            cyu: pcy,
        });

        // Detach u from its neighbours.
        for &(v, _, _) in &nbrs {
            nbr[v].remove(&u);
        }
        // Schur complement onto each neighbour: diagonal + rhs.
        for &(v, ox, oy) in &nbrs {
            cx[v] -= ox * pcx / pdx;
            cy[v] -= oy * pcy / pdy;
            dx[v] -= ox * ox / pdx;
            dy[v] -= oy * oy / pdy;
        }
        // Schur complement fill between neighbour pairs (series reduction).
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let (i, oxa, oya) = nbrs[a];
                let (j, oxb, oyb) = nbrs[b];
                let fx = -(oxa * oxb) / pdx;
                let fy = -(oya * oyb) / pdy;
                let e = nbr[i].entry(j).or_insert((0.0, 0.0));
                e.0 += fx;
                e.1 += fy;
                let e = nbr[j].entry(i).or_insert((0.0, 0.0));
                e.0 += fx;
                e.1 += fy;
            }
        }
        // Neighbours may now be peelable.
        for &(v, _, _) in &nbrs {
            if !eliminated[v] && nbr[v].len() <= 2 {
                queue.push_back(v);
            }
        }
    }

    let mut sol_x = vec![0.0; n];
    let mut sol_y = vec![0.0; n];

    // Solve the irreducible 2-core (degree ≥ 3 / ill-conditioned nodes) with
    // Jacobi-PCG. For a cave this is at most a few nodes per loop.
    let core: Vec<usize> = (0..n).filter(|&i| !eliminated[i]).collect();
    if !core.is_empty() {
        let k = core.len();
        let mut local: HashMap<usize, usize> = HashMap::with_capacity(k);
        for (li, &g) in core.iter().enumerate() {
            local.insert(g, li);
        }
        let cdx: Vec<f64> = core.iter().map(|&g| dx[g]).collect();
        let cdy: Vec<f64> = core.iter().map(|&g| dy[g]).collect();
        let crx: Vec<f64> = core.iter().map(|&g| cx[g]).collect();
        let cry: Vec<f64> = core.iter().map(|&g| cy[g]).collect();
        let mut cadj: Vec<Vec<(usize, f64, f64)>> = vec![Vec::new(); k];
        for (li, &g) in core.iter().enumerate() {
            for (&v, &(ox, oy)) in &nbr[g] {
                // Every surviving neighbour is itself a core node.
                if let Some(&lv) = local.get(&v) {
                    cadj[li].push((lv, ox, oy));
                }
            }
        }
        let cx0: Vec<f64> = core.iter().map(|&g| x0_x[g]).collect();
        let cy0: Vec<f64> = core.iter().map(|&g| x0_y[g]).collect();
        let (csx, report_x) = pcg_axis(k, &cdx, &crx, &cadj, 0, &cx0, control)?;
        axes[0] = report_x;
        let (csy, report_y) = pcg_axis(k, &cdy, &cry, &cadj, 1, &cy0, control)?;
        axes[1] = report_y;
        for (li, &g) in core.iter().enumerate() {
            sol_x[g] = csx[li];
            sol_y[g] = csy[li];
        }
    }

    // Back-substitute eliminated nodes in reverse elimination order.
    for rec in order.iter().rev() {
        super::check_cancelled(control)?;
        let mut sx = rec.cxu;
        let mut sy = rec.cyu;
        for &(v, ox, oy) in &rec.nbrs {
            sx -= ox * sol_x[v];
            sy -= oy * sol_y[v];
        }
        sol_x[rec.u] = sx / rec.dxu;
        sol_y[rec.u] = sy / rec.dyu;
    }

    Ok(KernelSolution {
        coordinates: [sol_x, sol_y],
        axes,
    })
}

/// Jacobi-preconditioned Conjugate Gradient for one axis of an SPD/PSD system
/// stored as diagonal + symmetric off-diagonal adjacency `(j, off_x, off_y)`
/// (`axis` selects the weight). Used only for the small reduced core.
///
/// Starts from the warm-start guess `x0` (the measured positions). This is not
/// just for convergence speed: a node that is unconstrained on this axis (e.g.
/// reachable only through a zero-length cardinal shot, so its weight is 0) has
/// an empty row, leaving that coordinate in the system's null space. CG never
/// moves a null-space direction, so it stays at `x0` - preserving the measured
/// position instead of snapping it to world 0.
fn pcg_axis(
    n: usize,
    diag: &[f64],
    rhs: &[f64],
    adj: &[Vec<(usize, f64, f64)>],
    axis: usize,
    x0: &[f64],
    control: &SolveControl,
) -> Result<(Vec<f64>, AxisDiagnostics), Error> {
    let off = |e: (usize, f64, f64)| if axis == 0 { e.1 } else { e.2 };
    let spmv = |x: &[f64]| -> Vec<f64> {
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] += diag[i] * x[i];
            for &edge in &adj[i] {
                y[i] += off(edge) * x[edge.0];
            }
        }
        y
    };
    let max_iter = n.min(10_000);
    let tol = 1e-12;
    let inv_diag: Vec<f64> = diag
        .iter()
        .map(|&d| if d.abs() > 1e-30 { 1.0 / d } else { 1.0 })
        .collect();

    super::check_cancelled(control)?;
    let mut report = AxisDiagnostics::default();
    let mut x = x0.to_vec();
    let ax0 = spmv(&x);
    let mut r: Vec<f64> = rhs.iter().zip(ax0.iter()).map(|(&b, &a)| b - a).collect();
    let mut z: Vec<f64> = r
        .iter()
        .zip(inv_diag.iter())
        .map(|(&ri, &mi)| ri * mi)
        .collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(z.iter()).map(|(&ri, &zi)| ri * zi).sum();

    // Converge relative to the *initial residual* ‖b - A·x0‖, not ‖b‖. With a
    // warm start that is the correct reference: return x0 only when it already
    // satisfies the system (residual ~0 - which includes pure null-space
    // directions), and otherwise let CG converge the range component while
    // leaving null-space directions at x0. Keying the early exit off ‖b‖ alone
    // would wrongly return a non-optimal x0 when Schur elimination cancels the
    // core RHS for a balanced loop even though the matrix rows are non-zero
    // (Codex #86).
    let r0_norm_sq: f64 = r.iter().map(|&v| v * v).sum();
    if r0_norm_sq < tol * tol {
        return Ok((x, report));
    }

    report.termination = Termination::IterationLimit;
    for iteration in 0..max_iter {
        super::check_cancelled(control)?;
        report.iterations = iteration + 1;
        let ap = spmv(&p);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(&pi, &api)| pi * api).sum();
        if pap.abs() < 1e-30 {
            report.termination = Termination::Stagnation;
            break;
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm_sq: f64 = r.iter().map(|&v| v * v).sum();
        if r_norm_sq < tol * tol * r0_norm_sq {
            report.termination = Termination::Converged;
            break;
        }
        for i in 0..n {
            z[i] = r[i] * inv_diag[i];
        }
        let rz_new: f64 = r.iter().zip(z.iter()).map(|(&ri, &zi)| ri * zi).sum();
        let beta = rz_new / rz;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }
    Ok((x, report))
}
