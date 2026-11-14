use super::*;
#[expect(
    clippy::too_many_arguments,
    reason = "preserve the independent pre-migration matrix cases"
)]
fn solve(
    n: usize,
    diag_x: &[f64],
    diag_y: &[f64],
    atb_x: &[f64],
    atb_y: &[f64],
    adj: &[Vec<(usize, f64, f64)>],
    x0_x: &[f64],
    x0_y: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(n, diag_x.len());
    let system = System {
        diag_x: diag_x.to_vec(),
        diag_y: diag_y.to_vec(),
        atb_x: atb_x.to_vec(),
        atb_y: atb_y.to_vec(),
        adj: adj.to_vec(),
    };
    let result =
        solver::solve_reduced_laplacian(&system, x0_x, x0_y, &SolveControl::default()).unwrap();
    let [x, y] = result.coordinates;
    (x, y)
}
/// Exactness guard for the degree-≤2 graph-reduction solver (#81): on a range of
/// random SPD weighted-Laplacian systems (chains, single cycles, and densely
/// connected cores), `solve_reduced_laplacian` must match a reference dense
/// Gaussian elimination on both axes.
#[test]
fn solve_reduced_laplacian_matches_dense() {
    // Deterministic xorshift RNG → reproducible [0,1) draws (no rand dep).
    fn next(state: &mut u64) -> f64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *state = x;
        ((x >> 11) as f64) / ((1u64 << 53) as f64)
    }

    // Dense SPD solve via Gaussian elimination with partial pivoting.
    fn dense_solve(mat: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
        let n = rhs.len();
        let mut a: Vec<Vec<f64>> = mat.to_vec();
        let mut b = rhs.to_vec();
        for col in 0..n {
            let mut piv = col;
            for r in (col + 1)..n {
                if a[r][col].abs() > a[piv][col].abs() {
                    piv = r;
                }
            }
            a.swap(col, piv);
            b.swap(col, piv);
            let d = a[col][col];
            let pivot_rhs = b[col];
            let (leading, trailing) = a.split_at_mut(col + 1);
            let pivot_row = &leading[col];
            for (row, rhs) in trailing.iter_mut().zip(&mut b[(col + 1)..]) {
                let f = row[col] / d;
                for (value, pivot_value) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                    *value -= f * pivot_value;
                }
                *rhs -= f * pivot_rhs;
            }
        }
        let mut x = vec![0.0; n];
        for col in (0..n).rev() {
            let mut s = b[col];
            for c in (col + 1)..n {
                s -= a[col][c] * x[c];
            }
            x[col] = s / a[col][col];
        }
        x
    }

    // (n, extra_random_edges) - 0 extra = pure chain (fully peels, empty core);
    // a few = small cycle core; many = a denser core that exercises the PCG path.
    let configs = [(2, 0), (6, 0), (8, 2), (20, 5), (40, 30), (60, 80)];
    let mut seed = 0x9E3779B97F4A7C15u64;

    for (n, extra) in configs {
        let mut adj: Vec<Vec<(usize, f64, f64)>> = vec![Vec::new(); n];
        let mut diag_x = vec![0.0; n];
        let mut diag_y = vec![0.0; n];

        let add_edge = |a: usize,
                        b: usize,
                        seed: &mut u64,
                        adj: &mut Vec<Vec<(usize, f64, f64)>>,
                        dx: &mut Vec<f64>,
                        dy: &mut Vec<f64>| {
            let wx = 0.5 + next(seed);
            let wy = 0.5 + next(seed);
            adj[a].push((b, -wx, -wy));
            adj[b].push((a, -wx, -wy));
            dx[a] += wx;
            dx[b] += wx;
            dy[a] += wy;
            dy[b] += wy;
        };

        // Spanning chain keeps the graph connected.
        for i in 0..n.saturating_sub(1) {
            add_edge(i, i + 1, &mut seed, &mut adj, &mut diag_x, &mut diag_y);
        }
        // Extra edges introduce loops (a non-trivial 2-core).
        for _ in 0..extra {
            let a = (next(&mut seed) * n as f64) as usize % n;
            let b = (next(&mut seed) * n as f64) as usize % n;
            if a != b {
                add_edge(a, b, &mut seed, &mut adj, &mut diag_x, &mut diag_y);
            }
        }
        // Grounding term on every diagonal → strictly SPD (mimics fixed-node edges).
        for i in 0..n {
            diag_x[i] += 0.25 + next(&mut seed);
            diag_y[i] += 0.25 + next(&mut seed);
        }
        let atb_x: Vec<f64> = (0..n).map(|_| next(&mut seed) * 4.0 - 2.0).collect();
        let atb_y: Vec<f64> = (0..n).map(|_| next(&mut seed) * 4.0 - 2.0).collect();

        // Dense reference matrices built from the SAME adj/diag (duplicate edges
        // summed identically to the solver).
        let mut ax = vec![vec![0.0; n]; n];
        let mut ay = vec![vec![0.0; n]; n];
        for i in 0..n {
            ax[i][i] = diag_x[i];
            ay[i][i] = diag_y[i];
            for &(j, ox, oy) in &adj[i] {
                ax[i][j] += ox;
                ay[i][j] += oy;
            }
        }
        let ref_x = dense_solve(&ax, &atb_x);
        let ref_y = dense_solve(&ay, &atb_y);

        // Strictly SPD (grounded) → unique solution, so a zero warm start matches
        // the dense reference exactly.
        let x0 = vec![0.0; n];
        let (sol_x, sol_y) = solve(n, &diag_x, &diag_y, &atb_x, &atb_y, &adj, &x0, &x0);

        for i in 0..n {
            assert!(
                (sol_x[i] - ref_x[i]).abs() < 1e-6,
                "n={n} extra={extra} x[{i}]: reduced={} dense={}",
                sol_x[i],
                ref_x[i]
            );
            assert!(
                (sol_y[i] - ref_y[i]).abs() < 1e-6,
                "n={n} extra={extra} y[{i}]: reduced={} dense={}",
                sol_y[i],
                ref_y[i]
            );
        }
    }
}

/// Regression for the warm-start fix (Codex #86): a node that is unconstrained on
/// one axis (zero pivot, zero-weight edge → null-space row) must keep its measured
/// position on that axis rather than snapping to world 0.
#[test]
fn solve_reduced_laplacian_preserves_unconstrained_axis_warm_start() {
    // Two nodes. Node 0 is grounded on both axes. Node 1 couples to node 0 on Y
    // only (the X off-diagonal weight is 0), and has zero X diagonal + zero X rhs
    // - i.e. nothing constrains its X. Its warm-start X (42) must survive.
    let n = 2;
    let diag_x = vec![1.0, 0.0];
    let diag_y = vec![2.0, 1.0];
    let atb_x = vec![5.0, 0.0];
    let atb_y = vec![0.0, 0.0];
    let adj = vec![vec![(1usize, 0.0, -1.0)], vec![(0usize, 0.0, -1.0)]];
    let x0_x = vec![0.0, 42.0];
    let x0_y = vec![0.0, 0.0];

    let (sol_x, _sol_y) = solve(n, &diag_x, &diag_y, &atb_x, &atb_y, &adj, &x0_x, &x0_y);

    assert!(
        (sol_x[1] - 42.0).abs() < 1e-9,
        "unconstrained-X node should retain its warm-start X (42), got {}",
        sol_x[1]
    );
    // The constrained node still solves normally.
    assert!(
        (sol_x[0] - 5.0).abs() < 1e-9,
        "node 0 X should be 5, got {}",
        sol_x[0]
    );
}

/// Regression for the early-exit fix (Codex #86 follow-up): when the reduced core
/// has a zero RHS but non-zero, full-rank matrix rows, the unique minimizer is 0 -
/// the warm start must NOT be returned just because `‖b‖ = 0`. Uses a K4 core (all
/// degree 3, so it survives degree-≤2 reduction and exercises the core PCG).
#[test]
fn solve_reduced_laplacian_zero_rhs_constrained_core_solves_not_warm_start() {
    let n = 4;
    // K4 + grounding: diagonal 4, every off-diagonal -1 → SPD, full rank.
    let diag = vec![4.0; 4];
    let mut adj: Vec<Vec<(usize, f64, f64)>> = vec![Vec::new(); n];
    for (i, row) in adj.iter_mut().enumerate() {
        for j in 0..n {
            if i != j {
                row.push((j, -1.0, -1.0));
            }
        }
    }
    let atb = vec![0.0; 4]; // zero RHS → unique solution is 0
    let x0_x = vec![10.0, -5.0, 3.0, 7.0];
    let x0_y = vec![-2.0, 8.0, -6.0, 1.0];

    let (sol_x, sol_y) = solve(n, &diag, &diag, &atb, &atb, &adj, &x0_x, &x0_y);

    for (i, v) in sol_x.iter().chain(sol_y.iter()).enumerate() {
        assert!(
            v.abs() < 1e-6,
            "zero-RHS constrained core must solve to 0, got {v} at {i}"
        );
    }
}

/// Codex #86 (3rd, null-space-through-peeled-nodes): a component that *floats* on
/// an axis (positive internal coupling but no grounding → 1-D translation null
/// space) must keep its measured positions, not jump, when its nodes are peeled by
/// degree-≤2 reduction before the core PCG sees them.
///
/// Mechanism: the surviving core node keeps its own measured value (warm start),
/// and back-substitution reproduces the measured *relative* offsets - so the whole
/// component lands exactly on its measured coordinates, matching the old
/// full-system CG. (In real surveys this case can't even arise: every
/// positive-length shot has positive weight on both axes, so floating components
/// require zero-length shots.)
#[test]
fn solve_reduced_laplacian_preserves_floating_component_nullspace() {
    let n = 3;
    // X: pure path Laplacian 0-1-2 (singular → floats). Y: grounded (nonsingular).
    let diag_x = vec![1.0, 2.0, 1.0];
    let diag_y = vec![2.0, 3.0, 2.0];
    let adj = vec![
        vec![(1usize, -1.0, -1.0)],
        vec![(0usize, -1.0, -1.0), (2usize, -1.0, -1.0)],
        vec![(1usize, -1.0, -1.0)],
    ];
    // X rhs consistent with measured x=(0,3,5): A_x·(0,3,5) = (-3,1,2).
    let atb_x = vec![-3.0, 1.0, 2.0];
    let atb_y = vec![0.0, 0.0, 0.0];
    let x0_x = vec![0.0, 3.0, 5.0];
    let x0_y = vec![0.0, 0.0, 0.0];

    let (sx, _sy) = solve(n, &diag_x, &diag_y, &atb_x, &atb_y, &adj, &x0_x, &x0_y);

    for (i, want) in [0.0_f64, 3.0, 5.0].iter().enumerate() {
        assert!(
            (sx[i] - want).abs() < 1e-6,
            "floating component node {i} jumped: got {}, want {want}",
            sx[i]
        );
    }
}

/// Regression for the scale-relative pivot floor (#89 / Codex). A *loop* that
/// floats on X with **varied** edge weights drives the last null-space pivot to
/// roundoff (~1e-12), not exactly 0 - so the old absolute 1e-30 floor eliminated
/// it and divided by ~0, shifting the component. The pivot must be judged relative
/// to its original diagonal so the node stays in the warm-started core and the
/// measured translation is preserved.
#[test]
fn solve_reduced_laplacian_floating_loop_varied_weights_no_shift() {
    // Triangle 0-1-2, pure (singular) Laplacian on X → floats; grounded on Y.
    // Varied off-diagonals so the final Schur pivot lands on roundoff, not 0.
    let n = 3;
    let diag_x = vec![3.3, 2.7, 4.0]; // = incident |off| sums (row-sum 0 ⇒ singular)
    let diag_y = vec![4.3, 3.7, 5.0]; // +1 grounding ⇒ nonsingular
    let adj = vec![
        vec![(1usize, -1.0, -1.0), (2usize, -2.3, -2.3)],
        vec![(0usize, -1.0, -1.0), (2usize, -1.7, -1.7)],
        vec![(0usize, -2.3, -2.3), (1usize, -1.7, -1.7)],
    ];
    // X rhs consistent with measured x=(0,2,5): atb_x = A_x·(0,2,5).
    let atb_x = vec![-13.5, -3.1, 16.6];
    let atb_y = vec![0.0, 0.0, 0.0];
    let x0_x = vec![0.0, 2.0, 5.0];
    let x0_y = vec![0.0, 0.0, 0.0];

    let (sx, _sy) = solve(n, &diag_x, &diag_y, &atb_x, &atb_y, &adj, &x0_x, &x0_y);

    for (i, want) in [0.0_f64, 2.0, 5.0].iter().enumerate() {
        assert!(
            (sx[i] - want).abs() < 1e-6,
            "floating loop node {i} shifted: got {}, want {want}",
            sx[i]
        );
    }
}
