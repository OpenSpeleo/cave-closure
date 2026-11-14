use cave_closure::horizontal::{HorizontalGroup, HorizontalObservation, HorizontalProblem};
use cave_closure::{Error, SolveControl, Termination};

fn problem() -> HorizontalProblem {
    HorizontalProblem {
        groups: vec![
            HorizontalGroup {
                initial: [100.0, 200.0],
                fixed: Some([100.0, 200.0]),
            },
            HorizontalGroup {
                initial: [110.0, 200.0],
                fixed: None,
            },
        ],
        observations: vec![observation(10.0, 1.0), observation(12.0, 1.0)],
        sigma_length_m: 0.1,
        sigma_azimuth_deg: 1.0,
    }
}
fn observation(length: f64, quality: f64) -> HorizontalObservation {
    HorizontalObservation {
        from: 0,
        to: 1,
        displacement: [length, 0.0],
        length_m: length,
        azimuth_rad: std::f64::consts::FRAC_PI_2,
        quality,
    }
}
fn solve(p: &HorizontalProblem) -> Vec<[f64; 2]> {
    let result = p.solve(&SolveControl::default()).unwrap();
    assert!(result.axes.iter().all(|a| a.residual_infinity < 1e-7));
    result.coordinates
}
#[test]
fn parallel_observations_have_analytical_weighted_mean() {
    let mut p = problem();
    assert_eq!(solve(&p)[0], [100.0, 200.0]);
    assert!((solve(&p)[1][0] - 111.0).abs() < 1e-10);
    p.observations[1].quality = 2.0;
    assert!((solve(&p)[1][0] - 110.4).abs() < 1e-10);
}
#[test]
fn quality_fallback_and_zero_sigmas_retain_source_behavior() {
    let mut p = problem();
    let expected = solve(&p);
    for quality in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        p.observations[1].quality = quality;
        assert_eq!(solve(&p), expected);
    }
    p.sigma_length_m = 0.0;
    p.sigma_azimuth_deg = 0.0;
    assert_eq!(solve(&p)[1], p.groups[1].initial);
}
#[test]
fn collapsed_and_fixed_observations_do_not_move_controls() {
    let mut p = problem();
    let expected = solve(&p);
    let mut internal = observation(999.0, 1.0);
    internal.to = 0;
    p.observations.push(internal);
    assert_eq!(solve(&p), expected);
    p.groups[1].fixed = Some([300.0, 400.0]);
    assert_eq!(solve(&p), vec![[100.0, 200.0], [300.0, 400.0]]);
}
#[test]
fn reversed_edges_and_disconnected_groups_preserve_geometry() {
    let mut p = problem();
    let expected = solve(&p);
    for shot in &mut p.observations {
        std::mem::swap(&mut shot.from, &mut shot.to);
        shot.displacement = shot.displacement.map(|v| -v);
    }
    assert_eq!(solve(&p), expected);
    p.groups.push(HorizontalGroup {
        initial: [42.0, -31.0],
        fixed: None,
    });
    assert_eq!(solve(&p)[2], [42.0, -31.0]);
}
#[test]
fn displacement_is_not_rederived_from_length_or_depth() {
    let mut p = problem();
    // The caller may clamp an impossible horizontal projection to zero while
    // retaining slope length for instrument weighting. No depth is interpreted here.
    for shot in &mut p.observations {
        shot.displacement = [0.0, 0.0];
    }
    let point = solve(&p)[1];
    assert!((point[0] - 100.0).abs() < 1e-10);
    assert!((point[1] - 200.0).abs() < 1e-10);
}
#[test]
fn empty_problem_and_cancelled_solve() {
    let mut p = problem();
    let control = SolveControl::default();
    control.cancel();
    assert!(matches!(
        p.solve(&control),
        Err(Error::Incomplete(Termination::Cancelled))
    ));
    p.groups.clear();
    p.observations.clear();
    assert!(solve(&p).is_empty());
}
#[test]
fn malformed_requests_return_errors_without_panicking() {
    let mut p = problem();
    p.observations[0].to = usize::MAX;
    assert!(p.solve(&SolveControl::default()).is_err());
    let mut p = problem();
    p.groups[0].initial[0] = f64::NAN;
    assert!(p.solve(&SolveControl::default()).is_err());
    let mut p = problem();
    p.observations[0].displacement[0] = f64::INFINITY;
    assert!(p.solve(&SolveControl::default()).is_err());
    let mut p = problem();
    p.sigma_length_m = f64::NAN;
    assert!(p.solve(&SolveControl::default()).is_err());
}
#[test]
fn concurrent_repeated_solves_are_deterministic_and_leave_input_intact() {
    let p = problem();
    let before = format!("{p:?}");
    let expected = solve(&p);
    std::thread::scope(|scope| {
        let threads: Vec<_> = (0..8).map(|_| scope.spawn(|| solve(&p))).collect();
        for thread in threads {
            assert_eq!(thread.join().unwrap(), expected);
        }
    });
    assert_eq!(before, format!("{p:?}"));
}
