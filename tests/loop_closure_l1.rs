//! Analytical, ownership and validation regressions for L1 loop closure.
use cave_closure::l1::*;
use cave_closure::*;

fn request() -> L1Input {
    L1Input {
        station_ids: vec!["A".into(), "B".into(), "C".into()],
        observations: vec![
            ("AB".into(), 0, 1, Some([5.0, 0.0, 0.0])),
            ("BC".into(), 1, 2, Some([15.0, 0.0, 0.0])),
        ],
        initial_positions: vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [24.0, 0.0, 0.0]],
        anchor_indices: vec![0, 2],
        options: L1Options {
            mode: "spatial".into(),
            length_floor_m: 0.1,
            minimum_anchor_count: 2,
            deduplicate_edges: true,
            max_iterations: None,
        },
    }
}
#[cfg(feature = "l1")]
fn solve(input: L1Input) -> SolveReport {
    prepare_l1_loop_closure(input)
        .unwrap()
        .solve(None, &SolveControl::default())
        .unwrap()
}
#[cfg(feature = "l1")]
fn point(report: &SolveReport, id: &str) -> Point {
    report
        .coordinates
        .iter()
        .find(|p| p.station.0.as_ref() == id)
        .unwrap()
        .value
}
#[cfg(feature = "l1")]
#[test]
fn l1_concentrates_the_error_on_the_longer_leg() {
    let r = solve(request());
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
    assert!((r.diagnostics.objective - 4.0 / 15.0).abs() < 1e-12);
    assert!(!r.diagnostics.calibrated_statistics_eligible);
}
#[cfg(feature = "l1")]
#[test]
fn all_control_coordinates_remain_exact() {
    let r = solve(request());
    assert_eq!(point(&r, "A"), [0.0, 0.0, 0.0]);
    assert_eq!(point(&r, "C"), [24.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn vertical_l1_adjustment_uses_the_same_relative_weights() {
    let mut p = request();
    p.observations[0].3 = Some([0.0, 0.0, 5.0]);
    p.observations[1].3 = Some([0.0, 0.0, 15.0]);
    p.initial_positions = vec![[0.0; 3], [0.0, 0.0, 5.0], [0.0, 0.0, 24.0]];
    assert_eq!(point(&solve(p), "B"), [0.0, 0.0, 5.0]);
}
#[cfg(feature = "l1")]
#[test]
fn horizontal_l1_preserves_recorded_depths() {
    let mut p = request();
    p.options.mode = "horizontal".into();
    p.initial_positions[1][2] = -3.0;
    p.initial_positions[2][2] = -8.0;
    let r = solve(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, -3.0]);
    assert_eq!(point(&r, "C")[2], -8.0);
}
#[cfg(feature = "l1")]
#[test]
fn reversing_every_observation_keeps_the_same_geometry() {
    let expected = solve(request());
    let mut p = request();
    for (_, u, v, delta) in &mut p.observations {
        std::mem::swap(u, v);
        *delta = delta.map(|p| p.map(|x| -x));
    }
    let r = solve(p);
    assert_eq!(point(&r, "B"), point(&expected, "B"));
}
#[cfg(feature = "l1")]
#[test]
fn large_global_offsets_do_not_change_the_solution() {
    let mut p = request();
    for point in &mut p.initial_positions {
        for a in point {
            *a += 1e6;
        }
    }
    assert_eq!(
        point(&solve(p), "B"),
        [1_000_005.0, 1_000_000.0, 1_000_000.0]
    );
}
#[cfg(feature = "l1")]
#[test]
fn zero_length_leg_uses_a_finite_weight_floor() {
    let mut p = request();
    p.observations[0].3 = Some([0.0; 3]);
    let r = solve(p);
    assert_eq!(point(&r, "B"), [0.0; 3]);
    assert!(r.diagnostics.objective.is_finite());
}
#[cfg(feature = "l1")]
#[test]
fn missing_measurement_uses_the_owned_initial_displacement() {
    let mut p = request();
    p.observations[1].3 = None;
    let r = solve(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
    assert!(r.diagnostics.objective < 1e-12);
}
#[cfg(feature = "l1")]
#[test]
fn missing_anchors_preserves_initial_positions() {
    let mut p = request();
    p.anchor_indices.clear();
    let r = solve(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
    assert_eq!(point(&r, "C"), [24.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn one_anchor_preserves_the_initial_guess_by_source_policy() {
    let mut p = request();
    p.anchor_indices.truncate(1);
    let r = solve(p);
    assert_eq!(point(&r, "C"), [24.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn configurable_anchor_threshold_allows_a_single_anchor_solve() {
    let mut p = request();
    p.anchor_indices.truncate(1);
    p.options.minimum_anchor_count = 1;
    let r = solve(p);
    assert_eq!(point(&r, "C"), [20.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn no_observations_retains_all_initial_positions() {
    let mut p = request();
    p.observations.clear();
    let r = solve(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn all_anchored_network_does_not_move() {
    let mut p = request();
    p.anchor_indices = vec![0, 1, 2];
    let r = solve(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn reverse_duplicates_are_suppressed_only_when_requested() {
    let mut p = request();
    p.observations
        .push(("reverse".into(), 1, 0, Some([-50.0, 0.0, 0.0])));
    let r = solve(p.clone());
    assert_eq!(r.observations.len(), 2);
    assert_eq!(point(&r, "B"), [5.0, 0.0, 0.0]);
    p.options.deduplicate_edges = false;
    assert_eq!(solve(p).observations.len(), 3);
}
#[cfg(feature = "l1")]
#[test]
fn self_observation_preserves_source_last_write_semantics() {
    let mut p = request();
    p.observations.push(("self".into(), 1, 1, Some([0.0; 3])));
    let r = solve(p);
    assert_eq!(point(&r, "B"), [12.0, 0.0, 0.0]);
}
#[cfg(feature = "l1")]
#[test]
fn preparation_is_owned_after_input_mutation() {
    let mut p = request();
    let prepared = prepare_l1_loop_closure(p.clone()).unwrap();
    p.initial_positions[2][0] = 1e9;
    let r = prepared.solve(None, &SolveControl::default()).unwrap();
    assert_eq!(point(&r, "C")[0], 24.0);
}
#[cfg(feature = "l1")]
#[test]
fn reusable_preparations_support_concurrent_solves() {
    let p = std::sync::Arc::new(prepare_l1_loop_closure(request()).unwrap());
    let workers: Vec<_> = (0..4)
        .map(|_| {
            let p = p.clone();
            std::thread::spawn(move || {
                p.solve(None, &SolveControl::default())
                    .unwrap()
                    .to_json()
                    .unwrap()
            })
        })
        .collect();
    let outputs: Vec<_> = workers.into_iter().map(|w| w.join().unwrap()).collect();
    assert!(outputs.windows(2).all(|v| v[0] == v[1]));
}
#[test]
fn cancellation_is_checked_before_starting_a_model() {
    let p = prepare_l1_loop_closure(request()).unwrap();
    let token = SolveControl::default();
    token.cancel();
    assert!(matches!(
        p.solve(None, &token),
        Err(Error::Incomplete(Termination::Cancelled))
    ));
}
#[test]
fn preparation_exposes_immutable_source_evidence() {
    let p = request();
    let owned = prepare_l1_loop_closure(p).unwrap();
    let roundtrip: L1Input = serde_json::from_str(&owned.input_json().unwrap()).unwrap();
    assert_eq!(roundtrip.station_ids, vec!["A", "B", "C"]);
}
#[cfg(not(feature = "l1"))]
#[test]
fn solving_requires_the_explicit_l1_feature() {
    let p = prepare_l1_loop_closure(request()).unwrap();
    assert!(matches!(
        p.solve(None, &SolveControl::default()),
        Err(Error::Validation(_))
    ));
}
macro_rules! invalid {
    ($name:ident,$change:expr) => {
        #[test]
        fn $name() {
            let mut p = request();
            ($change)(&mut p);
            assert!(prepare_l1_loop_closure(p).is_err());
        }
    };
}
invalid!(duplicate_station, |p: &mut L1Input| p.station_ids[1] =
    "A".into());
invalid!(empty_station, |p: &mut L1Input| p.station_ids[1].clear());
invalid!(missing_position, |p: &mut L1Input| {
    p.initial_positions.pop();
});
invalid!(nonfinite_position, |p: &mut L1Input| p.initial_positions
    [1][0] = f64::NAN);
invalid!(invalid_anchor, |p: &mut L1Input| p.anchor_indices[0] = 100);
invalid!(duplicate_anchor, |p: &mut L1Input| p.anchor_indices[1] = 0);
invalid!(invalid_endpoint, |p: &mut L1Input| p.observations[0].1 =
    100);
invalid!(duplicate_observation_id, |p: &mut L1Input| p.observations
    [1]
.0 = "AB".into());
invalid!(nonfinite_displacement, |p: &mut L1Input| p.observations
    [0]
.3 =
    Some([f64::NAN; 3]));
invalid!(overflowing_displacement_length, |p: &mut L1Input| p
    .observations[0]
    .3 =
    Some([1e308; 3]));
invalid!(unknown_mode, |p: &mut L1Input| p.options.mode =
    "typo".into());
invalid!(zero_length_floor, |p: &mut L1Input| p
    .options
    .length_floor_m =
    0.0);
invalid!(nonfinite_length_floor, |p: &mut L1Input| p
    .options
    .length_floor_m =
    f64::NAN);
invalid!(zero_anchor_threshold, |p: &mut L1Input| p
    .options
    .minimum_anchor_count =
    0);
invalid!(excessive_iteration_limit, |p: &mut L1Input| p
    .options
    .max_iterations =
    Some(usize::MAX));
