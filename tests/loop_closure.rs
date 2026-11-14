//! Native consumer-boundary regressions and independent analytical solutions.
use cave_closure::loop_closure::*;
use cave_closure::*;
use std::sync::Arc;

fn input() -> LoopClosureInput {
    LoopClosureInput {
        station_ids: vec!["A".into(), "B".into(), "C".into()],
        measurements: vec![
            ("AB".into(), 0, 1, 5.0, 90.0, 0.0, None, 1.0, false),
            ("BC".into(), 1, 2, 15.0, 90.0, 0.0, None, 1.0, false),
        ],
        initial_positions: vec![[0.0, 0.0, 0.0], [5.0, 0.0, -3.0], [24.0, 0.0, -8.0]],
        controls: vec![
            ("explicit:A".into(), 0, vec![0, 1], vec![0.0, 0.0]),
            ("explicit:C".into(), 2, vec![0, 1], vec![24.0, 0.0]),
        ],
        relationships: vec![],
        source_order: vec![0, 1, 2],
        options: LoopClosureOptions {
            mode: "horizontal".into(),
            algorithm: LoopClosureAlgorithm::LseUnreducedPcg,
            weighting: "independent_cartesian".into(),
            method: "least_squares".into(),
            sigma_m: 1.0,
            weight_cap: 10.0,
            sigma_length_m: 0.05,
            sigma_azimuth_deg: 1.0,
            average_azimuths: false,
            anchor_floating: false,
            preserve_centroids: false,
            fix_control_depths: false,
            max_iterations: 60_000,
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-5,
        },
    }
}
fn run(input: LoopClosureInput) -> SolveReport {
    prepare_loop_closure(input, &SolveControl::default())
        .unwrap()
        .solve(None, &SolveControl::default())
        .unwrap()
}
fn point(report: &SolveReport, id: &str) -> Point {
    report
        .coordinates
        .iter()
        .find(|p| p.station.0.as_ref() == id)
        .unwrap()
        .value
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-7, "{actual} != {expected}");
}
fn graph(reference: bool) -> LoopClosureInput {
    let mut p = input();
    p.options.method = if reference {
        "reference_cg"
    } else {
        "traverse_quality_cg"
    }
    .into();
    p.options.weighting = if reference {
        "inverse_length"
    } else {
        "traverse_quality_inverse_squared_length"
    }
    .into();
    p.options.absolute_tolerance = if reference { 0.001 } else { 1e-8 };
    p
}
#[test]
fn equal_variance_distributes_error_equally() {
    let r = run(input());
    near(point(&r, "B")[0], 7.0);
    near(r.diagnostics.objective, 8.0);
}
#[test]
fn inverse_length_distributes_error_by_length() {
    let mut p = input();
    p.options.weighting = "inverse_length".into();
    let r = run(p);
    near(point(&r, "B")[0], 6.0);
    near(r.diagnostics.objective, 0.8);
    assert!(!r.diagnostics.calibrated_statistics_eligible);
}
#[test]
fn reference_cg_solves_the_same_weighted_objective() {
    let r = run(graph(true));
    near(point(&r, "B")[0], 6.0);
    near(r.diagnostics.objective, 0.8);
}
#[test]
fn traverse_quality_uses_squared_lengths() {
    let r = run(graph(false));
    near(point(&r, "B")[0], 5.4);
}
#[test]
fn instrument_length_variance_controls_east_axis() {
    let mut p = input();
    p.options.weighting = "length_azimuth_diagonal".into();
    let r = run(p);
    near(point(&r, "B")[0], 7.0);
    assert!(!r.diagnostics.calibrated_statistics_eligible);
}
#[test]
fn instrument_quality_squares_the_sigma_multiplier() {
    let mut p = input();
    p.options.weighting = "length_azimuth_diagonal".into();
    p.measurements[0].7 = 2.0;
    near(point(&run(p), "B")[0], 8.2);
}
#[test]
fn invalid_section_quality_falls_back_to_one() {
    let mut p = input();
    p.options.weighting = "length_azimuth_diagonal".into();
    p.measurements[0].7 = f64::NAN;
    near(point(&run(p), "B")[0], 7.0);
}
#[test]
fn horizontal_adjustment_preserves_each_recorded_depth() {
    for p in [input(), graph(true), graph(false)] {
        let r = run(p);
        assert_eq!(point(&r, "B")[2], -3.0);
        assert_eq!(point(&r, "C")[2], -8.0);
    }
}
#[test]
fn spatial_gaussian_adjusts_vertical_measurements() {
    let mut p = input();
    p.options.mode = "spatial".into();
    p.options.fix_control_depths = true;
    p.initial_positions[2] = [0.0, 0.0, 12.0];
    p.controls[1].3 = vec![0.0, 0.0];
    for row in &mut p.measurements {
        row.3 = 5.0;
        row.4 = 0.0;
        row.5 = 5.0;
    }
    near(point(&run(p), "B")[2], 6.0);
}
#[test]
fn spatial_graph_adjusts_vertical_measurements() {
    let mut p = graph(false);
    p.options.mode = "spatial".into();
    p.options.fix_control_depths = true;
    p.initial_positions[2] = [0.0, 0.0, 12.0];
    p.controls[1].3 = vec![0.0, 0.0];
    for row in &mut p.measurements {
        row.3 = 5.0;
        row.4 = 0.0;
        row.5 = 5.0;
    }
    near(point(&run(p), "B")[2], 6.0);
}
#[test]
fn supplied_vertical_controls_override_recorded_control_depths() {
    let mut p = graph(false);
    p.options.mode = "spatial".into();
    p.options.fix_control_depths = true;
    p.controls.push(("z:A".into(), 0, vec![2], vec![2.0]));
    p.controls.push(("z:C".into(), 2, vec![2], vec![10.0]));
    let r = run(p);
    near(point(&r, "A")[2], 2.0);
    near(point(&r, "C")[2], 10.0);
}
#[test]
fn spatial_gaussian_requires_vertical_control() {
    let mut p = input();
    p.options.mode = "spatial".into();
    assert!(prepare_loop_closure(p, &SolveControl::default()).is_err());
}
#[test]
fn graph_with_one_anchor_returns_unadjusted_positions() {
    let mut p = graph(false);
    p.controls.truncate(1);
    let expected = p.initial_positions.clone();
    let r = run(p);
    for (i, id) in ["A", "B", "C"].into_iter().enumerate() {
        assert_eq!(point(&r, id), expected[i]);
    }
}
#[test]
fn graph_with_no_anchor_returns_unadjusted_positions() {
    let mut p = graph(false);
    p.controls.clear();
    let r = run(p);
    assert_eq!(point(&r, "B"), [5.0, 0.0, -3.0]);
}
#[test]
fn reference_cg_adjusts_with_only_one_anchor() {
    let mut p = graph(true);
    p.controls.truncate(1);
    let r = run(p);
    near(point(&r, "C")[0], 20.0);
}
#[test]
fn all_anchored_graph_preserves_controls() {
    let mut p = graph(false);
    p.controls
        .push(("explicit:B".into(), 1, vec![0, 1], vec![12.0, 7.0]));
    let r = run(p);
    assert_eq!(&point(&r, "B")[..2], &[12.0, 7.0]);
}
#[cfg(feature = "sparse-direct")]
#[test]
fn sparse_direct_fallback_matches_known_solution() {
    let p = prepare_loop_closure(graph(false), &SolveControl::default()).unwrap();
    let r = p.solve(Some(0), &SolveControl::default()).unwrap();
    near(point(&r, "B")[0], 5.4);
    assert!(
        r.diagnostics
            .warnings
            .iter()
            .any(|w| w.contains("fallback"))
    );
}
#[test]
fn reference_zero_iteration_exit_returns_initial_guess() {
    let p = prepare_loop_closure(graph(true), &SolveControl::default()).unwrap();
    let r = p.solve(Some(0), &SolveControl::default()).unwrap();
    near(point(&r, "B")[0], 5.0);
    assert!(
        r.diagnostics
            .warnings
            .iter()
            .any(|w| w.contains("current iterate"))
    );
}
#[test]
fn graph_weight_floor_keeps_zero_length_observations_finite() {
    let mut p = graph(false);
    p.measurements[0].3 = 0.0;
    let r = run(p);
    assert!(
        r.coordinates
            .iter()
            .flat_map(|p| p.value)
            .all(f64::is_finite)
    );
}
#[test]
fn reference_zero_length_cap_is_finite_and_configurable() {
    let mut p = graph(true);
    p.measurements[0].3 = 0.0;
    p.options.weight_cap = 2.0;
    near(point(&run(p), "B")[0], 9.0 / 31.0);
}
#[test]
fn station_equivalence_closes_a_loop() {
    let mut p = input();
    p.controls.truncate(1);
    p.relationships
        .push(("loop".into(), "equivalent".into(), 0, 2, [0.0; 3]));
    let r = run(p);
    near(point(&r, "C")[0], 0.0);
    near(point(&r, "B")[0], -5.0);
}
#[test]
fn graph_equivalence_merges_anchors_before_bypass() {
    let mut p = graph(false);
    p.controls.truncate(1);
    p.relationships
        .push(("loop".into(), "equivalent".into(), 0, 2, [0.0; 3]));
    let r = run(p);
    near(point(&r, "C")[0], 0.0);
    near(point(&r, "B")[0], 5.0);
}
#[test]
fn conflicting_geographic_controls_are_rejected_after_equivalence() {
    let mut p = graph(false);
    p.relationships
        .push(("loop".into(), "equivalent".into(), 0, 2, [0.0; 3]));
    assert!(prepare_loop_closure(p, &SolveControl::default()).is_err());
}
#[test]
fn closure_payload_does_not_become_a_measurement() {
    let mut p = input();
    p.controls.truncate(1);
    p.relationships.push((
        "loop".into(),
        "equivalent".into(),
        0,
        2,
        [900.0, 800.0, 700.0],
    ));
    let r = run(p);
    near(point(&r, "C")[0], 0.0);
    assert_eq!(r.observations.len(), 2);
}
#[test]
fn virtual_hard_offset_is_an_exact_constraint() {
    let mut p = input();
    p.measurements[0].8 = true;
    let r = run(p);
    near(point(&r, "B")[0], 5.0);
    assert_eq!(r.observations.len(), 1);
}
#[test]
fn source_compatible_graph_rejects_virtual_hard_offsets() {
    let mut p = graph(true);
    p.measurements[0].8 = true;
    assert!(prepare_loop_closure(p, &SolveControl::default()).is_err());
}
#[test]
fn reversed_measurement_orientation_preserves_geometry() {
    let expected = run(input());
    let mut p = input();
    for row in &mut p.measurements {
        std::mem::swap(&mut row.1, &mut row.2);
        row.4 += 180.0;
        row.5 = -row.5;
    }
    near(point(&run(p), "B")[0], point(&expected, "B")[0]);
}
#[test]
fn secondary_bearing_averages_across_north() {
    let mut p = input();
    p.controls.truncate(1);
    p.options.average_azimuths = true;
    p.measurements.truncate(1);
    p.measurements[0].3 = 10.0;
    p.measurements[0].4 = 350.0;
    p.measurements[0].6 = Some(10.0);
    p.station_ids.pop();
    p.initial_positions.pop();
    p.source_order.pop();
    let r = run(p);
    near(point(&r, "B")[0], 0.0);
    near(point(&r, "B")[1], 10.0);
}
#[test]
fn secondary_bearing_is_ignored_unless_requested() {
    let mut p = input();
    p.measurements[0].6 = Some(0.0);
    near(point(&run(p), "B")[0], 7.0);
}
#[test]
fn zero_instrument_axis_retains_measured_centroid() {
    let mut p = input();
    p.options.weighting = "length_azimuth_diagonal".into();
    p.controls.truncate(1);
    p.measurements.clear();
    p.station_ids = vec!["A".into(), "B".into()];
    p.initial_positions = vec![[0.0, 0.0, 0.0], [17.0, 23.0, -3.0]];
    p.source_order = vec![0, 1];
    p.options.anchor_floating = true;
    let r = run(p);
    near(point(&r, "B")[0], 17.0);
    near(point(&r, "B")[1], 23.0);
}
#[test]
fn floating_gaussian_component_requires_explicit_policy() {
    let mut p = input();
    p.controls.clear();
    assert!(prepare_loop_closure(p, &SolveControl::default()).is_err());
}
#[test]
fn floating_component_centroid_is_preserved_when_requested() {
    let mut p = input();
    p.controls.clear();
    p.options.preserve_centroids = true;
    let before = p.initial_positions.iter().map(|p| p[0]).sum::<f64>();
    let r = run(p);
    near(
        r.coordinates.iter().map(|p| p.value[0]).sum::<f64>(),
        before,
    );
}
#[test]
fn component_anchor_follows_source_order() {
    let mut p = input();
    p.controls.clear();
    p.options.anchor_floating = true;
    p.source_order = vec![1, 0, 2];
    let r = run(p);
    near(point(&r, "B")[0], 5.0);
    near(point(&r, "C")[0], 20.0);
}
#[test]
fn input_is_owned_after_preparation() {
    let mut p = input();
    let owned = prepare_loop_closure(p.clone(), &SolveControl::default()).unwrap();
    p.measurements[0].3 = 10000.0;
    near(
        point(&owned.solve(None, &SolveControl::default()).unwrap(), "B")[0],
        7.0,
    );
}
#[test]
fn shared_preparation_has_independent_concurrent_workspaces() {
    let p = Arc::new(prepare_loop_closure(input(), &SolveControl::default()).unwrap());
    let workers: Vec<_> = (0..8)
        .map(|_| {
            let p = Arc::clone(&p);
            std::thread::spawn(move || {
                p.solve(None, &SolveControl::default())
                    .unwrap()
                    .to_json()
                    .unwrap()
            })
        })
        .collect();
    let results: Vec<_> = workers.into_iter().map(|w| w.join().unwrap()).collect();
    assert!(results.windows(2).all(|v| v[0] == v[1]));
}
#[test]
fn cancellation_is_observed_by_all_adjustment_methods() {
    for p in [input(), graph(true), graph(false)] {
        let prepared = prepare_loop_closure(p, &SolveControl::default()).unwrap();
        let control = SolveControl::default();
        control.cancel();
        assert!(matches!(
            prepared.solve(None, &control),
            Err(Error::Incomplete(Termination::Cancelled))
        ));
    }
}
#[test]
fn scalar_instrument_solves_observe_cancellation() {
    let mut p = input();
    p.options.weighting = "length_azimuth_diagonal".into();
    let p = prepare_loop_closure(p, &SolveControl::default()).unwrap();
    let control = SolveControl::default();
    control.cancel();
    assert!(matches!(
        p.solve(None, &control),
        Err(Error::Incomplete(Termination::Cancelled))
    ));
}
#[test]
fn json_request_roundtrip_preserves_solution() {
    let p = input();
    let encoded = serde_json::to_string(&p).unwrap();
    let decoded = serde_json::from_str(&encoded).unwrap();
    near(point(&run(decoded), "B")[0], 7.0);
}
#[test]
fn json_rejects_unknown_configuration_fields() {
    let mut p = serde_json::to_value(input()).unwrap();
    p["options"]["typo"] = true.into();
    assert!(serde_json::from_value::<LoopClosureInput>(p).is_err());
}

macro_rules! invalid_request {
    ($name:ident,$change:expr) => {
        #[test]
        fn $name() {
            let mut p = input();
            ($change)(&mut p);
            assert!(prepare_loop_closure(p, &SolveControl::default()).is_err());
        }
    };
}
invalid_request!(duplicate_station_id, |p: &mut LoopClosureInput| p
    .station_ids[1] =
    "A".into());
invalid_request!(empty_station_id, |p: &mut LoopClosureInput| p.station_ids
    [1]
.clear());
invalid_request!(missing_initial_position, |p: &mut LoopClosureInput| {
    p.initial_positions.pop();
});
invalid_request!(nonfinite_initial_position, |p: &mut LoopClosureInput| p
    .initial_positions[1][0] =
    f64::NAN);
invalid_request!(duplicate_source_order_index, |p: &mut LoopClosureInput| p
    .source_order[1] =
    0);
invalid_request!(missing_source_order_index, |p: &mut LoopClosureInput| {
    p.source_order.pop();
});
invalid_request!(invalid_measurement_endpoint, |p: &mut LoopClosureInput| p
    .measurements[0]
    .2 =
    1000);
invalid_request!(duplicate_measurement_id, |p: &mut LoopClosureInput| p
    .measurements[1]
    .0 =
    "AB".into());
invalid_request!(negative_measured_length, |p: &mut LoopClosureInput| p
    .measurements[0]
    .3 =
    -1.0);
invalid_request!(
    length_shorter_than_depth_change,
    |p: &mut LoopClosureInput| p.measurements[0].5 = 6.0
);
invalid_request!(nonfinite_bearing, |p: &mut LoopClosureInput| p
    .measurements[0]
    .4 =
    f64::INFINITY);
invalid_request!(nonfinite_secondary_bearing, |p: &mut LoopClosureInput| p
    .measurements[0]
    .6 =
    Some(f64::NAN));
invalid_request!(control_unknown_station, |p: &mut LoopClosureInput| p
    .controls[0]
    .1 =
    1000);
invalid_request!(control_invalid_axis, |p: &mut LoopClosureInput| p
    .controls[0]
    .2[0] =
    3);
invalid_request!(
    control_mismatched_axis_values,
    |p: &mut LoopClosureInput| {
        p.controls[0].3.pop();
    }
);
invalid_request!(control_nonfinite_value, |p: &mut LoopClosureInput| p
    .controls[0]
    .3[0] =
    f64::INFINITY);
invalid_request!(relationship_unknown_station, |p: &mut LoopClosureInput| p
    .relationships
    .push(("loop".into(), "equivalent".into(), 0, 10, [0.0; 3])));
invalid_request!(relationship_unknown_kind, |p: &mut LoopClosureInput| p
    .relationships
    .push(("loop".into(), "typo".into(), 0, 1, [0.0; 3])));
invalid_request!(unknown_weighting, |p: &mut LoopClosureInput| p
    .options
    .weighting =
    "typo".into());
invalid_request!(unknown_method, |p: &mut LoopClosureInput| p
    .options
    .method =
    "typo".into());
invalid_request!(unknown_coordinate_mode, |p: &mut LoopClosureInput| p
    .options
    .mode =
    "typo".into());
invalid_request!(zero_precision_cap, |p: &mut LoopClosureInput| p
    .options
    .weight_cap =
    0.0);
invalid_request!(nonfinite_precision_cap, |p: &mut LoopClosureInput| p
    .options
    .weight_cap =
    f64::NAN);
invalid_request!(negative_length_sigma, |p: &mut LoopClosureInput| p
    .options
    .sigma_length_m =
    -1.0);
invalid_request!(zero_angle_sigma, |p: &mut LoopClosureInput| p
    .options
    .sigma_azimuth_deg =
    0.0);
invalid_request!(zero_absolute_tolerance, |p: &mut LoopClosureInput| p
    .options
    .absolute_tolerance =
    0.0);
invalid_request!(negative_relative_tolerance, |p: &mut LoopClosureInput| p
    .options
    .relative_tolerance =
    -1.0);
invalid_request!(
    zero_configured_iteration_limit,
    |p: &mut LoopClosureInput| p.options.max_iterations = 0
);
invalid_request!(
    traverse_weighting_requires_matching_method,
    |p: &mut LoopClosureInput| p.options.weighting =
        "traverse_quality_inverse_squared_length".into()
);
invalid_request!(
    instrument_weighting_rejects_spatial_mode,
    |p: &mut LoopClosureInput| {
        p.options.weighting = "length_azimuth_diagonal".into();
        p.options.mode = "spatial".into();
    }
);

#[test]
fn tiny_nonzero_rhs_must_satisfy_backward_error_before_early_exit() {
    let mut p = input();
    p.station_ids.pop();
    p.initial_positions.pop();
    p.source_order.pop();
    p.controls.truncate(1);
    p.measurements.truncate(1);
    p.measurements[0].3 = 1e-20;
    let report = run(p);
    report.require_converged().unwrap();
    let x = point(&report, "B")[0];
    assert!(x > 0.0);
    assert!(
        (x - 1e-20).abs() < 1e-30,
        "tiny measurement was dropped: {x}"
    );
    assert!(report.diagnostics.backward_error <= SolveConfig::default().backward_error_target);
}
