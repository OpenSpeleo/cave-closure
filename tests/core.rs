use cave_closure::*;
fn tree() -> NormalizedProblem {
    json::from_json(include_str!("fixtures/controlled_tree.json")).unwrap()
}
fn solve(p: &NormalizedProblem, reduce: bool) -> SolveReport {
    let mut problem = p.clone();
    if matches!(
        problem.algorithm,
        LoopClosureAlgorithm::LseReducedPcg | LoopClosureAlgorithm::LseUnreducedPcg
    ) {
        problem.algorithm = if reduce {
            LoopClosureAlgorithm::LseReducedPcg
        } else {
            LoopClosureAlgorithm::LseUnreducedPcg
        };
    }
    prepare(&problem, &PrepareOptions::default())
        .unwrap()
        .solve(
            &SolveConfig::default(),
            &mut Workspace::default(),
            &SolveControl::default(),
        )
        .unwrap()
}
#[test]
fn controlled_tree() {
    let report = solve(&tree(), true);
    report.require_converged().unwrap();
    assert!((report.coordinates[1].value[0] - 10.2).abs() < 1e-12);
    assert!((report.diagnostics.objective - 0.2).abs() < 1e-12);
    assert_eq!(report.diagnostics.cycle_rank, 0);
    assert_eq!(report.diagnostics.stochastic_redundancy, 3);
    assert_eq!(report.diagnostics.kernel_unknowns, 0);
}
#[test]
fn reduced_equals_unreduced() {
    let a = solve(&tree(), true);
    let b = solve(&tree(), false);
    for (a, b) in a.coordinates.iter().zip(b.coordinates) {
        assert_eq!(a.value, b.value);
    }
}

#[cfg(feature = "sparse-direct")]
#[test]
fn direct_solvers_preserve_exact_zero_stationarity() {
    let mut problem = tree();
    problem.stations.truncate(2);
    problem.controls.truncate(1);
    problem.observations = [("ab", "A", "B"), ("ba", "B", "A")]
        .into_iter()
        .map(|(id, from, to)| Observation::Displacement {
            id: id.into(),
            from: from.into(),
            to: to.into(),
            value: [8.660254037844386, -5., 0.],
            uncertainty: Uncertainty::Isotropic { sigma_m: 1. },
            reading_ids: vec![],
        })
        .collect();
    for backend in [
        LoopClosureAlgorithm::LseSparseQr,
        LoopClosureAlgorithm::LseSparseCholesky,
    ] {
        problem.algorithm = backend;
        let result = solve(&problem, false);
        result.require_converged().unwrap();
        assert_eq!(result.coordinates[1].value, [0.; 3]);
        assert_eq!(result.diagnostics.backward_error, 0.);
        assert!((result.diagnostics.objective - 200.).abs() < 1e-12);
    }
}
#[test]
fn self_loop_is_constant() {
    let mut p = tree();
    p.observations.push(Observation::Displacement {
        id: "check".into(),
        from: "B".into(),
        to: "B".into(),
        value: [2., 0., 0.],
        uncertainty: Uncertainty::Isotropic { sigma_m: 1. },
        reading_ids: vec![],
    });
    let r = solve(&p, true);
    assert!((r.coordinates[1].value[0] - 10.2).abs() < 1e-12);
    assert!((r.diagnostics.objective - 4.2).abs() < 1e-12);
    assert_eq!(r.diagnostics.stochastic_redundancy, 6);
}
#[test]
fn parallels_and_reversal_preserved() {
    let mut p = tree();
    let mut shot = p.observations[0].clone();
    let Observation::Displacement {
        id,
        from,
        to,
        value,
        ..
    } = &mut shot;
    *id = "reverse".into();
    std::mem::swap(from, to);
    for v in value {
        *v = -*v;
    }
    p.observations.push(shot);
    let r = solve(&p, true);
    assert_eq!(r.observations.len(), 3);
    assert!((r.coordinates[1].value[0] - (10. + 1. / 9.)).abs() < 1e-12);
}
#[test]
fn exact_offset_and_conflict() {
    let mut p = tree();
    p.relationships.push(Relationship::HardOffset {
        id: "protected".into(),
        from: "A".into(),
        to: "B".into(),
        value: [10., 0., 0.],
    });
    assert_eq!(solve(&p, true).coordinates[1].value, [10., 0., 0.]);
    p.relationships.push(Relationship::Equivalent {
        id: "bad".into(),
        from: "A".into(),
        to: "B".into(),
    });
    assert!(prepare(&p, &PrepareOptions::default()).is_err());
}
#[test]
fn datumless_tree_has_only_gauge() {
    let mut p = tree();
    p.controls.clear();
    let r = solve(&p, true);
    r.require_converged().unwrap();
    assert_eq!(r.diagnostics.artificial_gauges.len(), 3);
    assert_eq!(r.diagnostics.stochastic_redundancy, 0);
    assert!(r.diagnostics.variance_factor.is_none());
    assert!(
        prepare(
            &p,
            &PrepareOptions {
                require_control: true,
                ..Default::default()
            }
        )
        .is_err()
    );
}
#[test]
fn partial_and_soft_controls() {
    let mut p = tree();
    p.controls = vec![
        Control::Soft {
            id: "weak".into(),
            station: "A".into(),
            axes: vec![Axis::X],
            value: vec![100.],
            covariance: vec![vec![1e6]],
        },
        Control::Hard {
            id: "y".into(),
            station: "B".into(),
            axes: vec![Axis::Y],
            value: vec![3.],
        },
    ];
    let r = solve(&p, true);
    r.require_converged().unwrap();
    assert_eq!(r.diagnostics.artificial_gauges.len(), 1);
    assert!((r.coordinates[0].value[0] - 100.).abs() < 1e-6);
    assert_eq!(r.coordinates[1].value[1], 3.);
}
#[test]
fn shared_readings_rejected_or_reported() {
    let mut p = tree();
    for o in &mut p.observations {
        let Observation::Displacement { reading_ids, .. } = o;
        reading_ids.push("shared".into());
    }
    assert!(prepare(&p, &PrepareOptions::default()).is_err());
    let prepared = prepare(
        &p,
        &PrepareOptions {
            allow_shared_reading_approximation: true,
            ..Default::default()
        },
    )
    .unwrap();
    let r = prepared
        .solve(
            &Default::default(),
            &mut Workspace::default(),
            &Default::default(),
        )
        .unwrap();
    assert!(!r.diagnostics.calibrated_statistics_eligible);
    assert_eq!(r.diagnostics.warnings.len(), 1);
}
#[test]
fn horizontal_uses_marginal_covariance() {
    let mut p = tree();
    p.mode = LoopClosureMode::Horizontal {
        preserved_z: p.stations.iter().map(|s| (s.id.clone(), 42.)).collect(),
    };
    for c in &mut p.controls {
        let Control::Hard { axes, value, .. } = c else {
            unreachable!()
        };
        axes.pop();
        value.pop();
    }
    let Observation::Displacement { uncertainty, .. } = &mut p.observations[0];
    *uncertainty = Uncertainty::Covariance {
        matrix: [[1., 0., 0.8], [0., 1., 0.], [0.8, 0., 1.]],
    };
    let r = solve(&p, true);
    assert!((r.coordinates[1].value[0] - 10.2).abs() < 1e-12);
    assert_eq!(r.coordinates[1].value[2], 42.);
}
#[test]
fn rejects_invalid_covariances_and_ids() {
    let mut p = tree();
    p.stations.push(p.stations[0].clone());
    assert!(prepare(&p, &Default::default()).is_err());
    let mut p = tree();
    let Observation::Displacement { uncertainty, .. } = &mut p.observations[0];
    *uncertainty = Uncertainty::LengthProportional {
        length_m: 0.,
        variance_per_m: 1.,
    };
    assert!(prepare(&p, &Default::default()).is_err());
}
#[test]
fn iteration_limit_is_not_success() {
    let mut problem = tree();
    problem.algorithm = LoopClosureAlgorithm::LseUnreducedPcg;
    let r = prepare(&problem, &Default::default())
        .unwrap()
        .solve(
            &SolveConfig {
                max_iterations: 0,
                ..Default::default()
            },
            &mut Workspace::default(),
            &Default::default(),
        )
        .unwrap();
    assert_eq!(r.status, Termination::IterationLimit);
    assert!(r.require_converged().is_err());
}
#[test]
fn cancellation() {
    let c = SolveControl::default();
    c.cancel();
    assert!(
        prepare(&tree(), &Default::default())
            .unwrap()
            .solve(&Default::default(), &mut Workspace::default(), &c)
            .is_err()
    );
}
#[test]
fn strict_json() {
    let input = serde_json::to_string(&tree()).unwrap();
    let duplicate = input.replace("\"sigma_m\":1.0", "\"sigma_m\":1.0,\"sigma_m\":2.0");
    assert_ne!(duplicate, input);
    assert!(json::from_json(&duplicate).is_err());
    assert!(
        json::from_json(&input.replace("\"local_enu\"", "\"local_enu\",\"unknown\":true")).is_err()
    );
    assert!(json::from_json(&input.replace("{\"id\":\"A\"}", "[\"A\"]")).is_err());
}

#[test]
fn nonfinite_output_is_rejected() {
    let mut r = solve(&tree(), true);
    r.coordinates[0].value[0] = f64::NAN;
    assert!(r.to_json().is_err());
}
#[test]
fn normalized_polar_jacobian() {
    let (v, c) = normalize::polar(
        10.,
        0.,
        0.,
        [[0.01, 0., 0.], [0., 0.0001, 0.], [0., 0., 0.0004]],
    )
    .unwrap();
    assert_eq!(v, [0., 10., 0.]);
    assert!((c[0][0] - 0.01).abs() < 1e-14);
    assert!((c[2][2] - 0.04).abs() < 1e-14);
}
#[test]
fn underwater_invalid_geometry() {
    assert!(normalize::underwater(1., 0., 0., 2., [[0.; 4]; 4], 1e-6).is_err());
    assert!(normalize::underwater(1., 0., 0., 1., [[0.; 4]; 4], 1e-6).is_err());
}
#[test]
fn lost_grounding_is_failure() {
    let mut p = tree();
    p.stations.pop();
    p.observations.pop();
    p.controls = vec![Control::Soft {
        id: "weak".into(),
        station: "A".into(),
        axes: vec![Axis::X, Axis::Y, Axis::Z],
        value: vec![0.; 3],
        covariance: vec![vec![1e6, 0., 0.], vec![0., 1e6, 0.], vec![0., 0., 1e6]],
    }];
    let Observation::Displacement { uncertainty, .. } = &mut p.observations[0];
    *uncertainty = Uncertainty::Isotropic { sigma_m: 1e-6 };
    assert!(
        prepare(&p, &Default::default())
            .unwrap()
            .solve(
                &Default::default(),
                &mut Workspace::default(),
                &Default::default()
            )
            .is_err()
    );
}
#[test]
fn large_translation_keeps_local_accuracy() {
    let mut p = tree();
    for c in &mut p.controls {
        let Control::Hard { value, .. } = c else {
            unreachable!()
        };
        value[0] += 1e9;
    }
    let r = solve(&p, true);
    r.require_converged().unwrap();
    assert!((r.coordinates[1].value[0] - 1e9 - 10.2).abs() < 1e-7);
}
#[cfg(feature = "sparse-direct")]
#[test]
fn direct_backends_match_analytic_tree() {
    for backend in [
        LoopClosureAlgorithm::LseSparseQr,
        LoopClosureAlgorithm::LseSparseCholesky,
    ] {
        let mut p = tree();
        p.algorithm = backend;
        let r = solve(&p, true);
        r.require_converged().unwrap();
        assert!((r.coordinates[1].value[0] - 10.2).abs() < 1e-12);
    }
}

#[cfg(feature = "sparse-direct")]
#[test]
fn generated_noncommuting_blocks_with_partial_controls() {
    for seed in [13_u64, 71, 509] {
        let mut state = seed;
        let mut random = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 32) as u32
        };
        let mut p = tree();
        p.stations = (0..25)
            .map(|i| Station {
                id: StationId(format!("s{i:03}").into()),
            })
            .collect();
        p.controls = vec![
            Control::Hard {
                id: "origin".into(),
                station: "s000".into(),
                axes: vec![Axis::X, Axis::Y, Axis::Z],
                value: vec![0.; 3],
            },
            Control::Hard {
                id: "end".into(),
                station: "s024".into(),
                axes: vec![Axis::X, Axis::Z],
                value: vec![20., 3.],
            },
            Control::Hard {
                id: "partial".into(),
                station: "s012".into(),
                axes: vec![Axis::Y],
                value: vec![2.],
            },
        ];
        p.observations.clear();
        for i in 0..90 {
            let u = if i < 24 { i } else { random() as usize % 25 };
            let v = if i < 24 {
                i + 1
            } else {
                random() as usize % 25
            };
            let direction: [f64; 3] =
                std::array::from_fn(|_| f64::from(random() % 200) / 100. - 1.);
            let matrix = std::array::from_fn(|a| {
                std::array::from_fn(|b| {
                    direction[a] * direction[b]
                        + if a == b {
                            0.5 + f64::from(random() % 100) / 100.
                        } else {
                            0.
                        }
                })
            });
            let value = std::array::from_fn(|_| f64::from(random() % 100) / 10. - 5.);
            p.observations.push(Observation::Displacement {
                id: ObservationId(format!("o{i:03}").into()),
                from: p.stations[u].id.clone(),
                to: p.stations[v].id.clone(),
                value,
                uncertainty: Uncertainty::Covariance { matrix },
                reading_ids: vec![],
            });
        }
        p.algorithm = LoopClosureAlgorithm::LseSparseQr;
        let oracle = solve(&p, false);
        oracle.require_converged().unwrap();
        for backend in [
            LoopClosureAlgorithm::LseReducedPcg,
            LoopClosureAlgorithm::LseSparseCholesky,
        ] {
            p.algorithm = backend;
            for reduce in [false, true] {
                let actual = solve(&p, reduce);
                actual.require_converged().unwrap();
                for (a, b) in actual.coordinates.iter().zip(&oracle.coordinates) {
                    for axis in 0..3 {
                        assert!((a.value[axis] - b.value[axis]).abs() < 1e-9);
                    }
                }
                assert!((actual.diagnostics.objective - oracle.diagnostics.objective).abs() < 1e-9);
            }
        }
    }
}

#[test]
fn propagated_covariance_with_correlated_raw_readings() {
    let raw = [
        [0.04, 0.001, 0.0002],
        [0.001, 0.003, 0.0001],
        [0.0002, 0.0001, 0.002],
    ];
    let params = [12., 0.8, 0.4];
    let (_, covariance) = normalize::polar(params[0], params[1], params[2], raw).unwrap();
    let mut jacobian = [[0.; 3]; 3];
    for a in 0..3 {
        let mut plus = params;
        let mut minus = params;
        let step = 1e-5;
        plus[a] += step;
        minus[a] -= step;
        let (p, _) = normalize::polar(plus[0], plus[1], plus[2], raw).unwrap();
        let (m, _) = normalize::polar(minus[0], minus[1], minus[2], raw).unwrap();
        for i in 0..3 {
            jacobian[i][a] = (p[i] - m[i]) / (2. * step);
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let mut expected = 0.;
            for (a, row) in raw.iter().enumerate() {
                for (b, value) in row.iter().enumerate() {
                    expected += jacobian[i][a] * value * jacobian[j][b];
                }
            }
            assert!((covariance[i][j] - expected).abs() < 1e-8);
        }
    }
    let raw = [
        [0.04, 0.001, 0., 0.],
        [0.001, 0.002, 0., 0.],
        [0., 0., 0.01, 0.005],
        [0., 0., 0.005, 0.02],
    ];
    let (value, covariance) = normalize::underwater(12., 0.8, 2., 5., raw, 1e-6).unwrap();
    assert_eq!(value[2], -3.);
    assert!((covariance[2][2] - 0.02).abs() < 1e-14);
}
