use cave_closure::*;
use std::{hint::black_box, time::Instant};

fn main() {
    for size in [1000_usize, 10_000, 24_000] {
        for ring in [false, true] {
            let stations: Vec<_> = (0..size)
                .map(|i| Station {
                    id: StationId(format!("s{i:06}").into()),
                })
                .collect();
            let mut observations = Vec::new();
            for i in 0..size - 1 {
                observations.push(Observation::Displacement {
                    id: ObservationId(format!("o{i:06}").into()),
                    from: stations[i].id.clone(),
                    to: stations[i + 1].id.clone(),
                    value: [1., 0., 0.],
                    uncertainty: Uncertainty::Isotropic { sigma_m: 1. },
                    reading_ids: vec![],
                });
            }
            if ring {
                observations.push(Observation::Displacement {
                    id: "closing".into(),
                    from: stations[size - 1].id.clone(),
                    to: stations[0].id.clone(),
                    value: [-(size as f64 - 1.) + 0.2, 0., 0.],
                    uncertainty: Uncertainty::Isotropic { sigma_m: 1. },
                    reading_ids: vec![],
                });
            }
            let p = NormalizedProblem {
                schema_version: 1,
                frame: Frame {
                    kind: FrameKind::LocalEnu,
                    length_unit: LengthUnit::Metres,
                    provenance: vec!["synthetic chain/ring".into()],
                },
                mode: Mode::Spatial,
                stations,
                observations,
                controls: vec![],
                relationships: vec![],
                algorithm: Algorithm {
                    objective: Objective::Gaussian,
                    backend: LinearBackend::ReducedPcg,
                },
            };
            let mut preparation_times = Vec::new();
            let mut solve_times = Vec::new();
            for _ in 0..5 {
                let start = Instant::now();
                let prepared = prepare(black_box(&p), &Default::default()).unwrap();
                preparation_times.push(start.elapsed());
                let start = Instant::now();
                let result = prepared
                    .solve(
                        &Default::default(),
                        &mut Workspace::default(),
                        &Default::default(),
                    )
                    .unwrap();
                solve_times.push(start.elapsed());
                result.require_converged().unwrap();
                let expected = if ring { 0.04 / size as f64 } else { 0. };
                assert!((result.diagnostics.objective - expected).abs() < 1e-9);
                black_box(result);
            }
            preparation_times.sort();
            solve_times.sort();
            println!(
                "stations={size} topology={} repeats=5 preparation_median_ms={:.3} solve_median_ms={:.3}",
                if ring { "ring" } else { "chain" },
                preparation_times[2].as_secs_f64() * 1000.,
                solve_times[2].as_secs_f64() * 1000.
            );
        }
    }
}
