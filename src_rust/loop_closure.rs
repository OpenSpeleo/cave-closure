//! Owned loop-closure requests. Weighting, measurement vectors and gauges live here.
use crate::graph::{GraphConfig, GraphProblem};
use crate::*;
use std::collections::{BTreeMap, BTreeSet};

pub type MeasurementRow = (String, usize, usize, f64, f64, f64, Option<f64>, f64, bool);
pub type LoopClosureControl = (String, usize, Vec<usize>, Vec<f64>);
pub type LoopClosureRelationship = (String, String, usize, usize, Point);

#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoopClosureOptions {
    pub mode: String,
    pub algorithm: LoopClosureAlgorithm,
    pub weighting: String,
    pub method: String,
    pub sigma_m: f64,
    pub weight_cap: f64,
    pub sigma_length_m: f64,
    pub sigma_azimuth_deg: f64,
    pub average_azimuths: bool,
    pub anchor_floating: bool,
    pub preserve_centroids: bool,
    pub fix_control_depths: bool,
    pub max_iterations: usize,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
}
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoopClosureInput {
    pub station_ids: Vec<String>,
    pub measurements: Vec<MeasurementRow>,
    pub initial_positions: Vec<Point>,
    pub controls: Vec<LoopClosureControl>,
    pub relationships: Vec<LoopClosureRelationship>,
    pub source_order: Vec<usize>,
    pub options: LoopClosureOptions,
}
#[derive(Clone)]
pub enum LoopClosurePreparation {
    Standard(Box<PreparedProblem>),
    Axes {
        axes: Vec<PreparedProblem>,
        depths: Vec<(StationId, f64)>,
    },
    Graph(GraphProblem),
}
impl LoopClosurePreparation {
    pub fn solve(
        &self,
        max_iterations: Option<usize>,
        control: &SolveControl,
    ) -> Result<SolveReport, Error> {
        let config = SolveConfig {
            max_iterations: max_iterations.unwrap_or(10_000),
            ..Default::default()
        };
        match self {
            Self::Standard(p) => p.solve(&config, &mut Workspace::default(), control),
            Self::Graph(p) => p.solve(max_iterations, control),
            Self::Axes { axes, depths } => {
                let mut reports = axes
                    .iter()
                    .map(|p| p.solve(&config, &mut Workspace::default(), control))
                    .collect::<Result<Vec<_>, _>>()?;
                let y = reports.pop().unwrap();
                let mut result = reports.pop().unwrap();
                result.require_converged()?;
                y.require_converged()?;
                let y_values: BTreeMap<_, _> = y
                    .coordinates
                    .iter()
                    .map(|p| (&p.station, p.value[0]))
                    .collect();
                let z: BTreeMap<_, _> = depths.iter().map(|(id, v)| (id, *v)).collect();
                for point in &mut result.coordinates {
                    point.value = [point.value[0], y_values[&point.station], z[&point.station]];
                }
                let mut residuals = BTreeMap::<ObservationId, Point>::new();
                for observation in &result.observations {
                    residuals.entry(observation.id.clone()).or_default()[0] =
                        observation.residual[0];
                }
                for observation in &y.observations {
                    residuals.entry(observation.id.clone()).or_default()[1] =
                        observation.residual[0];
                }
                result.observations = residuals
                    .into_iter()
                    .map(|(id, residual)| ObservationResult {
                        id,
                        residual,
                        correction: residual,
                    })
                    .collect();
                result.diagnostics.objective += y.diagnostics.objective;
                result.diagnostics.calibrated_statistics_eligible = false;
                result.diagnostics.variance_factor = None;
                result
                    .diagnostics
                    .warnings
                    .push("Joint calibrated statistics are not reported".into());
                Ok(result)
            }
        }
    }
}
fn root(parents: &mut [usize], mut i: usize) -> usize {
    while parents[i] != i {
        parents[i] = parents[parents[i]];
        i = parents[i];
    }
    i
}
fn groups(n: usize, edges: impl IntoIterator<Item = (usize, usize)>) -> Vec<Vec<usize>> {
    let mut parents: Vec<_> = (0..n).collect();
    for (u, v) in edges {
        let (u, v) = (root(&mut parents, u), root(&mut parents, v));
        parents[u] = v;
    }
    let mut result = BTreeMap::<usize, Vec<usize>>::new();
    for i in 0..n {
        result.entry(root(&mut parents, i)).or_default().push(i);
    }
    result.into_values().collect()
}
fn axis(a: usize) -> Axis {
    [Axis::X, Axis::Y, Axis::Z][a]
}
fn control(row: &LoopClosureControl, ids: &[StationId]) -> Control {
    Control::Hard {
        id: row.0.clone(),
        station: ids[row.1].clone(),
        axes: row.2.iter().copied().map(axis).collect(),
        value: row.3.clone(),
    }
}
fn relationship(row: &LoopClosureRelationship, ids: &[StationId]) -> Relationship {
    if row.1 == "equivalent" {
        Relationship::Equivalent {
            id: RelationshipId(row.0.clone().into()),
            from: ids[row.2].clone(),
            to: ids[row.3].clone(),
        }
    } else {
        Relationship::HardOffset {
            id: RelationshipId(row.0.clone().into()),
            from: ids[row.2].clone(),
            to: ids[row.3].clone(),
            value: row.4,
        }
    }
}
fn problem(
    input: &LoopClosureInput,
    ids: &[StationId],
    observations: Vec<Observation>,
) -> NormalizedProblem {
    NormalizedProblem {
        schema_version: 1,
        frame: Frame {
            kind: FrameKind::LocalEnu,
            length_unit: LengthUnit::Metres,
            provenance: vec![
                "Owned measurement adjustment; weighting and preparation performed in Rust".into(),
            ],
        },
        mode: if input.options.mode == "horizontal" {
            LoopClosureMode::Horizontal {
                preserved_z: ids
                    .iter()
                    .cloned()
                    .zip(input.initial_positions.iter().map(|p| p[2]))
                    .collect(),
            }
        } else {
            LoopClosureMode::Spatial
        },
        stations: ids.iter().cloned().map(|id| Station { id }).collect(),
        observations,
        controls: input.controls.iter().map(|r| control(r, ids)).collect(),
        relationships: input
            .relationships
            .iter()
            .map(|r| relationship(r, ids))
            .collect(),
        algorithm: input.options.algorithm,
    }
}
fn prepare_problem(p: &NormalizedProblem, require_control: bool) -> Result<PreparedProblem, Error> {
    let options = PrepareOptions {
        require_control,
        ..Default::default()
    };
    // Retain the historical relative-weight interchange rounding contract.
    if p.observations
        .iter()
        .any(|o| matches!(o.parts().4, Uncertainty::Relative { .. }))
    {
        let source: NormalizedProblem = serde_json::from_str(&serde_json::to_string(p)?)?;
        prepare(&source, &options)
    } else {
        prepare(p, &options)
    }
}
/// Scaled compensated three-vector norm, retaining the source norm's rounding.
/// Uses error-free products/sums and one differential square-root correction.
fn measured_norm(value: Point) -> f64 {
    let maximum = value.iter().copied().map(f64::abs).fold(0.0, f64::max);
    if maximum == 0.0 || !maximum.is_finite() {
        return maximum;
    }
    let exponent = ((maximum.to_bits() >> 52) & 0x7ff) as i32 - 1022;
    if exponent < -1023 {
        return f64::MIN_POSITIVE * measured_norm(value.map(|x| x / f64::MIN_POSITIVE));
    }
    let scale = 2.0_f64.powi(-exponent);
    let (mut sum, mut fraction1, mut fraction2) = (1.0, 0.0, 0.0);
    for x in value {
        let x = x * scale;
        let product = x * x;
        let product_error = x.mul_add(x, -product);
        let combined = sum + product;
        let addition_error = (sum - combined) + product;
        sum = combined;
        fraction1 += product_error;
        fraction2 += addition_error;
    }
    let mut result = (sum - 1.0 + (fraction1 + fraction2)).sqrt();
    let product = -result * result;
    let product_error = (-result).mul_add(result, -product);
    let combined = sum + product;
    let addition_error = (sum - combined) + product;
    sum = combined;
    fraction1 += product_error;
    fraction2 += addition_error;
    let residual = sum - 1.0 + (fraction1 + fraction2);
    result += residual / (2.0 * result);
    result / scale
}

fn reading_square(value: f64) -> f64 {
    // Instrument propagation historically uses libm power, whose rounding is
    // distinct from a compiler-folded multiplication on some platforms.
    value.powf(std::hint::black_box(2.0))
}

fn validate(input: &LoopClosureInput) -> Result<(), Error> {
    let n = input.station_ids.len();
    let o = &input.options;
    if n == 0
        || n > 1_000_000
        || input.measurements.len() > 4_000_000
        || input.relationships.len() > 4_000_000
    {
        return Err(invalid("configured input resource limit exceeded"));
    }
    if input.initial_positions.len() != n
        || input.source_order.len() != n
        || input.source_order.iter().copied().collect::<BTreeSet<_>>() != (0..n).collect()
    {
        return Err(invalid("positions/source order must cover every station"));
    }
    if input.station_ids.iter().any(|id| id.is_empty())
        || input.station_ids.iter().collect::<BTreeSet<_>>().len() != n
        || !finite(input.initial_positions.iter().flatten().copied())
    {
        return Err(invalid("invalid station identities or initial positions"));
    }
    if !["horizontal", "spatial"].contains(&o.mode.as_str())
        || ![
            "independent_cartesian",
            "inverse_length",
            "length_azimuth_diagonal",
            "traverse_quality_inverse_squared_length",
        ]
        .contains(&o.weighting.as_str())
    {
        return Err(invalid("unknown adjustment mode or weighting"));
    }
    if !["least_squares", "reference_cg", "traverse_quality_cg"].contains(&o.method.as_str()) {
        return Err(invalid("unknown graph method"));
    }
    if !finite([
        o.sigma_m,
        o.weight_cap,
        o.sigma_length_m,
        o.sigma_azimuth_deg,
        o.absolute_tolerance,
        o.relative_tolerance,
    ]) || o.sigma_m <= 0.0
        || o.weight_cap <= 0.0
        || o.sigma_length_m <= 0.0
        || o.sigma_azimuth_deg <= 0.0
        || o.absolute_tolerance <= 0.0
        || o.relative_tolerance < 0.0
        || o.max_iterations == 0
    {
        return Err(invalid(
            "invalid adjustment uncertainty or stopping settings",
        ));
    }
    if o.weighting != "independent_cartesian" && o.sigma_m != 1.0 {
        return Err(invalid(
            "sigma_m only applies to independent Cartesian weighting",
        ));
    }
    if o.weighting == "length_azimuth_diagonal" && o.mode != "horizontal" {
        return Err(invalid(
            "Length/azimuth diagonal weighting is horizontal-only",
        ));
    }
    if (o.weighting == "traverse_quality_inverse_squared_length")
        != (o.method == "traverse_quality_cg")
    {
        return Err(invalid(
            "traverse-quality weighting requires its native graph method",
        ));
    }
    if o.method != "least_squares" && o.algorithm != LoopClosureAlgorithm::LseUnreducedPcg {
        return Err(invalid(
            "source-compatible graph methods require unreduced CG",
        ));
    }
    if o.method == "reference_cg" && (o.weighting != "inverse_length" || o.mode != "horizontal") {
        return Err(invalid(
            "reference CG requires horizontal inverse-length weights",
        ));
    }
    let mut identities = BTreeSet::new();
    for (id, u, v, length, bearing, vertical, secondary, _, offset) in &input.measurements {
        if *u >= n
            || *v >= n
            || id.is_empty()
            || !identities.insert(id)
            || !finite([*length, *bearing, *vertical])
            || secondary.is_some_and(|v| !v.is_finite())
            || *length < vertical.abs()
        {
            return Err(invalid("invalid measurement, station index or shot length"));
        }
        if *offset && o.method != "least_squares" {
            return Err(invalid(
                "graph loop closure requires measured VIRTUAL shots",
            ));
        }
    }
    for (_, u, axes, values) in &input.controls {
        if *u >= n
            || axes.is_empty()
            || axes.len() != values.len()
            || axes.iter().any(|a| *a > 2)
            || !finite(values.iter().copied())
        {
            return Err(invalid("invalid adjustment control"));
        }
    }
    for (_, kind, u, v, value) in &input.relationships {
        if *u >= n
            || *v >= n
            || !["equivalent", "hard_offset"].contains(&kind.as_str())
            || !finite(*value)
            || (kind == "hard_offset" && o.method != "least_squares")
        {
            return Err(invalid("invalid graph relationship"));
        }
    }
    Ok(())
}
fn add_centroid_controls(
    input: &LoopClosureInput,
    ids: &[StationId],
    p: &mut NormalizedProblem,
    floating: &[Vec<usize>],
    axis_index: Option<usize>,
    cancellation: &SolveControl,
) -> Result<(), Error> {
    if floating.is_empty() {
        return Ok(());
    }
    let report = prepare_problem(p, false)?.solve(
        &SolveConfig::default(),
        &mut Workspace::default(),
        cancellation,
    )?;
    report.require_converged()?;
    let coordinates: BTreeMap<_, _> = report
        .coordinates
        .iter()
        .map(|p| (&p.station, p.value))
        .collect();
    for group in floating {
        let first = group[0];
        let axes: Vec<_> = if axis_index.is_some() {
            vec![0]
        } else {
            vec![0, 1]
        };
        let values = axes
            .iter()
            .map(|&a| {
                coordinates[&ids[first]][a]
                    + group
                        .iter()
                        .map(|&i| {
                            input.initial_positions[i][axis_index.unwrap_or(a)]
                                - coordinates[&ids[i]][a]
                        })
                        .sum::<f64>()
                        / group.len() as f64
            })
            .collect();
        p.controls.push(Control::Hard {
            id: format!("centroid-gauge:{first}"),
            station: ids[first].clone(),
            axes: axes.into_iter().map(axis).collect(),
            value: values,
        });
    }
    Ok(())
}

pub fn prepare_loop_closure(
    mut input: LoopClosureInput,
    cancellation: &SolveControl,
) -> Result<LoopClosurePreparation, Error> {
    validate(&input)?;
    let ids: Vec<_> = input
        .station_ids
        .iter()
        .map(|s| StationId(s.clone().into()))
        .collect();
    let mut shots = Vec::new();
    let mut instrument = BTreeMap::new();
    for (id, u, v, length, bearing, vertical, secondary, quality, offset) in &input.measurements {
        let bearing = if input.options.average_azimuths {
            secondary.map_or(*bearing, |other| {
                let difference = (other - bearing + 540.0).rem_euclid(360.0) - 180.0;
                (bearing + difference / 2.0).rem_euclid(360.0)
            })
        } else {
            *bearing
        };
        let angle = bearing.to_radians();
        let horizontal = ((length - vertical.abs()) * (length + vertical.abs())).sqrt();
        let value = [
            horizontal * angle.sin(),
            horizontal * angle.cos(),
            *vertical,
        ];
        if *offset {
            input
                .relationships
                .push((id.clone(), "hard_offset".into(), *u, *v, value));
        } else {
            shots.push((ObservationId(id.clone().into()), *u, *v, value));
            let q = if quality.is_finite() && *quality > 0.0 {
                *quality
            } else {
                1.0
            };
            let sig_l = input.options.sigma_length_m * q;
            let sig_a = input.options.sigma_azimuth_deg.to_radians() * q;
            instrument.insert(
                id.clone(),
                [
                    reading_square(angle.sin() * sig_l)
                        + reading_square(length * angle.cos() * sig_a),
                    reading_square(angle.cos() * sig_l)
                        + reading_square(length * angle.sin() * sig_a),
                ],
            );
        }
    }
    let connected = groups(
        ids.len(),
        shots
            .iter()
            .map(|r| (r.1, r.2))
            .chain(input.relationships.iter().map(|r| (r.2, r.3))),
    );
    let mut fixed: BTreeSet<_> = input
        .controls
        .iter()
        .filter(|r| r.2.contains(&0))
        .map(|r| r.1)
        .collect();
    if input.options.anchor_floating {
        for group in &connected {
            if group.iter().all(|i| !fixed.contains(i)) {
                let first = *input
                    .source_order
                    .iter()
                    .find(|i| group.contains(i))
                    .unwrap();
                fixed.insert(first);
                input.controls.push((
                    format!("component:{first}"),
                    first,
                    vec![0, 1],
                    input.initial_positions[first][..2].to_vec(),
                ));
            }
        }
    }
    if input.options.fix_control_depths {
        for &station in &fixed {
            if !input
                .controls
                .iter()
                .any(|r| r.1 == station && r.2.contains(&2))
            {
                input.controls.push((
                    format!("depth:{station}"),
                    station,
                    vec![2],
                    vec![input.initial_positions[station][2]],
                ));
            }
        }
    }
    if input.options.method != "least_squares" {
        let mut roots: Vec<_> = (0..ids.len()).collect();
        for (_, kind, u, v, _) in &input.relationships {
            if kind != "equivalent" {
                return Err(invalid(
                    "graph loop closure requires measured VIRTUAL shots",
                ));
            }
            let (u, v) = (root(&mut roots, *u), root(&mut roots, *v));
            roots[v] = u;
        }
        let roots: Vec<_> = (0..ids.len()).map(|i| root(&mut roots, i)).collect();
        let mut positions: Vec<_> = roots.iter().map(|&i| input.initial_positions[i]).collect();
        let mut values = BTreeMap::<(usize, usize), f64>::new();
        for (_, station, axes, coordinates) in &input.controls {
            let station = roots[*station];
            for (&a, &value) in axes.iter().zip(coordinates) {
                if values
                    .insert((station, a), value)
                    .is_some_and(|old| (old - value).abs() > 1e-8)
                {
                    return Err(invalid("conflicting controls in one station group"));
                }
                positions[station][a] = value;
            }
        }
        let anchors: BTreeSet<_> = values
            .keys()
            .filter(|(_, a)| *a == 0)
            .map(|(i, _)| *i)
            .collect();
        let dimension = if input.options.mode == "horizontal" {
            2
        } else {
            3
        };
        if values.keys().any(|(i, _)| !anchors.contains(i))
            || anchors
                .iter()
                .any(|i| (0..dimension).any(|a| !values.contains_key(&(*i, a))))
        {
            return Err(invalid(
                "graph anchors require complete coordinate controls",
            ));
        }
        for (i, position) in positions.iter_mut().enumerate() {
            if input.options.mode == "horizontal" {
                position[2] = input.initial_positions[i][2];
            }
        }
        return Ok(LoopClosurePreparation::Graph(GraphProblem {
            ids,
            shots,
            positions,
            roots,
            anchors,
            config: GraphConfig {
                reference: input.options.method == "reference_cg",
                horizontal: input.options.mode == "horizontal",
                max_iterations: input.options.max_iterations,
                absolute_tolerance: input.options.absolute_tolerance,
                relative_tolerance: input.options.relative_tolerance,
                weight_cap: input.options.weight_cap,
            },
        }));
    }
    if input.options.weighting == "length_azimuth_diagonal" {
        return prepare_axes(&input, &ids, &shots, &instrument, cancellation);
    }
    let observations = shots
        .iter()
        .map(|(id, u, v, value)| Observation::Displacement {
            id: id.clone(),
            from: ids[*u].clone(),
            to: ids[*v].clone(),
            value: *value,
            uncertainty: if input.options.weighting == "inverse_length" {
                let length = measured_norm(*value);
                Uncertainty::Relative {
                    weight: if length == 0.0 {
                        input.options.weight_cap
                    } else {
                        input.options.weight_cap.min(1.0 / length)
                    },
                }
            } else {
                Uncertainty::Isotropic {
                    sigma_m: input.options.sigma_m,
                }
            },
            reading_ids: vec![],
        })
        .collect();
    let mut p = problem(&input, &ids, observations);
    if input.options.preserve_centroids {
        let floating: Vec<_> = connected
            .into_iter()
            .filter(|g| g.iter().all(|i| !fixed.contains(i)))
            .collect();
        add_centroid_controls(&input, &ids, &mut p, &floating, None, cancellation)?;
    }
    Ok(LoopClosurePreparation::Standard(Box::new(prepare_problem(
        &p, true,
    )?)))
}
fn prepare_axes(
    input: &LoopClosureInput,
    ids: &[StationId],
    shots: &[(ObservationId, usize, usize, Point)],
    variances: &BTreeMap<String, [f64; 2]>,
    cancellation: &SolveControl,
) -> Result<LoopClosurePreparation, Error> {
    let mut equivalent: Vec<_> = (0..ids.len()).collect();
    for (_, kind, u, v, _) in &input.relationships {
        if kind == "equivalent" {
            let (u, v) = (root(&mut equivalent, *u), root(&mut equivalent, *v));
            equivalent[u] = v;
        }
    }
    let mut controlled = BTreeSet::new();
    for (id, i, _, _) in &input.controls {
        if id.starts_with("explicit:") && !controlled.insert(root(&mut equivalent, *i)) {
            return Err(invalid("Multiple geographic controls in one station group"));
        }
    }
    let mut axes = vec![];
    for a in 0..2 {
        let active: Vec<_> = shots
            .iter()
            .filter(|(id, u, v, _)| {
                variances[id.0.as_ref()][a] > 0.0
                    && root(&mut equivalent, *u) != root(&mut equivalent, *v)
            })
            .collect();
        let observations = active
            .iter()
            .map(|(id, u, v, value)| Observation::Displacement {
                id: id.clone(),
                from: ids[*u].clone(),
                to: ids[*v].clone(),
                value: [value[a], 0.0, 0.0],
                uncertainty: Uncertainty::Covariance {
                    matrix: [
                        [variances[id.0.as_ref()][a], 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                },
                reading_ids: vec![],
            })
            .collect();
        let mut p = problem(input, ids, observations);
        p.mode = LoopClosureMode::Horizontal {
            preserved_z: ids.iter().cloned().map(|id| (id, 0.0)).collect(),
        };
        p.controls = (0..ids.len())
            .map(|i| Control::Hard {
                id: format!("inactive-y:{i}"),
                station: ids[i].clone(),
                axes: vec![Axis::Y],
                value: vec![0.0],
            })
            .collect();
        let mut fixed = BTreeSet::new();
        for (id, station, components, values) in &input.controls {
            if let Some(index) = components.iter().position(|c| *c == a) {
                fixed.insert(*station);
                p.controls.push(Control::Hard {
                    id: id.clone(),
                    station: ids[*station].clone(),
                    axes: vec![Axis::X],
                    value: vec![values[index]],
                });
            }
        }
        p.relationships = input
            .relationships
            .iter()
            .map(|r| {
                let mut row = r.clone();
                row.4 = [r.4[a], 0.0, 0.0];
                relationship(&row, ids)
            })
            .collect();
        let floating: Vec<_> = groups(
            ids.len(),
            active
                .iter()
                .map(|r| (r.1, r.2))
                .chain(input.relationships.iter().map(|r| (r.2, r.3))),
        )
        .into_iter()
        .filter(|g| g.iter().all(|i| !fixed.contains(i)))
        .collect();
        add_centroid_controls(input, ids, &mut p, &floating, Some(a), cancellation)?;
        axes.push(prepare_problem(&p, true)?);
    }
    Ok(LoopClosurePreparation::Axes {
        axes,
        depths: ids
            .iter()
            .cloned()
            .zip(input.initial_positions.iter().map(|p| p[2]))
            .collect(),
    })
}

#[cfg(test)]
mod numerical_contract_tests {
    use super::*;
    #[test]
    fn measured_norm_is_exact_for_a_pythagorean_triple() {
        assert_eq!(measured_norm([3.0, 4.0, 0.0]), 5.0);
    }
    #[test]
    fn measured_norm_retains_zero() {
        assert_eq!(measured_norm([0.0; 3]), 0.0);
    }
    #[test]
    fn measured_norm_preserves_subnormal_measurements() {
        let tiny = f64::from_bits(1);
        assert_eq!(measured_norm([tiny, 0.0, 0.0]), tiny);
    }
    #[test]
    fn measured_norm_avoids_intermediate_square_overflow() {
        let value = measured_norm([1e200, 1e200, 0.0]);
        assert!(value.is_finite());
        assert!((value / 1e200 - 2.0_f64.sqrt()).abs() < 1e-15);
    }
    #[test]
    fn measured_norm_avoids_intermediate_square_underflow() {
        let value = measured_norm([1e-200, 1e-200, 0.0]);
        assert!(value > 0.0);
        assert!((value / 1e-200 - 2.0_f64.sqrt()).abs() < 1e-15);
    }
}
