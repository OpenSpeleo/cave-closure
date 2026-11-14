//! Independent horizontal length/azimuth adjustment with measured-position warm starts.
//!
//! This compatibility model retains Stalagmight's scalar reduction and PCG rules.
//! Callers contract exact relationships, choose controls and supply world-frame
//! displacements. Group and observation order are significant for floating axes
//! and floating-point accumulation. No centroid gauges or variance floors are added.
mod solver;

use crate::{Error, SolveControl, Termination, finite, invalid};

/// One already-contracted station group. Fixed coordinates override the initial guess.
#[derive(Clone, Debug)]
pub struct HorizontalGroup {
    pub initial: [f64; 2],
    pub fixed: Option<[f64; 2]>,
}

/// An oriented observation. Displacement is interpreted by the caller, independently
/// of slope length, which is used only in uncertainty propagation.
#[derive(Clone, Debug)]
pub struct HorizontalObservation {
    pub from: usize,
    pub to: usize,
    pub displacement: [f64; 2],
    pub length_m: f64,
    pub azimuth_rad: f64,
    /// Multiplies both sigmas. Nonpositive or nonfinite quality falls back to one.
    pub quality: f64,
}

/// Ordered, owned problem with no consumer-specific station identities or depth model.
#[derive(Clone, Debug)]
pub struct HorizontalProblem {
    pub groups: Vec<HorizontalGroup>,
    pub observations: Vec<HorizontalObservation>,
    /// Finite sigmas; zero retains the historical omitted-axis behavior.
    pub sigma_length_m: f64,
    pub sigma_azimuth_deg: f64,
}

/// Kernel termination and the reconstructed original-system residual for one axis.
/// A finite historical iterate can be returned without convergence; inspect termination
/// when the caller requires a convergence guarantee.
#[derive(Clone, Copy, Debug)]
pub struct AxisDiagnostics {
    pub termination: Termination,
    pub iterations: usize,
    pub residual_infinity: f64,
}

impl Default for AxisDiagnostics {
    fn default() -> Self {
        Self {
            termination: Termination::Converged,
            iterations: 0,
            residual_infinity: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct HorizontalSolution {
    /// Coordinates in input group order, including exact supplied controls.
    pub coordinates: Vec<[f64; 2]>,
    pub axes: [AxisDiagnostics; 2],
}

struct System {
    diag_x: Vec<f64>,
    diag_y: Vec<f64>,
    atb_x: Vec<f64>,
    atb_y: Vec<f64>,
    adj: Vec<Vec<(usize, f64, f64)>>,
}

fn check_cancelled(control: &SolveControl) -> Result<(), Error> {
    if control.is_cancelled() {
        Err(Error::Incomplete(Termination::Cancelled))
    } else {
        Ok(())
    }
}

impl HorizontalProblem {
    /// Solve without modifying the request. Concurrent calls use independent workspaces.
    /// Finite incomplete iterates retain historical behavior; cancellation, invalid
    /// indexes and nonfinite reconstruction return errors rather than partial results.
    pub fn solve(&self, control: &SolveControl) -> Result<HorizontalSolution, Error> {
        check_cancelled(control)?;
        if !finite([self.sigma_length_m, self.sigma_azimuth_deg])
            || self
                .groups
                .iter()
                .any(|g| !finite(g.initial) || g.fixed.is_some_and(|p| !finite(p)))
        {
            return Err(invalid("nonfinite horizontal group or uncertainty"));
        }
        let mut unknown = Vec::with_capacity(self.groups.len());
        let mut x0 = [Vec::new(), Vec::new()];
        for group in &self.groups {
            unknown.push(if group.fixed.is_none() {
                let index = x0[0].len();
                x0[0].push(group.initial[0]);
                x0[1].push(group.initial[1]);
                Some(index)
            } else {
                None
            });
        }
        let n = x0[0].len();
        let mut system = System {
            diag_x: vec![0.0; n],
            diag_y: vec![0.0; n],
            atb_x: vec![0.0; n],
            atb_y: vec![0.0; n],
            adj: vec![Vec::new(); n],
        };
        for shot in &self.observations {
            check_cancelled(control)?;
            if shot.from >= self.groups.len()
                || shot.to >= self.groups.len()
                || !finite(shot.displacement)
                || !finite([shot.length_m, shot.azimuth_rad])
            {
                return Err(invalid("invalid horizontal observation"));
            }
            if shot.from == shot.to {
                continue;
            }
            let quality = if shot.quality.is_finite() && shot.quality > 0.0 {
                shot.quality
            } else {
                1.0
            };
            let sigma_len = self.sigma_length_m * quality;
            let sigma_az = self.sigma_azimuth_deg.to_radians() * quality;
            let az = shot.azimuth_rad;
            let len = shot.length_m;
            // Keep sqrt followed by squaring: removing that round trip changes
            // the source model's last bits and occasionally its stopping iteration.
            let sigma_dx =
                ((az.sin() * sigma_len).powi(2) + (len * az.cos() * sigma_az).powi(2)).sqrt();
            let sigma_dy =
                ((az.cos() * sigma_len).powi(2) + (len * az.sin() * sigma_az).powi(2)).sqrt();
            let wx = if sigma_dx > 0.0 {
                1.0 / (sigma_dx * sigma_dx)
            } else {
                0.0
            };
            let wy = if sigma_dy > 0.0 {
                1.0 / (sigma_dy * sigma_dy)
            } else {
                0.0
            };
            let i = unknown[shot.from];
            let j = unknown[shot.to];
            let rhs = |axis: usize| match (i, j) {
                (Some(_), Some(_)) => shot.displacement[axis],
                (Some(_), None) => {
                    shot.displacement[axis] - self.groups[shot.to].fixed.unwrap()[axis]
                }
                (None, Some(_)) => {
                    shot.displacement[axis] + self.groups[shot.from].fixed.unwrap()[axis]
                }
                (None, None) => 0.0,
            };
            for (w, b, diag, atb) in [
                (wx, rhs(0), &mut system.diag_x, &mut system.atb_x),
                (wy, rhs(1), &mut system.diag_y, &mut system.atb_y),
            ] {
                if w.is_finite() && w > 0.0 {
                    match (i, j) {
                        (Some(i), Some(j)) => {
                            diag[i] += w;
                            diag[j] += w;
                            atb[i] += -w * b;
                            atb[j] += w * b;
                        }
                        (Some(i), None) => {
                            diag[i] += w;
                            atb[i] += -w * b;
                        }
                        (None, Some(j)) => {
                            diag[j] += w;
                            atb[j] += w * b;
                        }
                        (None, None) => {}
                    }
                }
            }
            if let (Some(i), Some(j)) = (i, j) {
                for (from, to) in [(i, j), (j, i)] {
                    if let Some(edge) = system.adj[from].iter_mut().find(|e| e.0 == to) {
                        edge.1 -= wx;
                        edge.2 -= wy;
                    } else {
                        system.adj[from].push((to, -wx, -wy));
                    }
                }
            }
        }
        let mut solution = solver::solve_reduced_laplacian(&system, &x0[0], &x0[1], control)?;
        for (axis, diag, rhs) in [
            (0, &system.diag_x, &system.atb_x),
            (1, &system.diag_y, &system.atb_y),
        ] {
            let values = &solution.coordinates[axis];
            if !finite(values.iter().copied()) {
                return Err(Error::Numerical(
                    "nonfinite horizontal reconstruction".into(),
                ));
            }
            for i in 0..n {
                let mut ax = diag[i] * values[i];
                for &(j, wx, wy) in &system.adj[i] {
                    ax += if axis == 0 { wx } else { wy } * values[j];
                }
                let residual = (ax - rhs[i]).abs();
                if !residual.is_finite() {
                    return Err(Error::Numerical("nonfinite horizontal residual".into()));
                }
                solution.axes[axis].residual_infinity =
                    solution.axes[axis].residual_infinity.max(residual);
            }
        }
        check_cancelled(control)?;
        Ok(HorizontalSolution {
            coordinates: self
                .groups
                .iter()
                .zip(unknown)
                .map(|(group, index)| {
                    index.map_or_else(
                        || group.fixed.unwrap(),
                        |i| [solution.coordinates[0][i], solution.coordinates[1][i]],
                    )
                })
                .collect(),
            axes: solution.axes,
        })
    }
}

#[cfg(test)]
mod tests;
