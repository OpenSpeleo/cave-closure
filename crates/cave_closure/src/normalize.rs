//! First-order instrument normalization. All angles are radians, bearings clockwise from north.
use crate::*;
use nalgebra::{Matrix3, Matrix3x4, SMatrix};

fn array(m: Matrix3<f64>) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| m[(i, j)]))
}
pub fn polar(
    length: f64,
    bearing: f64,
    inclination: f64,
    covariance: [[f64; 3]; 3],
) -> Result<(Point, [[f64; 3]; 3]), Error> {
    if !finite([length, bearing, inclination])
        || length < 0.
        || inclination.abs() > std::f64::consts::FRAC_PI_2
    {
        return Err(invalid("invalid polar reading"));
    }
    let (sb, cb) = bearing.sin_cos();
    let (si, ci) = inclination.sin_cos();
    let j = Matrix3::new(
        ci * sb,
        length * ci * cb,
        -length * si * sb,
        ci * cb,
        -length * ci * sb,
        -length * si * cb,
        si,
        0.,
        length * ci,
    );
    let c = Matrix3::from_fn(|i, j| covariance[i][j]);
    super::prepare::precision(
        nalgebra::DMatrix::from_fn(3, 3, |i, j| c[(i, j)]),
        "raw polar covariance",
    )?;
    let propagated = j * c * j.transpose();
    let c = (propagated + propagated.transpose()) * 0.5;
    super::prepare::precision(
        nalgebra::DMatrix::from_fn(3, 3, |i, j| c[(i, j)]),
        "propagated polar covariance; supply declared placement noise for singular geometry",
    )?;
    Ok(([length * ci * sb, length * ci * cb, length * si], array(c)))
}
/// Depths are positive down. Shared input reading identities must still be tracked by the caller.
pub fn underwater(
    length: f64,
    bearing: f64,
    depth_u: f64,
    depth_v: f64,
    covariance: [[f64; 4]; 4],
    min_horizontal_ratio: f64,
) -> Result<(Point, [[f64; 3]; 3]), Error> {
    if !finite([length, bearing, depth_u, depth_v, min_horizontal_ratio])
        || length < 0.
        || !(0. ..1.).contains(&min_horizontal_ratio)
    {
        return Err(invalid(
            "invalid underwater reading or horizontal reliability policy",
        ));
    }
    let t = depth_v - depth_u;
    if length < t.abs() {
        return Err(invalid("length is shorter than absolute depth difference"));
    }
    let h = ((length - t.abs()) * (length + t.abs())).sqrt();
    if h == 0. || h <= length * min_horizontal_ratio {
        return Err(invalid(
            "near-vertical underwater normalization is singular or outside caller reliability policy",
        ));
    }
    let (sb, cb) = bearing.sin_cos();
    let j = Matrix3x4::new(
        length / h * sb,
        h * cb,
        t / h * sb,
        -t / h * sb,
        length / h * cb,
        -h * sb,
        t / h * cb,
        -t / h * cb,
        0.,
        0.,
        1.,
        -1.,
    );
    let c = SMatrix::<f64, 4, 4>::from_fn(|i, j| covariance[i][j]);
    super::prepare::precision(
        nalgebra::DMatrix::from_fn(4, 4, |i, j| c[(i, j)]),
        "raw underwater covariance",
    )?;
    let propagated = j * c * j.transpose();
    let c = (propagated + propagated.transpose()) * 0.5;
    super::prepare::precision(
        nalgebra::DMatrix::from_fn(3, 3, |i, j| c[(i, j)]),
        "propagated underwater covariance",
    )?;
    Ok(([h * sb, h * cb, -t], array(c)))
}
