//! Spatial diffusion solvers with von Neumann CFL numerical stability checking.

use crate::error::TissueError;
use crate::grid::{Grid1D, Grid2D};

/// Check the von Neumann / Courant-Friedrichs-Lewy (CFL) numerical stability condition for explicit diffusion.
pub fn check_diffusion_cfl(
    d_um2_s: f64,
    dt_s: f64,
    dx_um: f64,
    dim: usize,
) -> Result<(), TissueError> {
    let r = (d_um2_s * dt_s) / (dx_um * dx_um);
    let limit = if dim == 1 { 0.5 } else { 0.25 };

    if r > limit {
        Err(TissueError::CflStabilityViolation { value: r })
    } else {
        Ok(())
    }
}

/// Advance a 1D spatial concentration field by a single time step $\Delta t$:
///
/// $$\frac{\partial C}{\partial t} = D \nabla^2 C - k_{\text{deg}} C$$
pub fn step_diffusion_1d(
    grid: &mut Grid1D,
    d_um2_s: f64,
    k_deg_s: f64,
    dt_s: f64,
) -> Result<(), TissueError> {
    check_diffusion_cfl(d_um2_s, dt_s, grid.dx_um, 1)?;

    let laplacian = grid.compute_laplacian();
    for i in 0..grid.num_points {
        let delta = dt_s * (d_um2_s * laplacian[i] - k_deg_s * grid.values[i]);
        grid.values[i] = (grid.values[i] + delta).max(0.0);
    }
    Ok(())
}

/// Advance a 2D spatial concentration field by a single time step $\Delta t$.
pub fn step_diffusion_2d(
    grid: &mut Grid2D,
    d_um2_s: f64,
    k_deg_s: f64,
    dt_s: f64,
) -> Result<(), TissueError> {
    check_diffusion_cfl(d_um2_s, dt_s, grid.dx_um, 2)?;

    let laplacian = grid.compute_laplacian();
    for i in 0..grid.values.len() {
        let delta = dt_s * (d_um2_s * laplacian[i] - k_deg_s * grid.values[i]);
        grid.values[i] = (grid.values[i] + delta).max(0.0);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cfl_stability_check_rejects_unstable_parameters() {
        // D = 100, dt = 1.0, dx = 10.0 => r = 100 * 1 / 100 = 1.0 > 0.5 (Unstable!)
        let res = check_diffusion_cfl(100.0, 1.0, 10.0, 1);
        assert!(res.is_err());

        // D = 10, dt = 1.0, dx = 10.0 => r = 10 * 1 / 100 = 0.1 <= 0.5 (Stable!)
        let res_stable = check_diffusion_cfl(10.0, 1.0, 10.0, 1);
        assert!(res_stable.is_ok());
    }

    #[test]
    fn test_diffusion_step_smoothes_point_source() {
        let mut grid = Grid1D::new(11, 1.0, 0.0).unwrap();
        grid.values[5] = 100.0; // Central peak

        // Run 5 diffusion steps
        for _ in 0..5 {
            step_diffusion_1d(&mut grid, 0.1, 0.0, 0.5).unwrap();
        }

        // Peak should decrease and neighbors should gain concentration
        assert!(grid.values[5] < 100.0);
        assert!(grid.values[4] > 0.0);
        assert!(grid.values[6] > 0.0);
    }
}
