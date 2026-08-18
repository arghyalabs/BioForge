//! Turing reaction-diffusion pattern formation (Activator-Inhibitor / Gierer-Meinhardt).

use serde::{Deserialize, Serialize};

use crate::diffusion::check_diffusion_cfl;
use crate::error::TissueError;
use crate::grid::Grid1D;

/// Gierer-Meinhardt Activator-Inhibitor Reaction-Diffusion Turing Patterning System (1972).
///
/// Spontaneously generates spatial stripes, spots, and periodic pigmentation patterns
/// via local autocatalytic activation coupled to long-range lateral inhibition ($D_v \gg D_u$):
///
/// $$\frac{\partial u}{\partial t} = D_u \nabla^2 u + a - b u + \frac{u^2}{v}$$
/// $$\frac{\partial v}{\partial t} = D_v \nabla^2 v + c u^2 - d v$$
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TuringSystem {
    /// Activator diffusion coefficient $D_u$ in $\mu\text{m}^2/\text{s}$.
    pub d_activator: f64,
    /// Inhibitor diffusion coefficient $D_v$ in $\mu\text{m}^2/\text{s}$ ($D_v \gg D_u$).
    pub d_inhibitor: f64,
    /// Basal activator production rate ($a$).
    pub a: f64,
    /// Linear activator decay rate ($b$).
    pub b: f64,
    /// Inhibitor production coefficient ($c$).
    pub c: f64,
    /// Linear inhibitor decay rate ($d$).
    pub d: f64,
}

impl Default for TuringSystem {
    fn default() -> Self {
        Self {
            d_activator: 0.01,
            d_inhibitor: 0.20,
            a: 0.01,
            b: 0.10,
            c: 0.10,
            d: 0.10,
        }
    }
}

impl TuringSystem {
    /// Perform a single numerical reaction-diffusion step for activator $u$ and inhibitor $v$.
    pub fn step(
        &self,
        u_grid: &mut Grid1D,
        v_grid: &mut Grid1D,
        dt_s: f64,
    ) -> Result<(), TissueError> {
        check_diffusion_cfl(self.d_inhibitor.max(self.d_activator), dt_s, u_grid.dx_um, 1)?;

        let lap_u = u_grid.compute_laplacian();
        let lap_v = v_grid.compute_laplacian();
        let n = u_grid.num_points;

        for i in 0..n {
            let u = u_grid.values[i].max(1e-6);
            let v = v_grid.values[i].max(1e-6);

            // Reaction kinetics
            let r_u = self.a - self.b * u + (u * u) / v;
            let r_v = self.c * u * u - self.d * v;

            let delta_u = dt_s * (self.d_activator * lap_u[i] + r_u);
            let delta_v = dt_s * (self.d_inhibitor * lap_v[i] + r_v);

            u_grid.values[i] = (u + delta_u).max(1e-6);
            v_grid.values[i] = (v + delta_v).max(1e-6);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_turing_system_step_advancement() {
        let turing = TuringSystem::default();
        let mut u = Grid1D::new(20, 1.0, 1.0).unwrap();
        let mut v = Grid1D::new(20, 1.0, 1.0).unwrap();

        // Step 10 times
        for _ in 0..10 {
            turing.step(&mut u, &mut v, 0.1).unwrap();
        }

        // Field values remain positive and finite
        for &val in &u.values {
            assert!(val.is_finite() && val > 0.0);
        }
        for &val in &v.values {
            assert!(val.is_finite() && val > 0.0);
        }
    }
}
