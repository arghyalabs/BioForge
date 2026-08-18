//! # BioForge Tissue (`bioforge_tissue`)
//!
//! Tissue morphogenesis, spatial reaction-diffusion PDEs, and multicellular populations for BioForge.
//!
//! ## Scientific Architecture (Principle 3 & Principle 12)
//!
//! Models macroscopic multicellular systems across multiple spatial and developmental layers:
//! - **Spatial Grid & PDE Operators**: Finite-difference Laplacians in 1D and 2D with Neumann, Dirichlet, and Periodic boundaries.
//! - **Morphogen Gradients**: Source-Diffusion-Degradation (SDD) PDEs, exponential decay length ($\lambda = \sqrt{D/k_{\text{deg}}}$), and Wolpert's French Flag positional information.
//! - **Turing Pattern Formation**: Gierer-Meinhardt activator-inhibitor systems generating spontaneous pigmentation patterns.
//! - **Population Dynamics**: Verhulst logistic growth with carrying capacity and contact inhibition of proliferation.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod diffusion;
pub mod error;
pub mod grid;
pub mod morphogen;
pub mod population;
pub mod trajectory;
pub mod turing;

pub use diffusion::{check_diffusion_cfl, step_diffusion_1d, step_diffusion_2d};
pub use error::TissueError;
pub use grid::{BoundaryCondition, Grid1D, Grid2D};
pub use morphogen::{MorphogenGradient, WolpertFrenchFlag};
pub use population::CellPopulation;
pub use trajectory::SpatialTrajectory1D;
pub use turing::TuringSystem;

// ─── Morphogen Gradient Numerical Benchmark ────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full SDD Morphogen Gradient Numerical Benchmark:
    ///
    /// Simulates spatial diffusion & degradation from a localized anterior source ($x=0$)
    /// and verifies that the numerical PDE solution converges to the exact theoretical exponential gradient:
    ///
    /// $$C(x) = C_0 \exp\left( -\frac{x}{\lambda} \right), \quad \lambda = \sqrt{\frac{D}{k_{\text{deg}}}} = 100.0\,\mu\text{m}$$
    #[test]
    fn test_sdd_morphogen_gradient_numerical_convergence() {
        let bicoid = MorphogenGradient::new("Bicoid", 10.0, 0.001, 50.0);
        let domain_length_um = 400.0;
        let num_points = 41; // dx = 10.0 um

        // Simulate PDE for 5,000 seconds to reach steady state
        let grid = bicoid
            .simulate_to_steady_state_1d(domain_length_um, num_points, 5000.0)
            .unwrap();

        // Check numerical vs analytical values across interior points (x = 0 to 300 um)
        for i in 0..30 {
            let x_um = i as f64 * grid.dx_um;
            let numerical_c = grid.values[i];
            let analytical_c = bicoid.analytical_steady_state_at(x_um);

            let diff = (numerical_c - analytical_c).abs();
            // Max deviation < 1.0 nM (relative error < 2%)
            assert!(
                diff < 1.0,
                "at x = {:.1} um: expected C = {:.3} nM, got {:.3} nM (diff = {:.3} nM)",
                x_um,
                analytical_c,
                numerical_c,
                diff
            );
        }
    }
}
