//! Morphogen gradients, Source-Diffusion-Degradation (SDD) PDEs, and the Wolpert French Flag model.

use serde::{Deserialize, Serialize};

use crate::diffusion::step_diffusion_1d;
use crate::error::TissueError;
use crate::grid::{BoundaryCondition, Grid1D};

/// A morphogen signaling molecule governing spatial developmental pattern formation.
///
/// Models embryonic morphogen gradients (e.g. *Bicoid*, *Sonic Hedgehog*, *Dpp*, *Wnt*)
/// using the Source-Diffusion-Degradation (SDD) paradigm:
///
/// $$\frac{\partial C}{\partial t} = D \nabla^2 C - k_{\text{deg}} C$$
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MorphogenGradient {
    /// Morphogen name (e.g. "Bicoid", "Shh", "Dpp").
    pub name: String,
    /// Diffusion coefficient $D$ in $\mu\text{m}^2/\text{s}$.
    pub diffusion_coeff_um2_s: f64,
    /// Linear degradation rate constant $k_{\text{deg}}$ in $\text{s}^{-1}$.
    pub degradation_rate_s: f64,
    /// Fixed boundary source concentration $C_0 = C(0)$ in nanomolar ($\text{nM}$).
    pub source_boundary_conc_nM: f64,
}

impl MorphogenGradient {
    /// Create a new morphogen gradient specification.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        d_um2_s: f64,
        k_deg_s: f64,
        source_conc_nM: f64,
    ) -> Self {
        Self {
            name: name.into(),
            diffusion_coeff_um2_s: d_um2_s.max(1e-6),
            degradation_rate_s: k_deg_s.max(1e-9),
            source_boundary_conc_nM: source_conc_nM.max(0.0),
        }
    }

    /// Characteristic morphogen decay length scale $\lambda = \sqrt{\frac{D}{k_{\text{deg}}}}$ in micrometers ($\mu\text{m}$).
    #[must_use]
    pub fn decay_length_um(&self) -> f64 {
        (self.diffusion_coeff_um2_s / self.degradation_rate_s).sqrt()
    }

    /// Exact analytical steady-state concentration at physical position $x$:
    ///
    /// $$C(x) = C_0 \exp\left( -\frac{x}{\lambda} \right)$$
    #[must_use]
    pub fn analytical_steady_state_at(&self, x_um: f64) -> f64 {
        let lambda = self.decay_length_um();
        self.source_boundary_conc_nM * (-x_um.max(0.0) / lambda).exp()
    }

    /// Numerically integrate the 1D SDD PDE to steady state over a physical domain of length $L$ ($\mu\text{m}$).
    pub fn simulate_to_steady_state_1d(
        &self,
        domain_length_um: f64,
        num_points: usize,
        total_time_s: f64,
    ) -> Result<Grid1D, TissueError> {
        let dx = domain_length_um / (num_points - 1) as f64;
        let mut grid = Grid1D::new(num_points, dx, 0.0)?;
        grid.boundary_left = BoundaryCondition::FixedDirichlet(self.source_boundary_conc_nM);
        grid.boundary_right = BoundaryCondition::ZeroFluxNeumann;

        // Choose stable dt: r = D*dt/dx^2 = 0.45 <= 0.5
        let dt = 0.45 * (dx * dx) / self.diffusion_coeff_um2_s;
        let mut t = 0.0;

        while t < total_time_s {
            let step_dt = dt.min(total_time_s - t);
            step_diffusion_1d(&mut grid, self.diffusion_coeff_um2_s, self.degradation_rate_s, step_dt)?;
            // Keep fixed source at x=0
            grid.values[0] = self.source_boundary_conc_nM;
            t += step_dt;
        }

        Ok(grid)
    }
}

/// The French Flag Model of Positional Information (Lewis Wolpert, J. Theor. Biol. 1969).
///
/// Maps continuous morphogen concentration gradients into discrete cellular developmental cell fates.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WolpertFrenchFlag {
    /// Upper concentration threshold defining the Blue (anterior) fate boundary in $\text{nM}$.
    pub high_threshold_nM: f64,
    /// Lower concentration threshold defining the White (middle) / Red (posterior) fate boundary in $\text{nM}$.
    pub low_threshold_nM: f64,
}

impl Default for WolpertFrenchFlag {
    fn default() -> Self {
        Self {
            high_threshold_nM: 25.0,
            low_threshold_nM: 10.0,
        }
    }
}

impl WolpertFrenchFlag {
    /// Determine the discrete cellular developmental fate from local morphogen concentration.
    #[must_use]
    pub fn determine_fate(&self, conc_nM: f64) -> &'static str {
        if conc_nM >= self.high_threshold_nM {
            "Blue (Anterior / High-Threshold Fate)"
        } else if conc_nM >= self.low_threshold_nM {
            "White (Medial / Mid-Threshold Fate)"
        } else {
            "Red (Posterior / Low-Threshold Fate)"
        }
    }

    /// Calculate analytical spatial boundary positions $(x_{\text{high}}, x_{\text{low}})$ in $\mu\text{m}$.
    #[must_use]
    pub fn analytical_boundary_positions(&self, morphogen: &MorphogenGradient) -> (f64, f64) {
        let lambda = morphogen.decay_length_um();
        let c0 = morphogen.source_boundary_conc_nM;

        let x_high = if self.high_threshold_nM < c0 {
            lambda * (c0 / self.high_threshold_nM).ln()
        } else {
            0.0
        };

        let x_low = if self.low_threshold_nM < c0 {
            lambda * (c0 / self.low_threshold_nM).ln()
        } else {
            0.0
        };

        (x_high, x_low)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sdd_morphogen_analytical_decay_length() {
        // D = 10 um^2/s, k_deg = 0.001 s^-1 => lambda = sqrt(10 / 0.001) = sqrt(10000) = 100.0 um
        let bicoid = MorphogenGradient::new("Bicoid", 10.0, 0.001, 50.0);
        assert_eq!(bicoid.decay_length_um(), 100.0);

        // At x = 0 um, C(0) = 50.0 nM
        assert_eq!(bicoid.analytical_steady_state_at(0.0), 50.0);

        // At x = 100 um (1 lambda), C(100) = 50.0 / e = 18.394 nM
        let c_100 = bicoid.analytical_steady_state_at(100.0);
        assert!((c_100 - 50.0 / std::f64::consts::E).abs() < 1e-6);
    }

    #[test]
    fn test_french_flag_fate_and_boundary_positions() {
        let bicoid = MorphogenGradient::new("Bicoid", 10.0, 0.001, 50.0);
        let flag = WolpertFrenchFlag {
            high_threshold_nM: 25.0, // C0 / 2 => x = lambda * ln(2) = 100 * 0.69315 = 69.3 um
            low_threshold_nM: 10.0,  // C0 / 5 => x = lambda * ln(5) = 100 * 1.60944 = 160.9 um
        };

        let (x_high, x_low) = flag.analytical_boundary_positions(&bicoid);
        assert!((x_high - 69.315).abs() < 0.1);
        assert!((x_low - 160.944).abs() < 0.1);

        assert!(flag.determine_fate(30.0).starts_with("Blue"));
        assert!(flag.determine_fate(15.0).starts_with("White"));
        assert!(flag.determine_fate(5.0).starts_with("Red"));
    }
}
