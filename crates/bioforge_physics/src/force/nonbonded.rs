//! Non-bonded molecular interactions (Lennard-Jones 12-6 and Coulomb Electrostatics).

use bioforge_state::SimulationState;

use super::ForceField;
use crate::error::PhysicsError;

/// Universal biophysical electrostatic conversion factor:
///
/// $$f_{\text{coulomb}} = \frac{e^2 N_A}{4\pi \varepsilon_0 \times 10^{-7}} = 1389.35456\text{ kJ}\cdot\text{Å}/(\text{mol}\cdot e^2)$$
pub const COULOMB_CONSTANT: f64 = 1389.354_56;

/// Non-bonded force field evaluating Lennard-Jones (12-6) and Coulomb interactions.
///
/// $$U_{\text{nonbonded}} = \sum_{i < j} \left( 4\varepsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6 \right] + \frac{f_{\text{coulomb}}}{\varepsilon_r} \frac{q_i q_j}{r_{ij}} \right)$$
#[derive(Debug, Clone, PartialEq)]
pub struct NonBondedForce {
    /// Non-bonded spherical interaction cutoff in Ångströms ($\text{Å}$).
    pub cutoff: f64,
    /// Relative dielectric permittivity constant $\varepsilon_r$ (default $1.0$ in vacuum).
    pub dielectric_constant: f64,
    /// Whether Lennard-Jones van der Waals interactions are active.
    pub enable_lj: bool,
    /// Whether Coulomb electrostatic interactions are active.
    pub enable_coulomb: bool,
}

impl Default for NonBondedForce {
    fn default() -> Self {
        Self {
            cutoff: 10.0,
            dielectric_constant: 1.0,
            enable_lj: true,
            enable_coulomb: true,
        }
    }
}

impl NonBondedForce {
    /// Create a new non-bonded force field with default parameters (cutoff $10.0\text{ \AA}$, $\varepsilon_r = 1.0$).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the spherical cutoff distance in Ångströms.
    #[must_use]
    pub fn with_cutoff(mut self, cutoff: f64) -> Self {
        self.cutoff = cutoff.max(1.0);
        self
    }

    /// Set relative dielectric permittivity constant.
    #[must_use]
    pub fn with_dielectric(mut self, dielectric: f64) -> Self {
        self.dielectric_constant = dielectric.max(0.1);
        self
    }
}

impl ForceField for NonBondedForce {
    fn compute_forces(&self, state: &mut SimulationState) -> Result<f64, PhysicsError> {
        let n = state.num_atoms;
        if n < 2 {
            return Ok(0.0);
        }

        let exclusions = state.get_exclusions();
        let cutoff_sq = self.cutoff * self.cutoff;
        let box_size = state.box_size;
        let coulomb_prefactor = COULOMB_CONSTANT / self.dielectric_constant;

        let mut total_potential = 0.0;

        for i in 0..n {
            let pi = state.positions[i];
            let qi = state.charges[i];
            let (sigma_i, eps_i) = if i < state.vdw_params.len() {
                (state.vdw_params[i][0], state.vdw_params[i][1])
            } else {
                (3.0, 0.5)
            };

            for j in (i + 1)..n {
                // Skip 1-2 (bonded) and 1-3 (angle) exclusions
                if exclusions.contains(&(i, j)) {
                    continue;
                }

                let pj = state.positions[j];
                let mut dx = pi[0] - pj[0];
                let mut dy = pi[1] - pj[1];
                let mut dz = pi[2] - pj[2];

                // Periodic boundary conditions minimum image convention
                if let Some(box_l) = box_size {
                    dx -= box_l[0] * (dx / box_l[0]).round();
                    dy -= box_l[1] * (dy / box_l[1]).round();
                    dz -= box_l[2] * (dz / box_l[2]).round();
                }

                let r_sq = dx * dx + dy * dy + dz * dz;
                if r_sq > cutoff_sq {
                    continue;
                }
                if r_sq <= 1e-24 {
                    return Err(PhysicsError::ZeroDistanceBond {
                        atom1: i,
                        atom2: j,
                    });
                }

                let r = r_sq.sqrt();
                let inv_r = 1.0 / r;
                let mut f_factor = 0.0;

                // --- 1. Lennard-Jones (12-6) Potential ---
                if self.enable_lj {
                    let (sigma_j, eps_j) = if j < state.vdw_params.len() {
                        (state.vdw_params[j][0], state.vdw_params[j][1])
                    } else {
                        (3.0, 0.5)
                    };

                    // Lorentz-Berthelot combining rules
                    let sigma_ij = 0.5 * (sigma_i + sigma_j);
                    let eps_ij = (eps_i * eps_j).sqrt();

                    if eps_ij > 0.0 {
                        let sr = sigma_ij * inv_r;
                        let sr2 = sr * sr;
                        let sr6 = sr2 * sr2 * sr2;
                        let sr12 = sr6 * sr6;

                        // U_LJ = 4 * eps * (sr12 - sr6)
                        let u_lj = 4.0 * eps_ij * (sr12 - sr6);
                        total_potential += u_lj;

                        // F_mag_LJ = (24 * eps / r) * (2*sr12 - sr6)
                        // F_factor = F_mag / r = (24 * eps / r^2) * (2*sr12 - sr6)
                        let f_lj_factor = (24.0 * eps_ij * inv_r * inv_r) * (2.0 * sr12 - sr6);
                        f_factor += f_lj_factor;
                    }
                }

                // --- 2. Coulomb Electrostatics ---
                if self.enable_coulomb {
                    let qj = state.charges[j];
                    let q_product = qi * qj;

                    if q_product.abs() > 1e-12 {
                        // U_coul = (f_coulomb / eps_r) * (qi * qj / r)
                        let u_coul = coulomb_prefactor * q_product * inv_r;
                        total_potential += u_coul;

                        // F_coul_factor = (f_coulomb / eps_r) * (qi * qj / r^3)
                        let f_coul_factor = coulomb_prefactor * q_product * inv_r * inv_r * inv_r;
                        f_factor += f_coul_factor;
                    }
                }

                // Apply forces on atom i and Newton's 3rd law on atom j
                let fx = f_factor * dx;
                let fy = f_factor * dy;
                let fz = f_factor * dz;

                state.forces[i][0] += fx;
                state.forces[i][1] += fy;
                state.forces[i][2] += fz;

                state.forces[j][0] -= fx;
                state.forces[j][1] -= fy;
                state.forces[j][2] -= fz;
            }
        }

        Ok(total_potential)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::Element;

    #[test]
    fn test_lennard_jones_energy_minimum_and_zero_force() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.charges = vec![0.0, 0.0];
        let c = Element::from_symbol("C").unwrap();
        state.elements = vec![c, c];

        let sigma = 3.4; // A
        let epsilon = 1.0; // kJ/mol
        state.vdw_params = vec![[sigma, epsilon], [sigma, epsilon]];

        // Theoretical potential energy minimum is at r_min = 2^(1/6) * sigma
        let r_min = 2.0_f64.powf(1.0 / 6.0) * sigma;
        state.positions = vec![[0.0, 0.0, 0.0], [r_min, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];

        let nonbonded = NonBondedForce::new().with_cutoff(15.0);
        let u = nonbonded.compute_forces(&mut state).unwrap();

        // At r_min, potential energy must equal exactly -epsilon = -1.0 kJ/mol
        assert!(
            (u - (-epsilon)).abs() < 1e-6,
            "expected U = -{}, got U = {}",
            epsilon,
            u
        );

        // At r_min, net force on both atoms must be exactly 0.0 (equilibrium)
        assert!(state.forces[0][0].abs() < 1e-6);
        assert!(state.forces[1][0].abs() < 1e-6);
    }

    #[test]
    fn test_coulomb_exact_analytical_potential_and_force() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![22.990, 35.450]; // Na+ and Cl-
        state.charges = vec![1.0, -1.0]; // +1e and -1e
        state.elements = vec![
            Element::from_symbol("Na").unwrap(),
            Element::from_symbol("Cl").unwrap(),
        ];
        state.vdw_params = vec![[0.0, 0.0], [0.0, 0.0]]; // purely electrostatic test

        // Separation r = 10.0 A
        state.positions = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];

        let nonbonded = NonBondedForce::new()
            .with_cutoff(20.0)
            .with_dielectric(1.0);
        let u = nonbonded.compute_forces(&mut state).unwrap();

        // U_theoretical = 1389.35456 * (1.0 * -1.0) / 10.0 = -138.935456 kJ/mol
        let expected_u = -COULOMB_CONSTANT / 10.0;
        assert!(
            (u - expected_u).abs() < 1e-6,
            "expected U={}, got U={}",
            expected_u,
            u
        );

        // Attractive force acting on atom 0 (+x direction towards atom 1):
        // F_x = 1389.35456 * (1.0 * -1.0) / 10.0^3 * (-10.0) = +13.8935456 (kJ/mol)/A
        let expected_f = COULOMB_CONSTANT / 100.0;
        assert!((state.forces[0][0] - expected_f).abs() < 1e-6);
        assert!((state.forces[1][0] - (-expected_f)).abs() < 1e-6);
    }
}
