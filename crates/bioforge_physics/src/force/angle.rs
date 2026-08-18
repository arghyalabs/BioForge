//! Harmonic angle bending potential ($U = \frac{1}{2} k_\theta (\theta - \theta_0)^2$).

use bioforge_state::SimulationState;

use super::ForceField;
use crate::error::PhysicsError;

/// Harmonic valence angle bending force field.
///
/// For atom triplet $(i, j, k)$ with central vertex $j$:
///
/// $$U_{\text{angle}}(\theta_{ijk}) = \frac{1}{2} k_\theta (\theta_{ijk} - \theta_0)^2 \quad [\text{kJ/mol}]$$
///
/// Analytical 3-body forces strictly preserve translational invariance:
/// $$\vec{F}_i + \vec{F}_j + \vec{F}_k = \vec{0}$$
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HarmonicAngleForce;

impl HarmonicAngleForce {
    /// Create a new harmonic angle force evaluator.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl ForceField for HarmonicAngleForce {
    fn compute_forces(&self, state: &mut SimulationState) -> Result<f64, PhysicsError> {
        let mut total_potential = 0.0;
        let box_size = state.box_size;

        for angle in &state.angles {
            let i = angle.atom1;
            let j = angle.atom2; // Central vertex
            let k = angle.atom3;

            if i >= state.num_atoms || j >= state.num_atoms || k >= state.num_atoms {
                continue;
            }

            let pi = state.positions[i];
            let pj = state.positions[j];
            let pk = state.positions[k];

            // Vector r_ji = r_i - r_j
            let mut d_ji = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
            // Vector r_jk = r_k - r_j
            let mut d_jk = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];

            // Minimum image convention for PBC
            if let Some(box_l) = box_size {
                d_ji[0] -= box_l[0] * (d_ji[0] / box_l[0]).round();
                d_ji[1] -= box_l[1] * (d_ji[1] / box_l[1]).round();
                d_ji[2] -= box_l[2] * (d_ji[2] / box_l[2]).round();

                d_jk[0] -= box_l[0] * (d_jk[0] / box_l[0]).round();
                d_jk[1] -= box_l[1] * (d_jk[1] / box_l[1]).round();
                d_jk[2] -= box_l[2] * (d_jk[2] / box_l[2]).round();
            }

            let r_ji_sq = d_ji[0] * d_ji[0] + d_ji[1] * d_ji[1] + d_ji[2] * d_ji[2];
            let r_jk_sq = d_jk[0] * d_jk[0] + d_jk[1] * d_jk[1] + d_jk[2] * d_jk[2];

            if r_ji_sq <= 1e-24 || r_jk_sq <= 1e-24 {
                continue;
            }

            let r_ji = r_ji_sq.sqrt();
            let r_jk = r_jk_sq.sqrt();

            let dot = d_ji[0] * d_jk[0] + d_ji[1] * d_jk[1] + d_ji[2] * d_jk[2];
            let cos_theta = (dot / (r_ji * r_jk)).clamp(-1.0, 1.0);
            let theta = cos_theta.acos();

            let d_theta = theta - angle.theta0;

            // U = 0.5 * k_theta * (theta - theta0)^2 in kJ/mol
            let u_angle = 0.5 * angle.k_theta * d_theta * d_theta;
            total_potential += u_angle;

            let sin_theta = (1.0 - cos_theta * cos_theta).max(1e-12).sqrt();
            let factor = angle.k_theta * d_theta / sin_theta;

            // Gradient components
            let inv_r_ji = 1.0 / r_ji;
            let inv_r_jk = 1.0 / r_jk;
            let inv_prod = inv_r_ji * inv_r_jk;

            let inv_r_ji_sq = inv_r_ji * inv_r_ji;
            let inv_r_jk_sq = inv_r_jk * inv_r_jk;

            // Force on atom i: F_i = factor * (r_jk / (r_ji * r_jk) - (cos_theta / r_ji^2) * r_ji)
            let f_i_x = factor * (d_jk[0] * inv_prod - cos_theta * inv_r_ji_sq * d_ji[0]);
            let f_i_y = factor * (d_jk[1] * inv_prod - cos_theta * inv_r_ji_sq * d_ji[1]);
            let f_i_z = factor * (d_jk[2] * inv_prod - cos_theta * inv_r_ji_sq * d_ji[2]);

            // Force on atom k: F_k = factor * (r_ji / (r_ji * r_jk) - (cos_theta / r_jk^2) * r_jk)
            let f_k_x = factor * (d_ji[0] * inv_prod - cos_theta * inv_r_jk_sq * d_jk[0]);
            let f_k_y = factor * (d_ji[1] * inv_prod - cos_theta * inv_r_jk_sq * d_jk[1]);
            let f_k_z = factor * (d_ji[2] * inv_prod - cos_theta * inv_r_jk_sq * d_jk[2]);

            // Force on central atom j: F_j = -(F_i + F_k)
            let f_j_x = -(f_i_x + f_k_x);
            let f_j_y = -(f_i_y + f_k_y);
            let f_j_z = -(f_i_z + f_k_z);

            state.forces[i][0] += f_i_x;
            state.forces[i][1] += f_i_y;
            state.forces[i][2] += f_i_z;

            state.forces[j][0] += f_j_x;
            state.forces[j][1] += f_j_y;
            state.forces[j][2] += f_j_z;

            state.forces[k][0] += f_k_x;
            state.forces[k][1] += f_k_y;
            state.forces[k][2] += f_k_z;
        }

        Ok(total_potential)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::Element;
    use bioforge_state::StateAngle;

    #[test]
    fn test_harmonic_angle_force_and_3body_sum() {
        let mut state = SimulationState::empty();
        state.num_atoms = 3;
        state.masses = vec![1.008, 15.999, 1.008]; // H-O-H
        state.charges = vec![0.417, -0.834, 0.417];
        let h = Element::from_symbol("H").unwrap();
        let o = Element::from_symbol("O").unwrap();
        state.elements = vec![h, o, h];
        state.forces = vec![[0.0, 0.0, 0.0]; 3];

        // Central oxygen at origin (0,0,0)
        // H1 at (1.0, 0.0, 0.0), H2 at (0.0, 1.0, 0.0) -> current angle is 90 degrees (pi/2)
        state.positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

        // Equilibrium angle is 104.5 degrees = 1.82387 rad
        let theta0 = 104.5_f64.to_radians();
        let k_theta = 500.0; // (kJ/mol)/rad^2

        state.angles = vec![StateAngle {
            atom1: 0,
            atom2: 1, // vertex
            atom3: 2,
            theta0,
            k_theta,
        }];

        let force_field = HarmonicAngleForce::new();
        let u = force_field.compute_forces(&mut state).unwrap();

        // Current angle is 90 deg = pi/2 = 1.570796 rad
        // d_theta = 1.570796 - 1.82387 = -0.25307 rad
        // U = 0.5 * 500.0 * (-0.25307)^2 = 16.012 kJ/mol
        assert!(u > 15.0 && u < 17.0, "got U={}", u);

        // Sum of all 3-body forces must be exactly zero (Newton's third law)
        let total_fx = state.forces[0][0] + state.forces[1][0] + state.forces[2][0];
        let total_fy = state.forces[0][1] + state.forces[1][1] + state.forces[2][1];
        let total_fz = state.forces[0][2] + state.forces[1][2] + state.forces[2][2];

        assert!(total_fx.abs() < 1e-10, "total fx non-zero: {}", total_fx);
        assert!(total_fy.abs() < 1e-10, "total fy non-zero: {}", total_fy);
        assert!(total_fz.abs() < 1e-10, "total fz non-zero: {}", total_fz);
    }
}
