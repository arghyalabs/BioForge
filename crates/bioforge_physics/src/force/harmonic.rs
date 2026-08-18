//! Harmonic bond stretching force field ($U = \frac{1}{2} k_b (r - r_0)^2$).

use bioforge_state::SimulationState;

use super::ForceField;
use crate::error::PhysicsError;

/// Harmonic bond stretching force field based on Hooke's Law.
///
/// $$U_{\text{bond}}(r_{ij}) = \frac{1}{2} k_b (r_{ij} - r_0)^2 \quad [\text{kJ/mol}]$$
///
/// $$\vec{F}_i = -k_b (r_{ij} - r_0) \frac{\vec{r}_i - \vec{r}_j}{r_{ij}} \quad [(\text{kJ/mol})/\text{Å}]$$
/// $$\vec{F}_j = -\vec{F}_i$$
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HarmonicBondForce;

impl HarmonicBondForce {
    /// Create a new harmonic bond force evaluator.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl ForceField for HarmonicBondForce {
    fn compute_forces(&self, state: &mut SimulationState) -> Result<f64, PhysicsError> {
        let mut total_potential = 0.0;
        let box_size = state.box_size;

        for bond in &state.bonds {
            let i = bond.atom1;
            let j = bond.atom2;

            if i >= state.num_atoms || j >= state.num_atoms {
                continue;
            }

            let pi = state.positions[i];
            let pj = state.positions[j];

            let mut dx = pi[0] - pj[0];
            let mut dy = pi[1] - pj[1];
            let mut dz = pi[2] - pj[2];

            // Minimum image convention for periodic boundary conditions
            if let Some(box_l) = box_size {
                dx -= box_l[0] * (dx / box_l[0]).round();
                dy -= box_l[1] * (dy / box_l[1]).round();
                dz -= box_l[2] * (dz / box_l[2]).round();
            }

            let r_sq = dx * dx + dy * dy + dz * dz;
            if r_sq <= 1e-24 {
                return Err(PhysicsError::ZeroDistanceBond {
                    atom1: i,
                    atom2: j,
                });
            }

            let r = r_sq.sqrt();
            let dr = r - bond.r0;

            // U = 0.5 * kb * (r - r0)^2 in kJ/mol
            let u_bond = 0.5 * bond.kb * dr * dr;
            total_potential += u_bond;

            // F_mag = -kb * (r - r0)
            // F_vector_i = F_mag * (dr_vec / r)
            let f_factor = -bond.kb * dr / r;
            let fx = f_factor * dx;
            let fy = f_factor * dy;
            let fz = f_factor * dz;

            // Accumulate forces (Newton's third law: F_j = -F_i)
            state.forces[i][0] += fx;
            state.forces[i][1] += fy;
            state.forces[i][2] += fz;

            state.forces[j][0] -= fx;
            state.forces[j][1] -= fy;
            state.forces[j][2] -= fz;
        }

        Ok(total_potential)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::{BondOrder, Element};
    use bioforge_state::StateBond;

    #[test]
    fn test_harmonic_bond_force_and_energy() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.elements = vec![
            Element::from_symbol("C").unwrap(),
            Element::from_symbol("C").unwrap(),
        ];
        state.charges = vec![0.0, 0.0];
        state.velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];

        // Bond with r0 = 1.5 A, stretched to 1.6 A (dr = +0.1 A)
        // kb = 1000.0 (kJ/mol)/A^2
        state.positions = vec![[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]];
        state.bonds = vec![StateBond {
            atom1: 0,
            atom2: 1,
            r0: 1.5,
            kb: 1000.0,
            order: BondOrder::Single,
        }];

        let force_field = HarmonicBondForce::new();
        let u = force_field.compute_forces(&mut state).unwrap();

        // U = 0.5 * 1000.0 * 0.1^2 = 5.0 kJ/mol
        assert!((u - 5.0).abs() < 1e-6, "expected U=5.0, got U={}", u);

        // F acting on atom 0 (towards atom 1 in +x direction):
        // F_0 = -1000.0 * (1.6 - 1.5) * (0.0 - 1.6) / 1.6 = +100.0 (kJ/mol)/A
        assert!((state.forces[0][0] - 100.0).abs() < 1e-6);
        // F acting on atom 1 (towards atom 0 in -x direction):
        assert!((state.forces[1][0] - (-100.0)).abs() < 1e-6);
    }
}
