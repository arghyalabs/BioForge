//! Force field definitions and potential energy evaluations.

mod angle;
mod harmonic;
mod nonbonded;

pub use angle::HarmonicAngleForce;
pub use harmonic::HarmonicBondForce;
pub use nonbonded::{NonBondedForce, COULOMB_CONSTANT};

use bioforge_state::SimulationState;
use std::fmt::Debug;

use crate::error::PhysicsError;

/// Trait for molecular force field potential evaluations.
pub trait ForceField: Debug + Send + Sync {
    /// Evaluate forces acting on all atoms in `state` and return the total potential energy in $\text{kJ/mol}$.
    ///
    /// # Contract
    /// - Implementations **MUST accumulate (add)** forces into `state.forces` (i.e. `state.forces[i] += F_i`),
    ///   allowing multiple force fields to be composed together.
    /// - Returns total potential energy $U$ in $\text{kJ/mol}$.
    fn compute_forces(&self, state: &mut SimulationState) -> Result<f64, PhysicsError>;
}

/// A composite collection of multiple force fields (e.g., bonds + angles + electrostatics).
#[derive(Debug, Default)]
pub struct CompositeForceField {
    /// Individual force field components.
    pub components: Vec<Box<dyn ForceField>>,
}

impl CompositeForceField {
    /// Create an empty composite force field.
    #[must_use]
    pub fn new() -> Self {
        Self {
            components: Vec::new(),
        }
    }

    /// Add a force field component.
    pub fn add<F: ForceField + 'static>(&mut self, force: F) {
        self.components.push(Box::new(force));
    }

    /// Construct a standard full Molecular Mechanics force field:
    /// Harmonic Bonds + Harmonic Angles + Lennard-Jones (12-6) + Coulomb Electrostatics.
    #[must_use]
    pub fn standard_molecular_mechanics(dielectric_constant: f64, cutoff_angstrom: f64) -> Self {
        let mut ff = Self::new();
        ff.add(HarmonicBondForce::new());
        ff.add(HarmonicAngleForce::new());
        ff.add(
            NonBondedForce::new()
                .with_dielectric(dielectric_constant)
                .with_cutoff(cutoff_angstrom),
        );
        ff
    }
}

impl ForceField for CompositeForceField {
    fn compute_forces(&self, state: &mut SimulationState) -> Result<f64, PhysicsError> {
        let mut total_potential_energy = 0.0;
        for component in &self.components {
            total_potential_energy += component.compute_forces(state)?;
        }
        Ok(total_potential_energy)
    }
}
