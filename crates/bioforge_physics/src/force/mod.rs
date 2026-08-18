//! Force field definitions and potential energy evaluations.

mod harmonic;

pub use harmonic::HarmonicBondForce;

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
