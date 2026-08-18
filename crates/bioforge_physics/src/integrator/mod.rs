//! Numerical integration algorithms for advancing simulation state through time.

mod verlet;

pub use verlet::{VelocityVerlet, FORCE_TO_ACCELERATION_FACTOR};

use bioforge_state::SimulationState;
use std::fmt::Debug;

use crate::error::PhysicsError;
use crate::force::ForceField;

/// Trait for numerical solvers and time integrators.
pub trait Integrator: Debug + Send + Sync {
    /// Advance `state` forward by time step `dt` (in $\text{ps}$) using `force_field`.
    ///
    /// # Return
    /// Returns the instantaneous potential energy $U$ in $\text{kJ/mol}$ at the new coordinates.
    fn step(
        &mut self,
        state: &mut SimulationState,
        dt: f64,
        force_field: &dyn ForceField,
    ) -> Result<f64, PhysicsError>;
}
