//! Thermostat implementations for NVT ensemble temperature regulation.

mod berendsen;

pub use berendsen::BerendsenThermostat;

use bioforge_state::SimulationState;
use std::fmt::Debug;

/// Trait for temperature-coupling thermostats in NVT simulations.
pub trait Thermostat: Debug + Send + Sync {
    /// Apply temperature coupling to `state` over time step `dt` (in $\text{ps}$).
    fn apply(&mut self, state: &mut SimulationState, dt: f64);
}
