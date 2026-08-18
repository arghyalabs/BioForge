//! # BioForge State (`bioforge_state`)
//!
//! Simulation state and runtime representation for the BioForge multiscale simulation engine.
//!
//! ## Architectural Role (Principle 2)
//!
//! [`SimulationState`] is the **single source of truth** for all physical quantities during
//! simulation execution. The physical state owns coordinates, velocities, masses, and forces,
//! while observers (measurements, visualizers) have read-only access.
//!
//! ## Core Types
//!
//! - [`SimulationState`]: Complete physical snapshot of $N$ particles at time $t$.
//! - [`StateBond`]: Bond connectivity and harmonic force parameters.
//! - [`Trajectory`]: Ring buffer storing snapshots for analysis and export.
//! - [`TrajectoryFrame`]: Single historical time point.
//! - [`StateError`]: Errors relating to state construction and validation.

#![deny(unsafe_code)]

mod error;
mod state;
mod trajectory;

pub use error::StateError;
pub use state::{
    default_vdw_for_element, SimulationState, StateAngle, StateBond,
    DA_A2_PER_PS2_TO_KJ_PER_MOL, MOLAR_GAS_CONSTANT_R, THERMAL_VELOCITY_CONSTANT,
};
pub use trajectory::{Trajectory, TrajectoryFrame};
