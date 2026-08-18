//! Symplectic Velocity Verlet numerical integrator ($O(\Delta t^2)$).

use bioforge_state::SimulationState;

use super::Integrator;
use crate::error::PhysicsError;
use crate::force::ForceField;

/// Acceleration conversion constant:
///
/// $$1\text{ (kJ/mol)/Å} / (1\text{ Da}) = 100.0\text{ Å/ps}^2$$
pub const FORCE_TO_ACCELERATION_FACTOR: f64 = 100.0;

/// Symplectic Velocity Verlet numerical integrator.
///
/// Properties:
/// - $O(\Delta t^2)$ global accuracy
/// - Symplectic (preserves phase-space volume and constant-energy surface)
/// - Time-reversible
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct VelocityVerlet {
    initialized: bool,
}

impl VelocityVerlet {
    /// Create a new Velocity Verlet integrator.
    #[must_use]
    pub fn new() -> Self {
        Self { initialized: false }
    }
}

impl Integrator for VelocityVerlet {
    fn step(
        &mut self,
        state: &mut SimulationState,
        dt: f64,
        force_field: &dyn ForceField,
    ) -> Result<f64, PhysicsError> {
        if dt <= 0.0 || !dt.is_finite() {
            return Err(PhysicsError::InvalidTimeStep { dt });
        }
        let n = state.num_atoms;
        if n == 0 {
            return Ok(0.0);
        }

        // On the very first step, calculate initial forces if not already present
        if !self.initialized {
            state.zero_forces();
            let _ = force_field.compute_forces(state)?;
            self.initialized = true;
        }

        let half_dt = 0.5 * dt;
        let dt_sq_half = 0.5 * dt * dt;

        // Step 1: Update positions r(t + dt) and half-step velocities v(t + dt/2)
        for i in 0..n {
            let m = state.masses[i];
            if m <= 0.0 || !m.is_finite() {
                return Err(PhysicsError::NonPositiveMass {
                    atom_index: i,
                    mass: m,
                });
            }

            let inv_m = FORCE_TO_ACCELERATION_FACTOR / m;
            let ax = state.forces[i][0] * inv_m;
            let ay = state.forces[i][1] * inv_m;
            let az = state.forces[i][2] * inv_m;

            // r(t + dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
            state.positions[i][0] += state.velocities[i][0] * dt + ax * dt_sq_half;
            state.positions[i][1] += state.velocities[i][1] * dt + ay * dt_sq_half;
            state.positions[i][2] += state.velocities[i][2] * dt + az * dt_sq_half;

            // v(t + dt/2) = v(t) + 0.5*a(t)*dt
            state.velocities[i][0] += ax * half_dt;
            state.velocities[i][1] += ay * half_dt;
            state.velocities[i][2] += az * half_dt;

            // Check for numerical explosion
            let p = state.positions[i];
            if !p[0].is_finite() || !p[1].is_finite() || !p[2].is_finite() {
                return Err(PhysicsError::NumericalExplosion {
                    step: state.step,
                    atom_index: i,
                    x: p[0],
                    y: p[1],
                    z: p[2],
                });
            }
        }

        // Step 2: Compute forces F(t + dt) at new positions r(t + dt)
        state.zero_forces();
        let potential_energy = force_field.compute_forces(state)?;

        // Step 3: Complete velocity update v(t + dt) = v(t + dt/2) + 0.5*a(t + dt)*dt
        for i in 0..n {
            let m = state.masses[i];
            let inv_m = FORCE_TO_ACCELERATION_FACTOR / m;
            let ax_new = state.forces[i][0] * inv_m;
            let ay_new = state.forces[i][1] * inv_m;
            let az_new = state.forces[i][2] * inv_m;

            state.velocities[i][0] += ax_new * half_dt;
            state.velocities[i][1] += ay_new * half_dt;
            state.velocities[i][2] += az_new * half_dt;
        }

        // Advance simulation clock
        state.time += dt;
        state.step += 1;

        Ok(potential_energy)
    }
}
