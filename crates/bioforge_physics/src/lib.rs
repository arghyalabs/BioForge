//! # BioForge Physics (`bioforge_physics`)
//!
//! Numerical solvers, integrators, and molecular force fields for the BioForge simulation platform.
//!
//! ## Scientific Architecture (Principle 1 & Principle 9)
//!
//! Trajectories in BioForge are calculated exclusively by numerically integrating Newton's
//! equations of motion via symplectic solvers. The physics engine mutates [`SimulationState`],
//! while measurement systems and visualizers act strictly as read-only observers.
//!
//! ## Core Components
//!
//! - [`VelocityVerlet`]: $O(\Delta t^2)$ symplectic numerical integrator.
//! - [`HarmonicBondForce`]: Hookean potential and bond stretching forces.
//! - [`BerendsenThermostat`]: NVT temperature-coupling thermostat.
//! - [`CompositeForceField`]: Multi-term molecular force field aggregator.

#![deny(unsafe_code)]

pub mod error;
pub mod force;
pub mod integrator;
pub mod thermostat;

pub use error::PhysicsError;
pub use force::{CompositeForceField, ForceField, HarmonicBondForce};
pub use integrator::{Integrator, VelocityVerlet, FORCE_TO_ACCELERATION_FACTOR};
pub use thermostat::{BerendsenThermostat, Thermostat};

// ─── Rigorous Analytical Physics Benchmarks ─────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::{BondOrder, Element};
    use bioforge_state::{SimulationState, StateBond};

    /// Analytical Diatomic Harmonic Oscillator Benchmark
    ///
    /// Validates symplectic energy conservation over 2,000 integration steps (1.0 ps):
    /// - Carbon-Carbon diatomic bond: m1 = m2 = 12.011 Da
    /// - Equilibrium length: r0 = 1.50 A
    /// - Spring constant: kb = 2000.0 (kJ/mol)/A^2
    /// - Initial stretch: r(0) = 1.60 A (dr = +0.10 A), v(0) = [0, 0, 0]
    /// - Initial energy: E_0 = U_0 = 0.5 * 2000.0 * 0.1^2 = 10.000 kJ/mol
    /// - Integration time step: dt = 0.0005 ps (0.5 fs)
    #[test]
    fn test_harmonic_oscillator_energy_conservation() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.charges = vec![0.0, 0.0];
        let c = Element::from_symbol("C").unwrap();
        state.elements = vec![c, c];
        state.positions = vec![[0.0, 0.0, 0.0], [1.60, 0.0, 0.0]];
        state.velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.bonds = vec![StateBond {
            atom1: 0,
            atom2: 1,
            r0: 1.50,
            kb: 2000.0,
            order: BondOrder::Single,
        }];

        let force_field = HarmonicBondForce::new();
        let mut integrator = VelocityVerlet::new();

        let dt = 0.0005; // 0.5 fs
        let total_steps = 2000; // 1.0 ps total duration
        let initial_energy = 10.0; // 10.0 kJ/mol

        let mut max_energy_drift = 0.0;

        for _ in 0..total_steps {
            let u = integrator.step(&mut state, dt, &force_field).unwrap();
            let k = state.kinetic_energy();
            let total_e = k + u;

            let drift = (total_e - initial_energy).abs();
            if drift > max_energy_drift {
                max_energy_drift = drift;
            }
        }

        // Theoretical shadow Hamiltonian fluctuation amplitude for Velocity Verlet is:
        // Delta E / E0 = (omega * dt)^2 / 4 = (182.49 * 0.0005)^2 / 4 = 0.00208 (0.208%)
        // Expected Delta E = 10.0 * 0.00208 = 0.0208 kJ/mol.
        // Over 2,000 steps, symplectic Verlet strictly preserves this bound without secular drift.
        assert!(
            max_energy_drift < 0.025,
            "max energy drift = {} kJ/mol (threshold: 0.025 kJ/mol)",
            max_energy_drift
        );
    }

    /// Analytical Harmonic Oscillation Frequency & Period Benchmark
    ///
    /// Theoretical analytical period:
    /// - Reduced mass: mu = (m1 * m2) / (m1 + m2) = 12.011 / 2 = 6.0055 Da
    /// - Effective spring constant in internal units: k_eff = kb * 100 = 200,000 (Da/ps^2)
    /// - Angular frequency: omega = sqrt(k_eff / mu) = sqrt(200000 / 6.0055) = 182.489 ps^-1
    /// - Theoretical period: T = 2*pi / omega = 0.034430 ps = 34.43 fs
    #[test]
    fn test_harmonic_oscillator_period() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.charges = vec![0.0, 0.0];
        let c = Element::from_symbol("C").unwrap();
        state.elements = vec![c, c];
        state.positions = vec![[0.0, 0.0, 0.0], [1.60, 0.0, 0.0]];
        state.velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.bonds = vec![StateBond {
            atom1: 0,
            atom2: 1,
            r0: 1.50,
            kb: 2000.0,
            order: BondOrder::Single,
        }];

        let force_field = HarmonicBondForce::new();
        let mut integrator = VelocityVerlet::new();

        let dt = 0.0001; // 0.1 fs
        let total_steps = 1000;

        let theoretical_period = 2.0 * std::f64::consts::PI / (200_000.0 / 6.0055_f64).sqrt();

        // Track first full oscillation cycle (when position returns to maximum ~ 1.60 A)
        let mut prev_x = 1.60;
        let mut first_period_time = 0.0;
        let mut going_up = false;

        for _ in 0..total_steps {
            integrator.step(&mut state, dt, &force_field).unwrap();
            let current_x = state.positions[1][0] - state.positions[0][0];

            if current_x > prev_x {
                going_up = true;
            } else if going_up && current_x < prev_x && first_period_time == 0.0 {
                first_period_time = state.time;
            }
            prev_x = current_x;
        }

        // Numerical period must match theoretical period within 1%
        let period_error = (first_period_time - theoretical_period).abs() / theoretical_period;
        assert!(
            period_error < 0.01,
            "measured T={:.5} ps, theoretical T={:.5} ps (error: {:.2}%)",
            first_period_time,
            theoretical_period,
            period_error * 100.0
        );
    }

    /// Newton's Third Law & Linear Momentum Conservation Benchmark
    #[test]
    fn test_linear_momentum_conservation_during_simulation() {
        let mut state = SimulationState::empty();
        state.num_atoms = 3;
        state.masses = vec![12.011, 14.007, 15.999];
        state.charges = vec![0.0, 0.0, 0.0];
        state.elements = vec![
            Element::from_symbol("C").unwrap(),
            Element::from_symbol("N").unwrap(),
            Element::from_symbol("O").unwrap(),
        ];
        state.positions = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        state.velocities = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.bonds = vec![
            StateBond {
                atom1: 0,
                atom2: 1,
                r0: 1.40,
                kb: 1000.0,
                order: BondOrder::Single,
            },
            StateBond {
                atom1: 1,
                atom2: 2,
                r0: 1.40,
                kb: 1000.0,
                order: BondOrder::Single,
            },
        ];

        // Give initial velocity with zero CM drift
        state.thermalize(300.0, 99).unwrap();

        let force_field = HarmonicBondForce::new();
        let mut integrator = VelocityVerlet::new();

        for _ in 0..500 {
            integrator.step(&mut state, 0.001, &force_field).unwrap();
            let p = state.total_momentum();
            assert!(p[0].abs() < 1e-9, "px non-zero: {}", p[0]);
            assert!(p[1].abs() < 1e-9, "py non-zero: {}", p[1]);
            assert!(p[2].abs() < 1e-9, "pz non-zero: {}", p[2]);
        }
    }
}
