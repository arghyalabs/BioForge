//! # BioForge Measurement (`bioforge_measurement`)
//!
//! Scientific measurement and observer engine for the BioForge simulation platform.
//!
//! ## Scientific Architecture (Principle 2 & Principle 8)
//!
//! The measurement engine operates as a **pure read-only observer**. It extracts
//! physical metrics and observables from immutable references `&SimulationState`, guaranteeing
//! that observer directives in BioForge `.bio` source files cannot corrupt or alter the
//! underlying physical trajectory.
//!
//! ## Core Observables
//!
//! - [`DistanceObservable`]: Center-of-mass distance $d(A, B) = |\vec{R}_{\text{cm}}(A) - \vec{R}_{\text{cm}}(B)|$.
//! - [`RmsdObservable`]: Root Mean Square Deviation from reference coordinates.
//! - [`RadiusOfGyrationObservable`]: Spatial compactness $R_g$ of biomacromolecules.
//! - [`KineticEnergyObservable`], [`PotentialEnergyObservable`], [`TotalEnergyObservable`].
//! - [`TemperatureObservable`]: Instantaneous kinetic temperature.
//! - [`MeasurementEngine`]: Multi-observable time series recorder and CSV/JSON exporter.

#![deny(unsafe_code)]

pub mod engine;
pub mod error;
pub mod observable;

pub use engine::{MeasurementEngine, TimeSeries, TimeSeriesStatistics};
pub use error::MeasurementError;
pub use observable::{
    DistanceObservable, KineticEnergyObservable, Observable, PotentialEnergyObservable,
    RadiusOfGyrationObservable, RmsdObservable, TemperatureObservable, TotalEnergyObservable,
};

// ─── Trajectory Observation Integration Benchmark ───────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::{Element, Molecule};
    use bioforge_physics::{CompositeForceField, Integrator, VelocityVerlet};
    use bioforge_state::SimulationState;

    #[test]
    fn test_measurement_engine_observing_harmonic_trajectory() {
        use bioforge_biology::{Atom, Bond};

        let o = Element::from_symbol("O").unwrap();
        let h = Element::from_symbol("H").unwrap();
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(1, o, [0.0, 0.0, 0.0], "O"));
        mol.atoms.push(Atom::new(2, h, [1.0, 0.0, 0.0], "H1"));
        mol.atoms.push(Atom::new(3, h, [-0.333, 0.943, 0.0], "H2"));
        mol.bonds.push(Bond::single(1, 2));
        mol.bonds.push(Bond::single(1, 3));

        let mut state = SimulationState::from_molecule(&mol, None);
        state.thermalize(300.0, 42).unwrap();

        let mut engine = MeasurementEngine::new();
        engine.add_observable(DistanceObservable::pair("O-H1_dist", 0, 1));
        engine.add_observable(DistanceObservable::pair("O-H2_dist", 0, 2));
        engine.add_observable(RadiusOfGyrationObservable::new("water_Rg"));
        engine.add_observable(KineticEnergyObservable);
        engine.add_observable(PotentialEnergyObservable);
        engine.add_observable(TotalEnergyObservable);
        engine.add_observable(TemperatureObservable);

        let force_field = CompositeForceField::standard_molecular_mechanics(1.0, 10.0);
        let mut integrator = VelocityVerlet::new();

        // Simulate 200 steps and record measurements every 10 steps
        let dt = 0.0005; // 0.5 fs
        for step in 0..200 {
            let u = integrator.step(&mut state, dt, &force_field).unwrap();
            if step % 10 == 0 {
                engine.record_step(&state, u).unwrap();
            }
        }

        // 20 recorded data points
        assert_eq!(engine.time_series[0].data.len(), 20);

        let dist_stats = engine.statistics("O-H1_dist").unwrap();
        assert_eq!(dist_stats.count, 20);
        assert!(dist_stats.mean > 0.8 && dist_stats.mean < 1.2);

        let rg_stats = engine.statistics("water_Rg").unwrap();
        assert!(rg_stats.mean > 0.0);

        // Verify CSV export has all columns and valid numerical rows
        let csv = engine.export_csv();
        assert!(csv.contains("O-H1_dist"));
        assert!(csv.contains("water_Rg"));
        assert!(csv.contains("total_energy"));
    }
}
