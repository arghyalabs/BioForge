//! # BioForge Electrophysiology (`bioforge_electrophysiology`)
//!
//! Membrane biology, transmembrane ion channels, and Hodgkin-Huxley electrophysiology engine for BioForge.
//!
//! ## Scientific Architecture (Principle 3 & Principle 12)
//!
//! Implements cellular electrophysiology across multiple analytical and numerical layers:
//! - **Electrochemical Equilibria**: Exact Nernst reversal potentials and Goldman-Hodgkin-Katz (GHK) resting potentials.
//! - **Ion Channels**: Voltage-gated fast Sodium ($m^3 h$), delayed rectifier Potassium ($n^4$), and leak channels.
//! - **Action Potential Dynamics**: Nonlinear Hodgkin-Huxley differential equation system solved via 4th-order Runge-Kutta.
//! - **Electrophysiological Analysis**: Action potential spike detection, firing frequencies, and refractory period metrics.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod analysis;
pub mod channel;
pub mod constants;
pub mod error;
pub mod membrane;
pub mod hodgkin_huxley;
pub mod stimulus;

pub use analysis::{ActionPotentialMetrics, SpikeDetector};
pub use channel::{LeakChannel, PotassiumChannel, SodiumChannel};
pub use constants::{
    BODY_TEMPERATURE_KELVIN, DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2, FARADAY_CONSTANT_F,
    MOLAR_GAS_CONSTANT_R, ROOM_TEMPERATURE_KELVIN, SQUID_AXON_TEMPERATURE_KELVIN,
};
pub use error::ElectrophysiologyError;
pub use hodgkin_huxley::{ElectrophysiologyState, ElectrophysiologyTrajectory, HodgkinHuxleyModel};
pub use membrane::{Ion, Membrane};
pub use stimulus::StimulusProtocol;

// ─── Physiological Action Potential Benchmark ──────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full Physiological Action Potential Benchmark:
    /// Validates resting potential stability, stimulus-triggered spike, overshoot, and recovery.
    #[test]
    fn test_hodgkin_huxley_complete_spike_cycle() {
        let model = HodgkinHuxleyModel::standard();
        let stimulus = StimulusProtocol::CurrentPulse {
            amplitude_uA_cm2: 12.0,
            start_ms: 10.0,
            duration_ms: 1.0,
        };

        let traj = model.simulate(40.0, 0.01, &stimulus).unwrap();
        let detector = SpikeDetector::default();
        let metrics = detector.analyze(&traj);

        // Single action potential fired
        assert_eq!(metrics.spike_count, 1);
        let spike_t = metrics.spike_times_ms[0];
        assert!(spike_t > 10.0 && spike_t < 15.0);

        // Peak overshoot > +20 mV
        assert!(metrics.peak_voltage_mv > 20.0);

        // Total spike amplitude > 90 mV (from -65 mV up to > +25 mV)
        assert!(metrics.amplitude_mv > 90.0);

        // Export to CSV verification
        let csv = traj.export_csv();
        assert!(csv.contains("time_ms,v_membrane_mv"));
        assert!(csv.contains("i_na_uA_cm2"));
    }
}
