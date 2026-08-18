//! # BioForge Multiscale (`bioforge_multiscale`)
//!
//! Scale bridges, thermodynamic transformations, and sub-cycling time coordinators for BioForge.
//!
//! ## Scientific Architecture (Principle 12 & Principle 13)
//!
//! Implements explicit, bidirectional scale bridge interfaces connecting biological layers:
//! - **Molecular $\longleftrightarrow$ Reaction**: Eyring transition state theory ($k_{\text{cat}}$) & binding free energy ($\Delta G^\circ \to K_d$).
//! - **Reaction $\longleftrightarrow$ Electrophysiology**: ATP metabolic consumption to active ion pumps ($Na^+/K^+$-ATPase) & ligand-gated channel conductances.
//! - **Cell $\longleftrightarrow$ Tissue**: Intracellular gene translation to spatial tissue morphogen secretion ($S(\vec{x})$) and receptor feedback.
//! - **Multiscale Time Coordinator**: Multi-rate nested sub-cycling across disparate biological timescales.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod bridges;
pub mod coordinator;
pub mod error;

pub use bridges::{
    binding_kinetics_from_delta_g, dissociation_constant_from_delta_g,
    eyring_catalytic_rate_constant, AtpIonPumpBridge, LigandGatedChannelBridge,
    MorphogenEmissionBridge, MorphogenReceptorBridge,
};
pub use coordinator::MultiscaleCoordinator;
pub use error::MultiscaleError;

// ─── Multiscale Coupled Feedback Integration Benchmark ─────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full Multiscale Reaction-Electrophysiology Coupling Benchmark:
    ///
    /// An enzymatic catalytic reaction generates ligand $L$ over time ($S \to L$),
    /// which binds to a ligand-gated channel and progressively activates an inward depolarizing current:
    ///
    /// $$\text{Metabolic Reaction} \xrightarrow{[L](t)} \text{Scale Bridge} \xrightarrow{g(t)} \text{Electrophysiology Membrane Current}$$
    #[test]
    fn test_coupled_reaction_to_electrophysiology_scale_bridge() {
        let channel_bridge = LigandGatedChannelBridge {
            g_bar: 5.0,     // 5 mS/cm^2
            e_rev: 0.0,     // 0 mV
            kd_molar: 10.0e-6, // 10 uM
            hill_n: 1.0,
        };

        // At t = 0, [L] = 0 uM -> Conductance is 0 mS/cm^2, current is 0 uA/cm^2
        let g_0 = channel_bridge.conductance(0.0);
        let i_0 = channel_bridge.current(-65.0, 0.0);
        assert_eq!(g_0, 0.0);
        assert_eq!(i_0, 0.0);

        // At t = 1, enzymatic reaction produces [L] = 10 uM (Kd)
        let g_1 = channel_bridge.conductance(10.0e-6);
        let i_1 = channel_bridge.current(-65.0, 10.0e-6);
        assert_eq!(g_1, 2.5); // 0.5 * 5.0
        assert_eq!(i_1, 2.5 * -65.0); // -162.5 uA/cm^2 inward current

        // At t = 2, [L] saturates to 90 uM (9 * Kd) -> Conductance is 90% of g_bar (4.5 mS/cm^2)
        let g_2 = channel_bridge.conductance(90.0e-6);
        assert!((g_2 - 4.5).abs() < 1e-6);
    }
}
