//! # BioForge Digital Twin (`bioforge_twin`)
//!
//! Interactive simulation engine, multi-fidelity reporting, and programmable biological world for BioForge.
//!
//! ## Grand Synthesis (Phases 22–28)
//!
//! Unifies all 17 BioForge crates into a single programmable digital organism platform:
//! - **Phase 22**: Interactive simulation sessions with Play, Pause, Step, Rewind, and live perturbation steering.
//! - **Phases 23–26**: Multi-scale biological digital twins with real-time multi-channel telemetry streams.
//! - **Phase 27**: Multi-fidelity reporting documenting assumptions, approximations, and numerical error bounds.
//! - **Phase 28**: The Complete Programmable Biological World engine.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod controller;
pub mod error;
pub mod fidelity;
pub mod twin;
pub mod world;

pub use controller::{InteractiveSession, PlaybackState, SimulationCommand};
pub use error::TwinError;
pub use fidelity::{FidelityLevel, MultiFidelityReport, ScaleFidelitySpec};
pub use twin::{BiologicalDigitalTwin, ScaleLayerKind, ScaleLayerState, TelemetryFrame};
pub use world::BiologicalWorld;

// ─── Grand Phase 28 Programmable World Benchmark ───────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_cell::GeneExpression;
    use bioforge_electrophysiology::{HodgkinHuxleyModel, StimulusProtocol};
    use bioforge_tissue::morphogen::MorphogenGradient;
    use bioforge_tissue::Grid1D;

    /// Grand Phase 28 Multi-Scale Programmable World Benchmark:
    ///
    /// Executes a fully coupled multiscale biological world simultaneously advancing:
    /// - Membrane action potentials ($V_m$)
    /// - Intracellular mRNA and protein gene expression synthesis ($[p]$)
    /// - Spatial morphogen diffusion fields ($C(x)$)
    #[test]
    fn test_grand_phase_28_programmable_world_execution() {
        let mut world = BiologicalWorld::new("WholeOrganism-001", 0.0005); // 0.5 ms

        // Electrophysiology
        world.electrophysiology = Some(HodgkinHuxleyModel::default());
        world.stimulus = Some(StimulusProtocol::CurrentPulse {
            amplitude_uA_cm2: 15.0,
            start_ms: 1.0,
            duration_ms: 3.0,
        });

        // Central Dogma
        world.gene_expression = Some(GeneExpression::new("Insulin", 1.0, 0.001, 0.5, 0.0001));

        // Spatial Tissue Morphogenesis
        world.morphogen = Some(MorphogenGradient::new("SonicHedgehog", 20.0, 0.002, 100.0));
        world.tissue_grid = Some(Grid1D::new(51, 5.0, 0.0).unwrap());

        // Simulate 20 ms of whole-system multiscale dynamics (40 steps)
        let telemetry_log = world.run_duration(0.020).unwrap();

        assert_eq!(telemetry_log.len(), 40);
        assert!((world.current_time_s - 0.020).abs() < 1e-6);

        // Verify multi-scale state integrity
        let final_frame = telemetry_log.last().unwrap();
        assert_eq!(final_frame.layer_states.len(), 3);

        // mRNA was synthesized
        assert!(world.mrna_conc_nM > 0.0);
        // Protein was translated
        assert!(world.protein_conc_nM > 0.0);
        // Morphogen boundary maintained
        assert_eq!(world.tissue_grid.as_ref().unwrap().values[0], 100.0);
    }
}
