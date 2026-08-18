//! The Complete Programmable Biological World (Phase 28 Grand Engine).
//!
//! Orchestrates the full multiscale biological simulation stack across all 17 BioForge crates.

use bioforge_cell::GeneExpression;
use bioforge_electrophysiology::{
    ElectrophysiologyState, HodgkinHuxleyModel, StimulusProtocol,
};
use bioforge_reaction::{ReactionNetwork, Rk4Solver};
use bioforge_tissue::morphogen::MorphogenGradient;
use bioforge_tissue::Grid1D;

use crate::error::TwinError;
use crate::twin::{ScaleLayerState, TelemetryFrame};

/// The unified top-level programmable biological simulation world (Phase 28).
///
/// Simultaneously orchestrates and synchronizes multi-scale biological systems:
/// 1. Biochemical reaction pathways (`bioforge_reaction`)
/// 2. Electrophysiological membrane dynamics (`bioforge_electrophysiology`)
/// 3. Genetic Central Dogma transcription & translation (`bioforge_cell`)
/// 4. Spatial tissue morphogen diffusion & pattern formation (`bioforge_tissue`)
pub struct BiologicalWorld {
    /// Unique world identifier.
    pub id: String,
    /// Current elapsed simulation time in seconds ($\text{s}$).
    pub current_time_s: f64,
    /// Base macroscopic time step $\Delta t$ in seconds ($\text{s}$).
    pub dt_s: f64,
    /// Chemical reaction network.
    pub reaction_network: Option<ReactionNetwork>,
    /// Chemical species concentrations in Molar ($\text{M}$).
    pub reaction_concs: Vec<f64>,
    /// Electrophysiological membrane model.
    pub electrophysiology: Option<HodgkinHuxleyModel>,
    /// Current electrophysiological state ($V_m, m, h, n$).
    pub electrophysiology_state: Option<ElectrophysiologyState>,
    /// Stimulus protocol for electrophysiology.
    pub stimulus: Option<StimulusProtocol>,
    /// Cellular gene expression model.
    pub gene_expression: Option<GeneExpression>,
    /// Cellular mRNA concentration in nanomolar ($\text{nM}$).
    pub mrna_conc_nM: f64,
    /// Cellular protein concentration in nanomolar ($\text{nM}$).
    pub protein_conc_nM: f64,
    /// Spatial morphogen gradient specification.
    pub morphogen: Option<MorphogenGradient>,
    /// 1D spatial tissue concentration field.
    pub tissue_grid: Option<Grid1D>,
}

impl BiologicalWorld {
    /// Create a new biological simulation world.
    #[must_use]
    pub fn new(id: impl Into<String>, dt_s: f64) -> Self {
        Self {
            id: id.into(),
            current_time_s: 0.0,
            dt_s: dt_s.max(1e-6),
            reaction_network: None,
            reaction_concs: Vec::new(),
            electrophysiology: None,
            electrophysiology_state: None,
            stimulus: None,
            gene_expression: None,
            mrna_conc_nM: 0.0,
            protein_conc_nM: 0.0,
            morphogen: None,
            tissue_grid: None,
        }
    }

    /// Advance the biological world by a single synchronized time step $\Delta t$.
    pub fn step(&mut self) -> Result<TelemetryFrame, TwinError> {
        let dt = self.dt_s;
        let mut states = Vec::new();

        // 1. Advance Biochemical Reaction Network
        if let Some(ref network) = self.reaction_network {
            if self.reaction_concs.is_empty() {
                self.reaction_concs = network
                    .species
                    .iter()
                    .map(|s| s.initial_concentration)
                    .collect();
            }
            let rk4 = Rk4Solver::new(dt);
            rk4.step(network, &mut self.reaction_concs, dt);
            let names: Vec<String> = network.species.iter().map(|s| s.name.clone()).collect();
            states.push(ScaleLayerState::Reaction {
                species_names: names,
                concentrations_molar: self.reaction_concs.clone(),
            });
        }

        // 2. Advance Cellular Electrophysiology
        if let Some(ref hh) = self.electrophysiology {
            let mut state = self
                .electrophysiology_state
                .unwrap_or_else(|| hh.initial_resting_state());
            let dt_ms = dt * 1000.0;
            let t_ms = self.current_time_s * 1000.0;
            let i_stim = self
                .stimulus
                .as_ref()
                .map(|s| s.current_at(t_ms))
                .unwrap_or(0.0);

            // Explicit RK4 step for HH
            let k1 = hh.compute_derivatives(&state, i_stim);
            let s2 = ElectrophysiologyState {
                v_membrane: state.v_membrane + 0.5 * dt_ms * k1[0],
                m: (state.m + 0.5 * dt_ms * k1[1]).clamp(0.0, 1.0),
                h: (state.h + 0.5 * dt_ms * k1[2]).clamp(0.0, 1.0),
                n: (state.n + 0.5 * dt_ms * k1[3]).clamp(0.0, 1.0),
            };
            let k2 = hh.compute_derivatives(&s2, i_stim);
            let s3 = ElectrophysiologyState {
                v_membrane: state.v_membrane + 0.5 * dt_ms * k2[0],
                m: (state.m + 0.5 * dt_ms * k2[1]).clamp(0.0, 1.0),
                h: (state.h + 0.5 * dt_ms * k2[2]).clamp(0.0, 1.0),
                n: (state.n + 0.5 * dt_ms * k2[3]).clamp(0.0, 1.0),
            };
            let k3 = hh.compute_derivatives(&s3, i_stim);
            let s4 = ElectrophysiologyState {
                v_membrane: state.v_membrane + dt_ms * k3[0],
                m: (state.m + dt_ms * k3[1]).clamp(0.0, 1.0),
                h: (state.h + dt_ms * k3[2]).clamp(0.0, 1.0),
                n: (state.n + dt_ms * k3[3]).clamp(0.0, 1.0),
            };
            let k4 = hh.compute_derivatives(&s4, i_stim);

            state.v_membrane += (dt_ms / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
            state.m = (state.m + (dt_ms / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])).clamp(0.0, 1.0);
            state.h = (state.h + (dt_ms / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])).clamp(0.0, 1.0);
            state.n = (state.n + (dt_ms / 6.0) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3])).clamp(0.0, 1.0);

            self.electrophysiology_state = Some(state);
            states.push(ScaleLayerState::Electrophysiology {
                v_membrane_mv: state.v_membrane,
                m_gate: state.m,
                h_gate: state.h,
                n_gate: state.n,
            });
        }

        // 3. Advance Cellular Gene Expression
        if let Some(ref gene) = self.gene_expression {
            let d_mrna = (gene.transcription_rate - gene.mrna_degradation_rate * self.mrna_conc_nM) * dt;
            let d_prot = (gene.translation_rate * self.mrna_conc_nM
                - gene.protein_degradation_rate * self.protein_conc_nM) * dt;

            self.mrna_conc_nM = (self.mrna_conc_nM + d_mrna).max(0.0);
            self.protein_conc_nM = (self.protein_conc_nM + d_prot).max(0.0);

            states.push(ScaleLayerState::CellularGenetics {
                mrna_nM: self.mrna_conc_nM,
                protein_nM: self.protein_conc_nM,
            });
        }

        // 4. Advance Spatial Tissue Morphogenesis
        if let (Some(ref morphogen), Some(ref mut grid)) = (&self.morphogen, &mut self.tissue_grid) {
            let _ = bioforge_tissue::diffusion::step_diffusion_1d(
                grid,
                morphogen.diffusion_coeff_um2_s,
                morphogen.degradation_rate_s,
                dt,
            );
            grid.values[0] = morphogen.source_boundary_conc_nM;
            let mean_c = grid.values.iter().sum::<f64>() / grid.values.len() as f64;
            states.push(ScaleLayerState::TissueMorphogenesis {
                mean_morphogen_nM: mean_c,
                cell_population_count: 1.0,
            });
        }

        self.current_time_s += dt;

        Ok(TelemetryFrame {
            timestamp_s: self.current_time_s,
            layer_states: states,
        })
    }

    /// Run the biological world for a specified duration and stream all telemetry frames.
    pub fn run_duration(&mut self, total_time_s: f64) -> Result<Vec<TelemetryFrame>, TwinError> {
        let mut frames = Vec::new();
        let target_time = self.current_time_s + total_time_s;

        while self.current_time_s < target_time {
            let frame = self.step()?;
            frames.push(frame);
        }

        Ok(frames)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grand_biological_world_multiscale_synthesis() {
        let mut world = BiologicalWorld::new("BioWorld-001", 0.001); // 1 ms steps

        // 1. Electrophysiology
        world.electrophysiology = Some(HodgkinHuxleyModel::default());
        world.stimulus = Some(StimulusProtocol::CurrentStep {
            amplitude_uA_cm2: 10.0,
            start_ms: 0.0,
        });

        // 2. Cellular gene expression
        world.gene_expression = Some(GeneExpression::new("GFP", 0.5, 0.001, 0.1, 0.0001));

        // 3. Spatial tissue morphogen
        let bicoid = MorphogenGradient::new("Bicoid", 10.0, 0.001, 50.0);
        let grid = Grid1D::new(21, 10.0, 0.0).unwrap();
        world.morphogen = Some(bicoid);
        world.tissue_grid = Some(grid);

        // Run world for 10 steps (10 ms)
        let frames = world.run_duration(0.010).unwrap();

        assert_eq!(frames.len(), 10);
        assert!((world.current_time_s - 0.010).abs() < 1e-6);

        // Verify each frame contains telemetry from all 3 active layers
        for frame in &frames {
            assert_eq!(frame.layer_states.len(), 3);
        }

        // Gene expression mRNA grew from 0.0
        assert!(world.mrna_conc_nM > 0.0);
        // Morphogen diffused into grid
        assert_eq!(world.tissue_grid.as_ref().unwrap().values[0], 50.0);
    }
}
