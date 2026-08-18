//! Biological digital twin representations and multi-scale live telemetry feeds (Phases 23–26).

use serde::{Deserialize, Serialize};

/// Classification of biological scale layers comprising a digital twin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ScaleLayerKind {
    /// Atomic coordinates, forces, and molecular mechanics.
    MolecularMechanics,
    /// Chemical species concentrations and reaction kinetics.
    ReactionKinetics,
    /// Membrane voltages, ion channels, and action potentials.
    Electrophysiology,
    /// Genetic Central Dogma and gene regulatory networks.
    CellularGenetics,
    /// Spatial tissue morphogen fields and multicellular populations.
    TissueMorphogenesis,
}

/// Instantaneous snapshot of a specific scale layer in a biological digital twin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ScaleLayerState {
    /// Atomistic molecular state.
    Molecular {
        num_atoms: usize,
        potential_energy_kJ_mol: f64,
        temperature_k: f64,
    },
    /// Chemical reaction network state.
    Reaction {
        species_names: Vec<String>,
        concentrations_molar: Vec<f64>,
    },
    /// Electrophysiological membrane state.
    Electrophysiology {
        v_membrane_mv: f64,
        m_gate: f64,
        h_gate: f64,
        n_gate: f64,
    },
    /// Cellular gene expression state.
    CellularGenetics {
        mrna_nM: f64,
        protein_nM: f64,
    },
    /// Spatial tissue state.
    TissueMorphogenesis {
        mean_morphogen_nM: f64,
        cell_population_count: f64,
    },
}

/// A synchronized multi-scale telemetry data frame streamed from a digital twin.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct TelemetryFrame {
    /// Global physical timestamp in seconds ($\text{s}$).
    pub timestamp_s: f64,
    /// Snapshot states across all active biological scale layers.
    pub layer_states: Vec<ScaleLayerState>,
}

/// A complete, multi-scale computational biological digital twin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BiologicalDigitalTwin {
    /// Unique digital twin identifier.
    pub id: String,
    /// Human-readable entity name (e.g. "Pancreatic_Beta_Cell_Twin", "Neuron_Soma_Twin").
    pub name: String,
    /// Active biological scale layers.
    pub active_layers: Vec<ScaleLayerKind>,
    /// Live telemetry streaming history.
    pub history: Vec<TelemetryFrame>,
}

impl BiologicalDigitalTwin {
    /// Create a new biological digital twin.
    #[must_use]
    pub fn new(id: impl Into<String>, name: impl Into<String>, layers: Vec<ScaleLayerKind>) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            active_layers: layers,
            history: Vec::new(),
        }
    }

    /// Record a live multi-scale telemetry frame.
    pub fn record_telemetry(&mut self, frame: TelemetryFrame) {
        self.history.push(frame);
    }

    /// Retrieve the most recent telemetry frame.
    #[must_use]
    pub fn latest_telemetry(&self) -> Option<&TelemetryFrame> {
        self.history.last()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_biological_digital_twin_telemetry_streaming() {
        let layers = vec![
            ScaleLayerKind::Electrophysiology,
            ScaleLayerKind::CellularGenetics,
        ];
        let mut twin = BiologicalDigitalTwin::new("TWIN-001", "NeuronDigitalTwin", layers);

        let frame = TelemetryFrame {
            timestamp_s: 0.05,
            layer_states: vec![
                ScaleLayerState::Electrophysiology {
                    v_membrane_mv: -65.0,
                    m_gate: 0.05,
                    h_gate: 0.60,
                    n_gate: 0.32,
                },
                ScaleLayerState::CellularGenetics {
                    mrna_nM: 50.0,
                    protein_nM: 25000.0,
                },
            ],
        };

        twin.record_telemetry(frame);

        assert_eq!(twin.history.len(), 1);
        let latest = twin.latest_telemetry().unwrap();
        assert_eq!(latest.timestamp_s, 0.05);
        assert_eq!(latest.layer_states.len(), 2);
    }
}
