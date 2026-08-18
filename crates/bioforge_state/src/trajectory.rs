//! Trajectory recording and standard format exporters (XYZ, CSV).
//!
//! A [`Trajectory`] records physical state snapshots over time during
//! numerical simulation.

use crate::state::SimulationState;

/// A single snapshot frame in a physical simulation trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    /// Time in picoseconds ($\text{ps}$).
    pub time: f64,
    /// Integration step.
    pub step: u64,
    /// Atomic positions $[x, y, z]$ in Ångströms ($\text{Å}$).
    pub positions: Vec<[f64; 3]>,
    /// Instantaneous kinetic energy in $\text{kJ/mol}$.
    pub kinetic_energy: f64,
    /// Instantaneous potential energy in $\text{kJ/mol}$.
    pub potential_energy: f64,
    /// Instantaneous kinetic temperature in Kelvin ($\text{K}$).
    pub temperature: f64,
}

impl TrajectoryFrame {
    /// Total instantaneous energy ($E_{\text{total}} = E_{\text{kin}} + E_{\text{pot}}$) in $\text{kJ/mol}$.
    #[must_use]
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy + self.potential_energy
    }
}

/// A historical buffer of simulation snapshots.
#[derive(Debug, Clone, PartialEq)]
pub struct Trajectory {
    /// Recorded snapshot frames.
    pub frames: Vec<TrajectoryFrame>,
    /// Stride interval for recording frames (e.g. record every $N$ steps).
    pub sample_interval_steps: u64,
}

impl Trajectory {
    /// Create a new trajectory buffer.
    #[must_use]
    pub fn new(sample_interval_steps: u64) -> Self {
        Self {
            frames: Vec::new(),
            sample_interval_steps: sample_interval_steps.max(1),
        }
    }

    /// Record a snapshot from the current simulation state.
    pub fn record_frame(&mut self, state: &SimulationState, potential_energy: f64) {
        self.frames.push(TrajectoryFrame {
            time: state.time,
            step: state.step,
            positions: state.positions.clone(),
            kinetic_energy: state.kinetic_energy(),
            potential_energy,
            temperature: state.instantaneous_temperature(),
        });
    }

    /// Number of recorded frames.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether the trajectory contains zero frames.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Get reference to the latest recorded frame.
    #[must_use]
    pub fn last_frame(&self) -> Option<&TrajectoryFrame> {
        self.frames.last()
    }

    /// Get reference to frame at `index`.
    #[must_use]
    pub fn frame(&self, index: usize) -> Option<&TrajectoryFrame> {
        self.frames.get(index)
    }

    /// Export the trajectory into standard multi-frame XYZ format.
    ///
    /// Compatible with computational chemistry visualizers (VMD, PyMOL, Ovito).
    #[must_use]
    pub fn to_xyz(&self, state: &SimulationState) -> String {
        let mut out = String::new();
        let num_atoms = state.num_atoms;

        for frame in &self.frames {
            out.push_str(&format!("{}\n", num_atoms));
            out.push_str(&format!(
                "Time={:.3} ps, Step={}, T={:.1} K, E_tot={:.3} kJ/mol\n",
                frame.time,
                frame.step,
                frame.temperature,
                frame.total_energy()
            ));

            for i in 0..num_atoms {
                let symbol = if i < state.elements.len() {
                    state.elements[i].symbol
                } else {
                    "X"
                };
                let pos = frame.positions[i];
                out.push_str(&format!(
                    "{:<2} {:12.6} {:12.6} {:12.6}\n",
                    symbol, pos[0], pos[1], pos[2]
                ));
            }
        }

        out
    }

    /// Export thermodynamic time series data as CSV (time, step, T, E_kin, E_pot, E_total).
    #[must_use]
    pub fn to_csv(&self) -> String {
        let mut out = String::from("time_ps,step,temp_k,e_kin_kj_mol,e_pot_kj_mol,e_tot_kj_mol\n");
        for f in &self.frames {
            out.push_str(&format!(
                "{:.4},{},{:.2},{:.4},{:.4},{:.4}\n",
                f.time,
                f.step,
                f.temperature,
                f.kinetic_energy,
                f.potential_energy,
                f.total_energy()
            ));
        }
        out
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::Element;

    #[test]
    fn test_trajectory_recording_and_energy() {
        let mut state = SimulationState::empty();
        state.num_atoms = 1;
        state.masses = vec![12.011];
        state.velocities = vec![[1.0, 0.0, 0.0]];
        state.positions = vec![[0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0]];
        state.charges = vec![0.0];
        state.elements = vec![Element::from_symbol("C").unwrap()];

        let mut traj = Trajectory::new(1);
        assert!(traj.is_empty());

        traj.record_frame(&state, 10.0);
        assert_eq!(traj.len(), 1);

        let f = traj.last_frame().unwrap();
        assert_eq!(f.step, 0);
        assert!((f.potential_energy - 10.0).abs() < 1e-6);
        assert!((f.total_energy() - (f.kinetic_energy + 10.0)).abs() < 1e-6);
    }

    #[test]
    fn test_to_xyz_export() {
        let mut state = SimulationState::empty();
        state.num_atoms = 1;
        state.masses = vec![12.011];
        state.positions = vec![[1.5, 2.5, 3.5]];
        state.velocities = vec![[0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0]];
        state.charges = vec![0.0];
        state.elements = vec![Element::from_symbol("C").unwrap()];

        let mut traj = Trajectory::new(1);
        traj.record_frame(&state, 0.0);

        let xyz = traj.to_xyz(&state);
        assert!(xyz.starts_with("1\n"));
        assert!(xyz.contains("C "));
        assert!(xyz.contains("1.500000"));
        assert!(xyz.contains("2.500000"));
        assert!(xyz.contains("3.500000"));
    }

    #[test]
    fn test_to_csv_export() {
        let mut state = SimulationState::empty();
        state.num_atoms = 1;
        state.masses = vec![12.011];
        state.positions = vec![[0.0, 0.0, 0.0]];
        state.velocities = vec![[0.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0]];
        state.charges = vec![0.0];
        state.elements = vec![Element::from_symbol("C").unwrap()];

        let mut traj = Trajectory::new(1);
        traj.record_frame(&state, 5.5);

        let csv = traj.to_csv();
        assert!(csv.contains("time_ps,step,temp_k,e_kin_kj_mol,e_pot_kj_mol,e_tot_kj_mol"));
        assert!(csv.contains("5.5000"));
    }
}
