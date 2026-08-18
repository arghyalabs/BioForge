//! Interactive simulation controller with time-scrubbing, rewind, and live perturbation steering.

use serde::{Deserialize, Serialize};

use bioforge_hypothesis::Perturbation;

use crate::error::TwinError;

/// Live interactive simulation playback commands (Phase 22).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SimulationCommand {
    /// Resume forward execution.
    Play,
    /// Pause execution.
    Pause,
    /// Advance execution by $N$ discrete time steps.
    Step(usize),
    /// Rewind simulation state to a historical checkpoint step.
    Rewind(usize),
    /// Seek simulation state to a specific timestamp in seconds.
    SeekTime(f64),
    /// Dynamically inject an intervention/perturbation without restarting the simulation.
    InjectPerturbation(Perturbation),
}

/// Current playback state of the simulation session.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PlaybackState {
    /// Simulation is actively advancing.
    Running,
    /// Simulation is paused at the current state.
    Paused,
    /// Target duration has elapsed.
    Completed,
}

/// An interactive simulation session with historical checkpoint ring-buffering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InteractiveSession<T: Clone> {
    /// Current step number.
    pub current_step: usize,
    /// Current physical timestamp in seconds.
    pub current_time_s: f64,
    /// Current playback status.
    pub state: PlaybackState,
    /// Ring buffer of historical state checkpoints: `(step, timestamp_s, state)`.
    pub checkpoints: Vec<(usize, f64, T)>,
    /// Active perturbations dynamically steering the simulation.
    pub active_perturbations: Vec<Perturbation>,
}

impl<T: Clone> InteractiveSession<T> {
    /// Create a new interactive session initialized at $t=0$, step $0$.
    #[must_use]
    pub fn new(initial_state: T) -> Self {
        Self {
            current_step: 0,
            current_time_s: 0.0,
            state: PlaybackState::Paused,
            checkpoints: vec![(0, 0.0, initial_state)],
            active_perturbations: Vec::new(),
        }
    }

    /// Record a simulation checkpoint into the history buffer.
    pub fn record_checkpoint(&mut self, step: usize, time_s: f64, state: T) {
        self.current_step = step;
        self.current_time_s = time_s;
        self.checkpoints.push((step, time_s, state));
    }

    /// Execute an interactive playback command.
    pub fn execute_command(
        &mut self,
        cmd: SimulationCommand,
    ) -> Result<Option<T>, TwinError> {
        match cmd {
            SimulationCommand::Play => {
                self.state = PlaybackState::Running;
                Ok(None)
            }
            SimulationCommand::Pause => {
                self.state = PlaybackState::Paused;
                Ok(None)
            }
            SimulationCommand::Step(_n) => {
                self.state = PlaybackState::Paused;
                Ok(None)
            }
            SimulationCommand::Rewind(target_step) => {
                if let Some((step, time_s, state)) = self
                    .checkpoints
                    .iter()
                    .find(|(s, _, _)| *s == target_step)
                    .cloned()
                {
                    self.current_step = step;
                    self.current_time_s = time_s;
                    self.state = PlaybackState::Paused;
                    // Truncate forward history after rewind point
                    self.checkpoints.retain(|(s, _, _)| *s <= target_step);
                    Ok(Some(state))
                } else {
                    Err(TwinError::CheckpointNotFound { step: target_step })
                }
            }
            SimulationCommand::SeekTime(target_time_s) => {
                // Find closest historical checkpoint
                if let Some((step, time_s, state)) = self
                    .checkpoints
                    .iter()
                    .min_by(|(_, t1, _), (_, t2, _)| {
                        (t1 - target_time_s)
                            .abs()
                            .partial_cmp(&(t2 - target_time_s).abs())
                            .unwrap()
                    })
                    .cloned()
                {
                    self.current_step = step;
                    self.current_time_s = time_s;
                    self.state = PlaybackState::Paused;
                    Ok(Some(state))
                } else {
                    Err(TwinError::CheckpointNotFound { step: 0 })
                }
            }
            SimulationCommand::InjectPerturbation(pert) => {
                self.active_perturbations.push(pert);
                Ok(None)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interactive_session_rewind_and_perturbation() {
        let mut session = InteractiveSession::new(100.0); // State is f64

        // Advance 3 checkpoints
        session.record_checkpoint(1, 0.1, 90.0);
        session.record_checkpoint(2, 0.2, 80.0);
        session.record_checkpoint(3, 0.3, 70.0);

        assert_eq!(session.current_step, 3);
        assert_eq!(session.checkpoints.len(), 4);

        // Rewind to step 1
        let restored = session
            .execute_command(SimulationCommand::Rewind(1))
            .unwrap()
            .unwrap();

        assert_eq!(restored, 90.0);
        assert_eq!(session.current_step, 1);
        assert_eq!(session.current_time_s, 0.1);
        assert_eq!(session.checkpoints.len(), 2); // Step 0 and 1 preserved

        // Inject live drug perturbation
        let drug = Perturbation::DrugInhibition {
            drug_name: "Gleevec".to_string(),
            target_enzyme: "Abl".to_string(),
            concentration_nM: 100.0,
            ic50_nM: 25.0,
            hill_coeff: 1.0,
        };
        session
            .execute_command(SimulationCommand::InjectPerturbation(drug))
            .unwrap();

        assert_eq!(session.active_perturbations.len(), 1);
    }
}
