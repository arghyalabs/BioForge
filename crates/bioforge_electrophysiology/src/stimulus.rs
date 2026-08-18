//! Electrical stimulation protocols (Current pulses, steps, and ramps).

use serde::{Deserialize, Serialize};

/// Electrical current stimulation protocol injected into a cell membrane.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum StimulusProtocol {
    /// Zero injected current (passive recovery).
    None,

    /// A rectangular current pulse of specified amplitude and duration in $\mu\text{A/cm}^2$.
    CurrentPulse {
        /// Current density amplitude in $\mu\text{A/cm}^2$.
        amplitude_uA_cm2: f64,
        /// Pulse start time in milliseconds ($\text{ms}$).
        start_ms: f64,
        /// Pulse duration in milliseconds ($\text{ms}$).
        duration_ms: f64,
    },

    /// Constant current injection after start time.
    CurrentStep {
        /// Current density amplitude in $\mu\text{A/cm}^2$.
        amplitude_uA_cm2: f64,
        /// Step start time in milliseconds ($\text{ms}$).
        start_ms: f64,
    },

    /// Linear current ramp between two amplitudes.
    Ramp {
        start_amplitude: f64,
        end_amplitude: f64,
        start_ms: f64,
        duration_ms: f64,
    },
}

impl Default for StimulusProtocol {
    fn default() -> Self {
        Self::CurrentPulse {
            amplitude_uA_cm2: 10.0,
            start_ms: 5.0,
            duration_ms: 1.0,
        }
    }
}

impl StimulusProtocol {
    /// Injected stimulus current density $I_{\text{stim}}(t)$ in $\mu\text{A/cm}^2$ at time $t$ ($\text{ms}$).
    #[must_use]
    pub fn current_at(&self, t_ms: f64) -> f64 {
        match self {
            StimulusProtocol::None => 0.0,
            StimulusProtocol::CurrentPulse {
                amplitude_uA_cm2,
                start_ms,
                duration_ms,
            } => {
                if t_ms >= *start_ms && t_ms < (*start_ms + *duration_ms) {
                    *amplitude_uA_cm2
                } else {
                    0.0
                }
            }
            StimulusProtocol::CurrentStep {
                amplitude_uA_cm2,
                start_ms,
            } => {
                if t_ms >= *start_ms {
                    *amplitude_uA_cm2
                } else {
                    0.0
                }
            }
            StimulusProtocol::Ramp {
                start_amplitude,
                end_amplitude,
                start_ms,
                duration_ms,
            } => {
                if t_ms >= *start_ms && t_ms < (*start_ms + *duration_ms) {
                    let progress = (t_ms - *start_ms) / *duration_ms;
                    start_amplitude + (end_amplitude - start_amplitude) * progress
                } else {
                    0.0
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_current_pulse_timing() {
        let pulse = StimulusProtocol::CurrentPulse {
            amplitude_uA_cm2: 15.0,
            start_ms: 10.0,
            duration_ms: 2.0,
        };

        assert_eq!(pulse.current_at(5.0), 0.0);
        assert_eq!(pulse.current_at(10.0), 15.0);
        assert_eq!(pulse.current_at(11.5), 15.0);
        assert_eq!(pulse.current_at(12.0), 0.0);
        assert_eq!(pulse.current_at(20.0), 0.0);
    }
}
