//! Action potential analysis, spike detection, and firing frequency metrics.

use serde::{Deserialize, Serialize};

use crate::hodgkin_huxley::ElectrophysiologyTrajectory;

/// Quantitative metrics characterizing neuronal action potentials.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ActionPotentialMetrics {
    /// Total number of detected action potential spikes.
    pub spike_count: usize,
    /// Exact timestamps of peak action potentials in milliseconds ($\text{ms}$).
    pub spike_times_ms: Vec<f64>,
    /// Mean firing frequency in Hertz ($\text{Hz} = \text{spikes/s}$).
    pub mean_frequency_hz: f64,
    /// Maximum depolarization peak membrane potential in millivolts ($\text{mV}$).
    pub peak_voltage_mv: f64,
    /// Baseline resting potential in millivolts ($\text{mV}$).
    pub resting_voltage_mv: f64,
    /// Peak-to-trough spike amplitude in millivolts ($\text{mV}$).
    pub amplitude_mv: f64,
}

/// Action potential spike detector using positive voltage threshold crossings.
#[derive(Debug, Clone, PartialEq)]
pub struct SpikeDetector {
    /// Voltage threshold for spike classification (default $0.0\text{ mV}$).
    pub threshold_mv: f64,
    /// Minimum refractory time between consecutive spikes in milliseconds (default $2.0\text{ ms}$).
    pub min_interval_ms: f64,
}

impl Default for SpikeDetector {
    fn default() -> Self {
        Self {
            threshold_mv: 0.0,
            min_interval_ms: 2.0,
        }
    }
}

impl SpikeDetector {
    /// Create a new spike detector.
    #[must_use]
    pub fn new(threshold_mv: f64, min_interval_ms: f64) -> Self {
        Self {
            threshold_mv,
            min_interval_ms: min_interval_ms.max(0.1),
        }
    }

    /// Analyze an electrophysiology trajectory and extract action potential metrics.
    #[must_use]
    pub fn analyze(&self, traj: &ElectrophysiologyTrajectory) -> ActionPotentialMetrics {
        if traj.times_ms.is_empty() {
            return ActionPotentialMetrics {
                spike_count: 0,
                spike_times_ms: Vec::new(),
                mean_frequency_hz: 0.0,
                peak_voltage_mv: 0.0,
                resting_voltage_mv: 0.0,
                amplitude_mv: 0.0,
            };
        }

        let resting_voltage_mv = traj.v_membrane_mv[0];
        let mut peak_voltage_mv = resting_voltage_mv;
        let mut spike_times = Vec::new();

        let n = traj.times_ms.len();
        let mut i = 0;

        while i < n {
            let v = traj.v_membrane_mv[i];
            if v > peak_voltage_mv {
                peak_voltage_mv = v;
            }

            // Check if voltage crossed threshold from below
            if v >= self.threshold_mv {
                let t_cross = traj.times_ms[i];
                let can_add = match spike_times.last() {
                    Some(&last_t) => (t_cross - last_t) >= self.min_interval_ms,
                    None => true,
                };

                if can_add {
                    // Search forward for the local maximum peak of this spike
                    let mut local_max_v = v;
                    let mut local_max_t = t_cross;
                    let mut j = i;

                    while j < n && traj.times_ms[j] - t_cross < self.min_interval_ms {
                        let vj = traj.v_membrane_mv[j];
                        if vj > local_max_v {
                            local_max_v = vj;
                            local_max_t = traj.times_ms[j];
                        }
                        j += 1;
                    }

                    if local_max_v > peak_voltage_mv {
                        peak_voltage_mv = local_max_v;
                    }

                    spike_times.push(local_max_t);
                    i = j;
                    continue;
                }
            }
            i += 1;
        }

        let spike_count = spike_times.len();
        let total_time_sec = (traj.times_ms.last().copied().unwrap_or(0.0) - traj.times_ms[0]) / 1000.0;
        let mean_frequency_hz = if total_time_sec > 0.0 && spike_count > 0 {
            (spike_count as f64) / total_time_sec
        } else {
            0.0
        };

        let amplitude_mv = peak_voltage_mv - resting_voltage_mv;

        ActionPotentialMetrics {
            spike_count,
            spike_times_ms: spike_times,
            mean_frequency_hz,
            peak_voltage_mv,
            resting_voltage_mv,
            amplitude_mv,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hodgkin_huxley::HodgkinHuxleyModel;
    use crate::stimulus::StimulusProtocol;

    #[test]
    fn test_spike_detection_and_frequency() {
        let model = HodgkinHuxleyModel::standard();
        // 50 ms sustained current step -> triggers multiple action potentials
        let stimulus = StimulusProtocol::CurrentStep {
            amplitude_uA_cm2: 10.0,
            start_ms: 5.0,
        };

        let traj = model.simulate(50.0, 0.01, &stimulus).unwrap();
        let detector = SpikeDetector::default();
        let metrics = detector.analyze(&traj);

        // A 10 uA/cm^2 sustained step in Hodgkin-Huxley produces regular firing (approx 3 spikes in 45 ms ~ 60-70 Hz)
        assert!(metrics.spike_count >= 2, "got {} spikes", metrics.spike_count);
        assert!(metrics.mean_frequency_hz > 30.0);
        assert!(metrics.amplitude_mv > 80.0);
    }
}
