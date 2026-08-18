//! Numerical integrators and trajectory recordings for cellular simulations.

use serde::{Deserialize, Serialize};

use crate::error::CellError;

/// Recorded trajectory for cellular concentrations over time.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct CellSimulationTrajectory {
    /// Time points in seconds ($\text{s}$).
    pub times_s: Vec<f64>,
    /// Tracked species / protein names.
    pub names: Vec<String>,
    /// Concentrations over time: `values[species_idx][time_idx]` in $\text{nM}$.
    pub values: Vec<Vec<f64>>,
}

impl CellSimulationTrajectory {
    /// Create a new trajectory buffer.
    #[must_use]
    pub fn new(names: Vec<String>) -> Self {
        let n = names.len();
        Self {
            times_s: Vec::new(),
            names,
            values: vec![Vec::new(); n],
        }
    }

    /// Record a state frame.
    pub fn record(&mut self, time_s: f64, frame_values: &[f64]) {
        self.times_s.push(time_s);
        for (i, &v) in frame_values.iter().enumerate() {
            if i < self.values.len() {
                self.values[i].push(v);
            }
        }
    }

    /// Export recording to CSV format.
    #[must_use]
    pub fn export_csv(&self) -> String {
        let mut out = String::from("time_s");
        for name in &self.names {
            out.push_str(&format!(",\"[{}](nM)\"", name));
        }
        out.push('\n');

        let num_rows = self.times_s.len();
        for r in 0..num_rows {
            out.push_str(&format!("{:.4}", self.times_s[r]));
            for s_idx in 0..self.values.len() {
                let val = if r < self.values[s_idx].len() {
                    self.values[s_idx][r]
                } else {
                    0.0
                };
                out.push_str(&format!(",{}", val));
            }
            out.push('\n');
        }
        out
    }
}

/// Generic 4th-order Runge-Kutta integrator for coupled cellular systems.
pub fn integrate_rk4<F>(
    initial_state: &[f64],
    total_time_s: f64,
    dt_s: f64,
    names: Vec<String>,
    mut derivative_fn: F,
) -> Result<CellSimulationTrajectory, CellError>
where
    F: FnMut(&[f64], f64) -> Vec<f64>,
{
    let mut traj = CellSimulationTrajectory::new(names);
    let mut state = initial_state.to_vec();
    let n = state.len();

    let mut t = 0.0;
    let dt = dt_s.max(1e-6);
    let record_interval = dt * 10.0;
    let mut next_record_t = 0.0;

    traj.record(0.0, &state);

    while t < total_time_s - 1e-12 {
        let step_dt = dt.min(total_time_s - t);

        // k1
        let k1 = derivative_fn(&state, t);

        // k2
        let mut s2 = vec![0.0; n];
        for i in 0..n {
            s2[i] = (state[i] + 0.5 * step_dt * k1[i]).max(0.0);
        }
        let k2 = derivative_fn(&s2, t + 0.5 * step_dt);

        // k3
        let mut s3 = vec![0.0; n];
        for i in 0..n {
            s3[i] = (state[i] + 0.5 * step_dt * k2[i]).max(0.0);
        }
        let k3 = derivative_fn(&s3, t + 0.5 * step_dt);

        // k4
        let mut s4 = vec![0.0; n];
        for i in 0..n {
            s4[i] = (state[i] + step_dt * k3[i]).max(0.0);
        }
        let k4 = derivative_fn(&s4, t + step_dt);

        // Update
        for i in 0..n {
            state[i] = (state[i] + (step_dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])).max(0.0);
        }

        t += step_dt;

        if t >= next_record_t || t >= total_time_s - 1e-12 {
            traj.record(t, &state);
            next_record_t += record_interval;
        }
    }

    Ok(traj)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rk4_decay_integration() {
        // Simple decay: dy/dt = -0.1 * y
        let names = vec!["Y".to_string()];
        let traj = integrate_rk4(&[100.0], 10.0, 0.01, names, |state, _t| {
            vec![-0.1 * state[0]]
        })
        .unwrap();

        assert!(!traj.times_s.is_empty());
        let final_val = traj.values[0].last().copied().unwrap();
        // Theoretical: 100 * exp(-1.0) = 36.7879
        let expected = 100.0 * (-1.0_f64).exp();
        assert!((final_val - expected).abs() < 0.1);
    }
}
