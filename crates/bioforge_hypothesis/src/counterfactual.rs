//! Counterfactual simulation branching, differential trajectories, and fold-change analysis.

use serde::{Deserialize, Serialize};

use crate::error::HypothesisError;

/// Quantitative differential comparison between a Wild-Type baseline and a Counterfactual perturbed trajectory.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct DifferentialComparison {
    /// Time points in seconds.
    pub times_s: Vec<f64>,
    /// Wild-type baseline observable trajectory $O_{\text{WT}}(t)$.
    pub baseline_values: Vec<f64>,
    /// Counterfactual perturbed observable trajectory $O_{\text{pert}}(t)$.
    pub perturbed_values: Vec<f64>,
    /// Difference trajectory $\Delta O(t) = O_{\text{pert}}(t) - O_{\text{WT}}(t)$.
    pub delta_trajectory: Vec<f64>,
    /// Fold-change trajectory $\frac{O_{\text{pert}}(t)}{O_{\text{WT}}(t)}$.
    pub fold_change: Vec<f64>,
    /// Mean fold-change across entire trajectory.
    pub mean_fold_change: f64,
    /// Integrated absolute difference impact $\int_0^T |\Delta O(t)| \, dt$.
    pub integrated_impact: f64,
}

impl DifferentialComparison {
    /// Compute the differential comparison between baseline and perturbed simulation trajectories.
    pub fn compute(
        times_s: &[f64],
        baseline: &[f64],
        perturbed: &[f64],
    ) -> Result<Self, HypothesisError> {
        let n = times_s.len();
        if baseline.len() != n || perturbed.len() != n || n == 0 {
            return Err(HypothesisError::EvaluationError(
                "trajectory length mismatch in differential comparison".to_string(),
            ));
        }

        let mut delta_trajectory = Vec::with_capacity(n);
        let mut fold_change = Vec::with_capacity(n);
        let mut sum_fold = 0.0;
        let mut integrated_impact = 0.0;

        for i in 0..n {
            let b = baseline[i];
            let p = perturbed[i];
            let diff = p - b;
            delta_trajectory.push(diff);

            let fc = if b.abs() > 1e-12 { p / b } else { 1.0 };
            fold_change.push(fc);
            sum_fold += fc;

            if i > 0 {
                let dt = times_s[i] - times_s[i - 1];
                let avg_diff = 0.5 * (diff.abs() + delta_trajectory[i - 1].abs());
                integrated_impact += avg_diff * dt;
            }
        }

        let mean_fold_change = sum_fold / (n as f64);

        Ok(Self {
            times_s: times_s.to_vec(),
            baseline_values: baseline.to_vec(),
            perturbed_values: perturbed.to_vec(),
            delta_trajectory,
            fold_change,
            mean_fold_change,
            integrated_impact,
        })
    }

    /// Export differential comparison to CSV table.
    #[must_use]
    pub fn export_csv(&self) -> String {
        let mut out = String::from("time_s,baseline,perturbed,delta,fold_change\n");
        for i in 0..self.times_s.len() {
            out.push_str(&format!(
                "{:.4},{:.4},{:.4},{:.4},{:.4}\n",
                self.times_s[i],
                self.baseline_values[i],
                self.perturbed_values[i],
                self.delta_trajectory[i],
                self.fold_change[i]
            ));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_differential_trajectory_and_fold_change() {
        let times = vec![0.0, 1.0, 2.0];
        let wt = vec![10.0, 10.0, 10.0];      // 10 nM constant
        let ko = vec![10.0, 5.0, 0.0];        // Decaying under knockout

        let diff = DifferentialComparison::compute(&times, &wt, &ko).unwrap();

        assert_eq!(diff.delta_trajectory, vec![0.0, -5.0, -10.0]);
        assert_eq!(diff.fold_change, vec![1.0, 0.5, 0.0]);
        // Mean fold change: (1.0 + 0.5 + 0.0) / 3 = 0.5
        assert_eq!(diff.mean_fold_change, 0.5);
        // Integrated impact: trapezoid 1 (0 to 5) = 2.5, trapezoid 2 (5 to 10) = 7.5 => total = 10.0
        assert_eq!(diff.integrated_impact, 10.0);
    }
}
