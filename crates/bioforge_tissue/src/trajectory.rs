//! Spatial trajectory recordings and grid CSV export for tissue simulations.

use serde::{Deserialize, Serialize};

/// 1D Spatial profile recording over time.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct SpatialTrajectory1D {
    /// Recorded simulation timestamps in seconds.
    pub times_s: Vec<f64>,
    /// Spatial grid spacing $\Delta x$ in $\mu\text{m}$.
    pub dx_um: f64,
    /// Spatial profiles across time: `profiles[time_idx][spatial_x_idx]`.
    pub profiles: Vec<Vec<f64>>,
}

impl SpatialTrajectory1D {
    /// Create a new spatial trajectory buffer.
    #[must_use]
    pub fn new(dx_um: f64) -> Self {
        Self {
            times_s: Vec::new(),
            dx_um,
            profiles: Vec::new(),
        }
    }

    /// Record a spatial frame.
    pub fn record(&mut self, time_s: f64, values: &[f64]) {
        self.times_s.push(time_s);
        self.profiles.push(values.to_vec());
    }

    /// Export spatial history to CSV heatmap format.
    #[must_use]
    pub fn export_csv(&self) -> String {
        let mut out = String::from("time_s");
        if let Some(first) = self.profiles.first() {
            for i in 0..first.len() {
                let x_um = i as f64 * self.dx_um;
                out.push_str(&format!(",\"x_{:.1}um\"", x_um));
            }
        }
        out.push('\n');

        for (t_idx, &t) in self.times_s.iter().enumerate() {
            out.push_str(&format!("{:.4}", t));
            if let Some(row) = self.profiles.get(t_idx) {
                for &val in row {
                    out.push_str(&format!(",{}", val));
                }
            }
            out.push('\n');
        }

        out
    }
}
