//! Multi-fidelity dimensions, transparent approximation reporting, and error bounds (Principles 5 & 15).

use serde::{Deserialize, Serialize};

/// Computational modeling fidelity level across spatial and temporal dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FidelityLevel {
    /// Macroscopic coarse-grained continuum approximation (fastest).
    CoarseGrained,
    /// Deterministic reaction-diffusion continuum PDE / ODE.
    ContinuumApproximation,
    /// Discrete stochastic exact particle simulation (e.g. Gillespie SSA).
    DiscreteStochastic,
    /// Explicit all-atom molecular mechanics with symplectic integration (highest fidelity).
    AllAtomExplicit,
}

/// Specification of modeling fidelity and approximations for a single biological scale.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScaleFidelitySpec {
    /// Biological scale name.
    pub scale_name: String,
    /// Fidelity level used.
    pub fidelity_level: FidelityLevel,
    /// Spatial resolution in nanometers ($\text{nm}$).
    pub spatial_resolution_nm: f64,
    /// Temporal resolution (integration step) in seconds ($\text{s}$).
    pub temporal_resolution_s: f64,
    /// Exposed biological assumptions (Principle 4).
    pub assumptions: Vec<String>,
    /// Exposed mathematical approximations (Principle 5).
    pub approximations: Vec<String>,
    /// Expected numerical error upper bound in percentage ($\%$).
    pub expected_error_bound_percent: f64,
}

/// Multi-fidelity audit report documenting computational precision across all twin layers.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MultiFidelityReport {
    /// Digital twin identifier.
    pub twin_id: String,
    /// Fidelity specifications across each scale.
    pub scale_specs: Vec<ScaleFidelitySpec>,
}

impl MultiFidelityReport {
    /// Create a new multi-fidelity audit report.
    #[must_use]
    pub fn new(twin_id: impl Into<String>) -> Self {
        Self {
            twin_id: twin_id.into(),
            scale_specs: Vec::new(),
        }
    }

    /// Add a scale fidelity specification.
    pub fn add_spec(&mut self, spec: ScaleFidelitySpec) {
        self.scale_specs.push(spec);
    }

    /// Render a human-readable Markdown summary report.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut md = format!("# Multi-Fidelity Audit Report: {}\n\n", self.twin_id);
        md.push_str("| Biological Scale | Fidelity Level | Spatial Res | Temporal Step | Error Bound (±%) |\n");
        md.push_str("| :--- | :--- | :--- | :--- | :--- |\n");

        for s in &self.scale_specs {
            md.push_str(&format!(
                "| {} | {:?} | {:.2} nm | {:.2e} s | {:.2}% |\n",
                s.scale_name,
                s.fidelity_level,
                s.spatial_resolution_nm,
                s.temporal_resolution_s,
                s.expected_error_bound_percent
            ));
        }

        md
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_fidelity_report_generation() {
        let mut report = MultiFidelityReport::new("TWIN-001");

        report.add_spec(ScaleFidelitySpec {
            scale_name: "Molecular Dynamics".to_string(),
            fidelity_level: FidelityLevel::AllAtomExplicit,
            spatial_resolution_nm: 0.1,
            temporal_resolution_s: 0.5e-15,
            assumptions: vec!["CHARMM27 harmonic force field".to_string()],
            approximations: vec!["Non-bonded cutoff at 1.2 nm".to_string()],
            expected_error_bound_percent: 0.05,
        });

        report.add_spec(ScaleFidelitySpec {
            scale_name: "Cellular GRN".to_string(),
            fidelity_level: FidelityLevel::ContinuumApproximation,
            spatial_resolution_nm: 1000.0,
            temporal_resolution_s: 1.0,
            assumptions: vec!["Well-mixed cytoplasm".to_string()],
            approximations: vec!["Quasi-steady-state Hill promoter binding".to_string()],
            expected_error_bound_percent: 1.5,
        });

        let md = report.render_markdown();
        assert!(md.contains("TWIN-001"));
        assert!(md.contains("Molecular Dynamics"));
        assert!(md.contains("Cellular GRN"));
    }
}
