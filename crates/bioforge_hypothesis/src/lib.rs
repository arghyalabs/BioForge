//! # BioForge Hypothesis (`bioforge_hypothesis`)
//!
//! Scientific hypothesis formulation, counterfactual experiment branching, and provenance for BioForge.
//!
//! ## Scientific Architecture (Principles 4, 5, 8, 9, 10, 11)
//!
//! Implements rigorous scientific reasoning and counterfactual analysis:
//! - **Hypotheses & Causal Chains**: Exposes assumptions, approximations, and predicts testable observables.
//! - **Prediction Verification**: Statistical $Z$-score testing distinguishing empirical facts from predictions.
//! - **Counterfactual Experiment Branching**: Side-by-side Wild-Type vs Knockout / Drug Perturbation analysis.
//! - **Scientific Provenance & Reproducibility**: W3C PROV-O audit trails and SHA-256 cryptographic reproducibility receipts.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod counterfactual;
pub mod error;
pub mod hypothesis;
pub mod perturbation;
pub mod provenance;

pub use counterfactual::DifferentialComparison;
pub use error::HypothesisError;
pub use hypothesis::{CausalStep, EvidenceKind, Hypothesis, Prediction, PredictionStatus};
pub use perturbation::Perturbation;
pub use provenance::{ProvenanceRecord, ReproducibilityReceipt};

// ─── Counterfactual Gene Knockout Benchmark ────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full In-Silico Gene Knockout Counterfactual Benchmark:
    ///
    /// Evaluates Wild-Type baseline vs Gene Knockout counterfactual branch
    /// and verifies that the differential fold-change tracks the expected exponential decay to zero.
    #[test]
    fn test_counterfactual_gene_knockout_differential_trajectory() {
        let times = vec![0.0, 1000.0, 5000.0, 10000.0, 30000.0];

        // Wild-type maintained at 25,000 nM steady-state
        let wt = vec![25000.0, 25000.0, 25000.0, 25000.0, 25000.0];

        // Knockout: decay with delta_p = 0.0001 s^-1 (t_1/2 ~ 6931 s)
        let ko: Vec<f64> = times
            .iter()
            .map(|&t| 25000.0 * (-0.0001_f64 * t).exp())
            .collect();

        let diff = DifferentialComparison::compute(&times, &wt, &ko).unwrap();

        // At t = 0 -> fold-change is 1.0 (identical initial state)
        assert_eq!(diff.fold_change[0], 1.0);

        // At t = 30,000 s -> remaining protein is exp(-3) ~ 4.98%
        let final_fc = diff.fold_change.last().copied().unwrap();
        assert!((final_fc - (-3.0_f64).exp()).abs() < 0.001);

        // Delta trajectory is strictly negative (loss of function)
        for &delta in &diff.delta_trajectory[1..] {
            assert!(delta < 0.0);
        }

        // Generate CSV output
        let csv = diff.export_csv();
        assert!(csv.contains("time_s,baseline,perturbed,delta,fold_change"));
    }
}
