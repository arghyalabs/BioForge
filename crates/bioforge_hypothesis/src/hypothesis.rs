//! Scientific hypothesis, explicit assumptions, causal chains, and empirical validation (Principles 4, 5, 8, 9).

use serde::{Deserialize, Serialize};

/// Type of evidence supporting or testing a scientific claim (Principle 8).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum EvidenceKind {
    /// Direct empirical wet-lab experimental measurement.
    WetLabMeasurement {
        /// Digital Object Identifier or citation.
        doi: String,
        /// Number of experimental biological replicates $N$.
        sample_size: usize,
        /// Measured mean value.
        mean: f64,
        /// Standard error of the mean ($\text{SEM}$).
        std_error: f64,
    },

    /// Computational simulation prediction (distinct from empirical fact per Principle 9).
    SimulationPrediction {
        /// Computational model name.
        model_name: String,
        /// Predicted numerical value.
        value: f64,
    },

    /// Qualitative literature prior knowledge.
    LiteraturePrior {
        /// Citation reference.
        citation: String,
        /// Finding description.
        summary: String,
    },
}

/// Verification status of a hypothesis prediction against experimental evidence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PredictionStatus {
    /// Not yet tested against empirical data.
    Unverified,
    /// Verified by empirical data within statistical confidence ($Z \le 2.0$).
    Verified {
        /// Calculated standard score distance.
        z_score: f64,
    },
    /// Refuted by empirical data ($Z > 3.0$).
    Refuted {
        /// Calculated standard score distance.
        z_score: f64,
    },
    /// Inconclusive statistical test ($2.0 < Z \le 3.0$).
    Inconclusive {
        /// Calculated standard score distance.
        z_score: f64,
    },
}

/// A quantitative prediction derived from a mechanistic biological hypothesis.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Prediction {
    /// Prediction ID.
    pub id: usize,
    /// Physical observable name (e.g. "[Protein_GFP]", "action_potential_peak_mV").
    pub target_observable: String,
    /// Numerically predicted value.
    pub predicted_value: f64,
    /// Empirical verification status.
    pub status: PredictionStatus,
}

impl Prediction {
    /// Create a new unverified prediction.
    #[must_use]
    pub fn new(id: usize, target_observable: impl Into<String>, predicted_value: f64) -> Self {
        Self {
            id,
            target_observable: target_observable.into(),
            predicted_value,
            status: PredictionStatus::Unverified,
        }
    }

    /// Statistically test this prediction against empirical experimental data ($\mu_{\text{exp}} \pm \text{SEM}$).
    pub fn evaluate_against_experimental(&mut self, exp_mean: f64, exp_std_error: f64) {
        let sem = exp_std_error.max(1e-9);
        let z = (self.predicted_value - exp_mean).abs() / sem;

        self.status = if z <= 2.0 {
            PredictionStatus::Verified { z_score: z }
        } else if z > 3.0 {
            PredictionStatus::Refuted { z_score: z }
        } else {
            PredictionStatus::Inconclusive { z_score: z }
        };
    }
}

/// A mechanistic causal link in a biological hypothesis chain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CausalStep {
    /// Step order in causal sequence.
    pub step_index: usize,
    /// Direct upstream causal trigger (e.g. "Kinase Phosphorylation").
    pub cause: String,
    /// Direct downstream biological effect (e.g. "Promoter Activation").
    pub effect: String,
    /// Governing mathematical rate law or mechanism description.
    pub mechanism: String,
    /// Supporting empirical or computational evidence.
    pub evidence: Vec<EvidenceKind>,
}

/// A scientific biological hypothesis with explicit assumptions, causal chains, and predictions.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Hypothesis {
    /// Unique hypothesis identifier.
    pub id: String,
    /// Concise title.
    pub title: String,
    /// Detailed description.
    pub description: String,
    /// Explicit biological and physical assumptions (Scientific Principle 4).
    pub assumptions: Vec<String>,
    /// Explicit mathematical and numerical approximations (Scientific Principle 5).
    pub approximations: Vec<String>,
    /// Stepwise mechanistic causal chain.
    pub causal_chain: Vec<CausalStep>,
    /// Quantitative predictions derived from this hypothesis.
    pub predictions: Vec<Prediction>,
}

impl Hypothesis {
    /// Create a new biological hypothesis with exposed assumptions and approximations.
    #[must_use]
    pub fn new(id: impl Into<String>, title: impl Into<String>, description: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            title: title.into(),
            description: description.into(),
            assumptions: Vec::new(),
            approximations: Vec::new(),
            causal_chain: Vec::new(),
            predictions: Vec::new(),
        }
    }

    /// Add an explicit assumption to the hypothesis record.
    pub fn add_assumption(&mut self, assumption: impl Into<String>) {
        self.assumptions.push(assumption.into());
    }

    /// Add an explicit approximation to the hypothesis record.
    pub fn add_approximation(&mut self, approximation: impl Into<String>) {
        self.approximations.push(approximation.into());
    }

    /// Append a causal step to the hypothesis chain.
    pub fn add_causal_step(
        &mut self,
        cause: impl Into<String>,
        effect: impl Into<String>,
        mechanism: impl Into<String>,
    ) {
        let step_index = self.causal_chain.len();
        self.causal_chain.push(CausalStep {
            step_index,
            cause: cause.into(),
            effect: effect.into(),
            mechanism: mechanism.into(),
            evidence: Vec::new(),
        });
    }

    /// Add a testable prediction.
    pub fn add_prediction(&mut self, target_observable: impl Into<String>, predicted_value: f64) -> usize {
        let id = self.predictions.len();
        self.predictions.push(Prediction::new(id, target_observable, predicted_value));
        id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hypothesis_prediction_verification_z_score() {
        let mut hyp = Hypothesis::new(
            "HYP-001",
            "p53 Activation of p21",
            "p53 transcriptionally activates p21 cyclin inhibitor",
        );

        hyp.add_assumption("Cell is in G1 checkpoint arrest");
        hyp.add_approximation("Hill coefficient n = 2.0 quasi-steady state promoter binding");
        hyp.add_causal_step("p53 Phosphorylation", "p21 Transcription", "Hill Activation");

        let pred1 = hyp.add_prediction("[p21_mRNA]", 103.0); // 103 nM predicted
        let pred2 = hyp.add_prediction("[p21_mRNA]", 160.0); // 160 nM predicted (bad prediction)

        // Experimental ground truth: 100.0 nM +/- 2.0 nM SEM
        hyp.predictions[pred1].evaluate_against_experimental(100.0, 2.0);
        // Z = (103 - 100) / 2 = 1.5 <= 2.0 -> Verified
        assert!(matches!(
            hyp.predictions[pred1].status,
            PredictionStatus::Verified { .. }
        ));

        hyp.predictions[pred2].evaluate_against_experimental(100.0, 2.0);
        // Z = (160 - 100) / 2 = 30.0 > 3.0 -> Refuted
        assert!(matches!(
            hyp.predictions[pred2].status,
            PredictionStatus::Refuted { .. }
        ));
    }
}
