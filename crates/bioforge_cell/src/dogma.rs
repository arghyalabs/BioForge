//! Quantitative Central Dogma kinetics (transcription, translation, and half-life decay).

use serde::{Deserialize, Serialize};

/// Dynamic parameters governing the transcription and translation of a single gene.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeneExpression {
    /// Gene and protein descriptor (e.g. "GFP", "LacI", "p53").
    pub name: String,
    /// Basal transcription rate $k_{\text{tx}}$ in nanomolar per second ($\text{nM/s}$).
    pub transcription_rate: f64,
    /// mRNA degradation rate constant $\delta_m$ in $\text{s}^{-1}$.
    pub mrna_degradation_rate: f64,
    /// Translation rate constant $k_{\text{tl}}$ in $\text{s}^{-1}$ (protein molecules per mRNA per second).
    pub translation_rate: f64,
    /// Protein degradation rate constant $\delta_p$ in $\text{s}^{-1}$.
    pub protein_degradation_rate: f64,
}

impl Default for GeneExpression {
    fn default() -> Self {
        Self {
            name: "Gene".to_string(),
            transcription_rate: 0.1,        // 0.1 nM/s
            mrna_degradation_rate: 0.002,   // t_1/2 ~ 5.8 min
            translation_rate: 0.05,         // 0.05 protein/mRNA/s
            protein_degradation_rate: 0.0001, // t_1/2 ~ 1.9 hours
        }
    }
}

impl GeneExpression {
    /// Create a new gene expression kinetic unit.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        k_tx_nM_s: f64,
        delta_m_s: f64,
        k_tl_s: f64,
        delta_p_s: f64,
    ) -> Self {
        Self {
            name: name.into(),
            transcription_rate: k_tx_nM_s.max(0.0),
            mrna_degradation_rate: delta_m_s.max(1e-12),
            translation_rate: k_tl_s.max(0.0),
            protein_degradation_rate: delta_p_s.max(1e-12),
        }
    }

    /// Theoretical mRNA half-life $t_{1/2} = \frac{\ln(2)}{\delta_m}$ in seconds.
    #[must_use]
    pub fn mrna_half_life_seconds(&self) -> f64 {
        2.0_f64.ln() / self.mrna_degradation_rate
    }

    /// Theoretical protein half-life $t_{1/2} = \frac{\ln(2)}{\delta_p}$ in seconds.
    #[must_use]
    pub fn protein_half_life_seconds(&self) -> f64 {
        2.0_f64.ln() / self.protein_degradation_rate
    }

    /// Analytical steady-state mRNA concentration $[\text{mRNA}]^* = \frac{k_{\text{tx}}}{\delta_m}$ in $\text{nM}$.
    #[must_use]
    pub fn steady_state_mrna_nM(&self) -> f64 {
        self.transcription_rate / self.mrna_degradation_rate
    }

    /// Analytical steady-state protein concentration $[\text{Protein}]^* = \frac{k_{\text{tl}} k_{\text{tx}}}{\delta_m \delta_p}$ in $\text{nM}$.
    #[must_use]
    pub fn steady_state_protein_nM(&self) -> f64 {
        (self.translation_rate * self.transcription_rate)
            / (self.mrna_degradation_rate * self.protein_degradation_rate)
    }

    /// Evaluate time derivatives $(d[\text{mRNA}]/dt, d[\text{Protein}]/dt)$ given current concentrations.
    #[must_use]
    pub fn compute_derivatives(&self, mrna_nM: f64, protein_nM: f64) -> (f64, f64) {
        let d_mrna = self.transcription_rate - self.mrna_degradation_rate * mrna_nM.max(0.0);
        let d_protein = self.translation_rate * mrna_nM.max(0.0)
            - self.protein_degradation_rate * protein_nM.max(0.0);
        (d_mrna, d_protein)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_central_dogma_analytical_steady_states() {
        let gene = GeneExpression::new("GFP", 0.1, 0.002, 0.05, 0.0001);

        // [mRNA]* = 0.1 / 0.002 = 50.0 nM
        assert!((gene.steady_state_mrna_nM() - 50.0).abs() < 1e-9);

        // [Protein]* = (0.05 * 0.1) / (0.002 * 0.0001) = 0.005 / 2e-7 = 25,000 nM = 25.0 uM
        assert!((gene.steady_state_protein_nM() - 25_000.0).abs() < 1e-9);

        // At steady state, derivatives must be zero
        let (d_mrna, d_prot) = gene.compute_derivatives(50.0, 25_000.0);
        assert!(d_mrna.abs() < 1e-12);
        assert!(d_prot.abs() < 1e-12);
    }
}
