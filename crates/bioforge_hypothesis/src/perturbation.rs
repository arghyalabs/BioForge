//! Biological perturbations and pharmacological interventions for counterfactual simulations.

use serde::{Deserialize, Serialize};

/// A physical, genetic, or chemical perturbation applied to a biological system.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Perturbation {
    /// Complete genetic knockout ($k_{\text{tx}} \to 0$).
    GeneKnockout {
        /// Target gene name.
        gene_name: String,
    },

    /// Pharmacological small-molecule enzymatic inhibitor with Hill dose-response.
    ///
    /// Fractional remaining activity:
    /// $$\theta_{\text{act}}([I]) = \frac{1}{1 + \left(\frac{[I]}{\text{IC}_{50}}\right)^h}$$
    DrugInhibition {
        /// Drug chemical name (e.g. "Gleevec", "Aspirin", "Lapatinib").
        drug_name: String,
        /// Target protein/enzyme name.
        target_enzyme: String,
        /// Injected drug concentration in nanomolar ($\text{nM}$).
        concentration_nM: f64,
        /// Half-maximal inhibitory concentration $\text{IC}_{50}$ in nanomolar ($\text{nM}$).
        ic50_nM: f64,
        /// Hill slope coefficient $h$.
        hill_coeff: f64,
    },

    /// Point missense amino acid mutation in a protein sequence.
    MissenseMutation {
        /// Gene name.
        gene: String,
        /// 1-indexed amino acid residue position.
        codon_position: usize,
        /// Wild-type amino acid.
        wild_type_aa: char,
        /// Mutated amino acid.
        mutant_aa: char,
    },

    /// Environmental temperature shift.
    TemperatureShift {
        /// Target perturbed temperature in Kelvin ($\text{K}$).
        target_temp_k: f64,
    },
}

impl Perturbation {
    /// Calculate fractional remaining biological activity $\theta \in [0.0, 1.0]$.
    #[must_use]
    pub fn fractional_activity(&self) -> f64 {
        match self {
            Perturbation::GeneKnockout { .. } => 0.0,
            Perturbation::DrugInhibition {
                concentration_nM,
                ic50_nM,
                hill_coeff,
                ..
            } => {
                let i = concentration_nM.max(0.0);
                let ic50 = ic50_nM.max(1e-9);
                let ratio = (i / ic50).powf(*hill_coeff);
                1.0 / (1.0 + ratio)
            }
            _ => 1.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_drug_inhibition_dose_response_ic50() {
        let drug = Perturbation::DrugInhibition {
            drug_name: "Inhibitor-X".to_string(),
            target_enzyme: "Kinase-A".to_string(),
            concentration_nM: 50.0,
            ic50_nM: 50.0,
            hill_coeff: 1.0,
        };

        // When [I] = IC50 -> remaining activity is exactly 0.5 (50%)
        assert_eq!(drug.fractional_activity(), 0.5);

        // When [I] = 9 * IC50 (450 nM) -> remaining activity is 1 / (1 + 9) = 0.10 (10%)
        let drug_high = Perturbation::DrugInhibition {
            drug_name: "Inhibitor-X".to_string(),
            target_enzyme: "Kinase-A".to_string(),
            concentration_nM: 450.0,
            ic50_nM: 50.0,
            hill_coeff: 1.0,
        };
        assert!((drug_high.fractional_activity() - 0.10).abs() < 1e-6);

        // Knockout -> remaining activity is 0.0
        let ko = Perturbation::GeneKnockout {
            gene_name: "p53".to_string(),
        };
        assert_eq!(ko.fractional_activity(), 0.0);
    }
}
