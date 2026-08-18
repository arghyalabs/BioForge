//! Cellular Gene Expression to Spatial Tissue Morphogen Scale Bridge.

use serde::{Deserialize, Serialize};

/// Couples cellular protein translation to spatial tissue morphogen secretion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MorphogenEmissionBridge {
    /// Secretion turnover rate constant in $\text{s}^{-1}$.
    pub secretion_rate_s: f64,
}

impl Default for MorphogenEmissionBridge {
    fn default() -> Self {
        Self {
            secretion_rate_s: 0.01,
        }
    }
}

impl MorphogenEmissionBridge {
    /// Calculate spatial morphogen production source rate $S(\vec{x})$ in $\text{nM/s}$
    /// generated from intracellular protein concentration.
    #[must_use]
    pub fn compute_spatial_source_rate_nM_s(&self, intracellular_protein_nM: f64) -> f64 {
        self.secretion_rate_s * intracellular_protein_nM.max(0.0)
    }
}

/// Couples local spatial tissue morphogen concentration to intracellular transcription factor activation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MorphogenReceptorBridge {
    /// Maximum intracellular transcription factor concentration in $\text{nM}$.
    pub tf_max_nM: f64,
    /// Receptor-morphogen dissociation constant $K_d$ in nanomolar ($\text{nM}$).
    pub kd_nM: f64,
    /// Cooperativity exponent $n$.
    pub hill_n: f64,
}

impl Default for MorphogenReceptorBridge {
    fn default() -> Self {
        Self {
            tf_max_nM: 100.0,
            kd_nM: 10.0,
            hill_n: 1.0,
        }
    }
}

impl MorphogenReceptorBridge {
    /// Calculate active intracellular transcription factor concentration $[TF]$ in $\text{nM}$
    /// from local extracellular spatial morphogen concentration $C(\vec{x})$.
    #[must_use]
    pub fn compute_intracellular_tf_nM(&self, local_morphogen_nM: f64) -> f64 {
        let m = local_morphogen_nM.max(0.0);
        let m_n = m.powf(self.hill_n);
        let kd_n = self.kd_nM.powf(self.hill_n);

        if m_n + kd_n > 0.0 {
            self.tf_max_nM * (m_n / (m_n + kd_n))
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_morphogen_emission_and_receptor_activation() {
        let emitter = MorphogenEmissionBridge {
            secretion_rate_s: 0.05,
        };
        // 1000 nM protein -> 50 nM/s secretion into tissue grid
        assert_eq!(emitter.compute_spatial_source_rate_nM_s(1000.0), 50.0);

        let receptor = MorphogenReceptorBridge {
            tf_max_nM: 100.0,
            kd_nM: 10.0,
            hill_n: 1.0,
        };
        // At [M] = Kd (10 nM) -> [TF] = 0.5 * TF_max = 50 nM
        assert_eq!(receptor.compute_intracellular_tf_nM(10.0), 50.0);
    }
}
