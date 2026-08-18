//! Intracellular signal transduction cascades (MAPK / ERK ultrasensitive signaling).

use serde::{Deserialize, Serialize};

/// Classical 3-tier Mitogen-Activated Protein Kinase (MAPK / ERK) phosphorylation cascade.
///
/// Models the ultrasensitive sigmoidal response of cellular signaling (Huang & Ferrell, PNAS 1996).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapkCascade {
    /// Total MAPKKK (Raf) concentration in $\text{nM}$.
    pub total_mapkkk: f64,
    /// Total MAPKK (MEK) concentration in $\text{nM}$.
    pub total_mapkk: f64,
    /// Total MAPK (ERK) concentration in $\text{nM}$.
    pub total_mapk: f64,
    /// Activation rate for Tier 1.
    pub k_act_1: f64,
    /// Inactivation rate for Tier 1.
    pub k_inact_1: f64,
    /// Activation rate for Tier 2.
    pub k_act_2: f64,
    /// Inactivation rate for Tier 2.
    pub k_inact_2: f64,
    /// Activation rate for Tier 3.
    pub k_act_3: f64,
    /// Inactivation rate for Tier 3.
    pub k_inact_3: f64,
}

impl Default for MapkCascade {
    fn default() -> Self {
        Self {
            total_mapkkk: 100.0,
            total_mapkk: 300.0,
            total_mapk: 1000.0,
            k_act_1: 0.05,
            k_inact_1: 0.1,
            k_act_2: 0.02,
            k_inact_2: 0.05,
            k_act_3: 0.01,
            k_inact_3: 0.05,
        }
    }
}

impl MapkCascade {
    /// Compute derivatives for activated levels: `[mapkkk_act, mapkk_act, mapk_act]`.
    #[must_use]
    pub fn compute_derivatives(&self, state: &[f64; 3], input_signal: f64) -> [f64; 3] {
        let (k3_act, k2_act, k1_act) = (state[0], state[1], state[2]);

        let k3_inact = (self.total_mapkkk - k3_act).max(0.0);
        let k2_inact = (self.total_mapkk - k2_act).max(0.0);
        let k1_inact = (self.total_mapk - k1_act).max(0.0);

        // Tier 1: MAPKKK activation driven by input_signal (e.g. Ras-GTP)
        let dk3 = self.k_act_1 * input_signal * k3_inact - self.k_inact_1 * k3_act;

        // Tier 2: MAPKK activation driven by active MAPKKK
        let dk2 = self.k_act_2 * k3_act * k2_inact - self.k_inact_2 * k2_act;

        // Tier 3: MAPK activation driven by active MAPKK
        let dk1 = self.k_act_3 * k2_act * k1_inact - self.k_inact_3 * k1_act;

        [dk3, dk2, dk1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapk_signal_propagation() {
        let cascade = MapkCascade::default();
        // At rest with zero signal -> state remains at 0
        let deriv_zero = cascade.compute_derivatives(&[0.0, 0.0, 0.0], 0.0);
        assert_eq!(deriv_zero, [0.0, 0.0, 0.0]);

        // When upstream signal arrives -> Tier 1 activates first
        let deriv_active = cascade.compute_derivatives(&[0.0, 0.0, 0.0], 1.0);
        assert!(deriv_active[0] > 0.0);
    }
}
