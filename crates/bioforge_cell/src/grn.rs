//! Gene Regulatory Networks (Promoter Hill kinetics, Toggle Switch, and Repressilator).

use serde::{Deserialize, Serialize};

/// Regulation mechanism controlling transcriptional initiation at a gene promoter.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PromoterRegulation {
    /// Constant constitutive promoter.
    Constitutive {
        /// Basal transcription rate in $\text{nM/s}$.
        basal_rate: f64,
    },

    /// Transcriptional activation by a transcription factor (activator).
    ///
    /// $$v_{\text{tx}} = v_0 + v_{\max} \frac{[A]^n}{K_A^n + [A]^n}$$
    Activation {
        basal_rate: f64,
        max_rate: f64,
        activator_idx: usize,
        k_act: f64,
        hill_coeff: f64,
    },

    /// Transcriptional repression by a repressor protein.
    ///
    /// $$v_{\text{tx}} = v_0 + v_{\max} \frac{K_R^n}{K_R^n + [R]^n}$$
    Repression {
        basal_rate: f64,
        max_rate: f64,
        repressor_idx: usize,
        k_rep: f64,
        hill_coeff: f64,
    },
}

impl PromoterRegulation {
    /// Evaluate the instantaneous transcription rate in $\text{nM/s}$ given current protein concentrations.
    #[must_use]
    pub fn evaluate_rate(&self, protein_concs: &[f64]) -> f64 {
        match self {
            PromoterRegulation::Constitutive { basal_rate } => *basal_rate,
            PromoterRegulation::Activation {
                basal_rate,
                max_rate,
                activator_idx,
                k_act,
                hill_coeff,
            } => {
                let a = if *activator_idx < protein_concs.len() {
                    protein_concs[*activator_idx].max(0.0)
                } else {
                    0.0
                };
                let a_n = a.powf(*hill_coeff);
                let k_n = k_act.powf(*hill_coeff);
                basal_rate + max_rate * (a_n / (k_n + a_n))
            }
            PromoterRegulation::Repression {
                basal_rate,
                max_rate,
                repressor_idx,
                k_rep,
                hill_coeff,
            } => {
                let r = if *repressor_idx < protein_concs.len() {
                    protein_concs[*repressor_idx].max(0.0)
                } else {
                    0.0
                };
                let r_n = r.powf(*hill_coeff);
                let k_n = k_rep.powf(*hill_coeff);
                basal_rate + max_rate * (k_n / (k_n + r_n))
            }
        }
    }
}

/// The Genetic Toggle Switch (Gardner, Cantor, Collins, Nature 2000).
///
/// A synthetic, bistable gene regulatory network composed of two mutually inhibitory repressors.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ToggleSwitch {
    /// Effective synthesis rate for repressor 1 ($\alpha_1$).
    pub alpha1: f64,
    /// Effective synthesis rate for repressor 2 ($\alpha_2$).
    pub alpha2: f64,
    /// Cooperativity exponent for repressor 2 ($\beta$).
    pub beta: f64,
    /// Cooperativity exponent for repressor 1 ($\gamma$).
    pub gamma: f64,
}

impl Default for ToggleSwitch {
    fn default() -> Self {
        Self {
            alpha1: 156.25,
            alpha2: 15.6,
            beta: 2.5,
            gamma: 1.0,
        }
    }
}

impl ToggleSwitch {
    /// Compute state derivatives $(du/dt, dv/dt)$ given current dimensionless concentrations $(u, v)$.
    #[must_use]
    pub fn compute_derivatives(&self, u: f64, v: f64) -> (f64, f64) {
        let du_dt = (self.alpha1 / (1.0 + v.max(0.0).powf(self.beta))) - u;
        let dv_dt = (self.alpha2 / (1.0 + u.max(0.0).powf(self.gamma))) - v;
        (du_dt, dv_dt)
    }
}

/// The Repressilator (Elowitz & Leibler, Nature 2000).
///
/// A synthetic 3-gene cyclic transcriptional clock generating autonomous limit-cycle oscillations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Repressilator {
    /// Transcription rate under full de-repression ($\alpha$).
    pub alpha: f64,
    /// Basal transcription rate under full repression ($\alpha_0$).
    pub alpha0: f64,
    /// Ratio of protein decay rate to mRNA decay rate ($\beta = d_p / d_m$).
    pub beta: f64,
    /// Hill coefficient ($n$).
    pub n: f64,
}

impl Default for Repressilator {
    fn default() -> Self {
        Self {
            alpha: 216.0,
            alpha0: 0.216,
            beta: 0.2,
            n: 2.0,
        }
    }
}

impl Repressilator {
    /// Compute derivatives for 6-variable state: $[m_1, p_1, m_2, p_2, m_3, p_3]$.
    #[must_use]
    pub fn compute_derivatives(&self, state: &[f64; 6]) -> [f64; 6] {
        let (m1, p1, m2, p2, m3, p3) = (state[0], state[1], state[2], state[3], state[4], state[5]);

        // dm1/dt = -m1 + alpha / (1 + p3^n) + alpha0
        let dm1 = -m1 + (self.alpha / (1.0 + p3.max(0.0).powf(self.n))) + self.alpha0;
        // dp1/dt = -beta * (p1 - m1)
        let dp1 = -self.beta * (p1 - m1);

        // dm2/dt = -m2 + alpha / (1 + p1^n) + alpha0
        let dm2 = -m2 + (self.alpha / (1.0 + p1.max(0.0).powf(self.n))) + self.alpha0;
        // dp2/dt = -beta * (p2 - m2)
        let dp2 = -self.beta * (p2 - m2);

        // dm3/dt = -m3 + alpha / (1 + p2^n) + alpha0
        let dm3 = -m3 + (self.alpha / (1.0 + p2.max(0.0).powf(self.n))) + self.alpha0;
        // dp3/dt = -beta * (p3 - m3)
        let dp3 = -self.beta * (p3 - m3);

        [dm1, dp1, dm2, dp2, dm3, dp3]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_promoter_hill_repression() {
        let promoter = PromoterRegulation::Repression {
            basal_rate: 0.01,
            max_rate: 1.0,
            repressor_idx: 0,
            k_rep: 10.0,
            hill_coeff: 2.0,
        };

        // When [R] = 0 -> rate = basal + max = 1.01
        assert_eq!(promoter.evaluate_rate(&[0.0]), 1.01);

        // When [R] = K_rep (10.0) -> rate = basal + 0.5 * max = 0.51
        assert_eq!(promoter.evaluate_rate(&[10.0]), 0.51);

        // When [R] >> K_rep -> rate -> basal (0.01)
        let high_r = promoter.evaluate_rate(&[1000.0]);
        assert!((high_r - 0.01).abs() < 1e-4);
    }

    #[test]
    fn test_toggle_switch_bistable_derivatives() {
        let toggle = ToggleSwitch::default();
        let (du, dv) = toggle.compute_derivatives(1.0, 1.0);
        // Derivatives are well-defined and finite
        assert!(du.is_finite() && dv.is_finite());
    }
}
