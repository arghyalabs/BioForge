//! Biochemical kinetic rate laws (Mass-Action, Michaelis-Menten, Hill, and Inhibition models).

use serde::{Deserialize, Serialize};

use crate::species::AVOGADRO_CONSTANT;

/// Rate laws defining chemical reaction velocity $v(\vec{C})$ and stochastic propensity $a(\vec{X})$.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum RateLaw {
    /// Elementary mass-action kinetics:
    /// $$v = k_f \prod [S_i]^{\nu_i} - k_r \prod [P_j]^{\nu_j}$$
    MassAction {
        /// Forward reaction rate constant $k_f$.
        k_forward: f64,
        /// Reverse reaction rate constant $k_r$ ($0.0$ for irreversible).
        k_reverse: f64,
    },

    /// Michaelis-Menten enzyme kinetics:
    /// $$v = \frac{V_{\max} [S]}{K_m + [S]}$$
    MichaelisMenten {
        /// Maximum catalytic velocity $V_{\max} = k_{\text{cat}} [E]_{\text{tot}}\ (\text{M/s})$.
        vmax: f64,
        /// Michaelis constant $K_m\ (\text{M})$.
        km: f64,
        /// Index of substrate species $[S]$.
        substrate_idx: usize,
    },

    /// Competitive enzyme inhibition:
    /// $$v = \frac{V_{\max} [S]}{K_m \left(1 + \frac{[I]}{K_i}\right) + [S]}$$
    CompetitiveInhibition {
        vmax: f64,
        km: f64,
        ki: f64,
        substrate_idx: usize,
        inhibitor_idx: usize,
    },

    /// Non-competitive enzyme inhibition:
    /// $$v = \frac{V_{\max} [S]}{\left(K_m + [S]\right)\left(1 + \frac{[I]}{K_i}\right)}$$
    NonCompetitiveInhibition {
        vmax: f64,
        km: f64,
        ki: f64,
        substrate_idx: usize,
        inhibitor_idx: usize,
    },

    /// Uncompetitive enzyme inhibition:
    /// $$v = \frac{V_{\max} [S]}{K_m + [S]\left(1 + \frac{[I]}{K_i}\right)}$$
    UncompetitiveInhibition {
        vmax: f64,
        km: f64,
        ki: f64,
        substrate_idx: usize,
        inhibitor_idx: usize,
    },

    /// Allosteric Hill cooperative kinetics:
    /// $$v = \frac{V_{\max} [S]^n}{K_{0.5}^n + [S]^n}$$
    Hill {
        vmax: f64,
        k_half: f64,
        n: f64,
        substrate_idx: usize,
    },
}

impl RateLaw {
    /// Create an irreversible mass action rate law with forward rate constant $k$.
    #[must_use]
    pub fn mass_action_forward(k: f64) -> Self {
        Self::MassAction {
            k_forward: k.max(0.0),
            k_reverse: 0.0,
        }
    }

    /// Create a reversible mass action rate law with forward and reverse constants.
    #[must_use]
    pub fn mass_action_reversible(k_forward: f64, k_reverse: f64) -> Self {
        Self::MassAction {
            k_forward: k_forward.max(0.0),
            k_reverse: k_reverse.max(0.0),
        }
    }

    /// Create a Michaelis-Menten rate law.
    #[must_use]
    pub fn michaelis_menten(vmax: f64, km: f64, substrate_idx: usize) -> Self {
        Self::MichaelisMenten {
            vmax: vmax.max(0.0),
            km: km.max(1e-12),
            substrate_idx,
        }
    }

    /// Evaluate reaction velocity $v(\vec{C})$ in $\text{M/s}$ given molar concentration vector $\vec{C}$.
    #[must_use]
    pub fn evaluate_velocity(
        &self,
        concentrations: &[f64],
        reactants: &[(usize, f64)],
        products: &[(usize, f64)],
    ) -> f64 {
        match self {
            RateLaw::MassAction {
                k_forward,
                k_reverse,
            } => {
                let mut v_fwd = *k_forward;
                for &(idx, coeff) in reactants {
                    if idx < concentrations.len() {
                        let c = concentrations[idx].max(0.0);
                        v_fwd *= c.powf(coeff);
                    }
                }

                let mut v_rev = 0.0;
                if *k_reverse > 0.0 {
                    v_rev = *k_reverse;
                    for &(idx, coeff) in products {
                        if idx < concentrations.len() {
                            let c = concentrations[idx].max(0.0);
                            v_rev *= c.powf(coeff);
                        }
                    }
                }

                v_fwd - v_rev
            }
            RateLaw::MichaelisMenten {
                vmax,
                km,
                substrate_idx,
            } => {
                if *substrate_idx >= concentrations.len() {
                    return 0.0;
                }
                let s = concentrations[*substrate_idx].max(0.0);
                (vmax * s) / (km + s)
            }
            RateLaw::CompetitiveInhibition {
                vmax,
                km,
                ki,
                substrate_idx,
                inhibitor_idx,
            } => {
                if *substrate_idx >= concentrations.len() || *inhibitor_idx >= concentrations.len() {
                    return 0.0;
                }
                let s = concentrations[*substrate_idx].max(0.0);
                let i = concentrations[*inhibitor_idx].max(0.0);
                let km_apparent = km * (1.0 + i / ki.max(1e-12));
                (vmax * s) / (km_apparent + s)
            }
            RateLaw::NonCompetitiveInhibition {
                vmax,
                km,
                ki,
                substrate_idx,
                inhibitor_idx,
            } => {
                if *substrate_idx >= concentrations.len() || *inhibitor_idx >= concentrations.len() {
                    return 0.0;
                }
                let s = concentrations[*substrate_idx].max(0.0);
                let i = concentrations[*inhibitor_idx].max(0.0);
                (vmax * s) / ((km + s) * (1.0 + i / ki.max(1e-12)))
            }
            RateLaw::UncompetitiveInhibition {
                vmax,
                km,
                ki,
                substrate_idx,
                inhibitor_idx,
            } => {
                if *substrate_idx >= concentrations.len() || *inhibitor_idx >= concentrations.len() {
                    return 0.0;
                }
                let s = concentrations[*substrate_idx].max(0.0);
                let i = concentrations[*inhibitor_idx].max(0.0);
                (vmax * s) / (km + s * (1.0 + i / ki.max(1e-12)))
            }
            RateLaw::Hill {
                vmax,
                k_half,
                n,
                substrate_idx,
            } => {
                if *substrate_idx >= concentrations.len() {
                    return 0.0;
                }
                let s = concentrations[*substrate_idx].max(0.0);
                let s_n = s.powf(*n);
                let k_n = k_half.powf(*n);
                (vmax * s_n) / (k_n + s_n)
            }
        }
    }

    /// Evaluate discrete reaction propensity $a_j(\vec{X})$ for stochastic Gillespie SSA.
    #[must_use]
    pub fn evaluate_propensity(
        &self,
        counts: &[u64],
        volume_liters: f64,
        reactants: &[(usize, f64)],
    ) -> f64 {
        match self {
            RateLaw::MassAction { k_forward, .. } => {
                let mut h = 1.0;
                let mut order = 0;

                for &(idx, coeff) in reactants {
                    if idx >= counts.len() {
                        return 0.0;
                    }
                    let n = counts[idx];
                    let nu = coeff.round() as u64;
                    order += nu;

                    if nu == 1 {
                        h *= n as f64;
                    } else if nu == 2 {
                        if n < 2 {
                            return 0.0;
                        }
                        h *= (n * (n - 1)) as f64 * 0.5;
                    } else {
                        h *= (n as f64).powf(coeff);
                    }
                }

                // Volume scaling for mesoscopic rate constant c: c = k / (V * N_A)^(order - 1)
                let vol_na = volume_liters * AVOGADRO_CONSTANT;
                let c = if order <= 1 {
                    *k_forward
                } else {
                    *k_forward / vol_na.powi((order as i32) - 1)
                };

                (c * h).max(0.0)
            }
            RateLaw::MichaelisMenten {
                vmax,
                km,
                substrate_idx,
            } => {
                if *substrate_idx >= counts.len() {
                    return 0.0;
                }
                let n_s = counts[*substrate_idx] as f64;
                let vol_na = volume_liters * AVOGADRO_CONSTANT;
                let n_km = km * vol_na;

                // Propensity = (Vmax * V * N_A * n_S) / (n_Km + n_S)
                let v_molecules_per_sec = vmax * vol_na;
                (v_molecules_per_sec * n_s) / (n_km + n_s)
            }
            _ => {
                // Fallback for complex rate laws via macroscopic conversion
                let concentrations: Vec<f64> = counts
                    .iter()
                    .map(|&c| (c as f64) / (volume_liters * AVOGADRO_CONSTANT))
                    .collect();
                let v = self.evaluate_velocity(&concentrations, reactants, &[]);
                (v * volume_liters * AVOGADRO_CONSTANT).max(0.0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_michaelis_menten_half_maximal_rate() {
        let km = 50.0e-6; // 50 uM
        let vmax = 1.0e-6; // 1 uM/s
        let mm = RateLaw::michaelis_menten(vmax, km, 0);

        // At [S] = Km, rate v must be exactly 0.5 * Vmax
        let concs = vec![km];
        let v = mm.evaluate_velocity(&concs, &[(0, 1.0)], &[]);
        assert!((v - 0.5 * vmax).abs() < 1e-12, "expected 0.5*Vmax, got {}", v);
    }

    #[test]
    fn test_competitive_inhibition_rate() {
        let km = 50.0e-6;
        let vmax = 1.0e-6;
        let ki = 10.0e-6;
        let comp = RateLaw::CompetitiveInhibition {
            vmax,
            km,
            ki,
            substrate_idx: 0,
            inhibitor_idx: 1,
        };

        // When [I] = Ki, apparent Km doubles to 2*Km => at [S] = 2*Km, v = 0.5 * Vmax
        let concs = vec![2.0 * km, ki];
        let v = comp.evaluate_velocity(&concs, &[(0, 1.0)], &[]);
        assert!((v - 0.5 * vmax).abs() < 1e-12);
    }

    #[test]
    fn test_reversible_mass_action_equilibrium_rate() {
        // A <-> B with k_f = 2.0 s^-1, k_r = 1.0 s^-1 => K_eq = 2.0
        let k_fwd = 2.0;
        let k_rev = 1.0;
        let rev = RateLaw::mass_action_reversible(k_fwd, k_rev);

        // At equilibrium [A] = 1.0 M, [B] = 2.0 M => net velocity must be zero
        let concs = vec![1.0, 2.0];
        let v_net = rev.evaluate_velocity(&concs, &[(0, 1.0)], &[(1, 1.0)]);
        assert!(v_net.abs() < 1e-12, "expected net v = 0, got {}", v_net);
    }
}
