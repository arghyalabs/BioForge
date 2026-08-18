//! Voltage-gated ion channels and Hodgkin-Huxley gating particle kinetics ($m, h, n$).

use serde::{Deserialize, Serialize};

/// Fast voltage-gated Sodium ($\text{Na}^+$) channel ($m^3 h$ gating).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SodiumChannel {
    /// Maximum conductance in $\text{mS/cm}^2$ (default $120.0\text{ mS/cm}^2$).
    pub g_bar: f64,
    /// Reversal potential in $\text{mV}$ (default $+50.0\text{ mV}$).
    pub e_rev: f64,
    /// Activation gate particle ($m \in [0.0, 1.0]$).
    pub m: f64,
    /// Inactivation gate particle ($h \in [0.0, 1.0]$).
    pub h: f64,
}

impl Default for SodiumChannel {
    fn default() -> Self {
        Self {
            g_bar: 120.0,
            e_rev: 50.0,
            m: alpha_m(0.0) / (alpha_m(0.0) + beta_m(0.0)),
            h: alpha_h(0.0) / (alpha_h(0.0) + beta_h(0.0)),
        }
    }
}

impl SodiumChannel {
    /// Instantaneous conductance $g_{\text{Na}} = \bar{g}_{\text{Na}} m^3 h$ in $\text{mS/cm}^2$.
    #[must_use]
    pub fn conductance(&self) -> f64 {
        self.g_bar * self.m.powi(3) * self.h
    }

    /// Instantaneous ionic current $I_{\text{Na}} = g_{\text{Na}} (V_m - E_{\text{Na}})$ in $\mu\text{A/cm}^2$.
    #[must_use]
    pub fn current(&self, v_membrane_mv: f64) -> f64 {
        self.conductance() * (v_membrane_mv - self.e_rev)
    }
}

/// Delayed rectifier Potassium ($\text{K}^+$) channel ($n^4$ gating).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PotassiumChannel {
    /// Maximum conductance in $\text{mS/cm}^2$ (default $36.0\text{ mS/cm}^2$).
    pub g_bar: f64,
    /// Reversal potential in $\text{mV}$ (default $-77.0\text{ mV}$).
    pub e_rev: f64,
    /// Activation gate particle ($n \in [0.0, 1.0]$).
    pub n: f64,
}

impl Default for PotassiumChannel {
    fn default() -> Self {
        Self {
            g_bar: 36.0,
            e_rev: -77.0,
            n: alpha_n(0.0) / (alpha_n(0.0) + beta_n(0.0)),
        }
    }
}

impl PotassiumChannel {
    /// Instantaneous conductance $g_{\text{K}} = \bar{g}_{\text{K}} n^4$ in $\text{mS/cm}^2$.
    #[must_use]
    pub fn conductance(&self) -> f64 {
        self.g_bar * self.n.powi(4)
    }

    /// Instantaneous ionic current $I_{\text{K}} = g_{\text{K}} (V_m - E_{\text{K}})$ in $\mu\text{A/cm}^2$.
    #[must_use]
    pub fn current(&self, v_membrane_mv: f64) -> f64 {
        self.conductance() * (v_membrane_mv - self.e_rev)
    }
}

/// Non-specific passive leak channel ($g_L$).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LeakChannel {
    /// Leak conductance in $\text{mS/cm}^2$ (default $0.3\text{ mS/cm}^2$).
    pub g_leak: f64,
    /// Leak reversal potential in $\text{mV}$ (default $-54.4\text{ mV}$).
    pub e_leak: f64,
}

impl Default for LeakChannel {
    fn default() -> Self {
        Self {
            g_leak: 0.3,
            e_leak: -54.4,
        }
    }
}

impl LeakChannel {
    /// Instantaneous leak current $I_L = g_L (V_m - E_L)$ in $\mu\text{A/cm}^2$.
    #[must_use]
    pub fn current(&self, v_membrane_mv: f64) -> f64 {
        self.g_leak * (v_membrane_mv - self.e_leak)
    }
}

// ─── Hodgkin-Huxley Gating Rate Functions ──────────────────────────────────────

/// Rate $\alpha_m(V)$ for sodium activation particle $m$ in $\text{ms}^{-1}$.
///
/// $V$ is membrane potential in $\text{mV}$ relative to resting potential $-65\text{ mV}$ ($V = V_m + 65$).
#[must_use]
pub fn alpha_m(v: f64) -> f64 {
    let num = 0.1 * (25.0 - v);
    let den = ((25.0 - v) / 10.0).exp() - 1.0;
    if den.abs() < 1e-7 {
        1.0 // L'Hopital limit as v -> 25
    } else {
        num / den
    }
}

/// Rate $\beta_m(V)$ for sodium activation particle $m$ in $\text{ms}^{-1}$.
#[must_use]
pub fn beta_m(v: f64) -> f64 {
    4.0 * (-v / 18.0).exp()
}

/// Rate $\alpha_h(V)$ for sodium inactivation particle $h$ in $\text{ms}^{-1}$.
#[must_use]
pub fn alpha_h(v: f64) -> f64 {
    0.07 * (-v / 20.0).exp()
}

/// Rate $\beta_h(V)$ for sodium inactivation particle $h$ in $\text{ms}^{-1}$.
#[must_use]
pub fn beta_h(v: f64) -> f64 {
    1.0 / (((30.0 - v) / 10.0).exp() + 1.0)
}

/// Rate $\alpha_n(V)$ for potassium activation particle $n$ in $\text{ms}^{-1}$.
#[must_use]
pub fn alpha_n(v: f64) -> f64 {
    let num = 0.01 * (10.0 - v);
    let den = ((10.0 - v) / 10.0).exp() - 1.0;
    if den.abs() < 1e-7 {
        0.1 // L'Hopital limit as v -> 10
    } else {
        num / den
    }
}

/// Rate $\beta_n(V)$ for potassium activation particle $n$ in $\text{ms}^{-1}$.
#[must_use]
pub fn beta_n(v: f64) -> f64 {
    0.125 * (-v / 80.0).exp()
}

/// Steady-state gating value $x_\infty = \frac{\alpha_x}{\alpha_x + \beta_x}$.
#[must_use]
pub fn steady_state_gate(alpha: f64, beta: f64) -> f64 {
    alpha / (alpha + beta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resting_steady_state_gating_values() {
        // At rest (v = 0.0 mV relative to rest -65 mV)
        let m_inf = steady_state_gate(alpha_m(0.0), beta_m(0.0));
        let h_inf = steady_state_gate(alpha_h(0.0), beta_h(0.0));
        let n_inf = steady_state_gate(alpha_n(0.0), beta_n(0.0));

        // Theoretical Hodgkin-Huxley resting values:
        // m ~ 0.053 (mostly closed)
        // h ~ 0.596 (partially available)
        // n ~ 0.318 (partially open)
        assert!((m_inf - 0.0529).abs() < 0.005, "got m_inf={}", m_inf);
        assert!((h_inf - 0.5961).abs() < 0.005, "got h_inf={}", h_inf);
        assert!((n_inf - 0.3177).abs() < 0.005, "got n_inf={}", n_inf);
    }
}
