//! Reaction Kinetics to Electrophysiology Scale Bridge (ATP-driven pumps and Ligand-gated ion channels).

use serde::{Deserialize, Serialize};

/// Ligand-gated ion channel bridge (e.g. AMPA, NMDA, GABA receptors).
///
/// Converts macroscopic biochemical neurotransmitter concentrations $[L]$ into electrical conductances:
///
/// $$P_{\text{open}}([L]) = \frac{[L]^n}{K_d^n + [L]^n}, \quad g = \bar{g} \cdot P_{\text{open}}$$
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LigandGatedChannelBridge {
    /// Maximum channel conductance in $\text{mS/cm}^2$.
    pub g_bar: f64,
    /// Channel reversal potential in $\text{mV}$.
    pub e_rev: f64,
    /// Half-maximal dissociation constant $K_d$ in Molar ($\text{M}$).
    pub kd_molar: f64,
    /// Hill cooperativity coefficient $n$.
    pub hill_n: f64,
}

impl Default for LigandGatedChannelBridge {
    fn default() -> Self {
        Self {
            g_bar: 1.0,     // 1 mS/cm^2
            e_rev: 0.0,     // 0 mV (non-selective cation channel like AMPA)
            kd_molar: 1e-6, // 1 uM
            hill_n: 2.0,
        }
    }
}

impl LigandGatedChannelBridge {
    /// Calculate instantaneous open channel probability $P_{\text{open}} \in [0.0, 1.0]$.
    #[must_use]
    pub fn open_probability(&self, ligand_conc_molar: f64) -> f64 {
        let l = ligand_conc_molar.max(0.0);
        let l_n = l.powf(self.hill_n);
        let kd_n = self.kd_molar.powf(self.hill_n);
        if l_n + kd_n > 0.0 {
            l_n / (l_n + kd_n)
        } else {
            0.0
        }
    }

    /// Calculate instantaneous conductance in $\text{mS/cm}^2$.
    #[must_use]
    pub fn conductance(&self, ligand_conc_molar: f64) -> f64 {
        self.g_bar * self.open_probability(ligand_conc_molar)
    }

    /// Calculate instantaneous transmembrane current in $\mu\text{A/cm}^2$.
    #[must_use]
    pub fn current(&self, v_membrane_mv: f64, ligand_conc_molar: f64) -> f64 {
        self.conductance(ligand_conc_molar) * (v_membrane_mv - self.e_rev)
    }
}

/// Electrogenic ATP-driven Sodium-Potassium pump ($Na^+/K^+$-ATPase) bridge.
///
/// Couples metabolic ATP hydrolysis to active ionic current extrusion ($3 Na^+$ out, $2 K^+$ in):
///
/// $$I_{\text{pump}} = I_{\max} \cdot \left(\frac{[\text{ATP}]}{K_m^{\text{ATP}} + [\text{ATP}]}\right) \cdot \left(\frac{[Na^+]_{\text{in}}}{K_m^{Na} + [Na^+]_{\text{in}}}\right)^3 \cdot \left(\frac{[K^+]_{\text{out}}}{K_m^K + [K^+]_{\text{out}}}\right)^2$$
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtpIonPumpBridge {
    /// Maximum pump current density in $\mu\text{A/cm}^2$.
    pub i_max_uA_cm2: f64,
    /// Michaelis constant for ATP in millimolar ($\text{mM}$).
    pub km_atp_mM: f64,
    /// Michaelis constant for internal Sodium in millimolar ($\text{mM}$).
    pub km_na_mM: f64,
    /// Michaelis constant for external Potassium in millimolar ($\text{mM}$).
    pub km_k_mM: f64,
}

impl Default for AtpIonPumpBridge {
    fn default() -> Self {
        Self {
            i_max_uA_cm2: 2.0,
            km_atp_mM: 0.2,
            km_na_mM: 10.0,
            km_k_mM: 1.5,
        }
    }
}

impl AtpIonPumpBridge {
    /// Compute active outward net electrogenic pump current density in $\mu\text{A/cm}^2$.
    #[must_use]
    pub fn compute_pump_current(&self, atp_mM: f64, na_in_mM: f64, k_out_mM: f64) -> f64 {
        let f_atp = atp_mM.max(0.0) / (self.km_atp_mM + atp_mM.max(0.0));
        let f_na = (na_in_mM.max(0.0) / (self.km_na_mM + na_in_mM.max(0.0))).powi(3);
        let f_k = (k_out_mM.max(0.0) / (self.km_k_mM + k_out_mM.max(0.0))).powi(2);

        self.i_max_uA_cm2 * f_atp * f_na * f_k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ligand_gated_channel_half_maximal_conductance() {
        let channel = LigandGatedChannelBridge {
            g_bar: 2.0,
            e_rev: 0.0,
            kd_molar: 1.0e-6, // 1 uM
            hill_n: 1.0,
        };

        // When [L] = Kd -> open probability is exactly 0.5
        assert_eq!(channel.open_probability(1.0e-6), 0.5);
        assert_eq!(channel.conductance(1.0e-6), 1.0); // 0.5 * 2.0

        // At -65 mV -> current is 1.0 mS/cm^2 * (-65 mV - 0) = -65 uA/cm^2 (inward current)
        let i = channel.current(-65.0, 1.0e-6);
        assert_eq!(i, -65.0);
    }

    #[test]
    fn test_atp_pump_current_saturation() {
        let pump = AtpIonPumpBridge::default();
        // At saturated ATP (10 mM), Na+ (100 mM), K+ (20 mM) -> current approaches I_max (2.0)
        let i_sat = pump.compute_pump_current(10.0, 100.0, 20.0);
        assert!(i_sat > 1.2 && i_sat < 2.0);

        // When ATP is exhausted (0 mM) -> pump ceases completely
        let i_zero = pump.compute_pump_current(0.0, 100.0, 20.0);
        assert_eq!(i_zero, 0.0);
    }
}
