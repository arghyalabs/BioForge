//! Lipid bilayer membrane representations, ionic gradients, Nernst and GHK equations.

use serde::{Deserialize, Serialize};

use crate::constants::{
    BODY_TEMPERATURE_KELVIN, DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2, FARADAY_CONSTANT_F,
    MOLAR_GAS_CONSTANT_R, SQUID_AXON_TEMPERATURE_KELVIN,
};
use crate::error::ElectrophysiologyError;

/// An individual chemical ion species separated by a biological membrane.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ion {
    /// Chemical symbol (e.g. "K+", "Na+", "Cl-", "Ca2+").
    pub name: String,
    /// Electrical valence $z$ (e.g. $+1$ for $\text{K}^+$, $+2$ for $\text{Ca}^{2+}$, $-1$ for $\text{Cl}^-$).
    pub valence: i32,
    /// Intracellular concentration $[\text{ion}]_{\text{in}}$ in millimolar ($\text{mM} = \text{mol/m}^3$).
    pub conc_inside_mM: f64,
    /// Extracellular concentration $[\text{ion}]_{\text{out}}$ in millimolar ($\text{mM}$).
    pub conc_outside_mM: f64,
    /// Relative membrane permeability $P_{\text{ion}}$ (normalized to $P_{\text{K}} = 1.0$).
    pub relative_permeability: f64,
}

impl Ion {
    /// Create a new ion gradient.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        valence: i32,
        inside_mM: f64,
        outside_mM: f64,
        relative_permeability: f64,
    ) -> Self {
        Self {
            name: name.into(),
            valence,
            conc_inside_mM: inside_mM,
            conc_outside_mM: outside_mM,
            relative_permeability: relative_permeability.max(0.0),
        }
    }

    /// Standard mammalian neuronal Potassium ($\text{K}^+$): $[K^+]_{\text{in}} = 140\text{ mM}, [K^+]_{\text{out}} = 5\text{ mM}, P_K = 1.0$.
    #[must_use]
    pub fn potassium_mammalian() -> Self {
        Self::new("K+", 1, 140.0, 5.0, 1.0)
    }

    /// Standard mammalian neuronal Sodium ($\text{Na}^+$): $[Na^+]_{\text{in}} = 10\text{ mM}, [Na^+]_{\text{out}} = 145\text{ mM}, P_{\text{Na}} = 0.04$.
    #[must_use]
    pub fn sodium_mammalian() -> Self {
        Self::new("Na+", 1, 10.0, 145.0, 0.04)
    }

    /// Standard mammalian neuronal Chloride ($\text{Cl}^-$): $[Cl^-]_{\text{in}} = 10\text{ mM}, [Cl^-]_{\text{out}} = 110\text{ mM}, P_{\text{Cl}} = 0.45$.
    #[must_use]
    pub fn chloride_mammalian() -> Self {
        Self::new("Cl-", -1, 10.0, 110.0, 0.45)
    }

    /// Standard mammalian neuronal Calcium ($\text{Ca}^{2+}$): $[Ca^{2+}]_{\text{in}} = 0.0001\text{ mM}\ (100\text{ nM}), [Ca^{2+}]_{\text{out}} = 2.0\text{ mM}$.
    #[must_use]
    pub fn calcium_mammalian() -> Self {
        Self::new("Ca2+", 2, 0.0001, 2.0, 0.001)
    }

    /// Compute the analytical Nernst Equilibrium Potential in millivolts ($\text{mV}$):
    ///
    /// $$E_{\text{ion}} = \frac{R T}{z F} \ln\left(\frac{[\text{ion}]_{\text{out}}}{[\text{ion}]_{\text{in}}}\right) \times 1000 \quad [\text{mV}]$$
    pub fn nernst_potential(&self, temp_k: f64) -> Result<f64, ElectrophysiologyError> {
        if self.conc_inside_mM <= 0.0 || self.conc_outside_mM <= 0.0 {
            return Err(ElectrophysiologyError::InvalidIonConcentration {
                name: self.name.clone(),
                inside: self.conc_inside_mM,
                outside: self.conc_outside_mM,
            });
        }

        let z = self.valence as f64;
        let rt_over_zf = (MOLAR_GAS_CONSTANT_R * temp_k) / (z * FARADAY_CONSTANT_F);
        let ratio = self.conc_outside_mM / self.conc_inside_mM;

        // Convert Volts to millivolts (x 1000)
        Ok(rt_over_zf * ratio.ln() * 1000.0)
    }
}

/// Biological lipid bilayer membrane separating intra- and extracellular ionic compartments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Membrane {
    /// Specific membrane capacitance in $\mu\text{F/cm}^2$ (default $1.0\,\mu\text{F/cm}^2$).
    pub capacitance_uF_cm2: f64,
    /// Absolute temperature in Kelvin ($\text{K}$).
    pub temperature_k: f64,
    /// Active ionic species gradients.
    pub ions: Vec<Ion>,
}

impl Default for Membrane {
    fn default() -> Self {
        Self::standard_mammalian()
    }
}

impl Membrane {
    /// Construct a standard mammalian neuronal membrane at $37.0^\circ\text{C}$ ($310.15\text{ K}$).
    #[must_use]
    pub fn standard_mammalian() -> Self {
        Self {
            capacitance_uF_cm2: DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2,
            temperature_k: BODY_TEMPERATURE_KELVIN,
            ions: vec![
                Ion::potassium_mammalian(),
                Ion::sodium_mammalian(),
                Ion::chloride_mammalian(),
                Ion::calcium_mammalian(),
            ],
        }
    }

    /// Construct a classic Hodgkin-Huxley squid giant axon membrane at $6.3^\circ\text{C}$ ($279.45\text{ K}$).
    #[must_use]
    pub fn standard_squid_axon() -> Self {
        Self {
            capacitance_uF_cm2: DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2,
            temperature_k: SQUID_AXON_TEMPERATURE_KELVIN,
            ions: vec![
                Ion::new("K+", 1, 400.0, 20.0, 1.0),
                Ion::new("Na+", 1, 50.0, 440.0, 0.04),
                Ion::new("Cl-", -1, 40.0, 560.0, 0.45),
            ],
        }
    }

    /// Add an ion to the membrane.
    pub fn add_ion(&mut self, ion: Ion) {
        self.ions.push(ion);
    }

    /// Find an ion by name.
    #[must_use]
    pub fn find_ion(&self, name: &str) -> Option<&Ion> {
        self.ions.iter().find(|i| i.name == name)
    }

    /// Compute Nernst potential for a specific ion in this membrane.
    pub fn nernst_potential_for(&self, name: &str) -> Result<f64, ElectrophysiologyError> {
        self.find_ion(name)
            .ok_or_else(|| ElectrophysiologyError::IonNotFound {
                name: name.to_string(),
            })?
            .nernst_potential(self.temperature_k)
    }

    /// Compute the resting membrane potential using the Goldman-Hodgkin-Katz (GHK) voltage equation:
    ///
    /// $$V_m = \frac{R T}{F} \ln\left( \frac{\sum P_{\text{cat}} [\text{cat}]_{\text{out}} + \sum P_{\text{an}} [\text{an}]_{\text{in}}}{\sum P_{\text{cat}} [\text{cat}]_{\text{in}} + \sum P_{\text{an}} [\text{an}]_{\text{out}}} \right) \times 1000 \quad [\text{mV}]$$
    pub fn ghk_resting_potential(&self) -> Result<f64, ElectrophysiologyError> {
        let mut numerator = 0.0;
        let mut denominator = 0.0;

        for ion in &self.ions {
            if ion.conc_inside_mM <= 0.0 || ion.conc_outside_mM <= 0.0 {
                return Err(ElectrophysiologyError::InvalidIonConcentration {
                    name: ion.name.clone(),
                    inside: ion.conc_inside_mM,
                    outside: ion.conc_outside_mM,
                });
            }

            let p = ion.relative_permeability;
            if ion.valence == 1 {
                // Monovalent cation (K+, Na+)
                numerator += p * ion.conc_outside_mM;
                denominator += p * ion.conc_inside_mM;
            } else if ion.valence == -1 {
                // Monovalent anion (Cl-) -> reversed inside/outside
                numerator += p * ion.conc_inside_mM;
                denominator += p * ion.conc_outside_mM;
            }
        }

        if denominator <= 0.0 || numerator <= 0.0 {
            return Err(ElectrophysiologyError::SolverInstability {
                message: "GHK terms are non-positive".to_string(),
            });
        }

        let rt_over_f = (MOLAR_GAS_CONSTANT_R * self.temperature_k) / FARADAY_CONSTANT_F;
        Ok(rt_over_f * (numerator / denominator).ln() * 1000.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Analytical Nernst Potential Benchmarks at 37°C (310.15 K).
    #[test]
    fn test_nernst_potentials_at_body_temperature() {
        let k = Ion::potassium_mammalian();
        let na = Ion::sodium_mammalian();
        let cl = Ion::chloride_mammalian();
        let ca = Ion::calcium_mammalian();

        let temp = BODY_TEMPERATURE_KELVIN;

        // RT / F at 310.15 K = (8.314462618 * 310.15) / 96485.33212 = 26.726 mV
        // E_K = 26.726 * ln(5 / 140) = 26.726 * (-3.3322) = -89.06 mV
        let e_k = k.nernst_potential(temp).unwrap();
        assert!((e_k - (-89.06)).abs() < 0.1, "expected E_K ~ -89.06 mV, got {}", e_k);

        // E_Na = 26.726 * ln(145 / 10) = 26.726 * (2.67415) = +71.47 mV
        let e_na = na.nernst_potential(temp).unwrap();
        assert!((e_na - 71.47).abs() < 0.1, "expected E_Na ~ +71.47 mV, got {}", e_na);

        // E_Cl = -26.726 * ln(110 / 10) = -26.726 * 2.39789 = -64.09 mV
        let e_cl = cl.nernst_potential(temp).unwrap();
        assert!((e_cl - (-64.09)).abs() < 0.1, "expected E_Cl ~ -64.09 mV, got {}", e_cl);

        // E_Ca = (26.726 / 2) * ln(2.0 / 0.0001) = 13.363 * ln(20000) = 13.363 * 9.90348 = +132.34 mV
        let e_ca = ca.nernst_potential(temp).unwrap();
        assert!((e_ca - 132.34).abs() < 0.5, "expected E_Ca ~ +132.34 mV, got {}", e_ca);
    }

    /// Analytical Goldman-Hodgkin-Katz (GHK) Resting Potential Benchmark.
    #[test]
    fn test_ghk_resting_membrane_potential() {
        let membrane = Membrane::standard_mammalian();
        let v_rest = membrane.ghk_resting_potential().unwrap();

        // Standard mammalian neuronal resting potential is in the range -68 mV to -72 mV
        assert!(
            v_rest < -65.0 && v_rest > -75.0,
            "expected V_rest in [-75, -65] mV, got {} mV",
            v_rest
        );
    }
}
