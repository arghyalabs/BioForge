//! Chemical species representation, physical concentrations, and discrete particle counts.

use serde::{Deserialize, Serialize};

/// Exact Avogadro Constant ($N_A = 6.02214076 \times 10^{23}\text{ mol}^{-1}$, CODATA 2018).
pub const AVOGADRO_CONSTANT: f64 = 6.022_140_76e23;

/// Default cellular compartment volume: $1.0\text{ fL} = 1.0 \times 10^{-15}\text{ L}$ (typical mammalian/bacterial cell).
pub const DEFAULT_COMPARTMENT_VOLUME_LITERS: f64 = 1.0e-15;

/// A chemical or biological species participating in a reaction network.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Species {
    /// Unique integer identifier.
    pub id: usize,
    /// Human-readable chemical name (e.g., "ATP", "Glucose-6-Phosphate", "Hexokinase").
    pub name: String,
    /// Initial physical concentration in Molar ($\text{mol/L}$).
    pub initial_concentration: f64,
    /// Reaction compartment volume in Liters ($\text{L}$).
    pub compartment_volume: f64,
}

impl Species {
    /// Create a new chemical species with initial concentration in Molar ($\text{mol/L}$).
    #[must_use]
    pub fn new(id: usize, name: impl Into<String>, initial_concentration: f64) -> Self {
        Self {
            id,
            name: name.into(),
            initial_concentration: initial_concentration.max(0.0),
            compartment_volume: DEFAULT_COMPARTMENT_VOLUME_LITERS,
        }
    }

    /// Set custom compartment volume in Liters.
    #[must_use]
    pub fn with_volume(mut self, volume_liters: f64) -> Self {
        self.compartment_volume = volume_liters.max(1e-24);
        self
    }

    /// Convert molar concentration $[C]\ (\text{M})$ to discrete particle count $N = \lfloor [C] \cdot V \cdot N_A \rceil$.
    #[must_use]
    pub fn to_discrete_count(&self, concentration_molar: f64) -> u64 {
        Self::concentration_to_count(concentration_molar, self.compartment_volume)
    }

    /// Convert discrete particle count $N$ to molar concentration $[C] = \frac{N}{V \cdot N_A}\ (\text{M})$.
    #[must_use]
    pub fn to_molar_concentration(&self, count: u64) -> f64 {
        Self::count_to_concentration(count, self.compartment_volume)
    }

    /// Convert concentration in $\text{mol/L}$ to discrete count.
    #[must_use]
    pub fn concentration_to_count(concentration_molar: f64, volume_liters: f64) -> u64 {
        let moles = concentration_molar.max(0.0) * volume_liters;
        (moles * AVOGADRO_CONSTANT).round() as u64
    }

    /// Convert discrete particle count to concentration in $\text{mol/L}$.
    #[must_use]
    pub fn count_to_concentration(count: u64, volume_liters: f64) -> f64 {
        if volume_liters <= 0.0 {
            return 0.0;
        }
        (count as f64) / (volume_liters * AVOGADRO_CONSTANT)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_concentration_and_discrete_count_conversion() {
        let species = Species::new(0, "ATP", 1.0e-3); // 1 mM

        // In 1 fL (1e-15 L), 1 mM contains:
        // N = 1e-3 * 1e-15 * 6.02214076e23 = 602,214 molecules
        let count = species.to_discrete_count(1.0e-3);
        assert_eq!(count, 602_214);

        // Convert back to concentration
        let conc = species.to_molar_concentration(count);
        assert!((conc - 1.0e-3).abs() < 1e-8);
    }
}
