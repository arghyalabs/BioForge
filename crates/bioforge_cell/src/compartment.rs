//! Cellular organelles and trans-compartment passive & active transport.

use serde::{Deserialize, Serialize};

/// Biological organelle compartment classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OrganelleKind {
    /// Intracellular cytoplasm / cytosol.
    Cytoplasm,
    /// Cell nucleus.
    Nucleus,
    /// Mitochondrion matrix / intermembrane space.
    Mitochondria,
    /// Endoplasmic reticulum lumen.
    EndoplasmicReticulum,
    /// Extracellular microenvironment.
    Extracellular,
}

/// A physical cellular compartment with defined volume and boundary surface area.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Compartment {
    /// Unique compartment ID.
    pub id: usize,
    /// Organelle type.
    pub kind: OrganelleKind,
    /// Human-readable name (e.g. "Cytosol", "Nucleus").
    pub name: String,
    /// Volume in Liters ($\text{L}$).
    pub volume_liters: f64,
    /// Surface boundary area in square meters ($\text{m}^2$).
    pub surface_area_m2: f64,
}

impl Compartment {
    /// Standard mammalian cytosol ($V \approx 0.7\text{ fL} = 7 \times 10^{-16}\text{ L}$).
    #[must_use]
    pub fn cytoplasm(id: usize) -> Self {
        Self {
            id,
            kind: OrganelleKind::Cytoplasm,
            name: "Cytoplasm".to_string(),
            volume_liters: 7.0e-16,
            surface_area_m2: 4.0e-10,
        }
    }

    /// Standard mammalian nucleus ($V \approx 0.1\text{ fL} = 1 \times 10^{-16}\text{ L}$).
    #[must_use]
    pub fn nucleus(id: usize) -> Self {
        Self {
            id,
            kind: OrganelleKind::Nucleus,
            name: "Nucleus".to_string(),
            volume_liters: 1.0e-16,
            surface_area_m2: 1.0e-10,
        }
    }

    /// Standard mammalian mitochondrion ($V \approx 0.05\text{ fL} = 5 \times 10^{-17}\text{ L}$).
    #[must_use]
    pub fn mitochondria(id: usize) -> Self {
        Self {
            id,
            kind: OrganelleKind::Mitochondria,
            name: "Mitochondria".to_string(),
            volume_liters: 5.0e-17,
            surface_area_m2: 5.0e-11,
        }
    }
}

/// Trans-compartmental transport channel / pore between two organelles.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompartmentTransport {
    /// Source compartment index.
    pub from_compartment: usize,
    /// Destination compartment index.
    pub to_compartment: usize,
    /// Membrane permeability coefficient $P$ in meters per second ($\text{m/s}$).
    pub permeability_m_s: f64,
    /// Contact boundary surface area $A$ in square meters ($\text{m}^2$).
    pub contact_area_m2: f64,
}

impl CompartmentTransport {
    /// Calculate net transport flux in moles per second ($\text{mol/s}$) from concentration gradient:
    ///
    /// $$J = P \cdot A \cdot ([C_1] - [C_2]) \times 1000 \quad [\text{mol/s}]$$
    /// (where concentrations are in $\text{M} = \text{mol/L}$, so $\times 1000\text{ L/m}^3$).
    #[must_use]
    pub fn calculate_flux_mol_s(&self, conc_from_molar: f64, conc_to_molar: f64) -> f64 {
        let delta_c_mol_m3 = (conc_from_molar - conc_to_molar) * 1000.0;
        self.permeability_m_s * self.contact_area_m2 * delta_c_mol_m3
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compartment_creation_and_flux() {
        let cyto = Compartment::cytoplasm(0);
        let nuc = Compartment::nucleus(1);

        assert_eq!(cyto.kind, OrganelleKind::Cytoplasm);
        assert_eq!(nuc.kind, OrganelleKind::Nucleus);

        let transport = CompartmentTransport {
            from_compartment: cyto.id,
            to_compartment: nuc.id,
            permeability_m_s: 1.0e-6, // 1 um/s
            contact_area_m2: 1.0e-10, // 100 um^2
        };

        // When [C_cyto] = 1.0 mM, [C_nuc] = 0.0 mM -> positive forward flux
        let flux = transport.calculate_flux_mol_s(1.0e-3, 0.0);
        assert!(flux > 0.0);
    }
}
