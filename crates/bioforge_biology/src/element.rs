//! Chemical elements from the periodic table.
//!
//! Provides a lookup table of biologically relevant elements with
//! their atomic properties.

/// A chemical element from the periodic table.
///
/// Contains the physical properties needed for molecular simulation:
/// atomic mass (in Daltons) and van der Waals radius (in Ångströms).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    /// Chemical symbol (e.g., "C", "N", "O").
    pub symbol: &'static str,
    /// Full element name (e.g., "Carbon", "Nitrogen").
    pub name: &'static str,
    /// Atomic number (proton count).
    pub atomic_number: u8,
    /// Standard atomic mass in Daltons (Da).
    pub mass: f64,
    /// Van der Waals radius in Ångströms (Å).
    pub vdw_radius: f64,
    /// Single-bond covalent radius in Ångströms (Å) (Pyykkö 2009).
    pub covalent_radius: f64,
}

/// Built-in table of biologically relevant elements.
///
/// Masses are IUPAC 2021 standard atomic weights.
/// Van der Waals radii are from Bondi (1964) / CRC Handbook.
/// Covalent radii are from Pyykkö & Atsumi (2009).
const ELEMENTS: &[Element] = &[
    Element { symbol: "H",  name: "Hydrogen",   atomic_number: 1,  mass: 1.008,    vdw_radius: 1.20, covalent_radius: 0.31 },
    Element { symbol: "He", name: "Helium",     atomic_number: 2,  mass: 4.003,    vdw_radius: 1.40, covalent_radius: 0.28 },
    Element { symbol: "C",  name: "Carbon",     atomic_number: 6,  mass: 12.011,   vdw_radius: 1.70, covalent_radius: 0.76 },
    Element { symbol: "N",  name: "Nitrogen",   atomic_number: 7,  mass: 14.007,   vdw_radius: 1.55, covalent_radius: 0.71 },
    Element { symbol: "O",  name: "Oxygen",     atomic_number: 8,  mass: 15.999,   vdw_radius: 1.52, covalent_radius: 0.66 },
    Element { symbol: "F",  name: "Fluorine",   atomic_number: 9,  mass: 18.998,   vdw_radius: 1.47, covalent_radius: 0.57 },
    Element { symbol: "Na", name: "Sodium",     atomic_number: 11, mass: 22.990,   vdw_radius: 2.27, covalent_radius: 1.66 },
    Element { symbol: "Mg", name: "Magnesium",  atomic_number: 12, mass: 24.305,   vdw_radius: 1.73, covalent_radius: 1.41 },
    Element { symbol: "P",  name: "Phosphorus", atomic_number: 15, mass: 30.974,   vdw_radius: 1.80, covalent_radius: 1.07 },
    Element { symbol: "S",  name: "Sulfur",     atomic_number: 16, mass: 32.060,   vdw_radius: 1.80, covalent_radius: 1.05 },
    Element { symbol: "Cl", name: "Chlorine",   atomic_number: 17, mass: 35.450,   vdw_radius: 1.75, covalent_radius: 1.02 },
    Element { symbol: "K",  name: "Potassium",  atomic_number: 19, mass: 39.098,   vdw_radius: 2.75, covalent_radius: 2.03 },
    Element { symbol: "Ca", name: "Calcium",    atomic_number: 20, mass: 40.078,   vdw_radius: 2.31, covalent_radius: 1.76 },
    Element { symbol: "Mn", name: "Manganese",  atomic_number: 25, mass: 54.938,   vdw_radius: 2.05, covalent_radius: 1.39 },
    Element { symbol: "Fe", name: "Iron",       atomic_number: 26, mass: 55.845,   vdw_radius: 2.05, covalent_radius: 1.32 },
    Element { symbol: "Co", name: "Cobalt",     atomic_number: 27, mass: 58.933,   vdw_radius: 2.00, covalent_radius: 1.26 },
    Element { symbol: "Cu", name: "Copper",     atomic_number: 29, mass: 63.546,   vdw_radius: 1.40, covalent_radius: 1.32 },
    Element { symbol: "Zn", name: "Zinc",       atomic_number: 30, mass: 65.380,   vdw_radius: 1.39, covalent_radius: 1.22 },
    Element { symbol: "Se", name: "Selenium",   atomic_number: 34, mass: 78.971,   vdw_radius: 1.90, covalent_radius: 1.20 },
    Element { symbol: "Br", name: "Bromine",    atomic_number: 35, mass: 79.904,   vdw_radius: 1.85, covalent_radius: 1.20 },
    Element { symbol: "I",  name: "Iodine",     atomic_number: 53, mass: 126.904,  vdw_radius: 1.98, covalent_radius: 1.39 },
];

impl Element {
    /// Look up an element by its chemical symbol (case-sensitive).
    ///
    /// # Examples
    ///
    /// ```
    /// # use bioforge_biology::Element;
    /// let carbon = Element::from_symbol("C").unwrap();
    /// assert_eq!(carbon.atomic_number, 6);
    /// assert!((carbon.mass - 12.011).abs() < 0.001);
    /// ```
    #[must_use]
    pub fn from_symbol(symbol: &str) -> Option<Element> {
        ELEMENTS.iter().find(|e| e.symbol == symbol).copied()
    }

    /// Look up an element by atomic number.
    #[must_use]
    pub fn from_atomic_number(num: u8) -> Option<Element> {
        ELEMENTS.iter().find(|e| e.atomic_number == num).copied()
    }

    /// Get all available elements.
    #[must_use]
    pub fn all() -> &'static [Element] {
        ELEMENTS
    }
}

impl std::fmt::Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({})", self.name, self.symbol)
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_symbol_carbon() {
        let c = Element::from_symbol("C").unwrap();
        assert_eq!(c.symbol, "C");
        assert_eq!(c.name, "Carbon");
        assert_eq!(c.atomic_number, 6);
        assert!((c.mass - 12.011).abs() < 0.001);
        assert!((c.vdw_radius - 1.70).abs() < 0.01);
    }

    #[test]
    fn test_from_symbol_biologically_relevant() {
        // All key biological elements should be present
        for sym in &["H", "C", "N", "O", "S", "P", "Fe", "Ca", "Na", "K", "Cl", "Mg", "Zn"] {
            assert!(
                Element::from_symbol(sym).is_some(),
                "expected element '{}' to be in the table",
                sym
            );
        }
    }

    #[test]
    fn test_from_symbol_unknown() {
        assert!(Element::from_symbol("Xx").is_none());
        assert!(Element::from_symbol("Uuo").is_none());
    }

    #[test]
    fn test_from_atomic_number() {
        let oxygen = Element::from_atomic_number(8).unwrap();
        assert_eq!(oxygen.symbol, "O");
    }

    #[test]
    fn test_display() {
        let n = Element::from_symbol("N").unwrap();
        assert_eq!(format!("{}", n), "Nitrogen (N)");
    }

    #[test]
    fn test_hydrogen_is_lightest() {
        let h = Element::from_symbol("H").unwrap();
        assert!(h.mass < 2.0);
        assert_eq!(h.atomic_number, 1);
    }

    #[test]
    fn test_all_elements_have_positive_mass() {
        for elem in Element::all() {
            assert!(elem.mass > 0.0, "{} has non-positive mass", elem.symbol);
            assert!(elem.vdw_radius > 0.0, "{} has non-positive vdW radius", elem.symbol);
        }
    }
}
