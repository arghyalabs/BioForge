//! Atom representation.
//!
//! An [`Atom`] is a point-particle with position, mass, and charge.
//!
//! ## Assumptions (Scientific Principle 4)
//!
//! - **Point-particle approximation**: no electron cloud representation
//! - **Fixed partial charges**: no polarization model
//! - **Classical mechanics**: no quantum effects

use crate::element::Element;
use std::fmt;

/// A single atom with its physical properties.
///
/// Units (Scientific Principle 3):
/// - `position`: Ångströms (Å)
/// - `mass`: Daltons (Da)
/// - `charge`: elementary charges (e)
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Unique atom identifier (1-indexed, matching PDB convention).
    pub id: u32,
    /// The chemical element.
    pub element: Element,
    /// 3D position in Ångströms [x, y, z].
    pub position: [f64; 3],
    /// Atomic mass in Daltons. Defaults to element standard mass.
    pub mass: f64,
    /// Partial charge in elementary charges (e).
    pub charge: f64,
    /// Atom name from the structure file (e.g., "CA", "N", "CB").
    pub name: String,
    /// Residue name (e.g., "ALA", "GLY", "HOH").
    pub residue_name: Option<String>,
    /// Residue sequence number.
    pub residue_id: Option<i32>,
    /// Chain identifier (e.g., 'A', 'B').
    pub chain_id: Option<char>,
}

impl Atom {
    /// Create a new atom with minimal properties.
    ///
    /// Mass defaults to the element's standard atomic mass.
    /// Charge defaults to 0.0 (neutral).
    #[must_use]
    pub fn new(id: u32, element: Element, position: [f64; 3], name: impl Into<String>) -> Self {
        let mass = element.mass;
        Self {
            id,
            element,
            position,
            mass,
            charge: 0.0,
            name: name.into(),
            residue_name: None,
            residue_id: None,
            chain_id: None,
        }
    }

    /// Distance to another atom in Ångströms.
    #[must_use]
    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Atom({} {} [{:.3}, {:.3}, {:.3}])",
            self.id, self.element.symbol, self.position[0], self.position[1], self.position[2]
        )
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon() -> Element {
        Element::from_symbol("C").unwrap()
    }

    fn hydrogen() -> Element {
        Element::from_symbol("H").unwrap()
    }

    #[test]
    fn test_atom_construction() {
        let atom = Atom::new(1, carbon(), [1.0, 2.0, 3.0], "CA");
        assert_eq!(atom.id, 1);
        assert_eq!(atom.element.symbol, "C");
        assert_eq!(atom.name, "CA");
        assert!((atom.mass - 12.011).abs() < 0.001);
        assert!((atom.charge - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance() {
        let a = Atom::new(1, carbon(), [0.0, 0.0, 0.0], "C1");
        let b = Atom::new(2, carbon(), [3.0, 4.0, 0.0], "C2");
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_symmetric() {
        let a = Atom::new(1, hydrogen(), [1.0, 2.0, 3.0], "H1");
        let b = Atom::new(2, hydrogen(), [4.0, 6.0, 3.0], "H2");
        assert!((a.distance_to(&b) - b.distance_to(&a)).abs() < 1e-10);
    }

    #[test]
    fn test_display() {
        let atom = Atom::new(42, carbon(), [1.234, 5.678, 9.012], "CA");
        assert_eq!(format!("{}", atom), "Atom(42 C [1.234, 5.678, 9.012])");
    }
}
