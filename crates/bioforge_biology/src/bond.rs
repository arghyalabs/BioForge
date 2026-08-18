//! Chemical bonds between atoms.

use std::fmt;

/// A chemical bond between two atoms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Bond {
    /// ID of the first bonded atom.
    pub atom1: u32,
    /// ID of the second bonded atom.
    pub atom2: u32,
    /// The bond order (single, double, triple, aromatic).
    pub order: BondOrder,
}

/// The order (multiplicity) of a chemical bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    /// Single bond (σ bond).
    Single,
    /// Double bond (σ + π).
    Double,
    /// Triple bond (σ + 2π).
    Triple,
    /// Aromatic bond (delocalized, ~1.5 order).
    Aromatic,
}

impl Bond {
    /// Create a new single bond between two atoms.
    #[must_use]
    pub fn single(atom1: u32, atom2: u32) -> Self {
        Self {
            atom1,
            atom2,
            order: BondOrder::Single,
        }
    }

    /// Create a new bond with a specified order.
    #[must_use]
    pub fn new(atom1: u32, atom2: u32, order: BondOrder) -> Self {
        Self {
            atom1,
            atom2,
            order,
        }
    }

    /// Check if this bond involves a given atom ID.
    #[must_use]
    pub fn involves(&self, atom_id: u32) -> bool {
        self.atom1 == atom_id || self.atom2 == atom_id
    }

    /// Get the partner atom given one atom ID. Returns `None` if
    /// the given ID is not part of this bond.
    #[must_use]
    pub fn partner(&self, atom_id: u32) -> Option<u32> {
        if self.atom1 == atom_id {
            Some(self.atom2)
        } else if self.atom2 == atom_id {
            Some(self.atom1)
        } else {
            None
        }
    }
}

impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let order_str = match self.order {
            BondOrder::Single => "-",
            BondOrder::Double => "=",
            BondOrder::Triple => "≡",
            BondOrder::Aromatic => "~",
        };
        write!(f, "{}{}{}",  self.atom1, order_str, self.atom2)
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Single => write!(f, "single"),
            Self::Double => write!(f, "double"),
            Self::Triple => write!(f, "triple"),
            Self::Aromatic => write!(f, "aromatic"),
        }
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_bond() {
        let bond = Bond::single(1, 2);
        assert_eq!(bond.atom1, 1);
        assert_eq!(bond.atom2, 2);
        assert_eq!(bond.order, BondOrder::Single);
    }

    #[test]
    fn test_bond_involves() {
        let bond = Bond::single(3, 7);
        assert!(bond.involves(3));
        assert!(bond.involves(7));
        assert!(!bond.involves(1));
    }

    #[test]
    fn test_bond_partner() {
        let bond = Bond::single(3, 7);
        assert_eq!(bond.partner(3), Some(7));
        assert_eq!(bond.partner(7), Some(3));
        assert_eq!(bond.partner(1), None);
    }

    #[test]
    fn test_bond_display() {
        assert_eq!(format!("{}", Bond::single(1, 2)), "1-2");
        assert_eq!(format!("{}", Bond::new(1, 2, BondOrder::Double)), "1=2");
        assert_eq!(format!("{}", Bond::new(1, 2, BondOrder::Triple)), "1≡2");
        assert_eq!(format!("{}", Bond::new(1, 2, BondOrder::Aromatic)), "1~2");
    }
}
