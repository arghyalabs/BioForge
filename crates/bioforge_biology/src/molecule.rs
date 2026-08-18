//! Molecule: a collection of atoms and bonds.

use crate::atom::Atom;
use crate::bond::Bond;
use std::fmt;

/// A molecule: a collection of atoms and bonds.
///
/// This is the fundamental structural unit in BioForge.
/// Proteins, ligands, and other biological entities are all
/// represented as molecules internally.
///
/// ## Scientific Note
///
/// A `Molecule` represents static structural topology.
/// Dynamic properties (velocities, forces) are part of
/// [`SimulationState`] (Phase 5), not the molecule itself.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    /// Human-readable name (e.g., "1CRN", "aspirin").
    pub name: String,
    /// The atoms in this molecule.
    pub atoms: Vec<Atom>,
    /// The bonds in this molecule.
    pub bonds: Vec<Bond>,
}

impl Molecule {
    /// Create a new empty molecule with a name.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    /// Number of atoms in this molecule.
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds in this molecule.
    #[must_use]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Total mass in Daltons (Da).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }

    /// Center of mass in Ångströms [x, y, z].
    ///
    /// Returns `[0.0, 0.0, 0.0]` for an empty molecule.
    #[must_use]
    pub fn center_of_mass(&self) -> [f64; 3] {
        if self.atoms.is_empty() {
            return [0.0, 0.0, 0.0];
        }

        let total_mass = self.total_mass();
        if total_mass == 0.0 {
            return [0.0, 0.0, 0.0];
        }

        let mut com = [0.0, 0.0, 0.0];
        for atom in &self.atoms {
            com[0] += atom.position[0] * atom.mass;
            com[1] += atom.position[1] * atom.mass;
            com[2] += atom.position[2] * atom.mass;
        }
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
        com
    }

    /// Axis-aligned bounding box: (min_corner, max_corner) in Ångströms.
    ///
    /// Returns `None` for an empty molecule.
    #[must_use]
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.atoms.is_empty() {
            return None;
        }

        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];

        for atom in &self.atoms {
            for i in 0..3 {
                if atom.position[i] < min[i] {
                    min[i] = atom.position[i];
                }
                if atom.position[i] > max[i] {
                    max[i] = atom.position[i];
                }
            }
        }

        Some((min, max))
    }

    /// Get an atom by its ID.
    #[must_use]
    pub fn get_atom(&self, id: u32) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.id == id)
    }

    /// Get all unique chain IDs in this molecule.
    #[must_use]
    pub fn chain_ids(&self) -> Vec<char> {
        let mut ids: Vec<char> = self
            .atoms
            .iter()
            .filter_map(|a| a.chain_id)
            .collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    }

    /// Infer covalent bonds based on interatomic distances and covalent radii:
    /// $$d \le r_{\text{cov1}} + r_{\text{cov2}} + 0.45\text{ \AA} \quad (d \ge 0.4\text{ \AA})$$
    pub fn infer_bonds(&mut self) {
        let n = self.atoms.len();
        for i in 0..n {
            let p1 = self.atoms[i].position;
            let r1 = self.atoms[i].element.covalent_radius;
            let id1 = self.atoms[i].id;

            for j in (i + 1)..n {
                let p2 = self.atoms[j].position;
                let r2 = self.atoms[j].element.covalent_radius;
                let id2 = self.atoms[j].id;

                let dx = p1[0] - p2[0];
                let dy = p1[1] - p2[1];
                let dz = p1[2] - p2[2];
                let dist_sq = dx * dx + dy * dy + dz * dz;

                let max_bond_dist = r1 + r2 + 0.45;
                if dist_sq >= 0.16 && dist_sq <= max_bond_dist * max_bond_dist {
                    let exists = self.bonds.iter().any(|b| {
                        (b.atom1 == id1 && b.atom2 == id2) || (b.atom1 == id2 && b.atom2 == id1)
                    });
                    if !exists {
                        self.bonds.push(Bond::single(id1, id2));
                    }
                }
            }
        }
    }

    /// Get all unique residue names in this molecule.
    #[must_use]
    pub fn residue_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .atoms
            .iter()
            .filter_map(|a| a.residue_name.clone())
            .collect();
        names.sort();
        names.dedup();
        names
    }
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Molecule('{}', {} atoms, {} bonds, {:.1} Da)",
            self.name,
            self.atom_count(),
            self.bond_count(),
            self.total_mass()
        )
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    fn water() -> Molecule {
        let o = Element::from_symbol("O").unwrap();
        let h = Element::from_symbol("H").unwrap();

        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(1, o, [0.0, 0.0, 0.0], "O"));
        mol.atoms.push(Atom::new(2, h, [0.757, 0.586, 0.0], "H1"));
        mol.atoms.push(Atom::new(3, h, [-0.757, 0.586, 0.0], "H2"));
        mol.bonds.push(Bond::single(1, 2));
        mol.bonds.push(Bond::single(1, 3));
        mol
    }

    #[test]
    fn test_atom_count() {
        assert_eq!(water().atom_count(), 3);
    }

    #[test]
    fn test_bond_count() {
        assert_eq!(water().bond_count(), 2);
    }

    #[test]
    fn test_total_mass() {
        let mass = water().total_mass();
        // H2O: 15.999 + 2 * 1.008 = 18.015
        assert!((mass - 18.015).abs() < 0.01, "got {}", mass);
    }

    #[test]
    fn test_center_of_mass() {
        let com = water().center_of_mass();
        // Oxygen is much heavier, so COM should be near origin
        assert!(com[0].abs() < 0.1, "x: {}", com[0]);
        assert!(com[1] > 0.0, "y should be positive: {}", com[1]);
        assert!(com[2].abs() < 1e-10, "z: {}", com[2]);
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new("empty");
        assert_eq!(mol.atom_count(), 0);
        assert_eq!(mol.total_mass(), 0.0);
        assert_eq!(mol.center_of_mass(), [0.0, 0.0, 0.0]);
        assert!(mol.bounding_box().is_none());
    }

    #[test]
    fn test_bounding_box() {
        let (min, max) = water().bounding_box().unwrap();
        assert!(min[0] < 0.0); // H2 is at -0.757
        assert!(max[0] > 0.0); // H1 is at 0.757
    }

    #[test]
    fn test_get_atom() {
        let mol = water();
        assert_eq!(mol.get_atom(1).unwrap().element.symbol, "O");
        assert_eq!(mol.get_atom(2).unwrap().element.symbol, "H");
        assert!(mol.get_atom(99).is_none());
    }

    #[test]
    fn test_display() {
        let display = format!("{}", water());
        assert!(display.contains("water"));
        assert!(display.contains("3 atoms"));
        assert!(display.contains("2 bonds"));
    }
}
