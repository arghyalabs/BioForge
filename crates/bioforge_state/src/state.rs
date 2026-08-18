//! Simulation state representation and statistical mechanics calculations.
//!
//! # Scientific Architecture (Principle 2)
//!
//! [`SimulationState`] represents the complete physical state of the simulated
//! biological system at a specific instantaneous time $t$. It is the **single source of truth**
//! for all physical observables, measurements, and numerical integration steps.
//!
//! ## Internal Units (AKMA Standard)
//!
//! - **Length**: Ångströms ($\text{Å}$)
//! - **Time**: Picoseconds ($\text{ps}$)
//! - **Mass**: Daltons ($\text{Da}$)
//! - **Energy**: Kilojoules per mole ($\text{kJ/mol}$)
//! - **Velocity**: $\text{Å/ps}$ ($1\text{ Å/ps} = 100\text{ m/s}$)
//! - **Force**: $(\text{kJ/mol})/\text{Å}$ ($1\text{ (kJ/mol)/Å} = 1.660539 \times 10^{-11}\text{ N}$)
//! - **Charge**: Elementary charge ($e$)

use bioforge_biology::{BondOrder, Element, Molecule};
use std::fmt;

use crate::error::StateError;

/// Molar gas constant in $\text{kJ}/(\text{mol}\cdot\text{K})$.
pub const MOLAR_GAS_CONSTANT_R: f64 = 8.314_462_618_153_24e-3;

/// Conversion factor from $\text{Da}\cdot(\text{Å/ps})^2$ to $\text{kJ/mol}$.
///
/// $$1\text{ Da}\cdot(\text{Å/ps})^2 = (1.6605390666 \times 10^{-27}\text{ kg}) \times (100\text{ m/s})^2 = 1.6605390666 \times 10^{-23}\text{ J}$$
/// $$\text{Per mole: } 1.6605390666 \times 10^{-23}\text{ J} \times 6.02214076 \times 10^{23}\text{ mol}^{-1} = 0.01\text{ kJ/mol}$$
pub const DA_A2_PER_PS2_TO_KJ_PER_MOL: f64 = 0.01;

/// Thermal velocity variance constant: $\frac{R \cdot 100}{1\text{ Da}} = 0.8314462618\text{ (\AA/ps)}^2\cdot\text{Da}/\text{K}$.
pub const THERMAL_VELOCITY_CONSTANT: f64 = 0.831_446_261_815_324;

/// Bond topology inside the simulation state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateBond {
    /// 0-indexed atom index of first bonded partner.
    pub atom1: usize,
    /// 0-indexed atom index of second bonded partner.
    pub atom2: usize,
    /// Equilibrium bond length $r_0$ in Ångströms.
    pub r0: f64,
    /// Harmonic spring force constant $k_b$ in $(\text{kJ/mol})/\text{Å}^2$.
    pub kb: f64,
    /// Chemical bond order.
    pub order: BondOrder,
}

/// The complete, instantaneous physical state of a simulated biological system.
///
/// Uses contiguous Structure-of-Arrays (SoA) layout for high cache locality
/// and SIMD vectorization during force calculations and coordinate updates.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulationState {
    /// Current simulation time in picoseconds ($\text{ps}$).
    pub time: f64,
    /// Integration step counter.
    pub step: u64,
    /// Number of atoms in the system.
    pub num_atoms: usize,

    // --- Dynamic Cartesian Coordinates ---
    /// Positions $[x, y, z]$ in Ångströms ($\text{Å}$).
    pub positions: Vec<[f64; 3]>,
    /// Velocities $[v_x, v_y, v_z]$ in $\text{Å/ps}$.
    pub velocities: Vec<[f64; 3]>,
    /// Net forces $[F_x, F_y, F_z]$ in $(\text{kJ/mol})/\text{Å}$.
    pub forces: Vec<[f64; 3]>,

    // --- Static Atomic Properties ---
    /// Atomic masses in Daltons ($\text{Da}$).
    pub masses: Vec<f64>,
    /// Partial electrostatic charges in elementary charges ($e$).
    pub charges: Vec<f64>,
    /// Chemical element for each atom.
    pub elements: Vec<Element>,
    /// Atom names from the structure (e.g., "CA", "N", "O").
    pub atom_names: Vec<String>,

    // --- Structural Metadata ---
    /// Residue sequence numbers.
    pub residue_ids: Vec<Option<i32>>,
    /// Residue names (e.g., "ALA", "GLY").
    pub residue_names: Vec<Option<String>>,
    /// Chain identifiers (e.g., 'A', 'B').
    pub chain_ids: Vec<Option<char>>,

    // --- Topology & Boundary Conditions ---
    /// Bond connectivity with equilibrium parameters.
    pub bonds: Vec<StateBond>,
    /// Periodic boundary box dimensions $[L_x, L_y, L_z]$ in $\text{Å}$.
    pub box_size: Option<[f64; 3]>,
}

impl SimulationState {
    /// Create a new, empty simulation state.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            time: 0.0,
            step: 0,
            num_atoms: 0,
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            charges: Vec::new(),
            elements: Vec::new(),
            atom_names: Vec::new(),
            residue_ids: Vec::new(),
            residue_names: Vec::new(),
            chain_ids: Vec::new(),
            bonds: Vec::new(),
            box_size: None,
        }
    }

    /// Construct a `SimulationState` from a biological [`Molecule`].
    ///
    /// Initializes coordinates, masses, and metadata. Velocities and forces
    /// are initialized to zero.
    #[must_use]
    pub fn from_molecule(molecule: &Molecule, box_size: Option<[f64; 3]>) -> Self {
        Self::from_molecules(&[molecule.clone()], box_size)
    }

    /// Construct a `SimulationState` from multiple biological [`Molecule`]s.
    #[must_use]
    pub fn from_molecules(molecules: &[Molecule], box_size: Option<[f64; 3]>) -> Self {
        let total_atoms: usize = molecules.iter().map(|m| m.atom_count()).sum();

        let mut positions = Vec::with_capacity(total_atoms);
        let velocities = vec![[0.0, 0.0, 0.0]; total_atoms];
        let forces = vec![[0.0, 0.0, 0.0]; total_atoms];
        let mut masses = Vec::with_capacity(total_atoms);
        let mut charges = Vec::with_capacity(total_atoms);
        let mut elements = Vec::with_capacity(total_atoms);
        let mut atom_names = Vec::with_capacity(total_atoms);
        let mut residue_ids = Vec::with_capacity(total_atoms);
        let mut residue_names = Vec::with_capacity(total_atoms);
        let mut chain_ids = Vec::with_capacity(total_atoms);
        let mut bonds = Vec::new();

        let mut atom_offset = 0;
        for mol in molecules {
            for atom in &mol.atoms {
                positions.push(atom.position);
                masses.push(atom.mass);
                charges.push(atom.charge);
                elements.push(atom.element);
                atom_names.push(atom.name.clone());
                residue_ids.push(atom.residue_id);
                residue_names.push(atom.residue_name.clone());
                chain_ids.push(atom.chain_id);
            }

            for bond in &mol.bonds {
                let id1 = bond.atom1 as usize;
                let id2 = bond.atom2 as usize;

                // Only add bonds where both atoms are in range
                if id1 > 0 && id1 <= mol.atom_count() && id2 > 0 && id2 <= mol.atom_count() {
                    let idx1 = atom_offset + (id1 - 1);
                    let idx2 = atom_offset + (id2 - 1);

                    // Compute equilibrium distance from initial positions
                    let pos1 = mol.atoms[id1 - 1].position;
                    let pos2 = mol.atoms[id2 - 1].position;
                    let dx = pos1[0] - pos2[0];
                    let dy = pos1[1] - pos2[1];
                    let dz = pos1[2] - pos2[2];
                    let r0 = (dx * dx + dy * dy + dz * dz).sqrt();

                    // Standard default spring constant kb ~ 1250.0 kJ/(mol*A^2)
                    let kb = 1250.0;

                    bonds.push(StateBond {
                        atom1: idx1,
                        atom2: idx2,
                        r0,
                        kb,
                        order: bond.order,
                    });
                }
            }

            atom_offset += mol.atom_count();
        }

        Self {
            time: 0.0,
            step: 0,
            num_atoms: total_atoms,
            positions,
            velocities,
            forces,
            masses,
            charges,
            elements,
            atom_names,
            residue_ids,
            residue_names,
            chain_ids,
            bonds,
            box_size,
        }
    }

    /// Validate internal structural consistency (all parallel vectors have equal length).
    pub fn validate(&self) -> Result<(), StateError> {
        let n = self.num_atoms;
        if self.positions.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "positions",
                actual: self.positions.len(),
                expected: n,
            });
        }
        if self.velocities.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "velocities",
                actual: self.velocities.len(),
                expected: n,
            });
        }
        if self.forces.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "forces",
                actual: self.forces.len(),
                expected: n,
            });
        }
        if self.masses.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "masses",
                actual: self.masses.len(),
                expected: n,
            });
        }
        if self.charges.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "charges",
                actual: self.charges.len(),
                expected: n,
            });
        }
        if self.elements.len() != n {
            return Err(StateError::InconsistentDimensions {
                field: "elements",
                actual: self.elements.len(),
                expected: n,
            });
        }
        Ok(())
    }

    /// Total system mass in Daltons ($\text{Da}$).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }

    /// Center of mass $\vec{R}_{\text{cm}}$ in Ångströms.
    #[must_use]
    pub fn center_of_mass(&self) -> [f64; 3] {
        if self.num_atoms == 0 {
            return [0.0, 0.0, 0.0];
        }
        let total_m = self.total_mass();
        if total_m <= 0.0 {
            return [0.0, 0.0, 0.0];
        }

        let mut com = [0.0, 0.0, 0.0];
        for i in 0..self.num_atoms {
            let m = self.masses[i];
            com[0] += self.positions[i][0] * m;
            com[1] += self.positions[i][1] * m;
            com[2] += self.positions[i][2] * m;
        }
        com[0] /= total_m;
        com[1] /= total_m;
        com[2] /= total_m;
        com
    }

    /// Total linear momentum $\vec{P} = \sum m_i \vec{v}_i$ in $\text{Da}\cdot\text{Å/ps}$.
    #[must_use]
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut p = [0.0, 0.0, 0.0];
        for i in 0..self.num_atoms {
            let m = self.masses[i];
            p[0] += m * self.velocities[i][0];
            p[1] += m * self.velocities[i][1];
            p[2] += m * self.velocities[i][2];
        }
        p
    }

    /// Center of mass velocity $\vec{v}_{\text{cm}} = \frac{\sum m_i \vec{v}_i}{\sum m_i}$ in $\text{Å/ps}$.
    #[must_use]
    pub fn center_of_mass_velocity(&self) -> [f64; 3] {
        if self.num_atoms == 0 {
            return [0.0, 0.0, 0.0];
        }
        let total_m = self.total_mass();
        if total_m <= 0.0 {
            return [0.0, 0.0, 0.0];
        }
        let p = self.total_momentum();
        [p[0] / total_m, p[1] / total_m, p[2] / total_m]
    }

    /// Remove center-of-mass translational drift so $\vec{P}_{\text{cm}} = \vec{0}$.
    pub fn remove_cm_drift(&mut self) {
        if self.num_atoms == 0 {
            return;
        }
        let v_cm = self.center_of_mass_velocity();
        for v in &mut self.velocities {
            v[0] -= v_cm[0];
            v[1] -= v_cm[1];
            v[2] -= v_cm[2];
        }
    }

    /// Total kinetic energy $K$ in $\text{kJ/mol}$.
    ///
    /// $$K = 0.01 \times \left(\frac{1}{2} \sum_{i=1}^N m_i |\vec{v}_i|^2\right)$$
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let mut raw_energy = 0.0;
        for i in 0..self.num_atoms {
            let m = self.masses[i];
            let v = self.velocities[i];
            let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            raw_energy += 0.5 * m * v_sq;
        }
        raw_energy * DA_A2_PER_PS2_TO_KJ_PER_MOL
    }

    /// Instantaneous kinetic temperature $T$ in Kelvin ($\text{K}$).
    ///
    /// Uses equipartition theorem with $3N - 3$ translational degrees of freedom:
    /// $$T = \frac{2 K}{N_{\text{dof}} R}$$
    #[must_use]
    pub fn instantaneous_temperature(&self) -> f64 {
        if self.num_atoms == 0 {
            return 0.0;
        }
        let dof = if self.num_atoms == 1 {
            3.0
        } else {
            (3 * self.num_atoms - 3) as f64
        };

        let k = self.kinetic_energy();
        (2.0 * k) / (dof * MOLAR_GAS_CONSTANT_R)
    }

    /// Thermalize the system by assigning velocities sampled from a Maxwell–Boltzmann distribution.
    ///
    /// # Algorithm
    /// 1. Sample $v_{i, \alpha} \sim \mathcal{N}\left(0, \sigma_i^2\right)$ with $\sigma_i = \sqrt{\frac{R \cdot 100 \cdot T}{m_i}}$
    ///    using Box–Muller transformation.
    /// 2. Subtract center-of-mass momentum drift ($\vec{P} = \vec{0}$).
    /// 3. Rescale velocities by $\sqrt{T_{\text{target}} / T_{\text{instantaneous}}}$ to match target temperature exactly.
    ///
    /// # Errors
    /// Returns [`StateError::InvalidTemperature`] if $T \le 0.0$, or [`StateError::EmptyState`] if $N = 0$.
    pub fn thermalize(&mut self, target_temp_kelvin: f64, seed: u64) -> Result<(), StateError> {
        if target_temp_kelvin <= 0.0 {
            return Err(StateError::InvalidTemperature {
                temp_kelvin: target_temp_kelvin,
            });
        }
        if self.num_atoms == 0 {
            return Err(StateError::EmptyState);
        }

        let mut rng = XorShift128Plus::new(seed);

        for i in 0..self.num_atoms {
            let m = self.masses[i];
            if m <= 0.0 {
                continue;
            }
            let sigma = (THERMAL_VELOCITY_CONSTANT * target_temp_kelvin / m).sqrt();

            let (z0, z1) = rng.sample_gaussian();
            let (z2, _) = rng.sample_gaussian();

            self.velocities[i] = [z0 * sigma, z1 * sigma, z2 * sigma];
        }

        // Remove center of mass velocity drift
        self.remove_cm_drift();

        // Scale to exact target temperature
        let current_temp = self.instantaneous_temperature();
        if current_temp > 0.0 {
            let scale = (target_temp_kelvin / current_temp).sqrt();
            self.scale_velocities(scale);
        }

        Ok(())
    }

    /// Multiply all atomic velocities by a scalar factor.
    pub fn scale_velocities(&mut self, factor: f64) {
        for v in &mut self.velocities {
            v[0] *= factor;
            v[1] *= factor;
            v[2] *= factor;
        }
    }

    /// Reset all force vectors to zero.
    pub fn zero_forces(&mut self) {
        for f in &mut self.forces {
            *f = [0.0, 0.0, 0.0];
        }
    }

    /// Compute direct Euclidean distance between atom `i` and atom `j` in Ångströms.
    pub fn distance(&self, i: usize, j: usize) -> Result<f64, StateError> {
        if i >= self.num_atoms {
            return Err(StateError::AtomIndexOutOfBounds {
                index: i,
                num_atoms: self.num_atoms,
            });
        }
        if j >= self.num_atoms {
            return Err(StateError::AtomIndexOutOfBounds {
                index: j,
                num_atoms: self.num_atoms,
            });
        }

        let pi = self.positions[i];
        let pj = self.positions[j];
        let dx = pi[0] - pj[0];
        let dy = pi[1] - pj[1];
        let dz = pi[2] - pj[2];
        Ok((dx * dx + dy * dy + dz * dz).sqrt())
    }

    /// Compute distance between atom `i` and atom `j` considering periodic boundary conditions (PBC).
    pub fn distance_with_pbc(&self, i: usize, j: usize) -> Result<f64, StateError> {
        if i >= self.num_atoms {
            return Err(StateError::AtomIndexOutOfBounds {
                index: i,
                num_atoms: self.num_atoms,
            });
        }
        if j >= self.num_atoms {
            return Err(StateError::AtomIndexOutOfBounds {
                index: j,
                num_atoms: self.num_atoms,
            });
        }

        let pi = self.positions[i];
        let pj = self.positions[j];
        let mut dx = pi[0] - pj[0];
        let mut dy = pi[1] - pj[1];
        let mut dz = pi[2] - pj[2];

        if let Some(box_l) = self.box_size {
            dx -= box_l[0] * (dx / box_l[0]).round();
            dy -= box_l[1] * (dy / box_l[1]).round();
            dz -= box_l[2] * (dz / box_l[2]).round();
        }

        Ok((dx * dx + dy * dy + dz * dz).sqrt())
    }
}

impl fmt::Display for SimulationState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SimulationState(t={:.3} ps, step={}, {} atoms, T={:.1} K, E_kin={:.3} kJ/mol)",
            self.time,
            self.step,
            self.num_atoms,
            self.instantaneous_temperature(),
            self.kinetic_energy()
        )
    }
}

/// Deterministic, fast XorShift128+ PRNG with Box-Muller normal sampling.
struct XorShift128Plus {
    s0: u64,
    s1: u64,
}

impl XorShift128Plus {
    fn new(seed: u64) -> Self {
        // Initialize state avoiding all-zero state
        let s0 = if seed == 0 { 0x1234_5678_9ABC_DEF0 } else { seed };
        let s1 = s0.wrapping_mul(0x5851_F42D_4C95_7F2D).wrapping_add(1);
        Self { s0, s1 }
    }

    fn next_u64(&mut self) -> u64 {
        let mut s1 = self.s0;
        let s0 = self.s1;
        self.s0 = s0;
        s1 ^= s1 << 23;
        self.s1 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
        self.s1.wrapping_add(s0)
    }

    fn next_f64(&mut self) -> f64 {
        // Generate uniform f64 in (0, 1]
        let val = (self.next_u64() >> 11) as f64;
        (val + 1.0) / ((1u64 << 53) as f64 + 1.0)
    }

    /// Generate a pair of standard normal variables (mean 0, variance 1) using Box-Muller.
    fn sample_gaussian(&mut self) -> (f64, f64) {
        let u1 = self.next_f64();
        let u2 = self.next_f64();

        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;

        (r * theta.cos(), r * theta.sin())
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::pdb::parse_pdb;

    const TINY_PDB: &str = "\
ATOM      1  N   ALA A   1       1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       1.562   2.163   1.252  1.00  0.00           C
ATOM      4  O   ALA A   1       0.735   1.685   2.056  1.00  0.00           O
ATOM      5  CB  ALA A   1       3.529   1.388   0.000  1.00  0.00           C
END
";

    #[test]
    fn test_from_molecule_pdb() {
        let mol = parse_pdb(TINY_PDB, "alanine").unwrap();
        let state = SimulationState::from_molecule(&mol, Some([50.0, 50.0, 50.0]));

        assert_eq!(state.num_atoms, 5);
        assert_eq!(state.positions.len(), 5);
        assert_eq!(state.velocities.len(), 5);
        assert_eq!(state.forces.len(), 5);
        assert_eq!(state.masses.len(), 5);
        assert!(state.validate().is_ok());

        // Check exact nitrogen position from PDB
        assert!((state.positions[0][0] - 1.458).abs() < 1e-6);
        assert_eq!(state.atom_names[0], "N");
        assert_eq!(state.elements[0].symbol, "N");

        // Total mass check: N(14.007) + C(12.011)*3 + O(15.999) = 66.039
        assert!((state.total_mass() - 66.039).abs() < 0.01);
    }

    #[test]
    fn test_kinetic_energy_exact_unit_conversion() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        // Two carbon atoms (12.011 Da each) moving at 1.0 A/ps in opposite directions
        state.masses = vec![12.011, 12.011];
        state.velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        state.positions = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        state.forces = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        state.charges = vec![0.0, 0.0];
        let c = Element::from_symbol("C").unwrap();
        state.elements = vec![c, c];

        // E_kin = 0.01 * 0.5 * (12.011 * 1.0^2 + 12.011 * (-1.0)^2) = 0.01 * 12.011 = 0.12011 kJ/mol
        let ke = state.kinetic_energy();
        assert!((ke - 0.12011).abs() < 1e-6, "got kinetic energy: {}", ke);
    }

    #[test]
    fn test_thermalization_temperature_and_momentum() {
        let mol = parse_pdb(TINY_PDB, "alanine").unwrap();
        let mut state = SimulationState::from_molecule(&mol, None);

        let target_t = 310.0; // 310 K
        state.thermalize(target_t, 42).unwrap();

        // Temperature must match target temperature exactly
        let actual_t = state.instantaneous_temperature();
        assert!(
            (actual_t - target_t).abs() < 1e-6,
            "expected {} K, got {} K",
            target_t,
            actual_t
        );

        // Net linear momentum must be zero after CM drift removal
        let p = state.total_momentum();
        assert!(p[0].abs() < 1e-10, "px non-zero: {}", p[0]);
        assert!(p[1].abs() < 1e-10, "py non-zero: {}", p[1]);
        assert!(p[2].abs() < 1e-10, "pz non-zero: {}", p[2]);
    }

    #[test]
    fn test_thermalization_reproducibility_with_seed() {
        let mol = parse_pdb(TINY_PDB, "alanine").unwrap();
        let mut state1 = SimulationState::from_molecule(&mol, None);
        let mut state2 = SimulationState::from_molecule(&mol, None);

        state1.thermalize(300.0, 12345).unwrap();
        state2.thermalize(300.0, 12345).unwrap();

        assert_eq!(state1.velocities, state2.velocities);
    }

    #[test]
    fn test_distance_with_and_without_pbc() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.positions = vec![[1.0, 0.0, 0.0], [9.0, 0.0, 0.0]];
        state.box_size = Some([10.0, 10.0, 10.0]);

        // Direct Euclidean distance: 9.0 - 1.0 = 8.0 A
        let direct_d = state.distance(0, 1).unwrap();
        assert!((direct_d - 8.0).abs() < 1e-6);

        // With PBC on 10.0 A box: shortest distance is across boundary = 2.0 A
        let pbc_d = state.distance_with_pbc(0, 1).unwrap();
        assert!((pbc_d - 2.0).abs() < 1e-6);
    }
}
