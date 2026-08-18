//! Individual scientific observables (distance, RMSD, radius of gyration, energies, temperature).

use bioforge_state::SimulationState;
use std::fmt::Debug;

use crate::error::MeasurementError;

/// Trait for read-only scientific observables.
///
/// Per Scientific Principle 2, an `Observable` takes an immutable reference `&SimulationState`
/// and evaluates a scalar quantity without altering physical simulation state.
pub trait Observable: Debug + Send + Sync {
    /// Human-readable descriptor of the observable (e.g. "distance(receptor, drug)").
    fn name(&self) -> &str;

    /// Physical unit of the evaluated metric (e.g. "Å", "kJ/mol", "K").
    fn unit(&self) -> &str;

    /// Evaluate the observable against the current simulation state.
    fn evaluate(
        &self,
        state: &SimulationState,
        potential_energy: f64,
    ) -> Result<f64, MeasurementError>;
}

/// Center-of-mass distance between two atom selections: $d(A, B) = |\vec{R}_{\text{cm}}(A) - \vec{R}_{\text{cm}}(B)|$.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceObservable {
    pub name: String,
    pub group1: Vec<usize>,
    pub group2: Vec<usize>,
}

impl DistanceObservable {
    /// Create a new distance observable between single atom indices.
    #[must_use]
    pub fn pair(name: impl Into<String>, atom1: usize, atom2: usize) -> Self {
        Self {
            name: name.into(),
            group1: vec![atom1],
            group2: vec![atom2],
        }
    }

    /// Create a distance observable between two groups of atom indices.
    #[must_use]
    pub fn groups(name: impl Into<String>, group1: Vec<usize>, group2: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            group1,
            group2,
        }
    }

    fn compute_group_com(
        group: &[usize],
        state: &SimulationState,
    ) -> Result<[f64; 3], MeasurementError> {
        if group.is_empty() {
            return Err(MeasurementError::EmptySelection {
                name: "empty_group".to_string(),
            });
        }

        let mut total_m = 0.0;
        let mut com = [0.0, 0.0, 0.0];

        for &idx in group {
            if idx >= state.num_atoms {
                return Err(MeasurementError::AtomIndexOutOfBounds {
                    index: idx,
                    num_atoms: state.num_atoms,
                });
            }
            let m = state.masses[idx];
            let pos = state.positions[idx];
            total_m += m;
            com[0] += pos[0] * m;
            com[1] += pos[1] * m;
            com[2] += pos[2] * m;
        }

        if total_m <= 0.0 {
            total_m = 1.0;
        }
        Ok([com[0] / total_m, com[1] / total_m, com[2] / total_m])
    }
}

impl Observable for DistanceObservable {
    fn name(&self) -> &str {
        &self.name
    }

    fn unit(&self) -> &str {
        "Å"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        _potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        let com1 = Self::compute_group_com(&self.group1, state)?;
        let com2 = Self::compute_group_com(&self.group2, state)?;

        let dx = com1[0] - com2[0];
        let dy = com1[1] - com2[1];
        let dz = com1[2] - com2[2];

        Ok((dx * dx + dy * dy + dz * dz).sqrt())
    }
}

/// Root Mean Square Deviation (RMSD) relative to a reference coordinate frame:
///
/// $$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N |\vec{r}_i(t) - \vec{r}_i^{(0)}|^2}$$
#[derive(Debug, Clone, PartialEq)]
pub struct RmsdObservable {
    pub name: String,
    pub ref_positions: Vec<[f64; 3]>,
    pub atom_indices: Option<Vec<usize>>,
}

impl RmsdObservable {
    /// Create a new RMSD observable against reference positions.
    #[must_use]
    pub fn new(name: impl Into<String>, ref_positions: Vec<[f64; 3]>) -> Self {
        Self {
            name: name.into(),
            ref_positions,
            atom_indices: None,
        }
    }

    /// Create an RMSD observable restricted to specific atom indices.
    #[must_use]
    pub fn for_selection(
        name: impl Into<String>,
        ref_positions: Vec<[f64; 3]>,
        atom_indices: Vec<usize>,
    ) -> Self {
        Self {
            name: name.into(),
            ref_positions,
            atom_indices: Some(atom_indices),
        }
    }
}

impl Observable for RmsdObservable {
    fn name(&self) -> &str {
        &self.name
    }

    fn unit(&self) -> &str {
        "Å"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        _potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        let indices: Vec<usize> = match &self.atom_indices {
            Some(idxs) => idxs.clone(),
            None => (0..state.num_atoms).collect(),
        };

        if indices.is_empty() {
            return Ok(0.0);
        }
        if indices.len() != self.ref_positions.len() {
            return Err(MeasurementError::RmsdDimensionMismatch {
                ref_count: self.ref_positions.len(),
                state_count: indices.len(),
            });
        }

        let mut sum_sq_diff = 0.0;
        for (ref_idx, &state_idx) in indices.iter().enumerate() {
            if state_idx >= state.num_atoms {
                return Err(MeasurementError::AtomIndexOutOfBounds {
                    index: state_idx,
                    num_atoms: state.num_atoms,
                });
            }

            let pos = state.positions[state_idx];
            let ref_pos = self.ref_positions[ref_idx];

            let dx = pos[0] - ref_pos[0];
            let dy = pos[1] - ref_pos[1];
            let dz = pos[2] - ref_pos[2];

            sum_sq_diff += dx * dx + dy * dy + dz * dz;
        }

        Ok((sum_sq_diff / (indices.len() as f64)).sqrt())
    }
}

/// Radius of Gyration ($R_g$) measuring spatial compactness of a molecular structure:
///
/// $$R_g = \sqrt{\frac{\sum m_i |\vec{r}_i - \vec{R}_{\text{cm}}|^2}{\sum m_i}}$$
#[derive(Debug, Clone, PartialEq)]
pub struct RadiusOfGyrationObservable {
    pub name: String,
    pub atom_indices: Option<Vec<usize>>,
}

impl RadiusOfGyrationObservable {
    /// Create a new Radius of Gyration observable for the whole state.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            atom_indices: None,
        }
    }

    /// Create an $R_g$ observable for a subset of atoms.
    #[must_use]
    pub fn for_selection(name: impl Into<String>, atom_indices: Vec<usize>) -> Self {
        Self {
            name: name.into(),
            atom_indices: Some(atom_indices),
        }
    }
}

impl Observable for RadiusOfGyrationObservable {
    fn name(&self) -> &str {
        &self.name
    }

    fn unit(&self) -> &str {
        "Å"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        _potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        let indices: Vec<usize> = match &self.atom_indices {
            Some(idxs) => idxs.clone(),
            None => (0..state.num_atoms).collect(),
        };

        if indices.is_empty() {
            return Ok(0.0);
        }

        // Calculate center of mass for selection
        let mut total_m = 0.0;
        let mut com = [0.0, 0.0, 0.0];

        for &idx in &indices {
            if idx >= state.num_atoms {
                return Err(MeasurementError::AtomIndexOutOfBounds {
                    index: idx,
                    num_atoms: state.num_atoms,
                });
            }
            let m = state.masses[idx];
            let pos = state.positions[idx];
            total_m += m;
            com[0] += pos[0] * m;
            com[1] += pos[1] * m;
            com[2] += pos[2] * m;
        }

        if total_m <= 0.0 {
            return Ok(0.0);
        }

        com[0] /= total_m;
        com[1] /= total_m;
        com[2] /= total_m;

        let mut numerator = 0.0;
        for &idx in &indices {
            let m = state.masses[idx];
            let pos = state.positions[idx];

            let dx = pos[0] - com[0];
            let dy = pos[1] - com[1];
            let dz = pos[2] - com[2];

            numerator += m * (dx * dx + dy * dy + dz * dz);
        }

        Ok((numerator / total_m).sqrt())
    }
}

/// Instantaneous kinetic energy $K(t)$ in $\text{kJ/mol}$.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct KineticEnergyObservable;

impl Observable for KineticEnergyObservable {
    fn name(&self) -> &str {
        "kinetic_energy"
    }

    fn unit(&self) -> &str {
        "kJ/mol"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        _potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        Ok(state.kinetic_energy())
    }
}

/// Instantaneous potential energy $U(t)$ in $\text{kJ/mol}$.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PotentialEnergyObservable;

impl Observable for PotentialEnergyObservable {
    fn name(&self) -> &str {
        "potential_energy"
    }

    fn unit(&self) -> &str {
        "kJ/mol"
    }

    fn evaluate(
        &self,
        _state: &SimulationState,
        potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        Ok(potential_energy)
    }
}

/// Total physical energy $E_{\text{total}}(t) = K(t) + U(t)$ in $\text{kJ/mol}$.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TotalEnergyObservable;

impl Observable for TotalEnergyObservable {
    fn name(&self) -> &str {
        "total_energy"
    }

    fn unit(&self) -> &str {
        "kJ/mol"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        Ok(state.kinetic_energy() + potential_energy)
    }
}

/// Instantaneous kinetic temperature $T(t)$ in Kelvin ($\text{K}$).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TemperatureObservable;

impl Observable for TemperatureObservable {
    fn name(&self) -> &str {
        "temperature"
    }

    fn unit(&self) -> &str {
        "K"
    }

    fn evaluate(
        &self,
        state: &SimulationState,
        _potential_energy: f64,
    ) -> Result<f64, MeasurementError> {
        Ok(state.instantaneous_temperature())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::Element;

    #[test]
    fn test_distance_observable() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.positions = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];

        let obs = DistanceObservable::pair("C1-C2", 0, 1);
        let dist = obs.evaluate(&state, 0.0).unwrap();

        assert!((dist - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_rmsd_observable_unperturbed_and_shifted() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        let ref_pos = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        state.positions = ref_pos.clone();

        let obs = RmsdObservable::new("rmsd_test", ref_pos);

        // Unperturbed state must give RMSD = 0.0 A
        let rmsd_zero = obs.evaluate(&state, 0.0).unwrap();
        assert!(rmsd_zero.abs() < 1e-6);

        // Shift both particles by dx = 0.3 A
        state.positions = vec![[0.3, 0.0, 0.0], [1.8, 0.0, 0.0]];
        let rmsd_shifted = obs.evaluate(&state, 0.0).unwrap();
        assert!((rmsd_shifted - 0.3).abs() < 1e-6);
    }

    #[test]
    fn test_analytical_radius_of_gyration() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        // Two 12 Da carbon atoms separated by 4.0 A (at x = -2.0 A and x = +2.0 A)
        // Center of mass is origin (0,0,0)
        // Rg = sqrt( (12 * 2^2 + 12 * 2^2) / 24 ) = sqrt( (48 + 48) / 24 ) = sqrt(4) = 2.000 A
        state.masses = vec![12.011, 12.011];
        state.positions = vec![[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        state.elements = vec![
            Element::from_symbol("C").unwrap(),
            Element::from_symbol("C").unwrap(),
        ];

        let rg_obs = RadiusOfGyrationObservable::new("Rg");
        let rg = rg_obs.evaluate(&state, 0.0).unwrap();

        assert!((rg - 2.0).abs() < 1e-6, "expected Rg=2.0, got {}", rg);
    }
}
