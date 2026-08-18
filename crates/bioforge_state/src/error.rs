//! Error types for simulation state and trajectory management.

/// Errors from simulation state operations.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum StateError {
    /// The simulation state has no atoms.
    #[error("simulation state contains no atoms")]
    EmptyState,

    /// Atom index is out of bounds.
    #[error("atom index {index} is out of bounds (state contains {num_atoms} atoms)")]
    AtomIndexOutOfBounds {
        /// The invalid index requested.
        index: usize,
        /// Total number of atoms in the state.
        num_atoms: usize,
    },

    /// State arrays have mismatched lengths.
    #[error("inconsistent state dimension: {field} has length {actual}, expected {expected}")]
    InconsistentDimensions {
        /// The field name with mismatched length.
        field: &'static str,
        /// Actual length of the vector.
        actual: usize,
        /// Expected length (number of atoms).
        expected: usize,
    },

    /// Invalid temperature for thermalization.
    #[error("invalid thermalization temperature: {temp_kelvin} K (must be positive)")]
    InvalidTemperature {
        /// The invalid temperature value in Kelvin.
        temp_kelvin: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = StateError::AtomIndexOutOfBounds {
            index: 10,
            num_atoms: 5,
        };
        assert_eq!(
            err.to_string(),
            "atom index 10 is out of bounds (state contains 5 atoms)"
        );

        let dim_err = StateError::InconsistentDimensions {
            field: "velocities",
            actual: 3,
            expected: 5,
        };
        assert_eq!(
            dim_err.to_string(),
            "inconsistent state dimension: velocities has length 3, expected 5"
        );
    }
}
