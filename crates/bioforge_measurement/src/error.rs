//! Error types for scientific measurement evaluations.

/// Errors that can occur during observable evaluation.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MeasurementError {
    /// An atom index was out of bounds for the state.
    #[error("atom index {index} is out of bounds for state with {num_atoms} atoms")]
    AtomIndexOutOfBounds {
        /// Requested atom index.
        index: usize,
        /// Total number of atoms.
        num_atoms: usize,
    },

    /// An atom or group selection was empty.
    #[error("selection '{name}' is empty (no matching atoms found in state)")]
    EmptySelection {
        /// Selection name.
        name: String,
    },

    /// The reference structure for RMSD does not match the state particle count.
    #[error("RMSD reference coordinate count ({ref_count}) does not match state atom count ({state_count})")]
    RmsdDimensionMismatch {
        /// Number of reference coordinates.
        ref_count: usize,
        /// Number of atoms in current state.
        state_count: usize,
    },

    /// Serialization error during export.
    #[error("serialization error: {0}")]
    SerializationError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = MeasurementError::EmptySelection {
            name: "ligand".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "selection 'ligand' is empty (no matching atoms found in state)"
        );

        let rmsd_err = MeasurementError::RmsdDimensionMismatch {
            ref_count: 10,
            state_count: 5,
        };
        assert_eq!(
            rmsd_err.to_string(),
            "RMSD reference coordinate count (10) does not match state atom count (5)"
        );
    }
}
