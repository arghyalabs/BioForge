//! Error types for reaction kinetics and network solving.

/// Errors that can occur during reaction evaluation and numerical solving.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ReactionError {
    /// A species index was not found in the reaction network.
    #[error("species index {index} is out of bounds (network contains {num_species} species)")]
    SpeciesIndexOutOfBounds {
        /// Requested species index.
        index: usize,
        /// Total number of species in the network.
        num_species: usize,
    },

    /// A species with the specified name was not found.
    #[error("species '{name}' not found in reaction network")]
    SpeciesNotFound {
        /// Species name.
        name: String,
    },

    /// Reaction with the specified name was not found.
    #[error("reaction '{name}' not found in network")]
    ReactionNotFound {
        /// Reaction name.
        name: String,
    },

    /// Concentration became negative during numerical integration.
    #[error("negative concentration encountered for species '{name}': {value} M")]
    NegativeConcentration {
        /// Species name.
        name: String,
        /// Calculated negative value.
        value: f64,
    },

    /// Numerical instability or step failure in solver.
    #[error("numerical solver error: {0}")]
    SolverError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = ReactionError::SpeciesNotFound {
            name: "ATP".to_string(),
        };
        assert_eq!(err.to_string(), "species 'ATP' not found in reaction network");

        let neg_err = ReactionError::NegativeConcentration {
            name: "Glucose".to_string(),
            value: -0.05,
        };
        assert_eq!(
            neg_err.to_string(),
            "negative concentration encountered for species 'Glucose': -0.05 M"
        );
    }
}
