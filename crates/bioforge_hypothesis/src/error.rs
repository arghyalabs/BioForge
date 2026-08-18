//! Error types for hypothesis validation and counterfactual simulation.

/// Errors that can occur during hypothesis evaluation, counterfactual branching, and provenance checking.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum HypothesisError {
    /// An invalid perturbation was specified.
    #[error("invalid perturbation on target '{target}': {message}")]
    InvalidPerturbation {
        /// Target entity/gene/reaction.
        target: String,
        /// Description of failure.
        message: String,
    },

    /// Provenance hash or receipt mismatch (reproducibility failure).
    #[error("cryptographic reproducibility verification failed: expected hash {expected}, got {actual}")]
    ProvenanceMismatch {
        /// Expected SHA-256 hash.
        expected: String,
        /// Computed SHA-256 hash.
        actual: String,
    },

    /// Statistical hypothesis testing failure.
    #[error("hypothesis evaluation error: {0}")]
    EvaluationError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = HypothesisError::InvalidPerturbation {
            target: "p53".to_string(),
            message: "Target gene not found".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "invalid perturbation on target 'p53': Target gene not found"
        );
    }
}
