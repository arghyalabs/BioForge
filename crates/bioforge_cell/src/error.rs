//! Error types for cellular and genetic modeling.

/// Errors that can occur during genetic sequence processing, transcription, and cellular simulation.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CellError {
    /// Invalid nucleotide encountered in DNA or RNA sequence.
    #[error("invalid nucleotide '{nucleotide}' at position {position} (expected A, C, G, T/U)")]
    InvalidNucleotide {
        /// Character of invalid nucleotide.
        nucleotide: char,
        /// 0-indexed position in sequence.
        position: usize,
    },

    /// Sequence length is not a multiple of 3 for codon translation.
    #[error("sequence length {length} is not a multiple of 3 (incomplete codon)")]
    IncompleteCodonSequence {
        /// Length of sequence.
        length: usize,
    },

    /// A gene was not found in the regulatory network.
    #[error("gene '{name}' not found in cellular network")]
    GeneNotFound {
        /// Gene name.
        name: String,
    },

    /// A compartment was not found in the cell.
    #[error("compartment '{name}' not found in cell")]
    CompartmentNotFound {
        /// Compartment name.
        name: String,
    },

    /// Numerical instability in cellular ODE solver.
    #[error("cellular solver instability: {0}")]
    SolverError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = CellError::InvalidNucleotide {
            nucleotide: 'X',
            position: 4,
        };
        assert_eq!(
            err.to_string(),
            "invalid nucleotide 'X' at position 4 (expected A, C, G, T/U)"
        );

        let err2 = CellError::IncompleteCodonSequence { length: 11 };
        assert_eq!(
            err2.to_string(),
            "sequence length 11 is not a multiple of 3 (incomplete codon)"
        );
    }
}
