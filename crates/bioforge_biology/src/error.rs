//! Error types for biological structure operations.

/// Errors from biological structure operations.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum BiologyError {
    /// An element symbol was not recognized.
    #[error("unknown element symbol: '{symbol}'")]
    UnknownElement {
        /// The unrecognized element symbol.
        symbol: String,
    },

    /// An error occurred while parsing a PDB file.
    #[error("PDB parse error at line {line}: {message}")]
    PdbParseError {
        /// The 1-indexed line number where the error occurred.
        line: usize,
        /// A description of the parse error.
        message: String,
    },

    /// An I/O error occurred while reading a file.
    #[error("I/O error: {message}")]
    IoError {
        /// The error message.
        message: String,
    },
}

impl From<std::io::Error> for BiologyError {
    fn from(err: std::io::Error) -> Self {
        Self::IoError {
            message: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unknown_element_display() {
        let err = BiologyError::UnknownElement {
            symbol: "Xx".to_string(),
        };
        assert_eq!(err.to_string(), "unknown element symbol: 'Xx'");
    }

    #[test]
    fn test_pdb_error_display() {
        let err = BiologyError::PdbParseError {
            line: 42,
            message: "invalid coordinate".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "PDB parse error at line 42: invalid coordinate"
        );
    }
}
