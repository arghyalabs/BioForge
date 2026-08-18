//! Error types for dimensional operations.

/// Errors that can occur during dimensional operations on quantities.
///
/// These errors enforce the BioForge scientific principle that
/// "invalid dimensional operations fail at compile/runtime."
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum DimensionError {
    /// Attempted an operation on quantities with incompatible dimensions.
    ///
    /// For example, trying to add a length to a temperature:
    /// `5 nm + 310 K` → `Incompatible { op: "add", left: "Length", right: "Temperature" }`
    #[error("cannot {op} {left} and {right}: incompatible dimensions")]
    Incompatible {
        /// The operation that was attempted (e.g., "add", "subtract").
        op: String,
        /// Human-readable name of the left operand's dimension.
        left: String,
        /// Human-readable name of the right operand's dimension.
        right: String,
    },

    /// Attempted to look up a unit name that is not in the registry.
    #[error("unknown unit: '{name}'")]
    UnknownUnit {
        /// The unrecognized unit name.
        name: String,
    },

    /// Attempted to convert between units of incompatible dimensions.
    ///
    /// For example, converting nanometers to Kelvin.
    #[error("cannot convert from {from} to {to}: incompatible dimensions")]
    ConversionError {
        /// The source unit name.
        from: String,
        /// The target unit name.
        to: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_incompatible_display() {
        let err = DimensionError::Incompatible {
            op: "add".to_string(),
            left: "Length".to_string(),
            right: "Temperature".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "cannot add Length and Temperature: incompatible dimensions"
        );
    }

    #[test]
    fn test_unknown_unit_display() {
        let err = DimensionError::UnknownUnit {
            name: "furlongs".to_string(),
        };
        assert_eq!(err.to_string(), "unknown unit: 'furlongs'");
    }

    #[test]
    fn test_conversion_error_display() {
        let err = DimensionError::ConversionError {
            from: "nm".to_string(),
            to: "K".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "cannot convert from nm to K: incompatible dimensions"
        );
    }
}
