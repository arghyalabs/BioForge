//! Semantic error types for the HIR lowering pass.

use bioforge_diagnostics::Span;

/// Errors detected during semantic analysis (AST → HIR lowering).
///
/// Each error carries a source [`Span`] for precise diagnostic rendering
/// via ariadne.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SemanticError {
    /// A referenced entity was never declared in this experiment.
    #[error("undeclared entity '{name}'")]
    UndeclaredEntity {
        /// The undeclared entity name.
        name: String,
        /// Where the reference occurred.
        span: Span,
    },

    /// A property value has the wrong physical dimension.
    ///
    /// For example, `temperature = 310 nm` (Length instead of Temperature).
    #[error("invalid unit for '{property}': expected {expected} dimension, got {got}")]
    InvalidDimension {
        /// The property name (e.g., "temperature").
        property: String,
        /// The expected dimension name.
        expected: String,
        /// The actual dimension name found.
        got: String,
        /// Where the invalid value occurred.
        span: Span,
    },

    /// pH value is outside the valid 0–14 range.
    #[error("pH value {value} is out of range (must be 0.0–14.0)")]
    PhOutOfRange {
        /// The invalid pH value.
        value: f64,
        /// Where the value occurred.
        span: Span,
    },

    /// A required property is missing from a block.
    #[error("missing required property '{property}' in {block} block")]
    MissingRequiredProperty {
        /// The block type (e.g., "simulate").
        block: String,
        /// The required property name (e.g., "timestep").
        property: String,
        /// The span of the block.
        span: Span,
    },

    /// An entity name was declared more than once.
    #[error("duplicate entity name '{name}'")]
    DuplicateEntity {
        /// The duplicated name.
        name: String,
        /// Where the first declaration occurred.
        first: Span,
        /// Where the second (duplicate) declaration occurred.
        second: Span,
    },

    /// An unknown property appeared in a validated block.
    #[error("unknown property '{property}' in {block} block")]
    UnknownProperty {
        /// The block type.
        block: String,
        /// The unknown property name.
        property: String,
        /// Where the property occurred.
        span: Span,
    },

    /// A unit string could not be resolved.
    #[error("unknown unit '{unit}'")]
    UnknownUnit {
        /// The unrecognized unit string.
        unit: String,
        /// Where the unit occurred.
        span: Span,
    },
}

impl SemanticError {
    /// Get the primary source span for this error.
    #[must_use]
    pub fn span(&self) -> Span {
        match self {
            Self::UndeclaredEntity { span, .. }
            | Self::InvalidDimension { span, .. }
            | Self::PhOutOfRange { span, .. }
            | Self::MissingRequiredProperty { span, .. }
            | Self::UnknownProperty { span, .. }
            | Self::UnknownUnit { span, .. } => *span,
            Self::DuplicateEntity { second, .. } => *second,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = SemanticError::InvalidDimension {
            property: "temperature".to_string(),
            expected: "Temperature".to_string(),
            got: "Length".to_string(),
            span: Span::new(0, 10),
        };
        assert_eq!(
            err.to_string(),
            "invalid unit for 'temperature': expected Temperature dimension, got Length"
        );
    }

    #[test]
    fn test_ph_error_display() {
        let err = SemanticError::PhOutOfRange {
            value: 15.0,
            span: Span::new(0, 4),
        };
        assert_eq!(err.to_string(), "pH value 15 is out of range (must be 0.0–14.0)");
    }

    #[test]
    fn test_span_accessor() {
        let err = SemanticError::UndeclaredEntity {
            name: "foo".to_string(),
            span: Span::new(10, 13),
        };
        assert_eq!(err.span(), Span::new(10, 13));
    }
}
