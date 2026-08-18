//! Error types for multiscale coupling and scale bridges.

/// Errors that can occur during multiscale coupling, scale transformations, and nested time coordination.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MultiscaleError {
    /// Incompatible scales or missing bridge adapter.
    #[error("scale bridge error between '{from_scale}' and '{to_scale}': {message}")]
    ScaleBridgeError {
        /// Source biological scale.
        from_scale: String,
        /// Destination biological scale.
        to_scale: String,
        /// Description of failure.
        message: String,
    },

    /// An invalid thermodynamic parameter was provided (e.g. non-positive temperature).
    #[error("invalid thermodynamic parameter '{param}': {value} (must be > 0)")]
    InvalidThermodynamicParameter {
        /// Parameter name.
        param: String,
        /// Value.
        value: f64,
    },

    /// Multi-rate sub-cycling time mismatch.
    #[error("time sub-cycling mismatch: outer step ({outer_dt_s} s) is smaller than inner step ({inner_dt_s} s)")]
    InvalidSubcyclingTimestep {
        /// Outer scale timestep.
        outer_dt_s: f64,
        /// Inner scale timestep.
        inner_dt_s: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = MultiscaleError::ScaleBridgeError {
            from_scale: "Atomistic".to_string(),
            to_scale: "Reaction".to_string(),
            message: "Missing activation energy".to_string(),
        };
        assert_eq!(
            err.to_string(),
            "scale bridge error between 'Atomistic' and 'Reaction': Missing activation energy"
        );
    }
}
