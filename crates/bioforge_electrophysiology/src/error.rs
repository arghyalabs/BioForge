//! Error types for membrane biology and electrophysiology simulations.

/// Errors that can occur during electrophysiology modeling and action potential solving.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ElectrophysiologyError {
    /// An ion concentration was zero or negative (invalid for Nernst logarithm).
    #[error("invalid concentration for ion '{name}': inside={inside} mM, outside={outside} mM (must be > 0)")]
    InvalidIonConcentration {
        /// Ion name.
        name: String,
        /// Intracellular concentration in mM.
        inside: f64,
        /// Extracellular concentration in mM.
        outside: f64,
    },

    /// An ion was not found in the membrane configuration.
    #[error("ion '{name}' not found in membrane")]
    IonNotFound {
        /// Ion name.
        name: String,
    },

    /// Invalid membrane capacitance ($C_m \le 0$).
    #[error("invalid membrane capacitance: {value} uF/cm^2 (must be > 0)")]
    InvalidCapacitance {
        /// Capacitance value.
        value: f64,
    },

    /// Numerical instability in action potential integration.
    #[error("numerical solver instability: {message}")]
    SolverInstability {
        /// Description of instability.
        message: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = ElectrophysiologyError::InvalidIonConcentration {
            name: "Na+".to_string(),
            inside: -5.0,
            outside: 145.0,
        };
        assert_eq!(
            err.to_string(),
            "invalid concentration for ion 'Na+': inside=-5 mM, outside=145 mM (must be > 0)"
        );
    }
}
