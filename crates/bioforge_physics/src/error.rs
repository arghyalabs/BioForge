//! Error types for numerical physics solvers and force field calculations.

/// Errors that can occur during numerical integration and force evaluation.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum PhysicsError {
    /// Time step is invalid (zero, negative, or NaN).
    #[error("invalid integration time step dt={dt} ps (must be strictly positive and finite)")]
    InvalidTimeStep {
        /// The invalid time step.
        dt: f64,
    },

    /// Numerical explosion detected (coordinates or velocities became NaN or infinite).
    #[error("numerical instability / energy explosion at step {step}: atom {atom_index} position=[{x:.3}, {y:.3}, {z:.3}]")]
    NumericalExplosion {
        /// Integration step where instability was detected.
        step: u64,
        /// Index of the exploding atom.
        atom_index: usize,
        /// X coordinate.
        x: f64,
        /// Y coordinate.
        y: f64,
        /// Z coordinate.
        z: f64,
    },

    /// Atom mass is zero or negative, making acceleration calculation division by zero.
    #[error("atom {atom_index} has non-positive mass {mass} Da (cannot compute acceleration)")]
    NonPositiveMass {
        /// Index of the atom.
        atom_index: usize,
        /// The non-positive mass.
        mass: f64,
    },

    /// Two bonded atoms are at exactly the same coordinates (division by zero in unit vector).
    #[error("coincident bonded atoms: atom {atom1} and atom {atom2} have zero distance (cannot compute force direction)")]
    ZeroDistanceBond {
        /// Index of first atom.
        atom1: usize,
        /// Index of second atom.
        atom2: usize,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = PhysicsError::InvalidTimeStep { dt: -0.001 };
        assert_eq!(
            err.to_string(),
            "invalid integration time step dt=-0.001 ps (must be strictly positive and finite)"
        );

        let exp_err = PhysicsError::NumericalExplosion {
            step: 42,
            atom_index: 0,
            x: f64::NAN,
            y: 0.0,
            z: 0.0,
        };
        assert!(exp_err.to_string().contains("numerical instability"));
    }
}
