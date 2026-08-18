//! Error types for spatial tissue modeling and PDE solving.

/// Errors that can occur during spatial discretization, diffusion stability checks, and morphogenesis solving.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TissueError {
    /// Spatial grid dimensions were invalid.
    #[error("invalid grid dimensions: {width} x {height} (must be >= 1)")]
    InvalidGridDimensions {
        /// Grid width.
        width: usize,
        /// Grid height.
        height: usize,
    },

    /// Time step violates the von Neumann / CFL diffusion stability condition ($D \Delta t / \Delta x^2 \le 0.5$).
    #[error("diffusion CFL stability violation: D*dt/dx^2 = {value} > 0.5 (reduce dt or increase dx)")]
    CflStabilityViolation {
        /// Calculated CFL parameter.
        value: f64,
    },

    /// A coordinate was out of bounds on the grid.
    #[error("spatial coordinates ({x}, {y}) out of bounds for grid of size ({width}, {height})")]
    OutOfBounds {
        /// X coordinate.
        x: usize,
        /// Y coordinate.
        y: usize,
        /// Width.
        width: usize,
        /// Height.
        height: usize,
    },

    /// Numerical solver divergence or instability.
    #[error("tissue solver divergence: {0}")]
    SolverDivergence(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = TissueError::CflStabilityViolation { value: 0.62 };
        assert_eq!(
            err.to_string(),
            "diffusion CFL stability violation: D*dt/dx^2 = 0.62 > 0.5 (reduce dt or increase dx)"
        );
    }
}
