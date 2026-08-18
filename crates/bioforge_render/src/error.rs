//! Error types for 3D visualization and mesh generation.

/// Errors that can occur during 3D rendering and scene generation.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum RenderError {
    /// An atom index referenced in visual style was out of bounds.
    #[error("atom index {index} is out of bounds for state with {num_atoms} atoms")]
    AtomIndexOutOfBounds {
        /// Requested index.
        index: usize,
        /// Total atom count.
        num_atoms: usize,
    },

    /// Empty scene generated (zero renderable meshes).
    #[error("scene is empty (no visual geometry generated)")]
    EmptyScene,

    /// Mesh generation failed due to invalid dimensions.
    #[error("invalid geometry parameter: {parameter} = {value}")]
    InvalidGeometry {
        /// Parameter name.
        parameter: &'static str,
        /// Parameter value.
        value: f64,
    },

    /// I/O or export error.
    #[error("export error: {0}")]
    ExportError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = RenderError::InvalidGeometry {
            parameter: "radius",
            value: -1.0,
        };
        assert_eq!(err.to_string(), "invalid geometry parameter: radius = -1");
    }
}
