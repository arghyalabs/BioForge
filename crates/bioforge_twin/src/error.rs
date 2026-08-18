//! Error types for biological digital twin orchestration and interactive simulation.

/// Errors that can occur during interactive simulation control, digital twin execution, and world synthesis.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TwinError {
    /// Requested historical checkpoint or step index was not found in ring buffer.
    #[error("checkpoint for step {step} not found in historical playback buffer")]
    CheckpointNotFound {
        /// Requested step number.
        step: usize,
    },

    /// An invalid interactive simulation command was issued.
    #[error("invalid interactive simulation command: {0}")]
    InvalidCommand(String),

    /// A digital twin layer was not found.
    #[error("scale layer '{layer_name}' not found in biological digital twin '{twin_id}'")]
    LayerNotFound {
        /// Digital twin identifier.
        twin_id: String,
        /// Missing layer.
        layer_name: String,
    },

    /// World orchestration failure.
    #[error("biological world execution error: {0}")]
    WorldError(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = TwinError::CheckpointNotFound { step: 420 };
        assert_eq!(
            err.to_string(),
            "checkpoint for step 420 not found in historical playback buffer"
        );
    }
}
