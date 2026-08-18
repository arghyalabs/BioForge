//! Scientific provenance, audit trails, and cryptographic reproducibility receipts (Principles 10 & 11).

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// Scientific provenance audit trail record (W3C PROV-O compatible).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProvenanceRecord {
    /// Unique experiment ID.
    pub experiment_id: String,
    /// Simulation model name.
    pub model_name: String,
    /// Numerical solver algorithm used.
    pub solver_algorithm: String,
    /// Random number generator seed.
    pub random_seed: u64,
    /// Explicit list of all simulation parameters and initial values.
    pub parameters: Vec<(String, String)>,
    /// Git commit / build hash of the BioForge platform.
    pub build_version: String,
}

impl ProvenanceRecord {
    /// Create a new provenance audit record.
    #[must_use]
    pub fn new(
        experiment_id: impl Into<String>,
        model_name: impl Into<String>,
        solver: impl Into<String>,
        seed: u64,
    ) -> Self {
        Self {
            experiment_id: experiment_id.into(),
            model_name: model_name.into(),
            solver_algorithm: solver.into(),
            random_seed: seed,
            parameters: Vec::new(),
            build_version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }

    /// Add a key-value parameter to the provenance record.
    pub fn add_parameter(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.parameters.push((key.into(), value.into()));
    }
}

/// Cryptographic Reproducibility Receipt verifying deterministic simulation trajectories (Principle 10).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReproducibilityReceipt {
    /// Full provenance record.
    pub provenance: ProvenanceRecord,
    /// SHA-256 hash of the generated simulation trajectory.
    pub sha256_checksum: String,
}

impl ReproducibilityReceipt {
    /// Generate a cryptographic reproducibility receipt for a simulation output.
    #[must_use]
    pub fn generate(provenance: ProvenanceRecord, trajectory_bytes: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(trajectory_bytes);
        let result = hasher.finalize();
        let sha256_checksum = format!("{:x}", result);

        Self {
            provenance,
            sha256_checksum,
        }
    }

    /// Verify whether a reproduced trajectory matches the cryptographic checksum.
    #[must_use]
    pub fn verify(&self, reproduced_trajectory_bytes: &[u8]) -> bool {
        let mut hasher = Sha256::new();
        hasher.update(reproduced_trajectory_bytes);
        let result = hasher.finalize();
        let checksum = format!("{:x}", result);
        self.sha256_checksum == checksum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cryptographic_reproducibility_receipt() {
        let mut prov = ProvenanceRecord::new("EXP-101", "AlanineDynamics", "VelocityVerlet", 42);
        prov.add_parameter("timestep", "0.5 fs");
        prov.add_parameter("temperature", "310.0 K");

        let trajectory_bytes = b"Trajectory: step 0 to 1000 with exact coordinates...";
        let receipt = ReproducibilityReceipt::generate(prov, trajectory_bytes);

        // Checksum is 64 hex characters (SHA-256)
        assert_eq!(receipt.sha256_checksum.len(), 64);

        // Exact match passes verification
        assert!(receipt.verify(trajectory_bytes));

        // Corrupted / modified trajectory fails verification
        let corrupted_bytes = b"Trajectory: step 0 to 1000 with modified coordinates!";
        assert!(!receipt.verify(corrupted_bytes));
    }
}
