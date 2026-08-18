//! Numerical reaction kinetics solvers (Deterministic RK4 and Stochastic Gillespie SSA).

mod gillespie;
mod rk4;

pub use gillespie::GillespieSolver;
pub use rk4::Rk4Solver;

use crate::error::ReactionError;
use crate::reaction::ReactionNetwork;

/// Result trajectory produced by a reaction kinetics solver.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReactionTrajectory {
    /// Recorded simulation time points in seconds ($\text{s}$).
    pub times: Vec<f64>,
    /// Concentration values over time: `concentrations[species_idx][time_idx]` in Molar ($\text{mol/L}$).
    pub concentrations: Vec<Vec<f64>>,
}

impl ReactionTrajectory {
    /// Create a new empty trajectory for $N$ species.
    #[must_use]
    pub fn new(num_species: usize) -> Self {
        Self {
            times: Vec::new(),
            concentrations: vec![Vec::new(); num_species],
        }
    }

    /// Record a time point and concentration vector.
    pub fn record(&mut self, time: f64, concs: &[f64]) {
        self.times.push(time);
        for (i, &c) in concs.iter().enumerate() {
            if i < self.concentrations.len() {
                self.concentrations[i].push(c);
            }
        }
    }

    /// Total number of recorded time points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Whether trajectory is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// Export trajectory to multi-column CSV table.
    #[must_use]
    pub fn export_csv(&self, network: &ReactionNetwork) -> String {
        let mut out = String::from("time_s");
        for sp in &network.species {
            out.push_str(&format!(",\"[{}](M)\"", sp.name));
        }
        out.push('\n');

        let num_rows = self.times.len();
        for r in 0..num_rows {
            out.push_str(&format!("{:.6}", self.times[r]));
            for s_idx in 0..self.concentrations.len() {
                let val = if r < self.concentrations[s_idx].len() {
                    self.concentrations[s_idx][r]
                } else {
                    0.0
                };
                out.push_str(&format!(",{}", val));
            }
            out.push('\n');
        }

        out
    }
}

/// Trait for numerical reaction network solvers.
pub trait ReactionSolver: std::fmt::Debug {
    /// Solve the reaction network over `total_time` seconds.
    fn solve(
        &mut self,
        network: &ReactionNetwork,
        total_time: f64,
    ) -> Result<ReactionTrajectory, ReactionError>;
}
