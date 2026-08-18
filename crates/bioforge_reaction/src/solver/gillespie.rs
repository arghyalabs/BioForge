//! Exact Stochastic Simulation Algorithm (SSA / Gillespie Direct Method).

use super::{ReactionSolver, ReactionTrajectory};
use crate::error::ReactionError;
use crate::reaction::ReactionNetwork;

/// Gillespie Direct Method Stochastic Simulation Algorithm (SSA).
///
/// Simulates exact, discrete single-molecule reaction events for master equation systems.
#[derive(Debug, Clone, PartialEq)]
pub struct GillespieSolver {
    /// Deterministic PRNG seed.
    pub seed: u64,
    /// Sampling time interval for recorded trajectory points in seconds.
    pub record_interval: f64,
    /// Maximum allowed stochastic reaction events.
    pub max_steps: u64,
}

impl Default for GillespieSolver {
    fn default() -> Self {
        Self {
            seed: 42,
            record_interval: 1.0e-2, // 10 ms
            max_steps: 1_000_000,
        }
    }
}

impl GillespieSolver {
    /// Create a new Gillespie SSA solver with a deterministic seed.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        Self {
            seed,
            ..Default::default()
        }
    }

    /// Set trajectory recording interval in seconds.
    #[must_use]
    pub fn with_record_interval(mut self, interval: f64) -> Self {
        self.record_interval = interval.max(1e-6);
        self
    }
}

impl ReactionSolver for GillespieSolver {
    fn solve(
        &mut self,
        network: &ReactionNetwork,
        total_time: f64,
    ) -> Result<ReactionTrajectory, ReactionError> {
        let num_s = network.num_species();
        let num_r = network.num_reactions();
        if num_s == 0 || num_r == 0 {
            return Ok(ReactionTrajectory::new(num_s));
        }

        let mut traj = ReactionTrajectory::new(num_s);
        let mut counts = network.initial_counts();
        let mut rng = XorShift64::new(self.seed);

        let mut t = 0.0;
        let mut next_record_t = 0.0;
        let mut steps = 0;

        // Record initial state
        let initial_concs: Vec<f64> = (0..num_s)
            .map(|i| network.species[i].to_molar_concentration(counts[i]))
            .collect();
        traj.record(0.0, &initial_concs);

        while t < total_time && steps < self.max_steps {
            steps += 1;

            // 1. Calculate propensities a_j
            let propensities = network.compute_propensities(&counts);
            let a0: f64 = propensities.iter().sum();

            if a0 <= 1e-15 {
                // No reactions possible
                break;
            }

            // 2. Draw random numbers r1, r2 in (0, 1]
            let r1 = rng.next_f64().clamp(1e-15, 1.0);
            let r2 = rng.next_f64().clamp(0.0, 1.0);

            // 3. Time to next reaction: tau = (1/a0) * ln(1/r1)
            let tau = (1.0 / a0) * (1.0 / r1).ln();

            // 4. Record any intermediate time points passed during tau
            while next_record_t <= t + tau && next_record_t <= total_time {
                if next_record_t > t {
                    let concs: Vec<f64> = (0..num_s)
                        .map(|i| network.species[i].to_molar_concentration(counts[i]))
                        .collect();
                    traj.record(next_record_t, &concs);
                }
                next_record_t += self.record_interval;
            }

            t += tau;
            if t > total_time {
                break;
            }

            // 5. Select reaction mu
            let threshold = r2 * a0;
            let mut accum = 0.0;
            let mut selected_mu = num_r - 1;

            for (r_idx, &prop) in propensities.iter().enumerate() {
                accum += prop;
                if accum >= threshold {
                    selected_mu = r_idx;
                    break;
                }
            }

            // 6. Update molecule counts according to selected reaction
            let rxn = &network.reactions[selected_mu];
            let mut valid = true;

            for &(s_idx, coeff) in &rxn.reactants {
                let required = coeff.round() as u64;
                if counts[s_idx] < required {
                    valid = false;
                    break;
                }
            }

            if valid {
                for &(s_idx, coeff) in &rxn.reactants {
                    counts[s_idx] = counts[s_idx].saturating_sub(coeff.round() as u64);
                }
                for &(s_idx, coeff) in &rxn.products {
                    counts[s_idx] += coeff.round() as u64;
                }
            }
        }

        // Final state record
        let final_concs: Vec<f64> = (0..num_s)
            .map(|i| network.species[i].to_molar_concentration(counts[i]))
            .collect();
        traj.record(total_time.min(t), &final_concs);

        Ok(traj)
    }
}

/// Deterministic 64-bit XorShift PRNG.
#[derive(Debug, Clone, PartialEq)]
struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 0x853c49e6748fea9b } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() as f64) / (u64::MAX as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ratelaw::RateLaw;

    #[test]
    fn test_gillespie_stochastic_decay() {
        let mut net = ReactionNetwork::new();
        // 1,000 molecules of A in 1 fL
        let a = net.add_species("A", 1.660538e-6); // ~1,000 molecules
        let b = net.add_species("B", 0.0);

        net.add_reaction("decay", vec![(a, 1.0)], vec![(b, 1.0)], RateLaw::mass_action_forward(1.0));

        let mut ssa = GillespieSolver::new(12345);
        let traj = ssa.solve(&net, 3.0).unwrap();

        assert!(!traj.is_empty());
        let final_a = traj.concentrations[a].last().copied().unwrap();
        let final_b = traj.concentrations[b].last().copied().unwrap();

        // A should decay significantly over 3.0 seconds with k = 1.0 s^-1
        assert!(final_a < 1.0e-6);
        assert!(final_b > 1.0e-6);
    }
}
