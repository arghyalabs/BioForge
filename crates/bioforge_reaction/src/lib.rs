//! # BioForge Reaction (`bioforge_reaction`)
//!
//! Chemical reaction network and enzyme kinetics engine for the BioForge simulation platform.
//!
//! ## Scientific Architecture (Principle 3 & Principle 12)
//!
//! Models macroscopic chemical kinetics and biochemical pathways using dual mathematical formalisms:
//! - **Deterministic Continuous**: Ordinary Differential Equations (ODEs) integrated via Runge-Kutta 4 ([`Rk4Solver`]).
//! - **Stochastic Discrete**: Exact single-molecule master equation transitions via Gillespie SSA ([`GillespieSolver`]).
//!
//! ## Core Types
//!
//! - [`Species`]: Physical concentration $[C]\ (\text{M})$ and discrete particle count $N$.
//! - [`RateLaw`]: Mass-action, Michaelis-Menten, Hill cooperativity, and enzyme inhibition models.
//! - [`Reaction`]: Stoichiometric reaction equations.
//! - [`ReactionNetwork`]: Complete multi-reaction network with stoichiometry matrix $\mathbf{N}$.
//! - [`ReactionTrajectory`]: Time-series concentrations and CSV export.

#![deny(unsafe_code)]

pub mod error;
pub mod ratelaw;
pub mod reaction;
pub mod solver;
pub mod species;

pub use error::ReactionError;
pub use ratelaw::RateLaw;
pub use reaction::{Reaction, ReactionNetwork};
pub use solver::{GillespieSolver, ReactionSolver, ReactionTrajectory, Rk4Solver};
pub use species::{Species, AVOGADRO_CONSTANT, DEFAULT_COMPARTMENT_VOLUME_LITERS};

// ─── Elementary Enzyme Kinetics Benchmark ──────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full Elementary Michaelis-Menten Mechanism Benchmark:
    ///
    /// $$E + S \underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}} ES \xrightarrow{k_{\text{cat}}} E + P$$
    ///
    /// Validates that explicit elementary numerical ODE integration converges to the
    /// theoretical Briggs-Haldane quasi-steady state initial rate:
    ///
    /// $$K_m = \frac{k_{-1} + k_{\text{cat}}}{k_1}, \quad v_0 = \frac{k_{\text{cat}} [E]_{\text{tot}} [S]_0}{K_m + [S]_0}$$
    #[test]
    fn test_elementary_enzyme_mechanism_vs_briggs_haldane() {
        let mut net = ReactionNetwork::new();

        // Species
        let s = net.add_species("Substrate", 100.0e-6); // [S]_0 = 100 uM
        let e = net.add_species("Enzyme", 1.0e-6);      // [E]_tot = 1 uM
        let es = net.add_species("EnzymeSubstrate", 0.0);
        let p = net.add_species("Product", 0.0);

        // Kinetic Constants
        let k1 = 1.0e6;   // M^-1 s^-1 (association)
        let k_minus1 = 100.0; // s^-1 (dissociation)
        let kcat = 50.0;  // s^-1 (catalytic turnover)

        // Reaction 1: E + S <-> ES
        net.add_reaction(
            "Binding",
            vec![(e, 1.0), (s, 1.0)],
            vec![(es, 1.0)],
            RateLaw::mass_action_reversible(k1, k_minus1),
        );

        // Reaction 2: ES -> E + P
        net.add_reaction(
            "Catalysis",
            vec![(es, 1.0)],
            vec![(e, 1.0), (p, 1.0)],
            RateLaw::mass_action_forward(kcat),
        );

        // Theoretical Briggs-Haldane Km and initial rate v0
        let km = (k_minus1 + kcat) / k1; // (100 + 50) / 1e6 = 150 uM
        let vmax = kcat * 1.0e-6;        // 50 * 1 uM = 50 uM/s
        let expected_v0 = (vmax * 100.0e-6) / (km + 100.0e-6); // (50 * 100) / 250 = 20 uM/s = 2.0e-5 M/s

        let mut solver = Rk4Solver::new(1.0e-5); // dt = 10 us for fast pre-steady state resolution
        let traj = solver.solve(&net, 0.1).unwrap(); // 100 ms

        // Measure numerical initial rate of product formation: d[P]/dt
        let p_concs = &traj.concentrations[p];
        let times = &traj.times;
        let delta_p = p_concs.last().unwrap() - p_concs[times.len() / 2];
        let delta_t = times.last().unwrap() - times[times.len() / 2];
        let numerical_rate = delta_p / delta_t;

        // Verify numerical rate matches Briggs-Haldane quasi-steady state prediction within 2%
        let relative_error = (numerical_rate - expected_v0).abs() / expected_v0;
        assert!(
            relative_error < 0.02,
            "expected v0 = {:.3e} M/s, got {:.3e} M/s (rel error = {:.2}%)",
            expected_v0,
            numerical_rate,
            relative_error * 100.0
        );

        // Total Enzyme conservation: [E] + [ES] = 1 uM throughout
        let e_final = traj.concentrations[e].last().unwrap();
        let es_final = traj.concentrations[es].last().unwrap();
        assert!((e_final + es_final - 1.0e-6).abs() < 1e-9);
    }
}
