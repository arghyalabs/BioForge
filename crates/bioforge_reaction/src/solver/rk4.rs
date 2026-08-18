//! Deterministic 4th-order Runge-Kutta (RK4) ODE solver for reaction kinetics.

use super::{ReactionSolver, ReactionTrajectory};
use crate::error::ReactionError;
use crate::reaction::ReactionNetwork;

/// Classical 4th-order Runge-Kutta (RK4) numerical ODE integrator.
///
/// Integrates the macroscopic chemical master ODE system:
/// $$\frac{d\vec{C}}{dt} = \mathbf{N} \cdot \vec{v}(\vec{C})$$
#[derive(Debug, Clone, PartialEq)]
pub struct Rk4Solver {
    /// Numerical integration time step $\Delta t$ in seconds.
    pub dt: f64,
    /// Sampling time interval for recorded trajectory points in seconds.
    pub record_interval: f64,
}

impl Default for Rk4Solver {
    fn default() -> Self {
        Self {
            dt: 1.0e-4,             // 0.1 ms
            record_interval: 1.0e-2, // 10 ms
        }
    }
}

impl Rk4Solver {
    /// Create a new RK4 solver with specified time step $\Delta t$ in seconds.
    #[must_use]
    pub fn new(dt: f64) -> Self {
        Self {
            dt: dt.max(1e-9),
            record_interval: dt.max(1e-9) * 10.0,
        }
    }

    /// Set trajectory recording interval in seconds.
    #[must_use]
    pub fn with_record_interval(mut self, interval: f64) -> Self {
        self.record_interval = interval.max(self.dt);
        self
    }

    /// Advance concentrations by a single RK4 step of length `dt`.
    pub fn step(&self, network: &ReactionNetwork, concentrations: &mut [f64], dt: f64) {
        let n = concentrations.len();

        // k1 = f(C)
        let k1 = network.compute_derivatives(concentrations);

        // C2 = C + 0.5 * dt * k1
        let mut c2 = vec![0.0; n];
        for i in 0..n {
            c2[i] = (concentrations[i] + 0.5 * dt * k1[i]).max(0.0);
        }
        let k2 = network.compute_derivatives(&c2);

        // C3 = C + 0.5 * dt * k2
        let mut c3 = vec![0.0; n];
        for i in 0..n {
            c3[i] = (concentrations[i] + 0.5 * dt * k2[i]).max(0.0);
        }
        let k3 = network.compute_derivatives(&c3);

        // C4 = C + dt * k3
        let mut c4 = vec![0.0; n];
        for i in 0..n {
            c4[i] = (concentrations[i] + dt * k3[i]).max(0.0);
        }
        let k4 = network.compute_derivatives(&c4);

        // C_new = C + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..n {
            let delta = (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            concentrations[i] = (concentrations[i] + delta).max(0.0);
        }
    }
}

impl ReactionSolver for Rk4Solver {
    fn solve(
        &mut self,
        network: &ReactionNetwork,
        total_time: f64,
    ) -> Result<ReactionTrajectory, ReactionError> {
        let mut traj = ReactionTrajectory::new(network.num_species());
        let mut concs = network.initial_concentrations();
        let mut t = 0.0;
        let mut next_record_t = 0.0;

        // Record initial state at t=0
        traj.record(0.0, &concs);

        let dt = self.dt;
        let t_end = total_time.max(dt);

        while t < t_end - 1e-12 {
            let step_dt = dt.min(t_end - t);
            self.step(network, &mut concs, step_dt);
            t += step_dt;

            if t >= next_record_t || t >= t_end - 1e-12 {
                traj.record(t, &concs);
                next_record_t += self.record_interval;
            }
        }

        Ok(traj)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ratelaw::RateLaw;

    /// Analytical 1st-Order Exponential Decay Benchmark: A -> B with rate k.
    ///
    /// Analytical solution: [A](t) = [A]_0 * exp(-k * t)
    #[test]
    fn test_analytical_first_order_decay() {
        let mut net = ReactionNetwork::new();
        let a = net.add_species("A", 10.0e-3); // 10 mM
        let b = net.add_species("B", 0.0);

        let k = 0.5; // s^-1
        net.add_reaction("decay", vec![(a, 1.0)], vec![(b, 1.0)], RateLaw::mass_action_forward(k));

        let mut solver = Rk4Solver::new(0.001); // dt = 1 ms
        let traj = solver.solve(&net, 2.0).unwrap(); // 2.0 seconds

        // Expected at t = 2.0 s: [A](2) = 10.0e-3 * exp(-0.5 * 2.0) = 10.0e-3 * exp(-1.0) = 3.678794e-3 M
        let expected_a = 10.0e-3 * (-1.0_f64).exp();
        let final_a = traj.concentrations[a].last().copied().unwrap();

        assert!(
            (final_a - expected_a).abs() < 1e-6,
            "expected [A]={:.6e}, got {:.6e}",
            expected_a,
            final_a
        );

        // Mass conservation: [A] + [B] = 10 mM exactly
        let final_b = traj.concentrations[b].last().copied().unwrap();
        assert!((final_a + final_b - 10.0e-3).abs() < 1e-6);
    }

    /// Analytical Reversible Equilibrium Benchmark: A <-> B with k_f = 2.0, k_r = 1.0.
    ///
    /// Expected K_eq = 2.0 => [B]_eq / [A]_eq = 2.0
    #[test]
    fn test_analytical_reversible_equilibrium() {
        let mut net = ReactionNetwork::new();
        let a = net.add_species("A", 10.0e-3);
        let b = net.add_species("B", 0.0);

        net.add_reaction(
            "rev_rxn",
            vec![(a, 1.0)],
            vec![(b, 1.0)],
            RateLaw::mass_action_reversible(2.0, 1.0),
        );

        let mut solver = Rk4Solver::new(0.001);
        let traj = solver.solve(&net, 5.0).unwrap(); // 5 seconds to equilibrium

        let final_a = traj.concentrations[a].last().copied().unwrap();
        let final_b = traj.concentrations[b].last().copied().unwrap();

        // [A]_eq = 10.0 / 3.0 = 3.333 mM, [B]_eq = 20.0 / 3.0 = 6.667 mM
        let ratio = final_b / final_a;
        assert!(
            (ratio - 2.0).abs() < 1e-3,
            "expected K_eq = 2.0, got ratio = {}",
            ratio
        );
    }
}
