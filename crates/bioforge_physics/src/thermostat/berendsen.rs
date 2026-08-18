//! Berendsen weak-coupling thermostat for temperature regulation.

use bioforge_state::SimulationState;

use super::Thermostat;

/// Berendsen weak-coupling thermostat.
///
/// Rescales atomic velocities towards target temperature $T_0$ according to:
///
/// $$\lambda = \sqrt{1 + \frac{\Delta t}{\tau} \left(\frac{T_0}{T(t)} - 1\right)}$$
/// $$\vec{v}_i \leftarrow \lambda \vec{v}_i$$
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BerendsenThermostat {
    /// Target temperature in Kelvin ($\text{K}$).
    pub target_temperature: f64,
    /// Coupling time constant $\tau$ in picoseconds ($\text{ps}$).
    pub tau: f64,
}

impl BerendsenThermostat {
    /// Create a new Berendsen thermostat.
    ///
    /// # Arguments
    /// * `target_temperature` — Target temperature in $\text{K}$ (e.g., $310.0\text{ K}$).
    /// * `tau` — Coupling relaxation time constant in $\text{ps}$ (e.g., $0.1\text{ ps}$).
    #[must_use]
    pub fn new(target_temperature: f64, tau: f64) -> Self {
        Self {
            target_temperature: target_temperature.max(0.1),
            tau: tau.max(1e-4),
        }
    }
}

impl Thermostat for BerendsenThermostat {
    fn apply(&mut self, state: &mut SimulationState, dt: f64) {
        let current_t = state.instantaneous_temperature();
        if current_t <= 0.0 {
            return;
        }

        let ratio = self.target_temperature / current_t;
        let delta = (dt / self.tau) * (ratio - 1.0);

        // lambda^2 = 1 + (dt/tau) * (T0/T - 1)
        let lambda_sq = (1.0 + delta).max(0.01);
        let lambda = lambda_sq.sqrt();

        state.scale_velocities(lambda);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_biology::Element;

    #[test]
    fn test_berendsen_thermostat_cooling_and_heating() {
        let mut state = SimulationState::empty();
        state.num_atoms = 10;
        state.masses = vec![12.011; 10];
        state.positions = vec![[0.0, 0.0, 0.0]; 10];
        state.velocities = vec![[0.0, 0.0, 0.0]; 10];
        state.forces = vec![[0.0, 0.0, 0.0]; 10];
        state.charges = vec![0.0; 10];
        state.elements = vec![Element::from_symbol("C").unwrap(); 10];

        // Start at 100 K
        state.thermalize(100.0, 42).unwrap();
        assert!((state.instantaneous_temperature() - 100.0).abs() < 1e-6);

        // Heat to 300 K with tau = 0.1 ps, dt = 0.01 ps
        let mut thermostat = BerendsenThermostat::new(300.0, 0.1);
        for _ in 0..50 {
            thermostat.apply(&mut state, 0.01);
        }

        // Temperature should have relaxed towards 300 K
        let final_t = state.instantaneous_temperature();
        assert!(final_t > 280.0, "expected >280 K, got {}", final_t);
        assert!(final_t <= 300.1, "expected <=300.1 K, got {}", final_t);
    }
}
