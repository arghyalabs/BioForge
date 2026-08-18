//! Sub-cycling multi-rate time-stepping coordinator across nested biological scales.

use crate::error::MultiscaleError;

/// Multiscale time-stepping coordinator managing hierarchical sub-cycling.
///
/// Enables coarse-scale processes (gene regulation, tissue diffusion $\Delta t \sim 1\text{ s}$)
/// to coordinate fine-scale processes (reaction kinetics, electrophysiology $\Delta t \sim 1\text{ ms}$).
#[derive(Debug, Clone, PartialEq)]
pub struct MultiscaleCoordinator {
    /// Current global simulation time in seconds ($\text{s}$).
    pub current_time_s: f64,
    /// Total target duration in seconds ($\text{s}$).
    pub total_duration_s: f64,
    /// Macroscopic outer time step $\Delta t_{\text{outer}}$ in seconds ($\text{s}$).
    pub outer_dt_s: f64,
    /// Number of microscopic sub-cycling steps per outer step.
    pub inner_subcycles: usize,
}

impl MultiscaleCoordinator {
    /// Create a new multiscale coordinator.
    pub fn new(
        total_duration_s: f64,
        outer_dt_s: f64,
        inner_subcycles: usize,
    ) -> Result<Self, MultiscaleError> {
        if outer_dt_s <= 0.0 || total_duration_s <= 0.0 {
            return Err(MultiscaleError::InvalidThermodynamicParameter {
                param: "outer_dt_s".to_string(),
                value: outer_dt_s,
            });
        }
        if inner_subcycles == 0 {
            return Err(MultiscaleError::InvalidSubcyclingTimestep {
                outer_dt_s,
                inner_dt_s: 0.0,
            });
        }

        Ok(Self {
            current_time_s: 0.0,
            total_duration_s,
            outer_dt_s,
            inner_subcycles,
        })
    }

    /// Microscopic inner sub-cycling time step $\Delta t_{\text{inner}} = \frac{\Delta t_{\text{outer}}}{\text{subcycles}}$ in seconds.
    #[must_use]
    pub fn inner_dt_s(&self) -> f64 {
        self.outer_dt_s / (self.inner_subcycles as f64)
    }

    /// Advance one macroscopic outer step, executing `inner_subcycles` microscopic steps.
    pub fn step<FInner, FOuter>(
        &mut self,
        mut inner_step_fn: FInner,
        mut outer_step_fn: FOuter,
    ) -> Result<(), MultiscaleError>
    where
        FInner: FnMut(f64, f64) -> Result<(), MultiscaleError>,
        FOuter: FnMut(f64, f64) -> Result<(), MultiscaleError>,
    {
        let inner_dt = self.inner_dt_s();
        let mut t = self.current_time_s;

        // Execute inner microscopic sub-cycles
        for _ in 0..self.inner_subcycles {
            inner_step_fn(t, inner_dt)?;
            t += inner_dt;
        }

        // Execute outer macroscopic step
        outer_step_fn(self.current_time_s, self.outer_dt_s)?;
        self.current_time_s += self.outer_dt_s;

        Ok(())
    }

    /// Total number of outer steps to reach `total_duration_s`.
    #[must_use]
    pub fn total_outer_steps(&self) -> usize {
        (self.total_duration_s / self.outer_dt_s).ceil() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiscale_subcycling_coordination() {
        // Outer step: 1.0 s, Inner sub-cycles: 10 => Inner dt = 0.1 s
        let mut coord = MultiscaleCoordinator::new(5.0, 1.0, 10).unwrap();
        assert_eq!(coord.inner_dt_s(), 0.1);

        let mut inner_calls = 0;
        let mut outer_calls = 0;

        for _ in 0..5 {
            coord
                .step(
                    |_t, _dt| {
                        inner_calls += 1;
                        Ok(())
                    },
                    |_t, _dt| {
                        outer_calls += 1;
                        Ok(())
                    },
                )
                .unwrap();
        }

        assert_eq!(inner_calls, 50); // 5 * 10
        assert_eq!(outer_calls, 5);
        assert!((coord.current_time_s - 5.0).abs() < 1e-12);
    }
}
