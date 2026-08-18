//! Multicellular population dynamics, Verhulst logistic growth, and contact inhibition.

use serde::{Deserialize, Serialize};

/// Multicellular tissue population governed by logistic proliferation and contact inhibition.
///
/// Models macroscopic tissue growth and density-dependent growth arrest (Verhulst 1838):
///
/// $$\frac{dN}{dt} = r N \left(1 - \frac{N}{K}\right) - \delta_{\text{apop}} N$$
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CellPopulation {
    /// Initial cell count $N_0$.
    pub initial_count: f64,
    /// Environmental / tissue carrying capacity $K$ (contact inhibition limit).
    pub carrying_capacity: f64,
    /// Intrinsic cellular division rate $r$ in $\text{hour}^{-1}$.
    pub division_rate_per_hour: f64,
    /// Baseline apoptosis / cell death rate $\delta$ in $\text{hour}^{-1}$.
    pub apoptosis_rate_per_hour: f64,
}

impl CellPopulation {
    /// Create a new cell population model.
    #[must_use]
    pub fn new(
        initial_count: f64,
        carrying_capacity: f64,
        division_rate_per_hour: f64,
        apoptosis_rate_per_hour: f64,
    ) -> Self {
        Self {
            initial_count: initial_count.max(1.0),
            carrying_capacity: carrying_capacity.max(1.0),
            division_rate_per_hour: division_rate_per_hour.max(0.0),
            apoptosis_rate_per_hour: apoptosis_rate_per_hour.max(0.0),
        }
    }

    /// Net proliferation rate $r_{\text{net}} = r - \delta$ in $\text{hour}^{-1}$.
    #[must_use]
    pub fn net_growth_rate(&self) -> f64 {
        self.division_rate_per_hour - self.apoptosis_rate_per_hour
    }

    /// Exact analytical population size $N(t)$ at elapsed time $t$ in hours:
    ///
    /// $$N(t) = \frac{K N_0 e^{r t}}{K + N_0 (e^{r t} - 1)}$$
    #[must_use]
    pub fn analytical_count_at_hours(&self, t_hours: f64) -> f64 {
        let r = self.net_growth_rate();
        let k = self.carrying_capacity;
        let n0 = self.initial_count;

        if r.abs() < 1e-12 {
            return n0;
        }

        let exp_rt = (r * t_hours).exp();
        (k * n0 * exp_rt) / (k + n0 * (exp_rt - 1.0))
    }

    /// Time derivative $dN/dt$ given current population count $N$.
    #[must_use]
    pub fn compute_derivative(&self, count: f64) -> f64 {
        let n = count.max(0.0);
        let k = self.carrying_capacity;
        self.division_rate_per_hour * n * (1.0 - n / k) - self.apoptosis_rate_per_hour * n
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_verhulst_logistic_growth_analytical() {
        let pop = CellPopulation::new(100.0, 10_000.0, 0.1, 0.0);

        // At t = 0 -> N = 100
        assert_eq!(pop.analytical_count_at_hours(0.0), 100.0);

        // At large t -> N approaches K (10,000)
        let n_large = pop.analytical_count_at_hours(200.0);
        assert!((n_large - 10_000.0).abs() < 1.0);

        // At carrying capacity, derivative must be zero
        let d_at_k = pop.compute_derivative(10_000.0);
        assert!(d_at_k.abs() < 1e-12);
    }
}
