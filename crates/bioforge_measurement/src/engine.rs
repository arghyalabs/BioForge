//! Scientific measurement engine, statistical aggregators, and data exporters (CSV/JSON).

use bioforge_state::SimulationState;
use serde::{Deserialize, Serialize};

use crate::error::MeasurementError;
use crate::observable::Observable;

/// Statistical summary of a recorded scientific time series.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TimeSeriesStatistics {
    /// Arithmetic mean $\bar{x} = \frac{1}{M}\sum x_k$.
    pub mean: f64,
    /// Sample standard deviation $s = \sqrt{\frac{1}{M-1}\sum (x_k - \bar{x})^2}$.
    pub std_dev: f64,
    /// Minimum observed value.
    pub min: f64,
    /// Maximum observed value.
    pub max: f64,
    /// Number of recorded samples.
    pub count: usize,
}

/// Recorded historical time series data for a single observable.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TimeSeries {
    /// Observable name descriptor.
    pub name: String,
    /// Physical unit.
    pub unit: String,
    /// Vector of (time_ps, value) measurement pairs.
    pub data: Vec<(f64, f64)>,
}

impl TimeSeries {
    /// Create a new empty time series.
    #[must_use]
    pub fn new(name: impl Into<String>, unit: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            unit: unit.into(),
            data: Vec::new(),
        }
    }

    /// Record a single measurement point.
    pub fn record(&mut self, time_ps: f64, value: f64) {
        self.data.push((time_ps, value));
    }

    /// Compute summary statistics over the recorded time series.
    #[must_use]
    pub fn statistics(&self) -> Option<TimeSeriesStatistics> {
        let count = self.data.len();
        if count == 0 {
            return None;
        }

        let mut sum = 0.0;
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;

        for &(_t, v) in &self.data {
            sum += v;
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }

        let mean = sum / (count as f64);

        let std_dev = if count > 1 {
            let mut sum_sq_diff = 0.0;
            for &(_t, v) in &self.data {
                let diff = v - mean;
                sum_sq_diff += diff * diff;
            }
            (sum_sq_diff / ((count - 1) as f64)).sqrt()
        } else {
            0.0
        };

        Some(TimeSeriesStatistics {
            mean,
            std_dev,
            min,
            max,
            count,
        })
    }
}

/// Measurement engine that orchestrates multiple observers and records time series.
#[derive(Debug, Default)]
pub struct MeasurementEngine {
    /// Active observers.
    pub observables: Vec<Box<dyn Observable>>,
    /// Recorded time series per observer.
    pub time_series: Vec<TimeSeries>,
}

impl MeasurementEngine {
    /// Create a new, empty measurement engine.
    #[must_use]
    pub fn new() -> Self {
        Self {
            observables: Vec::new(),
            time_series: Vec::new(),
        }
    }

    /// Register a new observable.
    pub fn add_observable<O: Observable + 'static>(&mut self, observable: O) {
        let ts = TimeSeries::new(observable.name(), observable.unit());
        self.observables.push(Box::new(observable));
        self.time_series.push(ts);
    }

    /// Record a simulation step across all registered observables.
    pub fn record_step(
        &mut self,
        state: &SimulationState,
        potential_energy: f64,
    ) -> Result<(), MeasurementError> {
        let time = state.time;
        for (i, obs) in self.observables.iter().enumerate() {
            let val = obs.evaluate(state, potential_energy)?;
            self.time_series[i].record(time, val);
        }
        Ok(())
    }

    /// Get reference to a time series by name.
    #[must_use]
    pub fn get_time_series(&self, name: &str) -> Option<&TimeSeries> {
        self.time_series.iter().find(|ts| ts.name == name)
    }

    /// Get statistics for an observable by name.
    #[must_use]
    pub fn statistics(&self, name: &str) -> Option<TimeSeriesStatistics> {
        self.get_time_series(name).and_then(|ts| ts.statistics())
    }

    /// Export all recorded time series into a unified multi-column CSV table.
    #[must_use]
    pub fn export_csv(&self) -> String {
        if self.time_series.is_empty() {
            return String::from("time_ps\n");
        }

        let mut header = String::from("time_ps");
        for ts in &self.time_series {
            header.push_str(&format!(",\"{}[{}]\"", ts.name, ts.unit));
        }
        header.push('\n');

        let num_rows = self.time_series[0].data.len();
        let mut out = header;

        for r in 0..num_rows {
            let time = self.time_series[0].data[r].0;
            out.push_str(&format!("{:.4}", time));

            for ts in &self.time_series {
                let val = if r < ts.data.len() {
                    ts.data[r].1
                } else {
                    0.0
                };
                out.push_str(&format!(",{}", val));
            }
            out.push('\n');
        }

        out
    }

    /// Export all recorded time series and their statistics as structured JSON.
    pub fn export_json(&self) -> Result<String, MeasurementError> {
        #[derive(Serialize)]
        struct JsonExport<'a> {
            observables: Vec<JsonObservableEntry<'a>>,
        }

        #[derive(Serialize)]
        struct JsonObservableEntry<'a> {
            name: &'a str,
            unit: &'a str,
            statistics: Option<TimeSeriesStatistics>,
            data_points: usize,
            data: &'a Vec<(f64, f64)>,
        }

        let entries: Vec<JsonObservableEntry> = self
            .time_series
            .iter()
            .map(|ts| JsonObservableEntry {
                name: &ts.name,
                unit: &ts.unit,
                statistics: ts.statistics(),
                data_points: ts.data.len(),
                data: &ts.data,
            })
            .collect();

        let export = JsonExport {
            observables: entries,
        };

        serde_json::to_string_pretty(&export)
            .map_err(|e| MeasurementError::SerializationError(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observable::{DistanceObservable, KineticEnergyObservable, TemperatureObservable};

    #[test]
    fn test_time_series_statistics() {
        let mut ts = TimeSeries::new("test", "A");
        ts.record(0.0, 1.0);
        ts.record(0.1, 2.0);
        ts.record(0.2, 3.0);
        ts.record(0.3, 4.0);
        ts.record(0.4, 5.0);

        let stats = ts.statistics().unwrap();
        assert_eq!(stats.count, 5);
        assert!((stats.mean - 3.0).abs() < 1e-6);
        assert!((stats.min - 1.0).abs() < 1e-6);
        assert!((stats.max - 5.0).abs() < 1e-6);
        // Variance of [1,2,3,4,5] = (4+1+0+1+4)/4 = 2.5 => std_dev = sqrt(2.5) ~ 1.5811388
        assert!((stats.std_dev - 2.5_f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_measurement_engine_recording_and_export() {
        let mut state = SimulationState::empty();
        state.num_atoms = 2;
        state.masses = vec![12.011, 12.011];
        state.positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        state.velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];

        let mut engine = MeasurementEngine::new();
        engine.add_observable(DistanceObservable::pair("C1-C2", 0, 1));
        engine.add_observable(KineticEnergyObservable);
        engine.add_observable(TemperatureObservable);

        engine.record_step(&state, 0.0).unwrap();

        assert_eq!(engine.time_series.len(), 3);
        assert_eq!(engine.time_series[0].data.len(), 1);
        assert!((engine.time_series[0].data[0].1 - 3.0).abs() < 1e-6);

        let csv = engine.export_csv();
        assert!(csv.contains("time_ps,\"C1-C2[Å]\""));
        assert!(csv.contains("0.0000,3"));

        let json = engine.export_json().unwrap();
        assert!(json.contains("C1-C2"));
        assert!(json.contains("kinetic_energy"));
    }
}
