//! Hodgkin-Huxley numerical action potential integrator and electrophysiology trajectory recorder.

use serde::{Deserialize, Serialize};

use crate::channel::{
    alpha_h, alpha_m, alpha_n, beta_h, beta_m, beta_n, steady_state_gate, LeakChannel,
    PotassiumChannel, SodiumChannel,
};
use crate::constants::DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2;
use crate::error::ElectrophysiologyError;
use crate::stimulus::StimulusProtocol;

/// Instantaneous physical state of a cell membrane.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElectrophysiologyState {
    /// Transmembrane voltage $V_m$ in millivolts ($\text{mV}$).
    pub v_membrane: f64,
    /// Fast Sodium activation gate ($m \in [0.0, 1.0]$).
    pub m: f64,
    /// Fast Sodium inactivation gate ($h \in [0.0, 1.0]$).
    pub h: f64,
    /// Delayed rectifier Potassium activation gate ($n \in [0.0, 1.0]$).
    pub n: f64,
}

/// Recorded electrophysiology trajectory across time.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct ElectrophysiologyTrajectory {
    pub times_ms: Vec<f64>,
    pub v_membrane_mv: Vec<f64>,
    pub m_gates: Vec<f64>,
    pub h_gates: Vec<f64>,
    pub n_gates: Vec<f64>,
    pub i_na_uA_cm2: Vec<f64>,
    pub i_k_uA_cm2: Vec<f64>,
    pub i_stim_uA_cm2: Vec<f64>,
}

impl ElectrophysiologyTrajectory {
    /// Create empty trajectory buffer.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a time point.
    pub fn record(
        &mut self,
        time_ms: f64,
        state: &ElectrophysiologyState,
        i_na: f64,
        i_k: f64,
        i_stim: f64,
    ) {
        self.times_ms.push(time_ms);
        self.v_membrane_mv.push(state.v_membrane);
        self.m_gates.push(state.m);
        self.h_gates.push(state.h);
        self.n_gates.push(state.n);
        self.i_na_uA_cm2.push(i_na);
        self.i_k_uA_cm2.push(i_k);
        self.i_stim_uA_cm2.push(i_stim);
    }

    /// Export recording into multi-column CSV table.
    #[must_use]
    pub fn export_csv(&self) -> String {
        let mut out = String::from("time_ms,v_membrane_mv,m_gate,h_gate,n_gate,i_na_uA_cm2,i_k_uA_cm2,i_stim_uA_cm2\n");
        for i in 0..self.times_ms.len() {
            out.push_str(&format!(
                "{:.4},{:.3},{:.4},{:.4},{:.4},{:.3},{:.3},{:.3}\n",
                self.times_ms[i],
                self.v_membrane_mv[i],
                self.m_gates[i],
                self.h_gates[i],
                self.n_gates[i],
                self.i_na_uA_cm2[i],
                self.i_k_uA_cm2[i],
                self.i_stim_uA_cm2[i]
            ));
        }
        out
    }
}

/// Hodgkin-Huxley nonlinear 4-variable electrophysiological model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HodgkinHuxleyModel {
    /// Specific membrane capacitance in $\mu\text{F/cm}^2$.
    pub capacitance_uF_cm2: f64,
    /// Absolute resting potential in $\text{mV}$ (default $-65.0\text{ mV}$).
    pub v_rest_mv: f64,
    /// Fast voltage-gated Sodium channel.
    pub sodium_channel: SodiumChannel,
    /// Delayed rectifier Potassium channel.
    pub potassium_channel: PotassiumChannel,
    /// Passive leak channel.
    pub leak_channel: LeakChannel,
}

impl Default for HodgkinHuxleyModel {
    fn default() -> Self {
        Self {
            capacitance_uF_cm2: DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2,
            v_rest_mv: -65.0,
            sodium_channel: SodiumChannel::default(),
            potassium_channel: PotassiumChannel::default(),
            leak_channel: LeakChannel::default(),
        }
    }
}

impl HodgkinHuxleyModel {
    /// Create a standard Hodgkin-Huxley model.
    #[must_use]
    pub fn standard() -> Self {
        Self::default()
    }

    /// Initial steady-state resting state at $V_m = V_{\text{rest}}$.
    #[must_use]
    pub fn initial_resting_state(&self) -> ElectrophysiologyState {
        ElectrophysiologyState {
            v_membrane: self.v_rest_mv,
            m: steady_state_gate(alpha_m(0.0), beta_m(0.0)),
            h: steady_state_gate(alpha_h(0.0), beta_h(0.0)),
            n: steady_state_gate(alpha_n(0.0), beta_n(0.0)),
        }
    }

    /// Compute state derivatives $[dV/dt, dm/dt, dh/dt, dn/dt]$.
    #[must_use]
    pub fn compute_derivatives(
        &self,
        state: &ElectrophysiologyState,
        i_stim: f64,
    ) -> [f64; 4] {
        let v = state.v_membrane;
        // Relative depolarization from resting potential
        let v_rel = v - self.v_rest_mv;

        // Channel currents
        let i_na = self.sodium_channel.g_bar * state.m.powi(3) * state.h * (v - self.sodium_channel.e_rev);
        let i_k = self.potassium_channel.g_bar * state.n.powi(4) * (v - self.potassium_channel.e_rev);
        let i_l = self.leak_channel.g_leak * (v - self.leak_channel.e_leak);

        // dV/dt = (I_stim - I_Na - I_K - I_L) / C_m
        let dv_dt = (i_stim - i_na - i_k - i_l) / self.capacitance_uF_cm2;

        // dm/dt = alpha_m * (1 - m) - beta_m * m
        let dm_dt = alpha_m(v_rel) * (1.0 - state.m) - beta_m(v_rel) * state.m;

        // dh/dt = alpha_h * (1 - h) - beta_h * h
        let dh_dt = alpha_h(v_rel) * (1.0 - state.h) - beta_h(v_rel) * state.h;

        // dn/dt = alpha_n * (1 - n) - beta_n * n
        let dn_dt = alpha_n(v_rel) * (1.0 - state.n) - beta_n(v_rel) * state.n;

        [dv_dt, dm_dt, dh_dt, dn_dt]
    }

    /// Run numerical simulation of the Hodgkin-Huxley model over `total_time_ms` using 4th-order Runge-Kutta.
    pub fn simulate(
        &self,
        total_time_ms: f64,
        dt_ms: f64,
        stimulus: &StimulusProtocol,
    ) -> Result<ElectrophysiologyTrajectory, ElectrophysiologyError> {
        let mut traj = ElectrophysiologyTrajectory::new();
        let mut state = self.initial_resting_state();
        let dt = dt_ms.clamp(0.0001, 0.05);

        let mut t = 0.0;
        let mut next_record_t = 0.0;
        let record_interval = dt * 10.0;

        while t <= total_time_ms + 1e-9 {
            let i_stim = stimulus.current_at(t);
            let i_na = self.sodium_channel.g_bar * state.m.powi(3) * state.h * (state.v_membrane - self.sodium_channel.e_rev);
            let i_k = self.potassium_channel.g_bar * state.n.powi(4) * (state.v_membrane - self.potassium_channel.e_rev);

            if t >= next_record_t - 1e-9 {
                traj.record(t, &state, i_na, i_k, i_stim);
                next_record_t += record_interval;
            }

            // RK4 Step
            // k1 = f(t, y)
            let k1 = self.compute_derivatives(&state, i_stim);

            // y2 = y + 0.5 * dt * k1
            let state2 = ElectrophysiologyState {
                v_membrane: state.v_membrane + 0.5 * dt * k1[0],
                m: (state.m + 0.5 * dt * k1[1]).clamp(0.0, 1.0),
                h: (state.h + 0.5 * dt * k1[2]).clamp(0.0, 1.0),
                n: (state.n + 0.5 * dt * k1[3]).clamp(0.0, 1.0),
            };
            let k2 = self.compute_derivatives(&state2, stimulus.current_at(t + 0.5 * dt));

            // y3 = y + 0.5 * dt * k2
            let state3 = ElectrophysiologyState {
                v_membrane: state.v_membrane + 0.5 * dt * k2[0],
                m: (state.m + 0.5 * dt * k2[1]).clamp(0.0, 1.0),
                h: (state.h + 0.5 * dt * k2[2]).clamp(0.0, 1.0),
                n: (state.n + 0.5 * dt * k2[3]).clamp(0.0, 1.0),
            };
            let k3 = self.compute_derivatives(&state3, stimulus.current_at(t + 0.5 * dt));

            // y4 = y + dt * k3
            let state4 = ElectrophysiologyState {
                v_membrane: state.v_membrane + dt * k3[0],
                m: (state.m + dt * k3[1]).clamp(0.0, 1.0),
                h: (state.h + dt * k3[2]).clamp(0.0, 1.0),
                n: (state.n + dt * k3[3]).clamp(0.0, 1.0),
            };
            let k4 = self.compute_derivatives(&state4, stimulus.current_at(t + dt));

            // Update state
            state.v_membrane += (dt / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
            state.m = (state.m + (dt / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])).clamp(0.0, 1.0);
            state.h = (state.h + (dt / 6.0) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2])).clamp(0.0, 1.0);
            state.n = (state.n + (dt / 6.0) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3])).clamp(0.0, 1.0);

            if !state.v_membrane.is_finite() {
                return Err(ElectrophysiologyError::SolverInstability {
                    message: format!("voltage diverged to NaN/Inf at t = {:.2} ms", t),
                });
            }

            t += dt;
        }

        Ok(traj)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Action Potential Generation Benchmark: Injected 10 uA/cm^2 current pulse fires a full action potential.
    #[test]
    fn test_action_potential_firing_and_overshoot() {
        let model = HodgkinHuxleyModel::standard();
        let stimulus = StimulusProtocol::CurrentPulse {
            amplitude_uA_cm2: 10.0,
            start_ms: 5.0,
            duration_ms: 1.0,
        };

        let traj = model.simulate(25.0, 0.01, &stimulus).unwrap();

        // Check resting potential before stimulus
        assert!((traj.v_membrane_mv[0] - (-65.0)).abs() < 1.0);

        // Find peak voltage during action potential
        let mut max_v = f64::NEG_INFINITY;
        for &v in &traj.v_membrane_mv {
            if v > max_v {
                max_v = v;
            }
        }

        // Action potential must overshoot 0 mV into positive territory (> +20 mV)
        assert!(
            max_v > 20.0,
            "expected action potential peak > +20 mV, got peak = {:.2} mV",
            max_v
        );

        // Check hyperpolarizing afterpotential (undershoots resting -65 mV)
        let mut min_v_after = f64::INFINITY;
        for (i, &t) in traj.times_ms.iter().enumerate() {
            if t > 8.0 {
                let v = traj.v_membrane_mv[i];
                if v < min_v_after {
                    min_v_after = v;
                }
            }
        }
        assert!(
            min_v_after < -68.0,
            "expected hyperpolarizing undershoot < -68 mV, got {:.2} mV",
            min_v_after
        );
    }
}
