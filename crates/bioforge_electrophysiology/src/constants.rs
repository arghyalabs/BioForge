//! Fundamental physical and biophysical constants for electrophysiology.

/// Universal Faraday Constant $F = 96485.33212\text{ C/mol}$ (CODATA 2018).
pub const FARADAY_CONSTANT_F: f64 = 96485.332_12;

/// Molar Gas Constant $R = 8.314462618\text{ J}/(\text{mol}\cdot\text{K})$ (CODATA 2018).
pub const MOLAR_GAS_CONSTANT_R: f64 = 8.314_462_618;

/// Standard human physiological body temperature: $37.0^\circ\text{C} = 310.15\text{ K}$.
pub const BODY_TEMPERATURE_KELVIN: f64 = 310.15;

/// Standard room temperature: $20.0^\circ\text{C} = 293.15\text{ K}$.
pub const ROOM_TEMPERATURE_KELVIN: f64 = 293.15;

/// Classic Hodgkin-Huxley squid giant axon experimental temperature: $6.3^\circ\text{C} = 279.45\text{ K}$.
pub const SQUID_AXON_TEMPERATURE_KELVIN: f64 = 279.45;

/// Standard biological specific membrane capacitance: $1.0\,\mu\text{F/cm}^2$.
pub const DEFAULT_MEMBRANE_CAPACITANCE_UF_PER_CM2: f64 = 1.0;
