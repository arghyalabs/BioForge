//! Physical units and the unit registry.
//!
//! Each [`Unit`] maps a human-readable name (e.g., `"nm"`, `"fs"`) to its
//! physical [`Dimension`] and a conversion factor to the SI base unit.
//!
//! The [`UnitRegistry`] provides a single source of truth for all known
//! units, replacing the hardcoded `KNOWN_UNITS` list in the parser.

use crate::Dimension;
use std::collections::HashMap;

/// A physical unit with its dimension and conversion factor.
///
/// The `to_si` factor converts a value in this unit to the corresponding
/// SI base unit by multiplication:
///
/// ```text
/// value_in_si = value * unit.to_si
/// ```
///
/// # Examples
///
/// - `nm`: dimension = Length, to_si = 1e-9 (1 nm = 1e-9 m)
/// - `fs`: dimension = Time, to_si = 1e-15 (1 fs = 1e-15 s)
/// - `K`:  dimension = Temperature, to_si = 1.0 (Kelvin is the SI base)
#[derive(Debug, Clone, PartialEq)]
pub struct Unit {
    /// The unit symbol as written in BioForge source code (e.g., `"nm"`, `"fs"`).
    pub name: String,
    /// The physical dimension of this unit.
    pub dimension: Dimension,
    /// Multiplication factor to convert to SI base unit.
    pub to_si: f64,
    /// The name of the SI base unit (e.g., `"m"`, `"s"`, `"K"`).
    pub si_base: String,
}

/// Registry of all known units in the BioForge language.
///
/// This is the **single source of truth** for unit recognition. The parser
/// consults this registry to determine whether an identifier following a
/// number is a unit (forming a `Quantity`) or a regular identifier.
///
/// # Supported Unit Categories
///
/// | Category       | Units                         |
/// |----------------|-------------------------------|
/// | Length         | nm, um, mm, m, Å              |
/// | Time           | fs, ps, ns, us, ms, s         |
/// | Mass           | Da, kDa                       |
/// | Temperature    | K                             |
/// | Amount         | mol                           |
/// | Concentration  | nM, uM, mM, M                |
/// | Voltage        | mV, V                         |
/// | Pressure       | Pa, kPa, MPa, atm             |
/// | Frequency      | Hz                            |
#[derive(Debug)]
pub struct UnitRegistry {
    units: HashMap<String, Unit>,
}

impl UnitRegistry {
    /// Create a new registry pre-populated with all built-in BioForge units.
    #[must_use]
    pub fn new() -> Self {
        let mut registry = Self {
            units: HashMap::new(),
        };
        registry.register_builtin_units();
        registry
    }

    /// Check whether a given name is a known unit.
    ///
    /// This is used by the parser to decide whether an identifier following
    /// a number should be consumed as a unit token.
    #[must_use]
    pub fn is_known(&self, name: &str) -> bool {
        self.units.contains_key(name)
    }

    /// Look up a unit by name.
    #[must_use]
    pub fn resolve(&self, name: &str) -> Option<&Unit> {
        self.units.get(name)
    }

    /// Get all registered unit names.
    #[must_use]
    pub fn all_unit_names(&self) -> Vec<&str> {
        self.units.keys().map(|s| s.as_str()).collect()
    }

    /// Convert a value from one unit to another.
    ///
    /// Returns `None` if the units have incompatible dimensions.
    #[must_use]
    pub fn convert(&self, value: f64, from: &Unit, to: &Unit) -> Option<f64> {
        if !from.dimension.is_compatible(&to.dimension) {
            return None;
        }
        // value_in_from * from.to_si = value_in_si
        // value_in_si / to.to_si = value_in_to
        Some(value * from.to_si / to.to_si)
    }

    /// Register a single unit.
    fn register(&mut self, name: &str, dimension: Dimension, to_si: f64, si_base: &str) {
        self.units.insert(
            name.to_string(),
            Unit {
                name: name.to_string(),
                dimension,
                to_si,
                si_base: si_base.to_string(),
            },
        );
    }

    /// Register all built-in units.
    fn register_builtin_units(&mut self) {
        // ── Length ───────────────────────────────────────────────────────
        let length = Dimension::length();
        self.register("m", length, 1.0, "m");
        self.register("mm", length, 1e-3, "m");
        self.register("um", length, 1e-6, "m");     // micrometer
        self.register("nm", length, 1e-9, "m");
        // Å (Ångström) — internal unit for molecular biology
        // Note: The lexer does not yet support Å as a token; this is
        // registered for future use and programmatic access.

        // ── Time ────────────────────────────────────────────────────────
        let time = Dimension::time();
        self.register("s", time, 1.0, "s");
        self.register("ms", time, 1e-3, "s");
        self.register("us", time, 1e-6, "s");        // microsecond
        self.register("ns", time, 1e-9, "s");
        self.register("ps", time, 1e-12, "s");
        self.register("fs", time, 1e-15, "s");

        // ── Mass ────────────────────────────────────────────────────────
        let mass = Dimension::mass();
        self.register("Da", mass, 1.660_539_066_6e-27, "kg"); // dalton
        self.register("kDa", mass, 1.660_539_066_6e-24, "kg"); // kilodalton

        // ── Temperature ─────────────────────────────────────────────────
        let temp = Dimension::temperature();
        self.register("K", temp, 1.0, "K");

        // ── Amount of Substance ─────────────────────────────────────────
        let amount = Dimension::amount();
        self.register("mol", amount, 1.0, "mol");

        // ── Concentration (Amount / Volume = [N L⁻³]) ──────────────────
        let conc = Dimension::concentration();
        self.register("M", conc, 1000.0, "mol/m^3");   // 1 M = 1000 mol/m³
        self.register("mM", conc, 1.0, "mol/m^3");     // 1 mM = 1 mol/m³
        self.register("uM", conc, 1e-3, "mol/m^3");    // 1 µM = 0.001 mol/m³
        self.register("nM", conc, 1e-6, "mol/m^3");    // 1 nM = 1e-6 mol/m³

        // ── Voltage (Electric Potential) ────────────────────────────────
        let voltage = Dimension::voltage();
        self.register("V", voltage, 1.0, "V");
        self.register("mV", voltage, 1e-3, "V");

        // ── Pressure ────────────────────────────────────────────────────
        let pressure = Dimension::pressure();
        self.register("Pa", pressure, 1.0, "Pa");
        self.register("kPa", pressure, 1e3, "Pa");
        self.register("MPa", pressure, 1e6, "Pa");
        self.register("atm", pressure, 101_325.0, "Pa");

        // ── Frequency ───────────────────────────────────────────────────
        let freq = Dimension::frequency();
        self.register("Hz", freq, 1.0, "Hz");
    }
}

impl Default for UnitRegistry {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn registry() -> UnitRegistry {
        UnitRegistry::new()
    }

    #[test]
    fn test_all_known_units_resolve() {
        let reg = registry();
        let expected = [
            "m", "mm", "um", "nm",
            "s", "ms", "us", "ns", "ps", "fs",
            "Da", "kDa",
            "K",
            "mol",
            "M", "mM", "uM", "nM",
            "V", "mV",
            "Pa", "kPa", "MPa", "atm",
            "Hz",
        ];
        for name in expected {
            assert!(reg.is_known(name), "expected '{}' to be a known unit", name);
            assert!(reg.resolve(name).is_some(), "expected '{}' to resolve", name);
        }
    }

    #[test]
    fn test_unknown_units() {
        let reg = registry();
        assert!(!reg.is_known("furlongs"));
        assert!(!reg.is_known("cubits"));
        assert!(!reg.is_known("receptor"));
        assert!(reg.resolve("furlongs").is_none());
    }

    #[test]
    fn test_unit_dimensions() {
        let reg = registry();
        assert_eq!(reg.resolve("nm").unwrap().dimension, Dimension::length());
        assert_eq!(reg.resolve("fs").unwrap().dimension, Dimension::time());
        assert_eq!(reg.resolve("K").unwrap().dimension, Dimension::temperature());
        assert_eq!(reg.resolve("Da").unwrap().dimension, Dimension::mass());
        assert_eq!(reg.resolve("mol").unwrap().dimension, Dimension::amount());
        assert_eq!(reg.resolve("mM").unwrap().dimension, Dimension::concentration());
        assert_eq!(reg.resolve("mV").unwrap().dimension, Dimension::voltage());
        assert_eq!(reg.resolve("atm").unwrap().dimension, Dimension::pressure());
        assert_eq!(reg.resolve("Hz").unwrap().dimension, Dimension::frequency());
    }

    #[test]
    fn test_length_conversions() {
        let reg = registry();
        let nm = reg.resolve("nm").unwrap();
        let um = reg.resolve("um").unwrap();
        let mm = reg.resolve("mm").unwrap();

        // 1000 nm = 1 um
        let result = reg.convert(1000.0, nm, um).unwrap();
        assert!((result - 1.0).abs() < 1e-10, "got {}", result);

        // 1000 um = 1 mm
        let result = reg.convert(1000.0, um, mm).unwrap();
        assert!((result - 1.0).abs() < 1e-10, "got {}", result);
    }

    #[test]
    fn test_time_conversions() {
        let reg = registry();
        let fs = reg.resolve("fs").unwrap();
        let ps = reg.resolve("ps").unwrap();
        let ns = reg.resolve("ns").unwrap();

        // 1000 fs = 1 ps
        let result = reg.convert(1000.0, fs, ps).unwrap();
        assert!((result - 1.0).abs() < 1e-10, "got {}", result);

        // 1000 ps = 1 ns
        let result = reg.convert(1000.0, ps, ns).unwrap();
        assert!((result - 1.0).abs() < 1e-10, "got {}", result);
    }

    #[test]
    fn test_incompatible_conversion_returns_none() {
        let reg = registry();
        let nm = reg.resolve("nm").unwrap();
        let k = reg.resolve("K").unwrap();

        // Cannot convert nm to K
        assert!(reg.convert(1.0, nm, k).is_none());
    }

    #[test]
    fn test_pressure_conversion() {
        let reg = registry();
        let atm = reg.resolve("atm").unwrap();
        let pa = reg.resolve("Pa").unwrap();

        // 1 atm = 101325 Pa
        let result = reg.convert(1.0, atm, pa).unwrap();
        assert!((result - 101_325.0).abs() < 0.1, "got {}", result);
    }

    #[test]
    fn test_concentration_conversion() {
        let reg = registry();
        let m = reg.resolve("M").unwrap();
        let mm = reg.resolve("mM").unwrap();

        // 1 M = 1000 mM
        let result = reg.convert(1.0, m, mm).unwrap();
        assert!((result - 1000.0).abs() < 1e-10, "got {}", result);
    }
}
