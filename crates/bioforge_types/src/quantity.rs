//! Physical quantities: a value paired with a unit.
//!
//! [`Quantity`] is the core runtime type for all physical values in BioForge.
//! It enforces dimensional safety: you can add two lengths, but attempting
//! to add a length to a temperature produces a [`DimensionError`].
//!
//! # Examples
//!
//! ```text
//! 5 nm + 3 nm  →  8 nm       ✓ (same dimension)
//! 10 nm * 2.0  →  20 nm      ✓ (scalar multiplication)
//! 5 nm + 310 K →  ERROR      ✗ (incompatible dimensions)
//! ```

use crate::error::DimensionError;
use crate::unit::Unit;
use std::fmt;

/// A physical quantity: a numerical value with an associated unit.
///
/// Per the BioForge scientific principles:
/// > "All physical quantities MUST carry units. Naked floating-point numbers
/// > for physical values are strictly forbidden in the runtime."
///
/// # Dimensional Safety
///
/// - **Addition/Subtraction**: Both operands must have the same dimension.
///   The result is expressed in the left operand's unit.
/// - **Multiplication/Division**: Dimensions combine (exponents add/subtract).
/// - **Scalar operations**: Multiplying or dividing by a dimensionless number
///   preserves the original dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    /// The numerical value.
    pub value: f64,
    /// The unit (which carries its dimension).
    pub unit: Unit,
}

impl Quantity {
    /// Create a new quantity.
    #[must_use]
    pub fn new(value: f64, unit: Unit) -> Self {
        Self { value, unit }
    }

    /// Add two quantities. Both must have the same dimension.
    ///
    /// The result is expressed in the left operand's unit. If the right
    /// operand has a different unit of the same dimension, its value is
    /// converted first.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionError::Incompatible`] if the dimensions differ.
    pub fn add(&self, other: &Self) -> Result<Self, DimensionError> {
        if !self.unit.dimension.is_compatible(&other.unit.dimension) {
            return Err(DimensionError::Incompatible {
                op: "add".to_string(),
                left: self.unit.dimension.name(),
                right: other.unit.dimension.name(),
            });
        }

        // Convert other's value to self's unit
        let other_converted = other.value * other.unit.to_si / self.unit.to_si;
        Ok(Self::new(self.value + other_converted, self.unit.clone()))
    }

    /// Subtract two quantities. Both must have the same dimension.
    ///
    /// The result is expressed in the left operand's unit.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionError::Incompatible`] if the dimensions differ.
    pub fn sub(&self, other: &Self) -> Result<Self, DimensionError> {
        if !self.unit.dimension.is_compatible(&other.unit.dimension) {
            return Err(DimensionError::Incompatible {
                op: "subtract".to_string(),
                left: self.unit.dimension.name(),
                right: other.unit.dimension.name(),
            });
        }

        let other_converted = other.value * other.unit.to_si / self.unit.to_si;
        Ok(Self::new(self.value - other_converted, self.unit.clone()))
    }

    /// Multiply this quantity by a dimensionless scalar.
    #[must_use]
    pub fn scale(&self, factor: f64) -> Self {
        Self::new(self.value * factor, self.unit.clone())
    }

    /// Convert this quantity to a different unit of the same dimension.
    ///
    /// # Errors
    ///
    /// Returns [`DimensionError::ConversionError`] if the dimensions differ.
    pub fn convert_to(&self, target: &Unit) -> Result<Self, DimensionError> {
        if !self.unit.dimension.is_compatible(&target.dimension) {
            return Err(DimensionError::ConversionError {
                from: self.unit.name.clone(),
                to: target.name.clone(),
            });
        }

        let value_in_si = self.value * self.unit.to_si;
        let value_in_target = value_in_si / target.to_si;
        Ok(Self::new(value_in_target, target.clone()))
    }

    /// Get the value in SI base units.
    #[must_use]
    pub fn to_si_value(&self) -> f64 {
        self.value * self.unit.to_si
    }
}

impl fmt::Display for Quantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Display integers without decimal point
        if self.value == (self.value as i64) as f64 {
            write!(f, "{} {}", self.value as i64, self.unit.name)
        } else {
            write!(f, "{} {}", self.value, self.unit.name)
        }
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UnitRegistry;

    fn reg() -> UnitRegistry {
        UnitRegistry::new()
    }

    fn make_qty(value: f64, unit_name: &str) -> Quantity {
        let registry = reg();
        let unit = registry.resolve(unit_name).unwrap().clone();
        Quantity::new(value, unit)
    }

    #[test]
    fn test_add_same_unit() {
        let a = make_qty(5.0, "nm");
        let b = make_qty(3.0, "nm");
        let result = a.add(&b).unwrap();
        assert!((result.value - 8.0).abs() < 1e-10);
        assert_eq!(result.unit.name, "nm");
    }

    #[test]
    fn test_add_compatible_different_units() {
        // 1 um + 500 nm = 1.5 um
        let a = make_qty(1.0, "um");
        let b = make_qty(500.0, "nm");
        let result = a.add(&b).unwrap();
        assert!((result.value - 1.5).abs() < 1e-10, "got {}", result.value);
        assert_eq!(result.unit.name, "um");
    }

    #[test]
    fn test_add_incompatible_dimensions() {
        let length = make_qty(5.0, "nm");
        let temp = make_qty(310.0, "K");
        let err = length.add(&temp).unwrap_err();
        assert!(matches!(err, DimensionError::Incompatible { .. }));
        assert_eq!(
            err.to_string(),
            "cannot add Length and Temperature: incompatible dimensions"
        );
    }

    #[test]
    fn test_sub_same_unit() {
        let a = make_qty(10.0, "ps");
        let b = make_qty(3.0, "ps");
        let result = a.sub(&b).unwrap();
        assert!((result.value - 7.0).abs() < 1e-10);
        assert_eq!(result.unit.name, "ps");
    }

    #[test]
    fn test_sub_incompatible() {
        let a = make_qty(10.0, "ps");
        let b = make_qty(5.0, "nm");
        assert!(a.sub(&b).is_err());
    }

    #[test]
    fn test_scale() {
        let q = make_qty(5.0, "nm");
        let result = q.scale(3.0);
        assert!((result.value - 15.0).abs() < 1e-10);
        assert_eq!(result.unit.name, "nm");
    }

    #[test]
    fn test_convert_to_same_dimension() {
        let registry = reg();
        let q = make_qty(1000.0, "nm");
        let target = registry.resolve("um").unwrap();
        let result = q.convert_to(target).unwrap();
        assert!((result.value - 1.0).abs() < 1e-10, "got {}", result.value);
        assert_eq!(result.unit.name, "um");
    }

    #[test]
    fn test_convert_to_incompatible() {
        let registry = reg();
        let q = make_qty(100.0, "nm");
        let target = registry.resolve("K").unwrap();
        let err = q.convert_to(target).unwrap_err();
        assert!(matches!(err, DimensionError::ConversionError { .. }));
    }

    #[test]
    fn test_to_si_value() {
        let q = make_qty(5.0, "nm");
        let si = q.to_si_value();
        assert!((si - 5e-9).abs() < 1e-20, "got {}", si);
    }

    #[test]
    fn test_display_integer() {
        let q = make_qty(310.0, "K");
        assert_eq!(format!("{}", q), "310 K");
    }

    #[test]
    fn test_display_float() {
        let q = make_qty(7.4, "mM");
        assert_eq!(format!("{}", q), "7.4 mM");
    }

    #[test]
    fn test_time_conversion_chain() {
        // 1000 fs → ps → ns
        let registry = reg();
        let q = make_qty(1000.0, "fs");
        let ps_unit = registry.resolve("ps").unwrap();
        let ps = q.convert_to(ps_unit).unwrap();
        assert!((ps.value - 1.0).abs() < 1e-10);

        let ns_unit = registry.resolve("ns").unwrap();
        let ns = make_qty(1000.0, "ps").convert_to(ns_unit).unwrap();
        assert!((ns.value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_then_convert() {
        // 500 nm + 500 nm = 1000 nm, then convert to um
        let registry = reg();
        let a = make_qty(500.0, "nm");
        let b = make_qty(500.0, "nm");
        let sum = a.add(&b).unwrap();
        assert!((sum.value - 1000.0).abs() < 1e-10);

        let um_unit = registry.resolve("um").unwrap();
        let converted = sum.convert_to(um_unit).unwrap();
        assert!((converted.value - 1.0).abs() < 1e-10);
    }
}
