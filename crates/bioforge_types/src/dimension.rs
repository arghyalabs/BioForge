//! Physical dimensions using an exponent-vector representation.
//!
//! Every physical quantity has a dimension that is a product of powers of
//! the seven SI base dimensions. We represent this as an array of 7 integer
//! exponents: `[Length, Time, Mass, Temperature, Amount, Current, LuminousIntensity]`.
//!
//! # Examples
//!
//! | Quantity     | Exponents                |
//! |-------------|--------------------------|
//! | Length       | `[1, 0, 0, 0, 0, 0, 0]` |
//! | Velocity     | `[1,-1, 0, 0, 0, 0, 0]` |
//! | Force        | `[1,-2, 1, 0, 0, 0, 0]` |
//! | Energy       | `[2,-2, 1, 0, 0, 0, 0]` |
//! | Pressure     | `[-1,-2, 1, 0, 0, 0, 0]`|
//! | Dimensionless| `[0, 0, 0, 0, 0, 0, 0]` |

use std::fmt;

/// Index constants for the 7 SI base dimension exponents.
/// These document the array layout and are used in future phases
/// for direct index access.
#[allow(dead_code)]
const LENGTH: usize = 0;
#[allow(dead_code)]
const TIME: usize = 1;
#[allow(dead_code)]
const MASS: usize = 2;
#[allow(dead_code)]
const TEMPERATURE: usize = 3;
#[allow(dead_code)]
const AMOUNT: usize = 4;
#[allow(dead_code)]
const CURRENT: usize = 5;
#[allow(dead_code)]
const LUMINOUS_INTENSITY: usize = 6;

/// The names of the base dimensions, used for display.
const BASE_NAMES: [&str; 7] = [
    "Length",
    "Time",
    "Mass",
    "Temperature",
    "Amount",
    "Current",
    "LuminousIntensity",
];

/// Physical dimension represented as exponents of the 7 SI base dimensions.
///
/// Two dimensions are compatible (for addition/subtraction) if and only if
/// all 7 exponents are equal. Multiplication adds exponents; division
/// subtracts them.
///
/// # Scientific Principle
///
/// > "Units are first-class citizens. A number without a unit is
/// > scientifically meaningless."
/// >
/// > — BioForge Scientific Principles, Principle 3
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dimension {
    exponents: [i8; 7],
}

impl Dimension {
    /// Create a dimension from raw exponents.
    ///
    /// Order: `[Length, Time, Mass, Temperature, Amount, Current, LuminousIntensity]`
    #[must_use]
    pub const fn new(exponents: [i8; 7]) -> Self {
        Self { exponents }
    }

    /// The dimensionless quantity (all exponents zero).
    #[must_use]
    pub const fn dimensionless() -> Self {
        Self::new([0, 0, 0, 0, 0, 0, 0])
    }

    /// Length dimension: `[L]`
    #[must_use]
    pub const fn length() -> Self {
        Self::new([1, 0, 0, 0, 0, 0, 0])
    }

    /// Time dimension: `[T]`
    #[must_use]
    pub const fn time() -> Self {
        Self::new([0, 1, 0, 0, 0, 0, 0])
    }

    /// Mass dimension: `[M]`
    #[must_use]
    pub const fn mass() -> Self {
        Self::new([0, 0, 1, 0, 0, 0, 0])
    }

    /// Temperature dimension: `[Θ]`
    #[must_use]
    pub const fn temperature() -> Self {
        Self::new([0, 0, 0, 1, 0, 0, 0])
    }

    /// Amount of substance dimension: `[N]`
    #[must_use]
    pub const fn amount() -> Self {
        Self::new([0, 0, 0, 0, 1, 0, 0])
    }

    /// Electric current dimension: `[I]`
    #[must_use]
    pub const fn current() -> Self {
        Self::new([0, 0, 0, 0, 0, 1, 0])
    }

    /// Velocity: `[L T⁻¹]`
    #[must_use]
    pub const fn velocity() -> Self {
        Self::new([1, -1, 0, 0, 0, 0, 0])
    }

    /// Acceleration: `[L T⁻²]`
    #[must_use]
    pub const fn acceleration() -> Self {
        Self::new([1, -2, 0, 0, 0, 0, 0])
    }

    /// Force: `[M L T⁻²]`
    #[must_use]
    pub const fn force() -> Self {
        Self::new([1, -2, 1, 0, 0, 0, 0])
    }

    /// Energy: `[M L² T⁻²]`
    #[must_use]
    pub const fn energy() -> Self {
        Self::new([2, -2, 1, 0, 0, 0, 0])
    }

    /// Pressure: `[M L⁻¹ T⁻²]`
    #[must_use]
    pub const fn pressure() -> Self {
        Self::new([-1, -2, 1, 0, 0, 0, 0])
    }

    /// Voltage (electric potential): `[M L² T⁻³ I⁻¹]`
    #[must_use]
    pub const fn voltage() -> Self {
        Self::new([2, -3, 1, 0, 0, -1, 0])
    }

    /// Concentration: `[N L⁻³]` (amount per volume)
    #[must_use]
    pub const fn concentration() -> Self {
        Self::new([-3, 0, 0, 0, 1, 0, 0])
    }

    /// Frequency: `[T⁻¹]`
    #[must_use]
    pub const fn frequency() -> Self {
        Self::new([0, -1, 0, 0, 0, 0, 0])
    }

    /// Check if two dimensions are compatible (identical exponents).
    #[must_use]
    pub fn is_compatible(&self, other: &Self) -> bool {
        self.exponents == other.exponents
    }

    /// Check if this is the dimensionless dimension.
    #[must_use]
    pub fn is_dimensionless(&self) -> bool {
        self.exponents == [0, 0, 0, 0, 0, 0, 0]
    }

    /// Multiply two dimensions (add exponents).
    ///
    /// `[L] * [T⁻¹] = [L T⁻¹]` (velocity)
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = [0i8; 7];
        for i in 0..7 {
            result[i] = self.exponents[i] + other.exponents[i];
        }
        Self::new(result)
    }

    /// Divide two dimensions (subtract exponents).
    ///
    /// `[L] / [T] = [L T⁻¹]` (velocity)
    #[must_use]
    pub fn div(&self, other: &Self) -> Self {
        let mut result = [0i8; 7];
        for i in 0..7 {
            result[i] = self.exponents[i] - other.exponents[i];
        }
        Self::new(result)
    }

    /// Invert the dimension (negate all exponents).
    ///
    /// `[T]⁻¹ = [T⁻¹]`
    #[must_use]
    pub fn inverse(&self) -> Self {
        let mut result = [0i8; 7];
        for i in 0..7 {
            result[i] = -self.exponents[i];
        }
        Self::new(result)
    }

    /// Get the raw exponents array.
    #[must_use]
    pub const fn exponents(&self) -> &[i8; 7] {
        &self.exponents
    }

    /// Return a human-readable name for well-known dimensions,
    /// or a formula representation for derived dimensions.
    #[must_use]
    pub fn name(&self) -> String {
        // Check well-known dimensions first
        let known: &[(Self, &str)] = &[
            (Self::dimensionless(), "Dimensionless"),
            (Self::length(), "Length"),
            (Self::time(), "Time"),
            (Self::mass(), "Mass"),
            (Self::temperature(), "Temperature"),
            (Self::amount(), "Amount"),
            (Self::current(), "Current"),
            (Self::velocity(), "Velocity"),
            (Self::acceleration(), "Acceleration"),
            (Self::force(), "Force"),
            (Self::energy(), "Energy"),
            (Self::pressure(), "Pressure"),
            (Self::voltage(), "Voltage"),
            (Self::concentration(), "Concentration"),
            (Self::frequency(), "Frequency"),
        ];

        for (dim, name) in known {
            if self == dim {
                return name.to_string();
            }
        }

        // Build a formula representation for unknown derived dimensions
        self.formula()
    }

    /// Build a formula string like `"Length² × Mass / Time²"`.
    fn formula(&self) -> String {
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();

        for (i, &exp) in self.exponents.iter().enumerate() {
            if exp > 0 {
                if exp == 1 {
                    numerator.push(BASE_NAMES[i].to_string());
                } else {
                    numerator.push(format!("{}^{}", BASE_NAMES[i], exp));
                }
            } else if exp < 0 {
                let abs_exp = -exp;
                if abs_exp == 1 {
                    denominator.push(BASE_NAMES[i].to_string());
                } else {
                    denominator.push(format!("{}^{}", BASE_NAMES[i], abs_exp));
                }
            }
        }

        if numerator.is_empty() && denominator.is_empty() {
            return "Dimensionless".to_string();
        }

        let num_str = if numerator.is_empty() {
            "1".to_string()
        } else {
            numerator.join(" * ")
        };

        if denominator.is_empty() {
            num_str
        } else {
            format!("{} / {}", num_str, denominator.join(" * "))
        }
    }
}

impl fmt::Display for Dimension {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_dimensions() {
        assert_eq!(Dimension::length().exponents(), &[1, 0, 0, 0, 0, 0, 0]);
        assert_eq!(Dimension::time().exponents(), &[0, 1, 0, 0, 0, 0, 0]);
        assert_eq!(Dimension::mass().exponents(), &[0, 0, 1, 0, 0, 0, 0]);
        assert_eq!(Dimension::temperature().exponents(), &[0, 0, 0, 1, 0, 0, 0]);
        assert_eq!(Dimension::amount().exponents(), &[0, 0, 0, 0, 1, 0, 0]);
        assert_eq!(Dimension::current().exponents(), &[0, 0, 0, 0, 0, 1, 0]);
    }

    #[test]
    fn test_derived_dimensions() {
        assert_eq!(Dimension::velocity().exponents(), &[1, -1, 0, 0, 0, 0, 0]);
        assert_eq!(Dimension::force().exponents(), &[1, -2, 1, 0, 0, 0, 0]);
        assert_eq!(Dimension::energy().exponents(), &[2, -2, 1, 0, 0, 0, 0]);
        assert_eq!(Dimension::pressure().exponents(), &[-1, -2, 1, 0, 0, 0, 0]);
    }

    #[test]
    fn test_dimensionless() {
        assert!(Dimension::dimensionless().is_dimensionless());
        assert!(!Dimension::length().is_dimensionless());
    }

    #[test]
    fn test_compatibility() {
        assert!(Dimension::length().is_compatible(&Dimension::length()));
        assert!(!Dimension::length().is_compatible(&Dimension::time()));
        assert!(!Dimension::length().is_compatible(&Dimension::temperature()));
    }

    #[test]
    fn test_mul_gives_velocity() {
        // Length * Time⁻¹ = Velocity
        let velocity = Dimension::length().mul(&Dimension::time().inverse());
        assert_eq!(velocity, Dimension::velocity());
    }

    #[test]
    fn test_mul_gives_energy() {
        // Force * Length = Energy: [1,-2,1,0,0,0,0] * [1,0,0,0,0,0,0] = [2,-2,1,0,0,0,0]
        let energy = Dimension::force().mul(&Dimension::length());
        assert_eq!(energy, Dimension::energy());
    }

    #[test]
    fn test_div_gives_velocity() {
        // Length / Time = Velocity
        let velocity = Dimension::length().div(&Dimension::time());
        assert_eq!(velocity, Dimension::velocity());
    }

    #[test]
    fn test_div_gives_acceleration() {
        // Velocity / Time = Acceleration
        let accel = Dimension::velocity().div(&Dimension::time());
        assert_eq!(accel, Dimension::acceleration());
    }

    #[test]
    fn test_inverse() {
        let inv_time = Dimension::time().inverse();
        assert_eq!(inv_time, Dimension::frequency());
    }

    #[test]
    fn test_self_div_is_dimensionless() {
        let result = Dimension::length().div(&Dimension::length());
        assert!(result.is_dimensionless());
    }

    #[test]
    fn test_name_known_dimensions() {
        assert_eq!(Dimension::length().name(), "Length");
        assert_eq!(Dimension::time().name(), "Time");
        assert_eq!(Dimension::mass().name(), "Mass");
        assert_eq!(Dimension::temperature().name(), "Temperature");
        assert_eq!(Dimension::energy().name(), "Energy");
        assert_eq!(Dimension::force().name(), "Force");
        assert_eq!(Dimension::velocity().name(), "Velocity");
        assert_eq!(Dimension::pressure().name(), "Pressure");
        assert_eq!(Dimension::voltage().name(), "Voltage");
        assert_eq!(Dimension::concentration().name(), "Concentration");
        assert_eq!(Dimension::frequency().name(), "Frequency");
        assert_eq!(Dimension::dimensionless().name(), "Dimensionless");
    }

    #[test]
    fn test_name_unknown_derived() {
        // Length² (not a named dimension)
        let dim = Dimension::new([2, 0, 0, 0, 0, 0, 0]);
        assert_eq!(dim.name(), "Length^2");
    }

    #[test]
    fn test_display() {
        assert_eq!(format!("{}", Dimension::length()), "Length");
        assert_eq!(format!("{}", Dimension::energy()), "Energy");
    }

    #[test]
    fn test_force_equals_mass_times_acceleration() {
        let computed = Dimension::mass().mul(&Dimension::acceleration());
        assert_eq!(computed, Dimension::force());
    }
}
