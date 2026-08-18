//! # BioForge Types
//!
//! Physical dimensions, units, and quantities for the BioForge language.
//!
//! This crate implements the core type system for physical quantities,
//! ensuring that all numerical values in BioForge carry proper dimensional
//! information. Per the BioForge scientific principles:
//!
//! > "All physical quantities MUST carry units. Naked floating-point numbers
//! > for physical values are strictly forbidden in the runtime."
//!
//! ## Architecture
//!
//! - [`Dimension`] — Exponent-vector representation of physical dimensions
//!   (e.g., Length = \[1,0,0,0,0,0,0\], Energy = \[2,-2,1,0,0,0,0\])
//! - [`Unit`] — A named unit with its dimension and SI conversion factor
//! - [`UnitRegistry`] — Registry of all known units, used by the parser
//! - [`Quantity`] — A value paired with a unit, supporting dimensional arithmetic
//! - [`DimensionError`] — Errors for incompatible dimensional operations

#![deny(unsafe_code)]

mod dimension;
mod error;
mod quantity;
mod unit;

pub use dimension::Dimension;
pub use error::DimensionError;
pub use quantity::Quantity;
pub use unit::{Unit, UnitRegistry};
