//! # BioForge Biology
//!
//! Biological structure model for the BioForge language.
//!
//! This crate provides the foundational data types for representing
//! biological matter: atoms, bonds, molecules, and chemical elements.
//! It also includes a basic PDB file parser.
//!
//! ## Architecture
//!
//! This crate represents **structural data** — the topology and geometry
//! of biological systems. It is distinct from [`SimulationState`] (Phase 5),
//! which represents the dynamic state during a simulation.
//!
//! ```text
//! BioForge Source → Parser → AST → HIR → [Biology: Structure] → SimulationState
//! ```
//!
//! ## Internal Units
//!
//! Per SCIENTIFIC_PRINCIPLES.md, all quantities carry units:
//! - Positions: Ångströms (Å)
//! - Masses: Daltons (Da)
//! - Charges: Elementary charges (e)

#![deny(unsafe_code)]

mod atom;
mod bond;
mod element;
pub mod error;
mod molecule;
pub mod pdb;

pub use atom::Atom;
pub use bond::{Bond, BondOrder};
pub use element::Element;
pub use error::BiologyError;
pub use molecule::Molecule;
