//! # BioForge HIR (High-level Intermediate Representation)
//!
//! This crate provides semantic analysis for BioForge programs.
//! It lowers the syntactic AST into a typed, validated intermediate
//! representation where:
//!
//! - All unit strings are resolved to [`bioforge_types::Quantity`] values
//! - Environment properties are validated against expected dimensions
//! - Entity references are checked against declarations
//! - Simulation parameters are dimensionally verified
//!
//! ## Pipeline Position
//!
//! ```text
//! Source → Lexer → Parser → AST → [HIR Lowering] → BioIR → Runtime
//! ```
//!
//! The AST represents **what the user wrote**. The HIR represents
//! **what it means** — resolved, validated, and typed.

#![deny(unsafe_code)]

mod error;
mod hir;
mod lower;

pub use error::SemanticError;
pub use hir::*;
pub use lower::lower;
