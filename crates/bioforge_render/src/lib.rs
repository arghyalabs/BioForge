//! # BioForge Render (`bioforge_render`)
//!
//! 3D visual observation layer and molecular representation engine for BioForge.
//!
//! ## Scientific Architecture (Principle 2 & Invariant §9)
//!
//! The renderer operates as a **pure read-only observer**. It extracts 3D meshes, colors,
//! and transformation buffers from immutable references `&SimulationState`. The renderer
//! **NEVER owns state** and **NEVER determines biological behavior or physics**.
//!
//! ## Core Modules
//!
//! - [`Color`]: RGBA color model and standard CPK elemental palette.
//! - [`Mesh`]: 3D triangle mesh buffers and procedural generators (spheres, cylinders).
//! - [`Camera`]: 3D orbit camera with view and perspective projection matrix math.
//! - [`RenderStyle`]: Visual styles (Space-Filling, Ball-and-Stick, Backbone Trace).
//! - [`Scene`]: High-level scene representation with Wavefront `.obj` export.

#![deny(unsafe_code)]

pub mod camera;
pub mod color;
pub mod error;
pub mod mesh;
pub mod scene;
pub mod style;

pub use camera::Camera;
pub use color::{color_by_chain, cpk_color_for_element, Color};
pub use error::RenderError;
pub use mesh::{generate_cylinder, generate_sphere, generate_split_cylinder, Mesh, Vertex};
pub use scene::Scene;
pub use style::RenderStyle;
