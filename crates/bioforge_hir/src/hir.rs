//! Typed HIR (High-level Intermediate Representation) nodes.
//!
//! These nodes mirror the AST structure but with resolved types:
//! - Unit strings → [`Quantity`] values with real dimensions
//! - Entity references → validated against declarations
//! - Properties → validated against expected dimensions
//!
//! The HIR is independent of source syntax and ready for execution.

use bioforge_ast::EntityKind;
use bioforge_types::Quantity;
use std::fmt;

/// A fully validated, typed BioForge program.
#[derive(Debug, Clone, PartialEq)]
pub struct HirProgram {
    /// The validated experiments.
    pub experiments: Vec<HirExperiment>,
}

/// A validated experiment with resolved types.
#[derive(Debug, Clone, PartialEq)]
pub struct HirExperiment {
    /// Experiment name.
    pub name: String,
    /// Declared biological entities.
    pub entities: Vec<HirEntity>,
    /// Validated environment configuration (if present).
    pub environment: Option<HirEnvironment>,
    /// Validated simulation parameters (if present).
    pub simulation: Option<HirSimulation>,
    /// Measurement requests with validated references.
    pub measurements: Vec<HirMeasurement>,
    /// Visualization targets with validated references.
    pub visualizations: Vec<HirVisualization>,
}

/// A resolved biological entity declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct HirEntity {
    /// The biological entity kind (protein, ligand, atom, etc.).
    pub kind: EntityKind,
    /// The entity name (used for referencing).
    pub name: String,
    /// The resolved initializer expression.
    pub initializer: HirExpr,
}

/// Validated environment configuration.
///
/// Properties are validated against expected physical dimensions:
/// - `temperature` must be Temperature
/// - `pressure` must be Pressure
/// - `pH` must be a dimensionless number in [0.0, 14.0]
#[derive(Debug, Clone, PartialEq)]
pub struct HirEnvironment {
    /// Temperature setting (must have Temperature dimension).
    pub temperature: Option<Quantity>,
    /// pH value (dimensionless, range 0.0–14.0).
    pub ph: Option<f64>,
    /// Pressure setting (must have Pressure dimension).
    pub pressure: Option<Quantity>,
}

/// Validated simulation configuration.
///
/// Both `timestep` and `duration` must have Time dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct HirSimulation {
    /// Integration timestep (Time dimension).
    pub timestep: Quantity,
    /// Total simulation duration (Time dimension).
    pub duration: Quantity,
}

/// A resolved measurement request.
#[derive(Debug, Clone, PartialEq)]
pub struct HirMeasurement {
    /// The measurement function name (e.g., "distance", "energy").
    pub function: String,
    /// Arguments — entity names that have been validated to exist.
    pub args: Vec<String>,
}

/// A resolved visualization target.
#[derive(Debug, Clone, PartialEq)]
pub struct HirVisualization {
    /// The entity name to visualize (validated to exist).
    pub target: String,
}

/// Expression with resolved types.
///
/// Unlike AST `Expr`, `HirExpr` has all units resolved to real
/// [`Quantity`] values with dimensional information.
#[derive(Debug, Clone, PartialEq)]
pub enum HirExpr {
    /// A resolved physical quantity with dimension.
    Quantity(Quantity),
    /// A string value.
    String(String),
    /// A boolean value.
    Bool(bool),
    /// A bare number (dimensionless).
    Number(f64),
    /// A function call with resolved arguments.
    FunctionCall {
        /// Function name.
        name: String,
        /// Resolved arguments.
        args: Vec<HirExpr>,
    },
    /// A reference to a declared entity (validated to exist).
    EntityRef(String),
}

// ─── Display implementations for debugging ──────────────────────────────────

impl fmt::Display for HirExperiment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "experiment {} {{", self.name)?;
        for entity in &self.entities {
            writeln!(f, "  {} {} = {}", entity.kind, entity.name, entity.initializer)?;
        }
        if let Some(env) = &self.environment {
            writeln!(f, "  environment {{")?;
            if let Some(t) = &env.temperature {
                writeln!(f, "    temperature = {}", t)?;
            }
            if let Some(ph) = &env.ph {
                writeln!(f, "    pH = {}", ph)?;
            }
            if let Some(p) = &env.pressure {
                writeln!(f, "    pressure = {}", p)?;
            }
            writeln!(f, "  }}")?;
        }
        if let Some(sim) = &self.simulation {
            writeln!(f, "  simulate {{")?;
            writeln!(f, "    timestep = {}", sim.timestep)?;
            writeln!(f, "    duration = {}", sim.duration)?;
            writeln!(f, "  }}")?;
        }
        for m in &self.measurements {
            writeln!(f, "  measure: {}({})", m.function, m.args.join(", "))?;
        }
        for v in &self.visualizations {
            writeln!(f, "  visualize: {}", v.target)?;
        }
        write!(f, "}}")
    }
}

impl fmt::Display for HirExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Quantity(q) => write!(f, "{}", q),
            Self::String(s) => write!(f, "\"{}\"", s),
            Self::Bool(b) => write!(f, "{}", b),
            Self::Number(n) => {
                if *n == (*n as i64) as f64 {
                    write!(f, "{}", *n as i64)
                } else {
                    write!(f, "{}", n)
                }
            }
            Self::FunctionCall { name, args } => {
                let arg_strs: Vec<std::string::String> =
                    args.iter().map(|a| format!("{}", a)).collect();
                write!(f, "{}({})", name, arg_strs.join(", "))
            }
            Self::EntityRef(name) => write!(f, "{}", name),
        }
    }
}
