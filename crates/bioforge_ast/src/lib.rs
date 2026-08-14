//! # BioForge AST
//!
//! Abstract Syntax Tree definitions for the BioLang language.
//!
//! Every AST node carries source location information via [`Spanned<T>`].

use bioforge_diagnostics::Spanned;
use std::fmt;

// ─── Program ───────────────────────────────────────────────────────────────────

/// A complete BioLang program consisting of one or more experiments.
#[derive(Debug, Clone, PartialEq)]
pub struct Program {
    pub experiments: Vec<Spanned<ExperimentDecl>>,
}

// ─── Experiment ────────────────────────────────────────────────────────────────

/// An experiment declaration: `experiment Name { ... }`
#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentDecl {
    pub name: Spanned<String>,
    pub body: Vec<Spanned<Statement>>,
}

// ─── Statements ────────────────────────────────────────────────────────────────

/// Statements within an experiment body.
#[derive(Debug, Clone, PartialEq)]
pub enum Statement {
    /// Entity declaration: `protein name = expr`
    EntityDecl {
        kind: Spanned<EntityKind>,
        name: Spanned<String>,
        initializer: Spanned<Expr>,
    },
    /// Environment block: `environment { ... }` or `environment name { ... }`
    EnvironmentBlock {
        name: Option<Spanned<String>>,
        properties: Vec<Spanned<Property>>,
    },
    /// Simulate block: `simulate { ... }`
    SimulateBlock {
        properties: Vec<Spanned<Property>>,
    },
    /// Measure block: `measure { ... }`
    MeasureBlock {
        measurements: Vec<Spanned<Expr>>,
    },
    /// Visualize block: `visualize { ... }`
    VisualizeBlock {
        targets: Vec<Spanned<Expr>>,
    },
    /// Simple assignment: `name = expr`
    Assignment {
        name: Spanned<String>,
        value: Spanned<Expr>,
    },
}

// ─── Entity Kind ───────────────────────────────────────────────────────────────

/// Biological entity types supported in BioLang v0.1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityKind {
    Protein,
    Ligand,
    Ion,
    Molecule,
    Atom,
}

impl fmt::Display for EntityKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EntityKind::Protein => write!(f, "protein"),
            EntityKind::Ligand => write!(f, "ligand"),
            EntityKind::Ion => write!(f, "ion"),
            EntityKind::Molecule => write!(f, "molecule"),
            EntityKind::Atom => write!(f, "atom"),
        }
    }
}

// ─── Expressions ───────────────────────────────────────────────────────────────

/// Expressions in BioLang.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    /// Numeric literal: `310`, `7.4`
    NumberLiteral(f64),
    /// String literal: `"protein.pdb"`
    StringLiteral(String),
    /// Boolean literal: `true`, `false`
    BoolLiteral(bool),
    /// Quantity with unit: `310 K`, `1 fs`, `10 nm`
    Quantity {
        value: f64,
        unit: Spanned<String>,
    },
    /// Identifier reference: `receptor`, `drug`
    Identifier(String),
    /// Function call: `load_structure("protein.pdb")`
    FunctionCall {
        name: Spanned<String>,
        args: Vec<Spanned<Expr>>,
    },
    /// Member access: `receptor.chain_A`
    MemberAccess {
        object: Box<Spanned<Expr>>,
        member: Spanned<String>,
    },
    /// Binary operation: `a + b`
    BinaryOp {
        left: Box<Spanned<Expr>>,
        op: Spanned<BinOp>,
        right: Box<Spanned<Expr>>,
    },
}

// ─── Binary Operators ──────────────────────────────────────────────────────────

/// Binary operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
}

impl fmt::Display for BinOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BinOp::Add => write!(f, "+"),
            BinOp::Sub => write!(f, "-"),
            BinOp::Mul => write!(f, "*"),
            BinOp::Div => write!(f, "/"),
        }
    }
}

// ─── Property ──────────────────────────────────────────────────────────────────

/// A key-value property: `temperature = 310 K`
#[derive(Debug, Clone, PartialEq)]
pub struct Property {
    pub key: Spanned<String>,
    pub value: Spanned<Expr>,
}

// ─── Pretty Printer ────────────────────────────────────────────────────────────

impl Program {
    /// Produce a human-readable representation of the AST for debugging.
    pub fn pretty_print(&self) -> String {
        let mut out = String::new();
        for (i, exp) in self.experiments.iter().enumerate() {
            if i > 0 {
                out.push('\n');
            }
            pretty_experiment(&mut out, &exp.node, 0);
        }
        out
    }
}

fn pretty_experiment(out: &mut String, exp: &ExperimentDecl, indent: usize) {
    let pad = "  ".repeat(indent);
    out.push_str(&format!("{}experiment {} {{\n", pad, exp.name.node));
    for stmt in &exp.body {
        pretty_statement(out, &stmt.node, indent + 1);
    }
    out.push_str(&format!("{}}}\n", pad));
}

fn pretty_statement(out: &mut String, stmt: &Statement, indent: usize) {
    let pad = "  ".repeat(indent);
    match stmt {
        Statement::EntityDecl { kind, name, initializer } => {
            out.push_str(&format!(
                "{}{} {} = {}\n",
                pad, kind.node, name.node,
                pretty_expr(&initializer.node)
            ));
        }
        Statement::EnvironmentBlock { name, properties } => {
            if let Some(n) = name {
                out.push_str(&format!("{}environment {} {{\n", pad, n.node));
            } else {
                out.push_str(&format!("{}environment {{\n", pad));
            }
            for prop in properties {
                out.push_str(&format!(
                    "{}  {} = {}\n",
                    pad, prop.node.key.node,
                    pretty_expr(&prop.node.value.node)
                ));
            }
            out.push_str(&format!("{}}}\n", pad));
        }
        Statement::SimulateBlock { properties } => {
            out.push_str(&format!("{}simulate {{\n", pad));
            for prop in properties {
                out.push_str(&format!(
                    "{}  {} = {}\n",
                    pad, prop.node.key.node,
                    pretty_expr(&prop.node.value.node)
                ));
            }
            out.push_str(&format!("{}}}\n", pad));
        }
        Statement::MeasureBlock { measurements } => {
            out.push_str(&format!("{}measure {{\n", pad));
            for m in measurements {
                out.push_str(&format!("{}  {}\n", pad, pretty_expr(&m.node)));
            }
            out.push_str(&format!("{}}}\n", pad));
        }
        Statement::VisualizeBlock { targets } => {
            out.push_str(&format!("{}visualize {{\n", pad));
            for t in targets {
                out.push_str(&format!("{}  {}\n", pad, pretty_expr(&t.node)));
            }
            out.push_str(&format!("{}}}\n", pad));
        }
        Statement::Assignment { name, value } => {
            out.push_str(&format!(
                "{}{} = {}\n",
                pad, name.node,
                pretty_expr(&value.node)
            ));
        }
    }
}

fn pretty_expr(expr: &Expr) -> String {
    match expr {
        Expr::NumberLiteral(n) => {
            if *n == (*n as i64) as f64 {
                format!("{}", *n as i64)
            } else {
                format!("{}", n)
            }
        }
        Expr::StringLiteral(s) => format!("\"{}\"", s),
        Expr::BoolLiteral(b) => format!("{}", b),
        Expr::Quantity { value, unit } => {
            if *value == (*value as i64) as f64 {
                format!("{} {}", *value as i64, unit.node)
            } else {
                format!("{} {}", value, unit.node)
            }
        }
        Expr::Identifier(name) => name.clone(),
        Expr::FunctionCall { name, args } => {
            let arg_strs: Vec<String> = args.iter().map(|a| pretty_expr(&a.node)).collect();
            format!("{}({})", name.node, arg_strs.join(", "))
        }
        Expr::MemberAccess { object, member } => {
            format!("{}.{}", pretty_expr(&object.node), member.node)
        }
        Expr::BinaryOp { left, op, right } => {
            format!(
                "({} {} {})",
                pretty_expr(&left.node),
                op.node,
                pretty_expr(&right.node)
            )
        }
    }
}

// ─── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_diagnostics::Span;

    fn s(start: usize, end: usize) -> Span {
        Span::new(start, end)
    }

    #[test]
    fn test_entity_kind_display() {
        assert_eq!(EntityKind::Protein.to_string(), "protein");
        assert_eq!(EntityKind::Ligand.to_string(), "ligand");
        assert_eq!(EntityKind::Ion.to_string(), "ion");
        assert_eq!(EntityKind::Molecule.to_string(), "molecule");
        assert_eq!(EntityKind::Atom.to_string(), "atom");
    }

    #[test]
    fn test_binop_display() {
        assert_eq!(BinOp::Add.to_string(), "+");
        assert_eq!(BinOp::Sub.to_string(), "-");
        assert_eq!(BinOp::Mul.to_string(), "*");
        assert_eq!(BinOp::Div.to_string(), "/");
    }

    #[test]
    fn test_pretty_print_empty_program() {
        let program = Program { experiments: vec![] };
        assert_eq!(program.pretty_print(), "");
    }

    #[test]
    fn test_pretty_print_simple_experiment() {
        let program = Program {
            experiments: vec![Spanned::new(
                ExperimentDecl {
                    name: Spanned::new("Test".to_string(), s(0, 4)),
                    body: vec![Spanned::new(
                        Statement::EntityDecl {
                            kind: Spanned::new(EntityKind::Atom, s(0, 4)),
                            name: Spanned::new("h".to_string(), s(5, 6)),
                            initializer: Spanned::new(
                                Expr::FunctionCall {
                                    name: Spanned::new("atom".to_string(), s(9, 13)),
                                    args: vec![Spanned::new(
                                        Expr::StringLiteral("H".to_string()),
                                        s(14, 17),
                                    )],
                                },
                                s(9, 18),
                            ),
                        },
                        s(0, 18),
                    )],
                },
                s(0, 20),
            )],
        };
        let output = program.pretty_print();
        assert!(output.contains("experiment Test"));
        assert!(output.contains("atom h = atom(\"H\")"));
    }
}
