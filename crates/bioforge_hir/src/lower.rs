//! AST → HIR lowering pass.
//!
//! This module performs semantic analysis on the parsed AST, producing
//! a validated, typed HIR (High-level Intermediate Representation).
//!
//! ## Validations Performed
//!
//! 1. **Identifier resolution** — all entity references must be declared
//! 2. **Unit resolution** — unit strings → [`Quantity`] with real dimensions
//! 3. **Property validation** — environment/simulate properties checked
//! 4. **Duplicate detection** — entity names must be unique per experiment
//! 5. **Dimensional checking** — properties have correct physical dimensions

use std::collections::HashMap;

use bioforge_ast::{self as ast, Expr, Statement};
use bioforge_diagnostics::{Span, Spanned};
use bioforge_types::{Dimension, Quantity, UnitRegistry};

use crate::error::SemanticError;
use crate::hir::*;

/// Lower an AST program into a validated HIR program.
///
/// Returns the HIR program and a list of semantic errors. If there are
/// errors, the HIR may be partially constructed (best-effort).
///
/// # Examples
///
/// ```text
/// let (program, errors) = parse(source);
/// let (hir, semantic_errors) = lower(&program);
/// if semantic_errors.is_empty() {
///     // HIR is fully valid
/// }
/// ```
pub fn lower(program: &ast::Program) -> (HirProgram, Vec<SemanticError>) {
    let mut lowerer = Lowerer::new();
    let hir = lowerer.lower_program(program);
    (hir, lowerer.errors)
}

/// Internal lowering state.
struct Lowerer {
    /// Accumulated semantic errors.
    errors: Vec<SemanticError>,
    /// Unit registry for resolving unit strings.
    registry: UnitRegistry,
}

impl Lowerer {
    fn new() -> Self {
        Self {
            errors: Vec::new(),
            registry: UnitRegistry::new(),
        }
    }

    fn lower_program(&mut self, program: &ast::Program) -> HirProgram {
        let experiments = program
            .experiments
            .iter()
            .map(|exp| self.lower_experiment(&exp.node))
            .collect();
        HirProgram { experiments }
    }

    fn lower_experiment(&mut self, exp: &ast::ExperimentDecl) -> HirExperiment {
        // First pass: collect all entity declarations for reference checking
        let mut entity_spans: HashMap<String, Span> = HashMap::new();
        let mut entities = Vec::new();
        let mut environment = None;
        let mut simulation = None;
        let mut measurements = Vec::new();
        let mut visualizations = Vec::new();

        for stmt in &exp.body {
            match &stmt.node {
                Statement::EntityDecl {
                    kind,
                    name,
                    initializer,
                } => {
                    // Check for duplicates
                    if let Some(first_span) = entity_spans.get(&name.node) {
                        self.errors.push(SemanticError::DuplicateEntity {
                            name: name.node.clone(),
                            first: *first_span,
                            second: name.span,
                        });
                    } else {
                        entity_spans.insert(name.node.clone(), name.span);
                    }

                    entities.push(HirEntity {
                        kind: kind.node,
                        name: name.node.clone(),
                        initializer: self.lower_expr(&initializer.node, initializer.span),
                    });
                }
                Statement::EnvironmentBlock { properties, .. } => {
                    environment = Some(self.lower_environment(properties));
                }
                Statement::SimulateBlock { properties } => {
                    simulation = self.lower_simulation(properties, stmt.span);
                }
                Statement::MeasureBlock {
                    measurements: measure_exprs,
                } => {
                    for expr in measure_exprs {
                        if let Some(m) = self.lower_measurement(&expr.node, expr.span, &entity_spans)
                        {
                            measurements.push(m);
                        }
                    }
                }
                Statement::VisualizeBlock { targets } => {
                    for target in targets {
                        if let Some(v) =
                            self.lower_visualization(&target.node, target.span, &entity_spans)
                        {
                            visualizations.push(v);
                        }
                    }
                }
                Statement::Assignment { .. } => {
                    // Assignments are not yet lowered to HIR
                    // (future: variable bindings)
                }
            }
        }

        HirExperiment {
            name: exp.name.node.clone(),
            entities,
            environment,
            simulation,
            measurements,
            visualizations,
        }
    }

    // ─── Environment lowering ───────────────────────────────────────────

    fn lower_environment(
        &mut self,
        properties: &[Spanned<ast::Property>],
    ) -> HirEnvironment {
        let mut env = HirEnvironment {
            temperature: None,
            ph: None,
            pressure: None,
        };

        for prop in properties {
            let key = &prop.node.key.node;
            let value = &prop.node.value;

            match key.as_str() {
                "temperature" => {
                    if let Some(q) = self.resolve_quantity_property(
                        &value.node,
                        value.span,
                        "temperature",
                        &Dimension::temperature(),
                    ) {
                        env.temperature = Some(q);
                    }
                }
                "pH" => {
                    env.ph = self.resolve_ph(&value.node, value.span);
                }
                "pressure" => {
                    if let Some(q) = self.resolve_quantity_property(
                        &value.node,
                        value.span,
                        "pressure",
                        &Dimension::pressure(),
                    ) {
                        env.pressure = Some(q);
                    }
                }
                other => {
                    self.errors.push(SemanticError::UnknownProperty {
                        block: "environment".to_string(),
                        property: other.to_string(),
                        span: prop.node.key.span,
                    });
                }
            }
        }

        env
    }

    /// Resolve pH: must be a dimensionless number in [0.0, 14.0].
    fn resolve_ph(&mut self, expr: &Expr, span: Span) -> Option<f64> {
        match expr {
            Expr::NumberLiteral(n) => {
                if !(0.0..=14.0).contains(n) {
                    self.errors.push(SemanticError::PhOutOfRange {
                        value: *n,
                        span,
                    });
                    return None;
                }
                Some(*n)
            }
            Expr::Quantity { .. } => {
                self.errors.push(SemanticError::InvalidDimension {
                    property: "pH".to_string(),
                    expected: "dimensionless number".to_string(),
                    got: "quantity with unit".to_string(),
                    span,
                });
                None
            }
            _ => None, // Silently skip non-numeric (caught elsewhere)
        }
    }

    // ─── Simulation lowering ────────────────────────────────────────────

    fn lower_simulation(
        &mut self,
        properties: &[Spanned<ast::Property>],
        block_span: Span,
    ) -> Option<HirSimulation> {
        let mut timestep: Option<Quantity> = None;
        let mut duration: Option<Quantity> = None;

        for prop in properties {
            let key = &prop.node.key.node;
            let value = &prop.node.value;

            match key.as_str() {
                "timestep" => {
                    timestep = self.resolve_quantity_property(
                        &value.node,
                        value.span,
                        "timestep",
                        &Dimension::time(),
                    );
                }
                "duration" => {
                    duration = self.resolve_quantity_property(
                        &value.node,
                        value.span,
                        "duration",
                        &Dimension::time(),
                    );
                }
                other => {
                    self.errors.push(SemanticError::UnknownProperty {
                        block: "simulate".to_string(),
                        property: other.to_string(),
                        span: prop.node.key.span,
                    });
                }
            }
        }

        // Both timestep and duration are required
        let ts = match timestep {
            Some(t) => t,
            None => {
                self.errors.push(SemanticError::MissingRequiredProperty {
                    block: "simulate".to_string(),
                    property: "timestep".to_string(),
                    span: block_span,
                });
                return None;
            }
        };

        let dur = match duration {
            Some(d) => d,
            None => {
                self.errors.push(SemanticError::MissingRequiredProperty {
                    block: "simulate".to_string(),
                    property: "duration".to_string(),
                    span: block_span,
                });
                return None;
            }
        };

        Some(HirSimulation {
            timestep: ts,
            duration: dur,
        })
    }

    // ─── Measurement lowering ───────────────────────────────────────────

    fn lower_measurement(
        &mut self,
        expr: &Expr,
        span: Span,
        entities: &HashMap<String, Span>,
    ) -> Option<HirMeasurement> {
        if let Expr::FunctionCall { name, args } = expr {
            let mut arg_names = Vec::new();
            for arg in args {
                if let Expr::Identifier(ref id) = arg.node {
                    if !entities.contains_key(id) {
                        self.errors.push(SemanticError::UndeclaredEntity {
                            name: id.clone(),
                            span: arg.span,
                        });
                    }
                    arg_names.push(id.clone());
                }
            }
            Some(HirMeasurement {
                function: name.node.clone(),
                args: arg_names,
            })
        } else {
            // Non-function-call measurements are not supported yet
            let _ = span;
            None
        }
    }

    // ─── Visualization lowering ─────────────────────────────────────────

    fn lower_visualization(
        &mut self,
        expr: &Expr,
        span: Span,
        entities: &HashMap<String, Span>,
    ) -> Option<HirVisualization> {
        if let Expr::Identifier(ref name) = expr {
            if !entities.contains_key(name) {
                self.errors.push(SemanticError::UndeclaredEntity {
                    name: name.clone(),
                    span,
                });
            }
            Some(HirVisualization {
                target: name.clone(),
            })
        } else {
            None
        }
    }

    // ─── Expression lowering ────────────────────────────────────────────

    fn lower_expr(&mut self, expr: &Expr, span: Span) -> HirExpr {
        match expr {
            Expr::NumberLiteral(n) => HirExpr::Number(*n),
            Expr::StringLiteral(s) => HirExpr::String(s.clone()),
            Expr::BoolLiteral(b) => HirExpr::Bool(*b),
            Expr::Quantity { value, unit } => {
                if let Some(u) = self.registry.resolve(&unit.node) {
                    HirExpr::Quantity(Quantity::new(*value, u.clone()))
                } else {
                    self.errors.push(SemanticError::UnknownUnit {
                        unit: unit.node.clone(),
                        span: unit.span,
                    });
                    // Fallback to bare number
                    HirExpr::Number(*value)
                }
            }
            Expr::Identifier(name) => HirExpr::EntityRef(name.clone()),
            Expr::FunctionCall { name, args } => {
                let resolved_args: Vec<HirExpr> = args
                    .iter()
                    .map(|a| self.lower_expr(&a.node, a.span))
                    .collect();
                HirExpr::FunctionCall {
                    name: name.node.clone(),
                    args: resolved_args,
                }
            }
            Expr::MemberAccess { .. } => {
                // Member access not yet supported in HIR
                HirExpr::Number(0.0)
            }
            Expr::BinaryOp { .. } => {
                // Binary operations not yet supported in HIR
                // (future: evaluate constant expressions)
                let _ = span;
                HirExpr::Number(0.0)
            }
        }
    }

    // ─── Helpers ────────────────────────────────────────────────────────

    /// Resolve an expression as a Quantity and validate its dimension.
    fn resolve_quantity_property(
        &mut self,
        expr: &Expr,
        span: Span,
        property: &str,
        expected_dim: &Dimension,
    ) -> Option<Quantity> {
        match expr {
            Expr::Quantity { value, unit } => {
                if let Some(u) = self.registry.resolve(&unit.node) {
                    if !u.dimension.is_compatible(expected_dim) {
                        self.errors.push(SemanticError::InvalidDimension {
                            property: property.to_string(),
                            expected: expected_dim.name(),
                            got: u.dimension.name(),
                            span,
                        });
                        return None;
                    }
                    Some(Quantity::new(*value, u.clone()))
                } else {
                    self.errors.push(SemanticError::UnknownUnit {
                        unit: unit.node.clone(),
                        span: unit.span,
                    });
                    None
                }
            }
            _ => {
                // Non-quantity where quantity expected
                self.errors.push(SemanticError::InvalidDimension {
                    property: property.to_string(),
                    expected: expected_dim.name(),
                    got: "non-quantity value".to_string(),
                    span,
                });
                None
            }
        }
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bioforge_ast::*;
    use bioforge_diagnostics::Span;

    fn s(start: usize, end: usize) -> Span {
        Span::new(start, end)
    }

    fn sp<T>(node: T, start: usize, end: usize) -> Spanned<T> {
        Spanned::new(node, s(start, end))
    }

    /// Build a minimal valid program for testing.
    fn hello_bio_ast() -> Program {
        Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("HelloBiology".to_string(), 11, 23),
                    body: vec![
                        // protein receptor = load_structure("protein.pdb")
                        sp(
                            Statement::EntityDecl {
                                kind: sp(EntityKind::Protein, 30, 37),
                                name: sp("receptor".to_string(), 38, 46),
                                initializer: sp(
                                    Expr::FunctionCall {
                                        name: sp("load_structure".to_string(), 49, 63),
                                        args: vec![sp(
                                            Expr::StringLiteral("protein.pdb".to_string()),
                                            64,
                                            77,
                                        )],
                                    },
                                    49,
                                    78,
                                ),
                            },
                            30,
                            78,
                        ),
                        // ligand drug = load("drug.sdf")
                        sp(
                            Statement::EntityDecl {
                                kind: sp(EntityKind::Ligand, 83, 89),
                                name: sp("drug".to_string(), 90, 94),
                                initializer: sp(
                                    Expr::FunctionCall {
                                        name: sp("load".to_string(), 97, 101),
                                        args: vec![sp(
                                            Expr::StringLiteral("drug.sdf".to_string()),
                                            102,
                                            112,
                                        )],
                                    },
                                    97,
                                    113,
                                ),
                            },
                            83,
                            113,
                        ),
                        // environment { temperature = 310 K, pH = 7.4 }
                        sp(
                            Statement::EnvironmentBlock {
                                name: None,
                                properties: vec![
                                    sp(
                                        Property {
                                            key: sp("temperature".to_string(), 130, 141),
                                            value: sp(
                                                Expr::Quantity {
                                                    value: 310.0,
                                                    unit: sp("K".to_string(), 148, 149),
                                                },
                                                144,
                                                149,
                                            ),
                                        },
                                        130,
                                        149,
                                    ),
                                    sp(
                                        Property {
                                            key: sp("pH".to_string(), 154, 156),
                                            value: sp(Expr::NumberLiteral(7.4), 159, 162),
                                        },
                                        154,
                                        162,
                                    ),
                                ],
                            },
                            118,
                            168,
                        ),
                        // simulate { timestep = 1 fs, duration = 10 ps }
                        sp(
                            Statement::SimulateBlock {
                                properties: vec![
                                    sp(
                                        Property {
                                            key: sp("timestep".to_string(), 185, 193),
                                            value: sp(
                                                Expr::Quantity {
                                                    value: 1.0,
                                                    unit: sp("fs".to_string(), 198, 200),
                                                },
                                                196,
                                                200,
                                            ),
                                        },
                                        185,
                                        200,
                                    ),
                                    sp(
                                        Property {
                                            key: sp("duration".to_string(), 205, 213),
                                            value: sp(
                                                Expr::Quantity {
                                                    value: 10.0,
                                                    unit: sp("ps".to_string(), 219, 221),
                                                },
                                                216,
                                                221,
                                            ),
                                        },
                                        205,
                                        221,
                                    ),
                                ],
                            },
                            173,
                            227,
                        ),
                        // measure { distance(receptor, drug), energy(receptor, drug) }
                        sp(
                            Statement::MeasureBlock {
                                measurements: vec![
                                    sp(
                                        Expr::FunctionCall {
                                            name: sp("distance".to_string(), 244, 252),
                                            args: vec![
                                                sp(Expr::Identifier("receptor".to_string()), 253, 261),
                                                sp(Expr::Identifier("drug".to_string()), 263, 267),
                                            ],
                                        },
                                        244,
                                        268,
                                    ),
                                    sp(
                                        Expr::FunctionCall {
                                            name: sp("energy".to_string(), 273, 279),
                                            args: vec![
                                                sp(Expr::Identifier("receptor".to_string()), 280, 288),
                                                sp(Expr::Identifier("drug".to_string()), 290, 294),
                                            ],
                                        },
                                        273,
                                        295,
                                    ),
                                ],
                            },
                            232,
                            301,
                        ),
                        // visualize { receptor, drug }
                        sp(
                            Statement::VisualizeBlock {
                                targets: vec![
                                    sp(Expr::Identifier("receptor".to_string()), 318, 326),
                                    sp(Expr::Identifier("drug".to_string()), 331, 335),
                                ],
                            },
                            306,
                            341,
                        ),
                    ],
                },
                0,
                343,
            )],
        }
    }

    #[test]
    fn test_lower_valid_program() {
        let program = hello_bio_ast();
        let (hir, errors) = lower(&program);

        assert!(errors.is_empty(), "expected no errors, got: {:?}", errors);
        assert_eq!(hir.experiments.len(), 1);

        let exp = &hir.experiments[0];
        assert_eq!(exp.name, "HelloBiology");
        assert_eq!(exp.entities.len(), 2);
        assert_eq!(exp.entities[0].kind, EntityKind::Protein);
        assert_eq!(exp.entities[0].name, "receptor");
        assert_eq!(exp.entities[1].kind, EntityKind::Ligand);
        assert_eq!(exp.entities[1].name, "drug");
    }

    #[test]
    fn test_environment_validation() {
        let program = hello_bio_ast();
        let (hir, errors) = lower(&program);

        assert!(errors.is_empty());
        let env = hir.experiments[0].environment.as_ref().unwrap();
        assert!(env.temperature.is_some());
        let temp = env.temperature.as_ref().unwrap();
        assert!((temp.value - 310.0).abs() < 1e-10);
        assert_eq!(temp.unit.name, "K");
        assert!((env.ph.unwrap() - 7.4).abs() < 1e-10);
    }

    #[test]
    fn test_simulation_validation() {
        let program = hello_bio_ast();
        let (hir, errors) = lower(&program);

        assert!(errors.is_empty());
        let sim = hir.experiments[0].simulation.as_ref().unwrap();
        assert!((sim.timestep.value - 1.0).abs() < 1e-10);
        assert_eq!(sim.timestep.unit.name, "fs");
        assert!((sim.duration.value - 10.0).abs() < 1e-10);
        assert_eq!(sim.duration.unit.name, "ps");
    }

    #[test]
    fn test_measurement_validation() {
        let program = hello_bio_ast();
        let (hir, errors) = lower(&program);

        assert!(errors.is_empty());
        let exp = &hir.experiments[0];
        assert_eq!(exp.measurements.len(), 2);
        assert_eq!(exp.measurements[0].function, "distance");
        assert_eq!(exp.measurements[0].args, vec!["receptor", "drug"]);
        assert_eq!(exp.measurements[1].function, "energy");
    }

    #[test]
    fn test_visualization_validation() {
        let program = hello_bio_ast();
        let (hir, errors) = lower(&program);

        assert!(errors.is_empty());
        let exp = &hir.experiments[0];
        assert_eq!(exp.visualizations.len(), 2);
        assert_eq!(exp.visualizations[0].target, "receptor");
        assert_eq!(exp.visualizations[1].target, "drug");
    }

    #[test]
    fn test_wrong_dimension_temperature() {
        // temperature = 310 nm (Length instead of Temperature)
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![sp(
                        Statement::EnvironmentBlock {
                            name: None,
                            properties: vec![sp(
                                Property {
                                    key: sp("temperature".to_string(), 10, 21),
                                    value: sp(
                                        Expr::Quantity {
                                            value: 310.0,
                                            unit: sp("nm".to_string(), 28, 30),
                                        },
                                        24,
                                        30,
                                    ),
                                },
                                10,
                                30,
                            )],
                        },
                        5,
                        35,
                    )],
                },
                0,
                40,
            )],
        };

        let (_hir, errors) = lower(&program);
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            &errors[0],
            SemanticError::InvalidDimension {
                property,
                expected,
                got,
                ..
            } if property == "temperature" && expected == "Temperature" && got == "Length"
        ));
    }

    #[test]
    fn test_ph_out_of_range() {
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![sp(
                        Statement::EnvironmentBlock {
                            name: None,
                            properties: vec![sp(
                                Property {
                                    key: sp("pH".to_string(), 10, 12),
                                    value: sp(Expr::NumberLiteral(15.0), 15, 19),
                                },
                                10,
                                19,
                            )],
                        },
                        5,
                        25,
                    )],
                },
                0,
                30,
            )],
        };

        let (_hir, errors) = lower(&program);
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            &errors[0],
            SemanticError::PhOutOfRange { value, .. } if (*value - 15.0).abs() < 1e-10
        ));
    }

    #[test]
    fn test_undeclared_entity_in_visualize() {
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![sp(
                        Statement::VisualizeBlock {
                            targets: vec![sp(
                                Expr::Identifier("nonexistent".to_string()),
                                10,
                                21,
                            )],
                        },
                        5,
                        25,
                    )],
                },
                0,
                30,
            )],
        };

        let (_hir, errors) = lower(&program);
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            &errors[0],
            SemanticError::UndeclaredEntity { name, .. } if name == "nonexistent"
        ));
    }

    #[test]
    fn test_duplicate_entity() {
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![
                        sp(
                            Statement::EntityDecl {
                                kind: sp(EntityKind::Atom, 5, 9),
                                name: sp("h".to_string(), 10, 11),
                                initializer: sp(
                                    Expr::FunctionCall {
                                        name: sp("atom".to_string(), 14, 18),
                                        args: vec![sp(
                                            Expr::StringLiteral("H".to_string()),
                                            19,
                                            22,
                                        )],
                                    },
                                    14,
                                    23,
                                ),
                            },
                            5,
                            23,
                        ),
                        sp(
                            Statement::EntityDecl {
                                kind: sp(EntityKind::Atom, 28, 32),
                                name: sp("h".to_string(), 33, 34),
                                initializer: sp(
                                    Expr::FunctionCall {
                                        name: sp("atom".to_string(), 37, 41),
                                        args: vec![sp(
                                            Expr::StringLiteral("He".to_string()),
                                            42,
                                            46,
                                        )],
                                    },
                                    37,
                                    47,
                                ),
                            },
                            28,
                            47,
                        ),
                    ],
                },
                0,
                50,
            )],
        };

        let (_hir, errors) = lower(&program);
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            &errors[0],
            SemanticError::DuplicateEntity { name, .. } if name == "h"
        ));
    }

    #[test]
    fn test_missing_simulate_properties() {
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![sp(
                        Statement::SimulateBlock {
                            properties: vec![], // no timestep or duration
                        },
                        5,
                        20,
                    )],
                },
                0,
                25,
            )],
        };

        let (_hir, errors) = lower(&program);
        // Should get errors for both missing timestep and duration
        assert!(errors.len() >= 1);
        assert!(errors.iter().any(|e| matches!(
            e,
            SemanticError::MissingRequiredProperty { property, .. } if property == "timestep"
        )));
    }

    #[test]
    fn test_unknown_environment_property() {
        let program = Program {
            experiments: vec![sp(
                ExperimentDecl {
                    name: sp("Test".to_string(), 0, 4),
                    body: vec![sp(
                        Statement::EnvironmentBlock {
                            name: None,
                            properties: vec![sp(
                                Property {
                                    key: sp("gravity".to_string(), 10, 17),
                                    value: sp(Expr::NumberLiteral(9.8), 20, 23),
                                },
                                10,
                                23,
                            )],
                        },
                        5,
                        28,
                    )],
                },
                0,
                33,
            )],
        };

        let (_hir, errors) = lower(&program);
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            &errors[0],
            SemanticError::UnknownProperty { property, block, .. }
                if property == "gravity" && block == "environment"
        ));
    }

    #[test]
    fn test_hir_display() {
        let program = hello_bio_ast();
        let (hir, _) = lower(&program);
        let display = format!("{}", &hir.experiments[0]);
        assert!(display.contains("experiment HelloBiology"));
        assert!(display.contains("temperature = 310 K"));
        assert!(display.contains("timestep = 1 fs"));
    }
}
