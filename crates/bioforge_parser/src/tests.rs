//! Tests for the BioLang parser.

use crate::parse;
use bioforge_ast::*;

#[test]
fn test_parse_empty_experiment() {
    let source = "experiment Empty {}";
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    assert_eq!(program.experiments.len(), 1);
    assert_eq!(program.experiments[0].node.name.node, "Empty");
    assert!(program.experiments[0].node.body.is_empty());
}

#[test]
fn test_parse_entity_decl() {
    let source = r#"experiment Test {
    protein receptor = load_structure("protein.pdb")
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    assert_eq!(program.experiments.len(), 1);
    let body = &program.experiments[0].node.body;
    assert_eq!(body.len(), 1);
    match &body[0].node {
        Statement::EntityDecl { kind, name, initializer } => {
            assert_eq!(kind.node, EntityKind::Protein);
            assert_eq!(name.node, "receptor");
            match &initializer.node {
                Expr::FunctionCall { name, args } => {
                    assert_eq!(name.node, "load_structure");
                    assert_eq!(args.len(), 1);
                }
                _ => panic!("Expected FunctionCall"),
            }
        }
        _ => panic!("Expected EntityDecl"),
    }
}

#[test]
fn test_parse_environment_block() {
    let source = r#"experiment Test {
    environment {
        temperature = 310 K
        pH = 7.4
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    let body = &program.experiments[0].node.body;
    assert_eq!(body.len(), 1);
    match &body[0].node {
        Statement::EnvironmentBlock { name, properties } => {
            assert!(name.is_none());
            assert_eq!(properties.len(), 2);
            assert_eq!(properties[0].node.key.node, "temperature");
            match &properties[0].node.value.node {
                Expr::Quantity { value, unit } => {
                    assert_eq!(*value, 310.0);
                    assert_eq!(unit.node, "K");
                }
                _ => panic!("Expected Quantity"),
            }
            assert_eq!(properties[1].node.key.node, "pH");
            match &properties[1].node.value.node {
                Expr::NumberLiteral(n) => assert_eq!(*n, 7.4),
                _ => panic!("Expected NumberLiteral"),
            }
        }
        _ => panic!("Expected EnvironmentBlock"),
    }
}

#[test]
fn test_parse_simulate_block() {
    let source = r#"experiment Test {
    simulate {
        timestep = 1 fs
        duration = 10 ps
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    let body = &program.experiments[0].node.body;
    assert_eq!(body.len(), 1);
    match &body[0].node {
        Statement::SimulateBlock { properties } => {
            assert_eq!(properties.len(), 2);
            assert_eq!(properties[0].node.key.node, "timestep");
            assert_eq!(properties[1].node.key.node, "duration");
        }
        _ => panic!("Expected SimulateBlock"),
    }
}

#[test]
fn test_parse_measure_block() {
    let source = r#"experiment Test {
    measure {
        distance(receptor, drug)
        energy(receptor, drug)
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    let body = &program.experiments[0].node.body;
    assert_eq!(body.len(), 1);
    match &body[0].node {
        Statement::MeasureBlock { measurements } => {
            assert_eq!(measurements.len(), 2);
            match &measurements[0].node {
                Expr::FunctionCall { name, args } => {
                    assert_eq!(name.node, "distance");
                    assert_eq!(args.len(), 2);
                }
                _ => panic!("Expected FunctionCall"),
            }
        }
        _ => panic!("Expected MeasureBlock"),
    }
}

#[test]
fn test_parse_visualize_block() {
    let source = r#"experiment Test {
    visualize {
        receptor
        drug
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    let body = &program.experiments[0].node.body;
    assert_eq!(body.len(), 1);
    match &body[0].node {
        Statement::VisualizeBlock { targets } => {
            assert_eq!(targets.len(), 2);
        }
        _ => panic!("Expected VisualizeBlock"),
    }
}

#[test]
fn test_parse_quantity() {
    let source = r#"experiment Test {
    environment {
        temperature = 310 K
        timestep = 1 fs
        distance = 10.5 nm
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    let body = &program.experiments[0].node.body;
    match &body[0].node {
        Statement::EnvironmentBlock { properties, .. } => {
            match &properties[0].node.value.node {
                Expr::Quantity { value, unit } => {
                    assert_eq!(*value, 310.0);
                    assert_eq!(unit.node, "K");
                }
                _ => panic!("Expected Quantity"),
            }
            match &properties[1].node.value.node {
                Expr::Quantity { value, unit } => {
                    assert_eq!(*value, 1.0);
                    assert_eq!(unit.node, "fs");
                }
                _ => panic!("Expected Quantity"),
            }
            match &properties[2].node.value.node {
                Expr::Quantity { value, unit } => {
                    assert_eq!(*value, 10.5);
                    assert_eq!(unit.node, "nm");
                }
                _ => panic!("Expected Quantity"),
            }
        }
        _ => panic!("Expected EnvironmentBlock"),
    }
}

#[test]
fn test_parse_full_hello_bio() {
    let source = r#"experiment HelloBiology {
    protein receptor = load_structure("protein.pdb")
    ligand drug = load("drug.sdf")

    environment {
        temperature = 310 K
        pH = 7.4
    }

    simulate {
        timestep = 1 fs
        duration = 10 ps
    }

    measure {
        distance(receptor, drug)
        energy(receptor, drug)
    }

    visualize {
        receptor
        drug
    }
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    assert_eq!(program.experiments.len(), 1);
    let exp = &program.experiments[0].node;
    assert_eq!(exp.name.node, "HelloBiology");
    // Should have: protein, ligand, environment, simulate, measure, visualize = 6 statements
    assert_eq!(exp.body.len(), 6);
}

#[test]
fn test_error_recovery_missing_brace() {
    let source = r#"experiment Test {
    protein p = load("test.pdb")
"#;
    let (program, diagnostics) = parse(source);
    // Should still produce partial AST
    assert_eq!(program.experiments.len(), 1);
    // Should have at least one diagnostic about missing '}'
    assert!(!diagnostics.is_empty());
}

#[test]
fn test_snapshot_hello_bio_ast() {
    let source = r#"experiment HelloBiology {
    protein receptor = load_structure("protein.pdb")
    ligand drug = load("drug.sdf")
    environment {
        temperature = 310 K
        pH = 7.4
    }
    simulate {
        timestep = 1 fs
        duration = 10 ps
    }
    measure {
        distance(receptor, drug)
        energy(receptor, drug)
    }
    visualize {
        receptor
        drug
    }
}"#;
    let (program, _diagnostics) = parse(source);
    insta::assert_snapshot!(program.pretty_print());
}

#[test]
fn test_snapshot_simple_atom_ast() {
    let source = r#"experiment SimpleAtom {
    atom hydrogen = atom("H")
    environment {
        temperature = 300 K
    }
}"#;
    let (program, _diagnostics) = parse(source);
    insta::assert_snapshot!(program.pretty_print());
}

#[test]
fn test_multiple_experiments() {
    let source = r#"experiment A {
    protein p = load_structure("a.pdb")
}
experiment B {
    ligand l = load("b.sdf")
}"#;
    let (program, diagnostics) = parse(source);
    assert!(diagnostics.is_empty(), "Unexpected diagnostics: {:?}", diagnostics);
    assert_eq!(program.experiments.len(), 2);
    assert_eq!(program.experiments[0].node.name.node, "A");
    assert_eq!(program.experiments[1].node.name.node, "B");
}
