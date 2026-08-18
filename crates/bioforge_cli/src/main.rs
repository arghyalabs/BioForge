//! # BioForge CLI
//!
//! Command-line interface for the BioForge compiler and simulation platform.
//!
//! Usage:
//! ```text
//! bio parse <file.bio>    — Parse and display AST
//! bio check <file.bio>    — Parse, validate semantics, and check units
//! bio run <file.bio>      — Run end-to-end biological simulation pipeline
//! ```

#![deny(unsafe_code)]

use clap::{Parser, Subcommand};
use std::path::Path;
use std::process;

use bioforge_biology::pdb::parse_pdb_file;
use bioforge_biology::{Atom, Bond, Element, Molecule};
use bioforge_hir::{HirEntity, HirExperiment, HirExpr};
use bioforge_measurement::{
    DistanceObservable, KineticEnergyObservable, MeasurementEngine, PotentialEnergyObservable,
    RadiusOfGyrationObservable, TemperatureObservable, TotalEnergyObservable,
};
use bioforge_physics::{
    BerendsenThermostat, CompositeForceField, ForceField, Integrator, Thermostat, VelocityVerlet,
};
use bioforge_render::{RenderStyle, Scene};
use bioforge_state::{SimulationState, Trajectory};

#[derive(Parser)]
#[command(name = "bio")]
#[command(version)]
#[command(about = "BioForge — Biology-native programming language and simulation platform")]
#[command(
    long_about = "BioForge is a biology-native programming language and multiscale biological \
    simulation platform.\n\n\
    The BioForge language allows scientists to describe biological systems and mechanisms \
    in a computationally executable way.\n\n\
    Mission: Make biology computationally programmable."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Parse a BioForge source file and display the AST
    Parse {
        /// Path to the .bio source file
        file: String,
    },
    /// Check a BioForge source file for syntax, semantic, and unit errors
    Check {
        /// Path to the .bio source file
        file: String,
    },
    /// Run a BioForge simulation experiment end-to-end
    Run {
        /// Path to the .bio source file
        file: String,
    },
}

fn main() {
    let cli = Cli::parse();

    let exit_code = match cli.command {
        Commands::Parse { file } => cmd_parse(&file),
        Commands::Check { file } => cmd_check(&file),
        Commands::Run { file } => cmd_run(&file),
    };

    process::exit(exit_code);
}

/// Read a source file, returning its contents or printing an error.
fn read_source(path: &str) -> Result<String, i32> {
    match std::fs::read_to_string(path) {
        Ok(source) => Ok(source),
        Err(e) => {
            eprintln!("Error: Could not read file '{}': {}", path, e);
            Err(1)
        }
    }
}

/// Parse command: parse the file and display the AST.
fn cmd_parse(path: &str) -> i32 {
    let source = match read_source(path) {
        Ok(s) => s,
        Err(code) => return code,
    };

    let (program, diagnostics) = bioforge_parser::parse(&source);

    // Render any diagnostics
    if !diagnostics.is_empty() {
        bioforge_diagnostics::render_diagnostics(path, &source, &diagnostics);
    }

    let error_count = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Error)
        .count();

    // Print the AST
    println!("{}", program.pretty_print());

    // Print summary
    let exp_count = program.experiments.len();
    if error_count > 0 {
        eprintln!(
            "Parsed with {} error(s), {} experiment(s) found.",
            error_count, exp_count
        );
        1
    } else {
        eprintln!("Successfully parsed {} experiment(s).", exp_count);
        0
    }
}

/// Check command: parse the file, run semantic analysis, and report errors.
fn cmd_check(path: &str) -> i32 {
    let source = match read_source(path) {
        Ok(s) => s,
        Err(code) => return code,
    };

    let (program, diagnostics) = bioforge_parser::parse(&source);

    if !diagnostics.is_empty() {
        bioforge_diagnostics::render_diagnostics(path, &source, &diagnostics);
    }

    let parse_errors = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Error)
        .count();

    if parse_errors > 0 {
        eprintln!(
            "Check failed: {} parse error(s). Fix syntax errors first.",
            parse_errors
        );
        return 1;
    }

    // Semantic analysis — lower AST to HIR
    let (_hir, semantic_errors) = bioforge_hir::lower(&program);

    if !semantic_errors.is_empty() {
        let sem_diagnostics: Vec<bioforge_diagnostics::Diagnostic> = semantic_errors
            .iter()
            .map(|e| {
                bioforge_diagnostics::Diagnostic::error(e.to_string())
                    .with_label(e.span(), e.to_string())
            })
            .collect();

        bioforge_diagnostics::render_diagnostics(path, &source, &sem_diagnostics);

        eprintln!(
            "Check failed: {} semantic error(s) in {} experiment(s).",
            semantic_errors.len(),
            program.experiments.len()
        );
        return 1;
    }

    let exp_count = program.experiments.len();
    let warning_count = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Warning)
        .count();

    if warning_count > 0 {
        eprintln!(
            "Check passed with {} warning(s) in {} experiment(s).",
            warning_count, exp_count
        );
    } else {
        eprintln!(
            "Check passed: {} experiment(s), no errors. (syntax ✓, semantics ✓, units ✓)",
            exp_count
        );
    }
    0
}

/// Run command: parse, validate semantics, and execute the full simulation pipeline.
fn cmd_run(path: &str) -> i32 {
    let source = match read_source(path) {
        Ok(s) => s,
        Err(code) => return code,
    };

    let (program, diagnostics) = bioforge_parser::parse(&source);

    if !diagnostics.is_empty() {
        bioforge_diagnostics::render_diagnostics(path, &source, &diagnostics);
    }

    let parse_errors = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Error)
        .count();

    if parse_errors > 0 {
        eprintln!(
            "Cannot run: {} parse error(s). Fix syntax errors first.",
            parse_errors
        );
        return 1;
    }

    // Lower AST to HIR with full semantic & unit validation
    let (hir, semantic_errors) = bioforge_hir::lower(&program);

    if !semantic_errors.is_empty() {
        let sem_diagnostics: Vec<bioforge_diagnostics::Diagnostic> = semantic_errors
            .iter()
            .map(|e| {
                bioforge_diagnostics::Diagnostic::error(e.to_string())
                    .with_label(e.span(), e.to_string())
            })
            .collect();

        bioforge_diagnostics::render_diagnostics(path, &source, &sem_diagnostics);

        eprintln!(
            "Cannot run: {} semantic error(s). Fix validation errors first.",
            semantic_errors.len()
        );
        return 1;
    }

    let base_dir = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
    let output_dir = Path::new("output");
    let _ = std::fs::create_dir_all(output_dir);

    for exp in &hir.experiments {
        run_experiment(exp, base_dir, output_dir);
    }

    0
}

/// Execute a single compiled biological experiment end-to-end.
fn run_experiment(exp: &HirExperiment, base_dir: &Path, output_dir: &Path) {
    println!("\n╔══════════════════════════════════════════════════════════════╗");
    println!("║                BioForge Simulation Engine                    ║");
    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║ Experiment:   {:<47}║", exp.name);

    // 1. Biological Entity Ingestion
    let mut molecules = Vec::new();
    for entity in &exp.entities {
        let mol = resolve_entity_to_molecule(entity, base_dir);
        println!(
            "║ Entity:       {:<12} ({} atoms, {:.1} Da){:>14}║",
            entity.name,
            mol.atom_count(),
            mol.total_mass(),
            ""
        );
        molecules.push(mol);
    }

    // 2. Simulation State & Thermalization
    let mut state = SimulationState::from_molecules(&molecules, Some([50.0, 50.0, 50.0]));
    let temp_k = exp
        .environment
        .as_ref()
        .and_then(|e| e.temperature.as_ref())
        .map(|q| q.to_si_value())
        .unwrap_or(300.0);

    if let Err(e) = state.thermalize(temp_k, 42) {
        eprintln!("║ Error during thermalization: {:<32}║", e);
        return;
    }

    // 3. Physical Parameters Resolution
    let dt_ps = exp
        .simulation
        .as_ref()
        .map(|s| s.timestep.to_si_value() * 1e12)
        .unwrap_or(0.0005)
        .max(1e-5);

    let duration_ps = exp
        .simulation
        .as_ref()
        .map(|s| s.duration.to_si_value() * 1e12)
        .unwrap_or(1.0)
        .max(dt_ps);

    let total_steps = ((duration_ps / dt_ps).round() as u64).max(1);
    let stride = (total_steps / 100).max(1);

    println!(
        "║ Environment:  T = {:.1} K, pH = {:.1}{:>30}║",
        temp_k,
        exp.environment.as_ref().and_then(|e| e.ph).unwrap_or(7.0),
        ""
    );
    println!(
        "║ Simulation:   {} steps (dt = {:.2} fs, duration = {:.2} ps){:>6}║",
        total_steps,
        dt_ps * 1000.0,
        duration_ps,
        ""
    );
    println!("╠══════════════════════════════════════════════════════════════╣");

    // 4. Force Field, Numerical Integrator & Thermostat
    let force_field = CompositeForceField::standard_molecular_mechanics(1.0, 10.0);
    let mut integrator = VelocityVerlet::new();
    let mut thermostat = BerendsenThermostat::new(temp_k, 0.05); // tau = 50 fs

    // 5. Observables & Measurement Engine
    let mut measurement_engine = MeasurementEngine::new();
    let mut trajectory = Trajectory::new(stride);

    measurement_engine.add_observable(TotalEnergyObservable);
    measurement_engine.add_observable(KineticEnergyObservable);
    measurement_engine.add_observable(PotentialEnergyObservable);
    measurement_engine.add_observable(TemperatureObservable);
    measurement_engine.add_observable(RadiusOfGyrationObservable::new("radius_of_gyration"));

    if state.num_atoms >= 2 {
        measurement_engine.add_observable(DistanceObservable::pair(
            "end_to_end_dist",
            0,
            state.num_atoms - 1,
        ));
    }

    // Record t=0 initial state
    let u0 = force_field.compute_forces(&mut state).unwrap_or(0.0);
    trajectory.record_frame(&state, u0);
    let _ = measurement_engine.record_step(&state, u0);

    // 6. Numerical Integration Loop (with NVT Berendsen temperature regulation)
    for step in 1..=total_steps {
        let u = match integrator.step(&mut state, dt_ps, &force_field) {
            Ok(val) => val,
            Err(e) => {
                eprintln!("║ Simulation failed at step {}: {:<30}║", step, e);
                break;
            }
        };

        // Thermostat velocity rescaling
        let _ = thermostat.apply(&mut state, dt_ps);

        if step % stride == 0 || step == total_steps {
            trajectory.record_frame(&state, u);
            let _ = measurement_engine.record_step(&state, u);
        }
    }

    // 7. 3D Scene Generation
    let scene = Scene::from_state(&state, RenderStyle::ball_and_stick());

    // 8. Export Artifacts
    let exp_name = &exp.name;
    let traj_path = output_dir.join(format!("{}_trajectory.xyz", exp_name));
    let scene_path = output_dir.join(format!("{}_scene.obj", exp_name));
    let csv_path = output_dir.join(format!("{}_measurements.csv", exp_name));
    let json_path = output_dir.join(format!("{}_measurements.json", exp_name));

    let _ = std::fs::write(&traj_path, trajectory.to_xyz(&state));
    let _ = std::fs::write(&scene_path, scene.export_obj());
    let _ = std::fs::write(&csv_path, measurement_engine.export_csv());
    if let Ok(json) = measurement_engine.export_json() {
        let _ = std::fs::write(&json_path, json);
    }

    // 9. Display Rich Summary Dashboard
    println!("║ Observables Summary:                                         ║");
    if let Some(stats) = measurement_engine.statistics("total_energy") {
        println!(
            "║   • Total Energy:       {:.3} ± {:.3} kJ/mol{:>17}║",
            stats.mean, stats.std_dev, ""
        );
    }
    if let Some(stats) = measurement_engine.statistics("temperature") {
        println!(
            "║   • Temperature:        {:.1} ± {:.1} K{:>25}║",
            stats.mean, stats.std_dev, ""
        );
    }
    if let Some(stats) = measurement_engine.statistics("radius_of_gyration") {
        println!(
            "║   • Radius of Gyration: {:.3} ± {:.3} Å{:>22}║",
            stats.mean, stats.std_dev, ""
        );
    }

    println!("╠══════════════════════════════════════════════════════════════╣");
    println!("║ Generated Artifacts:                                         ║");
    println!("║   [XYZ]  {:<51}║", traj_path.display());
    println!("║   [OBJ]  {:<51}║", scene_path.display());
    println!("║   [CSV]  {:<51}║", csv_path.display());
    println!("║   [JSON] {:<51}║", json_path.display());
    println!("╚══════════════════════════════════════════════════════════════╝\n");
}

/// Resolve an entity initializer into a physical biological molecule.
fn resolve_entity_to_molecule(entity: &HirEntity, base_dir: &Path) -> Molecule {
    match &entity.initializer {
        HirExpr::FunctionCall { name, args }
            if name == "load_structure" || name == "load" || name == "pdb" =>
        {
            if let Some(HirExpr::String(file_path)) = args.first() {
                // Try candidate paths: relative to script, or relative to current working directory
                let candidate1 = base_dir.join(file_path);
                let candidate2 = Path::new(file_path).to_path_buf();

                let resolved_path = if candidate1.exists() {
                    Some(candidate1)
                } else if candidate2.exists() {
                    Some(candidate2)
                } else {
                    None
                };

                if let Some(path) = resolved_path {
                    if let Ok(mol) = parse_pdb_file(path.to_str().unwrap_or(file_path)) {
                        return mol;
                    }
                }
            }
        }
        _ => {}
    }

    // Default fallback: 3-atom water molecule
    let o = Element::from_symbol("O").unwrap();
    let h = Element::from_symbol("H").unwrap();
    let mut mol = Molecule::new(&entity.name);
    mol.atoms.push(Atom::new(1, o, [0.0, 0.0, 0.0], "O"));
    mol.atoms.push(Atom::new(2, h, [0.957, 0.0, 0.0], "H1"));
    mol.atoms.push(Atom::new(3, h, [-0.240, 0.927, 0.0], "H2"));
    mol.bonds.push(Bond::single(1, 2));
    mol.bonds.push(Bond::single(1, 3));
    mol
}
