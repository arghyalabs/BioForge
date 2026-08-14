//! # BioForge CLI
//!
//! Command-line interface for the BioForge compiler and simulation platform.
//!
//! Usage:
//! ```text
//! bio parse <file.bio>    — Parse and display AST
//! bio check <file.bio>    — Parse and report errors
//! bio run <file.bio>      — Run simulation (future)
//! ```

#![deny(unsafe_code)]

use clap::{Parser, Subcommand};
use std::process;

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
    /// Check a BioForge source file for errors without displaying the AST
    Check {
        /// Path to the .bio source file
        file: String,
    },
    /// Run a BioForge program (simulation runtime not yet implemented)
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
        eprintln!(
            "Successfully parsed {} experiment(s).",
            exp_count
        );
        0
    }
}

/// Check command: parse the file and report errors without showing AST.
fn cmd_check(path: &str) -> i32 {
    let source = match read_source(path) {
        Ok(s) => s,
        Err(code) => return code,
    };

    let (program, diagnostics) = bioforge_parser::parse(&source);

    if !diagnostics.is_empty() {
        bioforge_diagnostics::render_diagnostics(path, &source, &diagnostics);
    }

    let error_count = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Error)
        .count();
    let warning_count = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Warning)
        .count();
    let exp_count = program.experiments.len();

    if error_count > 0 {
        eprintln!(
            "Check failed: {} error(s), {} warning(s) in {} experiment(s).",
            error_count, warning_count, exp_count
        );
        1
    } else if warning_count > 0 {
        eprintln!(
            "Check passed with {} warning(s) in {} experiment(s).",
            warning_count, exp_count
        );
        0
    } else {
        eprintln!(
            "Check passed: {} experiment(s), no errors.",
            exp_count
        );
        0
    }
}

/// Run command: parse first, then explain that simulation is not yet available.
fn cmd_run(path: &str) -> i32 {
    let source = match read_source(path) {
        Ok(s) => s,
        Err(code) => return code,
    };

    let (program, diagnostics) = bioforge_parser::parse(&source);

    // Show any parse errors
    if !diagnostics.is_empty() {
        bioforge_diagnostics::render_diagnostics(path, &source, &diagnostics);
    }

    let error_count = diagnostics
        .iter()
        .filter(|d| d.severity == bioforge_diagnostics::DiagnosticSeverity::Error)
        .count();

    if error_count > 0 {
        eprintln!(
            "Cannot run: {} parse error(s) found. Fix syntax errors first.",
            error_count
        );
        return 1;
    }

    let exp_count = program.experiments.len();
    eprintln!("Parsed {} experiment(s) successfully.", exp_count);
    eprintln!();
    eprintln!("╔══════════════════════════════════════════════════════════════╗");
    eprintln!("║  Simulation runtime is not yet implemented.                 ║");
    eprintln!("║                                                             ║");
    eprintln!("║  The BioForge simulation engine (Phase 5+) will provide:    ║");
    eprintln!("║    • SimulationState construction                           ║");
    eprintln!("║    • Physics engine (Velocity Verlet integration)           ║");
    eprintln!("║    • Molecular structure loading (PDB, SDF, SMILES)         ║");
    eprintln!("║    • Trajectory generation                                  ║");
    eprintln!("║    • Measurement recording                                  ║");
    eprintln!("║    • 3D/4D visualization                                    ║");
    eprintln!("║                                                             ║");
    eprintln!("║  Current status: Phase 1 (Language Core) ✓                  ║");
    eprintln!("║  Next milestone: Phase 2 (Unit & Type System)               ║");
    eprintln!("╚══════════════════════════════════════════════════════════════╝");
    0
}
