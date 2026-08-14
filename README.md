# BioForge

## Project Overview

**Project Name**: BioForge  
**Language Name**: BioForge  

## Mission Statement

"A biology-native programming language and computational scientific environment for constructing high-fidelity, interactive, measurable, multiscale simulations of biological systems — from molecular mechanisms to cellular processes, tissues, organs, and physiological systems — in order to computationally investigate mechanisms, test hypotheses, explore counterfactual biological scenarios, and generate scientifically interpretable predictions."

**Deeper Mission**: Make biology computationally programmable.

## What This Is NOT

* **Not just a molecular viewer**: While visualization is critical, BioForge is primarily a platform for *simulating* and *measuring* biological systems, not just looking at static structures.
* **Not a Python wrapper**: BioForge is a completely novel, natively compiled language specifically designed for the unique abstractions of biology, rather than an API bolted onto an existing general-purpose language.
* **Not a game engine**: Scientific correctness, physical validity, and accurate measurements are always prioritized over visual beauty or real-time frame rates.
* **Not a general-purpose language**: You won't use BioForge to write a web server. It is focused entirely on biological modeling.

## The Two-Level Vision

1. **Immediate Vision (The First Vertical Slice)**: A complete pipeline from parsing BioForge code to simulating a protein-ligand interaction in an environment, performing measurements, and rendering the trajectory in 3D/4D. (BioForge → Parser → BioIR → Biological Structure → SimulationState → Physics → Measurements → Trajectory → 3D/4D Renderer)
2. **Long-Term Vision**: A truly programmable biological world spanning from atomic scales to full organismic scales.

## Core Philosophy

The fundamental pipeline of BioForge follows a strict scientific methodology:

**Biological Description** → **Physical Model** → **Simulation** → **Measurements** → **Analysis** → **Prediction**

1. Specify the system in biological terms.
2. The compiler maps this to the appropriate physical/mathematical models.
3. The runtime executes the simulation.
4. Measurements are taken directly from the physical state.
5. The results are analyzed and visualized.

## Architecture Diagram

```
                    BIOLOGICAL LANGUAGE
                            │
                            ▼
                       Lexer / Parser
                            │
                            ▼
                           AST
                            │
                            ▼
                     Semantic Analysis
                            │
                            ▼
                         BioIR
                            │
                            ▼
                  Biological Runtime
                            │
              ┌─────────────┼─────────────┐
              │             │             │
          Molecular      Cellular     Physiological
           Models         Models         Models
              │             │             │
              └─────────────┼─────────────┘
                            ▼
                    Multiscale State
                            │
              ┌─────────────┴─────────────┐
              ▼                           ▼
      Scientific Analysis            Renderer
              │                           │
              │                       3D / 4D
              │                           │
              └─────────────┬─────────────┘
                            ▼
                      BioStudio IDE
```

## Example BioForge Code

BioForge provides domain-specific primitives for experiments, structures, environments, and measurements:

```bio
experiment HelloBiology {
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
}
```

## Quick Start

To build the compiler and parse a BioForge file:

```bash
# Clone the repository
git clone <repository_url>
cd BioForge

# Build the project
cargo build

# Run the parser on an example file
cargo run --bin bio -- parse examples/hello.bio
```
*(Note: Execution and simulation runtimes are currently in development).*

## Repository Structure

```
BioForge/
├── Cargo.toml            # Workspace definition
├── README.md             # This file
├── ARCHITECTURE.md       # Technical architecture documentation
├── crates/               # Core Rust crates
│   ├── bioforge_cli/     # Command-line interface (`bio`)
│   ├── bioforge_lexer/   # Lexical analysis (source text to tokens)
│   ├── bioforge_parser/  # Syntactic analysis (tokens to AST)
│   ├── bioforge_ast/     # Abstract Syntax Tree definitions
│   └── bioforge_diagnostics/ # Error reporting and source mapping
├── examples/             # Example BioForge scripts
└── docs/                 # Additional project documentation
```

## Development Status

**Active Phase:** Phase 0 + 1 (Architecture + Language Core)

| Phase | Description | Status |
|---|---|---|
| **Phase 0** | Architecture and scaffolding | Active |
| **Phase 1** | Language Core (Lexer, Parser, AST) | Active |
| **Phase 2** | Semantic Analysis & BioIR | Planned |
| **Phase 3** | Runtime & Basic Simulation | Planned |
| **Phase 4** | Multiscale Physics Integration | Planned |
| **Phase 5** | BioStudio IDE & Advanced Rendering | Planned |

## Technology Stack

BioForge is built with a focus on performance, safety, and excellent developer experience:
* **Rust**: The core language for the compiler, runtime, and simulation engines.
* **logos**: Fast, state-machine based lexical analyzer.
* **ariadne**: Beautiful, compiler-grade error diagnostics.
* **insta**: Snapshot testing for parsing and AST verification.
* **clap**: Command-line argument parsing.

## License

This project is dual-licensed under either the **MIT license** or the **Apache License, Version 2.0**, at your option.

## Contributing

We welcome contributions! For AI agents and automated contributors, please refer to the `AGENTS.md` file for strict rules and context regarding architecture preservation and code style. Human contributors can open issues and pull requests following standard open-source workflows.
