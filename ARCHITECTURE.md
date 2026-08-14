# BioForge Technical Architecture

## 1. Overview

BioForge is uniquely positioned as both a **programming language (BioForge)** and a **scientific computing platform**. The architecture must support the parsing, semantic analysis, and execution of biological domain logic, while simultaneously providing a robust, highly performant engine for numerical simulation and multiscale physics.

This document details the architectural principles and design decisions driving the BioForge project.

## 2. Compiler Pipeline

BioForge utilizes a modern, multi-stage compiler pipeline designed to transform domain-specific biological descriptions into executable simulation plans.

```
Source Code → Lexer → Tokens → Parser → AST → Semantic Analysis → Type/Unit Checking → BioIR → Execution Plan → Runtime
```

1. **Source Code**: Raw `.bio` files written by the user.
2. **Lexer**: Breaks the raw text into a stream of meaningful tokens (keywords, identifiers, literals, operators).
3. **Parser**: Analyzes the token stream to produce an Abstract Syntax Tree (AST) representing the syntactic structure of the code.
4. **Semantic Analysis**: Resolves symbols, scopes, and ensures the biological meaning is coherent.
5. **Type/Unit Checking**: Validates data types and performs rigorous physical unit dimensional analysis (e.g., ensuring you cannot add Joules to Kelvin).
6. **BioIR (Biological Intermediate Representation)**: A lower-level, syntax-independent representation of the biological system and simulation plan.
7. **Execution Plan**: Translates BioIR into instructions for specific mathematical solvers and physics engines.
8. **Runtime**: The highly optimized engine that actually steps the simulation forward in time.

*(Currently implemented: Lexer → Parser → AST).*

## 3. Crate Architecture

The Rust workspace is divided into specialized, decoupled crates to ensure maintainability and fast compilation times.

```
bioforge_cli (binary frontend)
    │
    ├── bioforge_parser
    │   │
    │   ├── bioforge_lexer
    │   │   └── bioforge_diagnostics
    │   │
    │   ├── bioforge_ast
    │   │   └── bioforge_diagnostics
    │   │
    │   └── bioforge_diagnostics
    │
    └── bioforge_diagnostics
```

* `bioforge_diagnostics` is a core utility used across the pipeline for mapping errors back to original source code spans.

## 4. SimulationState Independence

**This is the most critical architectural principle of BioForge.**

```
Simulation
     ↓
SimulationState (The single source of scientific truth)
     ├──────────→ Analysis (Extracts data)
     ├──────────→ Measurements (Computes observables)
     └──────────→ Renderer (Visualizes state)
```

* The **SimulationState** is the absolute source of truth. It contains atomic coordinates, velocities, concentrations, topology, and all physical parameters.
* The **Renderer NEVER owns scientific state**. It acts purely as a read-only observer. It interpolates, shades, and visualizes the state for human consumption, but it cannot alter the physics.
* **Scientific correctness > Visual beauty**. If a rendering optimization compromises the accuracy of the displayed scientific state, it is rejected.

## 5. Solver Plurality

Biology spans dozens of orders of magnitude in space and time. No single mathematical framework can efficiently simulate everything from electron orbitals to whole-organism physiology. BioForge embraces **Solver Plurality**.

Different scales require specialized mathematical models:
* **Molecular dynamics (MD)**: For atomic/molecular resolution (Newtonian mechanics, force fields).
* **Brownian dynamics**: For mesoscale systems where solvent is implicit.
* **ODE systems (Ordinary Differential Equations)**: For reaction kinetics, metabolic pathways, and uniform concentrations.
* **PDE solvers (Partial Differential Equations)**: For diffusion, gradients, and spatial transport.
* **Stochastic simulations (SSA/Gillespie)**: For single-molecule regimes with high noise.
* **Agent-based models**: For cell populations and cellular behavior.
* **FEM (Finite Element Method)**: For tissue mechanics and macroscopic physics.
* **Electrical models (e.g., Hodgkin-Huxley)**: For neuronal action potentials.

The architecture provides a unified biological abstraction (BioIR) that orchestrates these underlying, specialized solvers.

## 6. Multiscale Architecture and Scale Bridges

To achieve true multiscale simulation, BioForge utilizes "Scale Bridges" that map information between different resolution domains.

```
Molecular Model (High spatial/temporal resolution, small domain)
      ↓
Scale Bridge (Aggregates parameters, passes signals, computes macroscopic concentrations/forces)
      ↓
Cellular Model (Medium resolution, larger domain)
      ↓
Scale Bridge (Maps cellular outputs to tissue properties)
      ↓
Tissue Model (Continuum mechanics)
      ↓
Scale Bridge
      ↓
Physiological Model (Systemic flows, ODE networks)
```
These bridges handle the rigorous mathematical coupling (e.g., passing average forces up, passing boundary conditions down).

## 7. Interactive Simulation

BioForge is designed to be interactive, not just a static batch-processing job.
* **Pause/Resume/Step**: The runtime supports pausing the simulation, advancing by discrete steps, and resuming.
* **Perturbations**: Users can inject perturbations (e.g., adding a drug, changing temperature, mutating a residue) during runtime.
* **Strict Interpretation**: All physical perturbations are interpreted *through the scientific model*. You cannot just "drag" an atom in the viewer and break physics; you must apply a virtual force, and the system responds according to the equations of motion.

## 8. Experiment Branching

Because biological systems are complex and chaotic, BioForge supports forkable simulation states. You can run a simulation up to a critical point (e.g., before cell division), fork the state into multiple branches, apply different conditions to each branch (counterfactual comparisons), and run them in parallel to observe diverging outcomes.

## 9. Fidelity Dimensions

A BioForge simulation is parameterized along several axes of fidelity:
* **Structural**: Coarse-grained vs. All-atom vs. Quantum.
* **Physical**: Implicit vs. Explicit solvent, fixed vs. flexible bonds.
* **Chemical**: Fixed protonation vs. dynamic pH, reactive vs. non-reactive force fields.
* **Biological**: Number of interacting pathways considered.
* **Spatial**: Point-particle vs. spatially resolved.
* **Temporal**: Timestep resolution.

## 10. Provenance Tracking

Every parameter, force field, structural starting point, and assumption must be traceable. BioForge tracks the *provenance* of its data. If a specific binding affinity is used in an ODE, the system should know which paper or database that number came from.

## 11. Reproducibility

Reproducibility is paramount. BioForge captures the full experiment state—the exact version of the compiler, the seed for random number generators, the starting coordinates, and the exact solver configurations. Running the same `.bio` file with the same environment must yield the exact same physical trajectory.

## 12. BioStudio Vision

While the core is a CLI and runtime, the ultimate vision is **BioStudio**: a tightly integrated IDE. It will combine a code editor for BioForge, a 3D/4D viewport for visualizing the multiscale world, an interactive timeline (scrubbing through simulation history), live measurement graphs, and scientific logs, all synchronized in real-time.

## 13. Future: BioIR

The Biological Intermediate Representation (BioIR) will be a critical layer. It decouples the BioForge syntax from the backend physics engines. BioIR represents a graph of biological entities, their states, and their interactions, allowing optimizations and transformations before the execution plan is handed to the specific numerical solvers.

## 14. Future: Unit System

In biology and physics, numbers without units are meaningless and dangerous. BioForge will implement units as first-class citizens in the type system. Dimensional analysis will occur at compile-time (or semantic analysis time). `310 K` is fundamentally a different type than `10 ps`. The compiler will statically prevent incorrect dimensional operations (e.g., adding a distance to an energy) and automatically handle necessary conversions (e.g., `1 ns + 500 ps = 1.5 ns`).
