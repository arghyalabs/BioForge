# BioForge: ROADMAP.md

This document outlines the comprehensive development roadmap for the BioForge project.

## 1. Current Status
**Phase 0 + Phase 1 (Architecture + Language Core) — Active**
We are currently establishing the foundational architecture, repository structure, and the core compiler components for BioLang.

## 2. Phase Table

### Phase 0 — Architecture & Documentation (CURRENT)
- Repository structure
- Documentation (`AGENTS.md`, `ROADMAP.md`)
- Scientific principles definition

### Phase 1 — Language Core (CURRENT)
- Lexer, Parser, AST, Diagnostics, CLI
- First command: `bio parse examples/hello.bio`

### Phase 2 — Unit & Type System
- Quantity `{ value, unit, dimension }`
- Support for distance, time, mass, temperature, energy, concentration, voltage, pressure, force
- Invalid dimensional operations fail at compile/runtime

### Phase 3 — BioIR (Intermediate Representation)
- `ProteinDeclaration`, `LigandDeclaration`, `EnvironmentDeclaration`, `SimulationDeclaration`
- AST → BioIR conversion
- Independent from source syntax

### Phase 4 — Structure Loading
- PDB, mmCIF, SDF, SMILES format support
- Using `pdbtbx` and `chematic-smiles` crates
- Internal types: `Structure`, `Atom`, `Bond`, `Residue`, `Chain`, `Molecule`, `Protein`, `Ligand`

### Phase 5 — SimulationState
- `SimulationState { time, entities, positions, velocities, forces, energies }`
- Serialization
- Independent from renderer

### Phase 6 — Basic Physics
- `F = ma`, position, velocity, mass, force, acceleration
- Velocity Verlet integrator
- Validation against harmonic oscillator, conservation laws

### Phase 7 — Molecular Physics
- Bond, angle, dihedral, van der Waals, electrostatic forces
- Validation benchmarks for each force type

### Phase 8 — Measurements
- `distance()`, `velocity()`, `force()`, `energy()`, `temperature()`
- Time series generation

### Phase 9 — 3D Renderer
- `wgpu`-based rendering
- Representations: Ball-and-stick, space-filling, cartoon, surface

### Phase 10 — 4D Timeline
- Play, pause, step, rewind, seek, time slider

### Phase 11 — Scientific Dashboard
- Synchronized 3D world, timeline, measurements, plots, logs

### Phase 12 — First Complete Vertical Slice
- Full pipeline: `example.bio` → compiler → AST → BioIR → Simulation → Trajectory → Measurements → 3D/4D
- Target command: `bio run examples/protein_ligand.bio`

### Phase 13 — Chemical Systems
- Reactions, concentrations, rates, equilibrium, diffusion

### Phase 14 — Membrane Biology
- Membrane, lipid, ion channel, ion transport, membrane potential

### Phase 15 — Neuroscience
- Neuron, synapse, ion channels, action potential, Hodgkin-Huxley

### Phase 16 — Cellular Systems
- Cell, organelle, transport, signaling, gene expression

### Phase 17 — Tissue & Physiology
- Cell populations, tissues, fluid dynamics, mechanical properties

### Phase 18 — Multiscale Coupling
- Scale bridges between molecular/cellular/tissue/organ/physiology

### Phase 19 — Hypothesis System
- `hypothesis`, `mechanism`, `prediction` as language primitives
- Causal chain representation

### Phase 20 — Counterfactual Simulation
- Experiment branching, perturbation comparison

### Phase 21 — Experiment Branching & Provenance
- Full experiment state capture, parameter provenance

### Phase 22 — Interactive Simulation
- Real-time pause/resume/perturb with model-consistent behavior

### Phase 23-27 — Biological Digital Twins & Advanced
- High-fidelity computational representations
- Fidelity dimensions and reporting

### Phase 28 — Full Programmable Biological World
- The complete vision: from atomic scale up to the organism scale, all computationally programmable and simulatable.

## 3. Success Criteria
- **Architecture Validation:** Strict adherence to module boundaries and scientific principles.
- **Test Coverage:** High coverage across units, integration, snapshot, and scientific validation tests.
- **Phase Completion:** Each phase is considered complete when its features are implemented, tested, validated, and integrated into the overarching pipeline.

## 4. First Vertical Slice Target
The most critical near-term target is completing Phase 12. This involves the end-to-end flow:
1. Parse `example.bio`
2. Generate AST
3. Convert to BioIR
4. Load structural data and initialize `SimulationState`
5. Run numerical simulation to generate trajectory
6. Record measurements
7. Visualize in the 3D/4D Renderer & Scientific Dashboard

## 5. Development Principles
- Architecture first.
- Scientific correctness is paramount.
- Validation must occur throughout every phase.
- Implementation must be iterative and incremental.
