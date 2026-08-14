# BioForge: AGENTS.md

Welcome to the BioForge project. This document is **mandatory reading** for any AI agent working on this codebase.

## 1. Project Purpose
BioForge is a biology-native programming language (BioForge) and multiscale biological simulation platform. The deeper mission of BioForge is to make biology computationally programmable, bridging the gap between biological description and computational simulation.

## 2. Architecture Overview
BioForge consists of a compiler pipeline and a decoupled simulation/rendering engine. The general flow is:
**BioForge → Parser → BioIR → Biological Structure → SimulationState → Physics → Measurements → Trajectory → 3D/4D Renderer**

The project is structured into multiple crates (modules) ensuring strict boundaries and unidirectionality. The first vertical slice connects the source code all the way through to a 3D/4D rendered trajectory.

## 3. The 15 Non-Negotiable Scientific Principles
As an agent working on this codebase, you must adhere strictly to these principles:
1. Scientific correctness is more important than visual appearance.
2. Simulation state is the source of truth (the renderer merely observes it).
3. Units are first-class citizens.
4. Models must expose assumptions.
5. Approximations must be visible.
6. Validation must be explicit.
7. Uncertainty must be represented where possible.
8. Experimental data must remain distinguishable from simulation results.
9. Predictions must not be presented as established biological facts.
10. Reproducibility must be supported.
11. Scientific provenance must be preserved.
12. Different biological scales may require different mathematical models.
13. Do not pretend one solver can accurately simulate all of biology.
14. Do not claim exact replication of biological reality.
15. Aim for progressively increasing fidelity.

## 4. Coding Conventions
- **Language:** Rust 2021 edition.
- **Safety:** No unsafe code (`#![deny(unsafe_code)]`).
- **Documentation:** All public APIs must have doc comments (`///`).
- **Traits:** Use `#[derive(Debug, Clone, PartialEq)]` on data types whenever possible.
- **Errors:** Use the `thiserror` crate for error types.
- **Diagnostics:** Use `ariadne` for diagnostic rendering (compiler errors/warnings).
- **Linting:** Follow all `clippy` warnings.
- **Serialization:** Use `serde` for serialization/deserialization where needed.

## 5. Module Boundaries & Dependency Rules
Crate dependencies must be strictly unidirectional. **NEVER create circular dependencies.**
- `diagnostics`: leaf crate, no internal dependencies.
- `ast`: depends only on `diagnostics`.
- `lexer`: depends only on `diagnostics`.
- `parser`: depends on `lexer`, `ast`, `diagnostics`.
- `cli`: depends on `parser`, `ast`, `diagnostics`.
- Future crates will follow this same unidirectional dependency rule.

## 6. Testing Requirements
- Every module must have unit tests.
- Use `insta` for snapshot testing of AST output and error messages.
- Place integration tests in the `tests/` directory.
- Every scientific model (future) must have explicit validation tests.

## 7. Units Policy
All physical quantities MUST carry units. **Naked floating-point numbers for physical values are strictly forbidden** in the runtime. Units must be explicitly defined and verified.

## 8. Numerical Precision Policy
Numerical methods, stability characteristics, and known limitations must be documented. We prioritize robust and well-understood numerical solvers.

## 9. Forbidden Architectural Changes (CRITICAL)
- `SimulationState` must NEVER be owned by the renderer.
- The renderer must NEVER independently determine biological behavior.
- No bypassing the scientific model pipeline (Scientific Model → Mathematical Model → Numerical Solver → SimulationState → Renderer).
- No unvalidated scientific claims in code or documentation.
- No arbitrary "magic" manipulation in simulations — all perturbations go through the model.
- Do not merge scales without explicit scale bridge interfaces.
- Do not force all biological systems into one universal solver.

## 10. AI Task Scoping Rules
- Every AI task must have a defined scope.
- Do not make uncontrolled architectural changes.
- Stay within the specified crate/module boundary.
- Do not modify code outside the task scope.
- Include tests for all new code.
- Do not add dependencies without justification.
- Read this `AGENTS.md` before modifying the repository.

## 11. API Stability
Public APIs require review before breaking changes. Design APIs with future extensibility and stability in mind.
