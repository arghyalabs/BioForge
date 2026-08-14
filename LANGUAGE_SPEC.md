# BioForge Language Specification (v0.1)

## 1. Language Overview

**BioForge** is a domain-specific programming language engineered specifically for describing biological systems, physical interactions, and scientific experiments. It serves as the primary interface for the BioForge platform.

## 2. Design Philosophy

Scientists should express biological systems, not generic computations. 

Traditional languages force scientists to translate their domain knowledge into arrays, loops, and pointers. BioForge flips this paradigm: the language primitives are the biological and physical concepts themselves. You declare what the system *is* and how it *behaves*, while the compiler maps this to the optimal mathematical solvers and data structures.

---

## 3. v0.1 Grammar (EBNF-Style)

```ebnf
Program         ::= ExperimentDecl*

ExperimentDecl  ::= 'experiment' IDENT '{' Statement* '}'

Statement       ::= EntityDecl 
                  | EnvironmentBlock 
                  | SimulateBlock 
                  | MeasureBlock 
                  | VisualizeBlock 
                  | Assignment

EntityDecl      ::= EntityKind IDENT '=' Expr

EntityKind      ::= 'protein' | 'ligand' | 'ion' | 'molecule' | 'atom'

EnvironmentBlock ::= 'environment' (IDENT)? '{' Property* '}'

SimulateBlock   ::= 'simulate' '{' Property* '}'

MeasureBlock    ::= 'measure' '{' MeasureExpr* '}'

VisualizeBlock  ::= 'visualize' '{' Expr* '}'

Property        ::= IDENT '=' Expr

MeasureExpr     ::= FunctionCall

Expr            ::= Primary (BinOp Primary)*

Primary         ::= Quantity | Literal | FunctionCall | Identifier | '(' Expr ')'

Quantity        ::= NUMBER UNIT

FunctionCall    ::= IDENT '(' (Expr (',' Expr)*)? ')'

BinOp           ::= '+' | '-' | '*' | '/'

Literal         ::= NUMBER | STRING | BOOL
```

---

## 4. Token Definitions

* **Keywords**: `experiment`, `protein`, `ligand`, `ion`, `molecule`, `atom`, `environment`, `simulate`, `measure`, `visualize`, `true`, `false`
* **Units**: `K`, `fs`, `ps`, `ns`, `nm`, `Å`, `mM`, `mV`, `atm`, `kJ/mol`, `kcal/mol`
* **Symbols**: `{ } ( ) [ ] = , . + - * /`
* **Numbers**: Standard integers (e.g., `310`) and floating-point values (e.g., `7.4`, `1.5e-3`)
* **Strings**: Double-quoted text (e.g., `"protein.pdb"`)
* **Identifiers**: Alphanumeric strings starting with a letter or underscore (e.g., `[a-zA-Z_][a-zA-Z0-9_]*`)
* **Comments**: `//` for single-line comments, `/* ... */` for block comments
* **Semicolons**: Not required to terminate statements

---

## 5. Semantic Blocks

### `experiment`
The root-level container for a scientific simulation. It groups the biological entities, the environmental conditions, the simulation parameters, and the observation directives.

### `environment`
Defines the boundary conditions and thermodynamic state of the system, such as temperature, pressure, solvent type, and ionic strength.

### `simulate`
Instructs the BioForge engine to advance the system through time. Properties inside this block dictate the temporal extent (e.g., `time`), integration step, and solver hints.

### `measure`
Defines the observables to extract from the `SimulationState`. Measurements do not affect the simulation; they only extract data (e.g., RMSD, energy, distances) into a trajectory or data file.

### `visualize`
Instructs the rendering engine on how to display the `SimulationState`. Like `measure`, visualization is strictly an observer.

---

## 6. Example Programs

### `hello.bio`
```bio
experiment HelloBiology {
    // A simple container holding water at body temperature
    environment {
        temperature = 310 K
        pressure = 1 atm
        solvent = "water"
    }

    simulate {
        time = 100 ns
    }

    measure {
        density()
    }
}
```

### `simple_atom.bio`
```bio
experiment SingleIonDynamics {
    ion sodium = Ion("Na+")
    
    environment {
        temperature = 298 K
        volume = 1000 nm3
    }
    
    simulate {
        time = 10 ps
        step = 2 fs
    }
    
    measure {
        kinetic_energy(sodium)
        position(sodium)
    }
    
    visualize {
        render(sodium, style="sphere")
    }
}
```

---

## 7. Future Language Primitives (Not in v0.1)

To fully capture biological complexity, future iterations of BioForge will introduce:
```
hypothesis, mechanism, prediction, perturbation,
observable, comparison, reaction, interaction,
force, field, concentration, signal, pathway,
constraint, scale, counterfactual
```

---

## 8. Fundamental Language Primitives (Long-Term Vision)

The ultimate architecture of BioForge will treat the following concepts as foundational types and blocks:
`ENTITY`, `STRUCTURE`, `STATE`, `FIELD`, `FORCE`, `INTERACTION`, `REACTION`, `CONSTRAINT`, `ENVIRONMENT`, `MECHANISM`, `SCALE`, `TIME`, `EXPERIMENT`, `OBSERVABLE`, `MEASUREMENT`, `HYPOTHESIS`, `PREDICTION`

---

## 9. Native Unit System (Future Phase 2)

A core feature planned for Phase 2 is the `Quantity` struct, which guarantees dimensional correctness.
Every physical value in BioForge will be internally represented as a `Quantity { value, unit, dimension }`. 
The compiler will track dimensions (e.g., `[Mass] * [Length]^2 / [Time]^2` for Energy) and reject mathematically invalid operations at compile-time.

---

## 10. Biological Type System (Future)

BioForge will feature a strictly typed hierarchy representing biological reality:
* `Atom` (Base particle with mass, charge, element)
* `Bond` (Interaction constraint between Atoms)
* `Molecule` (Graph of Atoms and Bonds)
* `Residue` (Subunit of biopolymers)
* `Protein` / `NucleicAcid` (Complex polymers with secondary/tertiary structures)
* `Membrane` / `Organelle` (Mesoscale structures) 
* `Cell` (Macro-container with distinct internal state)
