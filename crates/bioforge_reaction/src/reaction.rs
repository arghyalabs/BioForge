//! Reaction networks, stoichiometry matrices, and ODE derivative evaluations.

use serde::{Deserialize, Serialize};

use crate::error::ReactionError;
use crate::ratelaw::RateLaw;
use crate::species::Species;

/// A single biochemical reaction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reaction {
    /// Unique reaction index.
    pub id: usize,
    /// Human-readable reaction descriptor (e.g. "ATP_hydrolysis", "Phosphorylation").
    pub name: String,
    /// Reactant species indices and stoichiometric coefficients: `(species_idx, coeff)`.
    pub reactants: Vec<(usize, f64)>,
    /// Product species indices and stoichiometric coefficients: `(species_idx, coeff)`.
    pub products: Vec<(usize, f64)>,
    /// The kinetic rate law governing this reaction.
    pub rate_law: RateLaw,
}

/// A complete biochemical reaction network.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct ReactionNetwork {
    /// Species in the network.
    pub species: Vec<Species>,
    /// Reactions in the network.
    pub reactions: Vec<Reaction>,
}

impl ReactionNetwork {
    /// Create an empty reaction network.
    #[must_use]
    pub fn new() -> Self {
        Self {
            species: Vec::new(),
            reactions: Vec::new(),
        }
    }

    /// Add a species to the network and return its assigned integer index.
    pub fn add_species(&mut self, name: impl Into<String>, initial_concentration: f64) -> usize {
        let id = self.species.len();
        self.species.push(Species::new(id, name, initial_concentration));
        id
    }

    /// Add a reaction to the network.
    pub fn add_reaction(
        &mut self,
        name: impl Into<String>,
        reactants: Vec<(usize, f64)>,
        products: Vec<(usize, f64)>,
        rate_law: RateLaw,
    ) -> usize {
        let id = self.reactions.len();
        self.reactions.push(Reaction {
            id,
            name: name.into(),
            reactants,
            products,
            rate_law,
        });
        id
    }

    /// Number of species in the network.
    #[must_use]
    pub fn num_species(&self) -> usize {
        self.species.len()
    }

    /// Number of reactions in the network.
    #[must_use]
    pub fn num_reactions(&self) -> usize {
        self.reactions.len()
    }

    /// Build the stoichiometry matrix $\mathbf{N} \in \mathbb{R}^{S \times R}$ where $\mathbf{N}_{ij} = \nu_{ij}^+ - \nu_{ij}^-$.
    #[must_use]
    pub fn stoichiometry_matrix(&self) -> Vec<Vec<f64>> {
        let num_s = self.num_species();
        let num_r = self.num_reactions();
        let mut mat = vec![vec![0.0; num_r]; num_s];

        for (r_idx, rxn) in self.reactions.iter().enumerate() {
            // Subtract reactants
            for &(s_idx, coeff) in &rxn.reactants {
                if s_idx < num_s {
                    mat[s_idx][r_idx] -= coeff;
                }
            }
            // Add products
            for &(s_idx, coeff) in &rxn.products {
                if s_idx < num_s {
                    mat[s_idx][r_idx] += coeff;
                }
            }
        }

        mat
    }

    /// Compute the time derivatives of all species concentrations $\frac{d\vec{C}}{dt} = \mathbf{N} \cdot \vec{v}(\vec{C})$.
    #[must_use]
    pub fn compute_derivatives(&self, concentrations: &[f64]) -> Vec<f64> {
        let num_s = self.num_species();
        let mut dcdt = vec![0.0; num_s];

        for rxn in &self.reactions {
            let v = rxn
                .rate_law
                .evaluate_velocity(concentrations, &rxn.reactants, &rxn.products);

            for &(s_idx, coeff) in &rxn.reactants {
                if s_idx < num_s {
                    dcdt[s_idx] -= coeff * v;
                }
            }
            for &(s_idx, coeff) in &rxn.products {
                if s_idx < num_s {
                    dcdt[s_idx] += coeff * v;
                }
            }
        }

        dcdt
    }

    /// Compute stochastic reaction propensities $a_j(\vec{X})$ for all reactions in the network.
    #[must_use]
    pub fn compute_propensities(&self, counts: &[u64]) -> Vec<f64> {
        let volume = self
            .species
            .first()
            .map(|s| s.compartment_volume)
            .unwrap_or(1e-15);

        self.reactions
            .iter()
            .map(|rxn| rxn.rate_law.evaluate_propensity(counts, volume, &rxn.reactants))
            .collect()
    }

    /// Get initial physical concentration vector $\vec{C}_0$ in Molar ($\text{mol/L}$).
    #[must_use]
    pub fn initial_concentrations(&self) -> Vec<f64> {
        self.species.iter().map(|s| s.initial_concentration).collect()
    }

    /// Get initial discrete molecular counts vector $\vec{X}_0$.
    #[must_use]
    pub fn initial_counts(&self) -> Vec<u64> {
        self.species
            .iter()
            .map(|s| s.to_discrete_count(s.initial_concentration))
            .collect()
    }

    /// Find species index by name.
    pub fn find_species_index(&self, name: &str) -> Result<usize, ReactionError> {
        self.species
            .iter()
            .position(|s| s.name == name)
            .ok_or_else(|| ReactionError::SpeciesNotFound {
                name: name.to_string(),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reaction_network_stoichiometry_and_derivatives() {
        let mut net = ReactionNetwork::new();
        let a = net.add_species("A", 1.0); // 1.0 M
        let b = net.add_species("B", 0.0); // 0.0 M

        // Reaction 1: A -> B with k = 0.5 s^-1
        net.add_reaction(
            "A_to_B",
            vec![(a, 1.0)],
            vec![(b, 1.0)],
            RateLaw::mass_action_forward(0.5),
        );

        assert_eq!(net.num_species(), 2);
        assert_eq!(net.num_reactions(), 1);

        let mat = net.stoichiometry_matrix();
        // Row 0 (A): -1.0, Row 1 (B): +1.0
        assert_eq!(mat[0][0], -1.0);
        assert_eq!(mat[1][0], 1.0);

        // At [A]=1.0, d[A]/dt = -0.5 M/s, d[B]/dt = +0.5 M/s
        let dcdt = net.compute_derivatives(&[1.0, 0.0]);
        assert!((dcdt[0] - (-0.5)).abs() < 1e-12);
        assert!((dcdt[1] - 0.5).abs() < 1e-12);
    }
}
