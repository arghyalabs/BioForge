//! # BioForge Cell (`bioforge_cell`)
//!
//! Cellular systems, organelle compartments, Central Dogma genetics, and gene regulatory networks for BioForge.
//!
//! ## Scientific Architecture (Principle 3 & Principle 12)
//!
//! Models cellular biological processes across multiple scales:
//! - **Genetics & Central Dogma**: 64-codon universal genetic code translation, transcription ($k_{\text{tx}}$) and translation ($k_{\text{tl}}$) kinetics.
//! - **Gene Regulatory Networks**: Hill activation/repression promoters, synthetic bistable Toggle Switch, and the 3-gene Repressilator oscillator.
//! - **Organelle Compartments**: Cytosol, Nucleus, Mitochondria, ER, and trans-compartmental transport flux.
//! - **Signaling Cascades**: Multitier MAPK / ERK phosphorylation cascades.

#![deny(unsafe_code)]
#![allow(non_snake_case)]

pub mod codon;
pub mod compartment;
pub mod dogma;
pub mod error;
pub mod grn;
pub mod signaling;
pub mod solver;

pub use codon::{
    reverse_complement_dna, transcribe_dna_to_rna, translate_dna_to_protein, translate_rna_codon,
    translate_rna_to_protein,
};
pub use compartment::{Compartment, CompartmentTransport, OrganelleKind};
pub use dogma::GeneExpression;
pub use error::CellError;
pub use grn::{PromoterRegulation, Repressilator, ToggleSwitch};
pub use signaling::MapkCascade;
pub use solver::{integrate_rk4, CellSimulationTrajectory};

// ─── Central Dogma Dynamic Simulation Benchmark ────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Full Central Dogma Numerical Dynamic Simulation Benchmark:
    ///
    /// Simulates induction of gene expression from $t=0$ ($[\text{mRNA}]_0 = 0, [\text{Protein}]_0 = 0$)
    /// and verifies that the numerical trajectory asymptotically converges to the exact theoretical solution:
    ///
    /// $$p(t) = p^* \left[ 1 - \frac{\delta_m e^{-\delta_p t} - \delta_p e^{-\delta_m t}}{\delta_m - \delta_p} \right]$$
    #[test]
    fn test_central_dogma_dynamic_convergence_to_steady_state() {
        let gene = GeneExpression::new("GFP", 0.1, 0.002, 0.05, 0.0001);
        let names = vec!["GFP_mRNA".to_string(), "GFP_Protein".to_string()];

        let initial_state = [0.0, 0.0];
        let total_time = 30000.0; // 30,000 seconds (~8.3 hours, ~4.4 protein half-lives)
        let dt = 1.0;             // 1 second time step

        let traj = integrate_rk4(&initial_state, total_time, dt, names, |state, _t| {
            let (dm, dp) = gene.compute_derivatives(state[0], state[1]);
            vec![dm, dp]
        })
        .unwrap();

        let final_mrna = traj.values[0].last().copied().unwrap();
        let final_protein = traj.values[1].last().copied().unwrap();

        // mRNA reaches steady-state rapidly (t_1/2 ~ 5.8 min) -> must match 50.0 nM within 0.01%
        assert!(
            (final_mrna - 50.0).abs() < 0.01,
            "expected mRNA ~ 50.0 nM, got {:.3} nM",
            final_mrna
        );

        // Exact analytical transient solution for coupled 2-ODE linear system
        let dm = gene.mrna_degradation_rate;
        let dp = gene.protein_degradation_rate;
        let p_star = gene.steady_state_protein_nM();
        let expected_prot_at_t = p_star * (1.0 - (dm * (-dp * total_time).exp() - dp * (-dm * total_time).exp()) / (dm - dp));

        assert!(
            (final_protein - expected_prot_at_t).abs() < 1.0,
            "expected protein ~ {:.1} nM, got {:.1} nM",
            expected_prot_at_t,
            final_protein
        );
    }
}
