//! Molecular Mechanics to Reaction Kinetics Scale Bridge (Eyring Transition State & Binding Free Energies).

use crate::error::MultiscaleError;

/// Boltzmann Constant $k_B = 1.380649 \times 10^{-23}\text{ J/K}$ (CODATA 2018).
pub const BOLTZMANN_CONSTANT: f64 = 1.380_649e-23;

/// Planck Constant $h = 6.62607015 \times 10^{-34}\text{ J}\cdot\text{s}$ (CODATA 2018).
pub const PLANCK_CONSTANT: f64 = 6.626_070_15e-34;

/// Molar Gas Constant $R = 8.314462618\text{ J}/(\text{mol}\cdot\text{K})$.
pub const GAS_CONSTANT_R: f64 = 8.314_462_618;

/// Calculate macroscopic catalytic turnover rate constant $k_{\text{cat}}$ in $\text{s}^{-1}$
/// from atomistic activation free energy $\Delta G^\ddagger$ using Eyring-Polanyi Transition State Theory:
///
/// $$k_{\text{cat}} = \kappa \frac{k_B T}{h} \exp\left( -\frac{\Delta G^\ddagger}{R T} \right)$$
pub fn eyring_catalytic_rate_constant(
    delta_g_activation_kj_mol: f64,
    temp_k: f64,
    transmission_coeff: f64,
) -> Result<f64, MultiscaleError> {
    if temp_k <= 0.0 {
        return Err(MultiscaleError::InvalidThermodynamicParameter {
            param: "temperature_k".to_string(),
            value: temp_k,
        });
    }

    let kappa = transmission_coeff.clamp(0.0, 1.0);
    let kb_t_over_h = (BOLTZMANN_CONSTANT * temp_k) / PLANCK_CONSTANT;
    let exponent = -(delta_g_activation_kj_mol * 1000.0) / (GAS_CONSTANT_R * temp_k);

    Ok(kappa * kb_t_over_h * exponent.exp())
}

/// Calculate macroscopic thermodynamic dissociation constant $K_d$ in Molar ($\text{M}$)
/// from standard binding free energy $\Delta G_{\text{bind}}^\circ$ in $\text{kJ/mol}$:
///
/// $$K_d = \exp\left( \frac{\Delta G_{\text{bind}}^\circ}{R T} \right)$$
pub fn dissociation_constant_from_delta_g(
    delta_g_binding_kj_mol: f64,
    temp_k: f64,
) -> Result<f64, MultiscaleError> {
    if temp_k <= 0.0 {
        return Err(MultiscaleError::InvalidThermodynamicParameter {
            param: "temperature_k".to_string(),
            value: temp_k,
        });
    }

    let exponent = (delta_g_binding_kj_mol * 1000.0) / (GAS_CONSTANT_R * temp_k);
    Ok(exponent.exp())
}

/// Calculate complete forward and reverse binding kinetics $(K_d, k_{\text{off}})$ given $\Delta G_{\text{bind}}^\circ$ and association rate $k_{\text{on}}$:
///
/// $$k_{\text{off}} = k_{\text{on}} \cdot K_d$$
pub fn binding_kinetics_from_delta_g(
    delta_g_binding_kj_mol: f64,
    k_on_molar_s: f64,
    temp_k: f64,
) -> Result<(f64, f64), MultiscaleError> {
    let kd = dissociation_constant_from_delta_g(delta_g_binding_kj_mol, temp_k)?;
    let k_off = k_on_molar_s.max(0.0) * kd;
    Ok((kd, k_off))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eyring_transition_state_rate_at_body_temperature() {
        // Delta G# = 65.0 kJ/mol at 310.15 K (37°C), kappa = 1.0
        let kcat = eyring_catalytic_rate_constant(65.0, 310.15, 1.0).unwrap();
        // Theoretical: 73.03 s^-1
        assert!(
            (kcat - 73.03).abs() < 0.1,
            "expected kcat ~ 73.03 s^-1, got {}",
            kcat
        );
    }

    #[test]
    fn test_thermodynamic_binding_affinity_to_kd() {
        // Delta G_bind = -30.0 kJ/mol at 310.15 K => Kd = exp(-30000 / (8.314 * 310.15)) = 8.87 uM
        let kd = dissociation_constant_from_delta_g(-30.0, 310.15).unwrap();
        assert!((kd - 8.87e-6).abs() < 0.1e-6, "expected Kd ~ 8.87 uM, got {}", kd);

        // With k_on = 1e6 M^-1 s^-1 => k_off = 1e6 * 8.87e-6 = 8.87 s^-1
        let (kd_calc, k_off) = binding_kinetics_from_delta_g(-30.0, 1.0e6, 310.15).unwrap();
        assert_eq!(kd, kd_calc);
        assert!((k_off - 8.87).abs() < 0.1);
    }
}
