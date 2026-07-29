#!/usr/bin/env python3
"""Step 03: Effective Stress-Energy and Singularity-Theorem Evasion.

Compute the temporal-field stress-energy tensor T^(φ)_μν and evaluate
energy conditions to identify why standard FLRW singularity theorems
do not apply to the temporal-horizon geometry.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
    conformal_factor_A, DEFAULT_EPSILON_T, DEFAULT_Z_T, DEFAULT_N_T,
    H0_KM_S_MPC, OMEGA_M, OMEGA_L
)

STEP_ID = "step_03_effective_stress_energy"


def temporal_field_energy_density(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                                   z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute effective energy density ρ_φ of temporal field.
    
    For conformal coupling, the effective energy density from temporal shear:
    ρ_φ ≈ (ε_t / 8πG) * H² * (d ln A/d ln a)²
    """
    from th_common import e_z
    H = H0_KM_S_MPC * e_z(z)
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    
    # Derivative of ln A with respect to ln a
    # d ln A / d ln a = (dA/dz * dz/da) * (a/A)
    # Simplified: ~ ε_t * S(z) * (1 - n_t * (z/z_t)^n_t * ln(1+z))
    S = np.exp(-(z / z_t) ** n_t)
    dlnA_dlna = epsilon_t * S * (1 - n_t * (z / z_t)**n_t * np.log(1 + z))
    
    # Energy density (in units where 8πG = 1 for simplicity)
    rho_phi = H**2 * dlnA_dlna**2
    return rho_phi


def temporal_field_pressure(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                            z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute effective pressure p_φ of temporal field.
    
    For conformal coupling, equation of state is typically w ≈ -1 at late times.
    """
    rho_phi = temporal_field_energy_density(z, epsilon_t, z_t, n_t)
    
    # Equation of state parameter w_φ
    # For temporal shear, w ≈ -1 + small corrections
    w_phi = -1.0 + 0.1 * (1 + np.asarray(z))**(-1)  # Simplified
    
    p_phi = w_phi * rho_phi
    return p_phi


def equation_of_state(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                     z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute equation of state parameter w_φ = p_φ/ρ_φ."""
    rho_phi = temporal_field_energy_density(z, epsilon_t, z_t, n_t)
    p_phi = temporal_field_pressure(z, epsilon_t, z_t, n_t)
    
    w_phi = np.divide(p_phi, rho_phi, out=np.zeros_like(p_phi), where=rho_phi!=0)
    return w_phi


def energy_conditions(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                      z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Evaluate standard energy conditions for temporal field.
    
    Tests:
    - Null Energy Condition (NEC): ρ + p ≥ 0
    - Weak Energy Condition (WEC): ρ ≥ 0 and ρ + p ≥ 0
    - Dominant Energy Condition (DEC): ρ ≥ |p|
    - Strong Energy Condition (SEC): ρ + 3p ≥ 0 and ρ + p ≥ 0
    """
    rho = temporal_field_energy_density(z, epsilon_t, z_t, n_t)
    p = temporal_field_pressure(z, epsilon_t, z_t, n_t)
    
    rho_plus_p = rho + p
    rho_plus_3p = rho + 3 * p
    
    nec_satisfied = np.all(rho_plus_p >= -1e-10)
    wec_satisfied = np.all(rho >= -1e-10) and nec_satisfied
    dec_satisfied = np.all(rho >= np.abs(p) - 1e-10)
    sec_satisfied = np.all(rho_plus_3p >= -1e-10) and nec_satisfied
    
    return {
        'NEC': nec_satisfied,
        'WEC': wec_satisfied,
        'DEC': dec_satisfied,
        'SEC': sec_satisfied,
        'rho_plus_p_min': rounded(np.min(rho_plus_p), 10),
        'rho_plus_3p_min': rounded(np.min(rho_plus_3p), 10),
        'rho_min': rounded(np.min(rho), 10)
    }


def test_stress_energy_regularity(z_max: float = 1000.0, n_points: int = 1000,
                                 epsilon_t: float = DEFAULT_EPSILON_T,
                                 z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test stress-energy tensor regularity and energy conditions."""
    
    z_array = np.linspace(0, z_max, n_points)
    
    # Compute stress-energy components
    rho_phi = temporal_field_energy_density(z_array, epsilon_t, z_t, n_t)
    p_phi = temporal_field_pressure(z_array, epsilon_t, z_t, n_t)
    w_phi = equation_of_state(z_array, epsilon_t, z_t, n_t)
    A = conformal_factor_A(z_array, epsilon_t, z_t, n_t)
    
    # Check regularity
    rho_finite = np.all(np.isfinite(rho_phi))
    p_finite = np.all(np.isfinite(p_phi))
    w_finite = np.all(np.isfinite(w_phi))
    
    # Check behavior near horizon
    high_z_mask = z_array > z_max * 0.9
    rho_high_z = rho_phi[high_z_mask]
    p_high_z = p_phi[high_z_mask]
    
    rho_bounded = np.all(np.abs(rho_high_z) < 1e20) if len(rho_high_z) > 0 else True
    p_bounded = np.all(np.abs(p_high_z) < 1e20) if len(p_high_z) > 0 else True
    
    # Evaluate energy conditions
    ec_results = energy_conditions(z_array, epsilon_t, z_t, n_t)
    
    results = {
        'step': STEP_ID,
        'description': 'Effective stress-energy and energy condition analysis',
        'parameters': {
            'epsilon_t': epsilon_t,
            'z_t': z_t,
            'n_t': n_t,
            'z_max': z_max,
            'n_points': n_points
        },
        'stress_energy_components': {
            'rho_finite': rho_finite,
            'p_finite': p_finite,
            'w_finite': w_finite,
            'rho_bounded_near_horizon': rho_bounded,
            'p_bounded_near_horizon': p_bounded,
            'rho_max': rounded(np.max(rho_phi), 6) if rho_finite else None,
            'p_max': rounded(np.max(p_phi), 6) if p_finite else None,
            'w_at_z0': rounded(w_phi[0], 6) if w_finite else None,
            'w_at_zmax': rounded(w_phi[-1], 6) if w_finite else None
        },
        'energy_conditions': ec_results,
        'singularity_theorem_analysis': {
            'SEC_violated': not ec_results['SEC'],
            'NEC_satisfied': ec_results['NEC'],
            'standard_theorem_applies': ec_results['SEC'],  # SEC required for Hawking-Penrose
            'TEP_evades_singularity_theorem': not ec_results['SEC']
        },
        'interpretation': 'Temporal field violates SEC, evading standard singularity theorems' if 
                          not ec_results['SEC'] else
                          'Energy conditions satisfied; singularity theorems may apply'
    }
    
    # Generate CSV output
    csv_rows = []
    for i, z in enumerate(z_array):
        csv_rows.append({
            'z': rounded(z, 4),
            'A': rounded(A[i], 10),
            'rho_phi': rounded(rho_phi[i], 10) if np.isfinite(rho_phi[i]) else None,
            'p_phi': rounded(p_phi[i], 10) if np.isfinite(p_phi[i]) else None,
            'w_phi': rounded(w_phi[i], 6) if np.isfinite(w_phi[i]) else None
        })
    
    write_csv(step_csv_path(STEP_ID), csv_rows)
    
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Test with default parameters
    results = test_stress_energy_regularity()
    
    # Test with different epsilon_t values
    epsilon_tests = [0.0001, 0.001, 0.01]
    epsilon_results = []
    
    for eps in epsilon_tests:
        eps_result = test_stress_energy_regularity(epsilon_t=eps)
        epsilon_results.append({
            'epsilon_t': eps,
            'SEC_violated': eps_result['singularity_theorem_analysis']['SEC_violated'],
            'evades_singularity_theorem': eps_result['singularity_theorem_analysis']['TEP_evades_singularity_theorem']
        })
    
    results['epsilon_sensitivity'] = epsilon_results
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    print_status(f"SEC violated: {results['singularity_theorem_analysis']['SEC_violated']}", "INFO")
    print_status(f"Evades singularity theorem: {results['singularity_theorem_analysis']['TEP_evades_singularity_theorem']}", "INFO")
    
    return results


if __name__ == "__main__":
    run()
