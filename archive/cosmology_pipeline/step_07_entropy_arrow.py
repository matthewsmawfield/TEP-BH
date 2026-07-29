#!/usr/bin/env python3
"""Step 07: Entropy and Arrow of Time.

Track entropy in the temporal-horizon reconstruction to test whether
the temporal horizon is thermodynamically regular and does not require
an infinite-entropy initial singularity.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
    conformal_factor_A, DEFAULT_EPSILON_T, DEFAULT_Z_T, DEFAULT_N_T,
    OMEGA_M
)

STEP_ID = "step_07_entropy_arrow"


def radiation_entropy_density(z: np.ndarray) -> np.ndarray:
    """Compute radiation entropy density s_γ ∝ T³.
    
    For photons: s_γ = (4/3) ρ_γ / T ∝ T³
    """
    T_CMB = 2.725  # K
    T_z = T_CMB * (1 + z)
    s_gamma = 2.89e6 * (T_z / 2.725)**3  # Simplified units
    return s_gamma


def matter_entropy_density(z: np.ndarray) -> np.ndarray:
    """Compute matter entropy density (simplified).
    
    For non-relativistic matter, entropy is much smaller than radiation.
    """
    s_matter = 1e4 * (1 + z)**3  # Simplified
    return s_matter


def total_entropy_density(z: np.ndarray) -> np.ndarray:
    """Compute total entropy density."""
    s_gamma = radiation_entropy_density(z)
    s_matter = matter_entropy_density(z)
    return s_gamma + s_matter


def tep_entropy_density(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                       z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute entropy density in TEP matter frame.
    
    In conformal frame: s̃ = s / A³ (entropy scales with volume)
    """
    s_standard = total_entropy_density(z)
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    s_tep = s_standard / A**3
    return s_tep


def comoving_entropy_standard(z: np.ndarray) -> np.ndarray:
    """Standard comoving entropy: s a³ = constant."""
    s = total_entropy_density(z)
    a = 1.0 / (1 + z)
    s_comoving = s * a**3
    return s_comoving


def comoving_entropy_tep(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                         z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """TEP comoving entropy in the matter frame: s̃ (ã)^3 with ã = A a."""
    s_tilde = tep_entropy_density(z, epsilon_t, z_t, n_t)
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    a = 1.0 / (1.0 + z)
    a_tilde = A * a
    return s_tilde * a_tilde**3


def entropy_conservation_test(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                             z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test entropy conservation in standard and TEP frames."""
    
    # Standard comoving entropy
    s_comoving_std = comoving_entropy_standard(z)
    s_comoving_std_ratio = s_comoving_std / s_comoving_std[0]
    
    # TEP comoving entropy
    s_comoving_tep = comoving_entropy_tep(z, epsilon_t, z_t, n_t)
    s_comoving_tep_ratio = s_comoving_tep / s_comoving_tep[0]
    
    # Check conservation
    std_conserved = np.all(np.abs(s_comoving_std_ratio - 1.0) < 0.01)
    tep_conserved = np.all(np.abs(s_comoving_tep_ratio - 1.0) < 0.01)
    
    return {
        'standard_comoving_conserved': std_conserved,
        'tep_comoving_conserved': tep_conserved,
        's_comoving_std_range': (rounded(np.min(s_comoving_std_ratio), 6), 
                                  rounded(np.max(s_comoving_std_ratio), 6)),
        's_comoving_tep_range': (rounded(np.min(s_comoving_tep_ratio), 6),
                                 rounded(np.max(s_comoving_tep_ratio), 6))
    }


def arrow_of_time_indicator(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                           z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Arrow-of-time indicator using matter-frame comoving entropy.

    For adiabatic evolution, comoving entropy should be conserved (non-decreasing).
    """

    S_comoving = comoving_entropy_tep(z, epsilon_t, z_t, n_t)
    dS_dz = np.gradient(S_comoving, z)

    # Future direction is decreasing z; require S not to decrease toward the future
    # i.e. dS/dz should be <= small tolerance
    tol = 1e-6 * max(1.0, float(abs(S_comoving[0])))
    arrow_points_forward = bool(np.all(dS_dz <= tol))

    return {
        'arrow_points_forward': arrow_points_forward,
        'ds_dz_at_z0': rounded(float(dS_dz[0]), 6),
        'ds_dz_at_zmax': rounded(float(dS_dz[-1]), 6),
        'entropy_monotonic': arrow_points_forward
    }


def test_entropy_regularity(z_max: float = 1000.0, n_points: int = 1000,
                           epsilon_t: float = DEFAULT_EPSILON_T,
                           z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test entropy regularity at temporal horizon."""
    
    z_array = np.linspace(0, z_max, n_points)
    
    # Compute entropy densities
    s_std = total_entropy_density(z_array)
    s_tep = tep_entropy_density(z_array, epsilon_t, z_t, n_t)
    
    # Check for divergences
    s_std_finite = np.all(np.isfinite(s_std))
    s_tep_finite = np.all(np.isfinite(s_tep))
    
    # Check behavior near horizon
    high_z_mask = z_array > z_max * 0.9
    s_tep_high_z = s_tep[high_z_mask]
    
    s_tep_bounded = np.all(np.abs(s_tep_high_z) < 1e30) if len(s_tep_high_z) > 0 else True
    s_tep_positive = np.all(s_tep_high_z > 0) if len(s_tep_high_z) > 0 else True
    
    # Test entropy conservation
    conservation_test = entropy_conservation_test(z_array, epsilon_t, z_t, n_t)
    
    # Test arrow of time
    arrow_test = arrow_of_time_indicator(z_array, epsilon_t, z_t, n_t)
    
    results = {
        'step': STEP_ID,
        'description': 'Entropy and arrow-of-time regularity at temporal horizon',
        'parameters': {
            'epsilon_t': epsilon_t,
            'z_t': z_t,
            'n_t': n_t,
            'z_max': z_max,
            'n_points': n_points
        },
        'entropy_regularity': {
            's_standard_finite': s_std_finite,
            's_tep_finite': s_tep_finite,
            's_tep_bounded_near_horizon': s_tep_bounded,
            's_tep_positive_near_horizon': s_tep_positive,
            's_tep_max': rounded(np.max(s_tep), 6) if s_tep_finite else None,
            's_tep_min': rounded(np.min(s_tep), 6) if s_tep_finite else None
        },
        'entropy_conservation': conservation_test,
        'arrow_of_time': arrow_test,
        'thermodynamic_regularity': {
            'entropy_finite': s_tep_finite,
            'entropy_bounded': s_tep_bounded,
            'entropy_positive': s_tep_positive,
            'entropy_conserved': conservation_test['tep_comoving_conserved'],
            'arrow_of_time_valid': arrow_test['arrow_points_forward'],
            'all_regular': (s_tep_finite and s_tep_bounded and s_tep_positive and 
                           conservation_test['tep_comoving_conserved'] and 
                           arrow_test['arrow_points_forward'])
        },
        'interpretation': 'Temporal horizon is thermodynamically regular with finite entropy' if
                          (s_tep_finite and s_tep_bounded and s_tep_positive) else
                          'Entropy behavior requires further analysis'
    }
    
    # Generate CSV output
    csv_rows = []
    for i, z in enumerate(z_array):
        csv_rows.append({
            'z': rounded(z, 4),
            's_standard': rounded(s_std[i], 6),
            's_tep': rounded(s_tep[i], 6),
            's_comoving_std': rounded(conservation_test.get('s_comoving_std_range', (0, 0))[0], 6) if i == 0 else None,
            's_comoving_tep': rounded(conservation_test.get('s_comoving_tep_range', (0, 0))[0], 6) if i == 0 else None
        })
    
    write_csv(step_csv_path(STEP_ID), csv_rows)
    
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Test with default parameters
    results = test_entropy_regularity()
    
    # Test with different epsilon_t values
    epsilon_tests = [0.0001, 0.001, 0.01]
    epsilon_results = []
    
    for eps in epsilon_tests:
        eps_result = test_entropy_regularity(epsilon_t=eps)
        epsilon_results.append({
            'epsilon_t': eps,
            'all_regular': eps_result['thermodynamic_regularity']['all_regular'],
            'entropy_finite': eps_result['entropy_regularity']['s_tep_finite']
        })
    
    results['epsilon_sensitivity'] = epsilon_results
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    print_status(f"Thermodynamic regularity: {results['thermodynamic_regularity']['all_regular']}", "INFO")
    print_status(f"Arrow of time valid: {results['arrow_of_time']['arrow_points_forward']}", "INFO")
    
    return results


if __name__ == "__main__":
    run()
