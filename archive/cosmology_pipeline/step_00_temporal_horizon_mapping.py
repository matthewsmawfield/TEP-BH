#!/usr/bin/env python3
"""Step 00: Temporal Horizon Mapping - map a_eff → 0 to A(φ) → 0.

This step establishes the core temporal-horizon mapping: the apparent FLRW
Big Bang singularity (a_eff → 0) corresponds to the temporal-horizon limit
(A(φ) → 0) in the causal matter frame.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
    conformal_factor_A, effective_scale_factor, temporal_horizon_limit,
    DEFAULT_EPSILON_T, DEFAULT_Z_T, DEFAULT_N_T
)

STEP_ID = "step_00_temporal_horizon_mapping"


def compute_temporal_shear_integral(z_array: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                                     z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute the temporal shear integral ∫ Σ_∥ dℓ along light path.
    
    In the homogeneous limit, this reduces to: ε_t * ln(1+z) * S(z)
    where S(z) = exp(-(z/z_t)^n_t) is the screening function.
    """
    S = np.exp(-(z_array / z_t) ** n_t)
    integral = epsilon_t * np.log(1 + z_array) * S
    return integral


def test_horizon_mapping(z_max: float = 1000.0, n_points: int = 1000,
                        epsilon_t: float = DEFAULT_EPSILON_T,
                        z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test the temporal-horizon mapping across redshift range."""
    
    z_array = np.linspace(0, z_max, n_points)
    
    # Compute conformal factor A(z)
    A_z = conformal_factor_A(z_array, epsilon_t, z_t, n_t)
    
    # Compute effective scale factor a_eff(z)
    a_eff_z = effective_scale_factor(z_array, epsilon_t, z_t, n_t)
    
    # Compute temporal shear integral
    shear_integral = compute_temporal_shear_integral(z_array, epsilon_t, z_t, n_t)
    
    # Test horizon limit
    horizon_approach = temporal_horizon_limit(z_array, epsilon_t, z_t, n_t)
    
    # Find where a_eff → 0 (temporal horizon)
    horizon_indices = np.where(horizon_approach)[0]
    horizon_z = z_array[horizon_indices[0]] if len(horizon_indices) > 0 else z_max
    
    # Find where A → 0 (use a more realistic threshold)
    A_threshold = 0.01  # A < 0.01 considered near horizon
    A_indices = np.where(A_z < A_threshold)[0]
    A_horizon_z = z_array[A_indices[0]] if len(A_indices) > 0 else z_max
    
    results = {
        'step': STEP_ID,
        'description': 'Temporal horizon mapping: a_eff → 0 ↔ A(φ) → 0',
        'parameters': {
            'epsilon_t': epsilon_t,
            'z_t': z_t,
            'n_t': n_t,
            'z_max': z_max,
            'n_points': n_points
        },
        'horizon_mapping': {
            'a_eff_at_z0': rounded(a_eff_z[0], 10),
            'A_at_z0': rounded(A_z[0], 10),
            'a_eff_at_zmax': rounded(a_eff_z[-1], 15),
            'A_at_zmax': rounded(A_z[-1], 15),
            'shear_integral_at_zmax': rounded(shear_integral[-1], 6),
            'horizon_z_a_eff': rounded(horizon_z, 2) if horizon_z is not None else None,
            'horizon_z_A': rounded(A_horizon_z, 2) if A_horizon_z is not None else None,
            'horizons_coincide': abs(horizon_z - A_horizon_z) < 0.1 if (horizon_z is not None and A_horizon_z is not None) else None
        },
        'mapping_valid': (horizon_z is not None and A_horizon_z is not None and 
                         abs(horizon_z - A_horizon_z) < 0.1),
        'interpretation': 'a_eff → 0 and A → 0 occur at same redshift within tolerance' if 
                          (horizon_z is not None and A_horizon_z is not None and abs(horizon_z - A_horizon_z) < 0.1) 
                          else 'Horizon mapping requires further analysis'
    }
    
    # Generate CSV output for plotting
    csv_rows = []
    for i, z in enumerate(z_array):
        csv_rows.append({
            'z': rounded(z, 4),
            'A': rounded(A_z[i], 10),
            'a_eff': rounded(a_eff_z[i], 10),
            'shear_integral': rounded(shear_integral[i], 10),
            'horizon_approach': bool(horizon_approach[i])
        })
    
    write_csv(step_csv_path(STEP_ID), csv_rows)
    
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Test with default parameters
    results = test_horizon_mapping()
    
    # Test with different epsilon_t values to show parameter dependence
    epsilon_tests = [0.0001, 0.001, 0.01, 0.1]
    epsilon_results = []
    
    for eps in epsilon_tests:
        eps_result = test_horizon_mapping(epsilon_t=eps)
        epsilon_results.append({
            'epsilon_t': eps,
            'horizon_z': eps_result['horizon_mapping']['horizon_z_a_eff'],
            'mapping_valid': eps_result['mapping_valid']
        })
    
    results['epsilon_sensitivity'] = epsilon_results
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    print_status(f"Temporal horizon mapping: a_eff → 0 ↔ A → 0 at z ≈ {results['horizon_mapping']['horizon_z_a_eff']}", "INFO")
    
    return results


if __name__ == "__main__":
    run()
