#!/usr/bin/env python3
"""Step 01: Matter-Frame Geometry and Curvature Regularity.

Compute physical matter-frame curvature invariants to test whether the
causal matter frame remains regular at the temporal horizon (A → 0).
Tests: R̃, R̃_μν R̃^μν, R̃_μνρσ R̃^μνρσ, det(g̃_μν).
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
    conformal_factor_A, DEFAULT_EPSILON_T, DEFAULT_Z_T, DEFAULT_N_T,
    H0_KM_S_MPC, OMEGA_M, OMEGA_L, OMEGA_K, C_KM_S
)

STEP_ID = "step_01_matter_frame_curvature"


def flrw_scale_factor(z: np.ndarray) -> np.ndarray:
    """Standard FLRW scale factor a(z) = 1/(1+z)."""
    return 1.0 / (1 + np.asarray(z))


def flrw_hubble_parameter(z: np.ndarray) -> np.ndarray:
    """FLRW Hubble parameter H(z) = H0 * E(z)."""
    from th_common import e_z
    return H0_KM_S_MPC * e_z(z)


def flrw_ricci_scalar(z: np.ndarray) -> np.ndarray:
    """FLRW Ricci scalar R(z) for flat universe."""
    z = np.asarray(z)
    H = flrw_hubble_parameter(z)
    dH_dz = np.gradient(H, z)
    H_dot = -(1.0 + z) * H * dH_dz
    return 6.0 * (H_dot + 2.0 * H**2)


def matter_frame_metric(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                        z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute matter-frame metric g̃_μν = A²(φ) g_μν in pure conformal closure."""
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    # In FLRW, g_μν is diagonal with spatial components a²(t)
    # Matter frame: g̃_μν = A² * g_μν
    # For curvature invariants, we need the conformal factor
    return A

def disformal_lapse(z: np.ndarray, A: np.ndarray, power: float = 20.0, B0: float = 1.0) -> np.ndarray:
    z = np.asarray(z)
    A = np.asarray(A)
    A_safe = np.maximum(A, 1e-30)
    return np.sqrt(A_safe**2 + (B0**2) * (A_safe ** (-2.0 * power)))


def matter_frame_hubble(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                         z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T,
                         use_disformal: bool = True,
                         disformal_power: float = 20.0,
                         disformal_B0: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    z = np.asarray(z)
    H = flrw_hubble_parameter(z)
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    dlnA_dz = np.gradient(np.log(np.maximum(A, 1e-30)), z)
    dlnA_dt = -(1.0 + z) * H * dlnA_dz
    H_geom = H + dlnA_dt
    N = disformal_lapse(z, A, disformal_power, disformal_B0) if use_disformal else A
    H_tilde = H_geom / np.maximum(N, 1e-30)
    dHdz = np.gradient(H_tilde, z)
    Hdot_tilde = -(1.0 + z) * H * dHdz / np.maximum(N, 1e-30)
    a = flrw_scale_factor(z)
    a_tilde = A * a
    return H_tilde, Hdot_tilde, a_tilde, N



def matter_frame_ricci_scalar(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                              z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute matter-frame Ricci scalar from effective H̃ and Ḣ̃."""
    z = np.asarray(z)
    H_tilde, Hdot_tilde, _, _ = matter_frame_hubble(z, epsilon_t, z_t, n_t)
    return 6.0 * (Hdot_tilde + 2.0 * H_tilde**2)


def matter_frame_kretschmann(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                              z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute Kretschmann scalar in matter frame from H̃ and Ḣ̃."""
    z = np.asarray(z)
    H_tilde, Hdot_tilde, _, _ = matter_frame_hubble(z, epsilon_t, z_t, n_t)
    return 12.0 * (Hdot_tilde**2 + 3.0 * H_tilde**2 * Hdot_tilde + 3.0 * H_tilde**4)


def metric_determinant(z: np.ndarray, epsilon_t: float = DEFAULT_EPSILON_T,
                       z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> np.ndarray:
    """Compute metric determinant det(g̃_μν) including lapse regularization."""
    z = np.asarray(z)
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    _, _, a_tilde, N = matter_frame_hubble(z, epsilon_t, z_t, n_t)
    return - (np.maximum(N, 1e-30)**2) * (np.maximum(a_tilde, 1e-30)**6)


def test_curvature_regularity(z_max: float = 1000.0, n_points: int = 1000,
                             epsilon_t: float = DEFAULT_EPSILON_T,
                             z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test curvature invariants for regularity at temporal horizon."""
    
    z_array = np.linspace(0, z_max, n_points)
    
    # Compute curvature invariants
    R_tilde = matter_frame_ricci_scalar(z_array, epsilon_t, z_t, n_t)
    K_tilde = matter_frame_kretschmann(z_array, epsilon_t, z_t, n_t)
    det_g = metric_determinant(z_array, epsilon_t, z_t, n_t)
    A = conformal_factor_A(z_array, epsilon_t, z_t, n_t)
    
    # Check for divergences
    R_finite = np.all(np.isfinite(R_tilde))
    K_finite = np.all(np.isfinite(K_tilde))
    det_finite = np.all(np.isfinite(det_g))
    
    # Check behavior near horizon (high z)
    high_z_mask = z_array > z_max * 0.9
    R_high_z = R_tilde[high_z_mask]
    K_high_z = K_tilde[high_z_mask]
    det_high_z = det_g[high_z_mask]
    
    R_bounded = np.all(np.abs(R_high_z) < 1e10) if len(R_high_z) > 0 else True
    K_bounded = np.all(np.abs(K_high_z) < 1e20) if len(K_high_z) > 0 else True
    det_nonzero = np.all(np.abs(det_high_z) > 1e-20) if len(det_high_z) > 0 else True
    
    results = {
        'step': STEP_ID,
        'description': 'Matter-frame curvature regularity at temporal horizon',
        'parameters': {
            'epsilon_t': epsilon_t,
            'z_t': z_t,
            'n_t': n_t,
            'z_max': z_max,
            'n_points': n_points
        },
        'curvature_invariants': {
            'R_finite': R_finite,
            'K_finite': K_finite,
            'det_finite': det_finite,
            'R_bounded_near_horizon': R_bounded,
            'K_bounded_near_horizon': K_bounded,
            'det_nonzero_near_horizon': det_nonzero,
            'R_max': rounded(np.max(np.abs(R_tilde)), 6) if R_finite else None,
            'K_max': rounded(np.max(np.abs(K_tilde)), 6) if K_finite else None,
            'det_min': rounded(np.min(np.abs(det_g)), 15) if det_finite else None
        },
        'regularity_status': {
            'all_invariants_finite': R_finite and K_finite and det_finite,
            'all_invariants_bounded': R_bounded and K_bounded and det_nonzero,
            'matter_frame_regular': R_finite and K_finite and det_finite and R_bounded and K_bounded and det_nonzero
        },
        'interpretation': 'Matter-frame curvature invariants remain finite and bounded' if 
                          (R_finite and K_finite and det_finite and R_bounded and K_bounded and det_nonzero) else
                          'Matter-frame may require disformal/topological regularization'
    }
    
    # Generate CSV output
    csv_rows = []
    for i, z in enumerate(z_array):
        csv_rows.append({
            'z': rounded(z, 4),
            'A': rounded(A[i], 10),
            'R_tilde': rounded(R_tilde[i], 6) if np.isfinite(R_tilde[i]) else None,
            'K_tilde': rounded(K_tilde[i], 6) if np.isfinite(K_tilde[i]) else None,
            'det_g': rounded(det_g[i], 15) if np.isfinite(det_g[i]) else None
        })
    
    write_csv(step_csv_path(STEP_ID), csv_rows)
    
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Test with default parameters
    results = test_curvature_regularity()
    
    # Test with different epsilon_t values
    epsilon_tests = [0.0001, 0.001, 0.01]
    epsilon_results = []
    
    for eps in epsilon_tests:
        eps_result = test_curvature_regularity(epsilon_t=eps)
        epsilon_results.append({
            'epsilon_t': eps,
            'matter_frame_regular': eps_result['regularity_status']['matter_frame_regular']
        })
    
    results['epsilon_sensitivity'] = epsilon_results
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    print_status(f"Matter-frame regularity: {results['regularity_status']['matter_frame_regular']}", "INFO")
    
    return results


if __name__ == "__main__":
    run()
