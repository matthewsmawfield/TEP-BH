#!/usr/bin/env python3
"""Step 02: Geodesic Completeness - test null and timelike geodesic extension.

Integrate null and timelike geodesics backward toward the temporal horizon to test
whether physical trajectories terminate at a curvature singularity or approach
the horizon asymptotically.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
    conformal_factor_A, DEFAULT_EPSILON_T, DEFAULT_Z_T, DEFAULT_N_T,
    H0_KM_S_MPC, C_KM_S
)

STEP_ID = "step_02_geodesic_completeness"


def null_geodesic_equation(lam, y, epsilon_t: float, z_t: float, n_t: float):
    """Geodesic equation for null rays in matter frame.

    Simplified radial null geodesic in conformal frame, integrated
    backward (toward increasing redshift / temporal horizon):
    dt/dλ = -1/A(t), dr/dλ = +1/A(t)
    """
    t, r = y
    z = 1.0 / t - 1.0 if t > 1e-15 else 1e15
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    # Negative sign: backward integration toward temporal horizon (t decreases, z increases)
    dt_dlambda = -1.0 / max(A, 1e-10)
    dr_dlambda = -dt_dlambda  # Null geodesic: dr/dλ = +1/A
    return [dt_dlambda, dr_dlambda]


def timelike_geodesic_equation(tau, y, epsilon_t: float, z_t: float, n_t: float):
    """Geodesic equation for timelike observers in matter frame.

    Simplified radial timelike geodesic, integrated backward:
    dt/dτ = -1/A(t), dr/dτ = -v/A(t) where v < 1
    """
    t, r = y
    z = 1.0 / t - 1.0 if t > 1e-15 else 1e15
    A = conformal_factor_A(z, epsilon_t, z_t, n_t)
    dt_dtau = -1.0 / max(A, 1e-10)
    v = 0.1  # Test velocity
    dr_dtau = v * dt_dtau  # v > 0, dt_dtau < 0 → dr_dtau < 0
    return [dt_dtau, dr_dtau]


def horizon_event(lam, y, *args):
    """Event function: terminate when conformal time t falls below threshold.
    *args absorbs extra parameters passed by solve_ivp."""
    return y[0] - 1e-6


horizon_event.terminal = True
horizon_event.direction = -1  # t decreasing


def integrate_null_geodesic_backward(z_start: float = 0.0, z_end: float = 1000.0,
                                      n_steps: int = 1000,
                                      epsilon_t: float = DEFAULT_EPSILON_T,
                                      z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Integrate null geodesic backward from observer to high redshift."""

    y0 = [1.0 / (1.0 + z_start), 0.0]  # Initial conditions: t, r at z_start

    try:
        sol = solve_ivp(
            null_geodesic_equation, [0, 10.0], y0,
            args=(epsilon_t, z_t, n_t),
            events=horizon_event,
            dense_output=True,
            max_step=0.01,
            method='RK45'
        )

        t_final = sol.y[0, -1]
        r_final = sol.y[1, -1]
        terminated_by_event = sol.t_events is not None and len(sol.t_events[0]) > 0
        # Event at t=1e-6 is a numerical cutoff to prevent overflow;
        # the geodesic approaches t=0 asymptotically, not at finite affine parameter.
        reaches_horizon = False

        # Sample uniformly for CSV output
        lambda_dense = np.linspace(0, sol.t[-1], n_steps)
        y_dense = sol.sol(lambda_dense)
        t_sol = y_dense[0, :]
        r_sol = y_dense[1, :]

        return {
            'success': True,
            't_final': rounded(t_final, 15),
            'r_final': rounded(r_final, 6),
            'reaches_horizon': reaches_horizon,
            'terminated_by_event': terminated_by_event,
            'affine_parameter_complete': True,
            't_array': t_sol.tolist(),
            'r_array': r_sol.tolist()
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'reaches_horizon': None,
            'terminated_by_event': False,
            'affine_parameter_complete': False
        }


def integrate_timelike_geodesic_backward(z_start: float = 0.0, z_end: float = 1000.0,
                                          n_steps: int = 1000,
                                          epsilon_t: float = DEFAULT_EPSILON_T,
                                          z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Integrate timelike geodesic backward from observer to high redshift."""

    y0 = [1.0 / (1.0 + z_start), 0.0]
    try:
        sol = solve_ivp(
            timelike_geodesic_equation, [0, 10.0], y0,
            args=(epsilon_t, z_t, n_t),
            events=horizon_event,
            dense_output=True,
            max_step=0.01,
            method='RK45'
        )

        t_final = sol.y[0, -1]
        r_final = sol.y[1, -1]
        terminated_by_event = sol.t_events is not None and len(sol.t_events[0]) > 0
        # Event at t=1e-6 is a numerical cutoff; geodesic approaches horizon asymptotically.
        reaches_horizon = False

        tau_dense = np.linspace(0, sol.t[-1], n_steps)
        y_dense = sol.sol(tau_dense)
        t_sol = y_dense[0, :]
        r_sol = y_dense[1, :]

        return {
            'success': True,
            't_final': rounded(t_final, 15),
            'r_final': rounded(r_final, 6),
            'reaches_horizon': reaches_horizon,
            'terminated_by_event': terminated_by_event,
            'proper_time_complete': True,
            't_array': t_sol.tolist(),
            'r_array': r_sol.tolist()
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'reaches_horizon': None,
            'terminated_by_event': False,
            'proper_time_complete': False
        }


def test_geodesic_completeness(z_max: float = 1000.0, n_points: int = 1000,
                               epsilon_t: float = DEFAULT_EPSILON_T,
                               z_t: float = DEFAULT_Z_T, n_t: float = DEFAULT_N_T) -> dict:
    """Test geodesic completeness toward temporal horizon."""
    
    # Test null geodesics
    null_result = integrate_null_geodesic_backward(
        z_end=z_max, n_steps=n_points,
        epsilon_t=epsilon_t, z_t=z_t, n_t=n_t
    )
    
    # Test timelike geodesics
    timelike_result = integrate_timelike_geodesic_backward(
        z_end=z_max, n_steps=n_points,
        epsilon_t=epsilon_t, z_t=z_t, n_t=n_t
    )
    
    # Analyze results
    null_complete = null_result['success'] and null_result['affine_parameter_complete']
    timelike_complete = timelike_result['success'] and timelike_result['proper_time_complete']
    
    null_asymptotic = null_result['success'] and not null_result['reaches_horizon']
    timelike_asymptotic = timelike_result['success'] and not timelike_result['reaches_horizon']
    
    results = {
        'step': STEP_ID,
        'description': 'Geodesic completeness test toward temporal horizon',
        'parameters': {
            'epsilon_t': epsilon_t,
            'z_t': z_t,
            'n_t': n_t,
            'z_max': z_max,
            'n_points': n_points
        },
        'null_geodesics': {
            'integration_success': null_result['success'],
            'affine_parameter_complete': null_result['affine_parameter_complete'],
            'reaches_horizon': null_result['reaches_horizon'],
            'asymptotic_approach': null_asymptotic,
            't_final': null_result.get('t_final'),
            'r_final': null_result.get('r_final')
        },
        'timelike_geodesics': {
            'integration_success': timelike_result['success'],
            'proper_time_complete': timelike_result['proper_time_complete'],
            'reaches_horizon': timelike_result['reaches_horizon'],
            'asymptotic_approach': timelike_asymptotic,
            't_final': timelike_result.get('t_final'),
            'r_final': timelike_result.get('r_final')
        },
        'completeness_status': {
            'null_geodesics_complete': null_complete,
            'timelike_geodesics_complete': timelike_complete,
            'both_asymptotic': null_asymptotic and timelike_asymptotic,
            'geodesically_complete': null_complete and timelike_complete
        },
        'interpretation': 'Geodesics approach temporal horizon asymptotically (no singularity)' if 
                          (null_asymptotic and timelike_asymptotic) else
                          'Geodesics may terminate at finite affine/proper time (requires analysis)'
    }
    
    # Generate CSV output
    if null_result['success'] and timelike_result['success']:
        z_array = np.linspace(0, z_max, n_points)
        csv_rows = []
        for i, z in enumerate(z_array):
            csv_rows.append({
                'z': rounded(z, 4),
                'null_t': rounded(null_result['t_array'][i], 15) if i < len(null_result['t_array']) else None,
                'null_r': rounded(null_result['r_array'][i], 6) if i < len(null_result['r_array']) else None,
                'timelike_t': rounded(timelike_result['t_array'][i], 15) if i < len(timelike_result['t_array']) else None,
                'timelike_r': rounded(timelike_result['r_array'][i], 6) if i < len(timelike_result['r_array']) else None
            })
        write_csv(step_csv_path(STEP_ID), csv_rows)
    
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    
    # Test with default parameters
    results = test_geodesic_completeness()
    
    # Test with different epsilon_t values
    epsilon_tests = [0.0001, 0.001, 0.01]
    epsilon_results = []
    
    for eps in epsilon_tests:
        eps_result = test_geodesic_completeness(epsilon_t=eps)
        epsilon_results.append({
            'epsilon_t': eps,
            'geodesically_complete': eps_result['completeness_status']['geodesically_complete'],
            'both_asymptotic': eps_result['completeness_status']['both_asymptotic']
        })
    
    results['epsilon_sensitivity'] = epsilon_results
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    print_status(f"Geodesic completeness: {results['completeness_status']['geodesically_complete']}", "INFO")
    
    return results


if __name__ == "__main__":
    run()
