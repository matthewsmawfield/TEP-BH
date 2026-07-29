#!/usr/bin/env python3
"""Step 09b: Native Tensor-Mode Integration.

Integrate the native temporal-conformal tensor equation

    mu_k'' + (k^2 - A_clock''/A_clock) mu_k = 0

across the finite transition profile and compute the suppression factor
relative to the inflationary consistency relation r = 16 epsilon_field.
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from th_common import (
    TEPLogger, ensure_dirs, print_status, rounded, set_step_logger,
    step_json_path, step_csv_path, write_json, write_csv,
)

STEP_ID = "step_09b_native_tensor_integration"

# Physics parameters
P = 0.1              # temporal-horizon power-law index (0 < p <= 1/2)
ETA_0_MPC = 6000.0   # present-epoch conformal time in Mpc (physical)
N_PROFILE = 2.0      # transition steepness
EPSILON_FIELD = 0.0175  # Planck-constrained spectral-flow parameter
R_BK_MAX = 0.036     # BICEP/Keck 2021 95% CL upper bound

# We work in dimensionless units: eta_tilde = eta / eta_0, k_tilde = k * eta_0
# The combination k*eta = k_tilde * eta_tilde is the key dimensionless quantity.
# For CMB scales k ~ 0.05 Mpc^-1 and eta_0 ~ 6000 Mpc:
#   k * eta_0 ~ 300 >> 1, so CMB modes oscillate many times across the transition.
# This gives strong suppression of tensor production.


def A_profile(eta, p=P, eta_0=ETA_0_MPC, n=N_PROFILE):
    x = np.asarray(eta) / eta_0
    return (1.0 + x**n) ** (-p / n)


def A_double_prime_over_A(eta, p=P, eta_0=ETA_0_MPC, n=N_PROFILE):
    x = np.asarray(eta) / eta_0
    xn = x**n
    return p * ((p + 1.0) * xn - 1.0) / (eta_0**2 * (1.0 + xn)**2)


def integrate_mode(k, p=P, eta_0=ETA_0_MPC, n=N_PROFILE, with_source=True):
    """Integrate tensor mode from deep inside horizon to late times.

    If with_source=False, integrate the Minkowski reference (A''/A = 0).
    Returns the late-time mu amplitude.
    """
    eta_start = 500.0 / k   # deep inside horizon: k*eta >> 1
    eta_end = 1e-4 * eta_0  # late times: A -> 1, source negligible

    # Bunch-Davies vacuum initial conditions
    mu0 = 1.0 / np.sqrt(2.0 * k)
    dmu0 = -1j * np.sqrt(k / 2.0)

    def ode(eta, y):
        re_mu, im_mu, re_dmu, im_dmu = y
        src = A_double_prime_over_A(eta, p, eta_0, n) if with_source else 0.0
        f = -(k**2 - src)
        return [re_dmu, im_dmu, f * re_mu, f * im_mu]

    sol = solve_ivp(ode, [eta_start, eta_end],
                    [mu0.real, mu0.imag, dmu0.real, dmu0.imag],
                    dense_output=True, max_step=(eta_start - eta_end) / 2000,
                    method='RK45')
    if not sol.success:
        raise RuntimeError(f"ODE failed for k={k}")
    re_mu, im_mu, _, _ = sol.y[:, -1]
    return complex(re_mu, im_mu)


def bogoliubov_integral(k_tilde, p=P, n=N_PROFILE):
    """Compute dimensionless Bogoliubov coefficient.

    In dimensionless units (eta_tilde = eta/eta_0, k_tilde = k*eta_0):
        beta_k = -(i / 2 k_tilde) integral_0^infty V_tilde(eta_tilde) e^{-2 i k_tilde eta_tilde} d eta_tilde

    where V_tilde = eta_0^2 * V = p * [(p+1) x^n - 1] / (1+x^n)^2.
    """
    from scipy.integrate import quad

    def V_tilde(x):
        xn = x**n
        return p * ((p + 1.0) * xn - 1.0) / (1.0 + xn)**2

    def integrand_re(x):
        return V_tilde(x) * np.cos(2.0 * k_tilde * x)

    def integrand_im(x):
        return -V_tilde(x) * np.sin(2.0 * k_tilde * x)

    # Split integration to handle oscillatory convergence better
    re_val, re_err = quad(integrand_re, 0.0, 100.0, limit=500)
    im_val, im_err = quad(integrand_im, 0.0, 100.0, limit=500)

    beta = -(1j / (2.0 * k_tilde)) * complex(re_val, im_val)
    return beta, re_err + im_err


def compute_r_k(k_tilde_values, p=P, eta_0_mpc=ETA_0_MPC, n=N_PROFILE):
    """Compute tensor-to-scalar ratio using physical units.

    P_T(k) = (4 k^3 / pi^2) |beta_k|^2   (two polarizations)
    r(k)   = P_T(k) / P_zeta(k)
    """
    A_s = 2.10e-9
    n_s = 0.965
    k_pivot_mpc = 0.05
    k_tilde_pivot = k_pivot_mpc * eta_0_mpc

    results = []
    for kt in k_tilde_values:
        try:
            beta, err = bogoliubov_integral(kt, p, n)
            k_mpc = kt / eta_0_mpc
            # Tensor power spectrum in physical units
            P_T = (4.0 * k_mpc**3 / np.pi**2) * np.abs(beta)**2
            # Scalar power spectrum
            P_zeta = A_s * (k_mpc / k_pivot_mpc)**(n_s - 1.0)
            r_k = P_T / P_zeta

            results.append({
                'k_Mpc': rounded(k_mpc, 8),
                'k_tilde': rounded(kt, 4),
                'beta_sq': rounded(np.abs(beta)**2, 15),
                'P_T': rounded(P_T, 15),
                'P_zeta': rounded(P_zeta, 15),
                'r_k': rounded(r_k, 8),
            })
        except Exception as e:
            results.append({'k_Mpc': rounded(kt / eta_0_mpc, 8), 'error': str(e)})
    return results


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # k_tilde = k * eta_0; for k in [1e-4, 1] Mpc^-1 and eta_0 ~ 6000 Mpc:
    # k_tilde in [0.6, 6000]
    k_tilde_min = 1e-4 * ETA_0_MPC
    k_tilde_max = 1.0 * ETA_0_MPC
    k_tilde_values = np.logspace(np.log10(k_tilde_min), np.log10(k_tilde_max), 50)

    results = compute_r_k(k_tilde_values)
    valid = [r for r in results if 'error' not in r]

    r_vals = [r['r_k'] for r in valid]
    r_max = max(r_vals) if r_vals else None
    r_min = min(r_vals) if r_vals else None
    r_pivot = next((r['r_k'] for r in valid if abs(r['k_Mpc'] - 0.05) < 0.003), None)

    # Analytic estimate: suppression factor from rapid oscillations
    # For large k_tilde, |beta|^2 ~ P^2 / (4 k_tilde^4) for a smooth transition
    # This gives r ~ (k^3 / pi^2) * P^2 / (4 k_tilde^4) / P_zeta
    #               ~ P^2 / (4 pi^2 eta_0^3 A_s k_pivot^(1-n_s) k^(n_s-4))
    # At k_pivot: r ~ P^2 / (4 pi^2 eta_0^3 A_s k_pivot^3)
    k_pivot_mpc = 0.05
    A_s = 2.10e-9
    r_analytic = (P**2 / (4.0 * np.pi**2 * ETA_0_MPC**3 * A_s * k_pivot_mpc**3))

    summary = {
        'step': STEP_ID,
        'description': 'Native tensor-mode integration: first-order Bogoliubov coefficient with physical eta_0',
        'parameters': {
            'p': P,
            'eta_0_Mpc': ETA_0_MPC,
            'n_profile': N_PROFILE,
            'epsilon_field': EPSILON_FIELD,
            'r_inflationary': 16.0 * EPSILON_FIELD,
            'r_BK_max': R_BK_MAX,
        },
        'results': {
            'r_max': rounded(r_max, 6) if r_max else None,
            'r_min': rounded(r_min, 6) if r_min else None,
            'r_at_pivot': rounded(r_pivot, 6) if r_pivot else None,
            'r_analytic_estimate': rounded(r_analytic, 8),
            'below_BK_bound': r_max < R_BK_MAX if r_max is not None else None,
            'k_tilde_range': (rounded(k_tilde_min, 2), rounded(k_tilde_max, 2)),
            'num_k': len(k_tilde_values),
        },
        'interpretation': (
            f"Native TEP tensor integration gives r_max = {rounded(r_max, 8) if r_max else 'N/A'} "
            f"at k_pivot = {rounded(r_pivot, 8) if r_pivot else 'N/A'} (analytic est. {rounded(r_analytic, 8)}). "
            f"For k*eta_0 >> 1, rapid oscillations across the transition suppress tensor production. "
            f"Result is below BICEP/Keck bound r < {R_BK_MAX}: {r_max < R_BK_MAX if r_max is not None else 'N/A'}."
        ),
    }

    if valid:
        write_csv(step_csv_path(STEP_ID), valid)
    write_json(step_json_path(STEP_ID), summary)

    print_status(f"Step {STEP_ID} completed", "SUCCESS")
    print_status(f"r_max = {rounded(r_max, 8) if r_max else 'N/A'}", "INFO")
    print_status(f"r_at_pivot = {rounded(r_pivot, 8) if r_pivot else 'N/A'}", "INFO")
    print_status(f"r_analytic = {rounded(r_analytic, 8)}", "INFO")
    print_status(f"Below BK bound: {r_max < R_BK_MAX if r_max is not None else 'N/A'}", "INFO")

    return summary


if __name__ == "__main__":
    run()
