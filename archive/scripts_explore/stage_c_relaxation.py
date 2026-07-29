#!/usr/bin/env python3
"""Non-perturbative sGB solver — EF-regular formulation with horizon series.

Key insight: The scalar equation in EF-regular form has NO 1/f singularity:
  f·φ'' + (2f/r + f')·φ' + f·(N'/N)·φ' = -α·φ·G_GB/N²

At the horizon (f=0): f'·φ' = -α·φ·G_GB/N²  (FINITE!)

Strategy:
1. Start at r = r_h + ε with series expansion initial conditions
2. The horizon regularity condition gives φ'(r_h) from the scalar equation
3. Integrate outward (to r=∞) and inward (to r=0) from the horizon
4. The mass equation uses a regularized form near the horizon

For the mass and N equations, we use the fact that at the horizon,
the GB term cancels the scalar kinetic term (regularity condition),
making the 0/0 limit finite.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def gauss_bonnet_general(r, m, m_prime):
    """G_GB = 48(m - rm'/2)²/r⁶."""
    r_safe = np.maximum(np.abs(r), 1e-30)
    return 48.0 * (m - r_safe * m_prime / 2.0) ** 2 / r_safe ** 6


def sgb_equations_ef_regular(r, y, alpha):
    """EF-regular field equations for quadratic sGB coupling.

    State: y = [m, ln_N, phi, phi_prime]
    (Using ln_N instead of N for better numerics)

    The scalar equation is written in EF-regular form:
      f·φ'' + (2f/r + f')·φ' + f·(N'/N)·φ' = -α·φ·G_GB/N²

    The mass equation uses the regularized form:
      m' = r²·½(φ')²·N² / (2f)  [scalar kinetic only, regularized]

    The N equation uses:
      f·(N'/N) = r·½(φ')²·N²/2 - r·α·φ·G_GB·N²·φ'·r/(4f)
      [multiplied by f to remove 1/f from first term]
    """
    m, log_N, phi, phi_p = y
    r_safe = max(abs(r), 1e-30)

    N = np.exp(log_N)
    f = 1.0 - 2.0 * m / r_safe

    # Regularize f
    f_abs = abs(f)
    if f_abs < 1e-10:
        f_reg = 1e-10 * (1.0 if f >= 0 else -1.0)
    else:
        f_reg = f

    # First estimate m' (scalar kinetic only)
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 * N ** 2 / (2.0 * f_reg)
    G_GB = gauss_bonnet_general(r_safe, m, m_p_est)

    # Mass equation: m' = r²/(2f)·[½(φ')²·N² + α·φ·G_GB·N²·φ'·r/(2f)]
    # Near horizon, use regularized form: just scalar kinetic
    # (GB term is higher order and causes 1/f² singularity)
    rho_kin = 0.5 * phi_p ** 2 * N ** 2
    gb_term = alpha * phi * G_GB * N ** 2 * phi_p * r_safe / (2.0 * f_reg)

    # Blend: use full equation when |f| > 0.1, scalar-only when |f| < 0.01
    blend = min(1.0, max(0.0, (f_abs - 0.01) / 0.09))
    m_p = r_safe ** 2 / (2.0 * f_reg) * (rho_kin + blend * gb_term)

    # Recompute G_GB with full m'
    G_GB_full = gauss_bonnet_general(r_safe, m, m_p)

    # N equation (EF-regular form):
    # f·(N'/N) = r·½(φ')²·N²/2 - r·α·φ·G_GB·N²·φ'·r/(4f)
    # = r·N²·[½(φ')²/2 - α·φ·G_GB·φ'·r/(4f)]
    gb_N_term = alpha * phi * G_GB_full * phi_p * r_safe / (4.0 * f_reg)
    f_N_over_N = r_safe * N ** 2 * (rho_kin / 2.0 - blend * gb_N_term)
    N_over_N = f_N_over_N / f_reg  # This has 1/f but the numerator → 0 at horizon
    log_N_p = N_over_N

    # Scalar equation (EF-regular form):
    # f·φ'' + (2f/r + f')·φ' + f·(N'/N)·φ' = -α·φ·G_GB/N²
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe
    phi_pp = (-alpha * phi * G_GB_full / N ** 2
              - (2.0 * f_reg / r_safe + f_p) * phi_p
              - f_reg * N_over_N * phi_p) / f_reg

    return [m_p, log_N_p, phi_p, phi_pp]


def solve_from_horizon(r_h, alpha, m1, phi_h, N_h=1.0, direction='outward',
                       r_end=50.0, n_points=10000):
    """Solve from r_h + ε (or r_h - ε) to r_end.

    Parameters:
        r_h: horizon radius
        alpha: GB coupling
        m1: dm/dr at horizon (shooting parameter)
        phi_h: scalar at horizon
        N_h: lapse at horizon (gauge)
        direction: 'outward' (r > r_h) or 'inward' (r < r_h)
    """
    delta = 1e-3 * r_h

    if direction == 'outward':
        r_start = r_h + delta
        r_span = (r_start, r_end)
    else:
        r_start = r_h - delta
        r_span = (r_start, max(r_end, 1e-6))  # r_end should be small

    # Series expansion at r = r_h + δ
    M = r_h / 2.0
    G_GB_h = 48.0 * M ** 2 / r_h ** 6  # = 12/r_h⁴

    # f'(r_h) = (1 - 2m₁)/r_h
    f_p_h = (1.0 - 2.0 * m1) / r_h

    # Horizon regularity: φ'(r_h) = -α·φ_h·G_GB_h / (f'(r_h)·N_h²)
    phi_1 = -alpha * phi_h * G_GB_h / (f_p_h * N_h ** 2)

    # Initial conditions
    sign = 1.0 if direction == 'outward' else -1.0
    m_start = r_h / 2.0 + m1 * sign * delta
    phi_start = phi_h + phi_1 * sign * delta
    phi_p_start = phi_1
    log_N_start = np.log(N_h)

    y0 = [m_start, log_N_start, phi_start, phi_p_start]

    r_eval = np.logspace(np.log10(r_start), np.log10(r_span[1]), n_points)
    if direction == 'inward':
        # t_eval must be sorted in the direction of integration (decreasing)
        r_eval = r_eval[::-1]  # decreasing
        # Ensure sorted (descending)
        r_eval = np.sort(r_eval)[::-1]

    sol = solve_ivp(
        sgb_equations_ef_regular, r_span, y0,
        args=(alpha,),
        t_eval=r_eval,
        method='Radau',
        rtol=1e-8,
        atol=1e-10,
        max_step=0.5,
    )

    return sol


def compute_kretschmann(r, m, N):
    """Kretschmann scalar."""
    r_safe = np.maximum(np.abs(r), 1e-30)
    N_safe = np.where(np.abs(N) > 1e-30, N, np.sign(N) * 1e-30)

    f = 1.0 - 2.0 * m / r_safe
    f_safe = np.where(np.abs(f) > 1e-30, f, np.sign(f) * 1e-30)

    m_p = np.gradient(m, r)
    N_p = np.gradient(N, r)
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe
    f_pp = np.gradient(f_p, r)

    N_over_N = N_p / N_safe

    R_trtr = N_safe ** 2 * (f_pp / 2.0 + f_p * N_over_N + f_safe * N_over_N ** 2)
    R_tthtth = N_safe ** 2 * (f_safe * f_p / (2.0 * r_safe) + f_safe ** 2 * N_over_N / r_safe)
    R_rthrth = f_p / (2.0 * r_safe)
    R_thphthph = (1.0 - f_safe) / r_safe ** 2

    K = 4.0 * R_trtr ** 2 + 8.0 * R_tthtth ** 2 + 8.0 * R_rthrth ** 2 + 4.0 * R_thphthph ** 2

    return K


def scan_horizon_solutions():
    """Scan over parameters to find regular black hole solutions."""
    print("=" * 80)
    print("NON-PERTURBATIVE QUADRATIC sGB — HORIZON SERIES + EF-REGULAR")
    print("Coupling: f(φ) = (α/2)·φ²")
    print("=" * 80)

    M_target = 1.0
    r_h = 2.0 * M_target

    results = []

    for alpha in [0.5, 1.0, 2.0, 5.0]:
        for phi_h in [0.01, 0.1, 0.5, 1.0]:
            for m1 in [-0.1, 0.0, 0.1]:
                print(f"\n  α={alpha}, φ_h={phi_h}, m₁={m1}")

                # Solve outward
                sol_out = solve_from_horizon(r_h, alpha, m1, phi_h,
                                             direction='outward', r_end=30.0, n_points=5000)

                # Solve inward
                sol_in = solve_from_horizon(r_h, alpha, m1, phi_h,
                                            direction='inward', r_end=1e-5, n_points=5000)

                out_ok = sol_out.success and len(sol_out.t) > 10
                in_ok = sol_in.success and len(sol_in.t) > 10

                if not out_ok and not in_ok:
                    print(f"    Both failed: out={sol_out.message}, in={sol_in.message}")
                    continue

                # Combine
                if in_ok and out_ok:
                    r_in = sol_in.t[::-1]
                    m_in = sol_in.y[0][::-1]
                    logN_in = sol_in.y[1][::-1]
                    phi_in = sol_in.y[2][::-1]
                    phi_p_in = sol_in.y[3][::-1]

                    r_out = sol_out.t
                    m_out = sol_out.y[0]
                    logN_out = sol_out.y[1]
                    phi_out = sol_out.y[2]
                    phi_p_out = sol_out.y[3]

                    r = np.concatenate([r_in, r_out])
                    m = np.concatenate([m_in, m_out])
                    logN = np.concatenate([logN_in, logN_out])
                    phi = np.concatenate([phi_in, phi_out])
                    phi_p = np.concatenate([phi_p_in, phi_p_out])

                    sort_idx = np.argsort(r)
                    r = r[sort_idx]
                    m = m[sort_idx]
                    logN = logN[sort_idx]
                    phi = phi[sort_idx]
                    phi_p = phi_p[sort_idx]

                    N = np.exp(logN)
                    f = 1.0 - 2.0 * m / np.maximum(r, 1e-30)
                    K = compute_kretschmann(r, m, N)

                    M_adm = m[-1]
                    m_center = m[0]
                    f_center = f[0]
                    phi_center = phi[0]
                    K_center = K[0]

                    # Power law
                    mask = (r >= 1e-4) & (r <= 0.1) & np.isfinite(K) & (K > 0)
                    if np.sum(mask) > 10:
                        alpha_pow = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
                    else:
                        alpha_pow = np.nan

                    K_finite = np.isfinite(K_center) and K_center < 1e10
                    de_sitter = abs(f_center - 1.0) < 0.5
                    m_zero = abs(m_center) < 0.1

                    if K_finite and de_sitter and m_zero:
                        status = 'REGULAR BH ✓'
                    elif K_finite and de_sitter:
                        status = 'REGULAR (m≠0)'
                    elif K_finite:
                        status = 'K finite'
                    else:
                        status = 'SINGULAR'

                    print(f"    M_ADM={M_adm:.4f}, m(0)={m_center:.4e}, "
                          f"f(0)={f_center:.4f}, K(0)={K_center:.4e}, {status}")

                    if 'REGULAR BH' in status or 'REGULAR' in status:
                        print(f"    K power law: ~r^{{{alpha_pow:.3f}}}")
                        print(f"    Radial profile:")
                        for r_test in [1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 10.0, 30.0]:
                            idx = np.argmin(np.abs(r - r_test))
                            if idx < len(r):
                                K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
                                print(f"      r={r[idx]:.4e} m={m[idx]:.4e} f={f[idx]:.4f} "
                                      f"φ={phi[idx]:.4e} K={K_val:.4e}")

                    results.append({
                        'alpha': alpha, 'phi_h': phi_h, 'm1': m1,
                        'M_adm': float(M_adm),
                        'm_center': float(m_center),
                        'f_center': float(f_center),
                        'phi_center': float(phi_center),
                        'K_center': float(K_center) if np.isfinite(K_center) else None,
                        'K_power': float(alpha_pow) if np.isfinite(alpha_pow) else None,
                        'status': status,
                    })

    return results


def main():
    results = scan_horizon_solutions()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  {'α':>5} | {'φ_h':>6} | {'m₁':>6} | {'M_ADM':>8} | {'m(0)':>12} | "
          f"{'f(0)':>8} | {'K(0)':>14} | {'Status':>20}")
    print(f"  {'-'*5}-+-{'-'*6}-+-{'-'*6}-+-{'-'*8}-+-{'-'*12}-+-{'-'*8}-+-{'-'*14}-+-{'-'*20}")

    for r in results:
        K_str = f"{r['K_center']:.4e}" if r['K_center'] is not None else "—"
        print(f"  {r['alpha']:>5} | {r['phi_h']:>6} | {r['m1']:>6} | {r['M_adm']:>8.4f} | "
              f"{r['m_center']:>12.4e} | {r['f_center']:>8.4f} | {K_str:>14} | {r['status']:>20}")

    output_path = os.path.join(PROJECT_ROOT, 'results', 'stage_c_horizon_series.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_path}")


if __name__ == '__main__':
    main()
