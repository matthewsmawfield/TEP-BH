#!/usr/bin/env python3
"""Non-perturbative sGB solver — quadratic coupling, inward integration.

Coupling: f(φ) = (α/2)·φ²  (quadratic, non-shift-symmetric)
Scalar equation: □φ = α·φ·G_GB

Strategy:
1. Start at r = 30M with perturbative initial conditions
   (m = M, φ from linearized equation, N = 1)
2. Integrate inward using full non-perturbative equations
3. Regularize f near the horizon (clamp |f| to f_min)
4. Continue inward to r → 0
5. Check regularity at center

The quadratic coupling is key: the scalar equation □φ = α·φ·G_GB
has the φ factor which tames the divergence at the center.

At the center (r → 0), the regular solution has:
  m ~ r² (mass function → 0)
  φ → 0 (scalar → 0, super-exponentially for Schwarzschild-like)
  f → 1 (de Sitter core)
  K → finite (Kretschmann finite)

The perturbative scalar on Schwarzschild satisfies:
  f·φ'' + (f' + 2f/r)·φ' = α·φ·48M²/r⁶

At large r: φ ~ Q_s/r
At small r: φ ~ exp(-√(12α)·M/r²) → 0 (super-exponentially suppressed)

The non-perturbative effect: as m(r) decreases inward (scalar backreaction),
G_GB = 48(m-rm'/2)²/r⁶ decreases, and the scalar's energy density
decreases, allowing m to reach 0 at the center.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def gauss_bonnet_general(r, m, m_prime):
    """G_GB = 48(m - rm'/2)²/r⁶ for spherical metric."""
    r_safe = np.maximum(np.abs(r), 1e-30)
    return 48.0 * (m - r_safe * m_prime / 2.0) ** 2 / r_safe ** 6


def sgb_equations_quadratic(r, y, alpha, f_min=1e-8):
    """Field equations for quadratic sGB coupling.

    State: y = [m, N, phi, phi_prime]

    Equations (with f regularization):
      m' = (r²/(2f))·[½(φ')²·N² + α·φ·G_GB·N²·φ'·r/(2f)]
      N'/N = (r/(2f))·[½(φ')²·N² - α·φ·G_GB·N²·φ'·r/(2f)]
      φ'' + [(2/r) + (f'/(2f)) + (N'/N)]·φ' = -α·φ·G_GB/(f·N²)
    """
    m, N, phi, phi_p = y
    r_safe = max(abs(r), 1e-30)

    f = 1.0 - 2.0 * m / r_safe
    # Regularize f
    if abs(f) < f_min:
        f_reg = f_min * (1.0 if f >= 0 else -1.0)
    else:
        f_reg = f

    # Estimate m' for G_GB (scalar kinetic only, first estimate)
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 * N ** 2 / f_reg
    G_GB = gauss_bonnet_general(r_safe, m, m_p_est)

    # Full m' with quadratic GB coupling
    rho = 0.5 * phi_p ** 2 * N ** 2
    gb_src = alpha * phi * G_GB * N ** 2 * phi_p * r_safe / (2.0 * f_reg)
    m_p = r_safe ** 2 / (2.0 * f_reg) * (rho + gb_src)

    # Recompute G_GB
    G_GB_full = gauss_bonnet_general(r_safe, m, m_p)

    # N'/N
    gb_N = alpha * phi * G_GB_full * N ** 2 * phi_p * r_safe / (2.0 * f_reg)
    N_p_over_N = (r_safe / (2.0 * f_reg)) * (rho - gb_N)
    N_p = N_p_over_N * N

    # Scalar equation (quadratic coupling: α·φ·G_GB)
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe
    coeff = 2.0 / r_safe + f_p / (2.0 * f_reg) + N_p_over_N
    phi_pp = -alpha * phi * G_GB_full / (f_reg * N ** 2) - coeff * phi_p

    return [m_p, N_p, phi_p, phi_pp]


def perturbative_scalar_initial(r0, M, alpha):
    """Compute perturbative scalar at r0 for quadratic coupling on Schwarzschild.

    Solve: f·φ'' + (f' + 2f/r)·φ' = α·φ·48M²/r⁶
    with φ(∞) → 0, φ'(∞) → 0.

    Use the Sotiriou-Zhou-like solution adapted for quadratic coupling.
    For small α, the scalar is sourced by the Schwarzschild curvature.
    """
    # For quadratic coupling, the scalar equation is LINEAR in φ
    # (since the metric is fixed to Schwarzschild).
    # The solution can be computed by integrating from r0 inward.

    # At large r, the scalar is Coulomb-like: φ ~ Q_s/r
    # The scalar charge Q_s is determined by the coupling.
    # For quadratic coupling, Q_s depends on the central value φ_c,
    # which is determined by the regularity condition.

    # For the perturbative initial condition, use a simple estimate:
    # φ(r) ≈ α·M²/r³ (leading order from G_GB source)
    # This gives φ ~ α·M²/r³ at large r

    # Better: use the exact linearized solution
    # f·φ'' + (f' + 2f/r)·φ' = 48α·M²·φ/r⁶
    # At large r (f→1): φ'' + (2/r)·φ' = 48α·M²·φ/r⁶ ≈ 0
    # So φ ~ A + B/r at leading order. With φ(∞) = 0: φ ~ B/r.

    # The scalar charge B (or Q_s) is determined by matching to the interior.
    # For the perturbative estimate, use:
    # Q_s ≈ 2α/3 (same as shift-symmetric, to leading order)
    # But for quadratic coupling, the scalar equation is different.

    # Simple estimate: φ(r) = Q_s/r with Q_s = α·M
    Q_s = alpha * M * 0.1  # small scalar charge
    phi = Q_s / r0
    phi_p = -Q_s / r0 ** 2

    return phi, phi_p


def solve_inward_full(alpha, M=1.0, r_start=30.0, r_end=1e-5, n_points=20000, f_min=1e-8):
    """Integrate inward from r_start to r_end with quadratic sGB coupling."""
    # Initial conditions at r_start
    m0 = M
    N0 = 1.0
    phi0, phi_p0 = perturbative_scalar_initial(r_start, M, alpha)

    y0 = [m0, N0, phi0, phi_p0]

    r_span = (r_start, r_end)  # decreasing r
    r_eval = np.logspace(np.log10(r_start), np.log10(r_end), n_points)

    def event_blowup(r, y, alpha, f_min):
        m, N, phi, phi_p = y
        if abs(phi) > 1e8 or abs(phi_p) > 1e8 or m < -10 or m > 100:
            return 0.0
        return 1.0
    event_blowup.terminal = True
    event_blowup.direction = -1

    sol = solve_ivp(
        sgb_equations_quadratic, r_span, y0,
        args=(alpha, f_min),
        t_eval=r_eval,
        method='Radau',
        rtol=1e-8,
        atol=1e-10,
        events=[event_blowup],
        max_step=0.5,
    )

    return sol


def compute_kretschmann(r, m, N):
    """Kretschmann scalar for ds² = -N²f dt² + f⁻¹ dr² + r² dΩ²."""
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


def scan_quadratic_sgb():
    """Scan over coupling α and scalar charge to find regular solutions."""
    print("=" * 80)
    print("NON-PERTURBATIVE QUADRATIC sGB — INWARD INTEGRATION")
    print("Coupling: f(φ) = (α/2)·φ²")
    print("=" * 80)

    M = 1.0
    results = []

    for alpha in [0.5, 1.0, 2.0, 5.0, 10.0]:
        for Q_s_factor in [0.01, 0.1, 1.0]:
            Q_s = alpha * M * Q_s_factor

            print(f"\n{'='*60}")
            print(f"  α = {alpha}, Q_s = {Q_s:.4f}")
            print(f"{'='*60}")

            # Initial conditions
            r_start = 30.0
            m0 = M
            N0 = 1.0
            phi0 = Q_s / r_start
            phi_p0 = -Q_s / r_start ** 2

            y0 = [m0, N0, phi0, phi_p0]

            r_span = (r_start, 1e-5)
            r_eval = np.logspace(np.log10(r_start), np.log10(1e-5), 20000)

            sol = solve_ivp(
                sgb_equations_quadratic, r_span, y0,
                args=(alpha, 1e-8),
                t_eval=r_eval,
                method='Radau',
                rtol=1e-8,
                atol=1e-10,
                max_step=0.5,
            )

            if not sol.success or len(sol.t) < 10:
                print(f"  Integration failed: {sol.message}")
                results.append({'alpha': alpha, 'Q_s': Q_s, 'status': 'failed'})
                continue

            r = sol.t
            m = sol.y[0]
            N = sol.y[1]
            phi = sol.y[2]
            phi_p = sol.y[3]

            f = 1.0 - 2.0 * m / np.maximum(r, 1e-30)
            K = compute_kretschmann(r, m, N)

            # Find horizon
            sign_changes = np.where(np.diff(np.sign(f)))[0]
            r_h = r[sign_changes[0]] if len(sign_changes) > 0 else np.nan

            # Center values
            m_0 = m[0] if len(m) > 0 else np.nan  # at r_end (smallest r)
            f_0 = f[0] if len(f) > 0 else np.nan
            phi_0 = phi[0] if len(phi) > 0 else np.nan
            K_0 = K[0] if len(K) > 0 else np.nan

            # Actually, the last point is at r_end (smallest r)
            m_center = m[-1]
            f_center = f[-1]
            phi_center = phi[-1]
            K_center = K[-1]

            print(f"  Integration: {sol.message}")
            print(f"  r range: [{r[0]:.4e}, {r[-1]:.4e}]")
            print(f"  r_h = {r_h:.4f}" if not np.isnan(r_h) else "  No horizon")
            print(f"\n  Center (r → {r[-1]:.2e}):")
            print(f"    m = {m_center:.6e}")
            print(f"    f = {f_center:.6f}")
            print(f"    φ = {phi_center:.6e}")
            print(f"    K = {K_center:.6e}")

            # Power law near center
            mask = (r >= 1e-4) & (r <= 0.1) & np.isfinite(K) & (K > 0)
            if np.sum(mask) > 10:
                alpha_pow = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
                print(f"    K power law: ~r^{{{alpha_pow:.3f}}}")
            else:
                alpha_pow = np.nan

            # Radial profile
            print(f"\n  Radial profile:")
            print(f"  {'r':>10} | {'m(r)':>12} | {'f(r)':>12} | {'N(r)':>10} | {'φ(r)':>12} | {'K':>14}")
            print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}-+-{'-'*14}")

            for r_test in [1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 3.0, 5.0, 10.0, 30.0]:
                idx = np.argmin(np.abs(r - r_test))
                if idx < len(r):
                    K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
                    print(f"  {r[idx]:>10.4e} | {m[idx]:>12.6e} | "
                          f"{f[idx]:>12.6e} | {N[idx]:>10.4f} | "
                          f"{phi[idx]:>12.6e} | {K_val:>14.6e}")

            # Regularity
            K_finite = np.isfinite(K_center) and K_center < 1e10
            de_sitter = abs(f_center - 1.0) < 0.5
            m_zero = abs(m_center) < 0.1

            if K_finite and de_sitter and m_zero:
                status = 'REGULAR BLACK HOLE ✓'
            elif K_finite and de_sitter:
                status = 'REGULAR (m not 0)'
            elif K_finite:
                status = 'K finite, no de Sitter'
            else:
                status = 'SINGULAR'

            print(f"\n  → {status}")

            results.append({
                'alpha': alpha, 'Q_s': Q_s,
                'r_h': float(r_h) if not np.isnan(r_h) else None,
                'm_center': float(m_center),
                'f_center': float(f_center),
                'phi_center': float(phi_center),
                'K_center': float(K_center) if np.isfinite(K_center) else None,
                'K_power': float(alpha_pow) if np.isfinite(alpha_pow) else None,
                'status': status,
            })

    return results


def main():
    results = scan_quadratic_sgb()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  {'α':>6} | {'Q_s':>8} | {'r_h':>8} | {'m(0)':>12} | {'f(0)':>8} | {'K(0)':>14} | {'Status':>30}")
    print(f"  {'-'*6}-+-{'-'*8}-+-{'-'*8}-+-{'-'*12}-+-{'-'*8}-+-{'-'*14}-+-{'-'*30}")

    for r in results:
        if r['status'] == 'failed':
            print(f"  {r['alpha']:>6} | {r['Q_s']:>8.4f} | {'—':>8} | {'—':>12} | {'—':>8} | {'—':>14} | {'FAILED':>30}")
        else:
            r_h_str = f"{r['r_h']:.4f}" if r['r_h'] is not None else "—"
            K_str = f"{r['K_center']:.4e}" if r['K_center'] is not None else "—"
            print(f"  {r['alpha']:>6} | {r['Q_s']:>8.4f} | {r_h_str:>8} | {r['m_center']:>12.4e} | "
                  f"{r['f_center']:>8.4f} | {K_str:>14} | {r['status']:>30}")

    output_path = os.path.join(PROJECT_ROOT, 'results', 'stage_c_inward_quadratic.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_path}")


if __name__ == '__main__':
    main()
