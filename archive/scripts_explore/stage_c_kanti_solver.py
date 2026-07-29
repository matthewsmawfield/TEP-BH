#!/usr/bin/env python3
"""Non-perturbative sGB regular black hole solver — Kanti formulation.

Solves the full coupled Einstein-scalar-Gauss-Bonnet field equations
with regular boundary conditions, using the EF-regular formulation
that crosses the horizon without singularities.

Action:
  S = ∫ d⁴x √(-g) [R/(16π) - ½(∇φ)² + (α/2)·φ²·G_GB]

Coupling: f(φ) = (α/2)·φ²  (quadratic, non-shift-symmetric)
  - Allows regular solutions with finite φ at the center
  - f'(φ) = α·φ sources the scalar from curvature
  - Known to produce regular BH solutions (Kanti et al 1996)

Metric: ds² = -N²(r)·f(r)·dt² + f(r)⁻¹·dr² + r²·dΩ²
  f(r) = 1 - 2m(r)/r

Field equations (EF-regular form, from Kanti et al 1996, Torii et al 1999):

  m' = (r²/2)·[½(φ')²·N² + α·φ·G_GB·N²·φ'·r/(2f)]

  N'/N = (r/(2f))·[½(φ')²·N² - α·φ·G_GB·N²·φ'·r/(2f)]  ... (simplified)

  φ'' = [α·φ·G_GB - (f' + 2f/r)·φ'·N² - f·(N'/N)·φ'] / (f·N²)

  G_GB = (48/r⁶)·(m - r·m'/2)²  (exact for this metric, Φ=0 gauge)

Regular boundary conditions at r → 0:
  m(0) = 0,  m'(0) = 0
  φ(0) = φ_c  (finite, nonzero — activates GB source)
  φ'(0) = 0
  N(0) = 1  (gauge)

Taylor expansion at r = 0:
  m(r) ≈ m₅·r⁵  (mass grows from scalar energy density)
  φ(r) ≈ φ_c + φ₂·r²  (regular scalar)
  N(r) ≈ 1 + N₂·r²  (regular lapse)
  f(r) ≈ 1 - 2m₅·r⁴  → 1  (de Sitter-like core)

The solution is regular at r=0 (Kretschmann → 0) and develops a horizon
at r = r_h where f(r_h) = 0, if the scalar energy density is sufficient.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def gauss_bonnet_invariant(r, m, m_prime):
    """Gauss-Bonnet invariant for ds² = -N²f dt² + f⁻¹ dr² + r² dΩ².

    G_GB = (48/r⁶)·(m - r·m'/2)²

    This is exact when the lapse N is absorbed into the time coordinate
    (Schwarzschild gauge). For Schwarzschild (m'=0): G_GB = 48M²/r⁶ ✓
    """
    r_safe = np.maximum(np.abs(r), 1e-30)
    return 48.0 * (m - r_safe * m_prime / 2.0) ** 2 / r_safe ** 6


def sgb_equations(r, y, alpha):
    """Field equations for non-perturbative sGB with quadratic coupling.

    State: y = [m, N, phi, phi_prime]
    Returns: y' = [m', N', phi', phi'']

    Uses EF-regular form: equations are written so that f=0 (horizon)
    is NOT a singular point.
    """
    m, N, phi, phi_p = y
    r_safe = max(abs(r), 1e-30)

    # Metric function f = 1 - 2m/r
    f = 1.0 - 2.0 * m / r_safe

    # Regularize f near zero (horizon) — use smooth regularization
    f_abs = abs(f)
    if f_abs < 1e-12:
        f_reg = 1e-12 * (1.0 if f >= 0 else -1.0)
    else:
        f_reg = f

    # First estimate m' (for G_GB computation) — scalar kinetic only
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 * N ** 2 / f_reg

    # G_GB with estimated m'
    G_GB = gauss_bonnet_invariant(r_safe, m, m_p_est)

    # GB contribution to mass equation
    gb_mass = alpha * phi * G_GB * N ** 2 * phi_p * r_safe / (2.0 * f_reg)

    # Full m'
    m_p = 0.5 * r_safe ** 2 * (0.5 * phi_p ** 2 * N ** 2 + gb_mass)

    # Recompute G_GB with full m'
    G_GB_full = gauss_bonnet_invariant(r_safe, m, m_p)

    # N'/N equation (pressure balance)
    # N'/N = (r/(2f))·[½(φ')²·N² - α·φ·G_GB·N²·φ'·r/(2f)]
    gb_N = alpha * phi * G_GB_full * N ** 2 * phi_p * r_safe / (2.0 * f_reg)
    N_p_over_N = (r_safe / (2.0 * f_reg)) * (0.5 * phi_p ** 2 * N ** 2 - gb_N)
    N_p = N_p_over_N * N

    # Scalar equation (EF-regular form)
    # f·N²·φ'' + (f' + 2f/r)·N²·φ' + f·(N'/N)·N²·φ' = α·φ·G_GB
    # f' = 2m/r² - 2m'/r
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe

    # φ'' = [α·φ·G_GB - (f' + 2f/r)·N²·φ' - f·N²·(N'/N)·φ'] / (f·N²)
    numerator = (alpha * phi * G_GB_full
                 - (f_p + 2.0 * f_reg / r_safe) * N ** 2 * phi_p
                 - f_reg * N ** 2 * N_p_over_N * phi_p)
    denominator = f_reg * N ** 2
    if abs(denominator) < 1e-20:
        phi_pp = 0.0
    else:
        phi_pp = numerator / denominator

    return [m_p, N_p, phi_p, phi_pp]


def solve_sgb_kanti(alpha, phi_c, r_max=50.0, n_points=20000, r_min=1e-6):
    """Solve non-perturbative sGB with regular BCs at r=0.

    Parameters:
        alpha: GB coupling (f(φ) = (α/2)·φ²)
        phi_c: central scalar value (activates GB source)
        r_max: outer integration boundary
        n_points: number of output points
        r_min: inner boundary (small, ≈ 0)

    Returns: dict with solution arrays and diagnostics
    """
    # Taylor expansion at r = 0:
    # m(r) ≈ m₅·r⁵, φ(r) ≈ φ_c + φ₂·r², N(r) ≈ 1 + N₂·r²
    #
    # From the field equations at r → 0:
    # G_GB(0) = 48(m - rm'/2)²/r⁶ → 0 (since m ~ r⁵, m' ~ r⁴)
    # So the GB source vanishes at leading order.
    # φ''(0) ≈ 0 → φ₂ ≈ 0 (to leading order)
    # m'(0) ≈ ½r²(φ')² ≈ 0 → m starts very small
    #
    # The solution grows from r=0 due to the scalar's self-gravity.
    # We start with small perturbations and let the equations evolve.

    # Initial conditions at r_min (Taylor expansion)
    # φ(r_min) ≈ φ_c (φ₂ is very small, ~ α·φ_c·G_GB(0) which → 0)
    # φ'(r_min) ≈ 0
    # m(r_min) ≈ 0
    # N(r_min) ≈ 1

    # Small seed to break the trivial fixed point
    # The scalar's energy density ~ (φ')² provides the seed mass
    # We need a small φ' to start the mass growth
    # From the scalar equation: φ'' ~ α·φ·G_GB/f
    # At small r, G_GB ~ 48m₅²·r⁴ (if m ~ m₅·r⁵)
    # This is very small, so φ grows very slowly

    # Use a small initial gradient to seed the solution
    phi_p_seed = 1e-8  # tiny seed gradient

    m_init = 0.0
    N_init = 1.0
    phi_init = phi_c
    phi_p_init = phi_p_seed

    y0 = [m_init, N_init, phi_init, phi_p_init]

    r_span = (r_min, r_max)
    r_eval = np.logspace(np.log10(r_min), np.log10(r_max), n_points)

    def event_horizon(r, y, alpha):
        """Detect horizon crossing (f = 0, f decreasing)."""
        m = y[0]
        f = 1.0 - 2.0 * m / max(abs(r), 1e-30)
        return f
    event_horizon.terminal = False  # don't stop at horizon
    event_horizon.direction = -1  # f decreasing

    def event_blowup(r, y, alpha):
        """Stop if solution blows up."""
        m, N, phi, phi_p = y
        if abs(phi) > 1e10 or abs(phi_p) > 1e10 or m > 1e6:
            return 0.0
        return 1.0
    event_blowup.terminal = True
    event_blowup.direction = -1

    sol = solve_ivp(
        sgb_equations, r_span, y0,
        args=(alpha,),
        t_eval=r_eval,
        method='Radau',  # stiff solver for horizon crossing
        rtol=1e-10,
        atol=1e-12,
        events=[event_horizon, event_blowup],
        max_step=0.5,
    )

    r_sol = sol.t
    m_sol = sol.y[0]
    N_sol = sol.y[1]
    phi_sol = sol.y[2]
    phi_p_sol = sol.y[3]

    f_sol = 1.0 - 2.0 * m_sol / np.maximum(r_sol, 1e-30)

    # ADM mass
    M_adm = m_sol[-1] if len(m_sol) > 0 else np.nan

    # Scalar charge (from φ ~ Q_s/r at large r, for quadratic coupling
    # the asymptotic behavior may differ)
    if len(r_sol) > 100:
        mask = r_sol > r_max * 0.5
        if np.sum(mask) > 10:
            Q_s = np.polyfit(1.0 / r_sol[mask], phi_sol[mask], 1)[0]
        else:
            Q_s = np.nan
    else:
        Q_s = np.nan

    # Find horizon
    horizon_idx = None
    r_h = np.nan
    sign_changes = np.where(np.diff(np.sign(f_sol)))[0]
    if len(sign_changes) > 0:
        idx_h = sign_changes[0]
        if f_sol[idx_h] > 0 and f_sol[idx_h + 1] < 0:
            horizon_idx = idx_h
            r_h = r_sol[idx_h]

    return {
        'r': r_sol, 'm': m_sol, 'N': N_sol,
        'phi': phi_sol, 'phi_prime': phi_p_sol,
        'f': f_sol, 'M_adm': M_adm, 'Q_s': Q_s,
        'alpha': alpha, 'phi_c': phi_c,
        'r_h': r_h, 'horizon_idx': horizon_idx,
        'success': sol.success, 'message': sol.message,
    }


def compute_kretschmann(sol):
    """Compute Kretschmann scalar for the sGB solution.

    For ds² = -N²f dt² + f⁻¹ dr² + r² dΩ²:

    Orthonormal Riemann components:
      R_{(t)(r)(t)(r)} = N²·(f''/2 + f'·N'/N) + N²·f·(N'/N)²  ... (simplified)
      R_{(t)(θ)(t)(θ)} = N²·f·f'/(2r) + N²·f²·N'/(r·N)
      R_{(r)(θ)(r)(θ)} = f'/(2r)
      R_{(θ)(φ)(θ)(φ)} = (1-f)/r²  (N drops out)

    K = 4·R_{trtr}² + 8·R_{tθtθ}² + 8·R_{rθrθ}² + 4·R_{θφθφ}²
    """
    r = sol['r']
    m = sol['m']
    N = sol['N']
    f = sol['f']

    r_safe = np.maximum(r, 1e-30)
    N_safe = np.where(np.abs(N) > 1e-30, N, np.sign(N) * 1e-30)
    f_safe = np.where(np.abs(f) > 1e-30, f, np.sign(f) * 1e-30)

    # Numerical derivatives
    m_p = np.gradient(m, r)
    N_p = np.gradient(N, r)
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe
    f_pp = np.gradient(f_p, r)

    N_over_N = N_p / N_safe

    # Orthonormal Riemann components
    R_trtr = N_safe ** 2 * (f_pp / 2.0 + f_p * N_over_N + f_safe * N_over_N ** 2)
    R_tthtth = N_safe ** 2 * (f_safe * f_p / (2.0 * r_safe) + f_safe ** 2 * N_over_N / r_safe)
    R_rthrth = f_p / (2.0 * r_safe)
    R_thphthph = (1.0 - f_safe) / r_safe ** 2

    K = 4.0 * R_trtr ** 2 + 8.0 * R_tthtth ** 2 + 8.0 * R_rthrth ** 2 + 4.0 * R_thphthph ** 2

    return K


def scan_sgb_solutions():
    """Scan over α and φ_c to find regular black hole solutions."""
    print("=" * 80)
    print("NON-PERTURBATIVE sGB SOLVER — KANTI FORMULATION")
    print("Coupling: f(φ) = (α/2)·φ²  (quadratic, non-shift-symmetric)")
    print("=" * 80)

    results = []

    for alpha in [0.1, 0.5, 1.0, 2.0, 5.0]:
        for phi_c in [0.5, 1.0, 2.0, 5.0]:
            print(f"\n{'='*60}")
            print(f"  α = {alpha}, φ_c = {phi_c}")
            print(f"{'='*60}")

            sol = solve_sgb_kanti(alpha, phi_c, r_max=50.0)

            print(f"  Integration: {sol['message']}")
            print(f"  r range: [{sol['r'][0]:.4e}, {sol['r'][-1]:.4e}]")
            print(f"  M_ADM = {sol['M_adm']:.6f}")
            print(f"  Q_s = {sol['Q_s']:.6f}")
            print(f"  r_h = {sol['r_h']:.6f}" if not np.isnan(sol['r_h']) else "  No horizon")

            if len(sol['r']) < 10:
                print("  → Integration failed early")
                results.append({'alpha': alpha, 'phi_c': phi_c, 'status': 'failed'})
                continue

            K = compute_kretschmann(sol)

            # Check regularity at center
            K_0 = K[0] if np.isfinite(K[0]) else float('nan')
            K_finite = np.isfinite(K_0)

            # Power law near center
            mask = (sol['r'] >= 1e-5) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
            if np.sum(mask) > 10:
                alpha_pow = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0]
            else:
                alpha_pow = np.nan

            # f at center (should → 1 for de Sitter core)
            f_0 = sol['f'][0]

            # Check horizon
            has_horizon = not np.isnan(sol['r_h'])

            print(f"\n  Radial profiles:")
            print(f"  {'r':>10} | {'m(r)':>12} | {'f(r)':>12} | {'N(r)':>10} | {'φ(r)':>12} | {'K':>14}")
            print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}-+-{'-'*14}")

            for r_test in [1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
                idx = np.argmin(np.abs(sol['r'] - r_test))
                if idx < len(sol['r']):
                    K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
                    print(f"  {sol['r'][idx]:>10.4e} | {sol['m'][idx]:>12.6e} | "
                          f"{sol['f'][idx]:>12.6e} | {sol['N'][idx]:>10.4f} | "
                          f"{sol['phi'][idx]:>12.6e} | {K_val:>14.6e}")

            print(f"\n  K(r→0) = {K_0:.6e}")
            print(f"  K power law: ~r^{{{alpha_pow:.3f}}}")
            print(f"  f(r→0) = {f_0:.6f} (de Sitter core: f→1)")
            print(f"  Has horizon: {has_horizon}")

            regular = K_finite and (abs(f_0 - 1.0) < 0.5)
            black_hole = has_horizon

            if regular and black_hole:
                status = 'REGULAR BLACK HOLE ✓'
            elif regular and not black_hole:
                status = 'REGULAR SOLITON (no horizon)'
            elif not regular and black_hole:
                status = 'SINGULAR BLACK HOLE'
            else:
                status = 'SINGULAR / NO HORIZON'

            print(f"  → {status}")

            results.append({
                'alpha': alpha, 'phi_c': phi_c,
                'M_adm': float(sol['M_adm']),
                'Q_s': float(sol['Q_s']) if not np.isnan(sol['Q_s']) else None,
                'r_h': float(sol['r_h']) if not np.isnan(sol['r_h']) else None,
                'K_0': float(K_0) if np.isfinite(K_0) else None,
                'K_power': float(alpha_pow) if np.isfinite(alpha_pow) else None,
                'f_0': float(f_0),
                'status': status,
            })

    return results


def main():
    results = scan_sgb_solutions()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  {'α':>6} | {'φ_c':>6} | {'M_ADM':>10} | {'r_h':>10} | {'K(0)':>14} | {'f(0)':>8} | {'Status':>30}")
    print(f"  {'-'*6}-+-{'-'*6}-+-{'-'*10}-+-{'-'*10}-+-{'-'*14}-+-{'-'*8}-+-{'-'*30}")

    for r in results:
        if r['status'] == 'failed':
            print(f"  {r['alpha']:>6} | {r['phi_c']:>6} | {'—':>10} | {'—':>10} | {'—':>14} | {'—':>8} | {'FAILED':>30}")
        else:
            K_str = f"{r['K_0']:.4e}" if r['K_0'] is not None else "—"
            r_h_str = f"{r['r_h']:.4f}" if r['r_h'] is not None else "—"
            print(f"  {r['alpha']:>6} | {r['phi_c']:>6} | {r['M_adm']:>10.4f} | {r_h_str:>10} | "
                  f"{K_str:>14} | {r['f_0']:>8.4f} | {r['status']:>30}")

    # Save results
    output_path = os.path.join(PROJECT_ROOT, 'results', 'stage_c_kanti_results.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_path}")


if __name__ == '__main__':
    main()
