#!/usr/bin/env python3
"""Non-perturbative shift-symmetric sGB solver — horizon shooting.

For shift-symmetric sGB (f(φ) = α·φ), the scalar equation □φ = α·G_GB
is sourced by curvature. The scalar has a Coulomb profile φ ~ Q_s/r at
infinity, sourced by the Schwarzschild curvature G_GB = 48M²/r⁶.

Strategy: expand around the horizon r = r_h where f = 0, using the
horizon regularity condition, and integrate both outward (to r = ∞)
and inward (to r = 0).

Horizon regularity condition:
At f(r_h) = 0, m(r_h) = r_h/2. For m' to be finite, the source must
vanish: ½(φ')² + α·G_GB·φ'·r/2 = 0, giving:
  φ'(r_h) = -α·G_GB(r_h)·r_h

where G_GB(r_h) = 48M²/r_h⁶ (Schwarzschild value at the horizon).

Metric: ds² = -N²·f·dt² + f⁻¹·dr² + r²·dΩ²
  f = 1 - 2m(r)/r

Field equations (EF-regular form, from Kanti 1996):

  m' = (r²/(2f))·[½(φ')²·N² + α·G_GB·N²·φ'·r/(2f)]  ... has 1/f² in GB term

To avoid the 1/f singularity, we use the horizon expansion and integrate
in both directions from r_h.

The N equation:
  N'/N = (r/(2f))·[½(φ')²·N² - α·G_GB·N²·φ'·r/(2f)]

The scalar equation:
  φ'' + [(2/r) + (f'/(2f)) + (N'/N)]·φ' = -α·G_GB/(f·N²)

At the horizon, all these have 1/f singularities. The horizon expansion
resolves this: we expand in δ = r - r_h and determine the coefficients
from the regularity conditions.

Series expansion at r = r_h + δ:
  m(r) = r_h/2 + m₁·δ + m₂·δ² + ...
  φ(r) = φ_h + φ₁·δ + φ₂·δ² + ...
  N(r) = N_h + N₁·δ + N₂·δ² + ...

From f = 1 - 2m/r:
  f(r) = 1 - (r_h + 2m₁·δ + ...)/(r_h + δ)
       = 1 - (r_h/2 + m₁·δ)/(r_h + δ/2)  (to first order)
       ≈ (1 - 2m₁)·δ/r_h  (to leading order in δ)

So f ≈ (1 - 2m₁)·δ/r_h near the horizon. For f > 0 outside (δ > 0),
we need m₁ < 1/2.

The horizon regularity condition gives:
  φ₁ = φ'(r_h) = -α·G_GB(r_h)·r_h

where G_GB(r_h) = 48·(r_h/2)²/r_h⁶ = 12/r_h⁴ (using m = r_h/2, m' = 0
at leading order for Schwarzschild).

So: φ₁ = -12α/r_h³

The N equation at the horizon gives N_h (free parameter, set by gauge).
The m₁ parameter is determined by matching to the asymptotic mass M.

Shooting parameters: r_h (or M = r_h/2) and α (coupling).
The scalar charge Q_s is determined by the solution.
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def gauss_bonnet_schwarzschild(r, M):
    """G_GB for Schwarzschild: 48M²/r⁶."""
    return 48.0 * M ** 2 / np.maximum(r, 1e-30) ** 6


def gauss_bonnet_general(r, m, m_prime):
    """G_GB for general spherical metric: 48(m - rm'/2)²/r⁶."""
    r_safe = np.maximum(np.abs(r), 1e-30)
    return 48.0 * (m - r_safe * m_prime / 2.0) ** 2 / r_safe ** 6


def sgb_equations_outward(r, y, alpha):
    """Field equations for outward integration (r > r_h, f > 0).

    State: y = [m, N, phi, phi_prime]
    """
    m, N, phi, phi_p = y
    r_safe = max(abs(r), 1e-30)

    f = 1.0 - 2.0 * m / r_safe
    f_reg = f if abs(f) > 1e-14 else 1e-14 * (1.0 if f >= 0 else -1.0)

    # Estimate m' for G_GB
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 * N ** 2 / f_reg
    G_GB = gauss_bonnet_general(r_safe, m, m_p_est)

    # Full m'
    rho = 0.5 * phi_p ** 2 * N ** 2
    gb_src = alpha * G_GB * N ** 2 * phi_p * r_safe / (2.0 * f_reg)
    m_p = r_safe ** 2 / (2.0 * f_reg) * (rho + gb_src)

    # Recompute G_GB
    G_GB_full = gauss_bonnet_general(r_safe, m, m_p)

    # N'/N
    gb_N = alpha * G_GB_full * N ** 2 * phi_p * r_safe / (2.0 * f_reg)
    N_p_over_N = (r_safe / (2.0 * f_reg)) * (rho - gb_N)
    N_p = N_p_over_N * N

    # Scalar equation
    f_p = 2.0 * m / r_safe ** 2 - 2.0 * m_p / r_safe
    coeff = 2.0 / r_safe + f_p / (2.0 * f_reg) + N_p_over_N
    phi_pp = -alpha * G_GB_full / (f_reg * N ** 2) - coeff * phi_p

    return [m_p, N_p, phi_p, phi_pp]


def sgb_equations_inward(r, y, alpha):
    """Field equations for inward integration (r < r_h, f < 0).

    Same equations, but f < 0. The regularization handles the sign.
    """
    return sgb_equations_outward(r, y, alpha)


def horizon_expansion(r_h, alpha, N_h=1.0, phi_h=0.0):
    """Compute series coefficients at the horizon r = r_h.

    Returns initial conditions at r = r_h + δ (small δ) for outward
    integration, and at r = r_h - δ for inward integration.

    Horizon regularity:
      m(r_h) = r_h/2
      φ'(r_h) = -α·G_GB(r_h)·r_h = -12α/r_h³
      N(r_h) = N_h (gauge)

    The m₁ coefficient (dm/dr at r_h) is a free parameter determined
    by shooting to match the asymptotic mass. For Schwarzschild, m₁ = 0.
    For sGB, m₁ ≠ 0 due to the scalar backreaction.

    From the m' equation at the horizon (using L'Hôpital for the 1/f):
      m'(r_h) = m₁ = lim_{r→r_h} r²/(2f)·[½(φ')²·N² + α·G_GB·N²·φ'·r/(2f)]

    The numerator → 0 by the regularity condition, and f → 0, so we need
    to expand to next order. This gives m₁ as a function of the other
    parameters.

    For simplicity, we use m₁ as a shooting parameter and determine it
    by matching to M_ADM at infinity.
    """
    M = r_h / 2  # horizon mass

    # G_GB at horizon (Schwarzschild value)
    G_GB_h = gauss_bonnet_schwarzschild(r_h, M)

    # Horizon regularity: φ'(r_h) = -α·G_GB·r_h
    phi_1 = -alpha * G_GB_h * r_h

    # φ_h is a free parameter (scalar value at horizon)
    # For shift-symmetric sGB, φ_h is determined by the integration

    # N_h is gauge (set to 1)
    # m₁ is the shooting parameter

    return {
        'r_h': r_h, 'M': M, 'alpha': alpha,
        'phi_h': phi_h, 'phi_1': phi_1,
        'N_h': N_h, 'G_GB_h': G_GB_h,
    }


def solve_outward(r_h, alpha, m1, phi_h, N_h=1.0, r_max=50.0, n_points=5000):
    """Integrate outward from r_h + δ to r_max.

    Parameters:
        r_h: horizon radius
        alpha: GB coupling
        m1: dm/dr at horizon (shooting parameter)
        phi_h: scalar value at horizon
        N_h: lapse at horizon (gauge)
    """
    hz = horizon_expansion(r_h, alpha, N_h, phi_h)

    # Initial conditions at r = r_h + δ
    delta = 1e-4 * r_h
    r_start = r_h + delta

    # Series: m = r_h/2 + m1·δ, φ = φ_h + phi_1·δ, N = N_h + N1·δ
    m_start = r_h / 2.0 + m1 * delta
    phi_start = phi_h + hz['phi_1'] * delta
    phi_p_start = hz['phi_1']  # φ'(r_h)
    N_start = N_h  # N1 is higher order

    y0 = [m_start, N_start, phi_start, phi_p_start]

    r_span = (r_start, r_max)
    r_eval = np.logspace(np.log10(r_start), np.log10(r_max), n_points)

    sol = solve_ivp(
        sgb_equations_outward, r_span, y0,
        args=(alpha,),
        t_eval=r_eval,
        method='Radau',
        rtol=1e-10,
        atol=1e-12,
        max_step=0.5,
    )

    return sol


def solve_inward(r_h, alpha, m1, phi_h, N_h=1.0, r_min=1e-6, n_points=5000):
    """Integrate inward from r_h - δ to r_min.

    The equations are the same, but f < 0 inside the horizon.
    """
    hz = horizon_expansion(r_h, alpha, N_h, phi_h)

    delta = 1e-4 * r_h
    r_start = r_h - delta

    m_start = r_h / 2.0 + m1 * (-delta)  # m(r_h - δ)
    phi_start = phi_h + hz['phi_1'] * (-delta)
    phi_p_start = hz['phi_1']
    N_start = N_h

    y0 = [m_start, N_start, phi_start, phi_p_start]

    r_span = (r_start, r_min)  # decreasing r
    r_eval = np.logspace(np.log10(r_start), np.log10(r_min), n_points)

    sol = solve_ivp(
        sgb_equations_inward, r_span, y0,
        args=(alpha,),
        t_eval=r_eval,
        method='Radau',
        rtol=1e-10,
        atol=1e-12,
        max_step=0.1,
    )

    return sol


def compute_kretschmann(r, m, N):
    """Compute Kretschmann scalar for ds² = -N²f dt² + f⁻¹ dr² + r² dΩ²."""
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


def shoot_for_mass(r_h, alpha, M_target, phi_h=0.0, N_h=1.0):
    """Shoot over m1 to match M_ADM = M_target at infinity.

    Returns the solution that best matches the target mass.
    """
    def mass_residual(m1):
        sol = solve_outward(r_h, alpha, m1, phi_h, N_h, r_max=30.0, n_points=2000)
        if not sol.success or len(sol.t) < 10:
            return 1e6
        M_adm = sol.y[0][-1]
        return M_adm - M_target

    # Scan m1 to find the right value
    # For Schwarzschild: m1 = 0, M_ADM = r_h/2
    # For sGB: m1 < 0 (mass decreases inward due to scalar backreaction)
    print(f"  Shooting for M_target = {M_target}, r_h = {r_h}, α = {alpha}")

    m1_values = np.linspace(-0.5, 0.5, 21)
    residuals = []
    for m1 in m1_values:
        res = mass_residual(m1)
        residuals.append(res)
        print(f"    m1 = {m1:+.3f}, M_ADM - M_target = {res:+.6f}")

    # Find sign change
    residuals = np.array(residuals)
    sign_changes = np.where(np.diff(np.sign(residuals)))[0]

    if len(sign_changes) == 0:
        print(f"    No sign change found, using minimum |residual|")
        best_idx = np.argmin(np.abs(residuals))
        m1_best = m1_values[best_idx]
    else:
        # Linear interpolation
        idx = sign_changes[0]
        r1, r2 = residuals[idx], residuals[idx + 1]
        m1_1, m1_2 = m1_values[idx], m1_values[idx + 1]
        m1_best = m1_1 - r1 * (m1_2 - m1_1) / (r2 - r1)

    print(f"  Best m1 = {m1_best:.6f}")

    return m1_best


def solve_full(r_h, alpha, m1, phi_h=0.0, N_h=1.0):
    """Solve both outward and inward, combine into full solution."""
    sol_out = solve_outward(r_h, alpha, m1, phi_h, N_h, r_max=50.0, n_points=10000)
    sol_in = solve_inward(r_h, alpha, m1, phi_h, N_h, r_min=1e-6, n_points=10000)

    # Combine (inward first, then outward)
    r_in = sol_in.t[::-1]  # reverse to go from small r to r_h
    m_in = sol_in.y[0][::-1]
    N_in = sol_in.y[1][::-1]
    phi_in = sol_in.y[2][::-1]
    phi_p_in = sol_in.y[3][::-1]

    r_out = sol_out.t
    m_out = sol_out.y[0]
    N_out = sol_out.y[1]
    phi_out = sol_out.y[2]
    phi_p_out = sol_out.y[3]

    r = np.concatenate([r_in, r_out])
    m = np.concatenate([m_in, m_out])
    N = np.concatenate([N_in, N_out])
    phi = np.concatenate([phi_in, phi_out])
    phi_p = np.concatenate([phi_p_in, phi_p_out])

    # Sort by r
    sort_idx = np.argsort(r)
    r = r[sort_idx]
    m = m[sort_idx]
    N = N[sort_idx]
    phi = phi[sort_idx]
    phi_p = phi_p[sort_idx]

    f = 1.0 - 2.0 * m / np.maximum(r, 1e-30)
    K = compute_kretschmann(r, m, N)

    return {
        'r': r, 'm': m, 'N': N, 'phi': phi, 'phi_prime': phi_p,
        'f': f, 'K': K,
        'r_h': r_h, 'alpha': alpha, 'm1': m1, 'phi_h': phi_h,
        'M_adm': m[-1],
        'out_success': sol_out.success,
        'in_success': sol_in.success,
    }


def main():
    print("=" * 80)
    print("NON-PERTURBATIVE SHIFT-SYMMETRIC sGB — HORIZON SHOOTING")
    print("Coupling: f(φ) = α·φ  (shift-symmetric)")
    print("Strategy: Expand at horizon, integrate both directions")
    print("=" * 80)

    M_target = 1.0
    r_h = 2.0 * M_target  # Schwarzschild horizon

    results = []

    for alpha in [0.1, 0.5, 1.0, 2.0]:
        print(f"\n{'='*60}")
        print(f"  α = {alpha}, M_target = {M_target}, r_h = {r_h}")
        print(f"{'='*60}")

        # Shoot for m1 to match M_target
        m1_best = shoot_for_mass(r_h, alpha, M_target, phi_h=0.0)

        # Solve full problem
        sol = solve_full(r_h, alpha, m1_best, phi_h=0.0)

        print(f"\n  Outward success: {sol['out_success']}")
        print(f"  Inward success: {sol['in_success']}")
        print(f"  M_ADM = {sol['M_adm']:.6f}")

        # Check regularity at center
        K_0 = sol['K'][0]
        f_0 = sol['f'][0]
        m_0 = sol['m'][0]
        phi_0 = sol['phi'][0]

        print(f"\n  Center values (r → 0):")
        print(f"    m(0) = {m_0:.6e}")
        print(f"    f(0) = {f_0:.6f}")
        print(f"    φ(0) = {phi_0:.6f}")
        print(f"    K(0) = {K_0:.6e}")

        # Power law near center
        mask = (sol['r'] >= 1e-5) & (sol['r'] <= 0.1) & np.isfinite(sol['K']) & (sol['K'] > 0)
        if np.sum(mask) > 10:
            alpha_pow = np.polyfit(np.log(sol['r'][mask]), np.log(sol['K'][mask]), 1)[0]
            print(f"    K power law: ~r^{{{alpha_pow:.3f}}}")
        else:
            alpha_pow = np.nan
            print(f"    K power law: N/A")

        # Radial profile
        print(f"\n  Radial profile:")
        print(f"  {'r':>10} | {'m(r)':>12} | {'f(r)':>12} | {'N(r)':>10} | {'φ(r)':>12} | {'K':>14}")
        print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}-+-{'-'*14}")

        for r_test in [1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 3.0, 5.0, 10.0, 50.0]:
            idx = np.argmin(np.abs(sol['r'] - r_test))
            if idx < len(sol['r']):
                K_val = sol['K'][idx] if np.isfinite(sol['K'][idx]) else float('nan')
                print(f"  {sol['r'][idx]:>10.4e} | {sol['m'][idx]:>12.6e} | "
                      f"{sol['f'][idx]:>12.6e} | {sol['N'][idx]:>10.4f} | "
                      f"{sol['phi'][idx]:>12.6e} | {K_val:>14.6e}")

        # Regularity assessment
        K_finite = np.isfinite(K_0) and K_0 < 1e10
        de_sitter = abs(f_0 - 1.0) < 0.5
        m_zero = abs(m_0) < 0.01

        if K_finite and de_sitter and m_zero:
            status = 'REGULAR BLACK HOLE ✓'
        elif K_finite and de_sitter:
            status = 'REGULAR (m not quite 0)'
        elif K_finite:
            status = 'K finite but no de Sitter core'
        else:
            status = 'SINGULAR'

        print(f"\n  → {status}")

        results.append({
            'alpha': alpha, 'r_h': r_h, 'm1': m1_best,
            'M_adm': float(sol['M_adm']),
            'K_0': float(K_0) if np.isfinite(K_0) else None,
            'K_power': float(alpha_pow) if np.isfinite(alpha_pow) else None,
            'f_0': float(f_0),
            'm_0': float(m_0),
            'phi_0': float(phi_0),
            'status': status,
        })

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  {'α':>6} | {'m1':>10} | {'M_ADM':>10} | {'K(0)':>14} | {'f(0)':>8} | {'m(0)':>12} | {'Status':>30}")
    print(f"  {'-'*6}-+-{'-'*10}-+-{'-'*10}-+-{'-'*14}-+-{'-'*8}-+-{'-'*12}-+-{'-'*30}")

    for r in results:
        K_str = f"{r['K_0']:.4e}" if r['K_0'] is not None else "—"
        print(f"  {r['alpha']:>6} | {r['m1']:>10.6f} | {r['M_adm']:>10.4f} | "
              f"{K_str:>14} | {r['f_0']:>8.4f} | {r['m_0']:>12.4e} | {r['status']:>30}")

    # Save
    output_path = os.path.join(PROJECT_ROOT, 'results', 'stage_c_horizon_shooting.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_path}")


if __name__ == '__main__':
    main()
