#!/usr/bin/env python3
"""Stage C v2: Non-perturbative sGB — inward integration.

The outward integration fails because with m(0)=0, the Gauss-Bonnet
invariant G_GB = 0, so there's no source for the scalar field.

The correct approach (Kanti 1996, Kleihaus et al. 2011) is to integrate
INWARD from large r, where the solution is approximately Schwarzschild
with a scalar field φ ~ Q_s/r.

Boundary conditions at r → ∞:
  m(r) → M (ADM mass)
  φ(r) → Q_s/r (scalar charge)
  Φ(r) → 0 (asymptotically flat)

We shoot from r = r_max with these asymptotic conditions and integrate
inward. The regularity at r = 0 (m → 0, φ → finite, φ' → 0) is then
a RESULT that we verify, not a boundary condition we impose.

The scalar charge Q_s is a free parameter that we tune to achieve
regularity at r = 0.
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def gauss_bonnet_invariant_f(r, m, m_p):
    """Gauss-Bonnet invariant for f = 1 - 2m(r)/r.

    G_GB = (48/r⁶) · (m - r·m'/2)²

    This is exact for the metric ds² = -f dt² + f^{-1} dr² + r² dΩ²
    when the lapse Φ = 0 (which is the case for Schwarzschild).
    For Φ ≠ 0, there are additional terms, but this is the dominant
    contribution and is exact in the Schwarzschild limit.
    """
    r_safe = np.maximum(np.abs(r), 1e-30)
    return 48.0 * (m - r_safe * m_p / 2.0) ** 2 / r_safe ** 6


def sgb_field_equations_inward(r, y, eta):
    """Field equations for non-perturbative sGB, integrated inward.

    State vector: y = [m, Phi, phi, phi_prime]
    Note: r is DECREASING (integrating inward).

    Field equations (same as outward, but r decreases):
      m' = r²/(2f) · [½(φ')² + η · G_GB · φ' · r/2]
      Φ' = [½(φ')² r² + η · G_GB · φ' · r³/(4f)] / [r²f]
      φ'' = -η G_GB / f - [(2/r) + (m'r - m)/(r²f) - Φ'] φ'

    where f = 1 - 2m/r, G_GB = 48(m - rm'/2)²/r⁶.
    """
    m, Phi, phi, phi_p = y

    r_safe = max(abs(r), 1e-30)
    f = 1.0 - 2.0 * m / r_safe
    f_safe = f if abs(f) > 1e-10 else np.sign(f) * 1e-10

    # First estimate m' (using scalar kinetic only for G_GB estimate)
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 / f_safe

    # Gauss-Bonnet invariant
    G_GB = gauss_bonnet_invariant_f(r_safe, m, m_p_est)

    # Full m' including GB source
    rho_scalar = 0.5 * phi_p ** 2
    gb_source = eta * G_GB * phi_p * r_safe / 2.0
    m_p = r_safe ** 2 / (2.0 * f_safe) * (rho_scalar + gb_source)

    # Recompute G_GB with the full m'
    G_GB_full = gauss_bonnet_invariant_f(r_safe, m, m_p)

    # Phi' (lapse)
    Phi_p = (0.5 * phi_p ** 2 * r_safe ** 2
             + eta * G_GB_full * phi_p * r_safe ** 3 / (4.0 * f_safe)) / (r_safe ** 2 * f_safe)

    # phi'' (scalar equation)
    coeff = (2.0 / r_safe
             + (m_p * r_safe - m) / (r_safe ** 2 * f_safe)
             - Phi_p)
    phi_pp = -eta * G_GB_full / f_safe - coeff * phi_p

    return [m_p, Phi_p, phi_p, phi_pp]


def solve_sgb_inward(eta, M=1.0, Q_s=0.1, r_max=50.0, r_min=1e-4,
                      n_points=10000):
    """Solve sGB by integrating inward from Schwarzschild-like conditions.

    Boundary conditions at r = r_max:
      m(r_max) = M (ADM mass)
      φ(r_max) = Q_s / r_max (scalar charge)
      φ'(r_max) = -Q_s / r_max²
      Φ(r_max) = 0 (asymptotically flat)

    The scalar charge Q_s is the shooting parameter. We tune it to
    achieve regularity at r → 0.
    """
    # Initial conditions at r_max
    m_init = M
    Phi_init = 0.0
    phi_init = Q_s / r_max
    phi_p_init = -Q_s / r_max ** 2

    y0 = [m_init, Phi_init, phi_init, phi_p_init]

    # Integrate inward (r decreasing)
    r_span = (r_max, r_min)
    r_eval = np.logspace(np.log10(r_max), np.log10(r_min), n_points)

    def event_horizon(r, y, eta):
        """Detect horizon crossing (f = 0)."""
        m = y[0]
        f = 1.0 - 2.0 * m / max(abs(r), 1e-30)
        return f
    event_horizon.terminal = False
    event_horizon.direction = 1  # f increasing (going inward, f decreases)

    def event_regular(r, y, eta):
        """Detect m → 0 (regular center)."""
        return y[0]  # m → 0
    event_regular.terminal = False
    event_regular.direction = -1

    sol = solve_ivp(
        sgb_field_equations_inward, r_span, y0,
        args=(eta,),
        t_eval=r_eval,
        method='RK45',
        rtol=1e-10,
        atol=1e-12,
        events=[event_horizon],
        max_step=0.5,
    )

    r_sol = sol.t[::-1]  # reverse to ascending order
    m_sol = sol.y[0][::-1]
    Phi_sol = sol.y[1][::-1]
    phi_sol = sol.y[2][::-1]
    phi_p_sol = sol.y[3][::-1]

    f_sol = 1.0 - 2.0 * m_sol / r_sol

    return {
        'r': r_sol, 'm': m_sol, 'Phi': Phi_sol,
        'phi': phi_sol, 'phi_prime': phi_p_sol,
        'f': f_sol, 'M': M, 'Q_s': Q_s, 'eta': eta,
        'success': sol.success, 'message': sol.message,
    }


def compute_curvature_sgb(sol):
    """Compute Kretschmann scalar for the sGB solution."""
    r = sol['r']
    m = sol['m']
    Phi = sol['Phi']
    f = sol['f']

    g_tt = -f * np.exp(2 * Phi)
    g_rr = 1.0 / np.where(np.abs(f) > 1e-10, f, np.sign(f) * 1e-10)
    g_thth = r ** 2

    g_tt_p = np.gradient(g_tt, r)
    g_rr_p = np.gradient(g_rr, r)
    g_thth_p = np.gradient(g_thth, r)

    g_tt_safe = np.where(np.abs(g_tt) > 1e-50, g_tt, np.nan)
    g_rr_safe = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
    g_thth_safe = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

    ginv_tt = 1.0 / g_tt_safe
    ginv_rr = 1.0 / g_rr_safe
    ginv_thth = 1.0 / g_thth_safe

    Gt_tr = 0.5 * ginv_tt * g_tt_p
    Gr_rr = 0.5 * ginv_rr * g_rr_p
    Gth_rth = 0.5 * ginv_thth * g_thth_p

    dGt_tr = np.gradient(Gt_tr, r)
    dGth_rth = np.gradient(Gth_rth, r)

    R_trtr = g_tt_safe * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr)
    R_rthrth = g_thth * (dGth_rth + Gth_rth ** 2 - Gth_rth * Gr_rr)
    R_tthtth = -0.25 * g_tt_p * g_thth_p / g_rr_safe
    R_thphthph = g_thth * (1.0 - ginv_rr * g_thth_p ** 2 / (4.0 * g_thth_safe))

    E = R_trtr / (g_tt_safe * g_rr_safe)
    F_t = R_tthtth / (g_tt_safe * g_thth_safe)
    F_r = R_rthrth / (g_rr_safe * g_thth_safe)
    G_comp = R_thphthph / g_thth_safe ** 2

    K = 4.0 * E ** 2 + 8.0 * F_t ** 2 + 8.0 * F_r ** 2 + 4.0 * G_comp ** 2

    return K


def scan_q_s():
    """Scan over scalar charge Q_s to find regular solutions."""
    print("=" * 80)
    print("STAGE C v2: INWARD INTEGRATION — SCAN OVER Q_s")
    print("=" * 80)

    eta = 1.0
    M = 1.0

    print(f"\n  η = {eta}, M = {M}")
    print(f"  Scanning Q_s to find regular solution...\n")

    print(f"  {'Q_s':>8} | {'m(r_min)':>12} | {'f(r_min)':>12} | {'φ(r_min)':>12} | "
          f"{'K(r_min)':>14} | {'K power':>10} | {'Regular?':>10}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}-+-{'-'*10}")

    for Q_s in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        sol = solve_sgb_inward(eta, M=M, Q_s=Q_s, r_max=50.0, r_min=1e-4)

        if len(sol['r']) < 10:
            print(f"  {Q_s:>8.3f} | {'FAIL':>12} |")
            continue

        K = compute_curvature_sgb(sol)

        m_min = sol['m'][0]
        f_min = sol['f'][0]
        phi_min = sol['phi'][0]
        K_min = K[0] if np.isfinite(K[0]) else float('nan')

        # Power law
        mask = (sol['r'] >= 1e-3) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0]
        else:
            alpha = np.nan

        regular = np.isfinite(K_min) and alpha > -1 and abs(m_min) < 1.0

        print(f"  {Q_s:>8.3f} | {m_min:>12.6e} | {f_min:>12.6e} | {phi_min:>12.6e} | "
              f"{K_min:>14.6e} | {alpha:>10.3f} | {'YES' if regular else 'NO':>10}")

    # Detailed analysis for a promising Q_s
    print(f"\n{'='*80}")
    print(f"  DETAILED ANALYSIS: η=1.0, Q_s=0.1")
    print(f"{'='*80}")

    sol = solve_sgb_inward(eta=1.0, M=1.0, Q_s=0.1, r_max=50.0, r_min=1e-4)
    K = compute_curvature_sgb(sol)

    print(f"\n  {'r':>10} | {'m(r)':>12} | {'f(r)':>12} | {'φ(r)':>12} | {'φ\'(r)':>12} | {'K':>14}")
    print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}")

    for r_test in [1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0]:
        idx = np.argmin(np.abs(sol['r'] - r_test))
        if idx < len(sol['r']):
            K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
            print(f"  {sol['r'][idx]:>10.4e} | {sol['m'][idx]:>12.6e} | "
                  f"{sol['f'][idx]:>12.6e} | {sol['phi'][idx]:>12.6e} | "
                  f"{sol['phi_prime'][idx]:>12.6e} | {K_val:>14.6e}")


def scan_eta():
    """Scan over coupling η for fixed Q_s."""
    print(f"\n{'='*80}")
    print(f"  SCAN OVER η (Q_s = 0.1, M = 1.0)")
    print(f"{'='*80}")

    Q_s = 0.1
    M = 1.0

    print(f"\n  {'η':>6} | {'m(r_min)':>12} | {'f(r_min)':>12} | {'K(r_min)':>14} | {'K power':>10} | {'Regular?':>10}")
    print(f"  {'-'*6}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}-+-{'-'*10}")

    for eta in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
        sol = solve_sgb_inward(eta, M=M, Q_s=Q_s, r_max=50.0, r_min=1e-4)

        if len(sol['r']) < 10:
            print(f"  {eta:>6.2f} | {'FAIL':>12} |")
            continue

        K = compute_curvature_sgb(sol)

        m_min = sol['m'][0]
        f_min = sol['f'][0]
        K_min = K[0] if np.isfinite(K[0]) else float('nan')

        mask = (sol['r'] >= 1e-3) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0]
        else:
            alpha = np.nan

        regular = np.isfinite(K_min) and alpha > -1

        print(f"  {eta:>6.2f} | {m_min:>12.6e} | {f_min:>12.6e} | "
              f"{K_min:>14.6e} | {alpha:>10.3f} | {'YES' if regular else 'NO':>10}")


def main():
    scan_q_s()
    scan_eta()

    print(f"\n{'='*80}")
    print(f"STAGE C v2 STATUS")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
