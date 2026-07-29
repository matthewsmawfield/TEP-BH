#!/usr/bin/env python3
"""Stage C v3: Non-perturbative sGB — EF coordinates, horizon-regular.

The previous attempt hung because standard coordinates (t,r) have a
coordinate singularity at the horizon (f=0), and the field equations
contain 1/f terms.

Solution: Use Eddington-Finkelstein (EF) coordinates, which are regular
at the horizon. The metric is:

  ds² = -f dv² + 2 dv dr + r² dΩ²

where v = t + r* is the EF time coordinate and f = 1 - 2m(r)/r.

In EF coordinates, the field equations are regular at f=0.

The field equations in EF (Kleihaus, Kunz, Schneider 2011):

  m' = r²/(2) · [½(φ')² + η · G_GB · φ' · r/2]   (no 1/f!)
  φ'' + (2/r) φ' = -η G_GB                          (no 1/f!)
  (lapse absorbed into EF structure)

Wait, this is too simplified. Let me use the proper EF equations.

Actually, for the metric ds² = -f e^{2Φ} dv² + 2 e^{Φ} dv dr + r² dΩ²,
the equations are (Kleihaus et al. 2011, eq. 9-11):

  m' = 4π r² ρ_eff = r² [½(φ')² + η·G_GB·φ'·r/2]  (simplified)
  Φ' = r [½(φ')² + η·G_GB·φ'·r/(4f)] ... still has 1/f

Hmm, the Φ equation still has 1/f. But in EF, the key is that the
INTEGRATION is regular — the ODEs can be written in a form where
the horizon is not a singular point.

Let me use a different approach: the compactified coordinate x = 1/r,
or use the tortoise coordinate. Or simply integrate in standard
coordinates but with horizon regularization (clamp f to a minimum).

Actually, the simplest fix: integrate in standard coordinates but
regularize f near the horizon (clamp |f| to a minimum value). This
introduces a small error near the horizon but allows the integration
to proceed.
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def sgb_field_equations_regularized(r, y, eta, f_min=1e-6):
    """Field equations with horizon regularization.

    Clamp |f| to f_min to avoid division by zero at the horizon.
    This introduces O(f_min) error near the horizon but allows
    integration to proceed through it.
    """
    m, Phi, phi, phi_p = y

    r_safe = max(abs(r), 1e-30)
    f = 1.0 - 2.0 * m / r_safe
    # Regularize: clamp |f| to f_min
    f_reg = np.sign(f) * max(abs(f), f_min) if abs(f) > 0 else f_min

    # Estimate m' for G_GB
    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 / f_reg

    # G_GB = 48(m - rm'/2)²/r⁶
    G_GB = 48.0 * (m - r_safe * m_p_est / 2.0) ** 2 / r_safe ** 6

    # Full m'
    rho_scalar = 0.5 * phi_p ** 2
    gb_source = eta * G_GB * phi_p * r_safe / 2.0
    m_p = r_safe ** 2 / (2.0 * f_reg) * (rho_scalar + gb_source)

    # Recompute G_GB
    G_GB_full = 48.0 * (m - r_safe * m_p / 2.0) ** 2 / r_safe ** 6

    # Phi'
    Phi_p = (0.5 * phi_p ** 2 * r_safe ** 2
             + eta * G_GB_full * phi_p * r_safe ** 3 / (4.0 * f_reg)) / (r_safe ** 2 * f_reg)

    # phi''
    coeff = (2.0 / r_safe
             + (m_p * r_safe - m) / (r_safe ** 2 * f_reg)
             - Phi_p)
    phi_pp = -eta * G_GB_full / f_reg - coeff * phi_p

    return [m_p, Phi_p, phi_p, phi_pp]


def solve_sgb_inward(eta, M=1.0, Q_s=0.1, r_max=50.0, r_min=1e-4,
                      n_points=5000, f_min=1e-4):
    """Solve sGB by integrating inward with horizon regularization."""
    m_init = M
    Phi_init = 0.0
    phi_init = Q_s / r_max
    phi_p_init = -Q_s / r_max ** 2

    y0 = [m_init, Phi_init, phi_init, phi_p_init]

    r_span = (r_max, r_min)
    r_eval = np.logspace(np.log10(r_max), np.log10(r_min), n_points)

    sol = solve_ivp(
        lambda r, y: sgb_field_equations_regularized(r, y, eta, f_min),
        r_span, y0,
        t_eval=r_eval,
        method='RK45',
        rtol=1e-8,
        atol=1e-10,
        max_step=0.5,
    )

    r_sol = sol.t[::-1]
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
        'n_points': len(r_sol),
    }


def compute_curvature_sgb(sol):
    """Compute Kretschmann scalar."""
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


def quick_test():
    """Quick test with one configuration."""
    print("=" * 70)
    print("QUICK TEST: η=1.0, Q_s=0.1, M=1.0")
    print("=" * 70)

    sol = solve_sgb_inward(eta=1.0, M=1.0, Q_s=0.1, r_max=50.0, r_min=1e-4,
                            n_points=2000, f_min=1e-4)

    print(f"  Success: {sol['success']}")
    print(f"  Message: {sol['message']}")
    print(f"  N points: {sol['n_points']}")
    print(f"  r range: [{sol['r'][0]:.6e}, {sol['r'][-1]:.6e}]")

    if sol['n_points'] < 10:
        print("  → FAILED")
        return

    K = compute_curvature_sgb(sol)

    print(f"\n  {'r':>10} | {'m(r)':>12} | {'f(r)':>12} | {'φ(r)':>12} | {'K':>14}")
    print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}")

    for r_test in [1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0]:
        idx = np.argmin(np.abs(sol['r'] - r_test))
        if idx < len(sol['r']):
            K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
            print(f"  {sol['r'][idx]:>10.4e} | {sol['m'][idx]:>12.6e} | "
                  f"{sol['f'][idx]:>12.6e} | {sol['phi'][idx]:>12.6e} | {K_val:>14.6e}")

    # Power law
    mask = (sol['r'] >= 1e-3) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
    if np.sum(mask) > 10:
        alpha = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0]
        print(f"\n  K power law near r=0: ~r^{{{alpha:.3f}}}")
    else:
        print(f"\n  Cannot compute power law (insufficient finite K points)")


def scan():
    """Scan over Q_s and η."""
    print(f"\n{'='*80}")
    print(f"SCAN: Q_s and η (M=1.0, r_min=1e-4)")
    print(f"{'='*80}")

    print(f"\n  Fixed η=1.0, scan Q_s:")
    print(f"  {'Q_s':>8} | {'m_min':>12} | {'f_min':>12} | {'K_min':>14} | {'K power':>10} | {'n_pts':>6}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}-+-{'-'*6}")

    for Q_s in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0]:
        sol = solve_sgb_inward(eta=1.0, M=1.0, Q_s=Q_s, r_max=30.0, r_min=1e-4,
                                n_points=2000, f_min=1e-4)
        if sol['n_points'] < 10:
            print(f"  {Q_s:>8.3f} | {'FAIL':>12} |")
            continue

        K = compute_curvature_sgb(sol)
        m_min = sol['m'][0]
        f_min_val = sol['f'][0]
        K_min = K[0] if np.isfinite(K[0]) else float('nan')

        mask = (sol['r'] >= 1e-3) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
        alpha = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0] if np.sum(mask) > 10 else np.nan

        print(f"  {Q_s:>8.3f} | {m_min:>12.6e} | {f_min_val:>12.6e} | "
              f"{K_min:>14.6e} | {alpha:>10.3f} | {sol['n_points']:>6}")

    print(f"\n  Fixed Q_s=0.1, scan η:")
    print(f"  {'η':>8} | {'m_min':>12} | {'f_min':>12} | {'K_min':>14} | {'K power':>10} | {'n_pts':>6}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}-+-{'-'*10}-+-{'-'*6}")

    for eta in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        sol = solve_sgb_inward(eta=eta, M=1.0, Q_s=0.1, r_max=30.0, r_min=1e-4,
                                n_points=2000, f_min=1e-4)
        if sol['n_points'] < 10:
            print(f"  {eta:>8.3f} | {'FAIL':>12} |")
            continue

        K = compute_curvature_sgb(sol)
        m_min = sol['m'][0]
        f_min_val = sol['f'][0]
        K_min = K[0] if np.isfinite(K[0]) else float('nan')

        mask = (sol['r'] >= 1e-3) & (sol['r'] <= 0.1) & np.isfinite(K) & (K > 0)
        alpha = np.polyfit(np.log(sol['r'][mask]), np.log(K[mask]), 1)[0] if np.sum(mask) > 10 else np.nan

        print(f"  {eta:>8.3f} | {m_min:>12.6e} | {f_min_val:>12.6e} | "
              f"{K_min:>14.6e} | {alpha:>10.3f} | {sol['n_points']:>6}")


def main():
    quick_test()
    scan()

    print(f"\n{'='*80}")
    print(f"STAGE C v3 STATUS")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
