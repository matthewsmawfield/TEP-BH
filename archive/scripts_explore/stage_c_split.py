#!/usr/bin/env python3
"""Stage C v4: Non-perturbative sGB — split integration at horizon.

The ODE is stiff near the horizon (f→0). Solution: split the integration
into exterior (r > r_h) and interior (r < r_h) regions, using a stiff
solver (Radau) and matching at the horizon.

Also use compactified coordinate u = 1/r to improve resolution at large r
and avoid the horizon issue by integrating in the interior only from
just inside the horizon.
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


def sgb_eqs(r, y, eta, f_min=1e-3):
    """Field equations with regularization. Returns dy/dr."""
    m, Phi, phi, phi_p = y
    r_safe = max(abs(r), 1e-30)
    f = 1.0 - 2.0 * m / r_safe
    f_reg = np.sign(f) * max(abs(f), f_min) if f != 0 else f_min

    m_p_est = 0.5 * r_safe ** 2 * phi_p ** 2 / f_reg
    G_GB = 48.0 * (m - r_safe * m_p_est / 2.0) ** 2 / r_safe ** 6

    rho_s = 0.5 * phi_p ** 2
    gb_s = eta * G_GB * phi_p * r_safe / 2.0
    m_p = r_safe ** 2 / (2.0 * f_reg) * (rho_s + gb_s)

    G_GB_f = 48.0 * (m - r_safe * m_p / 2.0) ** 2 / r_safe ** 6
    Phi_p = (0.5 * phi_p ** 2 * r_safe ** 2
             + eta * G_GB_f * phi_p * r_safe ** 3 / (4.0 * f_reg)) / (r_safe ** 2 * f_reg)

    coeff = (2.0 / r_safe + (m_p * r_safe - m) / (r_safe ** 2 * f_reg) - Phi_p)
    phi_pp = -eta * G_GB_f / f_reg - coeff * phi_p

    return [m_p, Phi_p, phi_p, phi_pp]


def solve_exterior(eta, M, Q_s, r_max=30.0, r_h_plus=2.05, n=1000):
    """Integrate from r_max inward to just outside horizon."""
    y0 = [M, 0.0, Q_s / r_max, -Q_s / r_max ** 2]
    r_eval = np.linspace(r_max, r_h_plus, n)
    sol = solve_ivp(lambda r, y: sgb_eqs(r, y, eta), (r_max, r_h_plus), y0,
                    t_eval=r_eval, method='Radau', rtol=1e-6, atol=1e-8)
    return sol


def solve_interior(eta, y_horizon, r_h_minus=1.95, r_min=1e-3, n=2000):
    """Integrate from just inside horizon inward to r_min."""
    r_eval = np.linspace(r_h_minus, r_min, n)
    sol = solve_ivp(lambda r, y: sgb_eqs(r, y, eta), (r_h_minus, r_min), y_horizon,
                    t_eval=r_eval, method='Radau', rtol=1e-6, atol=1e-8)
    return sol


def compute_K(r, m, Phi):
    """Compute Kretschmann."""
    f = 1.0 - 2.0 * m / r
    g_tt = -f * np.exp(2 * Phi)
    g_rr = 1.0 / np.where(np.abs(f) > 1e-10, f, np.sign(f) * 1e-10)
    g_thth = r ** 2

    g_tt_p = np.gradient(g_tt, r)
    g_rr_p = np.gradient(g_rr, r)
    g_thth_p = np.gradient(g_thth, r)

    g_tt_s = np.where(np.abs(g_tt) > 1e-50, g_tt, np.nan)
    g_rr_s = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
    g_thth_s = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

    git = 1.0 / g_tt_s
    gir = 1.0 / g_rr_s
    gith = 1.0 / g_thth_s

    Gt_tr = 0.5 * git * g_tt_p
    Gr_rr = 0.5 * gir * g_rr_p
    Gth_rth = 0.5 * gith * g_thth_p

    dGt_tr = np.gradient(Gt_tr, r)
    dGth_rth = np.gradient(Gth_rth, r)

    R_trtr = g_tt_s * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr)
    R_rthrth = g_thth * (dGth_rth + Gth_rth ** 2 - Gth_rth * Gr_rr)
    R_tthtth = -0.25 * g_tt_p * g_thth_p / g_rr_s
    R_thphthph = g_thth * (1.0 - gir * g_thth_p ** 2 / (4.0 * g_thth_s))

    E = R_trtr / (g_tt_s * g_rr_s)
    F_t = R_tthtth / (g_tt_s * g_thth_s)
    F_r = R_rthrth / (g_rr_s * g_thth_s)
    G_c = R_thphthph / g_thth_s ** 2

    return 4.0 * E ** 2 + 8.0 * F_t ** 2 + 8.0 * F_r ** 2 + 4.0 * G_c ** 2


def main():
    print("=" * 70)
    print("STAGE C v4: SPLIT INTEGRATION (Radau solver)")
    print("=" * 70)

    for eta in [0.5, 1.0, 2.0, 5.0]:
        for Q_s in [0.05, 0.1, 0.3]:
            print(f"\n--- η={eta}, Q_s={Q_s} ---")

            # Exterior
            sol_ext = solve_exterior(eta, M=1.0, Q_s=Q_s, r_max=30.0, r_h_plus=2.05)
            print(f"  Exterior: {len(sol_ext.t)} pts, success={sol_ext.success}")

            if not sol_ext.success or len(sol_ext.t) < 5:
                print("  → Exterior failed")
                continue

            # Match at horizon (use last exterior point)
            y_horizon = sol_ext.y[:, -1].copy()
            print(f"  At r=2.05: m={y_horizon[0]:.6f}, φ={y_horizon[2]:.6f}")

            # Interior
            sol_int = solve_interior(eta, y_horizon, r_h_minus=1.95, r_min=1e-3)
            print(f"  Interior: {len(sol_int.t)} pts, success={sol_int.success}")

            if not sol_int.success or len(sol_int.t) < 5:
                print("  → Interior failed")
                continue

            # Combine (skip horizon gap)
            r_all = np.concatenate([sol_int.t[::-1], sol_ext.t[::-1]])
            m_all = np.concatenate([sol_int.y[0][::-1], sol_ext.y[0][::-1]])
            Phi_all = np.concatenate([sol_int.y[1][::-1], sol_ext.y[1][::-1]])
            phi_all = np.concatenate([sol_int.y[2][::-1], sol_ext.y[2][::-1]])

            K = compute_K(r_all, m_all, Phi_all)

            # Report
            print(f"\n  {'r':>10} | {'m':>12} | {'f':>12} | {'φ':>12} | {'K':>14}")
            print(f"  {'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*14}")
            for r_t in [1e-3, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0]:
                idx = np.argmin(np.abs(r_all - r_t))
                f_val = 1.0 - 2.0 * m_all[idx] / r_all[idx]
                K_v = K[idx] if np.isfinite(K[idx]) else float('nan')
                print(f"  {r_all[idx]:>10.4e} | {m_all[idx]:>12.6e} | {f_val:>12.6e} | "
                      f"{phi_all[idx]:>12.6e} | {K_v:>14.6e}")

            # Power law
            mask = (r_all >= 1e-3) & (r_all <= 0.5) & np.isfinite(K) & (K > 0)
            if np.sum(mask) > 10:
                alpha = np.polyfit(np.log(r_all[mask]), np.log(K[mask]), 1)[0]
                m_min = m_all[0]
                print(f"\n  K power law: ~r^{{{alpha:.3f}}}")
                print(f"  m(r_min) = {m_min:.6e}")
                print(f"  → {'REGULAR' if alpha > -1 and abs(m_min) < 0.5 else 'check'}")
            else:
                print(f"\n  Insufficient data for power law")


if __name__ == '__main__':
    main()
