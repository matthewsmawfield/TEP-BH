#!/usr/bin/env python3
"""CRITICAL CHECK: Does the Ricci scalar diverge for φ₀=1?

The conformal formula K = A^{-12} * K_Schw is only the LEADING TERM.
The full Kretschmann for g̃ = A² * g with non-constant A includes
derivative terms that may dominate.

Analytical result for R[g̃] with A = (r_h/r)^{φ₀}:
  R[g̃] = A^{-2} * (-6Fφ₀(1+φ₀)/r²) ~ r^{2φ₀-3}

  φ₀=1: R ~ r^{-1} → ∞  (DIVERGES!)
  φ₀=1.5: R ~ r^0 → const (finite)
  φ₀=2: R ~ r^1 → 0 (vanishes)

If R → ∞, then K cannot → 0 (for conformally flat: K = 2R_{μν}R^{μν} - R²/3).
This means the conformal formula A^{-12}*K_Schw is WRONG for non-constant A.

Let me compute the curvature in STANDARD Schwarzschild coordinates (t,r)
where g_rr = 1/F (doesn't vanish), avoiding the EF g_rr → 0 issue.
"""

import numpy as np

def compute_metric_standard_coords(r, params):
    """Compute matter metric in STANDARD Schwarzschild coordinates (t, r).

    Geometric: ds_g² = -F dt² + F^{-1} dr² + r² dΩ²
    Matter: g̃ = A² g + B ∇φ ∇φ, with φ = qt + ψ(r)

    Components:
      g̃_tt = -A²F + Bq²
      g̃_tr = Bqψ'  (mixed term from disformal)
      g̃_rr = A²/F + Bψ'²
      g̃_θθ = A²r²
    """
    r = np.asarray(r, dtype=float)
    M = params['M']
    q = params.get('q', 0.1)
    phi_0 = params.get('phi_0', 1.0)
    delta = params.get('delta', 0.05)
    B0 = params.get('B0', 1.0)
    sigma_B = params.get('sigma_B', 1.5)

    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)
    psi = phi_0 * np.log(r_safe / r_h) * S
    dpsi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    A = np.exp(-psi)  # beta_A = -1
    A2 = A ** 2

    abs_psi = np.abs(psi)
    B = B0 * abs_psi ** 2 / (1.0 + abs_psi ** 2) * np.exp(-(psi ** 4) / (2.0 * sigma_B ** 4))

    # Standard coordinates metric
    gtt = -A2 * F + B * q ** 2
    gtr = B * q * dpsi
    grr = A2 / F + B * dpsi ** 2  # NOTE: A²/F, not 0 like in EF!
    gthth = A2 * r ** 2

    # 2D (t,r) determinant
    det_2d = gtt * grr - gtr ** 2

    return {
        'r': r, 'F': F, 'psi': psi, 'dpsi': dpsi,
        'A': A, 'B': B, 'A2': A2,
        'gtt': gtt, 'gtr': gtr, 'grr': grr, 'gthth': gthth,
        'det_2d': det_2d,
    }


def compute_curvature_standard(r, metric, M=1.0):
    """Compute curvature in standard Schwarzschild coordinates.

    For spherically symmetric metric:
      ds² = g_tt dt² + 2 g_tr dt dr + g_rr dr² + g_θθ dΩ²

    Here g_rr = A²/F + Bψ'², which does NOT vanish in the deep interior
    (unlike EF where g_rr = Bψ'² → 0). This avoids the 0/0 numerical issue.
    """
    gtt = metric['gtt']
    gtr = metric['gtr']
    grr = metric['grr']
    gthth = metric['gthth']
    det_2d = metric['det_2d']

    # Derivatives
    gtt_p = np.gradient(gtt, r)
    gtr_p = np.gradient(gtr, r)
    grr_p = np.gradient(grr, r)
    gthth_p = np.gradient(gthth, r)

    gtt_pp = np.gradient(gtt_p, r)
    grr_pp = np.gradient(grr_p, r)
    gthth_pp = np.gradient(gthth_p, r)

    # Safe versions
    det_safe = np.where(np.abs(det_2d) > 1e-100, det_2d, np.nan)
    gthth_safe = np.where(np.abs(gthth) > 1e-100, gthth, np.nan)
    grr_safe = np.where(np.abs(grr) > 1e-100, grr, np.nan)

    # Inverse metric (2D block)
    ginv_tt = grr / det_safe
    ginv_tr = -gtr / det_safe
    ginv_rr = gtt / det_safe
    ginv_thth = 1.0 / gthth_safe

    # Christoffel symbols (only r-derivatives for static metric)
    # Γ^λ_μν = (1/2) g^{λσ} (∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
    # For static spherically symmetric, only ∂_r is nonzero

    Gt_tt = -0.5 * ginv_tr * gtt_p  # Γ^t_tt
    Gr_tt = -0.5 * ginv_rr * gtt_p  # Γ^r_tt
    Gt_tr = 0.5 * ginv_tt * gtt_p   # Γ^t_tr
    Gr_tr = 0.5 * ginv_tr * gtt_p   # Γ^r_tr
    Gt_rr = ginv_tt * gtr_p + 0.5 * ginv_tr * grr_p  # Γ^t_rr
    Gr_rr = ginv_tr * gtr_p + 0.5 * ginv_rr * grr_p  # Γ^r_rr
    Gt_thth = -0.5 * ginv_tr * gthth_p  # Γ^t_θθ
    Gr_thth = -0.5 * ginv_rr * gthth_p  # Γ^r_θθ
    Gth_rth = 0.5 * ginv_thth * gthth_p  # Γ^θ_rθ

    # Riemann component R_{trtr}
    R_trtr = -0.5 * gtt_pp + (
        gtt * (Gt_tt * Gt_rr - Gt_tr ** 2)
        + gtr * (Gt_tt * Gr_rr - Gt_tr * Gr_tr)
        + gtr * (Gr_tt * Gt_rr - Gr_tr * Gt_tr)
        + grr * (Gr_tt * Gr_rr - Gr_tr ** 2)
    )

    # Angular Riemann R_{rθrθ}
    R_rthrth = -0.5 * gthth_pp + gthth_p * Gr_thth

    # Ricci scalar
    R_2D = 2.0 * R_trtr / det_safe
    R_ang = -gthth_pp / gthth_safe + 0.5 * (gthth_p / gthth_safe) ** 2
    R_total = R_2D + 2.0 * R_ang / gthth_safe

    # Kretschmann
    K_term1 = 4.0 * R_trtr ** 2 / det_safe ** 2
    K_term2 = 4.0 * R_rthrth ** 2 / (gthth_safe * grr_safe)
    Kretschmann = K_term1 + K_term2

    # Ricci squared
    Ricci_tt = R_trtr * ginv_rr
    Ricci_rr = R_trtr * ginv_tt
    Ricci_tr = -R_trtr * ginv_tr
    Ricci_thth = R_ang
    Ricci_sq = (Ricci_tt ** 2 * ginv_tt + 2 * Ricci_tr ** 2 * ginv_tr
                + Ricci_rr ** 2 * ginv_rr + 2 * Ricci_thth ** 2 * ginv_thth)

    return {
        'Ricci_scalar': R_total,
        'Ricci_squared': Ricci_sq,
        'Kretschmann': Kretschmann,
        'R_trtr': R_trtr,
        'R_rthrth': R_rthrth,
    }


def test_curvature_standard_coords():
    """Test curvature in standard coordinates for various φ₀."""
    r = np.logspace(np.log10(1e-8), np.log10(50.0), 50000)

    print("=" * 80)
    print("CURVATURE IN STANDARD SCHWARZSCHILD COORDINATES")
    print("(Avoids EF g_rr → 0 numerical issue)")
    print("=" * 80)

    for phi_0 in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
        params = {'M': 1.0, 'q': 0.1, 'phi_0': phi_0, 'delta': 0.05,
                  'B0': 1.0, 'sigma_B': 1.5}
        metric = compute_metric_standard_coords(r, params)
        curvature = compute_curvature_standard(r, metric, M=1.0)

        K = curvature['Kretschmann']
        R = curvature['Ricci_scalar']
        R_sq = curvature['Ricci_squared']

        # Power laws near r=0
        n_check = 200
        r_inner = r[:n_check]

        results = {}
        for name, vals in [("K", K), ("R", R), ("R_sq", R_sq)]:
            mask = (np.abs(vals[:n_check]) > 0) & (r_inner > 0) & np.isfinite(vals[:n_check])
            if np.sum(mask) > 5:
                log_r = np.log(r_inner[mask])
                log_val = np.log(np.abs(vals[:n_check][mask]))
                alpha = np.polyfit(log_r, log_val, 1)[0]
                results[name] = alpha
            else:
                results[name] = np.nan

        # Areal radius
        areal = metric['A'] * r
        areal_inner = areal[0]
        areal_div = areal_inner > 1e3

        # Analytical predictions
        R_pred = 2 * phi_0 - 3  # R ~ r^{2φ₀-3}
        K_pred = 12 * phi_0 - 6  # Leading K term ~ r^{12φ₀-6} (but may not be dominant)

        print(f"\n  φ₀ = {phi_0:.2f}:")
        print(f"    Areal radius at r→0: {areal_inner:.4e} {'(DIVERGES)' if areal_div else '(finite)'}")
        print(f"    K power law: ~r^{{{results['K']:.2f}}} (leading term pred: r^{{{K_pred:.1f}}})")
        print(f"    R power law: ~r^{{{results['R']:.2f}}} (analytical pred: r^{{{R_pred:.1f}}})")
        print(f"    R² power law: ~r^{{{results['R_sq']:.2f}}}")

        # Verdict
        K_finite = results['K'] > 0 or (np.isfinite(results['K']) and results['K'] >= 0)
        R_finite = results['R'] >= 0
        areal_finite = not areal_div

        K_status = "VANISHES ✓" if results['K'] > 0 else ("DIVERGES ✗" if results['K'] < -1 else "finite")
        R_status = "VANISHES ✓" if results['R'] > 0 else ("DIVERGES ✗" if results['R'] < -1 else "finite")

        print(f"    K: {K_status}")
        print(f"    R: {R_status}")
        print(f"    Areal: {'finite ✓' if areal_finite else 'DIVERGES ✗'}")

        # Values at specific radii
        for r_test in [1.0, 0.1, 0.01, 0.001, 1e-6]:
            idx = np.argmin(np.abs(r - r_test))
            K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
            R_val = R[idx] if np.isfinite(R[idx]) else float('nan')
            print(f"      r={r_test:.0e}: K={K_val:.4e}, R={R_val:.4e}")

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY: Curvature power laws (standard coordinates)")
    print(f"{'='*80}")
    print(f"  {'φ₀':>6} | {'Areal':>12} | {'K':>12} | {'R':>12} | {'R²':>12} | {'All finite?':>12}")
    print(f"  {'-'*6}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}")

    for phi_0 in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
        params = {'M': 1.0, 'q': 0.1, 'phi_0': phi_0, 'delta': 0.05,
                  'B0': 1.0, 'sigma_B': 1.5}
        metric = compute_metric_standard_coords(r, params)
        curvature = compute_curvature_standard(r, metric, M=1.0)

        K = curvature['Kretschmann']
        R = curvature['Ricci_scalar']
        R_sq = curvature['Ricci_squared']
        areal = metric['A'] * r

        n_check = 200
        r_inner = r[:n_check]
        powers = {}
        for name, vals in [("K", K), ("R", R), ("R_sq", R_sq)]:
            mask = (np.abs(vals[:n_check]) > 0) & (r_inner > 0) & np.isfinite(vals[:n_check])
            if np.sum(mask) > 5:
                powers[name] = np.polyfit(np.log(r_inner[mask]), np.log(np.abs(vals[:n_check][mask])), 1)[0]
            else:
                powers[name] = np.nan

        areal_div = areal[0] > 1e3
        all_finite = (powers['K'] >= 0) and (powers['R'] >= 0) and not areal_div

        areal_str = "diverges" if areal_div else "finite"
        print(f"  {phi_0:>6.2f} | {areal_str:>12} | r^{{{powers['K']:>6.2f}}} | r^{{{powers['R']:>6.2f}}} | r^{{{powers['R_sq']:>6.2f}}} | {'YES' if all_finite else 'NO':>12}")


if __name__ == '__main__':
    test_curvature_standard_coords()
