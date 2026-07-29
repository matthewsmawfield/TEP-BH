#!/usr/bin/env python3
"""CRITICAL: Analytical Kretschmann for general φ₀ — the TRUE threshold.

The analytical computation for φ₀=1 shows K = 9/(M²r²) → ∞.
This contradicts the conformal formula A^{-12}*K_Schw ~ r^6 → 0.

The issue: the conformal formula only includes the LEADING Schwarzschild term.
The DERIVATIVE TERMS from the non-constant conformal factor DOMINATE for φ₀ < 3/2.

Analytical result (deep interior, B=0, standard coords):
  All three Kretschmann terms scale as r^{4φ₀-6}
  K_total ~ r^{4φ₀-6}

  φ₀ > 3/2: K → 0 (singularity removed) ✓
  φ₀ = 3/2: K → const
  φ₀ < 3/2: K → ∞ (singularity NOT removed) ✗

This means:
  - Finite areal radius (φ₀ ≤ 1) → K DIVERGES
  - Vanishing K (φ₀ > 3/2) → areal radius DIVERGES
  - MUTUALLY EXCLUSIVE

Verify this with SymPy for several φ₀ values.
"""

import sympy as sp

def compute_exact_kretschmann(phi_0_val):
    """Compute exact Kretschmann for specific φ₀ in deep interior."""

    r, M = sp.symbols('r M', positive=True)
    r_h = 2 * M

    # Deep interior: F ≈ -2M/r, A = (r_h/r)^{φ₀}
    A = (r_h / r)**phi_0_val
    A2 = A**2

    # Standard coordinates metric (B=0)
    g_tt = -A2 * (-2*M/r)  # = A² * 2M/r
    g_rr = A2 / (-2*M/r)   # = -A² * r/(2M)
    g_thth = A2 * r**2      # = r_h^{2φ₀} * r^{2-2φ₀}

    # Simplify
    g_tt = sp.simplify(g_tt)
    g_rr = sp.simplify(g_rr)
    g_thth = sp.simplify(g_thth)

    # Inverse metric
    ginv_tt = sp.simplify(1 / g_tt)
    ginv_rr = sp.simplify(1 / g_rr)
    ginv_thth = sp.simplify(1 / g_thth)

    # Derivatives
    g_tt_p = sp.diff(g_tt, r)
    g_rr_p = sp.diff(g_rr, r)
    g_thth_p = sp.diff(g_thth, r)

    # Christoffel symbols (diagonal metric, only r-derivatives)
    Gt_tr = sp.Rational(1, 2) * ginv_tt * g_tt_p
    Gr_tt = -sp.Rational(1, 2) * ginv_rr * g_tt_p
    Gr_rr = sp.Rational(1, 2) * ginv_rr * g_rr_p
    Gr_thth = -sp.Rational(1, 2) * ginv_rr * g_thth_p
    Gth_rth = sp.Rational(1, 2) * ginv_thth * g_thth_p

    # Riemann components
    # R^t_{rtr} = ∂_r(Γ^t_tr) + (Γ^t_tr)² - Γ^t_tr * Γ^r_rr
    dGt_tr = sp.diff(Gt_tr, r)
    R_trtr_upper = dGt_tr + Gt_tr**2 - Gt_tr * Gr_rr
    R_trtr = g_tt * R_trtr_upper

    # R^θ_{rθr} = ∂_r(Γ^θ_rθ) + (Γ^θ_rθ)² - Γ^θ_rθ * Γ^r_rr
    dGth_rth = sp.diff(Gth_rth, r)
    R_rthrth_upper = dGth_rth + Gth_rth**2 - Gth_rth * Gr_rr
    R_rthrth = g_thth * R_rthrth_upper

    # R^θ_{tθt} = -Γ^θ_{θr} * Γ^r_{tt}
    R_tthtth_upper = -Gth_rth * Gr_tt
    R_tthtth = g_thth * R_tthtth_upper

    # Kretschmann
    K = (4 * R_trtr**2 / (g_tt * g_rr)**2
         + 4 * R_rthrth**2 / (g_rr * g_thth)**2
         + 4 * R_tthtth**2 / (g_tt * g_thth)**2)

    K_simplified = sp.simplify(K)
    K_power = sp.simplify(K_simplified * r**(6 - 4*phi_0_val))  # should give const

    # Limit
    K_limit = sp.limit(K_simplified, r, 0)

    # Areal radius
    areal = A * r
    areal_limit = sp.limit(areal, r, 0)

    return {
        'phi_0': phi_0_val,
        'K': K_simplified,
        'K_power_law': f"r^{{{4*phi_0_val - 6}}}",
        'K_limit': K_limit,
        'areal_limit': areal_limit,
        'R_trtr': sp.simplify(R_trtr),
        'R_rthrth': sp.simplify(R_rthrth),
        'R_tthtth': sp.simplify(R_tthtth),
    }


def main():
    print("=" * 80)
    print("EXACT ANALYTICAL KRETSCHMANN FOR PURE CONFORMAL METRIC")
    print("g̃ = A² * g_Schw, A = (r_h/r)^{φ₀}, deep interior (B=0)")
    print("=" * 80)
    print()
    print("Formula: K ~ r^{4φ₀-6}")
    print("  φ₀ > 3/2: K → 0 (singularity removed)")
    print("  φ₀ = 3/2: K → const")
    print("  φ₀ < 3/2: K → ∞ (singularity NOT removed)")
    print()
    print("Conformal formula A^{-12}*K_Schw ~ r^{12φ₀-6} is WRONG for φ₀ < 3/2")
    print("  It only includes the leading Schwarzschild term, not derivative terms")
    print()

    results = []
    for phi_0_val in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]:
        res = compute_exact_kretschmann(phi_0_val)
        results.append(res)

        print(f"\n--- φ₀ = {phi_0_val} ---")
        print(f"  Areal radius → {res['areal_limit']} "
              f"{'(DIVERGES)' if res['areal_limit'] == sp.oo else '(finite)'}")
        print(f"  R_trtr = {res['R_trtr']}")
        print(f"  R_rθrθ = {res['R_rthrth']}")
        print(f"  R_tθtθ = {res['R_tthtth']}")
        print(f"  Kretschmann = {res['K']}")
        print(f"  K power law: {res['K_power_law']}")
        print(f"  K limit: {res['K_limit']}")
        print(f"  K → 0: {'YES ✓' if res['K_limit'] == 0 else 'NO ✗' if res['K_limit'] == sp.oo else 'CONST'}")

    # Summary table
    print("\n" + "=" * 80)
    print("SUMMARY: TRUE KRETSCHMANN BEHAVIOR vs CONFORMAL FORMULA")
    print("=" * 80)
    print(f"  {'φ₀':>6} | {'Areal':>12} | {'True K':>16} | {'Conf. K':>16} | {'K→0?':>8} | {'Areal finite?':>14}")
    print(f"  {'-'*6}-+-{'-'*12}-+-{'-'*16}-+-{'-'*16}-+-{'-'*8}-+-{'-'*14}")

    for res in results:
        phi_0 = res['phi_0']
        areal_str = "diverges" if res['areal_limit'] == sp.oo else "finite"
        true_K_power = f"r^{{{4*phi_0-6:.1f}}}"
        conf_K_power = f"r^{{{12*phi_0-6:.1f}}}"
        K_vanishes = res['K_limit'] == 0
        areal_finite = res['areal_limit'] != sp.oo

        print(f"  {phi_0:>6.2f} | {areal_str:>12} | {true_K_power:>16} | {conf_K_power:>16} | "
              f"{'YES' if K_vanishes else 'NO':>8} | {'YES' if areal_finite else 'NO':>14}")

    print()
    print("  TRUE threshold for K→0: φ₀ > 3/2 = 1.5")
    print("  CONFORMAL FORMULA threshold: φ₀ > 1/2 = 0.5 (WRONG!)")
    print()
    print("  FINITE AREAL RADIUS: φ₀ ≤ 1")
    print("  VANISHING K:         φ₀ > 3/2")
    print()
    print("  ★ THESE ARE MUTUALLY EXCLUSIVE ★")
    print("  No value of φ₀ gives both finite areal radius AND vanishing Kretschmann")
    print()
    print("  IMPLICATION:")
    print("  - The conformal prototype (φ₀=2) removes the singularity (K→0)")
    print("    but has diverging areal radius")
    print("  - The φ₀=1 branch has finite areal radius but K DIVERGES")
    print("  - The manuscript's conformal formula A^{-12}*K_Schw is incorrect")
    print("    for φ₀ < 3/2 — it underestimates the curvature")
    print("  - The pipeline has been using the wrong formula for curvature")


if __name__ == '__main__':
    main()
