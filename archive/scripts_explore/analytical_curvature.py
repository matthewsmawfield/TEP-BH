#!/usr/bin/env python3
"""Analytical verification of curvature for conformal metric g̃ = A² * g_Schw.

The numerical computations have been giving contradictory results:
- EF coordinates with blending: K ~ r^6 → 0 (vanishes) for φ₀=1
- Standard coordinates: K ~ r^{-1.6} → ∞ (diverges) for φ₀=1, with NEGATIVE values

Negative Kretschmann is impossible (K = sum of squares ≥ 0), so the standard
coordinates computation has a numerical error from the F→0 singularity at r=2M.

This script computes the curvature ANALYTICALLY for the pure conformal case
(B=0) to determine the true behavior.

For g̃ = Ω² * g in 4D, the Kretschmann transforms as:
  K[g̃] = Ω^{-12} * K[g] + derivative terms

The derivative terms involve:
  - (∇∇Ω)² terms: scale as Ω^{-12} * Ω² * (d²Ω/dr²)² ~ Ω^{-10} * (Ω/r²)²
  - (∇Ω)⁴ terms: scale as Ω^{-12} * (Ω/r)⁴ ~ Ω^{-8} / r⁴
  - Cross terms with Riemann: vanish for Schwarzschild (Weyl = trace-free Riemann)

For Ω = A = (r_h/r)^{φ₀} = C * r^{-φ₀}:
  Ω' = -φ₀ * Ω / r
  Ω'' = φ₀(φ₀+1) * Ω / r²

The derivative terms in K[g̃]:
  T1 = Ω^{-12} * K[g] ~ r^{12φ₀} * r^{-6} = r^{12φ₀-6}  (leading Schwarzschild)
  T2 = Ω^{-12} * (Ω'')² * (metric) ~ r^{12φ₀} * (r^{-φ₀-2})² * r² = r^{10φ₀-2}
       Wait, need to be more careful about metric contractions.

Let me compute this properly using the known formula for conformal Kretschmann.

For a conformal transformation g̃_{μν} = Ω² g_{μν} in n dimensions:
  K[g̃] = Ω^{-2n+2} * K[g] + ... (derivative terms)

In 4D (n=4): K[g̃] = Ω^{-6} * K[g] + ...

Wait, that's not right either. Let me look up the exact formula.

For g̃ = Ω² g in 4D:
  R̃_{μνρσ} = Ω² R_{μνρσ} - 2(δ_{[μ}^{[α} δ_{ν]}^{β]} g_{ρσ}) ∇_α∇_β Ω/Ω + ... (complex)

The Kretschmann is:
  K[g̃] = R̃_{μνρσ} R̃^{μνρσ}

For Schwarzschild (R_{μν} = 0, R = 0), the conformal Riemann is:
  R̃_{μνρσ} = Ω² [R_{μνρσ} - 2 g_{ρ[σ} ∇_{μ]} ∇_{ν]} (ln Ω) - ... ]

This is getting very complex. Let me use a different approach: compute the
Riemann tensor components directly for the specific metric.

For the spherically symmetric metric in EF coordinates:
  ds² = g_vv dv² + 2 g_vr dv dr + g_rr dr² + g_θθ dΩ²

With B=0 (pure conformal):
  g_vv = -A²F, g_vr = A², g_rr = 0, g_θθ = A²r²

The non-zero Riemann components for spherical symmetry are:
  R_{vrvr}, R_{vrθθ}, R_{θφθφ}, R_{rθrθ}

And K = 4 R_{vrvr}² / (det_2d)² + 4 R_{rθrθ}² / (g_θθ * g_rr) + ...

The problem: g_rr = 0, so the second term is 0/0.

But we can compute R_{rθrθ} directly:
R_{rθrθ} = -0.5 * g_θθ'' + g_θθ' * Γ^r_{θθ}

Γ^r_{θθ} = -0.5 * g^{rr} * g_θθ'

With g^{rr} = g_vv / det_2d = (-A²F) / (-A⁴) = F/A²

So Γ^r_{θθ} = -0.5 * (F/A²) * g_θθ'

For φ₀=1: g_θθ = A²r² = r_h² (constant), so g_θθ' = 0, g_θθ'' = 0
→ R_{rθrθ} = 0

So the angular Kretschmann term is 0² / 0 = 0/0, but the NUMERATOR is 0.
The limit is 0 (since R_{rθrθ} → 0 faster than sqrt(g_rr) → 0).

For the (v,r) component:
R_{vrvr} = -0.5 * g_vv'' + (Christoffel products)

g_vv = -A²F = -(r_h/r)² * (1-2M/r) = -r_h²/r² + 2Mr_h²/r³

g_vv' = 2r_h²/r³ - 6Mr_h²/r⁴
g_vv'' = -6r_h²/r⁴ + 24Mr_h²/r⁵

For the Christoffel products, we need:
Γ^v_vv, Γ^r_vv, Γ^v_vr, Γ^r_vr, Γ^v_rr, Γ^r_rr

With g^{vv} = g_rr/det = 0/(-A⁴) = 0
g^{vr} = -g_vr/det = -A²/(-A⁴) = 1/A²
g^{rr} = g_vv/det = (-A²F)/(-A⁴) = F/A²

Γ^v_vv = -0.5 * g^{vr} * g_vv' = -0.5 * (1/A²) * g_vv'
Γ^r_vv = -0.5 * g^{rr} * g_vv' = -0.5 * (F/A²) * g_vv'
Γ^v_vr = 0.5 * g^{vv} * g_vv' = 0 (since g^{vv} = 0)
Γ^r_vr = 0.5 * g^{vr} * g_vv' = 0.5 * (1/A²) * g_vv'
Γ^v_rr = g^{vv} * g_vr' + 0.5 * g^{vr} * g_rr' = 0 + 0.5 * (1/A²) * 0 = 0
         (g_vr = A², g_vr' = 2AA', g_rr = 0, g_rr' = 0)
         Wait: g_vr' = (A²)' = 2A*A' = 2A*(-A/r) = -2A²/r
         So Γ^v_rr = 0 * (-2A²/r) + 0.5 * (1/A²) * 0 = 0
Γ^r_rr = g^{vr} * g_vr' + 0.5 * g^{rr} * g_rr' = (1/A²)*(-2A²/r) + 0 = -2/r

Now R_{vrvr}:
R_{vrvr} = -0.5 * g_vv'' + g_vv*(Γ^v_vv*Γ^r_rr - Γ^v_vr²) + g_vr*(Γ^v_vv*Γ^r_vr - ...) + ...

This is getting very tedious. Let me just compute it numerically with HIGH PRECISION
using a proper coordinate system.

The best approach: use TORTOISE coordinate r* and double-null coordinates (u,v)
which are regular everywhere.

Or even simpler: compute in EF coordinates but use ANALYTICAL derivatives
instead of numerical finite differences.
"""

import numpy as np
import sympy as sp

def analytical_curvature_phi0_1():
    """Compute curvature analytically using SymPy for φ₀=1, B=0."""

    r, M, r_h = sp.symbols('r M r_h', positive=True)

    # Conformal factor: A = r_h/r (for φ₀=1, deep interior S=1)
    A = r_h / r
    A2 = A**2

    # Schwarzschild function
    F = 1 - 2*M/r

    # EF metric components (B=0, pure conformal)
    g_vv = -A2 * F
    g_vr = A2
    g_rr = sp.Integer(0)  # B=0
    g_thth = A2 * r**2  # = r_h² (constant for φ₀=1!)

    # 2D determinant
    det_2d = g_vv * g_rr - g_vr**2
    det_2d_simplified = sp.simplify(det_2d)
    print(f"det_2d = {det_2d_simplified}")

    # Inverse metric (2D block)
    # For [[a, b], [b, c]]: inv = (1/det) * [[c, -b], [-b, a]]
    # g^{vv} = g_rr/det = 0
    # g^{vr} = -g_vr/det
    # g^{rr} = g_vv/det
    ginv_vv = g_rr / det_2d
    ginv_vr = -g_vr / det_2d
    ginv_rr = g_vv / det_2d
    ginv_thth = 1 / g_thth

    ginv_vv = sp.simplify(ginv_vv)
    ginv_vr = sp.simplify(ginv_vr)
    ginv_rr = sp.simplify(ginv_rr)

    print(f"\ng^vv = {ginv_vv}")
    print(f"g^vr = {ginv_vr}")
    print(f"g^rr = {ginv_rr}")
    print(f"g^thth = {sp.simplify(ginv_thth)}")

    # First derivatives
    g_vv_p = sp.diff(g_vv, r)
    g_vr_p = sp.diff(g_vr, r)
    g_thth_p = sp.diff(g_thth, r)

    # Second derivatives
    g_vv_pp = sp.diff(g_vv_p, r)
    g_thth_pp = sp.diff(g_thth_p, r)

    print(f"\ng_vv = {sp.simplify(g_vv)}")
    print(f"g_vv' = {sp.simplify(g_vv_p)}")
    print(f"g_vv'' = {sp.simplify(g_vv_pp)}")
    print(f"g_thth = {sp.simplify(g_thth)}")
    print(f"g_thth' = {sp.simplify(g_thth_p)}")
    print(f"g_thth'' = {sp.simplify(g_thth_pp)}")

    # Christoffel symbols
    Gv_vv = -sp.Rational(1, 2) * ginv_vr * g_vv_p
    Gr_vv = -sp.Rational(1, 2) * ginv_rr * g_vv_p
    Gv_vr = sp.Rational(1, 2) * ginv_vv * g_vv_p  # = 0 since g^{vv}=0
    Gr_vr = sp.Rational(1, 2) * ginv_vr * g_vv_p
    Gv_rr = ginv_vv * g_vr_p + sp.Rational(1, 2) * ginv_vr * 0  # g_rr'=0
    Gr_rr = ginv_vr * g_vr_p + sp.Rational(1, 2) * ginv_rr * 0
    Gv_thth = -sp.Rational(1, 2) * ginv_vr * g_thth_p
    Gr_thth = -sp.Rational(1, 2) * ginv_rr * g_thth_p
    Gth_rth = sp.Rational(1, 2) * ginv_thth * g_thth_p

    # Simplify
    for name, val in [("Γ^v_vv", Gv_vv), ("Γ^r_vv", Gr_vv),
                       ("Γ^v_vr", Gv_vr), ("Γ^r_vr", Gr_vr),
                       ("Γ^v_rr", Gv_rr), ("Γ^r_rr", Gr_rr),
                       ("Γ^v_θθ", Gv_thth), ("Γ^r_θθ", Gr_thth),
                       ("Γ^θ_rθ", Gth_rth)]:
        print(f"  {name} = {sp.simplify(val)}")

    # Riemann component R_{vrvr}
    R_vrvr = (-sp.Rational(1, 2) * g_vv_pp
              + g_vv * (Gv_vv * Gr_rr - Gv_vr**2)
              + g_vr * (Gv_vv * Gr_vr - Gv_vr * Gr_vr)  # Note: Gr_vr = Gv_vr here? No
              + g_vr * (Gr_vv * Gv_rr - Gr_vr * Gv_vr)
              + g_rr * (Gr_vv * Gr_rr - Gr_vr**2))

    # Wait, I need to be more careful. The Riemann tensor formula for
    # a 2D metric with off-diagonal terms is:
    # R_{abcd} = ... (complex)
    # For our case with g_rr = 0, the formula simplifies.

    # Actually, for the metric ds² = g_vv dv² + 2 g_vr dv dr + g_θθ dΩ²
    # with g_rr = 0, the only nonzero Riemann components are:
    # R_{vrvr}, R_{vθvθ}, R_{rθrθ}, R_{θφθφ}

    # R_{vrvr} = -0.5 * g_vv'' + (products of Christoffels)
    # But with g_rr = 0, many terms simplify.

    # Let me use the standard formula:
    # R_{vrvr} = ∂_r Γ_{vv,r} - ∂_r Γ_{vr,v} + Γ^λ_{vv} Γ^r_{λr} - Γ^λ_{vr} Γ^r_{λv}
    # Wait, this isn't right either. Let me use:
    # R_{abcd} = g_{ae} R^e_{bcd}
    # R^a_{bcd} = ∂_c Γ^a_{bd} - ∂_d Γ^a_{bc} + Γ^a_{ce} Γ^e_{bd} - Γ^a_{de} Γ^e_{bc}

    # For R_{vrvr} = g_{vv} R^v_{rvr} + g_{vr} R^r_{rvr}
    # R^v_{rvr} = ∂_v Γ^v_{rr} - ∂_r Γ^v_{rv} + Γ^v_{vλ} Γ^λ_{rr} - Γ^v_{rλ} Γ^λ_{rv}
    # For static metric: ∂_v = 0
    # = -∂_r Γ^v_{rv} + Γ^v_{vv} Γ^v_{rr} + Γ^v_{vr} Γ^r_{rr} - Γ^v_{rv} Γ^v_{rv} - Γ^v_{rr} Γ^r_{rv}

    # This is getting very complex with SymPy. Let me try a different approach:
    # compute the curvature using the warped product decomposition.

    # For ds² = h_{AB} dx^A dx^B + R² dΩ² where A,B ∈ {v,r}:
    # R_{ABCD} = (1/2)(∂_B∂_D h_{AC} + ∂_A∂_C h_{BD} - ∂_B∂_C h_{AD} - ∂_A∂_D h_{BC}) + ...
    # R_{ArBr} = -R * ∂_A∂_B R / R + ... (warped product)

    # Actually, let me just compute the Kretschmann using the known formula
    # for spherically symmetric metrics.

    # For ds² = -e^{2Φ} dv² + 2 e^{Φ+Λ} dv dr + R² dΩ² (general form)
    # The Kretschmann involves:
    # K = 4(R_{vrvr})²/(det)² + 4(R_{rθrθ})²/(R² * g_rr) + 4(R_{vθvθ})²/(R² * g_vv) + ...

    # This is still complex. Let me just use SymPy to compute the full Riemann.

    # Actually, the simplest approach: compute in STANDARD Schwarzschild coordinates
    # but ONLY in the deep interior (r << 2M) where F ≈ -2M/r, avoiding the horizon.

    print("\n" + "=" * 60)
    print("ANALYTICAL COMPUTATION IN DEEP INTERIOR (r << 2M)")
    print("=" * 60)

    # In deep interior: F ≈ -2M/r, A = r_h/r, A² = r_h²/r²
    # Standard coordinates: ds² = A²(-F dt² + F^{-1} dr² + r² dΩ²)
    # = A² * (2M/r * dt² + r/(2M) * dr² + r² dΩ²)  (inside, F<0)
    # = (r_h²/r²) * (2M/r * dt² + r/(2M) * dr² + r² dΩ²)
    # = r_h² * (2M/r³ * dt² + 1/(2Mr) * dr² + dΩ²)

    # So: g_tt = 2M*r_h²/r³, g_rr = r_h²/(2Mr), g_θθ = r_h²
    # With r_h = 2M: g_tt = 8M³/r³, g_rr = 2M/r, g_θθ = 4M²

    # For a metric ds² = e^{2α(r)} dt² + e^{2β(r)} dr² + R² dΩ²:
    # K = 4[e^{-4β}(α'' + α'² - α'β' + 2α'/r - 2β'/r)²
    #      + e^{-4β}(α'/r - 1/(R²e^{2β}))² * R²
    #      + (1 - e^{-2β} R'²)² / R⁴]

    # Wait, this formula is for specific coordinate choices. Let me use the
    # general spherically symmetric formula.

    # For ds² = -N(r)² dt² + B(r)² dr² + R(r)² dΩ²:
    # (note: inside horizon, N² < 0, so -N² > 0, t is spacelike)

    # The Kretschmann scalar is:
    # K = 4/(N⁴B⁴) * [N''N - N'² + N'N'B'/B]²
    #   + 4/(N²B²R²) * [N'R'/(BR) - N'N/(B²R²)]²  ... (this isn't right)

    # Let me use the EXACT formula from Stephani's textbook:
    # For ds² = e^{2a(r)} dt² + e^{2b(r)} dr² + r̃² dΩ²:
    # K = 4*e^{-4b} * (a'' + a'² - a'b' + 2a'/r̃ - 2b'/r̃)² * (r̃'/r̃)²
    #   ... no, this depends on whether r̃ = r or r̃ = R(r).

    # OK let me just compute it directly with SymPy using the deep interior
    # approximation.

    r_sym = sp.Symbol('r', positive=True)
    M_sym = sp.Symbol('M', positive=True)
    r_h_sym = 2 * M_sym

    # Deep interior: F ≈ -2M/r
    F_approx = -2 * M_sym / r_sym

    # Standard coordinates metric (deep interior, B=0)
    A_approx = r_h_sym / r_sym
    A2_approx = A_approx**2

    g_tt = -A2_approx * F_approx  # = A² * 2M/r = 2M * r_h² / r³ = 8M³/r³
    g_rr = A2_approx / F_approx   # = A² / (-2M/r) = -A² * r/(2M) = -r_h²/(2Mr) = -2M/r
    g_thth = A2_approx * r_sym**2  # = r_h² = 4M²

    print(f"\nDeep interior metric (standard coords):")
    print(f"  g_tt = {sp.simplify(g_tt)} = {sp.simplify(g_tt).subs(M_sym, 1)}")
    print(f"  g_rr = {sp.simplify(g_rr)} = {sp.simplify(g_rr).subs(M_sym, 1)}")
    print(f"  g_θθ = {sp.simplify(g_thth)} = {sp.simplify(g_thth).subs(M_sym, 1)}")

    # Note: g_rr < 0 inside (r is timelike). For curvature, we use |g_rr|.
    # The 2D (t,r) determinant:
    det_2d = g_tt * g_rr  # no off-diagonal in standard coords
    print(f"  det_2d = {sp.simplify(det_2d)}")

    # Inverse metric
    ginv_tt = 1 / g_tt
    ginv_rr = 1 / g_rr
    ginv_thth = 1 / g_thth

    # Christoffel symbols (standard coords, only r-derivatives)
    # Γ^t_tr = 0.5 * g^{tt} * g_tt'
    # Γ^r_tt = -0.5 * g^{rr} * g_tt'
    # Γ^r_rr = 0.5 * g^{rr} * g_rr'
    # Γ^r_θθ = -0.5 * g^{rr} * g_thth'
    # Γ^θ_rθ = 0.5 * g^{θθ} * g_thth'

    g_tt_p = sp.diff(g_tt, r_sym)
    g_rr_p = sp.diff(g_rr, r_sym)
    g_thth_p = sp.diff(g_thth, r_sym)

    g_tt_pp = sp.diff(g_tt_p, r_sym)
    g_rr_pp = sp.diff(g_rr_p, r_sym)
    g_thth_pp = sp.diff(g_thth_p, r_sym)

    Gt_tr = sp.Rational(1, 2) * ginv_tt * g_tt_p
    Gr_tt = -sp.Rational(1, 2) * ginv_rr * g_tt_p
    Gr_rr = sp.Rational(1, 2) * ginv_rr * g_rr_p
    Gr_thth = -sp.Rational(1, 2) * ginv_rr * g_thth_p
    Gth_rth = sp.Rational(1, 2) * ginv_thth * g_thth_p

    print(f"\nChristoffel symbols:")
    for name, val in [("Γ^t_tr", Gt_tr), ("Γ^r_tt", Gr_tt),
                       ("Γ^r_rr", Gr_rr), ("Γ^r_θθ", Gr_thth),
                       ("Γ^θ_rθ", Gth_rth)]:
        simplified = sp.simplify(val)
        print(f"  {name} = {simplified}")

    # Riemann components
    # R_{trtr} = -0.5 * g_tt'' + g_tt * (Γ^t_tt*Γ^r_rr - ...) + ...
    # For standard coords (no off-diagonal):
    # R_{trtr} = -0.5 * g_tt'' + g_tt * Γ^r_tt * Γ^t_tr + g_rr * Γ^t_tr * Γ^r_tt
    # Wait, the general formula for diagonal metric:
    # R_{trtr} = -0.5 * g_tt'' + 0.5 * g_tt' * g_rr'/g_rr + 0.25 * (g_tt'/g_tt)² * g_tt - 0.25 * (g_tt'/g_rr)² * g_rr
    # This isn't right either. Let me use the component formula.

    # R^t_{rtr} = ∂_r Γ^t_rr - ∂_r Γ^t_tr + Γ^t_{rλ}Γ^λ_tr - Γ^t_{tλ}Γ^λ_rr
    # For diagonal metric: Γ^t_rr = 0, Γ^t_tr = 0.5 g^{tt} g_tt'
    # R^t_{rtr} = -∂_r(0.5 g^{tt} g_tt') + Γ^t_{rr}Γ^r_tr + Γ^t_{rt}Γ^t_tr - Γ^t_{tt}Γ^t_rr - Γ^t_{tr}Γ^r_rr
    # = -∂_r(0.5 g^{tt} g_tt') + 0 + 0 - 0 - 0.5 g^{tt} g_tt' * 0.5 g^{rr} g_rr'
    # Hmm, this is getting messy. Let me just use the standard result.

    # For a diagonal metric ds² = g_tt dt² + g_rr dr² + g_θθ dΩ²:
    # R_{trtr} = -0.5 * g_tt'' + (0.25/g_rr) * g_tt' * g_rr' + (0.25/g_tt) * g_tt'²
    # Wait, I need to be more careful.

    # R_{trtr} = g_{tt} R^t_{rtr}
    # R^t_{rtr} = ∂_r Γ^t_{tr} - ∂_t Γ^t_{rr} + Γ^t_{rα} Γ^α_{tr} - Γ^t_{tα} Γ^α_{rr}
    # Static: ∂_t = 0
    # = ∂_r Γ^t_{tr} + Γ^t_{rr} Γ^r_{tr} + Γ^t_{rt} Γ^t_{tr} - Γ^t_{tr} Γ^r_{rr} - Γ^t_{tt} Γ^t_{rr}
    # Γ^t_{rr} = 0 (diagonal), Γ^t_{tt} = 0 (diagonal, ∂_t g = 0)
    # = ∂_r Γ^t_{tr} + 0 + Γ^t_{tr}² - Γ^t_{tr} Γ^r_{rr} - 0
    # = ∂_r(0.5 g^{tt} g_tt') + (0.5 g^{tt} g_tt')² - 0.5 g^{tt} g_tt' * 0.5 g^{rr} g_rr'

    # Let me just compute numerically with SymPy
    Gt_tr_val = sp.Rational(1, 2) * ginv_tt * g_tt_p
    Gr_rr_val = sp.Rational(1, 2) * ginv_rr * g_rr_p

    dGt_tr_dr = sp.diff(Gt_tr_val, r_sym)

    R_trtr_upper = dGt_tr_dr + Gt_tr_val**2 - Gt_tr_val * Gr_rr_val
    R_trtr = g_tt * R_trtr_upper

    R_trtr_simplified = sp.simplify(R_trtr)
    print(f"\nR_trtr = {R_trtr_simplified}")

    # R_{rθrθ}
    # R^θ_{rθr} = ∂_r Γ^θ_{θr} - Γ^θ_{θλ} Γ^λ_{rr}
    # = ∂_r(0.5 g^{θθ} g_θθ') - Γ^θ_{θθ} Γ^θ_{rr} - Γ^θ_{θr} Γ^r_{rr}
    # Γ^θ_{θθ} = 0, Γ^θ_{rr} = 0
    # = ∂_r(0.5 g^{θθ} g_θθ') - 0.5 g^{θθ} g_θθ' * 0.5 g^{rr} g_rr'
    Gth_rth_val = sp.Rational(1, 2) * ginv_thth * g_thth_p
    dGth_rth_dr = sp.diff(Gth_rth_val, r_sym)
    R_rthrth_upper = dGth_rth_dr - Gth_rth_val * Gr_rr_val
    R_rthrth = g_thth * R_rthrth_upper

    R_rthrth_simplified = sp.simplify(R_rthrth)
    print(f"R_rθrθ = {R_rthrth_simplified}")

    # Kretschmann scalar
    # K = 4 R_{trtr}² / (g_tt * g_rr)² + 4 R_{rθrθ}² / (g_rr * g_θθ)² + 4 R_{tθtθ}² / (g_tt * g_θθ)²
    # Wait, for diagonal metric:
    # K = 4 R_{trtr}² / (g_tt² * g_rr²) + 4 R_{rθrθ}² / (g_rr² * g_θθ²) + 4 R_{tθtθ}² / (g_tt² * g_θθ²)
    # But this isn't right for the sign conventions. Let me use:
    # K = 4 R_{trtr}² / (det_2d)² + 4 R_{rθrθ}² / (g_rr * g_θθ)² + 4 R_{tθtθ}² / (g_tt * g_θθ)²

    # Actually for a diagonal metric with signature (-,+,+,+) or (+,-,+,+):
    # K = 4 R_{trtr}² / (g_tt * g_rr)² + 4 R_{rθrθ}² / (g_rr * g_θθ)² + 4 R_{tθtθ}² / (g_tt * g_θθ)²
    # But inside the horizon, g_tt > 0 and g_rr < 0, so g_tt * g_rr < 0.
    # The Kretschmann should use the absolute values or the metric determinant.

    # For the standard formula:
    # K = 4 R_{trtr}² / (g_tt * g_rr)² + 4 R_{rθrθ}² / (g_rr * g_θθ)² + 4 R_{tθtθ}² / (g_tt * g_θθ)²

    # R_{tθtθ}: similar computation
    # R^θ_{tθt} = ∂_θ Γ^θ_{tt} - ∂_t Γ^θ_{tθ} + ... = 0 (spherical symmetry, static)
    # Actually: R_{tθtθ} = g_{tt} * Γ^θ_{tθ} * Γ^t_{tθ} ... no
    # For diagonal static spherically symmetric:
    # R_{tθtθ} = -0.5 * g_tt * g_θθ' / (g_tt * r) * ... this is getting complex

    # Let me use the known formula for Schwarzschild-like metrics:
    # For ds² = -f(r) dt² + h(r) dr² + r² dΩ²:
    # K = 4(f''/f - f'²/(2fh) + f'h'/(2fh))² / (4fh) + 4(f'/(2rh))² / f + 4(h-1)²/r⁴
    # But our metric has R(r) = r_h (constant), not R(r) = r.

    # For ds² = g_tt dt² + g_rr dr² + R² dΩ² with R = const:
    # The angular Riemann components simplify since R' = 0.

    # R_{rθrθ} = -0.5 * R'' + R' * Γ^r_{θθ} = 0 (since R = const, R' = R'' = 0)
    # Wait, R = g_θθ^{1/2} = r_h = const, so R' = 0, R'' = 0.
    # But g_θθ = R² = r_h², and g_θθ' = 2R*R' = 0.
    # So Γ^θ_{rθ} = 0.5 * g^{θθ} * g_θθ' = 0.
    # And R_{rθrθ} = 0.

    print(f"\n  g_θθ' = {sp.simplify(g_thth_p)}")
    print(f"  → g_θθ is CONSTANT (R = r_h = const)")
    print("  → R_{rθrθ} = 0 (angular Riemann vanishes)")

    # R_{tθtθ}:
    # R^θ_{tθt} = ∂_θ Γ^θ_{tt} - ∂_t Γ^θ_{tθ} + Γ^θ_{θα}Γ^α_{tt} - Γ^θ_{tα}Γ^α_{tθ}
    # Static, spherically symmetric: Γ^θ_{tt} = 0, Γ^θ_{tθ} = 0
    # So R_{tθtθ} = 0.

    print("  → R_{tθtθ} = 0 (by symmetry)")

    # So the ONLY nonzero Riemann component is R_{trtr}!
    # K = 4 * R_{trtr}² / (g_tt * g_rr)²

    Kretschmann = 4 * R_trtr**2 / (g_tt * g_rr)**2
    K_simplified = sp.simplify(Kretschmann)
    print(f"\n  Kretschmann = 4 * R_trtr² / (g_tt * g_rr)²")
    print(f"  = {K_simplified}")

    # Power law: substitute r → ε and expand
    K_series = sp.series(K_simplified, r_sym, 0, n=3)
    print(f"\n  Series expansion near r=0:")
    print(f"  K = {K_series}")

    # Leading power
    K_leading = sp.simplify(K_simplified * r_sym**6)  # multiply by r^6 to get leading coeff
    print(f"\n  K * r^6 = {sp.simplify(K_leading)}")
    print(f"  → K ~ C / r^6 where C = {sp.simplify(K_leading).subs(r_sym, 1)}")

    # Check: does K → 0 or ∞?
    K_limit = sp.limit(K_simplified, r_sym, 0)
    print(f"\n  lim(r→0) K = {K_limit}")

    # Also compute for general φ₀
    print("\n" + "=" * 60)
    print("GENERAL φ₀ ANALYSIS")
    print("=" * 60)

    phi_0 = sp.Symbol('phi_0', positive=True)
    A_gen = (r_h_sym / r_sym)**phi_0
    A2_gen = A_gen**2
    F_gen = -2 * M_sym / r_sym  # deep interior

    g_tt_gen = -A2_gen * F_gen  # = A² * 2M/r
    g_rr_gen = A2_gen / F_gen   # = -A² * r/(2M)
    g_thth_gen = A2_gen * r_sym**2  # = r_h^{2φ₀} * r^{2-2φ₀}

    print(f"\n  g_tt = {sp.simplify(g_tt_gen)}")
    print(f"  g_rr = {sp.simplify(g_rr_gen)}")
    print(f"  g_θθ = {sp.simplify(g_thth_gen)}")

    # g_θθ' for general φ₀
    g_thth_gen_p = sp.diff(g_thth_gen, r_sym)
    print(f"  g_θθ' = {sp.simplify(g_thth_gen_p)}")
    print(f"  → g_θθ' = 0 only if φ₀ = 1 (constant areal radius)")

    # For φ₀ ≠ 1, R_{rθrθ} ≠ 0, and we need the full Kretschmann.
    # For φ₀ = 1, R_{rθrθ} = 0 and K = 4*R_{trtr}²/(g_tt*g_rr)².

    # Compute R_{trtr} for general φ₀
    ginv_tt_gen = 1 / g_tt_gen
    ginv_rr_gen = 1 / g_rr_gen

    g_tt_gen_p = sp.diff(g_tt_gen, r_sym)
    g_rr_gen_p = sp.diff(g_rr_gen, r_sym)

    Gt_tr_gen = sp.Rational(1, 2) * ginv_tt_gen * g_tt_gen_p
    Gr_rr_gen = sp.Rational(1, 2) * ginv_rr_gen * g_rr_gen_p

    dGt_tr_gen = sp.diff(Gt_tr_gen, r_sym)
    R_trtr_gen = g_tt_gen * (dGt_tr_gen + Gt_tr_gen**2 - Gt_tr_gen * Gr_rr_gen)
    R_trtr_gen_simp = sp.simplify(R_trtr_gen)

    # R_{rθrθ} for general φ₀
    ginv_thth_gen = 1 / g_thth_gen
    Gth_rth_gen = sp.Rational(1, 2) * ginv_thth_gen * g_thth_gen_p
    dGth_rth_gen = sp.diff(Gth_rth_gen, r_sym)
    R_rthrth_gen = g_thth_gen * (dGth_rth_gen - Gth_rth_gen * Gr_rr_gen)
    R_rthrth_gen_simp = sp.simplify(R_rthrth_gen)

    # R_{tθtθ} for general φ₀
    # R_{tθtθ} = g_θθ * Γ^t_{tθ} * ... actually for diagonal metric:
    # R_{tθtθ} = -0.5 * ∂_θ g_tt + ... = 0 by spherical symmetry
    # Wait, no. R_{tθtθ} involves the θ-derivative, which is zero for
    # spherical symmetry. But there's also the angular part.
    # Actually: R_{tθtθ} = g_{tt} * R^t_{θtθ}
    # R^t_{θtθ} = ∂_t Γ^t_{θθ} - ∂_θ Γ^t_{tθ} + Γ^t_{tα}Γ^α_{θθ} - Γ^t_{θα}Γ^α_{tθ}
    # = 0 - 0 + Γ^t_{tt}Γ^t_{θθ} + Γ^t_{tr}Γ^r_{θθ} - Γ^t_{θθ}Γ^θ_{tθ} - Γ^t_{θr}Γ^r_{tθ}
    # Γ^t_{θθ} = 0, Γ^t_{θr} = 0, Γ^θ_{tθ} = 0
    # = Γ^t_{tr}Γ^r_{θθ}
    # Γ^r_{θθ} = -0.5 * g^{rr} * g_θθ'
    # So R_{tθtθ} = g_tt * Γ^t_{tr} * (-0.5 * g^{rr} * g_θθ')
    # = g_tt * 0.5 * g^{tt} * g_tt' * (-0.5 * g^{rr} * g_θθ')
    # = -0.25 * g_tt' * g_θθ' / (g_tt * g_rr) * g_tt
    # = -0.25 * g_tt' * g_θθ' / g_rr

    R_tthtth_gen = -sp.Rational(1, 4) * g_tt_gen_p * g_thth_gen_p / g_rr_gen
    R_tthtth_gen_simp = sp.simplify(R_tthtth_gen)

    # Full Kretschmann
    K_gen = (4 * R_trtr_gen**2 / (g_tt_gen * g_rr_gen)**2
             + 4 * R_rthrth_gen**2 / (g_rr_gen * g_thth_gen)**2
             + 4 * R_tthtth_gen**2 / (g_tt_gen * g_thth_gen)**2)

    K_gen_simp = sp.simplify(K_gen)

    # Power law: find leading behavior as r → 0
    # K * r^p should give a finite limit for the right p
    for p in [0, 2, 4, 6, 8, 10, 12]:
        K_scaled = sp.simplify(K_gen_simp * r_sym**p)
        K_lim = sp.limit(K_scaled, r_sym, 0)
        if K_lim != 0 and K_lim != sp.oo:
            print(f"\n  K * r^{p} → {K_lim} (finite!)")
            print(f"  → K ~ {K_lim} / r^{p}")
            if p > 0:
                print(f"  → K VANISHES as r→0 ✓")
            elif p < 0:
                print(f"  → K DIVERGES as r→0 ✗")
            else:
                print(f"  → K → const as r→0")
            break
        elif K_lim == 0 and p > 0:
            print(f"  K * r^{p} → 0 (K vanishes faster than r^{-p})")
        elif K_lim == sp.oo:
            print(f"  K * r^{p} → ∞ (K diverges faster than r^{-p})")

    # Compute for specific φ₀ values
    print("\n" + "=" * 60)
    print("KRETSCHEMANN FOR SPECIFIC φ₀ VALUES")
    print("=" * 60)

    for phi_0_val in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
        K_specific = K_gen_simp.subs(phi_0, phi_0_val)
        K_specific = sp.simplify(K_specific)

        # Find leading power
        for p in range(-10, 20):
            K_scaled = sp.simplify(K_specific * r_sym**p)
            K_lim = sp.limit(K_scaled, r_sym, 0)
            if K_lim != 0 and K_lim != sp.oo and K_lim != -sp.oo:
                areal = (r_h_sym / r_sym)**phi_0_val * r_sym
                areal_lim = sp.limit(areal, r_sym, 0)
                print(f"\n  φ₀ = {phi_0_val}:")
                print(f"    Areal radius → {areal_lim} {'(DIVERGES)' if areal_lim == sp.oo else '(finite)'}")
                print(f"    K ~ {K_lim} * r^{{-{p}}} = {K_lim} / r^{p}")
                if p > 0:
                    print(f"    → K DIVERGES ✗")
                elif p < 0:
                    print(f"    → K VANISHES ✓ (as r^{{{-p}}})")
                else:
                    print(f"    → K → const")
                break


if __name__ == '__main__':
    analytical_curvature_phi0_1()
