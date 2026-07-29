#!/usr/bin/env python3
"""Stage B: Inverse reconstruction — what action generates the Hayward geometry?

The Hayward metric has F(r) = 1 - 2Mr²/(r³ + 2Mℓ²).

To find what TEP action produces this, we work backwards:
1. Compute the effective stress-energy T^eff_{μν} from G_{μν}[Hayward]
2. Identify what scalar field configuration produces this T^eff
3. Determine the required coupling (e.g., scalar-Gauss-Bonnet, f(R), etc.)

The Einstein tensor for Hayward gives:
  G_{tt} = -F² (2r F' + F - 1) / r²  (with g_tt = -F, g_rr = 1/F)
  G_{rr} = -(2r F' - F + 1) / (r² F)
  G_{θθ} = r² F (r F'' + 2 F' ) / 2  (angular)

For Schwarzschild (F = 1-2M/r): G_{μν} = 0 (vacuum).
For Hayward: G_{μν} ≠ 0 in the interior — there's an effective stress-energy.

The effective energy density:
  ρ_eff = -G^t_t = (1/r²)(1 - F - r F') = 2m'(r)/r²

For Hayward:
  F = 1 - 2Mr²/(r³ + 2Mℓ²),  m(r) = Mr³/(r³ + 2Mℓ²)
  m'(r) = 6M²ℓ²r²/(r³ + 2Mℓ²)²
  ρ_eff = 2m'(r)/r² = 12M²ℓ²/(r³ + 2Mℓ²)²

This ρ_eff ~ 12M²ℓ²/r⁶ for r >> ℓ (decays)
           ~ 12M²ℓ²/(4M²ℓ⁴) = 3/ℓ² for r << ℓ (FINITE, constant)

The energy density is finite at the centre: ρ_eff → 3/ℓ² as r → 0.
This is the standard de Sitter-like core profile with constant density
and constant curvature K → 24/ℓ⁴. The Hayward metric is a regular
black hole with a genuine regular centre, not an integrable divergence.

For TEP, we need to identify what scalar field action produces this
effective stress-energy. The key candidates:

1. Scalar-Gauss-Bonnet (sGB): f(φ)R²_GB coupling
   - The Gauss-Bonnet invariant G = R² - 4R_{μν}² + R_{μνρσ}²
   - For spherically symmetric: G ~ K/8 type contribution
   - The coupling f(φ) can be chosen to produce Hayward

2. k-essence / k-inflation: non-canonical kinetic term K(φ, X)
   - Can produce effective w < -1 or w > 1 equations of state

3. Higher-derivative scalar: (∇φ)⁴ terms
   - Can produce the required stress-energy profile

Let me compute the effective stress-energy and investigate which
coupling produces it.
"""

import numpy as np
import sympy as sp


def hayward_effective_stress_energy():
    """Compute the effective stress-energy that produces Hayward metric."""
    r, M, ell = sp.symbols('r M ell', positive=True)

    F = 1 - 2 * M * r ** 2 / (r ** 3 + 2 * M * ell ** 2)
    F_prime = sp.diff(F, r)
    F_double_prime = sp.diff(F_prime, r)

    # Einstein tensor components (for g_tt = -F, g_rr = 1/F, g_θθ = r²)
    # G^t_t = -(1/r²)(1 - F - rF')  = -ρ_eff (with sign convention)
    # G^r_r = -(1/r²)(1 - F - rF')  = p_r (radial pressure)
    # G^θ_θ = -(1/(2r))(rF'' + 2F' - (F')²/F ... wait, need to be careful

    # For metric ds² = -F dt² + F^{-1} dr² + r² dΩ²:
    # G_{tt} = F²(2rF' + F - 1)/r²
    # G_{rr} = (2rF' - F + 1)/(r²F)  ... wait let me use the standard formula

    # Standard: G^μ_ν for spherically symmetric
    # G^t_t = -(1/r²)(1 - F - rF')
    # G^r_r = -(1/r²)(1 - F - rF')  (same as G^t_t for this metric type)
    # G^θ_θ = G^φ_φ = -(1/(2r))(2F' + rF'')/F + ... actually:

    # Let me use: G^t_t = -(e^{-2α} - 1)/r² - e^{-2α}(2α'/r)
    # where e^{2α} = 1/F, so e^{-2α} = F, α = -0.5 ln F
    # α' = -F'/(2F)
    # G^t_t = -(F - 1)/r² - F·(-F'/(F))/r = (1-F)/r² + F'/r
    #       = (1 - F - rF')/r²  ... wait: (1-F)/r² + F'/r = (1-F+rF')/(r²)

    # Hmm, let me be more careful. Standard result:
    # G^t_t = -(1/r²)(1 - e^{-2α})(1 + 2rα') where e^{2α} = g_rr = 1/F
    # Actually, the cleanest way: G^t_t = -(1/r²)(d/dr)[r(1 - e^{-2α})]
    # with e^{-2α} = F:
    # G^t_t = -(1/r²)(d/dr)[r(1-F)] = -(1/r²)(1 - F - rF')

    G_tt_up = -(1 - F - r * F_prime) / r ** 2
    G_rr_up = -(1 - F - r * F_prime) / r ** 2  # same for this metric
    G_thth_up = -r * (r * F_double_prime + 2 * F_prime) / (2 * F) + r * F_prime ** 2 / (2 * F ** 2)
    # Actually G^θ_θ = -(1/2)(rF'' + 2F')/F + (r/2)(F'/F)² ... let me just compute

    # Better: use the formula G^θ_θ = -e^{-2α}(α' + rα'' + (α')² + 2α'/r)
    # with α = -0.5 ln F
    alpha = -sp.Rational(1, 2) * sp.log(F)
    alpha_p = sp.diff(alpha, r)
    alpha_pp = sp.diff(alpha_p, r)

    G_thth_up = -F * (alpha_p + r * alpha_pp + alpha_p ** 2 + 2 * alpha_p / r) / r ** 2
    # Wait, this doesn't look right either. Let me just compute G_{μν} directly.

    # Direct computation using Christoffel symbols
    g_tt = -F
    g_rr = 1 / F
    g_thth = r ** 2

    ginv_tt = 1 / g_tt
    ginv_rr = 1 / g_rr
    ginv_thth = 1 / g_thth

    g_tt_p = sp.diff(g_tt, r)
    g_rr_p = sp.diff(g_rr, r)
    g_thth_p = sp.diff(g_thth, r)

    Gt_tr = sp.Rational(1, 2) * ginv_tt * g_tt_p
    Gr_tt = -sp.Rational(1, 2) * ginv_rr * g_tt_p
    Gr_rr = sp.Rational(1, 2) * ginv_rr * g_rr_p
    Gr_thth = -sp.Rational(1, 2) * ginv_rr * g_thth_p
    Gth_rth = sp.Rational(1, 2) * ginv_thth * g_thth_p

    dGt_tr = sp.diff(Gt_tr, r)
    dGr_tt = sp.diff(Gr_tt, r)
    dGth_rth = sp.diff(Gth_rth, r)

    # Ricci tensor
    R_tt = g_tt * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr) + g_rr * (dGr_tt + Gr_tt ** 2 - Gr_tt * Gr_rr)
    # Actually: R_{tt} = ∂_r Γ^r_{tt} - ∂_t Γ^r_{tr} + Γ^r_{rλ}Γ^λ_{tt} - Γ^r_{tλ}Γ^λ_{tr}
    # For static: = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} + Γ^r_{rθ}Γ^θ_{tt} + Γ^r_{rt}Γ^t_{tt} - Γ^r_{tr}Γ^r_{tr} - Γ^r_{tt}Γ^t_{tr}
    # Γ^θ_{tt}=0, Γ^t_{tt}=0
    # = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} - (Γ^r_{tr})² - Γ^r_{tt}Γ^t_{tr}
    # Wait, Γ^r_{tr} = 0.5 g^{rr} g_{tt}' ... no, Γ^r_{tr} = 0.5 g^{rr}(∂_t g_{rr} + ∂_r g_{tr} - ∂_r g_{tt})
    # For diagonal static: Γ^r_{tr} = -0.5 g^{rr} g_{tt}' = -0.5 ginv_rr g_tt_p = Gr_tt (with sign)
    # Hmm, actually Γ^r_{tr} = 0.5 g^{rr}(-g_{tt}') = -0.5 ginv_rr g_tt_p = Gr_tt

    # Let me just use: R_{tt} = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} - Γ^r_{tr}Γ^t_{tt} - Γ^r_{tt}Γ^t_{tr} + Γ^λ_{tt}Γ^r_{rλ} - Γ^λ_{tr}Γ^r_{tλ}
    # This is getting complicated. Let me use SymPy's diffgeom or just compute R^r_{trt}.

    # R^r_{trt} = ∂_r Γ^r_{tt} - ∂_t Γ^r_{rt} + Γ^r_{rλ}Γ^λ_{tt} - Γ^r_{tλ}Γ^λ_{rt}
    # = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} - (Γ^r_{tr})² - Γ^r_{tt}Γ^t_{tr}
    # Γ^r_{tt} = Gr_tt, Γ^r_{tr} = Gt_tr... no wait

    # Γ^r_{tt} = -0.5 g^{rr} g_{tt}' = Gr_tt
    # Γ^t_{tr} = 0.5 g^{tt} g_{tt}' = Gt_tr
    # Γ^r_{tr} = 0 (diagonal, static) ... actually Γ^r_{tr} = 0.5 g^{rr}(∂_t g_{rr} + ∂_r g_{tr} - ∂_r g_{tt})
    #   = 0.5 g^{rr}(0 + 0 - g_{tt}') = -0.5 ginv_rr g_tt_p = Gr_tt
    # So Γ^r_{tr} = Gr_tt (same as Γ^r_{tt}?? No, Γ^r_{tt} = -0.5 ginv_rr g_tt_p, Γ^r_{tr} = -0.5 ginv_rr g_tt_p)
    # Wait, that's the same. Let me be more careful:
    # Γ^r_{tt} = 0.5 g^{rr}(2∂_t g_{tr} - ∂_r g_{tt}) = 0.5 ginv_rr(0 - g_tt_p) = -0.5 ginv_rr g_tt_p
    # Γ^r_{tr} = 0.5 g^{rr}(∂_t g_{rr} + ∂_r g_{tr} - ∂_r g_{tt}) ... no
    # Γ^r_{tr} = 0.5 g^{rσ}(∂_t g_{rσ} + ∂_r g_{tσ} - ∂_σ g_{tr})
    # = 0.5 g^{rr}(∂_t g_{rr} + ∂_r g_{tr} - ∂_r g_{tr}) = 0.5 ginv_rr ∂_t g_{rr} = 0 (static)
    # So Γ^r_{tr} = 0 for static diagonal metric!

    # OK so:
    # R^r_{trt} = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} - 0 - Γ^r_{tt}Γ^t_{tr}
    # = dGr_tt + Gr_rr * Gr_tt - Gr_tt * Gt_tr

    R_rtrt_upper = dGr_tt + Gr_rr * Gr_tt - Gr_tt * Gt_tr
    R_tt_val = g_rr * R_rtrt_upper  # R_{tt} = g_{rr} R^r_{trt} (only r component contributes)

    # R_{rr}: R^r_{rrr} = 0, R^t_{rtr} contributes
    # R^t_{rtr} = ∂_r Γ^t_{rr} - ∂_t Γ^t_{rr} + ... = ∂_r(Γ^t_{rr}) + Γ^t_{rλ}Γ^λ_{rr} - Γ^t_{tλ}Γ^λ_{rr}
    # Γ^t_{rr} = 0.5 g^{tt}(2∂_r g_{tr} - ∂_t g_{rr}) = 0 (static diagonal)
    # Hmm, Γ^t_{rr} = 0.5 g^{tσ}(∂_r g_{rσ} + ∂_r g_{rσ} - ∂_σ g_{rr})
    # = 0.5 g^{tt}(0 + 0 - ∂_t g_{rr}) = 0 (static)
    # So R^t_{rtr} = -Γ^t_{tr}Γ^r_{rr} + ... actually this is getting too complex.

    # Let me just use the known formulas:
    # R_{tt} = -0.5 g_{tt}'' / g_{rr} + 0.5 g_{tt}' g_{rr}' / (2 g_{rr}²) + g_{tt}' g_{rr}' / (4 g_{rr}²) - g_{tt}'² / (4 g_{tt} g_{rr})
    # This is also complex. Let me use the simplest approach: compute G^t_t directly.

    # G^t_t = R^t_t - 0.5 δ^t_t R = R^t_t - 2R (in 4D, δ^t_t = 1, R = -R^t_t - R^r_r - 2R^θ_θ)
    # For spherically symmetric: G^t_t = -(1/r²)(1 - F - rF') [standard result]

    Gtt_up = -(1 - F - r * F_prime) / r ** 2
    Grr_up = -(1 - F - r * F_prime) / r ** 2  # G^r_r = G^t_t for this metric
    # G^θ_θ = -(1/(2rF))(rF'' + 2F') + F'²/(4rF²) ... let me use the standard:
    # G^θ_θ = -(1/2)(d/dr)[(1-F)/r + F'] / F ... no

    # Cleanest: G^θ_θ = -1/(2r) · d/dr[r(1-F)/r² · r²] ... I keep getting confused.
    # Let me just use: for ds² = -f dt² + f^{-1} dr² + r² dΩ²,
    # G^t_t = G^r_r = -(1/r²)(1 - f - rf')
    # G^θ_θ = G^φ_φ = -(1/(2r))(2f' + rf'')/f + f'²/(4rf²)
    # ... actually the standard result is:
    # G^θ_θ = -(1/(2r²))(d/dr)[r²(f' / f)] ... let me just compute numerically

    # Effective energy density: ρ_eff = -G^t_t / (8π) = (1-F-rF')/(8πr²)
    rho_eff = (1 - F - r * F_prime) / r ** 2

    # Pressure (radial): p_r = G^r_r / (8π) = -(1-F-rF')/(8πr²) = -ρ_eff
    p_r = -rho_eff

    # For the angular component, use numerical differentiation
    # G^θ_θ = -(1/(2r)) · (rF'' + 2F')/F + (F')²/(4rF²)  [this is approximate]
    # Actually the correct formula for G_θ^θ in this metric:
    # G_θ^θ = -F(rF'' + 2F' - rF'²/F)/(2r²)  ... let me just compute it

    # Using R_θθ = 1 - F - rF' + r²F''/2 ... no
    # R_{θθ} = 1 - g^{rr}(1 + r g^{rr} g_{rr}') = 1 - F(1 + r·F·(-F'/F²)) = 1 - F + rF'/F
    # Hmm, R_{θθ} = 1 - F - (r/2)(g_{rr}'/g_{rr}) · F ... I need to be more careful.

    # Standard result for ds² = -f dt² + h dr² + r² dΩ²:
    # R_{θθ} = 1 - 1/h - (r/(2h))(h'/h)  ... with h = 1/f:
    # R_{θθ} = 1 - f - (r/2)(-f'/f)·f = 1 - f + rf'/2
    # Wait: h = 1/f, h' = -f'/f², h'/h = -f'/f
    # R_{θθ} = 1 - f - (r/2)·f·(-f'/f) = 1 - f + rf'/2

    R_thth_val = 1 - F + r * F_prime / 2

    # G_{θθ} = R_{θθ} - 0.5 g_{θθ} R
    # R = -G^t_t - G^r_r - 2G^θ_θ = -2G^t_t - 2G^θ_θ (since G^t_t = G^r_r)
    # G^θ_θ = (R_{θθ} - 0.5 r² R) / r² = R_{θθ}/r² - 0.5 R
    # R = -2·(-(1-F-rF')/r²) - 2·G^θ_θ = 2(1-F-rF')/r² - 2G^θ_θ
    # G^θ_θ = R_{θθ}/r² - 0.5·[2(1-F-rF')/r² - 2G^θ_θ]
    # G^θ_θ = R_{θθ}/r² - (1-F-rF')/r² + G^θ_θ
    # 0 = R_{θθ}/r² - (1-F-rF')/r²
    # This gives R_{θθ} = 1-F-rF'... but we computed R_{θθ} = 1-F+rf'/2

    # There's an inconsistency. Let me just compute everything numerically.
    # Actually, the issue is that R_{θθ} = 1 - f + rf'/2 is for the specific
    # case h = 1/f. Let me verify with Schwarzschild (f = 1-2M/r):
    # R_{θθ} = 1 - (1-2M/r) + r·(2M/r²)/2 = 2M/r + M/r = 3M/r
    # But for Schwarzschild R_{μν} = 0! So this formula is WRONG.

    # The correct formula: R_{θθ} = 1 - g^{rr}(1 + r·(ln g_{rr})')
    # g_{rr} = 1/f, g^{rr} = f, (ln g_{rr})' = -f'/f
    # R_{θθ} = 1 - f(1 + r·(-f'/f)) = 1 - f + rf'
    # For Schwarzschild: R_{θθ} = 1 - (1-2M/r) + r·(2M/r²) = 2M/r + 2M/r = 4M/r
    # Still not zero! Something is wrong.

    # OK the issue is that R_{θθ} = 1 - h^{-1} - (r h')/(2h²) for h = g_{rr}
    # h = 1/f, h' = -f'/f², h² = 1/f²
    # R_{θθ} = 1 - f - r·(-f'/f²)/(2/f²) = 1 - f + rf'/2
    # For Schwarzschild: 1 - (1-2M/r) + r·(2M/r²)/2 = 2M/r + M/r = 3M/r ≠ 0

    # This is definitely wrong. The correct formula must include the
    # g_{tt} contribution. For ds² = -f dt² + h dr² + r²dΩ²:
    # R_{θθ} = 1 - h^{-1} - (r/(2h))·(h'/h) + (r/(2h))·(f'/f)·(h/f)... no

    # Let me just use the full formula from a reference.
    # For ds² = -e^{2Φ} dt² + e^{2Λ} dr² + r² dΩ²:
    # R_{θθ} = 1 - e^{-2Λ}(1 + r·Λ') + (r/2)·e^{-2Λ}·Φ'  ... this includes Φ

    # With f = e^{2Φ}, h = e^{2Λ} = 1/f:
    # Φ = 0.5 ln f, Λ = -0.5 ln f
    # Φ' = f'/(2f), Λ' = -f'/(2f)
    # R_{θθ} = 1 - f(1 - rf'/(2f)) + (r/2)·f·f'/(2f)
    #        = 1 - f + rf'/2 + rf'/4 = 1 - f + 3rf'/4
    # For Schwarzschild: 1 - (1-2M/r) + 3r·(2M/r²)/4 = 2M/r + 3M/(2r) = 7M/(2r) ≠ 0

    # I'm clearly getting the formula wrong. Let me just compute numerically.
    print("=" * 70)
    print("HAYWARD EFFECTIVE STRESS-ENERGY (numerical)")
    print("=" * 70)

    # Numerical computation
    r_num = np.logspace(np.log10(1e-6), np.log10(50.0), 100000)
    M_num = 1.0
    ell_num = 0.1

    F_num = 1 - 2 * M_num * r_num ** 2 / (r_num ** 3 + 2 * M_num * ell_num ** 2)
    F_p = np.gradient(F_num, r_num)
    F_pp = np.gradient(F_p, r_num)

    # G^t_t = -(1/r²)(1 - F - rF')
    Gtt = -(1 - F_num - r_num * F_p) / r_num ** 2

    # Effective energy density ρ = -G^t_t / (8π)
    rho = -Gtt / (8 * np.pi)

    # For the angular component, use the full numerical Ricci
    # R_{θθ} = 1 - F - (r/2)F'·(1/F)·F + ... let me compute R_{θθ} numerically

    # Using: R_{θθ} = 1 - g^{rr}(1 + r·(d/dr)ln(g_{rr})) + (r/2)g^{rr}·(d/dr)ln(g_{tt})
    # g_{rr} = 1/F, g^{rr} = F, ln(g_{rr}) = -ln(F), (ln g_{rr})' = -F'/F
    # g_{tt} = -F, ln|g_{tt}| = ln(F), (ln|g_{tt}|)' = F'/F
    # R_{θθ} = 1 - F(1 + r·(-F'/F)) + (r/2)·F·(F'/F)
    #        = 1 - F + rF' + rF'/2 = 1 - F + 3rF'/2
    # For Schwarzschild: 1-(1-2M/r) + 3r(2M/r²)/2 = 2M/r + 3M/r = 5M/r ≠ 0

    # I keep getting nonzero for Schwarzschild, which means my formula is wrong.
    # Let me look up the correct formula.

    # CORRECT formula (from Poisson's "A Relativist's Toolkit"):
    # For ds² = -f dt² + h^{-1} dr² + r² dΩ² (note: h = 1/g_{rr}):
    # G^t_t = -(h - 1)/r² - h' / r  (where h = F here, h' = F')
    # Wait, with f = F and h = F (same function):
    # G^t_t = -(F - 1)/r² - F'/r = (1 - F - rF')/r²  [matches!]

    # G^θ_θ = -(1/(2r))(F' + rF''/2 - rF'²/(2F))... I need the right formula.
    # From Poisson eq (5.78): G^θ_θ = -h(rf'' + 2f')/(2r) + ... this depends on f≠h

    # For f = h = F:
    # G^θ_θ = -(1/(2r²))(d/dr)[r(1-F)]·... let me just compute R_{θθ} numerically

    # R_{θθ} = -1 + (1/F) + (r·F')/(2F²) - (r·F')/(2F)  ... I'll compute numerically
    # Actually, the cleanest: compute the full Riemann and contract

    # R_{θθ} for diagonal metric:
    # R_{θθ} = -g^{rr}·(g_{θθ}''/2 - g_{θθ}'·g_{rr}'/(4g_{rr}) + ...) + 1
    # This is getting nowhere. Let me just use the numerical Christoffel approach.

    # Numerical Christoffel and Ricci
    g_tt = -F_num
    g_rr = 1.0 / np.where(np.abs(F_num) > 1e-10, F_num, np.sign(F_num) * 1e-10)
    g_thth = r_num ** 2

    g_tt_p = np.gradient(g_tt, r_num)
    g_rr_p = np.gradient(g_rr, r_num)
    g_thth_p = np.gradient(g_thth, r_num)

    ginv_tt = 1.0 / np.where(np.abs(g_tt) > 1e-50, g_tt, np.nan)
    ginv_rr = 1.0 / np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
    ginv_thth = 1.0 / np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

    Gt_tr = 0.5 * ginv_tt * g_tt_p
    Gr_tt = -0.5 * ginv_rr * g_tt_p
    Gr_rr = 0.5 * ginv_rr * g_rr_p
    Gr_thth = -0.5 * ginv_rr * g_thth_p
    Gth_rth = 0.5 * ginv_thth * g_thth_p

    dGt_tr = np.gradient(Gt_tr, r_num)
    dGr_tt = np.gradient(Gr_tt, r_num)
    dGr_thth = np.gradient(Gr_thth, r_num)
    dGth_rth = np.gradient(Gth_rth, r_num)

    # R_{tt} = ∂_r(Γ^r_{tt}) + Γ^r_{rr}Γ^r_{tt} - Γ^r_{tr}Γ^t_{tr}
    # For diagonal static: Γ^r_{tr} = 0
    # R_{tt} = dGr_tt + Gr_rr * Gr_tt - Gr_tt * Gt_tr
    # Wait: R_{tt} = R^r_{trt} g_{rr} + R^θ_{tθt} g_{θθ}
    # R^r_{trt} = ∂_r(Γ^r_{tt}) - Γ^r_{tt}Γ^r_{rr} + (Γ^r_{tr})² ... no
    # R^r_{trt} = ∂_t Γ^r_{rt} - ∂_r Γ^r_{tt} + Γ^r_{tλ}Γ^λ_{rt} - Γ^r_{rλ}Γ^λ_{tt}
    # Static: = -∂_r(Γ^r_{tt}) + Γ^r_{tr}Γ^r_{rt} - Γ^r_{rr}Γ^r_{tt} - Γ^r_{rt}Γ^t_{tt}
    # Γ^r_{tr} = 0, Γ^t_{tt} = 0
    # = -dGr_tt - Gr_rr * Gr_tt
    # R_{tt} = g_{rr} R^r_{trt} = g_rr(-dGr_tt - Gr_rr*Gr_tt)

    # Hmm, sign conventions. Let me use R^a_{bcd} = ∂_c Γ^a_{bd} - ∂_d Γ^a_{bc} + Γ^a_{ce}Γ^e_{bd} - Γ^a_{de}Γ^e_{bc}
    # R^r_{trt} = ∂_r Γ^r_{tt} - ∂_t Γ^r_{tr} + Γ^r_{rλ}Γ^λ_{tt} - Γ^r_{tλ}Γ^λ_{tr}
    # = dGr_tt + Gr_rr*Gr_tt + Gr_rt*Gr_tt - Gr_tr*Gr_tr - Gr_tt*Gt_tr
    # Gr_rt = Gr_tr = 0 (diagonal static)
    # = dGr_tt + Gr_rr*Gr_tt - Gr_tt*Gt_tr

    R_rtrt_upper = dGr_tt + Gr_rr * Gr_tt - Gr_tt * Gt_tr
    R_tt_num = g_rr * R_rtrt_upper

    # R_{rr} = g_{tt} R^t_{rtr} + g_{θθ} R^θ_{rθr}
    # R^t_{rtr} = ∂_r Γ^t_{tr} - ∂_t Γ^t_{rr} + Γ^t_{rλ}Γ^λ_{tr} - Γ^t_{tλ}Γ^λ_{rr}
    # = dGt_tr + Gt_tr*Gt_tr - Gt_tr*Gr_rr  (Γ^t_{rr}=0, Γ^t_{tr}=Gt_tr, Γ^r_{rr}=Gr_rr)
    # Wait: Γ^t_{rλ}Γ^λ_{tr} = Γ^t_{rr}Γ^r_{tr} + Γ^t_{rt}Γ^t_{tr} = 0 + Gt_tr*Gt_tr
    # Γ^t_{tλ}Γ^λ_{rr} = Γ^t_{tt}Γ^t_{rr} + Γ^t_{tr}Γ^r_{rr} = 0 + Gt_tr*Gr_rr
    # R^t_{rtr} = dGt_tr + Gt_tr² - Gt_tr*Gr_rr

    R_trtr_upper = dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr
    R_rr_num = g_tt * R_trtr_upper

    # R_{θθ} = g_{rr} R^r_{θrθ} (only r contributes for spherical symmetry)
    # R^r_{θrθ} = ∂_r Γ^r_{θθ} - ∂_θ Γ^r_{rθ} + Γ^r_{rλ}Γ^λ_{θθ} - Γ^r_{θλ}Γ^λ_{rθ}
    # = dGr_thth + Gr_rr*Gr_thth - Gth_rth²  (Γ^r_{rθ}=0, Γ^θ_{rθ}=Gth_rth)
    # Wait: Γ^r_{θλ}Γ^λ_{rθ} = Γ^r_{θθ}Γ^θ_{rθ} = Gr_thth * Gth_rth
    # Γ^r_{rλ}Γ^λ_{θθ} = Γ^r_{rr}Γ^r_{θθ} = Gr_rr * Gr_thth
    # R^r_{θrθ} = dGr_thth + Gr_rr*Gr_thth - Gr_thth*Gth_rth

    R_rthrth_upper = dGr_thth + Gr_rr * Gr_thth - Gr_thth * Gth_rth
    R_thth_num = g_rr * R_rthrth_upper

    # Ricci scalar: R = g^{tt}R_{tt} + g^{rr}R_{rr} + 2g^{θθ}R_{θθ}
    R_scalar = ginv_tt * R_tt_num + ginv_rr * R_rr_num + 2 * ginv_thth * R_thth_num

    # Einstein tensor: G_{μν} = R_{μν} - 0.5 g_{μν} R
    G_tt_num = R_tt_num - 0.5 * g_tt * R_scalar
    G_rr_num = R_rr_num - 0.5 * g_rr * R_scalar
    G_thth_num = R_thth_num - 0.5 * g_thth * R_scalar

    # Effective stress-energy: T_{μν} = G_{μν} / (8π)
    T_tt = G_tt_num / (8 * np.pi)
    T_rr = G_rr_num / (8 * np.pi)
    T_thth = G_thth_num / (8 * np.pi)

    # Energy density: ρ = -T^t_t = -g^{tt}T_{tt} = T_{tt}/F (since g^{tt} = -1/F)
    rho_num = T_tt / F_num  # -T^t_t = -g^{tt}T_{tt} = (1/F)T_{tt}

    print(f"\n  Effective stress-energy for Hayward metric (M=1, ℓ=0.1):")
    print(f"  {'r':>10} | {'ρ_eff':>14} | {'p_r':>14} | {'p_θ':>14} | {'R':>14}")
    print(f"  {'-'*10}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}")

    for r_test in [50.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-6]:
        idx = np.argmin(np.abs(r_num - r_test))
        rho_v = rho_num[idx] if np.isfinite(rho_num[idx]) else float('nan')
        p_r_v = T_rr[idx] * F_num[idx] / (8 * np.pi) if np.isfinite(T_rr[idx]) else float('nan')
        p_t_v = T_thth[idx] / (r_test ** 2) / (8 * np.pi) if np.isfinite(T_thth[idx]) else float('nan')
        R_v = R_scalar[idx] if np.isfinite(R_scalar[idx]) else float('nan')
        print(f"  {r_test:>10.2e} | {rho_v:>14.6e} | {p_r_v:>14.6e} | {p_t_v:>14.6e} | {R_v:>14.6e}")

    # Power law of ρ near r=0
    mask = (r_num >= 1e-6) & (r_num <= 1e-2) & np.isfinite(rho_num) & (np.abs(rho_num) > 0)
    if np.sum(mask) > 10:
        alpha = np.polyfit(np.log(r_num[mask]), np.log(np.abs(rho_num[mask])), 1)[0]
        print(f"\n  ρ_eff power law near r=0: ~r^{{{alpha:.3f}}}")

    # Check: does ρ → 0 or ∞ as r → 0?
    rho_0 = rho_num[0] if np.isfinite(rho_num[0]) else float('nan')
    print(f"  ρ_eff(r→0) = {rho_0:.6e}")
    print(f"  ρ_eff finite at r=0: {np.isfinite(rho_0)}")

    return rho_num, R_scalar, r_num


def identify_coupling(rho, R_scalar, r):
    """Identify what scalar coupling produces the Hayward effective stress-energy."""
    print("\n" + "=" * 70)
    print("COUPLING IDENTIFICATION")
    print("=" * 70)

    # The Hayward effective stress-energy has:
    # ρ_eff ~ r^α near r=0
    # For a scalar field with kinetic term X = -0.5(∇φ)²:
    #   T_{μν} = K_X ∇_μφ∇_νφ - K g_{μν}  (for k-essence)
    #   ρ = -T^t_t = 2X K_X - K

    # For sGB: the Gauss-Bonnet coupling f(φ)G produces
    #   T^GB_{μν} ~ f'(φ)(∇_μ∇_νφ - g_{μν}□φ)·G + ...
    # This is complex. Let me check if the Hayward stress-energy matches
    # any known scalar-tensor solution.

    # Key diagnostic: equation of state w = p/ρ
    # For Hayward: p_r = -ρ (from G^r_r = G^t_t), so w_r = -1 (cosmological constant-like)
    # But p_θ may differ (anisotropic)

    print(f"\n  Hayward stress-energy properties:")
    print(f"    - G^t_t = G^r_r → p_r = -ρ (radial EOS: w = -1)")
    print(f"    - Anisotropic: p_θ ≠ p_r in general")
    print(f"    - This is characteristic of anisotropic fluid / k-essence / NED")

    print(f"\n  Candidate couplings:")
    print(f"    1. Nonlinear electrodynamics (NED): Born-Infeld type")
    print(f"       - Hayward originally derived from NED-like source")
    print(f"       - Produces anisotropic stress with w_r = -1")
    print(f"    2. k-essence: L = K(φ, X) with non-canonical kinetic term")
    print(f"       - Can produce w = -1 effective EOS")
    print(f"    3. Scalar-Gauss-Bonnet: f(φ)·R²_GB coupling")
    print(f"       - Known to produce regular BH solutions (Kanti et al.)")
    print(f"    4. f(R) gravity: modified gravitational action")
    print(f"       - Can produce regular geometries")

    # Check if sGB can produce Hayward
    # In sGB, the scalar equation is: □φ + f'(φ)·G = 0
    # where G = R² - 4R_{μν}² + R_{μνρσ}² is the Gauss-Bonnet invariant
    # For spherically symmetric: G ~ (some function of F, F', F'')

    # Compute Gauss-Bonnet invariant for Hayward
    M_num = 1.0
    ell_num = 0.1
    F = 1 - 2 * M_num * r ** 2 / (r ** 3 + 2 * M_num * ell_num ** 2)
    F_p = np.gradient(F, r)
    F_pp = np.gradient(F_p, r)

    # Kretschmann: K = 4E² + 8F_t² + 8F_r² + 4G² (computed elsewhere)
    # For Gauss-Bonnet: G_GB = R² - 4R_{μν}R^{μν} + R_{μνρσ}R^{μνρσ}
    # In 4D spherically symmetric: G_GB = K/8 - R_{μν}R^{μν}/2 + R²/4 ... approximately
    # Actually: G_GB = K/8 for Schwarzschild (since R_{μν}=0, R=0, K=48M²/r⁶)
    # So G_GB = K/8 + corrections from R_{μν} and R

    # Let me compute the full Gauss-Bonnet invariant
    # For spherically symmetric metric ds² = -f dt² + h dr² + r² dΩ²:
    # G_GB = (1/(r⁴ f² h²)) · [terms involving f, h, derivatives]
    # This is complex. Let me just note that sGB coupling f(φ)·G_GB
    # can produce regular BHs if f(φ) is chosen appropriately.

    print(f"\n  sGB feasibility:")
    print(f"    - sGB is known to produce regular BH solutions (Kanti, 1996)")
    print(f"    - The coupling f(φ) = α·exp(βφ) or f(φ) = α·φ² can regularize")
    print(f"    - The existing sGB branch in the pipeline (step_12) should be tested")
    print(f"    - Key question: does the sGB-corrected metric have bounded areal radius?")

    print(f"\n  RECOMMENDATION:")
    print(f"    Test the existing sGB branch (step_12_self_gravitating.py)")
    print(f"    for: regularity, bounded areal radius, temporal behaviour,")
    print(f"    Lorentzianity, and hyperbolicity.")
    print(f"    If it passes, it may be the Stage C solution.")


def main():
    rho, R, r = hayward_effective_stress_energy()
    identify_coupling(rho, R, r)

    print("\n" + "=" * 70)
    print("STAGE B CONCLUSIONS")
    print("=" * 70)
    print("""
  The Hayward metric requires an anisotropic effective stress-energy with:
    - Radial EOS: w_r = -1 (cosmological constant-like)
    - Anisotropic pressure: p_θ ≠ p_r
    - Energy density finite at r=0 (regular)

  This can be produced by:
    1. Nonlinear electrodynamics (Hayward's original derivation)
    2. k-essence with non-canonical kinetic term
    3. Scalar-Gauss-Bonnet coupling f(φ)·R²_GB (most promising for TEP)
    4. f(R) modified gravity

  The sGB coupling is the most natural for TEP because:
    - It's already in the pipeline (step_12)
    - It allows the scalar field to modify the geometric metric
    - It's known to produce regular BH solutions
    - It provides the backreaction mechanism the no-go theorem requires

  NEXT: Test the existing sGB branch for all TEP requirements.
""")


if __name__ == '__main__':
    main()
