#!/usr/bin/env python3
"""Stage A: TEP matter metric on a regular black-hole benchmark geometry.

The fixed-Schwarzschild no-go theorem (docs/fixed_schwarzschild_no_go_theorem.md)
proves that with g_{μν} fixed to Schwarzschild, curvature regularity (φ₀ > 3/2)
and bounded areal radius (φ₀ ≤ 1) are mutually exclusive.

The resolution: the temporal field must backreact on the geometric metric.
Before attempting the full coupled derivation (Stage C), we test whether
the intended temporal picture is mathematically possible once the geometric
singularity is removed.

STAGE A: Place the TEP matter metric on a HAYWARD regular black-hole background.

Hayward metric (2006):
  ds² = -(1 - 2Mr²/(r³ + 2Mℓ²)) dt² + (1 - 2Mr²/(r³ + 2Mℓ²))^{-1} dr²
        + r² dΩ²

  F_H(r) = 1 - 2Mr²/(r³ + 2Mℓ²)
         = (r³ + 2Mℓ² - 2Mr²) / (r³ + 2Mℓ²)
         = (r³ - 2Mr² + 2Mℓ²) / (r³ + 2Mℓ²)

Properties:
  - F_H(0) = 1 (Minkowski at center, no singularity)
  - F_H(r_h) = 0 (horizon, same as Schwarzschild for large r)
  - K_Hayward(0) = finite (curvature regular)
  - K_Hayward ~ 48M²/r⁶ for r >> ℓ (asymptotically Schwarzschild)
  - ℓ is the regularization length scale (typically ℓ ~ Planck length)

KEY ADVANTAGE: With a regular geometric metric, we can use BOUNDED A
(φ₀ ≤ 1, even A → const) and still have K → 0, because the geometric
metric itself is regular. The conformal factor no longer needs to do
all the work.

The temporal distortion can then come from the DISFORMAL term B∇φ∇φ
with a time-component scalar φ = qv + ψ(r), producing extreme temporal
shear without spatial enlargement.

This script:
1. Implements the Hayward metric
2. Verifies its curvature regularity
3. Places the TEP matter metric on it with BOUNDED A
4. Checks all gates: bounded areal radius, vanishing K, Lorentzianity
"""

import numpy as np
import sympy as sp


# =============================================================================
# HAYWARD METRIC
# =============================================================================

def hayward_F(r, M=1.0, ell=0.1):
    """Hayward metric function F(r) = 1 - 2Mr²/(r³ + 2Mℓ²).

    Properties:
      F(0) = 1 (regular center)
      F(r_h) = 0 (horizon)
      F → 1 - 2M/r for r >> ℓ (Schwarzschild)
    """
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-30)
    numerator = 2.0 * M * r_safe ** 2
    denominator = r_safe ** 3 + 2.0 * M * ell ** 2
    return 1.0 - numerator / denominator


def hayward_F_prime(r, M=1.0, ell=0.1):
    """Derivative dF/dr for Hayward metric."""
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-30)
    # F = 1 - 2Mr²/(r³ + 2Mℓ²)
    # dF/dr = -2M * [2r(r³+2Mℓ²) - r²·3r²] / (r³+2Mℓ²)²
    #       = -2M * [2r⁴ + 4Mℓ²r - 3r⁴] / (r³+2Mℓ²)²
    #       = -2M * [4Mℓ²r - r⁴] / (r³+2Mℓ²)²
    #       = 2M * r(r³ - 4Mℓ²) / (r³+2Mℓ²)²
    num = 2.0 * M * r_safe * (r_safe ** 3 - 4.0 * M * ell ** 2)
    den = (r_safe ** 3 + 2.0 * M * ell ** 2) ** 2
    return num / den


def find_hayward_horizon(M=1.0, ell=0.1):
    """Find the horizon radius r_h where F(r_h) = 0."""
    r = np.logspace(np.log10(ell * 0.01), np.log10(100 * M), 100000)
    F = hayward_F(r, M, ell)
    # Find sign changes
    sign_changes = np.where(np.diff(np.sign(F)))[0]
    if len(sign_changes) > 0:
        # Largest root (outer horizon)
        idx = sign_changes[-1]
        # Linear interpolation
        r1, r2 = r[idx], r[idx + 1]
        F1, F2 = F[idx], F[idx + 1]
        r_h = r1 - F1 * (r2 - r1) / (F2 - F1)
        return r_h
    return 2.0 * M  # fallback


# =============================================================================
# HAYWARD CURVATURE (analytical verification)
# =============================================================================

def hayward_kretschmann_analytical():
    """Verify Hayward metric is curvature-regular at r=0 using SymPy."""
    r, M, ell = sp.symbols('r M ell', positive=True)

    F = 1 - 2 * M * r ** 2 / (r ** 3 + 2 * M * ell ** 2)

    # Standard coordinates: g_tt = -F, g_rr = 1/F, g_θθ = r²
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
    Gr_rr = sp.Rational(1, 2) * ginv_rr * g_rr_p
    Gth_rth = sp.Rational(1, 2) * ginv_thth * g_thth_p

    dGt_tr = sp.diff(Gt_tr, r)
    dGth_rth = sp.diff(Gth_rth, r)

    R_trtr = g_tt * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr)
    R_rthrth = g_thth * (dGth_rth + Gth_rth ** 2 - Gth_rth * Gr_rr)
    R_tthtth = -sp.Rational(1, 4) * g_tt_p * g_thth_p / g_rr
    R_thphthph = g_thth * (1 - ginv_rr * g_thth_p ** 2 / (4 * g_thth))

    E = sp.simplify(R_trtr / (g_tt * g_rr))
    Ft = sp.simplify(R_tthtth / (g_tt * g_thth))
    Fr = sp.simplify(R_rthrth / (g_rr * g_thth))
    G_comp = sp.simplify(R_thphthph / g_thth ** 2)

    K = 4 * E ** 2 + 8 * Ft ** 2 + 8 * Fr ** 2 + 4 * G_comp ** 2
    K_simplified = sp.simplify(K)

    # Limit at r=0
    K_limit = sp.limit(K_simplified, r, 0)

    print("=" * 70)
    print("HAYWARD METRIC: ANALYTICAL KRETSCHMANN")
    print("=" * 70)
    print(f"  F = {sp.simplify(F)}")
    print(f"  K = {K_simplified}")
    print(f"  K(r→0) = {K_limit}")
    print(f"  K finite at r=0: {K_limit != sp.oo and K_limit != -sp.oo}")

    # Schwarzschild limit (ℓ → 0)
    K_schw_limit = sp.limit(K_simplified, ell, 0)
    print(f"  K(ℓ→0) = {sp.simplify(K_schw_limit)} (should be 48M²/r⁶)")

    return K_limit


# =============================================================================
# TEP MATTER METRIC ON HAYWARD BACKGROUND
# =============================================================================

def scalar_profile_bounded(r, M, r_h, phi_0=0.5, delta=0.05):
    """Bounded scalar field profile for use with regular background.

    With a regular geometric metric, we do NOT need A → ∞.
    We use a MODERATE φ₀ (e.g., 0.5) so that:
      - A = (r_h/r)^{φ₀} grows but stays bounded enough
      - The disformal term provides temporal distortion
      - The areal radius ρ = A·r stays finite

    For φ₀ = 0.5: ρ = A·r = (r_h/r)^{0.5} · r = (r_h·r)^{0.5} → 0 as r → 0
    For φ₀ = 0:   ρ = r → 0 (no conformal stretching at all)
    For φ₀ = 1:   ρ = r_h (constant, perfectly bounded)
    """
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-30)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    phi = phi_0 * np.log(r_safe / r_h) * S
    return phi


def scalar_gradient_bounded(r, M, r_h, phi_0=0.5, delta=0.05):
    """Gradient of bounded scalar profile."""
    r = np.asarray(r, dtype=float)
    r_safe = np.maximum(r, 1e-30)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)
    dphi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)
    return dphi


def compute_tep_hayward_metric(r, M=1.0, ell=0.1, phi_0=1.0, delta=0.05,
                                beta_A=-1.0, B0=1.0, sigma_B=1.5, q=0.1):
    """Compute TEP matter metric on Hayward background.

    g̃ = A² g_Hayward + B ∇φ ∇φ

    With φ = qv + ψ(r) (time-component scalar for temporal distortion):
      g̃_vv = -A²F + Bq²
      g̃_vr = A² + Bqψ'   (EF coordinates)
      g̃_rr = Bψ'²
      g̃_θθ = A²r²

    The KEY difference from the Schwarzschild case:
      - g_Hayward is regular at r=0 (F(0)=1, K_Hayward(0)=finite)
      - A can be BOUNDED (φ₀ ≤ 1) since the geometric metric is already regular
      - The disformal term with q≠0 provides temporal distortion
      - Areal radius ρ = A·r stays finite for φ₀ ≤ 1
    """
    r = np.asarray(r, dtype=float)
    r_h = find_hayward_horizon(M, ell)

    F = hayward_F(r, M, ell)
    phi = scalar_profile_bounded(r, M, r_h, phi_0, delta)
    dphi = scalar_gradient_bounded(r, M, r_h, phi_0, delta)

    A = np.exp(beta_A * phi)
    A2 = A ** 2

    # Disformal function (same quartic Gaussian bump)
    abs_phi = np.abs(phi)
    B = B0 * abs_phi ** 2 / (1.0 + abs_phi ** 2) * np.exp(-(phi ** 4) / (2.0 * sigma_B ** 4))

    # EF coordinates matter metric (with q≠0 for temporal distortion)
    gtilde_vv = -A2 * F + B * q ** 2
    gtilde_vr = A2 + B * q * dphi  # Note: +Bqψ' (not just A²)
    gtilde_rr = B * dphi ** 2
    gtilde_thth = A2 * r ** 2

    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2
    lorentzian = det_2d < 0

    # Areal radius
    areal_radius = A * r

    # Disformal ratio
    disformal_ratio = np.abs(F) * B * dphi ** 2 / np.where(A2 > 1e-30, A2, 1e-30)

    return {
        'r': r, 'F': F, 'phi': phi, 'dphi': dphi,
        'A': A, 'B': B, 'A2': A2,
        'gtilde_vv': gtilde_vv, 'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr, 'gtilde_thth': gtilde_thth,
        'det_2d': det_2d, 'lorentzian': lorentzian,
        'areal_radius': areal_radius,
        'disformal_ratio': disformal_ratio,
        'r_h': r_h, 'ell': ell, 'q': q,
    }


def compute_curvature_standard_coords_hayward(r, metric, M=1.0, ell=0.1):
    """Compute curvature of TEP matter metric on Hayward background.

    Uses standard coordinates (t, r) where the metric is:
      g̃_tt = -A²F + Bq²  (with q=0: -A²F)
      g̃_tr = Bqψ'         (with q=0: 0, diagonal)
      g̃_rr = A²/F + Bψ'²
      g̃_θθ = A²r²

    For q=0 (diagonal case), use the 4-term Kretschmann formula.
    For q≠0, we need the full off-diagonal computation.
    """
    r = np.asarray(r, dtype=float)
    A = metric['A']
    B = metric['B']
    A2 = metric['A2']
    F = metric['F']
    dphi = metric['dphi']
    q = metric['q']

    # Standard coordinates metric
    g_tt = -A2 * F + B * q ** 2
    g_tr = B * q * dphi  # off-diagonal
    g_rr = A2 / np.where(np.abs(F) > 1e-10, F, np.sign(F) * 1e-10) + B * dphi ** 2
    g_thth = A2 * r ** 2

    # For q=0 (diagonal), use the simple 4-term formula
    if abs(q) < 1e-10:
        return _curvature_diagonal(r, g_tt, g_rr, g_thth, M)
    else:
        return _curvature_offdiagonal(r, g_tt, g_tr, g_rr, g_thth, M)


def _curvature_diagonal(r, g_tt, g_rr, g_thth, M=1.0):
    """Compute curvature for diagonal spherically symmetric metric."""
    g_tt_safe = np.where(np.abs(g_tt) > 1e-50, g_tt, np.nan)
    g_rr_safe = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
    g_thth_safe = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

    ginv_tt = 1.0 / g_tt_safe
    ginv_rr = 1.0 / g_rr_safe
    ginv_thth = 1.0 / g_thth_safe

    g_tt_p = np.gradient(g_tt, r)
    g_rr_p = np.gradient(g_rr, r)
    g_thth_p = np.gradient(g_thth, r)

    Gt_tr = 0.5 * ginv_tt * g_tt_p
    Gr_rr = 0.5 * ginv_rr * g_rr_p
    Gth_rth = 0.5 * ginv_thth * g_thth_p

    dGt_tr = np.gradient(Gt_tr, r)
    dGth_rth = np.gradient(Gth_rth, r)

    R_trtr = g_tt * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr)
    R_rthrth = g_thth * (dGth_rth + Gth_rth ** 2 - Gth_rth * Gr_rr)
    R_tthtth = -0.25 * g_tt_p * g_thth_p / g_rr_safe
    R_thphthph = g_thth * (1.0 - ginv_rr * g_thth_p ** 2 / (4.0 * g_thth_safe))

    E = R_trtr / (g_tt_safe * g_rr_safe)
    F_t = R_tthtth / (g_tt_safe * g_thth_safe)
    F_r = R_rthrth / (g_rr_safe * g_thth_safe)
    G_comp = R_thphthph / g_thth_safe ** 2

    K = 4.0 * E ** 2 + 8.0 * F_t ** 2 + 8.0 * F_r ** 2 + 4.0 * G_comp ** 2

    g_thth_pp = np.gradient(g_thth_p, r)
    R_2D = 2.0 * R_trtr / (g_tt_safe * g_rr_safe)
    R_ang = -g_thth_pp / g_thth_safe + 0.5 * (g_thth_p / g_thth_safe) ** 2
    R_total = R_2D + 2.0 * R_ang / g_thth_safe

    return {'Kretschmann': K, 'Ricci_scalar': R_total,
            'R_trtr': R_trtr, 'R_rthrth': R_rthrth,
            'R_tthtth': R_tthtth, 'R_thphthph': R_thphthph}


def _curvature_offdiagonal(r, g_tt, g_tr, g_rr, g_thth, M=1.0):
    """Compute curvature for metric with off-diagonal g_tr term."""
    det_2d = g_tt * g_rr - g_tr ** 2
    det_2d_safe = np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    g_thth_safe = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)
    g_rr_safe = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)

    ginv_tt = g_rr / det_2d_safe
    ginv_tr = -g_tr / det_2d_safe
    ginv_rr = g_tt / det_2d_safe
    ginv_thth = 1.0 / g_thth_safe

    g_tt_p = np.gradient(g_tt, r)
    g_tr_p = np.gradient(g_tr, r)
    g_rr_p = np.gradient(g_rr, r)
    g_thth_p = np.gradient(g_thth, r)

    g_tt_pp = np.gradient(g_tt_p, r)
    g_thth_pp = np.gradient(g_thth_p, r)

    # Christoffel symbols
    Gt_tt = -0.5 * ginv_tr * g_tt_p
    Gr_tt = -0.5 * ginv_rr * g_tt_p
    Gt_tr = 0.5 * ginv_tt * g_tt_p
    Gr_tr = 0.5 * ginv_tr * g_tt_p
    Gt_rr = ginv_tt * g_tr_p + 0.5 * ginv_tr * g_rr_p
    Gr_rr = ginv_tr * g_tr_p + 0.5 * ginv_rr * g_rr_p
    Gt_thth = -0.5 * ginv_tr * g_thth_p
    Gr_thth = -0.5 * ginv_rr * g_thth_p
    Gth_rth = 0.5 * ginv_thth * g_thth_p

    # Riemann R_{trtr}
    R_trtr = (-0.5 * g_tt_pp
              + g_tt * (Gt_tt * Gr_rr - Gt_tr ** 2)
              + g_tr * (Gt_tt * Gr_rr - Gt_tr * Gr_tr)
              + g_tr * (Gr_tt * Gt_rr - Gr_tr * Gt_tr)
              + g_rr * (Gr_tt * Gr_rr - Gr_tr ** 2))

    R_rthrth = -0.5 * g_thth_pp + g_thth_p * Gr_thth

    # For off-diagonal, the simple 4-term formula needs modification
    # Use: K ~ 4*R_trtr²/det² + 4*R_rthrth²/(g_thth*g_rr) + angular terms
    K = (4.0 * R_trtr ** 2 / det_2d_safe ** 2
         + 4.0 * R_rthrth ** 2 / (g_thth_safe * g_rr_safe)
         + 4.0 * (g_thth * (1.0 - ginv_rr * g_thth_p ** 2 / (4.0 * g_thth_safe))) ** 2 / g_thth_safe ** 4)

    R_2D = 2.0 * R_trtr / det_2d_safe
    R_ang = -g_thth_pp / g_thth_safe + 0.5 * (g_thth_p / g_thth_safe) ** 2
    R_total = R_2D + 2.0 * R_ang / g_thth_safe

    return {'Kretschmann': K, 'Ricci_scalar': R_total,
            'R_trtr': R_trtr, 'R_rthrth': R_rthrth}


# =============================================================================
# GATE TESTS
# =============================================================================

def test_hayward_regular():
    """Test 1: Hayward metric itself is curvature-regular."""
    print("\n" + "=" * 70)
    print("TEST 1: Hayward metric curvature regularity")
    print("=" * 70)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    for ell in [0.01, 0.1, 0.5, 1.0]:
        F = hayward_F(r, M=1.0, ell=ell)
        r_h = find_hayward_horizon(M=1.0, ell=ell)

        # Compute Hayward curvature (A=1, B=0)
        g_tt = -F
        g_rr = 1.0 / np.where(np.abs(F) > 1e-10, F, np.sign(F) * 1e-10)
        g_thth = r ** 2

        curv = _curvature_diagonal(r, g_tt, g_rr, g_thth, M=1.0)
        K = curv['Kretschmann']

        K_at_0 = K[0] if np.isfinite(K[0]) else np.nan
        K_at_horizon = K[np.argmin(np.abs(r - r_h))]
        K_at_10 = K[np.argmin(np.abs(r - 10.0))]
        K_schw_10 = 48.0 / 10.0 ** 6

        print(f"\n  ℓ = {ell}:")
        print(f"    r_h = {r_h:.4f}")
        print(f"    K(r→0) = {K_at_0:.6e} (finite: {np.isfinite(K_at_0)})")
        print(f"    K(r_h) = {K_at_horizon:.6e}")
        print(f"    K(r=10) = {K_at_10:.6e} (Schw: {K_schw_10:.6e})")
        print(f"    K finite everywhere: {np.all(np.isfinite(K))}")


def test_tep_hayward_bounded():
    """Test 2: TEP on Hayward with bounded A — all gates."""
    print("\n" + "=" * 70)
    print("TEST 2: TEP matter metric on Hayward background")
    print("=" * 70)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    configs = [
        {'phi_0': 0.0, 'q': 0.0, 'ell': 0.1, 'label': 'A=1, q=0 (pure Hayward)'},
        {'phi_0': 0.5, 'q': 0.0, 'ell': 0.1, 'label': 'φ₀=0.5, q=0 (mild conformal)'},
        {'phi_0': 1.0, 'q': 0.0, 'ell': 0.1, 'label': 'φ₀=1.0, q=0 (bounded areal)'},
        {'phi_0': 0.5, 'q': 0.1, 'ell': 0.1, 'label': 'φ₀=0.5, q=0.1 (temporal disformal)'},
        {'phi_0': 1.0, 'q': 0.1, 'ell': 0.1, 'label': 'φ₀=1.0, q=0.1 (bounded + temporal)'},
        {'phi_0': 1.0, 'q': 0.5, 'ell': 0.1, 'label': 'φ₀=1.0, q=0.5 (strong temporal)'},
    ]

    print(f"\n  {'Config':>35} | {'Areal→0':>10} | {'K→0?':>6} | {'Lorentz?':>9} | {'K_max':>12}")
    print(f"  {'-'*35}-+-{'-'*10}-+-{'-'*6}-+-{'-'*9}-+-{'-'*12}")

    for cfg in configs:
        metric = compute_tep_hayward_metric(
            r, M=1.0, ell=cfg['ell'], phi_0=cfg['phi_0'],
            delta=0.05, beta_A=-1.0, B0=1.0, sigma_B=1.5, q=cfg['q'])

        curv = compute_curvature_standard_coords_hayward(r, metric, M=1.0, ell=cfg['ell'])
        K = curv['Kretschmann']

        areal = metric['areal_radius']
        areal_0 = areal[0] if np.isfinite(areal[0]) else np.nan
        lorentz_frac = np.mean(metric['lorentzian'])
        K_max = np.nanmax(K[np.isfinite(K)]) if np.any(np.isfinite(K)) else np.nan

        # K vanishes if K at small r << K at r=1
        K_small = K[np.argmin(np.abs(r - 1e-4))]
        K_ref = K[np.argmin(np.abs(r - 1.0))]
        K_vanishes = K_small < K_ref * 1e-3 if np.isfinite(K_small) and np.isfinite(K_ref) else False

        print(f"  {cfg['label']:>35} | {areal_0:>10.4f} | {'YES' if K_vanishes else 'NO':>6} | "
              f"{lorentz_frac*100:>7.1f}% | {K_max:>12.4e}")


def test_tep_hayward_detailed():
    """Test 3: Detailed analysis of best candidate configuration."""
    print("\n" + "=" * 70)
    print("TEST 3: Detailed analysis — φ₀=1.0, q=0.1 on Hayward (ℓ=0.1)")
    print("=" * 70)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)
    metric = compute_tep_hayward_metric(
        r, M=1.0, ell=0.1, phi_0=1.0, delta=0.05,
        beta_A=-1.0, B0=1.0, sigma_B=1.5, q=0.1)

    curv = compute_curvature_standard_coords_hayward(r, metric, M=1.0, ell=0.1)
    K = curv['Kretschmann']
    R = curv['Ricci_scalar']

    areal = metric['areal_radius']
    r_h = metric['r_h']

    print(f"\n  Hayward horizon: r_h = {r_h:.4f}")
    print(f"  Regularization scale: ℓ = {metric['ell']}")

    print(f"\n  {'r':>10} | {'A':>10} | {'B':>10} | {'ρ=A·r':>10} | {'K':>12} | {'det_2d':>12} | {'Lor?':>5}")
    print(f"  {'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*12}-+-{'-'*12}-+-{'-'*5}")

    for r_test in [50.0, 10.0, r_h, 1.0, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-6, 1e-8]:
        idx = np.argmin(np.abs(r - r_test))
        lor = metric['lorentzian'][idx]
        K_val = K[idx] if np.isfinite(K[idx]) else float('nan')
        print(f"  {r_test:>10.2e} | {metric['A'][idx]:>10.4e} | {metric['B'][idx]:>10.4e} | "
              f"{areal[idx]:>10.4e} | {K_val:>12.4e} | {metric['det_2d'][idx]:>12.4e} | {'Y' if lor else 'N':>5}")

    # Power law of K near r=0
    mask = (r >= 1e-6) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
    if np.sum(mask) > 10:
        alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
        print(f"\n  K power law (r in [1e-6, 1e-2]): K ~ r^{{{alpha:.3f}}}")

    # Areal radius behavior
    print(f"\n  Areal radius:")
    print(f"    ρ(r_h) = {areal[np.argmin(np.abs(r - r_h))]:.6f}")
    print(f"    ρ(0.1) = {areal[np.argmin(np.abs(r - 0.1))]:.6f}")
    print(f"    ρ(0.01) = {areal[np.argmin(np.abs(r - 0.01))]:.6f}")
    print(f"    ρ(r→0) = {areal[0]:.6f}")
    print(f"    ρ bounded: {areal[0] < 1e3}")

    # Lorentzian fraction
    lor_frac = np.mean(metric['lorentzian'])
    print(f"\n  Lorentzian fraction: {lor_frac*100:.2f}%")
    if lor_frac < 1.0:
        # Find where it fails
        fail_idx = np.where(~metric['lorentzian'])[0]
        if len(fail_idx) > 0:
            print(f"    First failure at r = {r[fail_idx[0]]:.6e}")
            print(f"    det_2d at failure = {metric['det_2d'][fail_idx[0]]:.6e}")


def main():
    print("STAGE A: TEP MATTER METRIC ON HAYWARD REGULAR BACKGROUND")
    print("Testing whether the temporal picture is possible with regular geometry")

    # First verify Hayward is regular
    hayward_kretschmann_analytical()

    # Test Hayward numerically
    test_hayward_regular()

    # Test TEP on Hayward — overview
    test_tep_hayward_bounded()

    # Detailed analysis of best candidate
    test_tep_hayward_detailed()


if __name__ == '__main__':
    main()
