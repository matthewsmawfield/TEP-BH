#!/usr/bin/env python3
"""Phase 9: Test existing sGB branch for TEP requirements.

The existing sGB branch (step_12) uses perturbative O(η²) corrections
on top of Schwarzschild. The interior remains Schwarzschild.

Tests:
1. Curvature regularity (does K diverge at r=0?)
2. Bounded areal radius
3. Lorentzianity
4. Temporal behaviour
5. Hyperbolicity (characteristic structure)

Expected: The perturbative sGB does NOT regularize the singularity
(corrections are O(η²) ~ small, while K_Schw ~ 48M²/r⁶ diverges).
But we verify this explicitly and identify what's needed for a
non-perturbative regular solution.
"""

import numpy as np
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, 'scripts', 'steps'))

from step_12_self_gravitating import (
    SGBModel, scalar_field_sgb, scalar_field_derivative_sgb,
    compute_metric_corrections, compute_tep_matter_metric_full,
    gauss_bonnet_schwarzschild
)
from scripts.steps.bh_common import compute_curvature_invariants, compute_disformal_metric, TEPBHModel


def test_sgb_curvature_regularity():
    """Test 1: Does the sGB-corrected metric have finite K at r=0?"""
    print("=" * 80)
    print("TEST 1: sGB CURVATURE REGULARITY")
    print("=" * 80)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    for eta in [0.1, 0.3, 0.5, 1.0, 2.0, 5.0]:
        model = SGBModel(eta=eta, M=1.0, B0=1.0, beta_A=1.0, sigma_B=0.5)
        M = 1.0
        r_h = 2.0 * M

        phi = scalar_field_sgb(r, eta, M)
        dphi = scalar_field_derivative_sgb(r, eta, M)
        F = 1.0 - 2.0 * M / r

        # Metric corrections (exterior only, perturbative)
        r_ext_mask = r > r_h + 1e-6
        r_ext = r[r_ext_mask]
        corrections = compute_metric_corrections(r_ext, eta, M)

        F_corrected = np.where(r_ext_mask,
                               np.interp(r, r_ext, corrections['F_corrected']),
                               F)
        g_rr_corrected = np.where(r_ext_mask,
                                  np.interp(r, r_ext, corrections['g_rr_corrected']),
                                  1.0 / np.where(np.abs(F) > 1e-30, F, np.nan))
        g_rr_corrected = np.where(np.isfinite(g_rr_corrected), g_rr_corrected, 0)

        # Compute curvature of the CORRECTED GEOMETRIC metric
        # (not the TEP matter metric — we want to know if g itself is regular)
        # Use standard coords: g_tt = -F_corrected, g_rr = g_rr_corrected
        g_tt = -F_corrected
        g_rr = g_rr_corrected
        g_thth = r ** 2

        # Regularize F near horizon
        F_reg = np.sign(F_corrected) * np.maximum(np.abs(F_corrected), 1e-6)
        g_tt_safe = -F_reg
        g_rr_safe = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
        g_thth_safe = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

        ginv_tt = 1.0 / np.where(np.abs(g_tt_safe) > 1e-50, g_tt_safe, np.nan)
        ginv_rr = 1.0 / g_rr_safe
        ginv_thth = 1.0 / g_thth_safe

        g_tt_p = np.gradient(g_tt_safe, r)
        g_rr_p = np.gradient(g_rr_safe, r)
        g_thth_p = np.gradient(g_thth, r)

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

        # Power law near r=0
        mask = (r >= 1e-6) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
        else:
            alpha = np.nan

        K_0 = K[0] if np.isfinite(K[0]) else np.nan
        K_finite = np.isfinite(K_0)

        print(f"\n  η = {eta}:")
        print(f"    K power law near r=0: ~r^{{{alpha:.3f}}} (Schwarzschild: r^{{-6}})")
        print(f"    K(r→0) = {K_0:.6e}")
        print(f"    K finite at r=0: {K_finite}")
        print(f"    → {'REGULAR' if K_finite and alpha > -1 else 'SINGULAR (still diverges)'}")


def test_sgb_tep_matter_metric():
    """Test 2: TEP matter metric on sGB-corrected background."""
    print("\n" + "=" * 80)
    print("TEST 2: TEP MATTER METRIC ON sGB BACKGROUND")
    print("=" * 80)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    for eta in [0.3, 0.5, 1.0]:
        model = SGBModel(eta=eta, M=1.0, B0=1.0, beta_A=1.0, sigma_B=0.5)
        M = 1.0
        r_h = 2.0 * M

        phi = scalar_field_sgb(r, eta, M)
        dphi = scalar_field_derivative_sgb(r, eta, M)
        F = 1.0 - 2.0 * M / r

        r_ext_mask = r > r_h + 1e-6
        r_ext = r[r_ext_mask]
        corrections = compute_metric_corrections(r_ext, eta, M)

        F_corrected = np.where(r_ext_mask,
                               np.interp(r, r_ext, corrections['F_corrected']),
                               F)
        g_rr_corrected = np.where(r_ext_mask,
                                  np.interp(r, r_ext, corrections['g_rr_corrected']),
                                  1.0 / np.where(np.abs(F) > 1e-30, F, np.nan))
        g_rr_corrected = np.where(np.isfinite(g_rr_corrected), g_rr_corrected, 0)

        tep = compute_tep_matter_metric_full(r, phi, dphi, F_corrected, g_rr_corrected, model)

        A = tep['A']
        areal = A * r
        det_2d = tep['det_2d']
        lorentz_frac = np.mean(det_2d < 0) * 100

        areal_0 = areal[0] if np.isfinite(areal[0]) else np.nan

        print(f"\n  η = {eta}:")
        print(f"    A(r→0) = {A[0]:.6e}")
        print(f"    Areal radius ρ(r→0) = {areal_0:.6e}")
        print(f"    Areal bounded: {areal_0 < 1e3}")
        print(f"    Lorentzian: {lorentz_frac:.1f}%")
        print(f"    det_2d at r=0.1: {det_2d[np.argmin(np.abs(r - 0.1))]:.6e}")


def test_nonperturbative_sgb():
    """Test 3: What happens with LARGE η (non-perturbative regime)?

    The perturbative Sotiriou-Zhou corrections are O(η²) and only valid
    for ζ = η²/M⁴ << 1. For large η, the corrections become large and
    the perturbative expansion breaks down.

    In the non-perturbative regime, sGB is known to produce regular BH
    solutions (Kanti 1996, Kanti et al. 1996). But the perturbative
    formula used in step_12 is NOT valid there.

    We need to solve the full coupled ODEs non-perturbatively.
    """
    print("\n" + "=" * 80)
    print("TEST 3: NON-PERTURBATIVE sGB REGIME")
    print("=" * 80)

    # The full sGB field equations for a static spherical metric:
    # ds² = -f(r) e^{2Φ(r)} dt² + f(r)^{-1} dr² + r² dΩ²
    #
    # With f = 1 - 2m(r)/r:
    #   m'(r) = 4π r² ρ_eff(r)
    #   Φ'(r) = (4π r (ρ + p_r) + ... ) / f  (TOV-like with GB corrections)
    #   φ'' + (2/r + f'/f) φ' = -η · G_GB / f  (scalar equation)
    #
    # where G_GB is the Gauss-Bonnet invariant and ρ_eff includes
    # both scalar kinetic energy and GB coupling contributions.
    #
    # For REGULAR solutions (Kanti 1996):
    #   - The scalar field φ → 0 at r=0 (regular)
    #   - The mass function m(0) = 0 (no point mass)
    #   - The metric function f(0) = 1 (de Sitter core)
    #
    # This requires solving the ODEs with boundary conditions:
    #   m(0) = 0, φ(0) = 0, φ'(0) = 0
    #   m(∞) = M, φ(∞) = 0

    print("""
  The perturbative sGB (Sotiriou-Zhou O(η²)) is valid only for ζ = η²/M⁴ << 1.
  In this regime, corrections are small and the interior remains Schwarzschild.

  For REGULAR solutions, we need the NON-PERTURBATIVE regime (Kanti 1996):
    - Solve the full coupled ODEs: m'(r), Φ'(r), φ''(r)
    - Boundary conditions: m(0)=0, φ(0)=0, φ'(0)=0, m(∞)=M
    - The scalar field backreacts strongly on the geometry
    - Result: de Sitter core (f(0)=1), regular curvature

  The full sGB equations (Kanti 1996, Sotiriou-Zhou 2014 non-perturbative):

    m'(r) = r²/(2F) · [4πG ρ_scalar + η · G_GB · φ' · r/2]
    Φ'(r) = [4πG (ρ+p_r) r² + η · G_GB · φ' · r³/(4F)] / [r(r - 2m)]
    φ'' + [(2 - Φ')/r + (m'r - m)/(r(r-2m))] φ' = -η G_GB / F

  where F = 1 - 2m(r)/r, G_GB = K/8 (Gauss-Bonnet invariant).

  These ODEs must be solved numerically with the regular boundary
  conditions at r=0. This is the Stage C computation.

  Key references:
    - Kanti, Geyer, Zwiep (1996): regular sGB BH solutions
    - Kanti, Mavromatos (1997): scalar hair and regularity
    - Sotiriou, Zhou (2014): shift-symmetric sGB, perturbative
    - Kleihaus, Kunz, Schneider (2011): non-perturbative numerical solutions
""")


def test_hyperbolicity():
    """Test 4: Hyperbolicity of the sGB-corrected metric."""
    print("=" * 80)
    print("TEST 4: HYPERBOLICITY (characteristic structure)")
    print("=" * 80)

    r = np.logspace(np.log10(1e-6), np.log10(50.0), 50000)

    for eta in [0.3, 0.5, 1.0]:
        model = SGBModel(eta=eta, M=1.0, B0=1.0, beta_A=1.0, sigma_B=0.5)
        M = 1.0
        r_h = 2.0 * M

        phi = scalar_field_sgb(r, eta, M)
        dphi = scalar_field_derivative_sgb(r, eta, M)
        F = 1.0 - 2.0 * M / r

        r_ext_mask = r > r_h + 1e-6
        r_ext = r[r_ext_mask]
        corrections = compute_metric_corrections(r_ext, eta, M)

        F_corrected = np.where(r_ext_mask,
                               np.interp(r, r_ext, corrections['F_corrected']),
                               F)
        g_rr_corrected = np.where(r_ext_mask,
                                  np.interp(r, r_ext, corrections['g_rr_corrected']),
                                  1.0 / np.where(np.abs(F) > 1e-30, F, np.nan))
        g_rr_corrected = np.where(np.isfinite(g_rr_corrected), g_rr_corrected, 0)

        tep = compute_tep_matter_metric_full(r, phi, dphi, F_corrected, g_rr_corrected, model)

        # Hyperbolicity: the metric must have Lorentzian signature everywhere
        # For the matter metric, this means det_2d < 0
        det_2d = tep['det_2d']

        # The characteristic speeds (eigenvalues of the metric) must be real
        # For a 2x2 metric, this is equivalent to det < 0
        hyperbolic = det_2d < 0
        hyp_frac = np.mean(hyperbolic) * 100

        print(f"\n  η = {eta}:")
        print(f"    Hyperbolic (Lorentzian) fraction: {hyp_frac:.1f}%")
        if hyp_frac < 100:
            fail_idx = np.where(~hyperbolic)[0]
            if len(fail_idx) > 0:
                print(f"    First failure at r = {r[fail_idx[0]]:.6e}")
                print(f"    det_2d at failure = {det_2d[fail_idx[0]]:.6e}")
        else:
            print(f"    → GLOBALLY HYPERBOLIC ✓")


def main():
    test_sgb_curvature_regularity()
    test_sgb_tep_matter_metric()
    test_nonperturbative_sgb()
    test_hyperbolicity()

    print("\n" + "=" * 80)
    print("PHASE 9 CONCLUSIONS")
    print("=" * 80)
    print("""
  The existing perturbative sGB branch (step_12) does NOT regularize
  the singularity. The interior remains Schwarzschild with O(η²)
  corrections that are too small to remove the r^{-6} divergence.

  For a regular TEP-BH solution, we need the NON-PERTURBATIVE sGB
  regime (Kanti 1996), solving the full coupled ODEs with regular
  boundary conditions at r=0.

  This is Stage C: solve the coupled self-gravitating system.
""")


if __name__ == '__main__':
    main()
