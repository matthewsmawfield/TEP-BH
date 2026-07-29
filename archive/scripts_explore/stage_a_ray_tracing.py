#!/usr/bin/env python3
"""Phase 6d: Two-way ray tracing on TEP-Hayward regular background.

Tests the optical properties of the Stage A configuration:
  - Inward photons: do they reach r=0 in finite affine parameter?
  - Outward photons: can they escape from the interior?
  - Temporal freeze: does dv/dr → 0 (asymptotic optical boundary)?
  - Redshift: how does the temporal distortion affect photon frequencies?

The Hayward regular background has F(0) = 1 (no singularity), so photons
can in principle reach r=0 and continue through. The TEP matter metric
adds temporal distortion via the disformal term.

Key questions:
  1. Are null geodesics complete (affine parameter → ∞ or finite)?
  2. Is there a "temporal horizon" where outgoing photons freeze?
  3. What is the redshift profile for infalling vs outgoing photons?
"""

import numpy as np
import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)

from scripts.explore.stage_a_hayward_benchmark import (
    hayward_F, find_hayward_horizon, scalar_profile_bounded,
    scalar_gradient_bounded, compute_tep_hayward_metric,
    compute_curvature_standard_coords_hayward
)


def null_geodesic_integrand(r, metric):
    """Compute the null geodesic integrand dv/dr for the matter metric.

    For EF coordinates with metric:
      g̃_vv dv² + 2 g̃_vr dv dr + g̃_rr dr² = 0

    Null condition: g̃_vv (dv/dr)² + 2 g̃_vr (dv/dr) + g̃_rr = 0

    Solutions:
      (dv/dr)_ingoing  = (-g̃_vr - sqrt(g̃_vr² - g̃_vv g̃_rr)) / g̃_vv
      (dv/dr)_outgoing = (-g̃_vr + sqrt(g̃_vr² - g̃_vv g̃_rr)) / g̃_vv

    The affine parameter integrand is related to dv/dr.
    For completeness, we check if ∫ dr / |dv/dr| converges or diverges.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']

    discriminant = gvr ** 2 - gvv * grr
    disc_safe = np.where(discriminant > 0, discriminant, np.nan)
    sqrt_disc = np.sqrt(disc_safe)

    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)

    dv_dr_in = (-gvr - sqrt_disc) / gvv_safe
    dv_dr_out = (-gvr + sqrt_disc) / gvv_safe

    return dv_dr_in, dv_dr_out


def compute_affine_parameter(r, dv_dr, r_start, r_end):
    """Compute affine parameter integral ∫ dr / (dv/dr) from r_start to r_end.

    If the integral diverges, the geodesic is complete (never reaches r_end).
    If it converges, the geodesic reaches r_end in finite affine parameter.
    """
    # Sort r in the direction of travel
    if r_end > r_start:
        mask = (r >= r_start) & (r <= r_end)
    else:
        mask = (r <= r_start) & (r >= r_end)

    r_seg = r[mask]
    dv_dr_seg = dv_dr[mask]

    if len(r_seg) < 2:
        return np.nan, False

    # Integrands: 1/|dv/dr| for affine parameter
    # (The exact relation depends on normalization, but convergence/divergence
    #  is what matters for completeness)
    integrand = 1.0 / np.abs(dv_dr_seg)

    # Trapezoidal integration
    dr = np.abs(np.diff(r_seg))
    integrand_mid = 0.5 * (integrand[:-1] + integrand[1:])
    integral = np.sum(dr * integrand_mid)

    # Check for divergence (if integrand grows too fast near r_end)
    finite = np.isfinite(integral)
    diverges = integral > 1e10 if finite else True

    return integral, diverges


def test_null_geodesics():
    """Test null geodesic structure on TEP-Hayward background."""
    print("=" * 80)
    print("PHASE 6d: NULL GEODESICS ON TEP-HAYWARD BACKGROUND")
    print("=" * 80)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    configs = [
        {'phi_0': 0.0, 'q': 0.0, 'label': 'Pure Hayward (no TEP)'},
        {'phi_0': 1.0, 'q': 0.0, 'label': 'φ₀=1.0, q=0 (bounded areal)'},
        {'phi_0': 1.0, 'q': 0.1, 'label': 'φ₀=1.0, q=0.1 (temporal)'},
        {'phi_0': 1.25, 'q': 0.0, 'label': 'φ₀=1.25, q=0 (K→0)'},
    ]

    for cfg in configs:
        print(f"\n{'='*80}")
        print(f"  Configuration: {cfg['label']}")
        print(f"{'='*80}")

        metric = compute_tep_hayward_metric(
            r, M=1.0, ell=0.1, phi_0=cfg['phi_0'], delta=0.05,
            beta_A=-1.0, B0=1.0, sigma_B=1.5, q=cfg['q'])

        r_h = metric['r_h']
        dv_dr_in, dv_dr_out = null_geodesic_integrand(r, metric)

        print(f"\n  Null geodesic slopes dv/dr:")
        print(f"  {'r':>10} | {'dv/dr_in':>14} | {'dv/dr_out':>14} | {'freeze?':>8}")
        print(f"  {'-'*10}-+-{'-'*14}-+-{'-'*14}-+-{'-'*8}")

        for r_test in [50.0, 10.0, r_h, 1.0, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-6, 1e-8]:
            idx = np.argmin(np.abs(r - r_test))
            dv_in = dv_dr_in[idx] if np.isfinite(dv_dr_in[idx]) else float('nan')
            dv_out = dv_dr_out[idx] if np.isfinite(dv_dr_out[idx]) else float('nan')
            freeze = np.abs(dv_out) < 1e-10 if np.isfinite(dv_out) else False
            print(f"  {r_test:>10.2e} | {dv_in:>14.6e} | {dv_out:>14.6e} | {'YES' if freeze else '':>8}")

        # Affine parameter: horizon to center (inward)
        lambda_in, div_in = compute_affine_parameter(r, dv_dr_in, r_h, 1e-6)
        print(f"\n  Inward affine parameter (r_h → 1e-6): {lambda_in:.4e} {'(DIVERGES = complete)' if div_in else '(FINITE)'}")

        # Affine parameter: center to horizon (outward)
        lambda_out, div_out = compute_affine_parameter(r, dv_dr_out, 1e-6, r_h)
        print(f"  Outward affine parameter (1e-6 → r_h): {lambda_out:.4e} {'(DIVERGES = complete)' if div_out else '(FINITE)'}")

        # Affine parameter: horizon to infinity (outward)
        lambda_esc, div_esc = compute_affine_parameter(r, dv_dr_out, r_h, 50.0)
        print(f"  Escape affine parameter (r_h → 50): {lambda_esc:.4e} {'(DIVERGES)' if div_esc else '(FINITE = can escape)'}")

        # Power law of outgoing dv/dr near r=0
        mask = (r >= 1e-6) & (r <= 1e-2) & np.isfinite(dv_dr_out) & (np.abs(dv_dr_out) > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(np.abs(dv_dr_out[mask])), 1)[0]
            print(f"  Outgoing dv/dr power law: ~r^{{{alpha:.3f}}}")
            if alpha > 0:
                print(f"    → dv/dr → 0 as r → 0 (TEMPORAL FREEZE — outgoing photons frozen)")
            elif alpha < 0:
                print(f"    → dv/dr → ∞ as r → 0 (photons can escape)")
            else:
                print(f"    → dv/dr → const (neutral)")


def test_redshift_profile():
    """Compute redshift profile for infalling and outgoing photons."""
    print("\n" + "=" * 80)
    print("REDSHIFT PROFILE ON TEP-HAYWARD BACKGROUND")
    print("=" * 80)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    for phi_0 in [1.0, 1.25]:
        for q in [0.0, 0.1]:
            metric = compute_tep_hayward_metric(
                r, M=1.0, ell=0.1, phi_0=phi_0, delta=0.05,
                beta_A=-1.0, B0=1.0, sigma_B=1.5, q=q)

            # The redshift factor for a static emitter at radius r
            # observed at infinity is related to sqrt(-g̃_vv / g̃_vv(∞))
            # For EF: g̃_vv(∞) → -1 (asymptotically flat)
            gvv = metric['gtilde_vv']
            gvv_inf = -1.0  # asymptotic value

            # Redshift 1+z = sqrt(gvv_inf / gvv) = sqrt(1/|gvv|)
            # (for gvv < 0 inside horizon-like region)
            redshift = np.sqrt(np.abs(gvv_inf) / np.abs(gvv))

            print(f"\n  φ₀={phi_0}, q={q}:")
            print(f"  {'r':>10} | {'g̃_vv':>14} | {'1+z':>14} | {'log(1+z)':>12}")
            print(f"  {'-'*10}-+-{'-'*14}-+-{'-'*14}-+-{'-'*12}")

            for r_test in [50.0, 10.0, metric['r_h'], 1.0, 0.5, 0.1, 0.01, 0.001, 1e-6]:
                idx = np.argmin(np.abs(r - r_test))
                gvv_val = gvv[idx]
                z_val = redshift[idx] if np.isfinite(redshift[idx]) else float('nan')
                log_z = np.log10(z_val) if z_val > 0 and np.isfinite(z_val) else float('nan')
                print(f"  {r_test:>10.2e} | {gvv_val:>14.6e} | {z_val:>14.6e} | {log_z:>12.4f}")


def test_geodesic_completeness_summary():
    """Summary of geodesic completeness for all configurations."""
    print("\n" + "=" * 80)
    print("GEODESIC COMPLETENESS SUMMARY")
    print("=" * 80)

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 100000)

    print(f"\n  {'Config':>30} | {'Null in':>10} | {'Null out':>10} | {'Freeze?':>8} | {'Esc?':>6}")
    print(f"  {'-'*30}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-{'-'*6}")

    configs = [
        {'phi_0': 0.0, 'q': 0.0, 'label': 'Pure Hayward'},
        {'phi_0': 0.5, 'q': 0.0, 'label': 'φ₀=0.5, q=0'},
        {'phi_0': 1.0, 'q': 0.0, 'label': 'φ₀=1.0, q=0'},
        {'phi_0': 1.0, 'q': 0.1, 'label': 'φ₀=1.0, q=0.1'},
        {'phi_0': 1.0, 'q': 0.5, 'label': 'φ₀=1.0, q=0.5'},
        {'phi_0': 1.25, 'q': 0.0, 'label': 'φ₀=1.25, q=0'},
        {'phi_0': 1.25, 'q': 0.1, 'label': 'φ₀=1.25, q=0.1'},
    ]

    for cfg in configs:
        metric = compute_tep_hayward_metric(
            r, M=1.0, ell=0.1, phi_0=cfg['phi_0'], delta=0.05,
            beta_A=-1.0, B0=1.0, sigma_B=1.5, q=cfg['q'])

        dv_dr_in, dv_dr_out = null_geodesic_integrand(r, metric)
        r_h = metric['r_h']

        # Inward: horizon to center
        lambda_in, div_in = compute_affine_parameter(r, dv_dr_in, r_h, 1e-6)
        # Outward: center to horizon
        lambda_out, div_out = compute_affine_parameter(r, dv_dr_out, 1e-6, r_h)
        # Escape: horizon to infinity
        lambda_esc, div_esc = compute_affine_parameter(r, dv_dr_out, r_h, 50.0)

        # Freeze: does dv/dr_out → 0 near r=0?
        idx_small = np.argmin(np.abs(r - 1e-4))
        freeze = np.abs(dv_dr_out[idx_small]) < 1e-10 if np.isfinite(dv_dr_out[idx_small]) else False

        in_status = "complete" if div_in else "finite"
        out_status = "complete" if div_out else "finite"
        esc_status = "YES" if not div_esc else "NO"

        print(f"  {cfg['label']:>30} | {in_status:>10} | {out_status:>10} | "
              f"{'YES' if freeze else 'NO':>8} | {esc_status:>6}")


def main():
    test_null_geodesics()
    test_redshift_profile()
    test_geodesic_completeness_summary()

    print("\n" + "=" * 80)
    print("RAY TRACING CONCLUSIONS")
    print("=" * 80)


if __name__ == '__main__':
    main()
