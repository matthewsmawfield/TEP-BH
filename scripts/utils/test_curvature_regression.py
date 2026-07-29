#!/usr/bin/env python3
"""Regression test: verify corrected curvature computation against exact formula.

Tests:
1. Schwarzschild (A=1, B=0): K = 48*M^2/r^6
2. Conformal phi_0=1 (deep interior): K ~ 9/(M^2*r^2) + 1/(4*M^4) -> infinity
3. Conformal phi_0=2 (deep interior): K ~ 39*r^2/(16*M^6) -> 0
4. Power law verification: K ~ r^{4*phi_0 - 6} for various phi_0
"""

import sys
import os
import numpy as np

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)

from scripts.steps.bh_common import TEPBHModel, compute_disformal_metric, compute_curvature_invariants


def test_schwarzschild():
    """Test: Schwarzschild metric (A=1, B=0) gives K = 48*M^2/r^6."""
    print("=" * 60)
    print("TEST 1: Schwarzschild K = 48*M^2/r^6")
    print("=" * 60)

    # Use phi_0=0 to get A=1, B=0 (Schwarzschild)
    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=0.001, delta=0.001, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(2.5), np.log10(50.0), 5000)  # exterior only
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']
    K_schw = curvature['Kretschmann_schwarzschild']

    # Check at several radii (exclude very large r where K is tiny and numerical noise dominates)
    max_ratio = 0
    for r_test in [3.0, 5.0, 10.0, 20.0, 30.0]:
        idx = np.argmin(np.abs(r - r_test))
        ratio = K[idx] / K_schw[idx] if K_schw[idx] > 0 else 0
        max_ratio = max(max_ratio, abs(ratio - 1.0))
        print(f"  r={r_test:.0f}: K={K[idx]:.6e}, K_Schw={K_schw[idx]:.6e}, ratio={ratio:.6f}")

    passed = max_ratio < 0.01  # 1% tolerance
    print(f"  Max deviation from 48M²/r⁶: {max_ratio:.4e}")
    print(f"  {'✓ PASS' if passed else '✗ FAIL'}")
    return passed


def test_conformal_phi0_1():
    """Test: phi_0=1 conformal gives K ~ 9/(M^2*r^2) in deep interior (DIVERGES)."""
    print("\n" + "=" * 60)
    print("TEST 2: phi_0=1 conformal, K ~ 9/(M²r²) -> DIVERGES")
    print("=" * 60)

    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=1.0, delta=0.05, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(1e-6), np.log10(1.5), 50000)  # deep interior
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']

    # Check power law: K ~ r^{4*1-6} = r^{-2}
    # Use range [1e-4, 1e-2] to avoid boundary effects at very small/large r
    mask = (r >= 1e-4) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
    if np.sum(mask) > 10:
        alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
    else:
        alpha = np.nan

    print(f"  Power law (r in [1e-4, 1e-2]): K ~ r^{{{alpha:.3f}}} (expected: r^{{-2}})")
    print(f"  K at r=0.1: {K[np.argmin(np.abs(r - 0.1))]:.6e}")
    print(f"  K at r=0.01: {K[np.argmin(np.abs(r - 0.01))]:.6e}")
    print(f"  K at r=0.001: {K[np.argmin(np.abs(r - 0.001))]:.6e}")

    # K should diverge (alpha < -1.5)
    passed = alpha < -1.5
    print(f"  K diverges: {'YES' if alpha < -1.5 else 'NO'}")
    print(f"  {'✓ PASS' if passed else '✗ FAIL'}")
    return passed


def test_conformal_phi0_2():
    """Test: phi_0=2 conformal gives K ~ r^2 -> 0 in deep interior (VANISHES)."""
    print("\n" + "=" * 60)
    print("TEST 3: phi_0=2 conformal, K ~ r^2 -> VANISHES")
    print("=" * 60)

    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=2.0, delta=0.05, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(1e-6), np.log10(1.5), 50000)
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']

    # Check power law: K ~ r^{4*2-6} = r^2
    # Use range [1e-4, 1e-2] to avoid boundary effects
    mask = (r >= 1e-4) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
    if np.sum(mask) > 10:
        alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
    else:
        alpha = np.nan

    print(f"  Power law (r in [1e-4, 1e-2]): K ~ r^{{{alpha:.3f}}} (expected: r^{{2}})")
    print(f"  K at r=0.1: {K[np.argmin(np.abs(r - 0.1))]:.6e}")
    print(f"  K at r=0.01: {K[np.argmin(np.abs(r - 0.01))]:.6e}")
    print(f"  K at r=0.001: {K[np.argmin(np.abs(r - 0.001))]:.6e}")

    # K should vanish (alpha > 0)
    passed = alpha > 0.5
    print(f"  K vanishes: {'YES' if alpha > 0.5 else 'NO'}")
    print(f"  {'✓ PASS' if passed else '✗ FAIL'}")
    return passed


def test_power_law_scan():
    """Test: power law K ~ r^{4*phi_0-6} for various phi_0."""
    print("\n" + "=" * 60)
    print("TEST 4: Power law K ~ r^{4*phi_0-6} for various phi_0")
    print("=" * 60)

    results = []
    for phi_0 in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]:
        model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=phi_0, delta=0.05, M=1.0, sigma_B=1.5)
        r = np.logspace(np.log10(1e-6), np.log10(1.5), 50000)
        metric = compute_disformal_metric(r, model)
        curvature = compute_curvature_invariants(r, metric, M=1.0)

        K = curvature['Kretschmann']
        # Use range [1e-4, 1e-2] for power law fit (avoids boundary effects)
        mask = (r >= 1e-4) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
        else:
            alpha = np.nan

        expected = 4 * phi_0 - 6
        diff = abs(alpha - expected)
        results.append((phi_0, alpha, expected, diff))

        areal = metric['A'] * r
        areal_div = areal[0] > 1e3

        print(f"  phi_0={phi_0:.2f}: K ~ r^{{{alpha:.2f}}}, expected r^{{{expected:.1f}}}, "
              f"diff={diff:.2f}, areal={'diverges' if areal_div else 'finite'}")

    # All power laws should match within tolerance
    max_diff = max(r[3] for r in results)
    passed = max_diff < 0.1  # tolerance for numerical effects
    print(f"\n  Max power law deviation: {max_diff:.4f}")
    print(f"  {'✓ PASS' if passed else '✗ FAIL'}")
    return passed


def test_no_go_theorem():
    """Test: no phi_0 gives both finite areal radius AND vanishing K."""
    print("\n" + "=" * 60)
    print("TEST 5: No-go theorem (finite areal + vanishing K = impossible)")
    print("=" * 60)

    found_both = False
    for phi_0 in np.arange(0.5, 3.01, 0.25):
        model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=phi_0, delta=0.05, M=1.0, sigma_B=1.5)
        r = np.logspace(np.log10(1e-6), np.log10(1.5), 50000)
        metric = compute_disformal_metric(r, model)
        curvature = compute_curvature_invariants(r, metric, M=1.0)

        K = curvature['Kretschmann']
        areal = metric['A'] * r

        # K vanishes if K at small r << K at r=0.5
        idx_small = np.argmin(np.abs(r - 1e-3))
        idx_ref = np.argmin(np.abs(r - 0.5))
        K_vanishes = K[idx_small] < K[idx_ref] * 1e-3
        areal_finite = areal[0] < 1e3

        if K_vanishes and areal_finite:
            found_both = True
            print(f"  phi_0={phi_0:.2f}: K_vanishes={K_vanishes}, areal_finite={areal_finite} *** BOTH! ***")
        else:
            status = ""
            if K_vanishes:
                status += "K->0 "
            if areal_finite:
                status += "areal_finite "
            if not status:
                status = "neither"
            print(f"  phi_0={phi_0:.2f}: K_vanishes={K_vanishes}, areal_finite={areal_finite} ({status})")

    passed = not found_both
    print(f"\n  No phi_0 achieves both: {'YES' if passed else 'NO'}")
    print(f"  {'✓ PASS (no-go confirmed)' if passed else '✗ FAIL (no-go violated!)'}")
    return passed


def main():
    print("CURVATURE REGRESSION TESTS")
    print("Verifying corrected compute_curvature_invariants against exact formula")
    print()

    tests = [
        test_schwarzschild,
        test_conformal_phi0_1,
        test_conformal_phi0_2,
        test_power_law_scan,
        test_no_go_theorem,
    ]

    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append(False)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for i, (test, result) in enumerate(zip(tests, results)):
        name = test.__name__.replace('test_', '').replace('_', ' ')
        print(f"  {i+1}. {name}: {'PASS' if result else 'FAIL'}")

    all_pass = all(results)
    print(f"\n  Overall: {'ALL PASS ✓' if all_pass else 'SOME FAIL ✗'}")
    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
