#!/usr/bin/env python3
"""Debug: check numerical curvature at specific radii vs analytical predictions."""

import sys
import os
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)

from scripts.steps.bh_common import TEPBHModel, compute_disformal_metric, compute_curvature_invariants


def check_phi0_1():
    """Check phi_0=1 curvature at specific radii."""
    print("=" * 70)
    print("phi_0=1: Numerical vs Analytical Kretschmann")
    print("Analytical (deep interior): K = 9/(M²r²) + 1/(4M⁴)")
    print("=" * 70)

    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=1.0, delta=0.05, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(1e-8), np.log10(1.5), 100000)
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']

    print(f"\n  {'r':>10} | {'K_numerical':>15} | {'K_analytical':>15} | {'ratio':>10}")
    print(f"  {'-'*10}-+-{'-'*15}-+-{'-'*15}-+-{'-'*10}")

    for r_test in [1.0, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
        idx = np.argmin(np.abs(r - r_test))
        K_num = K[idx]
        K_anal = 9.0 / (r_test ** 2) + 0.25  # deep interior approximation
        ratio = K_num / K_anal if K_anal > 0 else 0
        print(f"  {r_test:>10.1e} | {K_num:>15.6e} | {K_anal:>15.6e} | {ratio:>10.4f}")

    # Power law using different ranges
    print("\n  Power law fits (different ranges):")
    for r_min, r_max in [(1e-8, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 0.1), (1e-4, 1e-1)]:
        mask = (r >= r_min) & (r <= r_max) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
            print(f"    r in [{r_min:.0e}, {r_max:.0e}]: K ~ r^{{{alpha:.3f}}} (expected: r^{{-2}})")


def check_phi0_2():
    """Check phi_0=2 curvature at specific radii."""
    print("\n" + "=" * 70)
    print("phi_0=2: Numerical vs Analytical Kretschmann")
    print("Analytical (deep interior): K ~ 39*r²/(16*M⁶)")
    print("=" * 70)

    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=2.0, delta=0.05, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(1e-8), np.log10(1.5), 100000)
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']

    print(f"\n  {'r':>10} | {'K_numerical':>15} | {'K_analytical':>15} | {'ratio':>10}")
    print(f"  {'-'*10}-+-{'-'*15}-+-{'-'*15}-+-{'-'*10}")

    for r_test in [1.0, 0.5, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
        idx = np.argmin(np.abs(r - r_test))
        K_num = K[idx]
        K_anal = 39.0 / 16.0 * r_test ** 2  # deep interior approximation
        ratio = K_num / K_anal if K_anal > 0 else 0
        print(f"  {r_test:>10.1e} | {K_num:>15.6e} | {K_anal:>15.6e} | {ratio:>10.4f}")

    # Power law using different ranges
    print("\n  Power law fits (different ranges):")
    for r_min, r_max in [(1e-8, 1e-6), (1e-6, 1e-4), (1e-4, 1e-2), (1e-2, 0.1), (1e-4, 1e-1)]:
        mask = (r >= r_min) & (r <= r_max) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
            print(f"    r in [{r_min:.0e}, {r_max:.0e}]: K ~ r^{{{alpha:.3f}}} (expected: r^{{2}})")


def check_schwarzschild():
    """Check Schwarzschild curvature at specific radii."""
    print("\n" + "=" * 70)
    print("Schwarzschild: Numerical vs Analytical Kretschmann")
    print("Analytical: K = 48*M²/r⁶")
    print("=" * 70)

    # Use very small phi_0 to approximate Schwarzschild
    model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=0.001, delta=0.001, M=1.0, sigma_B=1.5)
    r = np.logspace(np.log10(2.5), np.log10(100.0), 50000)
    metric = compute_disformal_metric(r, model)
    curvature = compute_curvature_invariants(r, metric, M=1.0)

    K = curvature['Kretschmann']

    print(f"\n  {'r':>10} | {'K_numerical':>15} | {'K_analytical':>15} | {'ratio':>10}")
    print(f"  {'-'*10}-+-{'-'*15}-+-{'-'*15}-+-{'-'*10}")

    for r_test in [3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0]:
        idx = np.argmin(np.abs(r - r_test))
        K_num = K[idx]
        K_anal = 48.0 / r_test ** 6
        ratio = K_num / K_anal if K_anal > 0 else 0
        print(f"  {r_test:>10.1f} | {K_num:>15.6e} | {K_anal:>15.6e} | {ratio:>10.4f}")


def check_power_law_all_phi0():
    """Check power law for all phi_0 using appropriate range."""
    print("\n" + "=" * 70)
    print("Power law K ~ r^{4*phi_0-6} for various phi_0")
    print("Using range r in [1e-4, 1e-2] (deep interior, away from boundary effects)")
    print("=" * 70)

    print(f"\n  {'phi_0':>6} | {'alpha_numerical':>16} | {'alpha_expected':>16} | {'diff':>8} | {'areal':>12}")
    print(f"  {'-'*6}-+-{'-'*16}-+-{'-'*16}-+-{'-'*8}-+-{'-'*12}")

    for phi_0 in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5]:
        model = TEPBHModel(beta_A=-1.0, B0=0.0, phi_0=phi_0, delta=0.05, M=1.0, sigma_B=1.5)
        r = np.logspace(np.log10(1e-8), np.log10(1.5), 100000)
        metric = compute_disformal_metric(r, model)
        curvature = compute_curvature_invariants(r, metric, M=1.0)

        K = curvature['Kretschmann']
        areal = metric['A'] * r

        # Use range [1e-4, 1e-2] for power law fit
        mask = (r >= 1e-4) & (r <= 1e-2) & np.isfinite(K) & (K > 0)
        if np.sum(mask) > 10:
            alpha = np.polyfit(np.log(r[mask]), np.log(K[mask]), 1)[0]
        else:
            alpha = np.nan

        expected = 4 * phi_0 - 6
        diff = abs(alpha - expected)
        areal_str = "diverges" if areal[0] > 1e3 else "finite"

        print(f"  {phi_0:>6.2f} | {alpha:>16.3f} | {expected:>16.1f} | {diff:>8.3f} | {areal_str:>12}")


if __name__ == '__main__':
    check_schwarzschild()
    check_phi0_1()
    check_phi0_2()
    check_power_law_all_phi0()
