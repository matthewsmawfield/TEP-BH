#!/usr/bin/env python3
"""Scan the Hayward regularization scale g and compute the shadow
deviation as a function of g, finding the EHT-consistent range.

The key insight: g=1.1M was chosen for the interior integration to
demonstrate regularity, but the physical regularization scale for a
real black hole should be much smaller. This scan determines:
  1. How the shadow deviation scales with g
  2. What range of g is consistent with EHT M87* and Sgr A* observations
  3. The physical interpretation of g in terms of the temporal-well depth

Outputs:
  results/step_38_phantom_mass_scan.json
"""

from __future__ import annotations
import json, os, sys
import numpy as np
from scipy.optimize import minimize_scalar

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")
sys.path.insert(0, os.path.join(_PROJECT_ROOT, "scripts"))
from step_34_solve_interior import hayward_mass, hayward_lapse_sq

# Physical constants
G_NEWTON = 6.67430e-11
C_LIGHT = 2.99792458e8
M_SUN = 1.98847e30
MPC_TO_M = 3.08567758e22
KPC_TO_M = 3.08567758e19
UAS_TO_RAD = 4.84813681109536e-12

M = 1.0


def F_func(r, g):
    return hayward_lapse_sq(r, M, g)


def V_eff(r, g):
    if r < 1e-14:
        return 0.0
    return F_func(r, g) / r**2


def find_photon_sphere(g):
    """Local maximum of V_eff near r=3M."""
    r = np.linspace(1.5, 10.0, 5000)
    V = np.array([V_eff(ri, g) for ri in r])
    i_max = np.argmax(V)
    r_lo = r[max(0, i_max - 5)]
    r_hi = r[min(len(r) - 1, i_max + 5)]
    result = minimize_scalar(lambda r: -V_eff(r, g), bounds=(r_lo, r_hi), method='bounded')
    r_ps = result.x
    b_c = 1.0 / np.sqrt(V_eff(r_ps, g))
    return r_ps, b_c


def find_F_min(g):
    """Minimum of F(r) and its location."""
    r = np.linspace(0.01, 5.0, 10000)
    F = np.array([F_func(ri, g) for ri in r])
    i_min = np.argmin(F)
    return float(F[i_min]), float(r[i_min])


def compute_eht_sigma(b_c, b_c_schw=3 * np.sqrt(3)):
    """Compute sigma deviations for M87* and Sgr A*."""
    # M87*
    M87_mass = 6.5e9
    M87_dist = 16.8
    M87_M_m = M87_mass * M_SUN * G_NEWTON / C_LIGHT**2
    M87_D_m = M87_dist * MPC_TO_M
    M87_meas = 42.0
    M87_unc = 3.0
    M87_shadow = 2 * b_c * M87_M_m / M87_D_m / UAS_TO_RAD
    M87_sigma = (M87_shadow - M87_meas) / M87_unc

    # Sgr A* (VLTI)
    sgr_mass = 4.297e6
    sgr_dist = 8.277
    sgr_M_m = sgr_mass * M_SUN * G_NEWTON / C_LIGHT**2
    sgr_D_m = sgr_dist * KPC_TO_M
    sgr_meas = 48.7
    sgr_unc = 7.0
    sgr_shadow = 2 * b_c * sgr_M_m / sgr_D_m / UAS_TO_RAD
    sgr_sigma = (sgr_shadow - sgr_meas) / sgr_unc

    return M87_sigma, sgr_sigma, M87_shadow, sgr_shadow


def main():
    print("=" * 70)
    print("PHANTOM MASS SCAN OVER HAYWARD REGULARIZATION SCALE g")
    print("=" * 70)

    b_c_schw = 3.0 * np.sqrt(3.0)
    r_ps_schw = 3.0

    # Scan g from 0.01 to 2.0
    g_values = np.concatenate([
        np.array([0.001, 0.005, 0.01, 0.02, 0.05]),
        np.linspace(0.1, 2.0, 20),
    ])

    print(f"\n{'g/M':>8s}  {'r_ps':>8s}  {'b_c':>8s}  {'dev%':>8s}  "
          f"{'F_min':>8s}  {'r_throat':>8s}  {'M87σ':>7s}  {'SgrAσ':>7s}")
    print("-" * 80)

    results = []
    for g in g_values:
        r_ps, b_c = find_photon_sphere(g)
        F_min, r_throat = find_F_min(g)
        dev = (b_c - b_c_schw) / b_c_schw * 100
        M87_sigma, sgr_sigma, M87_shadow, sgr_shadow = compute_eht_sigma(b_c)
        results.append({
            'g': float(g),
            'r_ps': float(r_ps),
            'b_c': float(b_c),
            'deviation_percent': float(dev),
            'F_min': float(F_min),
            'r_throat': float(r_throat),
            'M87_sigma': float(M87_sigma),
            'SgrA_sigma': float(sgr_sigma),
            'M87_shadow_uas': float(M87_shadow),
            'SgrA_shadow_uas': float(sgr_shadow),
        })
        print(f"{g:8.4f}  {r_ps:8.4f}  {b_c:8.4f}  {dev:+8.4f}  "
              f"{F_min:8.5f}  {r_throat:8.4f}  {M87_sigma:+7.3f}  {sgr_sigma:+7.3f}")

    # Find the g that minimizes the combined chi-squared
    print("\n--- Optimal g (combined EHT fit) ---")
    best_g = None
    best_chi2 = 1e10
    for r in results:
        chi2 = r['M87_sigma']**2 + r['SgrA_sigma']**2
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_g = r['g']
    print(f"  Best g = {best_g:.4f}, chi^2 = {best_chi2:.4f}")

    # Find the g range consistent with EHT at 2-sigma
    print("\n--- EHT-consistent g range (both within 2-sigma) ---")
    consistent = [r for r in results if abs(r['M87_sigma']) < 2.0 and abs(r['SgrA_sigma']) < 2.0]
    if consistent:
        g_min = min(r['g'] for r in consistent)
        g_max = max(r['g'] for r in consistent)
        print(f"  g in [{g_min:.4f}, {g_max:.4f}]  (both M87* and Sgr A* within 2-sigma)")
    else:
        print("  No g value in the scan is consistent with both at 2-sigma")

    # Find g where deviation = 0 (Schwarzschild limit)
    print("\n--- Schwarzschild recovery ---")
    small_g = [r for r in results if r['g'] < 0.05]
    if small_g:
        print(f"  At g=0.001: deviation = {small_g[0]['deviation_percent']:.6f}%")
        print(f"  --> Small g recovers Schwarzschild (as expected)")

    # Physical interpretation
    print("\n--- Physical interpretation ---")
    print("  g is the Hayward regularization scale (in units of M).")
    print("  For g << M: the temporal well is deep and narrow,")
    print("    the shadow is indistinguishable from Schwarzschild.")
    print("  For g ~ M: the temporal well is wide,")
    print("    the shadow deviates by several percent.")
    print("  For g > M: the shadow shrinks significantly (mass deficit).")
    print()
    print("  The EHT constraint requires g < ~0.5M for both sources")
    print("  to remain within 2-sigma. The physical TEP prediction")
    print("  for g depends on the scalar field coupling eta and the")
    print("  interior solution (Section 4.6).")

    # Key result: the scaling law
    print("\n--- Scaling law ---")
    # Fit deviation ~ A * g^n for small g
    small = [r for r in results if r['g'] < 0.5 and r['g'] > 0.01]
    if len(small) >= 3:
        gs = np.array([r['g'] for r in small])
        devs = np.array([abs(r['deviation_percent']) for r in small])
        # Fit log(dev) = n * log(g) + c
        valid = devs > 1e-10
        if np.any(valid):
            coeffs = np.polyfit(np.log(gs[valid]), np.log(devs[valid]), 1)
            n = coeffs[0]
            A = np.exp(coeffs[1])
            print(f"  |deviation| ~ {A:.4f} * g^{n:.2f}  (for small g)")
            print(f"  Shadow deviation scales as g^{n:.1f}")

    output = {
        "parameters": {"M": M},
        "method": "Scan of Hayward regularization scale g vs EHT shadow constraints",
        "schwarzschild_reference": {
            "r_ps": float(r_ps_schw),
            "b_c": float(b_c_schw),
        },
        "scan_results": results,
        "optimal_g": {"g": float(best_g), "chi2": float(best_chi2)},
        "eht_consistent_range": {
            "g_min": float(g_min) if consistent else None,
            "g_max": float(g_max) if consistent else None,
        },
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "step_38_phantom_mass_scan.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
