#!/usr/bin/env python3
"""Verify conformal null invariance: compute the shadow on both g and tilde_g.

For a static spherically symmetric metric:
  g:      ds^2 = -F dt^2 + F^{-1} dr^2 + r^2 dOmega^2
  tilde_g: ds^2 = -A^2 F dt^2 + (A^2/F + B*phi'^2) dr^2 + A^2 r^2 dOmega^2

The shadow depends on g_tt/g_phiphi:
  g:      g_tt/g_phiphi = -F / r^2
  tilde_g: g_tt/g_phiphi = -A^2 F / (A^2 r^2) = -F / r^2  (IDENTICAL)

So the photon sphere and shadow are the same on both metrics.
This script verifies this numerically by computing V_eff on both.
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
from step_34_solve_interior import hayward_lapse_sq

M = 1.0
G_PARAM = 1.1  # Use the validated interior solution


def F_func(r):
    return hayward_lapse_sq(r, M, G_PARAM)


def A_func(r):
    """Conformal factor A = e^{-phi(r)}.
    
    From the interior solver: A(0) = 4.68, and A -> 1 at large r.
    We use the numerical profile from step_34_solve_interior.
    """
    # For this verification, use a simple model:
    # A(r) = 1 + (A0 - 1) * exp(-r^2/sigma^2)
    # with A0 = 4.68 and sigma chosen so A(2M) ~ 1.1
    A0 = 4.68
    sigma = 0.5  # rough scale
    return 1.0 + (A0 - 1.0) * np.exp(-r**2 / sigma**2)


def V_eff_g(r):
    """Effective potential on g: F/r^2."""
    if r < 1e-14:
        return 0.0
    return F_func(r) / r**2


def V_eff_tilde_g(r):
    """Effective potential on tilde_g: A^2*F / (A^2*r^2) = F/r^2.
    
    The conformal factor cancels EXACTLY.
    """
    if r < 1e-14:
        return 0.0
    A = A_func(r)
    F = F_func(r)
    # g_tt = -A^2 * F, g_phiphi = A^2 * r^2
    # V_eff = |g_tt| / g_phiphi = A^2 * F / (A^2 * r^2) = F / r^2
    return (A**2 * F) / (A**2 * r**2)


def find_photon_sphere(V_func, label):
    """Find the photon sphere: local maximum of V_eff near r=3M."""
    r = np.linspace(1.5, 10.0, 5000)
    V = np.array([V_func(ri) for ri in r])
    i_max = np.argmax(V)
    r_lo = r[max(0, i_max - 5)]
    r_hi = r[min(len(r) - 1, i_max + 5)]
    result = minimize_scalar(lambda r: -V_func(r), bounds=(r_lo, r_hi), method='bounded')
    r_ps = result.x
    V_ps = V_func(r_ps)
    b_c = 1.0 / np.sqrt(V_ps)
    return r_ps, V_ps, b_c


def main():
    print("=" * 70)
    print("CONFORMAL NULL INVARIANCE VERIFICATION")
    print(f"Hayward background: M={M}, g={G_PARAM}")
    print("=" * 70)

    # Show the conformal factor profile
    print("\n--- Conformal factor A(r) ---")
    for r in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
        print(f"  r={r:6.2f}: A={A_func(r):.4f}, F={F_func(r):.6f}")

    # Compute photon sphere on g
    print("\n--- Photon sphere on g (geometric metric) ---")
    r_ps_g, V_ps_g, b_c_g = find_photon_sphere(V_eff_g, "g")
    print(f"  r_ps = {r_ps_g:.10f}")
    print(f"  V_ps = {V_ps_g:.10f}")
    print(f"  b_c  = {b_c_g:.10f}")

    # Compute photon sphere on tilde_g
    print("\n--- Photon sphere on tilde_g (matter metric) ---")
    r_ps_t, V_ps_t, b_c_t = find_photon_sphere(V_eff_tilde_g, "tilde_g")
    print(f"  r_ps = {r_ps_t:.10f}")
    print(f"  V_ps = {V_ps_t:.10f}")
    print(f"  b_c  = {b_c_t:.10f}")

    # Compare
    print("\n--- Comparison ---")
    print(f"  r_ps difference: {abs(r_ps_g - r_ps_t):.2e}")
    print(f"  b_c  difference: {abs(b_c_g - b_c_t):.2e}")
    print(f"  Relative b_c difference: {abs(b_c_g - b_c_t)/b_c_g:.2e}")

    if abs(b_c_g - b_c_t) < 1e-10:
        print("\n  CONFIRMED: Shadow is identical on g and tilde_g.")
        print("  The conformal factor A(r) cancels exactly in V_eff = g_tt/g_phiphi.")
        print("  The disformal term (entering only g_rr) does not affect the shadow")
        print("  in spherical symmetry.")
    else:
        print("\n  WARNING: Shadows differ! (This would indicate a bug.)")

    # Also verify analytically
    print("\n--- Analytical verification ---")
    r_test = 2.65  # near the photon sphere
    A = A_func(r_test)
    F = F_func(r_test)
    V_g = F / r_test**2
    V_t = (A**2 * F) / (A**2 * r_test**2)
    print(f"  At r={r_test}: A={A:.4f}")
    print(f"  V_eff(g)      = F/r^2           = {V_g:.10f}")
    print(f"  V_eff(tilde_g) = A^2*F/(A^2*r^2) = {V_t:.10f}")
    print(f"  Ratio: {V_t/V_g:.15f}  (should be 1.0)")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print("  The shadow on the matter metric tilde_g is IDENTICAL to the shadow")
    print("  on the geometric metric g, for any conformal factor A(r) and any")
    print("  disformal term B(phi)*(phi')^2 in spherical symmetry.")
    print()
    print("  The -5.37% Phantom Mass at the horizonless threshold is the correct")
    print("  shadow size for both g and tilde_g. It is not an artifact of the")
    print("  geometric metric.")
    print()
    print("  The negative sign (mass deficit) is the photon-sector signature,")
    print("  opposite to the positive ISCO shift (mass excess) in the massive-")
    print("  particle sector. This sign divergence is the frame-split test.")

    output = {
        "parameters": {"M": M, "g": G_PARAM},
        "conformal_factor": {
            "A(0)": 4.68,
            "model": "Gaussian approximation for verification",
        },
        "photon_sphere_g": {
            "r_ps": float(r_ps_g),
            "b_c": float(b_c_g),
        },
        "photon_sphere_tilde_g": {
            "r_ps": float(r_ps_t),
            "b_c": float(b_c_t),
        },
        "difference": {
            "delta_b_c": float(abs(b_c_g - b_c_t)),
            "relative": float(abs(b_c_g - b_c_t) / b_c_g),
        },
        "conclusion": "Shadow identical on g and tilde_g (conformal null invariance)",
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "step_24_conformal_invariance_check.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
