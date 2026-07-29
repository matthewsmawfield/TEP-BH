#!/usr/bin/env python3
"""Track 1: Full-geometry ray-tracer for the solved temporal well.

Integrates null geodesics through the full Hayward + sGB interior
(g=1.1, eta=-0.1) and computes the shadow boundary, comparing against
the Schwarzschild benchmark and EHT M87* / Sgr A* observations.

Key physics:
  - The Hayward background has NO horizon: F_min = 0.038 at r ~ 1.39
  - Photons below the Schwarzschild critical impact parameter (b_c = 3*sqrt(3)*M)
    would be captured by a Schwarzschild BH, but on the temporal well they
    transit through the throat and may emerge
  - This modifies the shadow boundary -> the "Phantom Mass" signature

The shadow boundary is determined by the critical impact parameter:
  b_c = r_ps / sqrt(F(r_ps))
where r_ps is the photon sphere (maximum of V_eff = F/r^2).

For the full geometry, we also integrate geodesics to check which
impact parameters transit through the throat vs. scatter back.

Outputs:
  results/step_37_phantom_mass_raytrace.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")

sys.path.insert(0, os.path.join(_PROJECT_ROOT, "scripts"))
from step_34_solve_interior import hayward_mass, hayward_N, hayward_lapse_sq

# Physical constants
G_NEWTON = 6.67430e-11
C_LIGHT = 2.99792458e8
M_SUN = 1.98847e30
MPC_TO_M = 3.08567758e22
KPC_TO_M = 3.08567758e19
UAS_TO_RAD = 4.84813681109536e-12

# Parameters
M = 1.0
G_PARAM = 1.1
ETA = -0.1


# ---------------------------------------------------------------------------
# Background geometry
# ---------------------------------------------------------------------------
def F_func(r):
    """F(r) = N^2(r) = 1 - 2m(r)/r on Hayward background."""
    return hayward_lapse_sq(r, M, G_PARAM)


def N_func(r):
    return hayward_N(r, M, G_PARAM)


def m_func(r):
    return hayward_mass(r, M, G_PARAM)


# ---------------------------------------------------------------------------
# Effective potential and photon sphere
# ---------------------------------------------------------------------------
def V_eff_null(r):
    """Effective potential for null geodesics: V = F(r) / r^2.

    For a metric ds^2 = -N^2 dt^2 + N^{-2} dr^2 + r^2 dOmega^2,
    the radial equation is (dr/dlambda)^2 = E^2 - L^2 * F / r^2.
    The photon sphere is at the maximum of F/r^2.
    """
    if r < 1e-14:
        return 0.0
    return F_func(r) / r**2


def find_photon_sphere():
    """Find the photon sphere: LOCAL maximum of V_eff = F/r^2 near r=3M.

    Note: V_eff = F/r^2 diverges as r -> 0 (since F -> 1 and r^2 -> 0),
    creating an inner potential barrier. The photon sphere is the LOCAL
    maximum in the exterior region (near r=3M for small g), which is the
    unstable circular orbit that determines the shadow boundary.
    """
    # Search for the local maximum near r=3M (exterior region)
    r = np.linspace(1.5, 10.0, 5000)
    V = np.array([V_eff_null(ri) for ri in r])
    i_max = np.argmax(V)
    r_lo = r[max(0, i_max - 5)]
    r_hi = r[min(len(r) - 1, i_max + 5)]
    from scipy.optimize import minimize_scalar
    result = minimize_scalar(lambda r: -V_eff_null(r), bounds=(r_lo, r_hi), method='bounded')
    r_ps = result.x
    V_ps = V_eff_null(r_ps)
    b_c = 1.0 / np.sqrt(V_ps)  # critical impact parameter
    return r_ps, V_ps, b_c


def find_photon_sphere_schwarzschild():
    """Schwarzschild photon sphere for comparison."""
    r_ps = 3.0 * M
    F_ps = 1.0 - 2.0 * M / r_ps
    V_ps = F_ps / r_ps**2
    b_c = 1.0 / np.sqrt(V_ps)
    return r_ps, V_ps, b_c


# ---------------------------------------------------------------------------
# Null geodesic integration
# ---------------------------------------------------------------------------
def geodesic_ode(r, y, b):
    """Null geodesic ODE in the equatorial plane.

    State: y = [phi, u] where u = 1/r (we use the standard approach)
    Actually, let's use (r, phi) with the orbit equation.

    For null geodesics: (dr/dlambda)^2 = E^2 - L^2 * F / r^2
    With b = L/E: (dr/dlambda)^2 = E^2 * (1 - b^2 * F / r^2)

    dphi/dlambda = L / r^2 = b*E / r^2
    dr/dlambda = +/- E * sqrt(1 - b^2 * F / r^2)

    So: dphi/dr = (b / r^2) / sqrt(1 - b^2 * F / r^2)
    """
    phi, r_curr = y
    if r_curr < 1e-14:
        return [0.0, 0.0]

    F = F_func(r_curr)
    discriminant = 1.0 - b**2 * F / r_curr**2

    if discriminant < 0:
        # Turning point - can't proceed further inward
        return [0.0, 0.0]

    dphi_dr = b / (r_curr**2 * np.sqrt(discriminant))
    return [dphi_dr, 1.0]  # phi as function of r, r increases


def integrate_geodesic(b, r_start=50.0, r_end=0.01, direction='inward'):
    """Integrate a null geodesic with impact parameter b.

    Returns (r_turnaround, phi_total, transited) where:
    - r_turnaround: minimum r reached (if geodesic turns back)
    - phi_total: total deflection angle
    - transited: True if geodesic passed through the throat (r < r_throat)
    """
    r_throat = 1.39  # approximate throat radius

    # Start from far away, incoming
    # dphi/dr = b / (r^2 * sqrt(1 - b^2 * F / r^2))
    # We integrate from r_start inward

    def ode_inward(r, y):
        phi = y[0]
        if r < 1e-14:
            return [0.0]
        F = F_func(r)
        disc = 1.0 - b**2 * F / r**2
        if disc < 0:
            return [0.0]  # stop
        dphi_dr = b / (r**2 * np.sqrt(disc))
        return [dphi_dr]

    # Find the turning point first
    r_grid = np.linspace(r_start, r_end, 5000)
    F_grid = np.array([F_func(r) for r in r_grid])
    disc_grid = 1.0 - b**2 * F_grid / r_grid**2

    if np.all(disc_grid > 0):
        # Geodesic transits through the throat!
        # Integrate all the way to r_end
        sol = solve_ivp(ode_inward, [r_start, r_end], [0.0],
                        method='DOP853', rtol=1e-12, atol=1e-14,
                        max_step=0.1, dense_output=True)
        if sol.success:
            return r_end, sol.y[0, -1], True
        else:
            return r_end, 0.0, True
    else:
        # Find the turning point (where disc = 0)
        i_turn = np.where(disc_grid < 0)[0][0]
        r_turn = r_grid[i_turn]

        # Integrate to just before the turning point
        r_stop = r_turn * 1.001
        sol = solve_ivp(ode_inward, [r_start, r_stop], [0.0],
                        method='DOP853', rtol=1e-12, atol=1e-14,
                        max_step=0.1)
        if sol.success:
            return r_turn, sol.y[0, -1], False
        else:
            return r_turn, 0.0, False


def compute_shadow_boundary():
    """Compute the shadow boundary by scanning impact parameters.

    For b > b_c: geodesic scatters back (turning point outside throat)
    For b < b_c: on Schwarzschild, geodesic falls into horizon
                 on temporal well, geodesic may transit through throat

    The shadow boundary is the critical b where the behavior changes.
    On the temporal well, this is modified by the throat.
    """
    r_ps, V_ps, b_c = find_photon_sphere()
    r_ps_s, V_ps_s, b_c_s = find_photon_sphere_schwarzschild()

    print(f"  Photon sphere (Hayward): r_ps = {r_ps:.6f}, b_c = {b_c:.6f}")
    print(f"  Photon sphere (Schwarzschild): r_ps = {r_ps_s:.6f}, b_c = {b_c_s:.6f}")
    print(f"  Deviation: {(b_c - b_c_s) / b_c_s * 100:.6f}%")

    # Scan impact parameters around b_c
    b_values = np.linspace(0.9 * b_c_s, 1.2 * b_c_s, 50)
    results = []
    for b in b_values:
        r_min, phi, transited = integrate_geodesic(b)
        results.append({
            'b': b,
            'r_min': r_min,
            'phi': phi,
            'transited': transited,
        })
        status = "TRANSITED" if transited else f"turned at r={r_min:.3f}"
        print(f"  b={b:.4f}: {status}")

    return r_ps, b_c, r_ps_s, b_c_s, results


# ---------------------------------------------------------------------------
# Phantom Mass calculation
# ---------------------------------------------------------------------------
def compute_phantom_mass(b_c, b_c_schw):
    """Compute the Phantom Mass signature.

    If the shadow is larger/smaller than Schwarzschild for the same
    asymptotic mass M, the observer infers a different mass:
      M_inferred / M_true = b_c / b_c_schw

    The "Phantom Mass" is the apparent mass excess:
      Delta_M / M = (b_c / b_c_schw) - 1
    """
    mass_ratio = b_c / b_c_schw
    phantom_mass_fraction = mass_ratio - 1.0
    return mass_ratio, phantom_mass_fraction


# ---------------------------------------------------------------------------
# EHT comparison
# ---------------------------------------------------------------------------
def compute_eht_comparison(b_c, b_c_schw):
    """Compare shadow sizes against EHT M87* and Sgr A* observations."""
    # M87* parameters
    M87_mass_msun = 6.5e9
    M87_dist_mpc = 16.8
    M87_M_m = M87_mass_msun * M_SUN * G_NEWTON / C_LIGHT**2
    M87_D_m = M87_dist_mpc * MPC_TO_M
    M87_shadow_meas = 42.0  # uas
    M87_shadow_unc = 3.0  # uas

    # Sgr A* parameters (VLTI distance/mass)
    sgr_mass_msun = 4.297e6
    sgr_dist_kpc = 8.277
    sgr_M_m = sgr_mass_msun * M_SUN * G_NEWTON / C_LIGHT**2
    sgr_D_m = sgr_dist_kpc * KPC_TO_M
    sgr_shadow_meas = 48.7  # uas
    sgr_shadow_unc = 7.0  # uas

    # Shadow angular diameter: theta = 2 * b_c * M_phys / D
    M87_shadow_tep = 2 * b_c * M87_M_m / M87_D_m / UAS_TO_RAD
    M87_shadow_schw = 2 * b_c_schw * M87_M_m / M87_D_m / UAS_TO_RAD

    sgr_shadow_tep = 2 * b_c * sgr_M_m / sgr_D_m / UAS_TO_RAD
    sgr_shadow_schw = 2 * b_c_schw * sgr_M_m / sgr_D_m / UAS_TO_RAD

    # Sigma deviations
    M87_sigma = (M87_shadow_tep - M87_shadow_meas) / M87_shadow_unc
    sgr_sigma = (sgr_shadow_tep - sgr_shadow_meas) / sgr_shadow_unc

    M87_sigma_schw = (M87_shadow_schw - M87_shadow_meas) / M87_shadow_unc
    sgr_sigma_schw = (sgr_shadow_schw - sgr_shadow_meas) / sgr_shadow_unc

    return {
        'M87': {
            'mass_Msun': M87_mass_msun,
            'distance_Mpc': M87_dist_mpc,
            'shadow_tep_uas': float(M87_shadow_tep),
            'shadow_schw_uas': float(M87_shadow_schw),
            'shadow_measured_uas': M87_shadow_meas,
            'shadow_uncertainty_uas': M87_shadow_unc,
            'sigma_tep': float(M87_sigma),
            'sigma_schw': float(M87_sigma_schw),
            'tep_deviation_percent': float((M87_shadow_tep - M87_shadow_schw) / M87_shadow_schw * 100),
        },
        'SgrA': {
            'mass_Msun': sgr_mass_msun,
            'distance_kpc': sgr_dist_kpc,
            'shadow_tep_uas': float(sgr_shadow_tep),
            'shadow_schw_uas': float(sgr_shadow_schw),
            'shadow_measured_uas': sgr_shadow_meas,
            'shadow_uncertainty_uas': sgr_shadow_unc,
            'sigma_tep': float(sgr_sigma),
            'sigma_schw': float(sgr_sigma_schw),
            'tep_deviation_percent': float((sgr_shadow_tep - sgr_shadow_schw) / sgr_shadow_schw * 100),
        },
    }


# ---------------------------------------------------------------------------
# Effective potential profile
# ---------------------------------------------------------------------------
def compute_V_eff_profile():
    """Compute V_eff = F/r^2 profile for the full geometry."""
    r = np.concatenate([
        np.linspace(0.01, 1.0, 500),
        np.linspace(1.0, 5.0, 500),
        np.linspace(5.0, 20.0, 500),
    ])
    V = np.array([V_eff_null(ri) for ri in r])
    F = np.array([F_func(ri) for ri in r])
    return r.tolist(), V.tolist(), F.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("TRACK 1: PHANTOM MASS RAY-TRACING ON THE SOLVED TEMPORAL WELL")
    print(f"Parameters: M={M}, g={G_PARAM}, eta={ETA}")
    print("=" * 70)

    # 1. Geometry check
    print("\n--- 1. Geometry ---")
    r_check = np.linspace(0.01, 10, 1000)
    F_check = np.array([F_func(r) for r in r_check])
    print(f"  F_min = {F_check.min():.6f} at r = {r_check[np.argmin(F_check)]:.4f}")
    print(f"  F(2M=2.0) = {F_func(2.0):.6f}")
    print(f"  F(3M=3.0) = {F_func(3.0):.6f}")
    print(f"  No horizon: F > 0 everywhere")

    # 2. Photon sphere and shadow boundary
    print("\n--- 2. Photon sphere and shadow boundary ---")
    r_ps, b_c, r_ps_s, b_c_s, geodesic_results = compute_shadow_boundary()

    # 3. Phantom Mass
    print("\n--- 3. Phantom Mass ---")
    mass_ratio, phantom_fraction = compute_phantom_mass(b_c, b_c_s)
    print(f"  M_inferred / M_true = {mass_ratio:.8f}")
    print(f"  Phantom Mass fraction = {phantom_fraction * 100:.6f}%")
    if phantom_fraction > 0:
        print(f"  --> Shadow is LARGER than Schwarzschild -> apparent mass EXCESS")
    else:
        print(f"  --> Shadow is SMALLER than Schwarzschild -> apparent mass DEFICIT")

    # 4. EHT comparison
    print("\n--- 4. EHT observational comparison ---")
    eht = compute_eht_comparison(b_c, b_c_s)
    for source in ['M87', 'SgrA']:
        d = eht[source]
        print(f"\n  {source}:")
        print(f"    Shadow (TEP):     {d['shadow_tep_uas']:.2f} uas")
        print(f"    Shadow (Schw):    {d['shadow_schw_uas']:.2f} uas")
        print(f"    Shadow (measured):{d['shadow_measured_uas']:.2f} +/- {d['shadow_uncertainty_uas']:.2f} uas")
        print(f"    sigma (TEP):      {d['sigma_tep']:.3f}")
        print(f"    sigma (Schw):     {d['sigma_schw']:.3f}")
        print(f"    TEP deviation:    {d['tep_deviation_percent']:.6f}%")

    # 5. V_eff profile
    print("\n--- 5. Effective potential profile ---")
    r_prof, V_prof, F_prof = compute_V_eff_profile()
    V_max = max(V_prof)
    r_at_max = r_prof[V_prof.index(V_max)]
    print(f"  V_eff maximum: {V_max:.6f} at r = {r_at_max:.4f}")
    print(f"  Schwarzschild V_eff max: {V_eff_null(3.0):.6f} at r = 3.0")

    # 6. Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Photon sphere: r_ps = {r_ps:.6f} (Schwarzschild: {r_ps_s:.6f})")
    print(f"  Shadow boundary: b_c = {b_c:.6f} (Schwarzschild: {b_c_s:.6f})")
    print(f"  Phantom Mass: {phantom_fraction * 100:+.6f}%")
    print(f"  M87* sigma: {eht['M87']['sigma_tep']:.3f} (Schw: {eht['M87']['sigma_schw']:.3f})")
    print(f"  Sgr A* sigma: {eht['SgrA']['sigma_tep']:.3f} (Schw: {eht['SgrA']['sigma_schw']:.3f})")

    # Count transiting geodesics
    n_transit = sum(1 for g in geodesic_results if g['transited'])
    n_total = len(geodesic_results)
    print(f"  Geodesics transiting through throat: {n_transit}/{n_total}")

    output = {
        "parameters": {"M": M, "g": G_PARAM, "eta": ETA},
        "method": "Full-geometry null geodesic integration on Hayward + sGB background",
        "geometry": {
            "F_min": float(F_check.min()),
            "r_F_min": float(r_check[np.argmin(F_check)]),
            "has_horizon": False,
        },
        "photon_sphere": {
            "r_ps_hayward": float(r_ps),
            "r_ps_schwarzschild": float(r_ps_s),
            "b_c_hayward": float(b_c),
            "b_c_schwarzschild": float(b_c_s),
            "deviation_percent": float((b_c - b_c_s) / b_c_s * 100),
        },
        "phantom_mass": {
            "mass_ratio": float(mass_ratio),
            "phantom_fraction": float(phantom_fraction),
            "phantom_percent": float(phantom_fraction * 100),
            "interpretation": "apparent mass excess" if phantom_fraction > 0 else "apparent mass deficit",
        },
        "eht_comparison": eht,
        "geodesic_scan": [
            {"b": g["b"], "r_min": g["r_min"], "transited": g["transited"]}
            for g in geodesic_results
        ],
        "V_eff_profile": {
            "r": r_prof[:200],  # subsample for JSON
            "V": V_prof[:200],
            "F": F_prof[:200],
        },
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "step_37_phantom_mass_raytrace.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
