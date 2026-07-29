#!/usr/bin/env python3
"""Step 04: Observer redshift and clock-rate profile.

Computes the gravitational redshift for static observers at different radii
in the TEP matter metric, for both Schwarzschild and Hayward backgrounds.

For a static observer at radius r in a static metric ds^2 = -F dt^2 + ...,
the proper time is d tau = sqrt(F) dt. The redshift relative to infinity is:

  z(r) = 1/sqrt(F_tilde(r)) - 1

where F_tilde = A^2 * F is the tt-component of the TEP matter metric.

For the TEP metric with A = (r_h/r)^{phi0}:
  z(r) = (r/r_h)^{phi0} / sqrt(F(r)) - 1

Key results:
  - At the horizon (r = r_h): z -> infinity (standard gravitational redshift)
  - On Schwarzschild inside the horizon: static observers don't exist (F < 0)
  - On Hayward at the centre (r -> 0): F -> 1, A -> infinity, z -> -1
    (the conformal divergence produces a BLUESHIFT, not a redshift, at the
    centre — the infinite redshift is at the horizon, not the centre)

The "temporal freeze" (infinite redshift of outgoing light) is a coordinate
effect in Eddington-Finkelstein time v, not a static-observer redshift.

Outputs:
  results/step_18_observer_redshift.json
  results/step_18_observer_redshift.csv
"""

from __future__ import annotations

import csv
import json
import os
import sys

import sympy as sp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from scripts.utils.logger import TEPLogger  # noqa: E402

r, M, ell, phi0, r_h = sp.symbols("r M ell phi0 r_h", positive=True)
RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")


def observer_redshift(F, A_factor):
    """Gravitational redshift z(r) = 1/sqrt(A^2 * F) - 1 for a static observer.

    Valid only where F > 0 (outside the horizon for Schwarzschild).
    Inside the horizon, static observers don't exist.
    """
    F_tilde = A_factor**2 * F
    z = 1 / sp.sqrt(F_tilde) - 1
    return sp.simplify(z)


def clock_rate_ratio(F, A_factor):
    """Clock rate ratio d(tau_local)/d(tau_infinity) = sqrt(A^2 * F).

    This is the ratio of proper time for a static observer at r to coordinate
    time (which equals proper time at infinity).
    """
    F_tilde = A_factor**2 * F
    return sp.sqrt(F_tilde)


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    logger = TEPLogger("step_18_observer_redshift")
    logger.info("Step 04: Observer Redshift and Clock-Rate Profile")

    # --- Symbolic computation ---
    F_schw = 1 - 2 * M / r
    F_hay = 1 - 2 * M * r**2 / (r**3 + 2 * M * ell**2)

    A1 = (r_h / r) ** 1  # phi0 = 1
    A2 = (r_h / r) ** 2  # phi0 = 2

    # Redshift formulas (symbolic)
    z_schw_phi0_1 = observer_redshift(F_schw, A1)
    z_schw_phi0_2 = observer_redshift(F_schw, A2)
    z_hay_phi0_1 = observer_redshift(F_hay, A1)

    # Clock rate ratios
    clock_schw_phi0_1 = clock_rate_ratio(F_schw, A1)
    clock_hay_phi0_1 = clock_rate_ratio(F_hay, A1)

    print("=== Symbolic Redshift Formulas ===")
    print(f"  z_Schw(phi0=1) = {z_schw_phi0_1}")
    print(f"  z_Schw(phi0=2) = {z_schw_phi0_2}")
    print(f"  z_Hay(phi0=1)  = {z_hay_phi0_1}")
    print()
    print("=== Clock Rate Ratios ===")
    print(f"  dtau/dt_Schw(phi0=1) = {clock_schw_phi0_1}")
    print(f"  dtau/dt_Hay(phi0=1)  = {clock_hay_phi0_1}")

    # --- Limits ---
    print("\n=== Limits ===")
    # Schwarzschild horizon: z -> infinity
    z_schw_horizon = sp.limit(z_schw_phi0_1, r, 2 * M, "+")
    print(f"  z_Schw(phi0=1) at horizon (r->2M+) = {z_schw_horizon}")

    # Hayward centre: F -> 1, A -> infinity
    z_hay_centre = sp.limit(z_hay_phi0_1, r, 0)
    print(f"  z_Hay(phi0=1) at centre (r->0) = {z_hay_centre}")

    # Hayward horizon
    # r_h for Hayward is slightly less than 2M
    # At the horizon F -> 0, so z -> infinity
    z_hay_horizon = sp.limit(z_hay_phi0_1, r, 2 * M, "+")
    print(f"  z_Hay(phi0=1) at r->2M+ = {z_hay_horizon}")

    # Clock rate at Hayward centre
    clock_hay_centre = sp.limit(clock_hay_phi0_1, r, 0)
    print(f"  dtau/dt_Hay(phi0=1) at centre = {clock_hay_centre}")

    # --- Numerical radial scan ---
    M_num = 1.0
    ell_num = 0.5
    r_h_num = 2.0 * M_num  # approximate horizon

    # Radial grid (outside horizon for Schwarzschild, full range for Hayward)
    r_vals_ext = np.linspace(2.01, 50.0, 200)  # exterior
    r_vals_int_hay = np.linspace(0.001, 1.99, 200)  # interior (Hayward only)

    z_schw_func = sp.lambdify((r, M, r_h), z_schw_phi0_1, "numpy")
    z_hay_func = sp.lambdify((r, M, ell, r_h), z_hay_phi0_1, "numpy")
    clock_hay_func = sp.lambdify((r, M, ell, r_h), clock_hay_phi0_1, "numpy")

    # --- Key finding ---
    print("\n=== Key Finding ===")
    print("On the Hayward background with phi0=1:")
    print(f"  z(centre) = {z_hay_centre}  (BLUESHIFT, not infinite redshift!)")
    print(f"  z(horizon) -> infinity  (standard gravitational redshift at horizon)")
    print(f"  dtau/dt(centre) = {clock_hay_centre}  (clock rate -> infinity, not zero!)")
    print()
    print("The 'temporal freeze' (infinite redshift) occurs at the HORIZON,")
    print("not at the centre. At the centre, the conformal divergence A -> infinity")
    print("makes the local clock run FASTER (dtau/dt -> infinity), producing a")
    print("blueshift z -> -1, not a redshift.")
    print()
    print("The manuscript's claim that 'redshift becomes infinite as r -> 0'")
    print("is correct for the Schwarzschild background (horizon at r=2M), but")
    print("for the Hayward regular benchmark the infinite redshift is at the")
    print("horizon, and the centre has a blueshift.")

    # --- Write outputs ---
    output = {
        "step": "04_observer_redshift",
        "description": "Observer redshift and clock-rate profile",
        "formulas": {
            "z_Schwarzschild_phi0_1": str(z_schw_phi0_1),
            "z_Schwarzschild_phi0_2": str(z_schw_phi0_2),
            "z_Hayward_phi0_1": str(z_hay_phi0_1),
            "clock_rate_Schwarzschild_phi0_1": str(clock_schw_phi0_1),
            "clock_rate_Hayward_phi0_1": str(clock_hay_phi0_1),
        },
        "limits": {
            "z_Schwarzschild_horizon": str(z_schw_horizon),
            "z_Hayward_centre": str(z_hay_centre),
            "z_Hayward_horizon_approx": str(z_hay_horizon),
            "clock_rate_Hayward_centre": str(clock_hay_centre),
        },
        "key_finding": (
            "The infinite redshift (temporal freeze) occurs at the HORIZON, not at "
            "the centre. On the Hayward regular benchmark with phi0=1, the centre "
            "has z -> -1 (blueshift) and dtau/dt -> infinity (clock runs faster), "
            "because the conformal divergence A -> infinity dominates over F -> 1. "
            "The manuscript's 'redshift becomes infinite as r -> 0' is correct for "
            "the Schwarzschild horizon but not for the Hayward centre."
        ),
    }

    out_json = os.path.join(RESULTS_DIR, "step_18_observer_redshift.json")
    with open(out_json, "w") as fh:
        json.dump(output, fh, indent=2, default=str)
    print(f"\nWrote {out_json}")

    # CSV: radial scan
    out_csv = os.path.join(RESULTS_DIR, "step_18_observer_redshift.csv")
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["r", "z_Schwarzschild_phi0_1", "z_Hayward_phi0_1", "clock_rate_Hayward_phi0_1"])
        # Exterior
        for rv in r_vals_ext:
            try:
                zs = float(z_schw_func(rv, M_num, r_h_num))
            except (ValueError, ZeroDivisionError, OverflowError):
                zs = float("nan")
            try:
                zh = float(z_hay_func(rv, M_num, ell_num, r_h_num))
            except (ValueError, ZeroDivisionError, OverflowError):
                zh = float("nan")
            try:
                ch = float(clock_hay_func(rv, M_num, ell_num, r_h_num))
            except (ValueError, ZeroDivisionError, OverflowError):
                ch = float("nan")
            writer.writerow([rv, zs, zh, ch])
        # Interior (Hayward only)
        for rv in r_vals_int_hay:
            try:
                zh = float(z_hay_func(rv, M_num, ell_num, r_h_num))
            except (ValueError, ZeroDivisionError, OverflowError):
                zh = float("nan")
            try:
                ch = float(clock_hay_func(rv, M_num, ell_num, r_h_num))
            except (ValueError, ZeroDivisionError, OverflowError):
                ch = float("nan")
            writer.writerow([rv, float("nan"), zh, ch])
    print(f"Wrote {out_csv}")

    logger.info("Step 04 complete")
    return output


if __name__ == "__main__":
    main()
