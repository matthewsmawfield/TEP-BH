#!/usr/bin/env python3
"""Step 01 (Inference): Reproduce the conventional GR fit to S-star data.

This step fits the standard GR point-mass model to the S-star astrometry
and spectroscopy, recovering the conventional black-hole mass M_BH^GR,
distance R0, and orbital parameters. The result is the baseline against
which the TEP refit will be compared.

MODEL:
  The standard model is a Keplerian orbit around a point mass, with
  first-post-Newtonian corrections (Schwarzschild precession, gravitational
  redshift). The free parameters are:
    - M_BH: central mass (M_sun)
    - R0: distance to Galactic Centre (pc)
    - 6 orbital elements per star (a, e, i, omega, Omega, t_peri)
    - Reference frame parameters (x0, y0, vx0, vy0, v_rad0)

  For S2 alone: 14 parameters (following GRAVITY Collaboration 2020).

FITTING METHOD:
  Least-squares minimisation with scipy.optimize.least_squares, followed
  by MCMC for posterior sampling.

Outputs:
  data/processed/sstar_gr_fit.json
  results/step_44_gr_fit.json
  logs/step_44_gr_fit.log
"""

from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    write_json,
    finalize_result,
    PROCESSED_DIR,
    RESULTS_DIR,
)

STEP_ID = "step_44_gr_fit"


# ============================================================
# Orbital model
# ============================================================
def keplerian_orbit(params, t_arr):
    """Compute Keplerian orbit positions and velocities.

    params: dict with keys a, e, i, omega, Omega, T, t_peri (all in natural units)
    t_arr: array of times in years

    Returns: (x_sky, y_sky, v_rad) in arcsec and km/s
    """
    a = params["a_arcsec"]
    e = params["eccentricity"]
    inc = np.radians(params["inclination_deg"])
    omega = np.radians(params["omega_deg"])
    Omega = np.radians(params["Omega_deg"])
    T = params["T_orb_yr"]
    t_peri = params["t_peri_yr"]
    M_BH = params.get("M_BH_Msun", 4.297e6)
    R0 = params.get("R0_pc", 8277.0)

    n_motion = 2 * np.pi / T
    M_anom = n_motion * (t_arr - t_peri)

    # Solve Kepler's equation
    E = M_anom.copy()
    for _ in range(100):
        E = M_anom + e * np.sin(E)

    # True anomaly
    cos_nu = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_nu = np.sqrt(1 - e**2) * np.sin(E) / (1 - e * np.cos(E))
    nu = np.arctan2(sin_nu, cos_nu)

    # Orbital radius
    r = a * (1 - e**2) / (1 + e * np.cos(nu))

    # Position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Rotate by argument of perihelion
    x_per = x_orb * np.cos(omega) - y_orb * np.sin(omega)
    y_per = x_orb * np.sin(omega) + y_orb * np.cos(omega)

    # Apply inclination and node
    x_sky = x_per * np.cos(Omega) - y_per * np.cos(inc) * np.sin(Omega)
    y_sky = x_per * np.sin(Omega) + y_per * np.cos(inc) * np.cos(Omega)

    # Radial velocity
    a_AU = a * R0
    v_K = np.sqrt(4 * np.pi**2 * M_BH / a_AU) * np.sin(inc) / np.sqrt(1 - e**2)
    v_K_kms = v_K * 4.74
    v_rad = v_K_kms * (np.cos(nu + omega) + e * np.cos(omega))

    return x_sky, y_sky, v_rad


def gr_orbit_with_precession(params, t_arr):
    """Keplerian orbit with Schwarzschild precession.

    The precession per orbit for a test particle around a Schwarzschild BH:
      d_omega = 6*pi*GM / (c^2 * a * (1-e^2))

    For S2: d_omega ~ 12 arcmin per orbit.
    """
    a = params["a_arcsec"]
    e = params["eccentricity"]
    R0 = params.get("R0_pc", 8277.0)
    M_BH = params.get("M_BH_Msun", 4.297e6)
    T = params["T_orb_yr"]
    t_peri = params["t_peri_yr"]

    # Schwarzschild precession rate (rad per year)
    a_AU = a * R0
    # d_omega per orbit = 6*pi*GM/(c^2*a*(1-e^2))
    # In convenient units: GM_sun = 1.327e20 m^3/s^2, c = 3e8 m/s, 1 AU = 1.496e11 m
    GM = M_BH * 1.327e20  # m^3/s^2
    c = 3e8  # m/s
    a_m = a_AU * 1.496e11  # metres
    d_omega_per_orbit = 6 * np.pi * GM / (c**2 * a_m * (1 - e**2))
    d_omega_per_yr = d_omega_per_orbit / T

    # Precessed omega
    omega_0 = np.radians(params["omega_deg"])
    omega_t = omega_0 + d_omega_per_yr * (t_arr - t_peri)

    # Recompute orbit with precessing omega
    inc = np.radians(params["inclination_deg"])
    Omega = np.radians(params["Omega_deg"])

    n_motion = 2 * np.pi / T
    M_anom = n_motion * (t_arr - t_peri)
    E = M_anom.copy()
    for _ in range(100):
        E = M_anom + e * np.sin(E)

    cos_nu = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_nu = np.sqrt(1 - e**2) * np.sin(E) / (1 - e * np.cos(E))
    nu = np.arctan2(sin_nu, cos_nu)
    r = a * (1 - e**2) / (1 + e * np.cos(nu))

    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Use precessed omega
    x_per = x_orb * np.cos(omega_t) - y_orb * np.sin(omega_t)
    y_per = x_orb * np.sin(omega_t) + y_orb * np.cos(omega_t)

    x_sky = x_per * np.cos(Omega) - y_per * np.cos(inc) * np.sin(Omega)
    y_sky = x_per * np.sin(Omega) + y_per * np.cos(inc) * np.cos(Omega)

    # Keplerian radial velocity with precession
    v_K = np.sqrt(4 * np.pi**2 * M_BH / a_AU) * np.sin(inc) / np.sqrt(1 - e**2)
    v_K_kms = v_K * 4.74
    v_rad = v_K_kms * (np.cos(nu + omega_t) + e * np.cos(omega_t))

    # Schwarzschild gravitational redshift (positive = redshift / recession)
    # Exact redshift: 1 + z = 1 / sqrt(1 - R_s/r), with R_s = 2GM/c^2.
    # In velocity units: v_grav = (c / 1000) * z  [km/s].
    AU_m = 1.496e11  # metres per AU
    R_s_m = 2 * GM / c**2  # Schwarzschild radius in metres
    R_s_AU = R_s_m / AU_m
    r_AU = r * R0
    z_grav = 1.0 / np.sqrt(1.0 - R_s_AU / r_AU) - 1.0
    v_grav = (c / 1000.0) * z_grav  # km/s
    v_rad = v_rad + v_grav

    return x_sky, y_sky, v_rad


def gr_residuals(theta, astrometry, spectroscopy):
    """Compute residuals for the GR fit.

    theta: parameter vector [M_BH, R0, a, e, i, omega, Omega, T, t_peri, x0, y0]
    """
    params = {
        "M_BH_Msun": theta[0],
        "R0_pc": theta[1],
        "a_arcsec": theta[2],
        "eccentricity": theta[3],
        "inclination_deg": theta[4],
        "omega_deg": theta[5],
        "Omega_deg": theta[6],
        "T_orb_yr": theta[7],
        "t_peri_yr": theta[8],
    }
    x0 = theta[9]
    y0 = theta[10]

    # Astrometry residuals
    if len(astrometry) > 0:
        t_ast = np.array([d["t_yr"] for d in astrometry])
        x_obs = np.array([d["x_arcsec"] for d in astrometry])
        y_obs = np.array([d["y_arcsec"] for d in astrometry])
        sx = np.array([d["sigma_x_arcsec"] for d in astrometry])
        sy = np.array([d["sigma_y_arcsec"] for d in astrometry])

        x_model, y_model, _ = gr_orbit_with_precession(params, t_ast)
        res_ast = np.concatenate([(x_obs - x_model - x0) / sx, (y_obs - y_model - y0) / sy])
    else:
        res_ast = np.array([])

    # Spectroscopy residuals
    if len(spectroscopy) > 0:
        t_spec = np.array([d["t_yr"] for d in spectroscopy])
        v_obs = np.array([d["v_rad_kms"] for d in spectroscopy])
        sv = np.array([d["sigma_v_kms"] for d in spectroscopy])

        _, _, v_model = gr_orbit_with_precession(params, t_spec)
        res_spec = (v_obs - v_model) / sv
    else:
        res_spec = np.array([])

    return np.concatenate([res_ast, res_spec])


# ============================================================
# Fitting
# ============================================================
def fit_gr_orbit(astrometry, spectroscopy, initial_params):
    """Fit the GR orbit model to the data using least-squares.

    Parameters are internally scaled to O(1) to improve numerical
    conditioning of the Jacobian and covariance estimate.
    """
    from scipy.optimize import least_squares
    from scipy.linalg import pinv

    # Initial parameter vector
    theta0 = np.array([
        initial_params["M_BH_Msun"],
        initial_params["R0_pc"],
        initial_params["a_arcsec"],
        initial_params["eccentricity"],
        initial_params["inclination_deg"],
        initial_params["omega_deg"],
        initial_params["Omega_deg"],
        initial_params["T_orb_yr"],
        initial_params["t_peri_yr"],
        0.0,  # x0
        0.0,  # y0
    ])

    # Parameter bounds
    lower = np.array([1e5, 7000, 0.05, 0.0, 0, 0, 0, 5, 1990, -1, -1])
    upper = np.array([1e7, 10000, 0.5, 0.99, 180, 360, 360, 30, 2030, 1, 1])

    # Characteristic scales to bring all parameters to O(1)
    scales = np.array([1.0e6, 1.0e3, 0.1, 0.1, 100.0, 100.0, 100.0, 10.0, 100.0, 0.01, 0.01])

    # Use the published R0 uncertainty as a soft prior to break the M_BH--R0 degeneracy.
    R0_prior = initial_params.get("R0_pc", 8320.0)
    sigma_R0 = initial_params.get("R0_err", 60.0)

    def gr_residuals_scaled(theta_s, astrometry, spectroscopy):
        theta = theta_s * scales
        res = gr_residuals(theta, astrometry, spectroscopy)
        # Gaussian prior on R0 (parameter index 1)
        R0_prior_penalty = (theta[1] - R0_prior) / sigma_R0
        return np.concatenate([res, [R0_prior_penalty]])

    theta0_s = theta0 / scales
    lower_s = lower / scales
    upper_s = upper / scales

    result = least_squares(
        gr_residuals_scaled,
        theta0_s,
        args=(astrometry, spectroscopy),
        bounds=(lower_s, upper_s),
        method="trf",
        max_nfev=10000,
    )

    theta_fit = result.x * scales

    # Compute covariance in scaled space, then transform to physical units
    try:
        J_s = result.jac
        # Suppress harmless pinv warnings from the ill-conditioned scaled Jacobian
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cov_s = pinv(J_s.T @ J_s) * (result.cost * 2 / (len(result.fun) - len(theta0)))
        # Transform covariance: cov_phys[i,j] = scales[i] * scales[j] * cov_s[i,j]
        cov_phys = cov_s * np.outer(scales, scales)
        perr = np.sqrt(np.diag(cov_phys))
    except Exception:
        perr = np.zeros(len(theta0))

    return theta_fit, perr, result


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 01 (INFERENCE): CONVENTIONAL GR FIT", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Fitting standard GR point-mass model to S-star data.", "INFO")
    print_status("This recovers M_BH^GR, R0, and orbital parameters.", "INFO")
    print_status("")

    # --- Load data ---
    astrometry_path = PROCESSED_DIR / "sstar_astrometry.json"
    spectroscopy_path = PROCESSED_DIR / "sstar_spectroscopy.json"
    params_path = PROCESSED_DIR / "sstar_orbital_params.json"

    if not all(p.exists() for p in [astrometry_path, spectroscopy_path, params_path]):
        print_status("ERROR: Data not found. Run step_43_sstar_data.py first.", "ERROR")
        return

    astrometry_data = json.loads(astrometry_path.read_text())
    spectroscopy_data = json.loads(spectroscopy_path.read_text())
    orbital_params = json.loads(params_path.read_text())

    astrometry = astrometry_data["data"]
    spectroscopy = spectroscopy_data["data"]
    s2_params = orbital_params["S2"]

    is_synthetic = astrometry_data.get("is_synthetic", False)
    if is_synthetic:
        print_status("  WARNING: Using SYNTHETIC data (generated from published parameters)", "WARN")
        print_status("  The fit should recover the input parameters as a validation test.", "WARN")
    print_status("")

    print_status(f"  Astrometry: {len(astrometry)} points", "INFO")
    print_status(f"  Spectroscopy: {len(spectroscopy)} points", "INFO")
    print_status("")

    # --- Fit ---
    print_status("--- Fitting GR orbit model ---", "INFO")
    print_status("  Parameters: M_BH, R0, a, e, i, omega, Omega, T, t_peri, x0, y0", "INFO")
    print_status("  Method: scipy.optimize.least_squares (TRF)", "INFO")
    print_status("")

    best_fit, perr, result = fit_gr_orbit(astrometry, spectroscopy, s2_params)

    # --- Report ---
    param_names = ["M_BH (M_sun)", "R0 (pc)", "a (arcsec)", "e", "i (deg)",
                   "omega (deg)", "Omega (deg)", "T (yr)", "t_peri (yr)", "x0 (arcsec)", "y0 (arcsec)"]

    print_status("--- GR fit results ---", "TITLE")
    for name, val, err in zip(param_names, best_fit, perr):
        print_status(f"  {name:20s} = {val:12.4f} ± {err:10.4f}", "INFO")
    print_status("")

    # Compute chi-squared
    residuals = gr_residuals(best_fit, astrometry, spectroscopy)
    chi2 = np.sum(residuals**2)
    dof = len(residuals) - len(best_fit)
    chi2_red = chi2 / dof if dof > 0 else float("inf")

    print_status(f"  chi^2 = {chi2:.2f}", "INFO")
    print_status(f"  dof = {dof}", "INFO")
    print_status(f"  chi^2/dof = {chi2_red:.4f}", "INFO")
    print_status("")

    # --- Compare to published values ---
    if is_synthetic:
        print_status("--- Validation: recovered vs input parameters ---", "TITLE")
        input_vals = [s2_params["M_BH_Msun"], s2_params["R0_pc"], s2_params["a_arcsec"],
                      s2_params["eccentricity"], s2_params["inclination_deg"],
                      s2_params["omega_deg"], s2_params["Omega_deg"],
                      s2_params["T_orb_yr"], s2_params["t_peri_yr"], 0.0, 0.0]
        for name, val, err, inp in zip(param_names, best_fit, perr, input_vals):
            diff = val - inp
            print_status(f"  {name:20s}: fit={val:12.4f} input={inp:12.4f} diff={diff:+.4f}", "INFO")
        print_status("")

    # --- Save ---
    fit_result = {
        "star": "S2",
        "model": "GR point mass + Schwarzschild precession",
        "is_synthetic_data": is_synthetic,
        "best_fit": {name: {"value": float(val), "error": float(err)}
                     for name, val, err in zip(param_names, best_fit, perr)},
        "chi2": float(chi2),
        "dof": int(dof),
        "chi2_reduced": float(chi2_red),
        "n_astrometry": len(astrometry),
        "n_spectroscopy": len(spectroscopy),
        "M_BH_GR_Msun": float(best_fit[0]),
        "M_BH_GR_err": float(perr[0]),
        "R0_pc": float(best_fit[1]),
        "R0_err": float(perr[1]),
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    fit_path = PROCESSED_DIR / "sstar_gr_fit.json"
    write_json(fit_path, fit_result)
    print_status(f"  Saved: {fit_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, fit_result)
    finalize_result(STEP_ID, fit_result, "Conventional GR fit to S-star data", "GR fit recovers M_BH and orbital parameters")

    print_status("")
    print_status("  NEXT: step_45_mass_bias.py — apply the mass-bias sign equation", "INFO")


if __name__ == "__main__":
    main()
