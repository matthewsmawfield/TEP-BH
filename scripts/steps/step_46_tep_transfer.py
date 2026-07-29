#!/usr/bin/env python3
"""Step 03 (Inference): TEP transfer-function fit to S-star data.

Fit the TEP model (with temporal transfer functions) to the S-star data,
allowing the transfer factors S_a, T_P, D_dyn to be free parameters
constrained by the data. This is the non-isochronous refit.

The TEP model modifies the standard GR orbit by:
  1. Temporal transfer: observed periods are stretched by T_P(r)
  2. Spatial calibration: inferred distances are scaled by S_a(r)
  3. Dynamical modification: orbital dynamics modified by D_dyn(r)

The fit determines whether the data prefer TEP over GR, and what
the inferred Phantom Mass is.

MODEL:
  The TEP orbit model has additional parameters beyond the GR fit:
    - eta_TEP: TEP coupling strength (replaces eta as a free parameter)
    - A_c: central lapse (deep-interior regularisation, not relevant for S2)
    - The transfer functions are computed from the TEP temporal-well geometry

  For S2 (weak field), the TEP corrections are perturbative:
    delta_x / x ~ eta * GM/(rc^2) ~ 10^-4 * eta
    delta_v / v ~ eta * GM/(rc^2) ~ 10^-4 * eta

LIKELIHOOD:
  chi^2_TEP = sum_k [(x_obs - x_TEP)^2 / sigma_x^2 + (y_obs - y_TEP)^2 / sigma_y^2]
            + sum_k [(v_obs - v_TEP)^2 / sigma_v^2]

  The TEP model reduces to GR when eta_TEP = 0.

BAYES FACTOR:
  The comparison between GR (eta=0) and TEP (eta free) is quantified
  by the Bayes factor B = L_TEP / L_GR, computed via the Bayesian
  Information Criterion (BIC) as a first approximation.

Outputs:
  data/processed/sstar_tep_fit.json
  results/step_46_tep_transfer.json
  logs/step_46_tep_transfer.log
"""

from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

# Import from step_01
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps" / "inference"))
from step_44_gr_fit import gr_orbit_with_precession
from step_45_mass_bias import (
    temporal_transfer_factor,
    spatial_calibration_factor,
    dynamical_modification_factor,
    G_SI, c_si, M_sun_kg, AU_m,
)

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    write_json,
    finalize_result,
    tep_s2_conformal_factor,
    PROCESSED_DIR,
    RESULTS_DIR,
)

STEP_ID = "step_46_tep_transfer"


# ============================================================
# TEP orbit model
# ============================================================
def tep_orbit(params, t_arr):
    """Compute TEP-corrected orbit positions and velocities.

    The GR orbit (from gr_orbit_with_precession) ALREADY includes all
    GR effects: gravitational lensing, gravitational redshift, and
    Schwarzschild precession. The synthetic data is generated from
    this GR orbit.

    The TEP correction is the ADDITIONAL effect from the conformal
    factor A(r) = exp(+eta/(3*r_rs)).  This is the mass-inflation branch
    (negative alpha_GB), where A > 1 for eta > 0.  Spatial scales are
    magnified relative to GR, producing positive Phantom Mass.

      S_a^TEP / S_a^GR = A(r)   (conformal spatial magnification)
      T_P^TEP / T_P^GR = 1/A(r) (conformal temporal transfer)

    So the TEP-predicted observables are:
      x_TEP = A(r) * x_GR
      y_TEP = A(r) * y_GR
      v_TEP = v_GR / A(r)

    At eta=0: A=1, so TEP = GR exactly.
    """
    # Get GR orbit (includes all GR effects)
    x_gr, y_gr, v_gr = gr_orbit_with_precession(params, t_arr)

    # Compute radius at each time
    r = np.sqrt(x_gr**2 + y_gr**2)

    # Convert to Schwarzschild radii
    M_BH = params.get("M_BH_Msun", 4.297e6)
    R0 = params.get("R0_pc", 8277.0)
    eta = params.get("eta_TEP", 0.0)

    R_s_m = 2 * G_SI * M_BH * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m
    R_s_arcsec = R_s_AU / R0

    r_rs = r / R_s_arcsec  # radius in R_s units
    r_rs = np.maximum(r_rs, 1.0)  # avoid singularity

    # TEP-specific conformal factor: A(r) = exp(+eta/(3*r_rs))
    # Mass-inflation branch (alpha_GB < 0).  At eta=0, A=1 and TEP = GR.
    A = np.exp(eta / (3.0 * r_rs))

    # Apply TEP corrections (perturbative, on top of GR)
    x_tep = A * x_gr
    y_tep = A * y_gr
    v_tep = v_gr / A  # scaling of observed velocities by the conformal factor

    return x_tep, y_tep, v_tep


def tep_residuals(theta, astrometry, spectroscopy):
    """Compute residuals for the TEP fit.

    theta: [M_BH, R0, a, e, i, omega, Omega, T, t_peri, x0, y0, eta_TEP]
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
        "eta_TEP": theta[11],
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

        x_model, y_model, _ = tep_orbit(params, t_ast)
        res_ast = np.concatenate([(x_obs - x_model - x0) / sx, (y_obs - y_model - y0) / sy])
    else:
        res_ast = np.array([])

    # Spectroscopy residuals
    if len(spectroscopy) > 0:
        t_spec = np.array([d["t_yr"] for d in spectroscopy])
        v_obs = np.array([d["v_rad_kms"] for d in spectroscopy])
        sv = np.array([d["sigma_v_kms"] for d in spectroscopy])

        _, _, v_model = tep_orbit(params, t_spec)
        res_spec = (v_obs - v_model) / sv
    else:
        res_spec = np.array([])

    return np.concatenate([res_ast, res_spec])


# ============================================================
# Fitting
# ============================================================
def fit_tep_orbit(astrometry, spectroscopy, gr_fit_params):
    """Fit the TEP orbit model to the data.

    Parameters are internally scaled to O(1) and a Gaussian prior on R0
    is added to break the M_BH--R0 degeneracy, matching step_01.
    """
    from scipy.optimize import least_squares
    from scipy.linalg import pinv

    theta0 = np.array([
        gr_fit_params["M_BH_Msun"],
        gr_fit_params["R0_pc"],
        gr_fit_params["a_arcsec"],
        gr_fit_params["eccentricity"],
        gr_fit_params["inclination_deg"],
        gr_fit_params["omega_deg"],
        gr_fit_params["Omega_deg"],
        gr_fit_params["T_orb_yr"],
        gr_fit_params["t_peri_yr"],
        0.0,  # x0
        0.0,  # y0
        0.1,  # eta_TEP initial guess
    ])

    # eta_TEP is restricted to [0, 1] because TEP requires the
    # mass-inflation branch (alpha_GB < 0, A > 1, M_phantom^T > 0).
    lower = np.array([1e5, 7000, 0.05, 0.0, 0, 0, 0, 5, 1990, -1, -1, 0.0])
    upper = np.array([1e7, 10000, 0.5, 0.99, 180, 360, 360, 30, 2030, 1, 1, 1.0])

    # Characteristic scales for the 12 parameters (same as step_01 plus eta_TEP scale)
    scales = np.array([1.0e6, 1.0e3, 0.1, 0.1, 100.0, 100.0, 100.0, 10.0, 100.0, 0.01, 0.01, 1.0])

    # Soft Gaussian priors on M_BH and R0 from the GR fit values/errors.
    # These break the M_BH--eta and M_BH--R0 degeneracies so the data can
    # isolate the TEP signal rather than absorb it into the mass scale.
    M_BH_prior = gr_fit_params.get("M_BH_Msun", 4.297e6)
    sigma_M_BH = gr_fit_params.get("M_BH_err", 1.0e5)
    R0_prior = gr_fit_params.get("R0_pc", 8320.0)
    sigma_R0 = gr_fit_params.get("R0_err", 60.0)

    def tep_residuals_scaled(theta_s, astrometry, spectroscopy):
        theta = theta_s * scales
        res = tep_residuals(theta, astrometry, spectroscopy)
        M_BH_prior_penalty = (theta[0] - M_BH_prior) / sigma_M_BH
        R0_prior_penalty = (theta[1] - R0_prior) / sigma_R0
        return np.concatenate([res, [M_BH_prior_penalty, R0_prior_penalty]])

    theta0_s = theta0 / scales
    lower_s = lower / scales
    upper_s = upper / scales

    # Multi-start optimization: try non-negative eta initial values
    best_result = None
    best_cost = np.inf
    eta_starts = [0.0, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    for eta_start in eta_starts:
        theta0_try = theta0_s.copy()
        theta0_try[11] = eta_start
        try:
            result_try = least_squares(
                tep_residuals_scaled,
                theta0_try,
                args=(astrometry, spectroscopy),
                bounds=(lower_s, upper_s),
                method="trf",
                max_nfev=10000,
            )
            if result_try.cost < best_cost:
                best_cost = result_try.cost
                best_result = result_try
        except Exception as exc:
            print(f"DEBUG least_squares eta_start={eta_start}: {exc}")
            continue

    result = best_result

    theta_fit = result.x * scales

    try:
        J_s = result.jac
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cov_s = pinv(J_s.T @ J_s) * (result.cost * 2 / (len(result.fun) - len(theta0)))
        cov_phys = cov_s * np.outer(scales, scales)
        perr = np.sqrt(np.diag(cov_phys))
        # Clamp reported errors to the physical parameter range (especially for
        # parameters sitting on a bound, where the linear covariance can blow up).
        perr = np.minimum(perr, (upper - lower) * scales)
    except Exception:
        perr = np.zeros(len(theta0))

    return theta_fit, perr, result


# ============================================================
# Bayes factor (BIC approximation)
# ============================================================
def compute_bic(chi2, n_params, n_data):
    """Bayesian Information Criterion."""
    return chi2 + n_params * np.log(n_data)


def bayes_factor_bic(bic_simple, bic_complex):
    """Bayes factor from BIC difference.

    B = exp(-Delta_BIC / 2)
    Delta_BIC = BIC_complex - BIC_simple
    B > 1 favours the complex model
    """
    delta_bic = bic_complex - bic_simple
    log_b = -delta_bic / 2
    return np.exp(log_b), delta_bic


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 03 (INFERENCE): TEP TRANSFER-FUNCTION FIT", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Fitting TEP model (with transfer functions) to S-star data.", "INFO")
    print_status("Additional parameter: eta_TEP (TEP coupling strength)", "INFO")
    print_status("")

    # --- Load data and GR fit ---
    astrometry_path = PROCESSED_DIR / "sstar_astrometry.json"
    spectroscopy_path = PROCESSED_DIR / "sstar_spectroscopy.json"
    gr_fit_path = PROCESSED_DIR / "sstar_gr_fit.json"

    if not all(p.exists() for p in [astrometry_path, spectroscopy_path, gr_fit_path]):
        print_status("ERROR: Data or GR fit not found. Run steps 00-01 first.", "ERROR")
        return

    astrometry_data = json.loads(astrometry_path.read_text())
    spectroscopy_data = json.loads(spectroscopy_path.read_text())
    gr_fit = json.loads(gr_fit_path.read_text())

    astrometry = astrometry_data["data"]
    spectroscopy = spectroscopy_data["data"]

    # Extract GR fit parameters
    gr_params = {
        "M_BH_Msun": gr_fit["best_fit"]["M_BH (M_sun)"]["value"],
        "M_BH_err": gr_fit["best_fit"]["M_BH (M_sun)"]["error"],
        "R0_pc": gr_fit["best_fit"]["R0 (pc)"]["value"],
        "R0_err": gr_fit["best_fit"]["R0 (pc)"]["error"],
        "a_arcsec": gr_fit["best_fit"]["a (arcsec)"]["value"],
        "eccentricity": gr_fit["best_fit"]["e"]["value"],
        "inclination_deg": gr_fit["best_fit"]["i (deg)"]["value"],
        "omega_deg": gr_fit["best_fit"]["omega (deg)"]["value"],
        "Omega_deg": gr_fit["best_fit"]["Omega (deg)"]["value"],
        "T_orb_yr": gr_fit["best_fit"]["T (yr)"]["value"],
        "t_peri_yr": gr_fit["best_fit"]["t_peri (yr)"]["value"],
    }

    is_synthetic = astrometry_data.get("is_synthetic", False)
    print_status(f"  Data: {'SYNTHETIC' if is_synthetic else 'REAL'}", "INFO")
    print_status(f"  Astrometry: {len(astrometry)} points", "INFO")
    print_status(f"  Spectroscopy: {len(spectroscopy)} points", "INFO")
    print_status("")

    # --- Fit TEP model ---
    print_status("--- Fitting TEP orbit model ---", "INFO")
    print_status("  Parameters: M_BH, R0, a, e, i, omega, Omega, T, t_peri, x0, y0, eta_TEP", "INFO")
    print_status("  Method: scipy.optimize.least_squares (TRF)", "INFO")
    print_status("")

    best_fit, perr, result = fit_tep_orbit(astrometry, spectroscopy, gr_params)

    param_names = ["M_BH (M_sun)", "R0 (pc)", "a (arcsec)", "e", "i (deg)",
                   "omega (deg)", "Omega (deg)", "T (yr)", "t_peri (yr)",
                   "x0 (arcsec)", "y0 (arcsec)", "eta_TEP"]

    print_status("--- TEP fit results ---", "TITLE")
    for name, val, err in zip(param_names, best_fit, perr):
        print_status(f"  {name:20s} = {val:12.6f} ± {err:10.6f}", "INFO")
    print_status("")

    # Compute chi-squared
    residuals = tep_residuals(best_fit, astrometry, spectroscopy)
    chi2_tep = np.sum(residuals**2)
    dof_tep = len(residuals) - len(best_fit)
    chi2_red_tep = chi2_tep / dof_tep if dof_tep > 0 else float("inf")

    print_status(f"  chi^2_TEP = {chi2_tep:.2f}", "INFO")
    print_status(f"  dof = {dof_tep}", "INFO")
    print_status(f"  chi^2/dof = {chi2_red_tep:.4f}", "INFO")
    print_status("")

    # --- Compare to GR ---
    chi2_gr = gr_fit["chi2"]
    n_params_gr = 11
    n_params_tep = 12
    n_data = len(residuals)

    bic_gr = compute_bic(chi2_gr, n_params_gr, n_data)
    bic_tep = compute_bic(chi2_tep, n_params_tep, n_data)

    bf, delta_bic = bayes_factor_bic(bic_gr, bic_tep)

    print_status("--- GR vs TEP comparison ---", "TITLE")
    print_status(f"  chi^2_GR  = {chi2_gr:.2f}  (n_params = {n_params_gr})", "INFO")
    print_status(f"  chi^2_TEP = {chi2_tep:.2f}  (n_params = {n_params_tep})", "INFO")
    print_status(f"  Delta chi^2 = {chi2_gr - chi2_tep:.2f}", "INFO")
    print_status(f"  BIC_GR  = {bic_gr:.2f}", "INFO")
    print_status(f"  BIC_TEP = {bic_tep:.2f}", "INFO")
    print_status(f"  Delta BIC = {delta_bic:.2f}", "INFO")
    print_status(f"  Bayes factor B(TEP/GR) = {bf:.4e}", "INFO")
    print_status("")

    if bf > 1:
        print_status("  -> TEP model is FAVOURED over GR (B > 1)", "INFO")
    else:
        print_status("  -> GR model is FAVOURED over TEP (B < 1)", "INFO")
        print_status(f"  -> eta_TEP is consistent with 0 within {perr[11]:.4f}", "INFO")

    print_status("")

    # --- Compute Phantom Mass ---
    eta_best = best_fit[11]
    M_BH_tep = best_fit[0]
    M_BH_gr = gr_params["M_BH_Msun"]

    # M_phantom = M_GR - M_TEP (residual definition)
    M_phantom = M_BH_gr - M_BH_tep
    M_phantom_frac = M_phantom / M_BH_gr

    print_status("--- Phantom Mass estimate ---", "TITLE")
    print_status(f"  M_BH^GR  = {M_BH_gr:.4e} M_sun", "INFO")
    print_status(f"  M_BH^TEP = {M_BH_tep:.4e} M_sun", "INFO")
    print_status(f"  M_phantom^T = M_GR - M_TEP = {M_phantom:.4e} M_sun", "INFO")
    print_status(f"  M_phantom fraction = {M_phantom_frac:.6e}", "INFO")
    print_status("")

    # --- Save ---
    fit_result = {
        "step": STEP_ID,
        "description": "TEP transfer-function fit to S-star data",
        "is_synthetic_data": is_synthetic,
        "best_fit": {name: {"value": float(val), "error": float(err)}
                     for name, val, err in zip(param_names, best_fit, perr)},
        "chi2": float(chi2_tep),
        "dof": int(dof_tep),
        "chi2_reduced": float(chi2_red_tep),
        "comparison": {
            "chi2_GR": float(chi2_gr),
            "chi2_TEP": float(chi2_tep),
            "delta_chi2": float(chi2_gr - chi2_tep),
            "BIC_GR": float(bic_gr),
            "BIC_TEP": float(bic_tep),
            "delta_BIC": float(delta_bic),
            "bayes_factor": float(bf),
            "n_params_GR": n_params_gr,
            "n_params_TEP": n_params_tep,
        },
        "phantom_mass": {
            "M_GR_Msun": float(M_BH_gr),
            "M_TEP_Msun": float(M_BH_tep),
            "M_phantom_Msun": float(M_phantom),
            "M_phantom_fraction": float(M_phantom_frac),
            "eta_TEP": float(eta_best),
            "eta_TEP_err": float(perr[11]),
        },
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    fit_path = PROCESSED_DIR / "sstar_tep_fit.json"
    write_json(fit_path, fit_result)
    print_status(f"  Saved: {fit_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, fit_result)
    finalize_result(STEP_ID, fit_result, "TEP transfer-function fit", "TEP model fitted with eta_TEP as free parameter")

    print_status("")
    print_status("  NEXT: step_47_joint_forward.py — joint forward-model", "INFO")


if __name__ == "__main__":
    main()
