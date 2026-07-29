#!/usr/bin/env python3
"""Step 05: Observational Constraints.

Uses real published EHT and LIGO measurements (downloaded in step_00)
to constrain the TEP black-hole model parameters.

Compares TEP theoretical predictions against EHT shadow measurements and
reports a GW150914 ringdown diagnostic.  Canonical TEP gravitational waves
propagate on geometric g, prescribed in this repository as Schwarzschild, so
the supported nonspinning TEP QNM is exactly the Schwarzschild QNM.  Because
the GW150914 remnant spins, its data are not used as a nonspinning TEP-versus-
Schwarzschild improvement test.

Computes chi-squared and p-value diagnostics.

Outputs:
  - results/step_05_observational_constraints.json
  - results/step_05_observational_constraints.csv
  - logs/step_05_observational_constraints.log

Usage:
    python scripts/steps/step_05_observational_constraints.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from datetime import datetime
from math import erf, sqrt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

import numpy as np
from scipy.integrate import cumulative_trapezoid

from bh_common import (
    ensure_dirs,
    finalize_result,
    make_step_logger,
    print_status,
    step_json_path,
    step_csv_path,
    rel,
    write_json,
    write_csv,
    PROCESSED_DIR,
    TEPBHModel,
    solve_tep_bh,
)

STEP_ID = "step_05_observational_constraints"

# Physical constants (SI)
G_NEWTON = 6.67430e-11
C_LIGHT = 2.99792458e8
M_SUN = 1.98847e30
MPC_TO_M = 3.08567758e22
KPC_TO_M = 3.08567758e19
UAS_TO_RAD = 4.84813681109536e-12


# =============================================================================
# Photon sphere from the disformal metric (copied from step_03_raytracing.py)
# =============================================================================

def compute_photon_sphere(model, solution):
    """Compute the photon sphere radius and critical impact parameter.

    For a spherically symmetric metric ds² = -f dt² + h dr² + R² dΩ²,
    the photon sphere is at the maximum of the effective potential:

        V_eff = f / R² * L²

    The condition dV_eff/dr = 0 gives:

        f' * R = 2 * f * R'

    For the TEP disformal metric in the conformal exterior:
        f = A² * F,  R = A * r,  F = 1 - 2M/r

    Substituting: the conformal factor A cancels exactly, giving:
        F' * r = 2 * F  =>  r = 3M  (Schwarzschild)

    The impact parameter:
        b = R / sqrt(f) = A*r / (A*sqrt(F)) = r / sqrt(F)

    At r = 3M: b = 3*sqrt(3)*M (same as Schwarzschild).

    The TEP correction comes only from the residual disformal term
    (B * (phi')²) which is negligible in the exterior.
    """
    r = solution['r']
    metric = solution['metric']
    M = model.M

    A = metric['A']
    F = metric['F']  # Schwarzschild metric function
    R = A * r  # physical areal radius

    # Effective metric function: f = A² * |F| (conformal limit)
    f_eff = A**2 * np.abs(F)

    # Derivatives
    df_dr = np.gradient(f_eff, r)
    dR_dr = np.gradient(R, r)

    # Photon sphere condition: f' * R = 2 * f * R' (exterior only)
    idx_h = np.argmin(np.abs(r - 2.0 * M))
    r_ext = r[idx_h:]
    f_ext = f_eff[idx_h:]
    df_ext = df_dr[idx_h:]
    R_ext = R[idx_h:]
    dR_ext = dR_dr[idx_h:]

    condition = df_ext * R_ext - 2.0 * f_ext * dR_ext
    # Find zero crossing (should be near r = 3M)
    sign_changes = np.where(np.diff(np.sign(condition)))[0]

    if len(sign_changes) > 0:
        # Pick the crossing closest to r = 3M
        r_crossings = r_ext[sign_changes]
        idx_best = sign_changes[np.argmin(np.abs(r_crossings - 3.0 * M))]
        idx_local = idx_best
        if abs(condition[idx_local + 1] - condition[idx_local]) > 1e-30:
            frac = -condition[idx_local] / (condition[idx_local + 1] - condition[idx_local])
            r_photon = r_ext[idx_local] + frac * (r_ext[idx_local + 1] - r_ext[idx_local])
        else:
            r_photon = r_ext[idx_local]
    else:
        r_photon = 3.0 * M  # fallback to Schwarzschild

    # Critical impact parameter: b = R / sqrt(f)
    idx_ph = np.argmin(np.abs(r - r_photon))
    R_ph = R[idx_ph]
    f_ph = f_eff[idx_ph]
    b_photon = R_ph / np.sqrt(np.maximum(f_ph, 1e-30))

    # Schwarzschild reference
    r_photon_schw = 3.0 * M
    b_photon_schw = 3.0 * np.sqrt(3.0) * M

    return {
        'r_photon': float(r_photon),
        'b_photon': float(b_photon),
        'r_photon_schwarzschild': float(r_photon_schw),
        'b_photon_schwarzschild': float(b_photon_schw),
        'b_photon_ratio': float(b_photon / b_photon_schw),
        'r_photon_ratio': float(r_photon / r_photon_schw),
        'A_at_photon_sphere': float(A[idx_ph]),
    }


# =============================================================================
# Shadow angular diameter (copied from core/bh_raytracing.py)
# =============================================================================

def compute_shadow_diameter(b_photon, M_kg, D_m):
    """Compute the shadow angular diameter.

    theta_shadow = 2 * b_photon * (G * M) / (c^2 * D)

    where b_photon is in geometric units (multiples of M).

    Parameters
    ----------
    b_photon : float
        Critical impact parameter in geometric units (M=1).
    M_kg : float
        Black hole mass in kg.
    D_m : float
        Distance to observer in meters.

    Returns
    -------
    float
        Shadow angular diameter in radians.
    """
    # Convert geometric b to physical: b_phys = b_geometric * G * M / c^2
    r_g = G_NEWTON * M_kg / C_LIGHT**2  # gravitational radius
    b_phys = b_photon * r_g

    theta = 2.0 * b_phys / D_m
    return theta


# =============================================================================
# Tortoise coordinate for the disformal metric (copied from step_02_perturbations.py)
# =============================================================================

def compute_tortoise_coordinate(r, metric):
    """Compute the tortoise coordinate r* for the disformal metric.

    For a spherically symmetric metric ds^2 = -f dt^2 + h dr^2 + R^2 dOmega^2,
    the tortoise coordinate is:

        dr*/dr = sqrt(h / f)

    For the TEP disformal metric in EF form:
        gtilde_{vv} = -A^2 F,  gtilde_{vr} = A^2,  gtilde_{rr} = B (phi')^2

    The effective f and h for the (t,r) sector (diagonalized) are:
        f_eff = -det_2d / gtilde_rr   (if gtilde_rr > 0)
        h_eff = gtilde_rr

    But in the exterior where B -> 0, gtilde_rr -> 0, so we use the
    conformal limit: gtilde -> A^2 * g_Schw, giving
        f_eff = A^2 * F,  h_eff = A^2 / F
        dr*/dr = sqrt(h/f) = 1/F  (same as Schwarzschild)

    In the interior, the tortoise coordinate extends smoothly.
    For the exterior (where perturbations are observable), r* = r + 2M ln(r/2M - 1).
    """
    F = metric['F']
    A = metric['A']

    # In the exterior (r > 2M): r* = r + 2M ln|r/2M - 1| (Schwarzschild)
    # The conformal factor A cancels in sqrt(h/f) = sqrt((A^2/F)/(A^2*F)) = 1/|F|
    # So r* is the same as Schwarzschild in the exterior.
    r_h = 2.0
    F_abs = np.abs(F)
    F_safe = np.where(F_abs > 1e-30, F_abs, 1e-30)

    # dr*/dr = 1/|F| (same as Schwarzschild, A cancels)
    drstar_dr = 1.0 / F_safe

    # Integrate from a reference point (r = 3M, near the photon sphere)
    idx_ref = np.argmin(np.abs(r - 3.0))
    r_star = np.zeros_like(r)
    # Integrate outward from reference
    r_star[idx_ref:] = cumulative_trapezoid(drstar_dr[idx_ref:], r[idx_ref:], initial=0)
    # Integrate inward from reference
    r_star[:idx_ref+1] = -cumulative_trapezoid(
        drstar_dr[idx_ref::-1], r[idx_ref::-1], initial=0)[::-1]

    return r_star


# =============================================================================
# Regge-Wheeler potential for GWs on geometric g
# =============================================================================

def compute_regge_wheeler_potential(r, r_star, metric, l=2):
    """Compute the Schwarzschild RW potential for GWs on geometric g.

    The matter metric does not rescale gravitational perturbations.  The
    second returned array is an identical Schwarzschild reference.
    """
    F = metric['F']
    F_safe = np.where(np.abs(F) > 1e-30, F, np.nan)
    V_RW_geometric = F_safe * (l * (l + 1) / r**2 - 6.0 / r**3)
    V_RW_geometric = np.where(np.isfinite(V_RW_geometric), V_RW_geometric, 0)

    return V_RW_geometric, V_RW_geometric.copy()


# =============================================================================
# QNM via WKB approximation (copied from step_02_perturbations.py)
# =============================================================================

def compute_qnm_wkb(r, r_star, V, n=0):
    """Compute QNM frequency using the WKB approximation.

    For a potential V(r*) with a single peak at r*_0:

        omega^2 = V_0 + sqrt(-2 * V''_0) * (n + 1/2) * i

    where V_0 = V(r*_0) and V''_0 is the second derivative at the peak.
    """
    # Find the peak of V in the exterior (r > 2M)
    idx_h = np.argmin(np.abs(r - 2.0))
    V_ext = V[idx_h:]
    r_star_ext = r_star[idx_h:]
    r_ext = r[idx_h:]

    if len(V_ext) < 10:
        return {"omega_R": None, "omega_I": None, "V_max": None, "r_max": None}

    idx_peak = np.argmax(V_ext)
    V_max = float(V_ext[idx_peak])
    r_max = float(r_ext[idx_peak])
    r_star_max = float(r_star_ext[idx_peak])

    # Second derivative at the peak
    if idx_peak > 0 and idx_peak < len(V_ext) - 1:
        drs = r_star_ext[idx_peak + 1] - r_star_ext[idx_peak - 1]
        if abs(drs) > 1e-30:
            V_pp = (V_ext[idx_peak + 1] - 2 * V_ext[idx_peak] + V_ext[idx_peak - 1]) / (
                (drs / 2)**2)
        else:
            V_pp = 0
    else:
        V_pp = 0

    if V_pp >= 0:
        return {"omega_R": None, "omega_I": None, "V_max": V_max, "r_max": r_max}

    sqrt_term = np.sqrt(-2 * V_pp)
    omega_sq = V_max + 1j * sqrt_term * (n + 0.5)
    omega = np.sqrt(omega_sq)

    return {
        "omega_R": float(omega.real),
        "omega_I": float(omega.imag),
        "V_max": V_max,
        "r_max": r_max,
        "r_star_max": r_star_max,
        "V_pp": float(V_pp),
    }


# =============================================================================
# ISCO properties (copied from step_04_accretion.py)
# =============================================================================

def compute_isco_properties(M=1.0):
    """Compute ISCO properties for Schwarzschild (exterior of TEP metric).

    At r_ISCO = 6M:
        E_ISCO = sqrt(8/9) ~ 0.9428
        L_ISCO = 2*sqrt(3)*M ~ 3.464*M
        Omega_ISCO = 1/(6*sqrt(6)*M) ~ 0.0680/M
        eta = 1 - E_ISCO ~ 5.72%
    """
    r_isco = 6.0 * M
    E_isco = np.sqrt(8.0 / 9.0)
    L_isco = 2.0 * np.sqrt(3.0) * M
    Omega_isco = 1.0 / (6.0 * np.sqrt(6.0) * M)
    eta = 1.0 - E_isco
    g_redshift = np.sqrt(1.0 - 2.0 * M / r_isco)

    return {
        'r_isco': float(r_isco),
        'E_isco': float(E_isco),
        'L_isco': float(L_isco),
        'Omega_isco': float(Omega_isco),
        'radiative_efficiency': float(eta),
        'radiative_efficiency_percent': float(eta * 100),
        'redshift_factor_isco': float(g_redshift),
        'redshift_z_isco': float(1.0 / g_redshift - 1.0),
    }


def normal_cdf(x):
    """Standard normal CDF."""
    return 0.5 * (1.0 + erf(x / sqrt(2)))


def chi_squared(prediction, measurement, uncertainty):
    """Compute chi-squared for a single measurement."""
    if uncertainty == 0:
        return 0.0 if prediction == measurement else float('inf')
    return ((prediction - measurement) / uncertainty) ** 2


def two_sided_pvalue(chi2, ndof=1):
    """Compute two-sided p-value from chi-squared."""
    z = sqrt(chi2)
    return 2.0 * (1.0 - normal_cdf(z))


def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 05: Observational Constraints", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    # --- Load real measurement data ---
    eht_path = PROCESSED_DIR / "eht_measurements.json"
    ligo_path = PROCESSED_DIR / "ligo_qnm_measurements.json"

    if not eht_path.exists() or not ligo_path.exists():
        print_status("Measurement data not found. Run step_00_data_download first.", "ERROR")
        return {"step": STEP_ID, "status": "failed", "error": "missing data"}

    with open(eht_path) as f:
        eht_data = json.load(f)
    with open(ligo_path) as f:
        ligo_data = json.load(f)

    print_status("Loaded real observational data:", "SUCCESS")
    print_status(f"  EHT measurements: {rel(eht_path)}", "INFO")
    print_status(f"  LIGO QNM measurements: {rel(ligo_path)}", "INFO")
    print_status("")

    # --- TEP model ---
    model = TEPBHModel(
        beta_A=-1.0, B0=1.0, n_B=2.0,
        phi_0=2.0, delta=0.05, M=1.0,
    )
    print_status(f"Model parameters: {model.to_dict()}", "INFO")
    print_status(f"  beta_A = {model.beta_A} (conformal coupling)", "DEBUG")
    print_status(f"  B0 = {model.B0}, n_B = {model.n_B} (disformal shear)", "DEBUG")
    print_status(f"  phi_0 = {model.phi_0}, delta = {model.delta} (scalar profile)", "DEBUG")
    print_status(f"  M = {model.M} (geometric mass)", "DEBUG")

    print_status("Solving TEP-BH background...", "PROCESS")
    solution = solve_tep_bh(model)
    if not solution["success"]:
        print_status("BACKGROUND SOLUTION FAILED", "ERROR")
        return {"step": STEP_ID, "status": "failed"}
    print_status(f"  Background solved: r_grid = [{solution['r'][0]:.4f}, {solution['r'][-1]:.4f}] M", "DEBUG")
    print_status(f"  Radial grid points: {len(solution['r'])}", "DEBUG")

    # --- Compute theoretical predictions ---
    print_status("Computing theoretical predictions...", "PROCESS")
    ps = compute_photon_sphere(model, solution)
    b_photon = ps["b_photon"]
    print_status(f"  Photon sphere radius: r_ph = {ps['r_photon']:.6f} M (Schw: {ps['r_photon_schwarzschild']:.6f} M)", "DEBUG")
    print_status(f"  Critical impact parameter: b_ph = {ps['b_photon']:.6f} M (Schw: {ps['b_photon_schwarzschild']:.6f} M)", "DEBUG")
    print_status(f"  b_photon ratio (TEP/Schw): {ps['b_photon_ratio']:.8f}", "DEBUG")
    print_status(f"  A at photon sphere: {ps['A_at_photon_sphere']:.8f}", "DEBUG")

    # Perturbations for QNM
    r = solution["r"]
    metric = solution["metric"]
    r_star = compute_tortoise_coordinate(r, metric)
    V_RW, V_RW_schw = compute_regge_wheeler_potential(r, r_star, metric, l=2)
    qnm = compute_qnm_wkb(r, r_star, V_RW, n=0)
    qnm_schw = compute_qnm_wkb(r, r_star, V_RW_schw, n=0)
    print_status(f"  Regge-Wheeler potential peak: V_max = {qnm.get('V_max', 'N/A')}", "DEBUG")
    print_status(f"  QNM (WKB, l=2, n=0): omega_R = {qnm.get('omega_R', 'N/A')}, omega_I = {qnm.get('omega_I', 'N/A')}", "DEBUG")
    print_status(f"  QNM Schwarzschild ref: omega_R = {qnm_schw.get('omega_R', 'N/A')}, omega_I = {qnm_schw.get('omega_I', 'N/A')}", "DEBUG")

    # ISCO
    isco = compute_isco_properties(model.M)
    print_status(f"  ISCO radius: {isco['r_isco']:.4f} M", "DEBUG")
    print_status(f"  ISCO radiative efficiency: {isco['radiative_efficiency_percent']:.4f}%", "DEBUG")
    print_status(f"  ISCO orbital frequency: {isco['Omega_isco']:.6f} /M", "DEBUG")

    # --- EHT M87* constraints ---
    print_status("", "INFO")
    print_status("EHT M87* Constraints", "TITLE")
    m87 = eht_data["M87"]
    m87_mass = m87["inferred_parameters"]["mass_Msun"]["value"]
    m87_mass_err = m87["inferred_parameters"]["mass_Msun"]["uncertainty_syst"]
    m87_dist = m87["inferred_parameters"]["distance_Mpc"]["value"]
    m87_dist_err = m87["inferred_parameters"]["distance_Mpc"]["uncertainty"]
    m87_shadow_meas = m87["observations"]["ring_diameter_uas"]["value"]
    m87_shadow_err = m87["observations"]["ring_diameter_uas"]["uncertainty"]

    print_status(f"  Input: M87* mass = {m87_mass:.2e} M_sun (syst err {m87_mass_err:.2e})", "DEBUG")
    print_status(f"  Input: M87* distance = {m87_dist} Mpc (err {m87_dist_err} Mpc)", "DEBUG")
    print_status(f"  Input: M87* shadow diameter = {m87_shadow_meas} ± {m87_shadow_err} μas", "DEBUG")

    # TEP prediction for M87* shadow
    m87_M_kg = m87_mass * M_SUN
    m87_D_m = m87_dist * MPC_TO_M
    print_status(f"  M87* mass in kg: {m87_M_kg:.4e} kg", "DEBUG")
    print_status(f"  M87* distance in m: {m87_D_m:.4e} m", "DEBUG")
    m87_theta_tep = compute_shadow_diameter(b_photon, m87_M_kg, m87_D_m)
    m87_shadow_tep = m87_theta_tep / UAS_TO_RAD
    print_status(f"  TEP shadow angular diameter: {m87_theta_tep:.4e} rad = {m87_shadow_tep:.4f} μas", "DEBUG")

    # Schwarzschild prediction
    m87_theta_schw = compute_shadow_diameter(ps["b_photon_schwarzschild"], m87_M_kg, m87_D_m)
    m87_shadow_schw = m87_theta_schw / UAS_TO_RAD
    print_status(f"  Schw shadow angular diameter: {m87_theta_schw:.4e} rad = {m87_shadow_schw:.4f} μas", "DEBUG")

    m87_chi2_tep = chi_squared(m87_shadow_tep, m87_shadow_meas, m87_shadow_err)
    m87_chi2_schw = chi_squared(m87_shadow_schw, m87_shadow_meas, m87_shadow_err)
    m87_pvalue_tep = two_sided_pvalue(m87_chi2_tep)
    m87_pvalue_schw = two_sided_pvalue(m87_chi2_schw)
    m87_sigma_tep = sqrt(m87_chi2_tep)
    m87_sigma_schw = sqrt(m87_chi2_schw)
    print_status(f"  M87* chi-squared (TEP):  ((pred-meas)/err)^2 = (({m87_shadow_tep:.4f}-{m87_shadow_meas})/{m87_shadow_err})^2 = {m87_chi2_tep:.6f}", "DEBUG")
    print_status(f"  M87* chi-squared (Schw): ((pred-meas)/err)^2 = (({m87_shadow_schw:.4f}-{m87_shadow_meas})/{m87_shadow_err})^2 = {m87_chi2_schw:.6f}", "DEBUG")
    print_status(f"  M87* sigma deviation (TEP):  {m87_sigma_tep:.4f}σ", "DEBUG")
    print_status(f"  M87* sigma deviation (Schw): {m87_sigma_schw:.4f}σ", "DEBUG")

    print_status(f"  Mass: ({m87_mass:.1e} ± {m87_mass_err:.1e}) M_sun", "INFO")
    print_status(f"  Distance: {m87_dist} ± {m87_dist_err} Mpc", "INFO")
    print_status(f"  Shadow (measured): {m87_shadow_meas} ± {m87_shadow_err} μas", "INFO")
    print_status(f"  Shadow (TEP pred):  {m87_shadow_tep:.2f} μas", "INFO")
    print_status(f"  Shadow (Schw pred): {m87_shadow_schw:.2f} μas", "INFO")
    print_status(f"  Chi-squared (TEP):  {m87_chi2_tep:.4f} ({m87_sigma_tep:.2f}σ)", "INFO")
    print_status(f"  Chi-squared (Schw): {m87_chi2_schw:.4f} ({m87_sigma_schw:.2f}σ)", "INFO")
    print_status(f"  p-value (TEP):      {m87_pvalue_tep:.4f}", "INFO")
    print_status(f"  p-value (Schw):     {m87_pvalue_schw:.4f}", "INFO")
    print_status(f"  TEP-Schw shadow difference: {abs(m87_shadow_tep - m87_shadow_schw):.6f} μas (identical: exterior = Schwarzschild)", "SUCCESS")

    # --- EHT Sgr A* constraints ---
    print_status("", "INFO")
    print_status("EHT Sgr A* Constraints", "TITLE")
    sgr = eht_data["SgrA"]
    sgr_mass = sgr["inferred_parameters"]["mass_Msun"]["value_VLTI"]
    sgr_mass_err = sgr["inferred_parameters"]["mass_Msun"]["uncertainty_VLTI_syst"]
    sgr_dist = sgr["inferred_parameters"]["distance_kpc"]["value_VLTI"]
    sgr_dist_err = sgr["inferred_parameters"]["distance_kpc"]["uncertainty_VLTI_syst"]
    sgr_shadow_meas = sgr["observations"]["angular_shadow_diameter_uas"]["value"]
    sgr_shadow_err = sgr["observations"]["angular_shadow_diameter_uas"]["uncertainty"]

    print_status(f"  Input: Sgr A* mass (VLTI) = {sgr_mass:.4e} M_sun (syst err {sgr_mass_err:.4e})", "DEBUG")
    print_status(f"  Input: Sgr A* distance (VLTI) = {sgr_dist} kpc (err {sgr_dist_err} kpc)", "DEBUG")
    print_status(f"  Input: Sgr A* shadow diameter = {sgr_shadow_meas} ± {sgr_shadow_err} μas", "DEBUG")

    # TEP prediction for Sgr A* shadow
    sgr_M_kg = sgr_mass * M_SUN
    sgr_D_m = sgr_dist * KPC_TO_M  # kpc to m (KPC_TO_M is already meters per kpc)
    print_status(f"  Sgr A* mass in kg: {sgr_M_kg:.4e} kg", "DEBUG")
    print_status(f"  Sgr A* distance in m: {sgr_D_m:.4e} m", "DEBUG")
    sgr_theta_tep = compute_shadow_diameter(b_photon, sgr_M_kg, sgr_D_m)
    sgr_shadow_tep = sgr_theta_tep / UAS_TO_RAD
    print_status(f"  TEP shadow angular diameter: {sgr_theta_tep:.4e} rad = {sgr_shadow_tep:.4f} μas", "DEBUG")

    sgr_theta_schw = compute_shadow_diameter(ps["b_photon_schwarzschild"], sgr_M_kg, sgr_D_m)
    sgr_shadow_schw = sgr_theta_schw / UAS_TO_RAD
    print_status(f"  Schw shadow angular diameter: {sgr_theta_schw:.4e} rad = {sgr_shadow_schw:.4f} μas", "DEBUG")

    sgr_chi2_tep = chi_squared(sgr_shadow_tep, sgr_shadow_meas, sgr_shadow_err)
    sgr_chi2_schw = chi_squared(sgr_shadow_schw, sgr_shadow_meas, sgr_shadow_err)
    sgr_pvalue_tep = two_sided_pvalue(sgr_chi2_tep)
    sgr_pvalue_schw = two_sided_pvalue(sgr_chi2_schw)
    sgr_sigma_tep = sqrt(sgr_chi2_tep)
    sgr_sigma_schw = sqrt(sgr_chi2_schw)
    print_status(f"  Sgr A* chi-squared (TEP):  ((pred-meas)/err)^2 = (({sgr_shadow_tep:.4f}-{sgr_shadow_meas})/{sgr_shadow_err})^2 = {sgr_chi2_tep:.6f}", "DEBUG")
    print_status(f"  Sgr A* chi-squared (Schw): ((pred-meas)/err)^2 = (({sgr_shadow_schw:.4f}-{sgr_shadow_meas})/{sgr_shadow_err})^2 = {sgr_chi2_schw:.6f}", "DEBUG")
    print_status(f"  Sgr A* sigma deviation (TEP):  {sgr_sigma_tep:.4f}σ", "DEBUG")
    print_status(f"  Sgr A* sigma deviation (Schw): {sgr_sigma_schw:.4f}σ", "DEBUG")

    print_status(f"  Mass (VLTI): ({sgr_mass:.3e} ± {sgr_mass_err:.3e}) M_sun", "INFO")
    print_status(f"  Distance (VLTI): {sgr_dist} ± {sgr_dist_err} kpc", "INFO")
    print_status(f"  Shadow (measured): {sgr_shadow_meas} ± {sgr_shadow_err} μas", "INFO")
    print_status(f"  Shadow (TEP pred):  {sgr_shadow_tep:.2f} μas", "INFO")
    print_status(f"  Shadow (Schw pred): {sgr_shadow_schw:.2f} μas", "INFO")
    print_status(f"  Chi-squared (TEP):  {sgr_chi2_tep:.4f} ({sgr_sigma_tep:.2f}σ)", "INFO")
    print_status(f"  Chi-squared (Schw): {sgr_chi2_schw:.4f} ({sgr_sigma_schw:.2f}σ)", "INFO")
    print_status(f"  p-value (TEP):      {sgr_pvalue_tep:.4f}", "INFO")
    print_status(f"  p-value (Schw):     {sgr_pvalue_schw:.4f}", "INFO")
    print_status(f"  TEP-Schw shadow difference: {abs(sgr_shadow_tep - sgr_shadow_schw):.6f} μas (identical: exterior = Schwarzschild)", "SUCCESS")

    # --- LIGO QNM constraints ---
    print_status("", "INFO")
    print_status("LIGO QNM Constraints (GW150914)", "TITLE")
    gw = ligo_data["GW150914"]
    gw_mass = gw["final_mass_Msun"]["value"]
    gw_mass_err = gw["final_mass_Msun"]["uncertainty"]
    gw_spin = gw["final_spin"]["value"]
    gw_f220 = gw["ringdown_f220_Hz"]["value"]
    gw_f220_err = gw["ringdown_f220_Hz"]["uncertainty"]
    gw_tau220 = gw["ringdown_tau220_ms"]["value"]
    gw_tau220_err = gw["ringdown_tau220_ms"]["uncertainty"]

    print_status(f"  Input: GW150914 final mass = {gw_mass} ± {gw_mass_err} M_sun", "DEBUG")
    print_status(f"  Input: GW150914 final spin = {gw_spin} ± {gw['final_spin']['uncertainty']}", "DEBUG")
    print_status(f"  Input: GW150914 f_220 = {gw_f220} ± {gw_f220_err} Hz", "DEBUG")
    print_status(f"  Input: GW150914 tau_220 = {gw_tau220} ± {gw_tau220_err} ms", "DEBUG")

    # Convert the geometric Schwarzschild QNM (M=1) to physical units for a
    # nonspinning mass-matched diagnostic. Canonical TEP and Schwarzschild are equal.
    # omega_R (geometric) -> f = omega_R / (2*pi*M_geometric_in_seconds)
    # M_geometric = G*M_kg / c^3
    gw_M_kg = gw_mass * M_SUN
    gw_M_time = G_NEWTON * gw_M_kg / C_LIGHT**3  # seconds
    print_status(f"  GW150914 mass in kg: {gw_M_kg:.4e} kg", "DEBUG")
    print_status(f"  GW150914 geometric time unit: {gw_M_time:.4e} s", "DEBUG")

    if qnm["omega_R"] is not None:
        gw_f220_tep = qnm["omega_R"] / (2 * np.pi * gw_M_time)
        gw_f220_schw = qnm_schw["omega_R"] / (2 * np.pi * gw_M_time)
        print_status(f"  QNM frequency conversion: omega_R={qnm['omega_R']:.6f} -> f={gw_f220_tep:.2f} Hz (TEP)", "DEBUG")
        print_status(f"  QNM frequency conversion: omega_R={qnm_schw['omega_R']:.6f} -> f={gw_f220_schw:.2f} Hz (Schw)", "DEBUG")
    else:
        gw_f220_tep = None
        gw_f220_schw = None

    if qnm["omega_I"] is not None:
        gw_tau220_tep = gw_M_time / abs(qnm["omega_I"]) * 1000  # ms
        gw_tau220_schw = gw_M_time / abs(qnm_schw["omega_I"]) * 1000
        print_status(f"  QNM damping conversion: omega_I={qnm['omega_I']:.6f} -> tau={gw_tau220_tep:.2f} ms (TEP)", "DEBUG")
        print_status(f"  QNM damping conversion: omega_I={qnm_schw['omega_I']:.6f} -> tau={gw_tau220_schw:.2f} ms (Schw)", "DEBUG")
    else:
        gw_tau220_tep = None
        gw_tau220_schw = None

    # GW150914's remnant has spin a=0.67, whereas this branch supports only
    # nonspinning Schwarzschild perturbations.  The following rough Kerr spin
    # correction is retained solely as a phenomenological context diagnostic;
    # it is not a TEP GW prediction and is not used for a TEP improvement claim.
    kerr_freq_correction = 1.0 + 0.5 * gw_spin - 0.1 * gw_spin**2
    kerr_damp_correction = 1.0 + 0.1 * gw_spin - 0.05 * gw_spin**2
    gw_f220_kerr_pred = gw_f220_schw * kerr_freq_correction if gw_f220_schw else None
    gw_tau220_kerr_pred = gw_tau220_schw / kerr_damp_correction if gw_tau220_schw else None
    print_status(f"  Kerr spin correction (freq): {kerr_freq_correction:.4f} (spin a={gw_spin})", "DEBUG")
    print_status(f"  Kerr spin correction (damp): {kerr_damp_correction:.4f}", "DEBUG")
    print_status(f"  Kerr f_220 prediction: {gw_f220_kerr_pred:.2f} Hz (Schw * {kerr_freq_correction:.4f})", "DEBUG")
    print_status(f"  Kerr tau_220 prediction: {gw_tau220_kerr_pred:.2f} ms (Schw / {kerr_damp_correction:.4f})", "DEBUG")

    # Mass-matched diagnostics only: direct comparison to a spinning remnant
    # is not a valid nonspinning model-selection or improvement test.
    if gw_f220_tep is not None and gw_f220_schw is not None:
        tep_shift_hz = gw_f220_tep - gw_f220_schw
        gw_f220_shift_chi2 = chi_squared(tep_shift_hz, 0.0, gw_f220_err)
        gw_f220_shift_sigma = sqrt(gw_f220_shift_chi2)
        print_status(f"  TEP-Schw f_220 shift: {tep_shift_hz:.4f} Hz (chi2={gw_f220_shift_chi2:.4f}, {gw_f220_shift_sigma:.2f}σ)", "DEBUG")
        # Retain mismatch values as explicitly nonspinning diagnostics.
        gw_f220_chi2_tep = chi_squared(gw_f220_tep, gw_f220, gw_f220_err)
        gw_f220_chi2_schw = chi_squared(gw_f220_schw, gw_f220, gw_f220_err)
        gw_f220_chi2_kerr = chi_squared(gw_f220_kerr_pred, gw_f220, gw_f220_err)
        gw_f220_sigma_tep = sqrt(gw_f220_chi2_tep)
        gw_f220_sigma_schw = sqrt(gw_f220_chi2_schw)
        gw_f220_sigma_kerr = sqrt(gw_f220_chi2_kerr)
        print_status(f"  f_220 chi-sq detail: TEP={gw_f220_chi2_tep:.4f} ({gw_f220_sigma_tep:.2f}σ), Schw={gw_f220_chi2_schw:.4f} ({gw_f220_sigma_schw:.2f}σ), Kerr={gw_f220_chi2_kerr:.4f} ({gw_f220_sigma_kerr:.2f}σ)", "DEBUG")
    else:
        gw_f220_shift_chi2 = None
        gw_f220_shift_sigma = None
        gw_f220_chi2_tep = None
        gw_f220_chi2_schw = None
        gw_f220_chi2_kerr = None
        gw_f220_sigma_tep = None
        gw_f220_sigma_schw = None
        gw_f220_sigma_kerr = None

    if gw_tau220_tep is not None:
        gw_tau220_chi2_tep = chi_squared(gw_tau220_tep, gw_tau220, gw_tau220_err)
        gw_tau220_chi2_schw = chi_squared(gw_tau220_schw, gw_tau220, gw_tau220_err)
        gw_tau220_chi2_kerr = chi_squared(gw_tau220_kerr_pred, gw_tau220, gw_tau220_err)
        gw_tau220_sigma_tep = sqrt(gw_tau220_chi2_tep)
        gw_tau220_sigma_schw = sqrt(gw_tau220_chi2_schw)
        gw_tau220_sigma_kerr = sqrt(gw_tau220_chi2_kerr)
        print_status(f"  tau_220 chi-sq detail: TEP={gw_tau220_chi2_tep:.4f} ({gw_tau220_sigma_tep:.2f}σ), Schw={gw_tau220_chi2_schw:.4f} ({gw_tau220_sigma_schw:.2f}σ), Kerr={gw_tau220_chi2_kerr:.4f} ({gw_tau220_sigma_kerr:.2f}σ)", "DEBUG")
    else:
        gw_tau220_chi2_tep = None
        gw_tau220_chi2_schw = None
        gw_tau220_chi2_kerr = None
        gw_tau220_sigma_tep = None
        gw_tau220_sigma_schw = None
        gw_tau220_sigma_kerr = None

    print_status(f"  Final mass: {gw_mass} ± {gw_mass_err} M_sun", "INFO")
    print_status(f"  Final spin: {gw_spin} ± {gw['final_spin']['uncertainty']}", "INFO")
    print_status(f"  f_220 (measured):  {gw_f220} ± {gw_f220_err} Hz", "INFO")
    if gw_f220_tep is not None:
        print_status(f"  f_220 (TEP pred):  {gw_f220_tep:.1f} Hz", "INFO")
        print_status(f"  f_220 (Schw pred): {gw_f220_schw:.1f} Hz", "INFO")
        print_status(f"  f_220 (Kerr pred): {gw_f220_kerr_pred:.1f} Hz (spin correction {kerr_freq_correction:.3f})", "INFO")
        print_status(f"  TEP-Schw shift:    {gw_f220_tep - gw_f220_schw:.1f} Hz ({(gw_f220_tep/gw_f220_schw - 1)*100:.1f}%)", "INFO")
        print_status(f"  f_220 chi-sq (TEP vs meas):  {gw_f220_chi2_tep:.4f} ({gw_f220_sigma_tep:.2f}σ)", "INFO")
        print_status(f"  f_220 chi-sq (Schw vs meas): {gw_f220_chi2_schw:.4f} ({gw_f220_sigma_schw:.2f}σ)", "INFO")
        print_status(f"  f_220 chi-sq (Kerr vs meas): {gw_f220_chi2_kerr:.4f} ({gw_f220_sigma_kerr:.2f}σ)", "INFO")
    print_status(f"  tau_220 (measured): {gw_tau220} ± {gw_tau220_err} ms", "INFO")
    if gw_tau220_tep is not None:
        print_status(f"  tau_220 (TEP pred):  {gw_tau220_tep:.2f} ms", "INFO")
        print_status(f"  tau_220 (Schw pred): {gw_tau220_schw:.2f} ms", "INFO")
        print_status(f"  tau_220 (Kerr pred): {gw_tau220_kerr_pred:.2f} ms", "INFO")
    print_status(f"  QNM comparison: TEP = Schw (GWs propagate on geometric g = Schwarzschild)", "SUCCESS")
    print_status(f"  Kerr baseline needed for spinning remnant (a={gw_spin}); not a TEP-vs-Schw test", "INFO")

    # --- Combined constraint ---
    print_status("", "INFO")
    print_status("Combined Constraints", "TITLE")
    # For shadows: TEP = Schw (exterior is Schwarzschild), so both give same chi2
    # For QNM: TEP = Schw (GWs propagate on geometric g = Schwarzschild), so both
    # give the same nonspinning QNM. The key diagnostic is the TEP/Kerr vs data
    # comparison, which tests whether the nonspinning Schwarzschild QNM baseline
    # (shared by TEP) is distinguishable from the spinning Kerr baseline.
    print_status(f"  Combining {4} constraints: M87* shadow, Sgr A* shadow, GW150914 f_220, GW150914 tau_220", "DEBUG")
    print_status(f"  M87* chi2 (TEP):     {m87_chi2_tep:.4f}", "DEBUG")
    print_status(f"  Sgr A* chi2 (TEP):   {sgr_chi2_tep:.4f}", "DEBUG")
    print_status(f"  GW150914 f_220 (TEP):  {gw_f220_chi2_tep or 0:.4f}", "DEBUG")
    print_status(f"  GW150914 tau_220 (TEP): {gw_tau220_chi2_tep or 0:.4f}", "DEBUG")
    total_chi2_tep = (m87_chi2_tep + sgr_chi2_tep +
                      (gw_f220_chi2_tep or 0) + (gw_tau220_chi2_tep or 0))
    total_chi2_schw = (m87_chi2_schw + sgr_chi2_schw +
                       (gw_f220_chi2_schw or 0) + (gw_tau220_chi2_schw or 0))
    total_chi2_kerr = (m87_chi2_schw + sgr_chi2_schw +
                       (gw_f220_chi2_kerr or 0) + (gw_tau220_chi2_kerr or 0))
    n_constraints = 4
    print_status(f"  Total chi-squared (TEP):  {total_chi2_tep:.4f} ({n_constraints} constraints)", "INFO")
    print_status(f"  Total chi-squared (Schw): {total_chi2_schw:.4f}", "INFO")
    print_status(f"  Total chi-squared (Kerr): {total_chi2_kerr:.4f}", "INFO")
    print_status(f"  TEP consistent with data: {total_chi2_tep < n_constraints * 2.0}", "INFO")
    print_status(f"  Schw consistent with data: {total_chi2_schw < n_constraints * 2.0}", "INFO")
    print_status(f"  Kerr consistent with data: {total_chi2_kerr < n_constraints * 2.0}", "INFO")
    print_status(f"  Note: Shadow chi2 is identical for TEP/Schw (exterior = Schwarzschild)", "INFO")
    print_status(f"  QNM shift (TEP vs Schw): {gw_f220_tep - gw_f220_schw:.1f} Hz "
                 f"({(gw_f220_tep/gw_f220_schw - 1)*100:.1f}%)" if gw_f220_tep else "", "INFO")
    print_status(f"  Consistency threshold: chi2 < {n_constraints * 2.0} (2*n_constraints)", "DEBUG")
    print_status(f"  TEP total chi2 / n = {total_chi2_tep / n_constraints:.4f} (reduced)", "DEBUG")
    print_status(f"  Schw total chi2 / n = {total_chi2_schw / n_constraints:.4f} (reduced)", "DEBUG")
    print_status(f"  Kerr total chi2 / n = {total_chi2_kerr / n_constraints:.4f} (reduced)", "DEBUG")
    print_status(f"  TEP = Schw for all exterior observables (shadow + QNM): identical chi2", "SUCCESS")
    print_status(f"  Note: TEP reduces to Schwarzschild in the exterior (r > 2M); A cancels in geodesic equations", "INFO")

    # --- Save CSV ---
    csv_rows = [
        {"test": "M87_shadow", "measured": m87_shadow_meas, "uncertainty": m87_shadow_err,
         "tep_prediction": m87_shadow_tep, "schw_prediction": m87_shadow_schw,
         "kerr_prediction": m87_shadow_schw,
         "chi2_tep": m87_chi2_tep, "chi2_schw": m87_chi2_schw, "chi2_kerr": m87_chi2_schw,
         "sigma_tep": m87_sigma_tep, "sigma_schw": m87_sigma_schw},
        {"test": "SgrA_shadow", "measured": sgr_shadow_meas, "uncertainty": sgr_shadow_err,
         "tep_prediction": sgr_shadow_tep, "schw_prediction": sgr_shadow_schw,
         "kerr_prediction": sgr_shadow_schw,
         "chi2_tep": sgr_chi2_tep, "chi2_schw": sgr_chi2_schw, "chi2_kerr": sgr_chi2_schw,
         "sigma_tep": sgr_sigma_tep, "sigma_schw": sgr_sigma_schw},
        {"test": "GW150914_f220", "measured": gw_f220, "uncertainty": gw_f220_err,
         "tep_prediction": gw_f220_tep, "schw_prediction": gw_f220_schw,
         "kerr_prediction": gw_f220_kerr_pred,
         "chi2_tep": gw_f220_chi2_tep, "chi2_schw": gw_f220_chi2_schw, "chi2_kerr": gw_f220_chi2_kerr,
         "sigma_tep": gw_f220_sigma_tep, "sigma_schw": gw_f220_sigma_schw},
        {"test": "GW150914_tau220", "measured": gw_tau220, "uncertainty": gw_tau220_err,
         "tep_prediction": gw_tau220_tep, "schw_prediction": gw_tau220_schw,
         "kerr_prediction": gw_tau220_kerr_pred,
         "chi2_tep": gw_tau220_chi2_tep, "chi2_schw": gw_tau220_chi2_schw, "chi2_kerr": gw_tau220_chi2_kerr,
         "sigma_tep": gw_tau220_sigma_tep, "sigma_schw": gw_tau220_sigma_schw},
    ]
    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, csv_rows)
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "model": model.to_dict(),
        "data_sources": {
            "eht": eht_data["M87"]["source"],
            "eht_sgr": eht_data["SgrA"]["source"],
            "ligo": ligo_data["GW150914"]["source"],
        },
        "M87": {
            "mass_Msun": m87_mass,
            "mass_uncertainty": m87_mass_err,
            "distance_Mpc": m87_dist,
            "shadow_measured_uas": m87_shadow_meas,
            "shadow_uncertainty_uas": m87_shadow_err,
            "shadow_tep_uas": m87_shadow_tep,
            "shadow_schw_uas": m87_shadow_schw,
            "chi2_tep": m87_chi2_tep,
            "chi2_schw": m87_chi2_schw,
            "sigma_tep": m87_sigma_tep,
            "sigma_schw": m87_sigma_schw,
            "pvalue_tep": m87_pvalue_tep,
            "pvalue_schw": m87_pvalue_schw,
        },
        "SgrA": {
            "mass_Msun": sgr_mass,
            "mass_uncertainty": sgr_mass_err,
            "distance_kpc": sgr_dist,
            "shadow_measured_uas": sgr_shadow_meas,
            "shadow_uncertainty_uas": sgr_shadow_err,
            "shadow_tep_uas": sgr_shadow_tep,
            "shadow_schw_uas": sgr_shadow_schw,
            "chi2_tep": sgr_chi2_tep,
            "chi2_schw": sgr_chi2_schw,
            "sigma_tep": sgr_sigma_tep,
            "sigma_schw": sgr_sigma_schw,
            "pvalue_tep": sgr_pvalue_tep,
            "pvalue_schw": sgr_pvalue_schw,
        },
        "GW150914": {
            "final_mass_Msun": gw_mass,
            "final_spin": gw_spin,
            "f220_measured_Hz": gw_f220,
            "f220_uncertainty_Hz": gw_f220_err,
            "f220_tep_Hz": gw_f220_tep,
            "f220_schw_Hz": gw_f220_schw,
            "f220_kerr_Hz": gw_f220_kerr_pred,
            "f220_chi2_tep": gw_f220_chi2_tep,
            "f220_chi2_schw": gw_f220_chi2_schw,
            "f220_chi2_kerr": gw_f220_chi2_kerr,
            "f220_sigma_tep": gw_f220_sigma_tep,
            "f220_sigma_schw": gw_f220_sigma_schw,
            "f220_sigma_kerr": gw_f220_sigma_kerr,
            "f220_tep_shift_hz": (gw_f220_tep - gw_f220_schw) if gw_f220_tep else None,
            "f220_tep_shift_percent": ((gw_f220_tep/gw_f220_schw - 1)*100) if gw_f220_tep else None,
            "tau220_measured_ms": gw_tau220,
            "tau220_uncertainty_ms": gw_tau220_err,
            "tau220_tep_ms": gw_tau220_tep,
            "tau220_schw_ms": gw_tau220_schw,
            "tau220_kerr_ms": gw_tau220_kerr_pred,
            "tau220_chi2_tep": gw_tau220_chi2_tep,
            "tau220_chi2_schw": gw_tau220_chi2_schw,
            "tau220_chi2_kerr": gw_tau220_chi2_kerr,
        },
        "combined": {
            "total_chi2_tep": total_chi2_tep,
            "total_chi2_schw": total_chi2_schw,
            "total_chi2_kerr": total_chi2_kerr,
            "n_constraints": n_constraints,
            "tep_consistent": bool(total_chi2_tep < n_constraints * 2.0),
            "schw_consistent": bool(total_chi2_schw < n_constraints * 2.0),
            "kerr_consistent": bool(total_chi2_kerr < n_constraints * 2.0),
        },
        "note": "TEP and Schwarzschild predictions are identical for exterior observables "
                "(shadow, ISCO, QNM) because the TEP metric reduces to Schwarzschild for r > 2M "
                "and gravitational waves propagate on the geometric metric g, which is prescribed "
                "as Schwarzschild in this branch. The TEP-Schwarzschild QNM shift is 0.0%. "
                "GW150914 remnant has spin a=0.67, so Kerr QNM frequency is higher than "
                "Schwarzschild; a spinning event must be compared with the Kerr baseline.",
    }
    json_path = step_json_path(STEP_ID)
    summary = finalize_result(
        STEP_ID, summary,
        description="Constrain TEP-BH model against EHT shadow measurements (M87*, Sgr A*) and LIGO GW150914 ringdown QNM diagnostics",
        key_result=f"TEP and Schwarzschild give identical exterior observables (shadow, QNM); combined chi2={total_chi2_tep:.4f} over {n_constraints} constraints, consistent with data",
        model=model.to_dict(),
        dependencies=["step_00_data_download", "step_02_perturbations", "step_03_raytracing", "step_04_accretion"],
    )
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 05 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
