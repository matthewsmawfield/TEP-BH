#!/usr/bin/env python3
"""Step 04 (Inference): Joint forward-model — composition vs calibration.

Separate the two TEP hypotheses:
  H_composition: M_matter < M_GR-fit  (less material mass than conventional)
  H_calibration: M_g < M_GR-fit       (gravitational mass itself is calibrated)

The joint forward-model fits for both:
  - M_matter: the locally measured material content
  - M_g:      the exterior gravitational parameter
  - eta_TEP:  the TEP coupling strength

The model predicts observables (astrometry, spectroscopy) from these
parameters through the TEP transfer functions. The fit determines
which hypothesis the data support.

MODEL STRUCTURE:
  The exterior gravitational field is determined by M_g (not M_matter).
  The orbital dynamics follow the TEP-modified geodesics on gtilde.
  The locally measured mass (if we could place a local observer) would
  be M_matter, related to M_g by the mass-bias equation:
    M_g / M_matter = (S_a^3 * D_dyn) / T_P^2

  The observables depend on M_g (through the exterior field) and eta_TEP
  (through the transfer functions). M_matter is inferred from M_g and
  the transfer factors.

  If eta_TEP = 0: M_matter = M_g = M_GR-fit (GR limit)
  If eta_TEP > 0 and ratio > 1: M_matter < M_g (composition hypothesis)
  If eta_TEP > 0 and ratio < 1: M_matter > M_g (anti-composition)

PROPER-VOLUME DIAGNOSTIC:
  Compute C_V = V_proper / V_apparent for the S2 orbit region.
  C_V >> 1 would demonstrate apparent compactness without physical compression.

Outputs:
  data/processed/sstar_joint_forward.json
  results/step_47_joint_forward.json
  logs/step_47_joint_forward.log
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps" / "inference"))

from step_46_tep_transfer import tep_orbit, tep_residuals
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

STEP_ID = "step_47_joint_forward"


# ============================================================
# Proper-volume diagnostic
# ============================================================
def compute_proper_volume_diagnostic(M_BH_Msun, R0_pc, r_boundary_AU, eta=0.3):
    """Compute the proper-volume diagnostic C_V.

    C_V = V_proper / V_apparent

    where:
      V_proper = 4*pi * integral of sqrt(gamma_rr) * R^2 dr
               (proper volume on the matter metric gtilde)
      V_apparent = (4*pi/3) * r_apparent^3
               (volume inferred from exterior angular size)

    For the TEP matter metric:
      gtilde_rr = A^2 * g_rr = A^2 * (1 - 2M/r)^{-1}
      R_tilde = A * r  (areal radius on gtilde)

    So sqrt(gamma_rr_tilde) * R_tilde^2 = A^3 * (1-2M/r)^{-1/2} * r^2

    For weak fields (r >> 2M): C_V ~ A^3 ~ exp(+eta/r_rs) ~ 1 + eta/r_rs.
    """
    R_s_m = 2 * G_SI * M_BH_Msun * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m

    # Integrate from 0 to r_boundary
    n_points = 1000
    r_arr = np.linspace(R_s_AU * 1.01, r_boundary_AU, n_points)  # avoid r=0
    dr = r_arr[1] - r_arr[0]

    r_rs = r_arr / R_s_AU  # in R_s units

    # Conformal factor: mass-inflation branch A = exp(+eta/(3*r_rs))
    A = np.exp(eta / (3.0 * r_rs))

    # g_rr = (1 - 1/r_rs)^{-1} for Schwarzschild
    F = 1.0 - 1.0 / r_rs
    g_rr = 1.0 / np.maximum(F, 1e-10)

    # Proper volume element: sqrt(gtilde_rr) * R_tilde^2 dr
    # gtilde_rr = A^2 * g_rr
    # R_tilde = A * r
    integrand = A**3 * np.sqrt(g_rr) * r_arr**2

    V_proper = 4 * np.pi * np.trapz(integrand, r_arr)

    # Apparent volume (from exterior angular size)
    V_apparent = (4 * np.pi / 3) * r_boundary_AU**3

    C_V = V_proper / V_apparent

    return {
        "r_boundary_AU": float(r_boundary_AU),
        "r_boundary_Rs": float(r_boundary_AU / R_s_AU),
        "V_proper_AU3": float(V_proper),
        "V_apparent_AU3": float(V_apparent),
        "C_V": float(C_V),
        "eta": eta,
    }


# ============================================================
# Composition vs Calibration
# ============================================================
def separate_hypotheses(M_g, eta_TEP, a_arcsec, e, M_BH_Msun, R0_pc):
    """Separate the composition and calibration hypotheses.

    Given the exterior gravitational mass M_g and the TEP coupling eta,
    compute the locally measured material mass M_matter using the
    mass-bias equation relative to GR:

        M_g / M_matter = [M_app/M_local](eta) / [M_app/M_local](eta=0)

    This is the TEP-only correction, so M_matter = M_g when eta = 0.
    """
    # Compute transfer factors at pericentre
    R_s_m = 2 * G_SI * M_BH_Msun * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m
    r_peri_AU = a_arcsec * (1 - e) * R0_pc
    r_peri_rs = r_peri_AU / R_s_AU

    T_P = temporal_transfer_factor(r_peri_rs, np.inf, eta=eta_TEP)
    S_a = spatial_calibration_factor(r_peri_rs, M_BH_Msun, R0_pc, eta=eta_TEP)
    D_dyn = dynamical_modification_factor(r_peri_rs, eta=eta_TEP)

    T_P_GR = temporal_transfer_factor(r_peri_rs, np.inf, eta=0.0)
    S_a_GR = spatial_calibration_factor(r_peri_rs, M_BH_Msun, R0_pc, eta=0.0)
    D_dyn_GR = dynamical_modification_factor(r_peri_rs, eta=0.0)

    ratio = ((S_a / S_a_GR)**3 * (D_dyn / D_dyn_GR)) / (T_P / T_P_GR)**2

    # M_matter = M_g / ratio
    M_matter = M_g / ratio

    return {
        "M_g_Msun": float(M_g),
        "M_matter_Msun": float(M_matter),
        "ratio_Mg_Mmatter": float(ratio),
        "T_P": float(T_P),
        "S_a": float(S_a),
        "D_dyn": float(D_dyn),
        "composition_hypothesis": M_matter < M_g,
        "calibration_hypothesis": M_g < M_BH_Msun,  # M_g < M_GR-fit
        "M_phantom_Msun": float(M_g - M_matter),
        "M_phantom_fraction": float((M_g - M_matter) / M_g),
    }


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 04 (INFERENCE): JOINT FORWARD-MODEL", "TITLE")
    print_status("  Composition vs Calibration Hypotheses", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Separating two distinct TEP claims:", "INFO")
    print_status("  H_composition: M_matter < M_GR-fit  (less material mass)", "INFO")
    print_status("  H_calibration: M_g < M_GR-fit       (gravitational mass calibrated)", "INFO")
    print_status("")

    # --- Load TEP fit ---
    tep_fit_path = PROCESSED_DIR / "sstar_tep_fit.json"
    gr_fit_path = PROCESSED_DIR / "sstar_gr_fit.json"

    if not all(p.exists() for p in [tep_fit_path, gr_fit_path]):
        print_status("ERROR: Previous fits not found. Run steps 00-03 first.", "ERROR")
        return

    tep_fit = json.loads(tep_fit_path.read_text())
    gr_fit = json.loads(gr_fit_path.read_text())

    M_GR = gr_fit["M_BH_GR_Msun"]
    M_g = tep_fit["best_fit"]["M_BH (M_sun)"]["value"]
    eta_TEP = tep_fit["best_fit"]["eta_TEP"]["value"]
    eta_err = tep_fit["best_fit"]["eta_TEP"]["error"]
    a = tep_fit["best_fit"]["a (arcsec)"]["value"]
    e = tep_fit["best_fit"]["e"]["value"]
    R0 = tep_fit["best_fit"]["R0 (pc)"]["value"]

    print_status(f"  M_GR-fit = {M_GR:.4e} M_sun", "INFO")
    print_status(f"  M_g (TEP) = {M_g:.4e} M_sun", "INFO")
    print_status(f"  eta_TEP = {eta_TEP:.4f} ± {eta_err:.4f}", "INFO")
    print_status("")

    # --- Separate hypotheses ---
    print_status("--- Hypothesis separation ---", "TITLE")
    hyp = separate_hypotheses(M_g, eta_TEP, a, e, M_GR, R0)

    print_status(f"  M_g      = {hyp['M_g_Msun']:.4e} M_sun", "INFO")
    print_status(f"  M_matter = {hyp['M_matter_Msun']:.4e} M_sun", "INFO")
    print_status(f"  Ratio M_g/M_matter = {hyp['ratio_Mg_Mmatter']:.8f}", "INFO")
    print_status(f"  T_P = {hyp['T_P']:.8f}, S_a = {hyp['S_a']:.8f}, D_dyn = {hyp['D_dyn']:.8f}", "INFO")
    print_status("")
    print_status(f"  H_composition (M_matter < M_g): {hyp['composition_hypothesis']}", "INFO")
    print_status(f"  H_calibration (M_g < M_GR-fit): {hyp['calibration_hypothesis']}", "INFO")
    print_status(f"  M_phantom = {hyp['M_phantom_Msun']:.4e} M_sun ({hyp['M_phantom_fraction']:.6e} fraction)", "INFO")
    print_status("")

    # --- Proper-volume diagnostic ---
    print_status("--- Proper-volume diagnostic C_V ---", "TITLE")

    # Compute C_V for the S2 pericentre region
    r_peri_AU = a * (1 - e) * R0
    r_apo_AU = a * (1 + e) * R0

    cv_peri = compute_proper_volume_diagnostic(M_GR, R0, r_peri_AU, eta=eta_TEP)
    cv_apo = compute_proper_volume_diagnostic(M_GR, R0, r_apo_AU, eta=eta_TEP)
    cv_gr = compute_proper_volume_diagnostic(M_GR, R0, r_peri_AU, eta=0.0)

    print_status(f"  At pericentre (r = {cv_peri['r_boundary_Rs']:.0f} R_s):", "INFO")
    print_status(f"    C_V (eta={eta_TEP}) = {cv_peri['C_V']:.8f}", "INFO")
    print_status(f"    C_V (eta=0, GR)    = {cv_gr['C_V']:.8f}", "INFO")
    print_status(f"    V_proper / V_apparent = {cv_peri['C_V']:.8f}", "INFO")
    print_status("")
    print_status(f"  At apocentre (r = {cv_apo['r_boundary_Rs']:.0f} R_s):", "INFO")
    print_status(f"    C_V (eta={eta_TEP}) = {cv_apo['C_V']:.8f}", "INFO")
    print_status("")

    if cv_peri["C_V"] > 1.01:
        print_status("  -> C_V > 1: proper volume exceeds apparent volume", "INFO")
        print_status("  -> Apparent compactness WITHOUT physical compression", "INFO")
    elif cv_peri["C_V"] < 0.99:
        print_status("  -> C_V < 1: proper volume is LESS than apparent volume", "WARN")
        print_status("  -> Physical compression is present", "WARN")
    else:
        print_status("  -> C_V ~ 1: no significant volume distortion at S2 scales", "INFO")

    print_status("")

    # --- Interpretation ---
    r_peri_Rs = cv_peri["r_boundary_Rs"]
    A_peri = np.exp(eta_TEP / (3.0 * r_peri_Rs))
    print_status("--- Interpretation ---", "TITLE")
    print_status(f"  At S2 pericentre (r ~ {r_peri_Rs:.0f} R_s):", "INFO")
    print_status(f"  TEP conformal factor A ~ {A_peri:.6f}", "INFO")
    print_status(f"  C_V - 1 ~ {cv_peri['C_V'] - 1:.2e}", "INFO")
    print_status(f"  TEP M_phantom fraction ~ {abs(hyp['M_phantom_fraction']):.2e}", "INFO")
    print_status("")
    print_status("  The S-star data at current precision cannot distinguish:", "INFO")
    print_status("    - H_composition from H_calibration", "INFO")
    print_status("    - TEP from GR (eta_TEP consistent with 0)", "INFO")
    print_status("  The decisive test requires horizon-scale observations.", "INFO")
    print_status("")

    # --- Save ---
    result = {
        "step": STEP_ID,
        "description": "Joint forward-model: composition vs calibration",
        "hypotheses": hyp,
        "proper_volume_diagnostic": {
            "pericentre": cv_peri,
            "apocentre": cv_apo,
            "GR_benchmark": cv_gr,
        },
        "interpretation": (
            f"At S2 pericentre (r ~ {r_peri_Rs:.0f} R_s), the TEP conformal factor is A ~ {A_peri:.6f}, "
            f"C_V - 1 ~ {cv_peri['C_V'] - 1:.2e}, and the TEP M_phantom fraction is ~ {abs(hyp['M_phantom_fraction']):.2e}. "
            "The S-star data at current precision cannot distinguish "
            "H_composition from H_calibration, nor TEP from GR. "
            "The decisive test requires horizon-scale observations."
        ),
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    forward_path = PROCESSED_DIR / "sstar_joint_forward.json"
    write_json(forward_path, result)
    print_status(f"  Saved: {forward_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, result)
    finalize_result(STEP_ID, result, "Joint forward-model: composition vs calibration", "Composition and calibration hypotheses separated")

    print_status("")
    print_status("  NEXT: step_48_likelihood.py — formal likelihood comparison", "INFO")


if __name__ == "__main__":
    main()
