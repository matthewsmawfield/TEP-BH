#!/usr/bin/env python3
"""Step 11: EHT Polarization — TEP Exterior-Null Result.

TEP does not predict observable birefringence in the exterior of a
black hole. The conformal factor A(phi) = exp(beta_A * phi) approaches
unity in the weak-field exterior (phi -> 0), and the bounded
disformal function B(phi) = B0 * |phi|^n / (1 + |phi|^n) approaches zero as
phi -> 0. Consequently the disformal metric contribution and its coupling
to electromagnetic propagation vanish asymptotically in the exterior, and
photon polarization is unaffected by this model at leading order.
The disformal term only becomes dynamically relevant inside the horizon,
where phi grows large (negative) and B saturates at B0.

This step records the EHT polarimetric measurements for M87* and
Sgr A* as observational data, and honestly states that TEP does not
explain the observed low polarization fraction. The low polarization
(~4-7%) versus the standard GR+MHD prediction (~10-15%) is attributed
in the literature to Faraday rotation internal to the emission region,
not to a TEP birefringence signal.

Data sources:
  - EHT 2021, ApJ 910, L13 (M87* polarization paper VII)
  - EHT 2021, ApJ 910, L12 (M87* polarization paper VI)
  - EHT 2024, ApJ 964, L26 (M87* updated polarimetry)
  - Akiyama et al. 2015 (ray-tracing polarimetry methods)

Outputs:
  - data/processed/eht_m87_polarization.json
  - results/step_11_eht_polarization.json
  - results/step_11_eht_polarization.csv
  - logs/step_11_eht_polarization.log
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

import numpy as np

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
)

STEP_ID = "step_11_eht_polarization"

# M87* parameters (EHT 2019 Paper VI)
M87_M_MSUN = 6.5e9
M87_DIST_MPC = 16.8
M87_INCLINATION_DEG = 17.0  # degrees

# EHT polarization observations (Paper VII, 2021)
# The observed polarization fraction and EVPA pattern
EHT_POL_DATA = {
    "m87_star": {
        "observed_ring_polarization_fraction": {
            "value": 0.04,  # ~4% linear polarization
            "uncertainty": 0.01,
            "citation": "EHT 2021, ApJ 910, L13",
        },
        "evpa_swirl_angle_deg": {
            "value": 180.0,  # total EVPA rotation around ring
            "uncertainty": 30.0,
            "citation": "EHT 2021, ApJ 910, L13",
        },
        "magnetic_field_topology": "predominantly normal (poloidal)",
        "citations": [
            "EHT 2021, ApJ 910, L13 (Paper VII: Polarimetric Performance)",
            "EHT 2021, ApJ 910, L12 (Paper VI: Polarimetric Imaging)",
        ],
    },
    "sgr_a_star": {
        "observed_ring_polarization_fraction": {
            "value": 0.07,  # ~7% linear polarization
            "uncertainty": 0.02,
            "citation": "EHT 2024, ApJ 964, L26",
        },
        "evpa_swirl_angle_deg": {
            "value": 90.0,
            "uncertainty": 30.0,
            "citation": "EHT 2024, ApJ 964, L26",
        },
        "magnetic_field_topology": "mixed",
        "citations": [
            "EHT 2024, ApJ 964, L26 (Sgr A* Polarimetry)",
        ],
    },
}


def tep_exterior_polarization_null() -> dict:
    """Document the TEP exterior-null result for photon polarization.

    Under TEP, the conformal factor A(phi) = exp(beta_A * phi) and the
    bounded disformal function B(phi) = B0 * |phi|^n / (1 + |phi|^n)
    govern the coupling between the scalar field and the electromagnetic
    field. In the exterior of the black hole, the scalar field phi is
    small (weak-field regime), so:

      A(phi) ~ exp(beta_A * phi) ~ 1 + beta_A * phi  ->  1
      B(phi) ~ B0 * |phi|^n  ->  0

    The disformal coupling B * (dphi/dr)^2 therefore vanishes in the
    exterior, and the effective metric for electromagnetic propagation
    reduces to the conformal (Schwarzschild-like) metric with A ~ 1.
    There is no birefringence: both polarization modes propagate with
    the same phase velocity, and the net linear polarization is
    unaffected by TEP corrections in the exterior.

    The disformal term only becomes dynamically relevant inside the
    horizon, where phi grows large in magnitude and B -> B0. However, photons
    observed by the EHT escape from the exterior (r > r_h), so any
    interior disformal effect is not imprinted on the observed
    polarization signal.

    Returns
    -------
    dict
        Summary of the exterior-null result.
    """
    # In the exterior (weak field), phi -> 0, so:
    #   A -> exp(0) = 1
    #   B -> 0
    # The disformal birefringence signal is identically zero.
    A_exterior = 1.0
    B_exterior = 0.0
    delta_phase_total = 0.0  # no birefringence in exterior
    pol_suppression = 1.0  # no suppression from TEP

    return {
        "A_exterior": A_exterior,
        "B_exterior": B_exterior,
        "delta_phase_total_rad": delta_phase_total,
        "pol_suppression_factor": pol_suppression,
        "tep_predicts_birefringence": False,
        "explanation": (
            "In the exterior, A(phi) ~ 1 and B(phi) ~ 0, so the disformal "
            "coupling vanishes and photon polarization is unaffected. "
            "TEP does not predict observable birefringence in the exterior."
        ),
    }


def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 11: EHT Polarization — TEP Exterior-Null Result", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    # --- Save compiled EHT polarization data ---
    eht_path = PROCESSED_DIR / "eht_m87_polarization.json"
    write_json(eht_path, EHT_POL_DATA)
    print_status(f"EHT polarization data saved to {rel(eht_path)}", "SUCCESS")

    # --- M87 source parameters ---
    print_status("", "INFO")
    print_status("M87* Source Parameters", "TITLE")
    print_status(f"  Mass: M = {M87_M_MSUN:.2e} M_sun", "INFO")
    print_status(f"  Distance: D = {M87_DIST_MPC} Mpc", "INFO")
    print_status(f"  Inclination: i = {M87_INCLINATION_DEG} deg", "INFO")
    print_status(f"  [DEBUG] Schwarzschild shadow diameter ~ 2*sqrt(27)*M*km_per_Msun/(D*km_per_pc)*uas_per_rad", "DEBUG")

    # --- Polarization model details ---
    print_status("", "INFO")
    print_status("Polarization Model Channels", "TITLE")
    print_status("  Linear polarization: Stokes Q, U — primary EHT observable", "INFO")
    print_status("  Circular polarization: Stokes V — sub-percent, not resolved by EHT 2017", "INFO")
    print_status("  Birefringence: differential phase between polarization modes", "INFO")
    print_status("  [DEBUG] TEP birefringence source: disformal coupling B*(dphi)^2", "DEBUG")
    print_status("  [DEBUG] In exterior: phi -> 0 => B -> 0 => birefringence vanishes", "DEBUG")

    # --- Observed polarization summary ---
    print_status("", "INFO")
    print_status("Observed EHT Polarization", "TITLE")
    for source, data in EHT_POL_DATA.items():
        pol_frac = data["observed_ring_polarization_fraction"]["value"]
        pol_err = data["observed_ring_polarization_fraction"]["uncertainty"]
        evpa = data["evpa_swirl_angle_deg"]["value"]
        evpa_err = data["evpa_swirl_angle_deg"]["uncertainty"]
        print_status(f"  {source}:", "INFO")
        print_status(f"    Linear polarization: {pol_frac*100:.1f}% +/- {pol_err*100:.1f}%", "INFO")
        print_status(f"    EVPA swirl: {evpa:.0f} +/- {evpa_err:.0f} deg", "INFO")
        print_status(f"    B-field topology: {data['magnetic_field_topology']}", "INFO")
        print_status(
            f"    [DEBUG] {source}: pol_frac={pol_frac:.4f} +/- {pol_err:.4f}, "
            f"EVPA={evpa:.1f} deg, citation={data['observed_ring_polarization_fraction']['citation']}",
            "DEBUG",
        )

    # --- Standard GR+MHD prediction ---
    print_status("", "INFO")
    print_status("Standard GR+MHD Prediction", "TITLE")
    print_status("  Polarization fraction: 10-15% (for magnetically arrested disk)", "INFO")
    print_status("  EVPA: smooth azimuthal pattern", "INFO")
    print_status("  Issue: Observed polarization (4-7%) is LOWER than predicted", "INFO")
    print_status("  Standard explanation: Faraday rotation internal to disk", "INFO")

    # --- TEP exterior-null analysis ---
    print_status("", "INFO")
    print_status("TEP Exterior-Null Analysis", "TITLE")
    print_status("  Under TEP, the conformal factor A(phi) = exp(beta_A * phi)", "INFO")
    print_status("  and the bounded disformal function B(phi) govern the coupling", "INFO")
    print_status("  between the scalar field and the EM field.", "INFO")
    print_status("")
    print_status("  In the exterior (weak field, phi -> 0):", "INFO")
    print_status("    A(phi) ~ exp(beta_A * phi) -> 1", "INFO")
    print_status("    B(phi) ~ B0 * |phi|^n -> 0", "INFO")
    print_status("  The disformal coupling B * (dphi/dr)^2 vanishes in the exterior,", "INFO")
    print_status("  so there is no birefringence: both polarization modes propagate", "INFO")
    print_status("  with the same phase velocity. Photon polarization is unaffected.", "INFO")
    print_status("")
    print_status("  The disformal term only matters inside the horizon, where", "INFO")
    print_status("  phi grows large in magnitude and B -> B0. Photons observed by the EHT", "INFO")
    print_status("  escape from the exterior (r > r_h), so any interior disformal", "INFO")
    print_status("  effect is not imprinted on the observed polarization signal.", "INFO")
    print_status("")

    null_result = tep_exterior_polarization_null()

    print_status("", "INFO")
    print_status("TEP Birefringence Null-Result Computation", "TITLE")
    print_status(f"  A_exterior = {null_result['A_exterior']:.6f} (unity => no conformal phase shift)", "INFO")
    print_status(f"  B_exterior = {null_result['B_exterior']:.6f} (zero => no disformal birefringence)", "INFO")
    print_status(f"  delta_phase_total = {null_result['delta_phase_total_rad']:.6f} rad (identically zero)", "INFO")
    print_status(f"  pol_suppression_factor = {null_result['pol_suppression_factor']:.6f} (unity => no suppression)", "INFO")
    print_status(f"  tep_predicts_birefringence = {null_result['tep_predicts_birefringence']}", "SUCCESS")
    print_status(
        f"  [DEBUG] Null result: A->1, B->0 in exterior => disformal coupling "
        f"B*(dphi/dr)^2 -> 0, both polarization modes propagate identically",
        "DEBUG",
    )

    csv_rows = []
    for source, data in EHT_POL_DATA.items():
        pol_obs = data["observed_ring_polarization_fraction"]["value"]
        pol_pred_std = 0.10  # standard GR+MHD prediction
        print_status(f"  {source}:", "INFO")
        print_status(f"    Observed polarization: {pol_obs*100:.1f}%", "INFO")
        print_status(f"    Standard GR+MHD prediction: {pol_pred_std*100:.1f}%", "INFO")
        print_status(f"    TEP birefringence prediction: NONE (exterior-null)", "INFO")
        print_status(f"    TEP does not explain the low polarization fraction", "INFO")
        print_status(
            f"    [DEBUG] {source}: observed={pol_obs:.4f} vs predicted={pol_pred_std:.4f} "
            f"=> deficit={(pol_pred_std - pol_obs)*100:.1f} pp (Faraday rotation, not TEP)",
            "DEBUG",
        )

        csv_rows.append({
            "source": source,
            "mass_msun": M87_M_MSUN if "m87" in source else 4.0e6,
            "pol_observed": pol_obs,
            "pol_observed_uncertainty": data["observed_ring_polarization_fraction"]["uncertainty"],
            "pol_standard_predicted": 0.10,
            "tep_predicts_birefringence": False,
            "tep_pol_suppression_factor": null_result["pol_suppression_factor"],
            "tep_delta_phase_rad": null_result["delta_phase_total_rad"],
            "tep_explains_low_pol": False,
        })

    # --- Key finding ---
    print_status("", "INFO")
    print_status("Key Finding", "TITLE")
    print_status("  TEP does not predict observable birefringence in the exterior.", "INFO")
    print_status("  The conformal factor A ~ 1 and the disformal term B ~ 0 in", "INFO")
    print_status("  the exterior, so photon polarization is unaffected.", "INFO")
    print_status("  The observed low polarization fraction (~4-7%) is not", "INFO")
    print_status("  explained by TEP. The standard explanation (Faraday", "INFO")
    print_status("  rotation internal to the disk) remains the leading", "INFO")
    print_status("  interpretation in the literature.", "INFO")

    # --- Save CSV ---
    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, csv_rows)
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "data_sources": {
            "m87_star": EHT_POL_DATA["m87_star"]["citations"],
            "sgr_a_star": EHT_POL_DATA["sgr_a_star"]["citations"],
        },
        "m87_parameters": {
            "mass_msun": M87_M_MSUN,
            "distance_mpc": M87_DIST_MPC,
            "inclination_deg": M87_INCLINATION_DEG,
        },
        "tep_model": {
            "conformal_factor": "A(phi) = exp(beta_A * phi)",
            "disformal_function": "B(phi) = B0 * |phi|^n / (1 + |phi|^n)",
            "exterior_behavior": "A -> 1, B -> 0 as phi -> 0 in the weak-field exterior",
            "predicts_birefringence": False,
            "birefringence_regime": "interior only (r < r_h), not observable by EHT",
        },
        "tep_exterior_null_result": null_result,
        "results": {
            "m87_star": {
                "pol_observed": csv_rows[0]["pol_observed"],
                "tep_predicts_birefringence": csv_rows[0]["tep_predicts_birefringence"],
                "tep_explains_low_pol": csv_rows[0]["tep_explains_low_pol"],
            },
            "sgr_a_star": {
                "pol_observed": csv_rows[1]["pol_observed"],
                "tep_predicts_birefringence": csv_rows[1]["tep_predicts_birefringence"],
                "tep_explains_low_pol": csv_rows[1]["tep_explains_low_pol"],
            },
        },
        "conclusion": (
            "TEP does not predict observable birefringence in the exterior. "
            "The conformal factor A ~ 1 and the disformal term B ~ 0 in the "
            "exterior, so photon polarization is unaffected. The observed low "
            "polarization fraction is not explained by TEP."
        ),
    }
    json_path = step_json_path(STEP_ID)
    summary = finalize_result(
        STEP_ID, summary,
        description=(
            "Record EHT polarimetric measurements for M87* and Sgr A* and "
            "document the TEP exterior-null result: the conformal factor A~1 "
            "and disformal function B~0 in the exterior, so TEP predicts no "
            "observable birefringence."
        ),
        key_result=(
            "TEP does not predict observable birefringence in the exterior; "
            "the observed low polarization fraction (~4-7%) is not explained "
            "by TEP and is attributed to Faraday rotation in standard GR+MHD."
        ),
        dependencies=["step_00_data_download"],
    )
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 11 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
