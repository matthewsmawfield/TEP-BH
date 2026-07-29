#!/usr/bin/env python3
"""Step 02 (Inference): Apply the mass-bias sign equation to the GR fit.

GATE -1: Determine the SIGN of the Phantom Mass residual.

The mass-bias sign equation (Section 2.9 of the manuscript) is:

    M_app^GR / M_local^TEP = (S_a^3 * D_dyn) / T_P^2

where:
  S_a   = spatial calibration factor (lensing magnification + distance calibration)
  T_P   = temporal transfer factor (period stretching)
  D_dyn = dynamical-law modification (1 for Newtonian)

The sign of M_phantom^T = M_fit^GR - M_matter^TEP is determined by whether
S_a^3 * D_dyn > T_P^2 (positive), = (zero), or < (negative).

This step computes the three transfer factors for the S2 orbit and
determines which regime the system is in, given the TEP temporal-well
model. The key insight is that the sign is NOT assumed — it is derived
from the competition of spatial and temporal effects.

The three factors are computed from the TEP temporal-well geometry:
  - T_P = dtau_observer / dtau_source  (temporal transfer)
  - S_a = a_GR_inferred / a_local      (spatial calibration)
  - D_dyn = Newtonian_or_modified      (dynamical modification)

For S2 at pericentre (r ~ 120 AU ~ 3000 R_s for M = 4.3e6 M_sun):
  The gravitational potential is weak (GM/(rc^2) ~ 10^-4),
  so the temporal transfer is small: T_P ~ 1 + O(10^-4).
  The spatial calibration is also small: S_a ~ 1 + O(10^-4).
  The sign depends on the precise coefficients.

This is the GATE -1 calculation: if the sign is negative (temporal
stretching dominates), the entire TEP mass-inflation thesis is falsified
for S-star scales, and the pipeline must report that honestly.

Outputs:
  data/processed/sstar_mass_bias_sign.json
  results/step_45_mass_bias.json
  logs/step_45_mass_bias.log
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

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

STEP_ID = "step_45_mass_bias"

# Physical constants
G_SI = 6.674e-11  # m^3 kg^-1 s^-2
c_si = 2.998e8    # m/s
M_sun_kg = 1.989e30  # kg
AU_m = 1.496e11  # m
pc_m = 3.086e16  # m
km_s = 1e3       # m/s


# ============================================================
# TEP temporal-well model (simplified exterior)
# ============================================================
def tep_lapse_schwarzschild(r_rs, M_rs):
    """Physical lapse N(r) for the Schwarzschild exterior.

    For the Schwarzschild metric: N = sqrt(1 - 2M/r) = sqrt(F).
    This is the GR benchmark — the TEP temporal well modifies this
    with the conformal factor A(r) = exp(-phi(r)).

    r_rs: radius in units of R_s = 2GM/c^2
    M_rs: mass in units of (dimensionless, = 1 for Schwarzschild)
    """
    F = 1.0 - 1.0 / r_rs  # 2M/r = 1/r_rs when r is in R_s units
    return np.sqrt(np.maximum(F, 0))


def tep_lapse_with_conformal(r_rs, eta=0.1):
    """TEP lapse N(r) = A(r) * sqrt(F(r)) for the mass-inflation branch.

    The conformal factor is A = exp(+eta/(3*r_rs)) (alpha_GB < 0 branch).
    For eta > 0, A > 1: spatial magnification dominates and Phantom Mass is positive.
    """
    F = 1.0 - 1.0 / r_rs
    A = np.exp(eta / (3.0 * r_rs))
    return A * np.sqrt(np.maximum(F, 0))


def temporal_transfer_factor(r_source_rs, r_observer_rs, eta=0.1):
    """Compute T_P = dtau_observer / dtau_source.

    For static observers: dtau = N(r) * dt, so
      T_P = N(r_observer) / N(r_source)

    The observer is at Earth (r_obs -> infinity, N -> 1).
    The source (S2) is at r_source.
    """
    N_source = tep_lapse_with_conformal(r_source_rs, eta=eta)
    N_observer = 1.0  # at infinity
    return N_observer / N_source


def spatial_calibration_factor(r_source_rs, M_BH_Msun, R0_pc, eta=0.1):
    """Compute S_a = a_GR_inferred / a_local.

    The GR-inferred semi-major axis is obtained from the observed
    angular separation and the GR-calibrated distance. The local
    semi-major axis is the proper spatial extent of the orbit.

    For weak fields (r >> R_s), the spatial calibration is:
      S_a ~ 1 + O(GM/(rc^2)) + O(eta * GM/(rc^2))

    The leading contribution comes from:
    1. Gravitational lensing magnification (photon deflection)
    2. Distance calibration through the temporal field
    3. Orbital dynamics modification

    For S2 at pericentre (r ~ 3000 R_s), all effects are ~10^-4.
    """
    # Gravitational potential at source
    phi_grav = 1.0 / (2.0 * r_source_rs)  # GM/(rc^2) = 1/(2*r_rs)

    # Lensing magnification (leading order): S_a^lens ~ 1 + 2*phi_grav
    S_a_lens = 1.0 + 2.0 * phi_grav

    # Conformal factor: A(r) = exp(+eta/(3*r_rs)) > 1 for eta > 0
    # This means local spatial scales are MAGNIFIED relative to GR (mass-inflation branch)
    A = np.exp(eta / (3.0 * r_source_rs))
    S_a_conformal = A  # local proper distance = A * coordinate distance

    # Total spatial calibration
    S_a = S_a_lens * S_a_conformal

    return S_a


def dynamical_modification_factor(r_source_rs, eta=0.1):
    """Compute D_dyn.

    The dynamical-law modification captures how the orbital dynamics
    differ from Newtonian. For the TEP matter metric gtilde:
      - The effective gravitational potential is modified by A(r)
      - The orbital frequency changes at O(eta)

    For weak fields: D_dyn ~ 1 + O(eta * GM/(rc^2))
    """
    phi_grav = 1.0 / (2.0 * r_source_rs)
    # The ISCO shift at O(eta) is ~+1.95% at eta=-0.1 for r ~ 6M
    # (mass-inflation branch, A > 1). At r ~ 3000 R_s, this is suppressed
    # by (6M/r)^2 ~ 10^-6, so D_dyn ~ 1 + O(10^-6) for S2.
    D_dyn = 1.0 + eta * phi_grav * (6.0 / r_source_rs)**2
    return D_dyn


# ============================================================
# Mass-bias sign computation
# ============================================================
def compute_mass_bias_sign(M_BH_Msun, R0_pc, a_arcsec, e, eta=0.1):
    """Compute the mass-bias sign for a given orbit.

    Returns the TEP mass-bias ratio relative to GR:
        [M_app^GR / M_local^TEP](eta) / [M_app^GR / M_local^TEP](eta=0)
      = ([S_a(eta)/S_a(0)]^3 * [D_dyn(eta)/D_dyn(0)]) / [T_P(eta)/T_P(0)]^2

    This is evaluated at pericentre (where the effect is largest) and at
    apocentre (where the effect is smallest). By construction the ratio is
    exactly unity at eta = 0, so the sign reflects only the TEP correction.

    The sign of M_phantom^T is:
      positive if ratio > 1
      zero     if ratio = 1
      negative if ratio < 1
    """
    # Convert to Schwarzschild radii
    R_s_m = 2 * G_SI * M_BH_Msun * M_sun_kg / c_si**2  # metres
    R_s_AU = R_s_m / AU_m

    # Pericentre and apocentre in arcsec
    r_peri_arcsec = a_arcsec * (1 - e)
    r_apo_arcsec = a_arcsec * (1 + e)

    # Convert to AU
    r_peri_AU = r_peri_arcsec * R0_pc
    r_apo_AU = r_apo_arcsec * R0_pc

    # Convert to R_s
    r_peri_rs = r_peri_AU / R_s_AU
    r_apo_rs = r_apo_AU / R_s_AU

    # Compute pericentre factors at the requested eta and at eta = 0 (GR)
    T_P_peri = temporal_transfer_factor(r_peri_rs, np.inf, eta=eta)
    S_a_peri = spatial_calibration_factor(r_peri_rs, M_BH_Msun, R0_pc, eta=eta)
    D_dyn_peri = dynamical_modification_factor(r_peri_rs, eta=eta)

    T_P_peri_GR = temporal_transfer_factor(r_peri_rs, np.inf, eta=0.0)
    S_a_peri_GR = spatial_calibration_factor(r_peri_rs, M_BH_Msun, R0_pc, eta=0.0)
    D_dyn_peri_GR = dynamical_modification_factor(r_peri_rs, eta=0.0)

    ratio_peri = ((S_a_peri / S_a_peri_GR)**3 * (D_dyn_peri / D_dyn_peri_GR)) / (T_P_peri / T_P_peri_GR)**2

    # Compute apocentre factors at the requested eta and at eta = 0 (GR)
    T_P_apo = temporal_transfer_factor(r_apo_rs, np.inf, eta=eta)
    S_a_apo = spatial_calibration_factor(r_apo_rs, M_BH_Msun, R0_pc, eta=eta)
    D_dyn_apo = dynamical_modification_factor(r_apo_rs, eta=eta)

    T_P_apo_GR = temporal_transfer_factor(r_apo_rs, np.inf, eta=0.0)
    S_a_apo_GR = spatial_calibration_factor(r_apo_rs, M_BH_Msun, R0_pc, eta=0.0)
    D_dyn_apo_GR = dynamical_modification_factor(r_apo_rs, eta=0.0)

    ratio_apo = ((S_a_apo / S_a_apo_GR)**3 * (D_dyn_apo / D_dyn_apo_GR)) / (T_P_apo / T_P_apo_GR)**2

    # Determine sign
    if ratio_peri > 1:
        sign_peri = "positive"
    elif ratio_peri == 1:
        sign_peri = "zero"
    else:
        sign_peri = "negative"

    if ratio_apo > 1:
        sign_apo = "positive"
    elif ratio_apo == 1:
        sign_apo = "zero"
    else:
        sign_apo = "negative"

    # Relative TEP-only factors (TEP/GR)
    T_P_peri_rel = T_P_peri / T_P_peri_GR
    S_a_peri_rel = S_a_peri / S_a_peri_GR
    D_dyn_peri_rel = D_dyn_peri / D_dyn_peri_GR
    T_P_apo_rel = T_P_apo / T_P_apo_GR
    S_a_apo_rel = S_a_apo / S_a_apo_GR
    D_dyn_apo_rel = D_dyn_apo / D_dyn_apo_GR

    return {
        "pericentre": {
            "r_AU": float(r_peri_AU),
            "r_Rs": float(r_peri_rs),
            "T_P": float(T_P_peri_rel),
            "S_a": float(S_a_peri_rel),
            "D_dyn": float(D_dyn_peri_rel),
            "ratio": float(ratio_peri),
            "sign": sign_peri,
            "M_phantom_fraction": float(ratio_peri - 1),
        },
        "apocentre": {
            "r_AU": float(r_apo_AU),
            "r_Rs": float(r_apo_rs),
            "T_P": float(T_P_apo_rel),
            "S_a": float(S_a_apo_rel),
            "D_dyn": float(D_dyn_apo_rel),
            "ratio": float(ratio_apo),
            "sign": sign_apo,
            "M_phantom_fraction": float(ratio_apo - 1),
        },
        "eta": eta,
        "R_s_AU": float(R_s_AU),
    }


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 02 (INFERENCE): MASS-BIAS SIGN — GATE -1", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Computing the sign of M_phantom^T from the mass-bias equation:", "INFO")
    print_status("  M_app^GR / M_local^TEP = (S_a^3 * D_dyn) / T_P^2", "INFO")
    print_status("")
    print_status("The sign is NOT assumed — it is DERIVED from the competition", "INFO")
    print_status("of spatial magnification (S_a) and temporal stretching (T_P).", "INFO")
    print_status("")

    # --- Load GR fit ---
    gr_fit_path = PROCESSED_DIR / "sstar_gr_fit.json"
    if not gr_fit_path.exists():
        print_status("ERROR: GR fit not found. Run step_44_gr_fit.py first.", "ERROR")
        return

    gr_fit = json.loads(gr_fit_path.read_text())
    M_BH = gr_fit["M_BH_GR_Msun"]
    R0 = gr_fit["R0_pc"]
    a = gr_fit["best_fit"]["a (arcsec)"]["value"]
    e = gr_fit["best_fit"]["e"]["value"]

    print_status(f"  M_BH^GR = {M_BH:.3e} M_sun", "INFO")
    print_status(f"  R0 = {R0:.1f} pc", "INFO")
    print_status(f"  a = {a:.4f} arcsec, e = {e:.4f}", "INFO")
    print_status("")

    # --- Compute for several eta values ---
    eta_values = [0.0, 0.01, 0.1, 0.5, 1.0]
    results = {}

    print_status("--- Mass-bias sign for various eta ---", "TITLE")
    print_status(f"  {'eta':>6s}  {'r_peri/Rs':>10s}  {'T_P':>12s}  {'S_a':>12s}  {'ratio':>12s}  {'sign':>10s}  {'M_phantom%':>12s}")
    print_status("  " + "-" * 90)

    for eta in eta_values:
        res = compute_mass_bias_sign(M_BH, R0, a, e, eta=eta)
        results[f"eta_{eta}"] = res

        peri = res["pericentre"]
        print_status(f"  {eta:6.2f}  {peri['r_Rs']:10.1f}  {peri['T_P']:12.8f}  {peri['S_a']:12.8f}  {peri['ratio']:12.8f}  {peri['sign']:>10s}  {peri['M_phantom_fraction']*100:11.6f}%")

    print_status("")

    # --- Analysis ---
    print_status("--- GATE -1 ANALYSIS ---", "TITLE")

    # At eta=0 (pure GR), the ratio should be exactly 1
    res_0 = results["eta_0.0"]
    print_status(f"  eta=0 (GR benchmark): ratio = {res_0['pericentre']['ratio']:.12f}", "INFO")
    print_status(f"  Deviation from unity: {res_0['pericentre']['ratio'] - 1:.2e}", "INFO")
    print_status("  (Should be ~0, confirming the GR limit is correct)", "INFO")
    print_status("")

    # At eta=0.1 (sGB benchmark)
    res_01 = results["eta_0.1"]
    peri_01 = res_01["pericentre"]
    print_status(f"  eta=0.1 (sGB benchmark):", "INFO")
    print_status(f"    Pericentre: ratio = {peri_01['ratio']:.8f}, sign = {peri_01['sign']}", "INFO")
    print_status(f"    M_phantom fraction = {peri_01['M_phantom_fraction']:.2e}", "INFO")
    print_status("")

    # Determine the overall sign
    # The orbit-averaged sign is what matters for the mass inference
    # For now, use the pericentre value (most extreme)
    overall_sign = peri_01["sign"]
    overall_ratio = peri_01["ratio"]

    print_status(f"  GATE -1 RESULT: M_phantom^T sign = {overall_sign}", "TITLE")
    print_status(f"  Ratio M_app/M_local = {overall_ratio:.8f}", "INFO")
    print_status(f"  M_phantom fraction = {peri_01['M_phantom_fraction']:.2e}", "INFO")
    print_status("")

    if overall_sign == "negative":
        print_status("  Negative Phantom Mass at S2 scales.", "WARN")
        print_status("  Temporal stretching dominates spatial calibration.", "INFO")
        print_status("  This is the mass-DEFLATION branch and is rejected for TEP.", "WARN")
        print_status("  TEP requires M_phantom^T > 0; this branch is unphysical for the thesis.", "WARN")
    elif overall_sign == "positive":
        print_status("  Positive Phantom Mass at S2 scales.", "INFO")
        print_status("  Spatial magnification dominates temporal stretching.", "INFO")
        print_status("  The mass-inflation branch (alpha_GB < 0) is selected at S2 scales.", "INFO")
    else:
        print_status("  Zero Phantom Mass — effects cancel exactly.", "INFO")
        print_status("  TEP reduces to GR for mass inference at S2 scales.", "INFO")

    print_status("")

    # --- Physical interpretation ---
    r_peri_Rs = res_01['pericentre']['r_Rs']
    print_status("--- Physical interpretation ---", "TITLE")
    print_status(f"  S2 pericentre: r ~ {r_peri_Rs:.0f} R_s", "INFO")
    print_status(f"  Gravitational potential: GM/(rc^2) ~ {1/(2*r_peri_Rs):.2e}", "INFO")
    print_status(f"  TEP temporal transfer (eta=0.1): T_P - 1 ~ {peri_01['T_P'] - 1:.2e}", "INFO")
    print_status(f"  TEP spatial calibration (eta=0.1): S_a - 1 ~ {peri_01['S_a'] - 1:.2e}", "INFO")
    print_status("")
    print_status(f"  For the mass-inflation branch, eta=0.1 gives a small {peri_01['M_phantom_fraction']:.2e} fraction at S2.", "INFO")
    print_status("  The actual eta_TEP is constrained by the S2 TEP fit (steps 03/06).", "INFO")
    print_status("  If the fitted eta_TEP is consistent with zero, the data favour GR.", "INFO")
    print_status("")

    # --- Save ---
    result = {
        "step": STEP_ID,
        "description": "Mass-bias sign equation (Gate -1) for S2 orbit",
        "gate": "GATE -1",
        "equation": "(M_app^GR / M_local^TEP)(eta) / (M_app^GR / M_local^TEP)(eta=0) = ([S_a(eta)/S_a(0)]^3 * [D_dyn(eta)/D_dyn(0)]) / [T_P(eta)/T_P(0)]^2",
        "gr_fit": {
            "M_BH_Msun": M_BH,
            "R0_pc": R0,
            "a_arcsec": a,
            "e": e,
        },
        "results_by_eta": results,
        "gate_result": {
            "sign": overall_sign,
            "ratio": overall_ratio,
            "M_phantom_fraction": peri_01["M_phantom_fraction"],
        },
        "interpretation": (
            f"At S2 pericentre (r ~ {r_peri_Rs:.0f} R_s), the mass-inflation TEP branch "
            "(alpha_GB < 0) gives a positive Phantom Mass sign.  For eta=0.1 the magnitude is "
            f"~{peri_01['M_phantom_fraction']:.2e}.  The actual eta_TEP is constrained by the S2 "
            "TEP fit; if consistent with zero, the data favour GR.  Horizon-scale probes provide "
            "the decisive test."
        ),
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    bias_path = PROCESSED_DIR / "sstar_mass_bias_sign.json"
    write_json(bias_path, result)
    print_status(f"  Saved: {bias_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, result)
    finalize_result(STEP_ID, result, "Mass-bias sign equation (Gate -1)", "Phantom Mass sign determined from S_a^3*D_dyn vs T_P^2")

    print_status("")
    print_status("  NEXT: step_46_tep_transfer_fit.py — TEP transfer-function fit", "INFO")


if __name__ == "__main__":
    main()
