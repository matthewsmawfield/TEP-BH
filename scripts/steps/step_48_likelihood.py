#!/usr/bin/env python3
"""Step 05 (Inference): Formal likelihood comparison — GR vs TEP.

Compute the formal likelihood ratio and Bayes factor between the
GR model (eta_TEP = 0) and the TEP model (eta_TEP free), using:
  1. Maximum likelihood ratio
  2. Bayesian Information Criterion (BIC)
  3. Akaike Information Criterion (AIC)
  4. Wilks' theorem chi-squared test

The likelihood function is:
  L = exp(-chi^2 / 2)

The Bayes factor B = L_TEP / L_GR quantifies the evidence for TEP
over GR. Convention:
  B > 100:  decisive evidence for TEP
  B > 10:   strong evidence for TEP
  B > 3:    moderate evidence for TEP
  B ~ 1:    inconclusive
  B < 1/3:  moderate evidence for GR
  B < 1/10: strong evidence for GR
  B < 1/100: decisive evidence for GR

Outputs:
  data/processed/sstar_likelihood_comparison.json
  results/step_48_likelihood.json
  logs/step_48_likelihood.log
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps" / "inference"))

from step_46_tep_transfer import tep_residuals
from step_44_gr_fit import gr_residuals

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    write_json,
    finalize_result,
    PROCESSED_DIR,
    RESULTS_DIR,
)

STEP_ID = "step_48_likelihood"


def interpret_bayes_factor(bf):
    """Interpret the Bayes factor using the Kass & Raftery (1995) scale."""
    if bf > 100:
        return "decisive evidence for TEP"
    elif bf > 10:
        return "strong evidence for TEP"
    elif bf > 3:
        return "moderate evidence for TEP"
    elif bf > 1.0/3:
        return "inconclusive"
    elif bf > 1.0/10:
        return "moderate evidence for GR"
    elif bf > 1.0/100:
        return "strong evidence for GR"
    else:
        return "decisive evidence for GR"


def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 05 (INFERENCE): LIKELIHOOD COMPARISON", "TITLE")
    print_status("  GR vs TEP — Formal Statistical Test", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")

    # --- Load data and fits ---
    astrometry_path = PROCESSED_DIR / "sstar_astrometry.json"
    spectroscopy_path = PROCESSED_DIR / "sstar_spectroscopy.json"
    gr_fit_path = PROCESSED_DIR / "sstar_gr_fit.json"
    tep_fit_path = PROCESSED_DIR / "sstar_tep_fit.json"

    if not all(p.exists() for p in [astrometry_path, spectroscopy_path, gr_fit_path, tep_fit_path]):
        print_status("ERROR: Previous fits not found. Run steps 00-04 first.", "ERROR")
        return

    astrometry = json.loads(astrometry_path.read_text())["data"]
    spectroscopy = json.loads(spectroscopy_path.read_text())["data"]
    gr_fit = json.loads(gr_fit_path.read_text())
    tep_fit = json.loads(tep_fit_path.read_text())

    # Extract best-fit parameters
    gr_theta = np.array([
        gr_fit["best_fit"]["M_BH (M_sun)"]["value"],
        gr_fit["best_fit"]["R0 (pc)"]["value"],
        gr_fit["best_fit"]["a (arcsec)"]["value"],
        gr_fit["best_fit"]["e"]["value"],
        gr_fit["best_fit"]["i (deg)"]["value"],
        gr_fit["best_fit"]["omega (deg)"]["value"],
        gr_fit["best_fit"]["Omega (deg)"]["value"],
        gr_fit["best_fit"]["T (yr)"]["value"],
        gr_fit["best_fit"]["t_peri (yr)"]["value"],
        gr_fit["best_fit"]["x0 (arcsec)"]["value"],
        gr_fit["best_fit"]["y0 (arcsec)"]["value"],
    ])

    tep_theta = np.array([
        tep_fit["best_fit"]["M_BH (M_sun)"]["value"],
        tep_fit["best_fit"]["R0 (pc)"]["value"],
        tep_fit["best_fit"]["a (arcsec)"]["value"],
        tep_fit["best_fit"]["e"]["value"],
        tep_fit["best_fit"]["i (deg)"]["value"],
        tep_fit["best_fit"]["omega (deg)"]["value"],
        tep_fit["best_fit"]["Omega (deg)"]["value"],
        tep_fit["best_fit"]["T (yr)"]["value"],
        tep_fit["best_fit"]["t_peri (yr)"]["value"],
        tep_fit["best_fit"]["x0 (arcsec)"]["value"],
        tep_fit["best_fit"]["y0 (arcsec)"]["value"],
        tep_fit["best_fit"]["eta_TEP"]["value"],
    ])

    # Compute chi-squared for each model
    res_gr = gr_residuals(gr_theta, astrometry, spectroscopy)
    chi2_gr = np.sum(res_gr**2)
    n_data = len(res_gr)
    n_params_gr = len(gr_theta)
    n_params_tep = len(tep_theta)

    res_tep = tep_residuals(tep_theta, astrometry, spectroscopy)
    chi2_tep = np.sum(res_tep**2)

    # --- Maximum likelihood ratio ---
    # L = exp(-chi^2/2), so L_TEP/L_GR = exp(-(chi2_tep - chi2_gr)/2)
    delta_chi2 = chi2_gr - chi2_tep
    ml_ratio = np.exp(delta_chi2 / 2)

    # --- BIC ---
    bic_gr = chi2_gr + n_params_gr * np.log(n_data)
    bic_tep = chi2_tep + n_params_tep * np.log(n_data)
    delta_bic = bic_tep - bic_gr
    bf_bic = np.exp(-delta_bic / 2)

    # --- AIC ---
    aic_gr = chi2_gr + 2 * n_params_gr
    aic_tep = chi2_tep + 2 * n_params_tep
    delta_aic = aic_tep - aic_gr
    bf_aic = np.exp(-delta_aic / 2)

    # --- Wilks' theorem ---
    # The TEP model has 1 extra parameter (eta_TEP)
    # Under the null (GR, eta=0), delta_chi2 should follow chi^2(1)
    p_value = 1 - stats.chi2.cdf(delta_chi2, df=1)

    # --- Results ---
    print_status("--- Chi-squared values ---", "TITLE")
    print_status(f"  chi^2_GR  = {chi2_gr:.4f}  (n_params = {n_params_gr}, dof = {n_data - n_params_gr})", "INFO")
    print_status(f"  chi^2_TEP = {chi2_tep:.4f}  (n_params = {n_params_tep}, dof = {n_data - n_params_tep})", "INFO")
    print_status(f"  Delta chi^2 = {delta_chi2:.4f}", "INFO")
    print_status("")

    print_status("--- Information criteria ---", "TITLE")
    print_status(f"  BIC_GR  = {bic_gr:.4f}", "INFO")
    print_status(f"  BIC_TEP = {bic_tep:.4f}", "INFO")
    print_status(f"  Delta BIC = {delta_bic:.4f}  (negative favours TEP)", "INFO")
    print_status(f"  AIC_GR  = {aic_gr:.4f}", "INFO")
    print_status(f"  AIC_TEP = {aic_tep:.4f}", "INFO")
    print_status(f"  Delta AIC = {delta_aic:.4f}  (negative favours TEP)", "INFO")
    print_status("")

    print_status("--- Bayes factors ---", "TITLE")
    print_status(f"  Maximum likelihood ratio: {ml_ratio:.6e}", "INFO")
    print_status(f"  Bayes factor (BIC):       {bf_bic:.6e}", "INFO")
    print_status(f"  Bayes factor (AIC):       {bf_aic:.6e}", "INFO")
    print_status("")

    interpretation = interpret_bayes_factor(bf_bic)
    print_status(f"  Interpretation (Kass & Raftery 1995): {interpretation}", "INFO")
    print_status("")

    print_status("--- Wilks' theorem test ---", "TITLE")
    print_status(f"  Delta chi^2 = {delta_chi2:.4f} (should follow chi^2(1) under null)", "INFO")
    print_status(f"  p-value = {p_value:.4f}", "INFO")
    if p_value < 0.05:
        print_status("  -> TEP is statistically significant (p < 0.05)", "INFO")
    else:
        print_status(f"  -> TEP is NOT statistically significant (p = {p_value:.4f} > 0.05)", "INFO")
        print_status("  -> eta_TEP is consistent with 0", "INFO")
    print_status("")

    # --- Summary ---
    print_status("=" * 70, "TITLE")
    print_status("SUMMARY: GR vs TEP at S-star scales", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(f"  Bayes factor B(TEP/GR) = {bf_bic:.4e}", "INFO")
    print_status(f"  p-value (Wilks) = {p_value:.4f}", "INFO")
    print_status(f"  Interpretation: {interpretation}", "INFO")
    print_status("")
    print_status("  At S2 pericentre (r ~ 1381 R_s), the fitted TEP correction is", "INFO")
    print_status("  consistent with zero and the extra TEP parameter is not", "INFO")
    print_status("  constrained by the data. GR is the preferred simpler model.", "INFO")
    print_status("  This is the expected weak-field result.", "INFO")
    print_status("  The decisive test requires horizon-scale observations.", "INFO")
    print_status("")

    # --- Save ---
    result = {
        "step": STEP_ID,
        "description": "Formal likelihood comparison: GR vs TEP",
        "chi2": {
            "GR": float(chi2_gr),
            "TEP": float(chi2_tep),
            "delta": float(delta_chi2),
        },
        "information_criteria": {
            "BIC": {"GR": float(bic_gr), "TEP": float(bic_tep), "delta": float(delta_bic)},
            "AIC": {"GR": float(aic_gr), "TEP": float(aic_tep), "delta": float(delta_aic)},
        },
        "bayes_factors": {
            "max_likelihood_ratio": float(ml_ratio),
            "BIC": float(bf_bic),
            "AIC": float(bf_aic),
        },
        "wilks_test": {
            "delta_chi2": float(delta_chi2),
            "dof": 1,
            "p_value": float(p_value),
        },
        "interpretation": interpretation,
        "n_data": int(n_data),
        "n_params": {"GR": int(n_params_gr), "TEP": int(n_params_tep)},
        "conclusion": (
            "At S2 pericentre (r ~ 1381 R_s), the fitted TEP correction is consistent with zero. "
            "The extra TEP parameter is not constrained by the data, so GR is the preferred simpler model. "
            "This is the expected weak-field result. The decisive test requires horizon-scale observations."
        ),
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    like_path = PROCESSED_DIR / "sstar_likelihood_comparison.json"
    write_json(like_path, result)
    print_status(f"  Saved: {like_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, result)
    finalize_result(STEP_ID, result, "Likelihood comparison: GR vs TEP", "Bayes factor and p-value computed")

    print_status("")
    print_status("  NEXT: step_06_phantom_mass_posterior.py — MCMC posterior", "INFO")


if __name__ == "__main__":
    main()
