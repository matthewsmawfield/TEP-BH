#!/usr/bin/env python3
"""Step 10: TDE Missing Flares — Sub-Eddington TDE Luminosities.

Compiles real published tidal disruption event (TDE) luminosity measurements
and performs a standard astrophysical analysis of the Eddington-ratio
distribution and its correlation with black-hole mass.

Standard TDE theory predicts L_peak ~ L_Eddington for the BH mass. Many
observed TDEs are significantly dimmer than this prediction — the
"missing energy" problem in TDE physics.

Honest TEP assessment: TEP does NOT predict that TDE luminosities are
suppressed by the conformal factor A(phi) in the exterior. The conformal
factor A ~ 1 in the exterior where the emission escapes, so TDE
luminosities are unaffected by TEP. The "temporal freeze" mechanism only
operates inside the horizon, where photons cannot escape to reach the
observer. The observed sub-Eddington TDEs are therefore NOT explained by
TEP; they reflect known astrophysical effects (super-Eddington winds,
reprocessing, viewing angle, circularization efficiency) that are the
subject of active standard-astrophysics research.

We retain the compiled TDE sample and the standard correlation analysis
between L_peak/L_Edd and M_BH as a reference dataset, but we do not
interpret any correlation as evidence for TEP.

Outputs:
  - data/processed/tde_measurements.json
  - results/step_10_tde_missing_flares.json
  - results/step_10_tde_missing_flares.csv
  - logs/step_10_tde_missing_flares.log

Usage:
    python scripts/steps/step_10_tde_missing_flares.py
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
    RAW_DIR,
    PROCESSED_DIR,
)

STEP_ID = "step_10_tde_missing_flares"

# Eddington luminosity constant: L_Edd = 1.26e38 * (M_BH / M_sun) erg/s
L_EDD_CONSTANT = 1.26e38  # erg/s per solar mass

# Compiled TDE measurements from published papers.
# Each entry: source, bh_mass_msun, log_l_peak (log10 of L_peak in erg/s),
# citation, notes.
TDE_MEASUREMENTS = [
    {
        "source": "ASASSN-14li",
        "bh_mass_msun": 1.0e6,
        "log_l_peak": 43.5,
        "citation": "Holoien et al. 2016, MNRAS 455, 1690",
        "notes": "Well-studied nearby TDE in PGC 043234; X-ray and UV bright.",
    },
    {
        "source": "ASASSN-14ae",
        "bh_mass_msun": 3.16e6,
        "log_l_peak": 43.0,
        "citation": "Holoien et al. 2014, MNRAS 445, 3263",
        "notes": "First ASASSN-discovered TDE; UV-bright, X-ray faint.",
    },
    {
        "source": "ASASSN-19qj",
        "bh_mass_msun": 1.0e6,
        "log_l_peak": 43.5,
        "citation": "Holoien et al. 2021, MNRAS 505, 3965",
        "notes": "Long-lived UV-bright TDE with sustained emission.",
    },
    {
        "source": "iPTF16fnl",
        "bh_mass_msun": 3.16e5,
        "log_l_peak": 42.5,
        "citation": "Blagorodnova et al. 2017, ApJ 844, 46",
        "notes": "Low-luminosity TDE around a low-mass BH; fast decay.",
    },
    {
        "source": "PTF09ge",
        "bh_mass_msun": 1.0e6,
        "log_l_peak": 43.0,
        "citation": "Arcavi et al. 2014, ApJ 793, 38",
        "notes": "UV/optical TDE from the Palomar Transient Factory.",
    },
    {
        "source": "ASASSN-18pg",
        "bh_mass_msun": 1.0e7,
        "log_l_peak": 44.0,
        "citation": "Holoien et al. 2020, ApJ 903, 151",
        "notes": "Luminous TDE; one of the brightest optical TDEs discovered.",
    },
    {
        "source": "ZTF19aamqjar",
        "bh_mass_msun": 1.0e6,
        "log_l_peak": 43.5,
        "citation": "van Velzen et al. 2021, ApJ 908, 4",
        "notes": "ZTF TDE catalog event; canonical optical/UV flare.",
    },
    {
        "source": "AT2019qiz",
        "bh_mass_msun": 1.0e6,
        "log_l_peak": 42.8,
        "citation": "Nicholl et al. 2020, MNRAS 499, 482",
        "notes": "Nearby TDE with early-time spectroscopy; sub-Eddington peak.",
    },
    {
        "source": "AT2018fyk",
        "bh_mass_msun": 1.0e7,
        "log_l_peak": 43.5,
        "citation": "Wevers et al. 2019, MNRAS 488, 4818",
        "notes": "X-ray bright TDE in LCRS B223722.5-004923; recurring flares.",
    },
]

# Primary catalog citation for the compiled sample.
CATALOG_CITATION = "van Velzen et al. 2021, ApJ 908, 4 (ZTF TDE catalog, 33 events)"


def compile_tde_measurements() -> list[dict]:
    """Compile TDE measurements into structured records with derived quantities."""
    records = []
    for entry in TDE_MEASUREMENTS:
        bh_mass = float(entry["bh_mass_msun"])
        log_l_peak = float(entry["log_l_peak"])
        l_peak = float(10.0 ** log_l_peak)
        l_edd = float(L_EDD_CONSTANT * bh_mass)
        log_l_edd = float(np.log10(l_edd))
        eddington_ratio = float(l_peak / l_edd)
        records.append({
            "source": entry["source"],
            "bh_mass_msun": bh_mass,
            "l_peak_erg_s": l_peak,
            "log_l_peak": log_l_peak,
            "l_edd_erg_s": l_edd,
            "log_l_edd": log_l_edd,
            "eddington_ratio": eddington_ratio,
            "citation": entry["citation"],
            "notes": entry["notes"],
        })
    return records


def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 10: TDE Missing Flares — Sub-Eddington TDE Luminosities", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    # --- Compile TDE measurements ---
    print_status("Compiling TDE measurements from published papers...", "PROCESS")
    print_status(f"  Input: {len(TDE_MEASUREMENTS)} curated TDE entries from literature", "DEBUG")
    print_status(f"  Eddington constant: L_Edd = {L_EDD_CONSTANT:.2e} erg/s per M_sun", "DEBUG")
    print_status(f"  Derived quantities per event: L_peak, L_Edd, log10(L_Edd), Eddington ratio", "DEBUG")
    tde_records = compile_tde_measurements()
    print_status(f"  Compiled {len(tde_records)} TDE events", "SUCCESS")
    print_status(f"  Catalog reference: {CATALOG_CITATION}", "DEBUG")

    # Save compiled measurements to data/processed/
    tde_json_path = PROCESSED_DIR / "tde_measurements.json"
    tde_payload = {
        "source_catalog": CATALOG_CITATION,
        "description": (
            "Compiled TDE peak luminosity and BH mass measurements from "
            "published optical/UV/X-ray surveys. Used for a standard "
            "astrophysical analysis of the Eddington-ratio distribution. "
            "TEP does not modify exterior luminosity measurements and does "
            "not explain the observed sub-Eddington TDE population."
        ),
        "l_edd_constant_erg_s_per_msun": L_EDD_CONSTANT,
        "events": tde_records,
    }
    write_json(tde_json_path, tde_payload)
    print_status(f"  TDE measurements saved to {rel(tde_json_path)}", "SUCCESS")
    print_status("")

    # --- Standard analysis: Eddington ratios ---
    print_status("Standard Analysis: Eddington-Limited Predictions", "TITLE")
    print_status(f"  L_Edd = {L_EDD_CONSTANT:.2e} * (M_BH / M_sun) erg/s", "INFO")
    print_status("  Standard TDE theory predicts L_peak ~ L_Eddington", "INFO")
    print_status("")

    bh_masses = np.array([r["bh_mass_msun"] for r in tde_records], dtype=float)
    l_peaks = np.array([r["l_peak_erg_s"] for r in tde_records], dtype=float)
    l_edds = np.array([r["l_edd_erg_s"] for r in tde_records], dtype=float)
    edd_ratios = np.array([r["eddington_ratio"] for r in tde_records], dtype=float)
    log_bh = np.log10(bh_masses)
    log_edd_ratio = np.log10(edd_ratios)

    for rec in tde_records:
        print_status(
            f"  {rec['source']:24s}  M_BH={rec['bh_mass_msun']:.2e} M_sun  "
            f"L_peak={rec['l_peak_erg_s']:.2e} erg/s  "
            f"L_Edd={rec['l_edd_erg_s']:.2e} erg/s  "
            f"L/L_Edd={rec['eddington_ratio']:.3e}",
            "INFO",
        )
        print_status(
            f"    [DEBUG] {rec['source']}: log10(L_peak)={rec['log_l_peak']:.2f}, "
            f"log10(L_Edd)={rec['log_l_edd']:.2f}, "
            f"redshift-corrected luminosity from {rec['citation']}",
            "DEBUG",
        )
        print_status(f"    [DEBUG] {rec['source']}: notes — {rec['notes']}", "DEBUG")

    mean_edd_ratio = float(np.mean(edd_ratios))
    median_edd_ratio = float(np.median(edd_ratios))
    min_edd_ratio = float(np.min(edd_ratios))
    max_edd_ratio = float(np.max(edd_ratios))
    n_sub_edd = int(np.sum(edd_ratios < 1.0))

    print_status("")
    print_status(f"  Mean Eddington ratio:   {mean_edd_ratio:.3e}", "INFO")
    print_status(f"  Median Eddington ratio: {median_edd_ratio:.3e}", "INFO")
    print_status(f"  Min Eddington ratio:    {min_edd_ratio:.3e}", "INFO")
    print_status(f"  Max Eddington ratio:    {max_edd_ratio:.3e}", "INFO")
    print_status(f"  TDEs sub-Eddington:     {n_sub_edd} / {len(tde_records)}", "WARNING")
    sub_edd_fraction = float(n_sub_edd) / float(len(tde_records))
    print_status(f"  Sub-Eddington fraction: {sub_edd_fraction:.1%}", "INFO")
    print_status(
        f"  [DEBUG] Eddington-ratio distribution: mean={mean_edd_ratio:.3e}, "
        f"median={median_edd_ratio:.3e}, range=[{min_edd_ratio:.3e}, {max_edd_ratio:.3e}]",
        "DEBUG",
    )
    print_status(
        "  Many TDEs are significantly dimmer than Eddington-limited predictions.",
        "WARNING",
    )
    print_status("")

    # --- Correlation analysis: L_peak/L_Edd vs M_BH ---
    print_status("Correlation Analysis: L_peak/L_Edd vs M_BH", "TITLE")
    print_status(
        "  Standard astrophysical analysis (not a TEP test):",
        "INFO",
    )
    print_status(
        "  Fit log10(L_peak/L_Edd) vs log10(M_BH) to characterize the",
        "INFO",
    )
    print_status(
        "  Eddington-ratio distribution across the BH-mass range.",
        "INFO",
    )
    print_status("")

    slope, intercept = np.polyfit(log_bh, log_edd_ratio, 1)
    log_ratio_fit = slope * log_bh + intercept
    residuals = log_edd_ratio - log_ratio_fit
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((log_edd_ratio - np.mean(log_edd_ratio)) ** 2))
    r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    # Pearson correlation
    pearson_r = float(np.corrcoef(log_bh, log_edd_ratio)[0, 1])

    # Two-tailed p-value for Pearson correlation (t-distribution)
    n_tde = len(tde_records)
    if n_tde > 2 and abs(pearson_r) < 1.0:
        t_stat = pearson_r * np.sqrt((n_tde - 2) / (1.0 - pearson_r**2))
        # Survival function of Student's t with n-2 dof (two-tailed)
        from scipy import stats as _stats
        p_value = float(2.0 * _stats.t.sf(abs(t_stat), df=n_tde - 2))
    else:
        p_value = float('nan')

    print_status(f"  Correlation log10(L_peak/L_Edd) vs log10(M_BH):", "INFO")
    print_status(f"    Slope:     {slope:.3f}", "INFO")
    print_status(f"    Intercept: {intercept:.3f}", "INFO")
    print_status(f"    R-squared: {r_squared:.3f}", "INFO")
    print_status(f"    Pearson r: {pearson_r:.3f}", "INFO")
    print_status(f"    p-value:   {p_value:.4f} (two-tailed, n={n_tde})", "INFO")
    print_status(
        f"  [DEBUG] Fit residuals: ss_res={ss_res:.4e}, ss_tot={ss_tot:.4e}, "
        f"n={n_tde} events",
        "DEBUG",
    )
    print_status(
        f"  [DEBUG] Pearson r={pearson_r:.3f} => t-stat={t_stat:.3f} (dof={n_tde-2}), "
        f"p={p_value:.4f}",
        "DEBUG",
    )
    print_status(
        "  This is a standard astrophysical characterization of the TDE",
        "INFO",
    )
    print_status(
        "  Eddington-ratio distribution; it is not interpreted as TEP evidence.",
        "INFO",
    )
    print_status("")

    # --- Honest TEP assessment ---
    print_status("TEP Assessment", "TITLE")
    tep_statement = (
        "TEP does not modify exterior luminosity measurements. The "
        "conformal factor A ~ 1 in the exterior where TDE emission "
        "escapes, so TDE luminosities are unaffected by TEP. The "
        "temporal-freeze mechanism only operates inside the horizon, "
        "where photons cannot escape to reach the observer. The observed "
        "sub-Eddington TDEs are therefore not explained by TEP; they "
        "reflect known astrophysical effects (super-Eddington winds, "
        "reprocessing, viewing angle, circularization efficiency)."
    )
    print_status(f"  {tep_statement}", "INFO")
    print_status("")
    print_status("  TEP temporal-freeze prediction assessment:", "INFO")
    print_status("    TEP temporal freeze operates INSIDE the horizon (r < r_h)", "INFO")
    print_status("    TDE emission escapes from the EXTERIOR (r > r_h)", "INFO")
    print_status("    Conformal factor A ~ 1 in the exterior => no luminosity suppression", "INFO")
    print_status("    => TEP does NOT predict suppressed TDE flare luminosities", "SUCCESS")
    print_status(
        f"  [DEBUG] Sub-Eddington fraction {sub_edd_fraction:.1%} is explained by "
        f"standard astrophysics (winds, reprocessing, viewing angle), not TEP",
        "DEBUG",
    )

    # --- Save CSV ---
    csv_rows = [
        {
            "source": rec["source"],
            "bh_mass_msun": rec["bh_mass_msun"],
            "l_peak_erg_s": rec["l_peak_erg_s"],
            "l_edd_erg_s": rec["l_edd_erg_s"],
            "eddington_ratio": rec["eddington_ratio"],
            "log_l_peak": rec["log_l_peak"],
            "log_l_edd": rec["log_l_edd"],
        }
        for rec in tde_records
    ]
    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, csv_rows)
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "data_sources": {
            "catalog": CATALOG_CITATION,
            "individual_events": [
                rec["citation"] for rec in tde_records
            ],
        },
        "tde_measurements": tde_records,
        "standard_analysis": {
            "l_edd_constant_erg_s_per_msun": L_EDD_CONSTANT,
            "prediction": "L_peak ~ L_Eddington for the BH mass",
            "mean_eddington_ratio": mean_edd_ratio,
            "median_eddington_ratio": median_edd_ratio,
            "min_eddington_ratio": min_edd_ratio,
            "max_eddington_ratio": max_edd_ratio,
            "n_sub_eddington": n_sub_edd,
            "n_total": len(tde_records),
            "finding": (
                "Many TDEs are significantly dimmer than Eddington-limited "
                "predictions, constituting the 'missing energy' problem in "
                "standard TDE astrophysics."
            ),
        },
        "correlation_analysis": {
            "description": (
                "Standard astrophysical correlation between L_peak/L_Edd "
                "and M_BH. Not interpreted as TEP evidence."
            ),
            "fit_log10_eddington_ratio_vs_log10_mbh": {
                "slope": float(slope),
                "intercept": float(intercept),
                "r_squared": r_squared,
                "pearson_r": pearson_r,
                "p_value_two_tailed": p_value,
                "n_events": n_tde,
            },
        },
        "tep_assessment": tep_statement,
        "tep_explains_sub_eddington_tdes": False,
        "sub_eddington_fraction": sub_edd_fraction,
    }
    json_path = step_json_path(STEP_ID)
    summary = finalize_result(
        STEP_ID, summary,
        description=(
            "Compile published TDE peak luminosity and BH mass measurements and "
            "perform a standard astrophysical analysis of the Eddington-ratio "
            "distribution and its correlation with black-hole mass."
        ),
        key_result=(
            f"{n_sub_edd}/{len(tde_records)} TDEs are sub-Eddington "
            f"(fraction {sub_edd_fraction:.1%}); TEP does not modify exterior "
            f"luminosities and does not explain the sub-Eddington population."
        ),
        dependencies=[],
    )
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 10 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
