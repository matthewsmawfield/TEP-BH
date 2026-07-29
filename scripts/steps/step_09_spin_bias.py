#!/usr/bin/env python3
"""Step 09: Spin Bias Analysis.

Compiles published black-hole spin measurements from X-ray reflection
and continuum-fitting spectroscopy and analyses the observed clustering
of spin values near maximal Kerr.

TEP does not modify exterior spectral measurements. The conformal factor
A approx 1 in the exterior and the disformal term B approx 0, so spin
measurements derived from X-ray reflection or continuum-fitting are
unaffected by TEP. The observed spin clustering is not explained by TEP.

Outputs:
  - data/processed/bh_spin_measurements.json
  - results/step_09_spin_bias.json
  - results/step_09_spin_bias.csv
  - logs/step_09_spin_bias.log

Usage:
    python scripts/steps/step_09_spin_bias.py
"""

from __future__ import annotations

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

try:
    from scipy import stats as sp_stats
    SCIPY_AVAILABLE = True
except ImportError:  # pragma: no cover - scipy is a soft dependency
    sp_stats = None
    SCIPY_AVAILABLE = False

STEP_ID = "step_09_spin_bias"


# ---------------------------------------------------------------------------
# Published X-ray reflection spin measurements
# ---------------------------------------------------------------------------
# Each entry: source, spin (a/M), 1-sigma uncertainty, method, environment,
# and the bibliographic reference.  Values are taken from the continuum-fitting
# and X-ray reflection spectroscopy literature as compiled in
# McClintock, Narayan & Steiner (2014) and Reynolds (2021).
SPIN_MEASUREMENTS = [
    {
        "source": "GRS 1915+105",
        "spin_measured": 0.98,
        "spin_uncertainty": 0.02,
        "method": "continuum-fitting",
        "environment": "dense",
        "reference": "McClintock et al. 2006",
    },
    {
        "source": "Cygnus X-1",
        "spin_measured": 0.97,
        "spin_uncertainty": 0.02,
        "method": "continuum-fitting",
        "environment": "dense",
        "reference": "Gou et al. 2011, ApJ 742, 85",
    },
    {
        "source": "LMC X-1",
        "spin_measured": 0.92,
        "spin_uncertainty": 0.05,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Gou et al. 2009",
    },
    {
        "source": "M33 X-7",
        "spin_measured": 0.84,
        "spin_uncertainty": 0.05,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Liu et al. 2008",
    },
    {
        "source": "4U 1543-47",
        "spin_measured": 0.80,
        "spin_uncertainty": 0.05,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Park et al. 2004",
    },
    {
        "source": "GRO J1655-40",
        "spin_measured": 0.70,
        "spin_uncertainty": 0.05,
        "method": "continuum-fitting",
        "environment": "dense",
        "reference": "Shafee et al. 2006",
    },
    {
        "source": "A0620-00",
        "spin_measured": 0.12,
        "spin_uncertainty": 0.19,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Gou et al. 2010",
    },
    {
        "source": "XTE J1550-564",
        "spin_measured": 0.34,
        "spin_uncertainty": 0.15,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Steiner et al. 2011",
    },
    {
        "source": "H 1743-322",
        "spin_measured": 0.20,
        "spin_uncertainty": 0.30,
        "method": "continuum-fitting",
        "environment": "isolated",
        "reference": "Steiner et al. 2011",
    },
    {
        "source": "1E 1743.1-2843",
        "spin_measured": 0.95,
        "spin_uncertainty": 0.05,
        "method": "X-ray reflection",
        "environment": "dense",
        "reference": "Porquet et al. 2018",
    },
    {
        "source": "M87*",
        "spin_measured": 0.90,
        "spin_uncertainty": 0.10,
        "method": "EHT shadow",
        "environment": "dense",
        "reference": "EHT 2019",
    },
    {
        "source": "Sgr A*",
        "spin_measured": 0.50,
        "spin_uncertainty": 0.20,
        "method": "EHT shadow",
        "environment": "dense",
        "reference": "EHT 2022",
    },
]

DATA_SOURCES = {
    "mcclintock_2014": {
        "citation": "McClintock, Narayan & Steiner 2014, arXiv:1303.1583 (Space Sci Rev)",
        "description": "Review of black-hole spin via continuum fitting and reflection spectroscopy",
    },
    "reynolds_2021": {
        "citation": "Reynolds 2021, arXiv:2104.10300",
        "description": "Review of X-ray reflection spin measurements and systematics",
    },
    "gou_2011": {
        "citation": "Gou et al. 2011, ApJ 742, 85",
        "description": "Cygnus X-1 spin via continuum fitting",
    },
    "eht_2019": {
        "citation": "Event Horizon Telescope 2019",
        "description": "M87* shadow and spin constraint",
    },
    "eht_2022": {
        "citation": "Event Horizon Telescope 2022",
        "description": "Sgr A* shadow and spin constraint",
    },
}


def descriptive_stats(spins: list[float]) -> dict:
    """Return mean, median, std, min, max for a list of spins."""
    arr = np.asarray(spins, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "fraction_gt_0.9": float(np.mean(arr > 0.9)),
    }


def ks_test_uniform(spins: list[float]) -> dict:
    """One-sample KS test against Uniform[0, 1].

    Returns the statistic and p-value.  Falls back to None if scipy is
    unavailable.
    """
    if not SCIPY_AVAILABLE:
        return {"ks_statistic": None, "p_value": None, "scipy_available": False}
    arr = np.asarray(spins, dtype=float)
    result = sp_stats.kstest(arr, "uniform", args=(0.0, 1.0))
    return {
        "ks_statistic": float(result.statistic),
        "p_value": float(result.pvalue),
        "scipy_available": True,
    }


def compile_spin_data() -> list[dict]:
    """Compile the published spin measurements into processed records."""
    records = []
    for m in SPIN_MEASUREMENTS:
        records.append({
            "source": m["source"],
            "spin_measured": m["spin_measured"],
            "spin_uncertainty": m["spin_uncertainty"],
            "method": m["method"],
            "environment": m["environment"],
            "reference": m["reference"],
        })
    return records


def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 09: Spin Bias Analysis", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status(f"scipy available: {SCIPY_AVAILABLE}", "INFO")
    print_status("")

    # --- Compile published spin measurements ---
    print_status("Compiling published X-ray spin measurements...", "PROCESS")
    records = compile_spin_data()

    spin_json_path = PROCESSED_DIR / "bh_spin_measurements.json"
    spin_payload = {
        "description": "Compiled black-hole spin measurements from X-ray reflection "
                       "and continuum-fitting literature",
        "n_sources": len(records),
        "sources": records,
        "citations": DATA_SOURCES,
    }
    write_json(spin_json_path, spin_payload)
    print_status(f"Spin measurements saved to {rel(spin_json_path)}", "SUCCESS")
    print_status(f"  {len(records)} sources compiled", "INFO")
    print_status("")

    # --- Standard (measured) distribution ---
    print_status("Standard (Measured) Spin Distribution", "TITLE")
    measured_spins = [r["spin_measured"] for r in records]
    std_stats = descriptive_stats(measured_spins)
    std_ks = ks_test_uniform(measured_spins)
    n_near_max = int(np.sum(np.asarray(measured_spins) > 0.9))
    n_over_maximal = int(np.sum(np.asarray(measured_spins) > 0.95))

    print_status(f"  N sources:            {std_stats['n']}", "INFO")
    print_status(f"  Mean a/M:             {std_stats['mean']:.3f}", "INFO")
    print_status(f"  Median a/M:           {std_stats['median']:.3f}", "INFO")
    print_status(f"  Std dev:              {std_stats['std']:.3f}", "INFO")
    print_status(f"  Min / Max:            {std_stats['min']:.2f} / {std_stats['max']:.2f}", "INFO")
    print_status(f"  Fraction with a/M>0.9: {std_stats['fraction_gt_0.9']:.3f} "
                 f"({n_near_max}/{std_stats['n']})", "INFO")
    print_status(f"  Over-maximal spin fraction (a/M>0.95): {n_over_maximal}/{std_stats['n']} = {n_over_maximal/std_stats['n']:.3f}", "INFO")
    print_status(f"  All measured spins: {sorted(measured_spins)}", "DEBUG")
    print_status(f"  Spin distribution: {n_near_max}/{std_stats['n']} near maximal (>0.9), {n_over_maximal}/{std_stats['n']} over-maximal (>0.95)", "DEBUG")
    if std_ks["ks_statistic"] is not None:
        print_status(f"  KS statistic vs U[0,1]: {std_ks['ks_statistic']:.4f}", "INFO")
        print_status(f"  KS p-value:              {std_ks['p_value']:.4f}", "INFO")
        print_status(f"  KS critical value (alpha=0.05, n={std_stats['n']}): {1.36/np.sqrt(std_stats['n']):.4f}", "DEBUG")
        print_status(f"  KS test: statistic={std_ks['ks_statistic']:.4f} vs critical={1.36/np.sqrt(std_stats['n']):.4f}", "DEBUG")
        if std_ks["p_value"] < 0.05:
            print_status("  Distribution INCONSISTENT with uniform [0,1] (p < 0.05)", "WARNING")
            print_status(f"  Spin clustering is statistically significant (p={std_ks['p_value']:.4f})", "DEBUG")
        else:
            print_status("  Distribution consistent with uniform [0,1]", "INFO")
    else:
        print_status("  scipy not available; KS test skipped", "WARNING")
    print_status("")

    # --- TEP assessment ---
    print_status("TEP Assessment", "TITLE")
    print_status("  TEP does not modify exterior spectral measurements.", "INFO")
    print_status("  The conformal factor A ~ 1 in the exterior and the disformal", "INFO")
    print_status("  term B ~ 0, so spin measurements from X-ray reflection and", "INFO")
    print_status("  continuum-fitting are unaffected by TEP.", "INFO")
    print_status("  The observed spin clustering is not explained by TEP.", "INFO")
    print_status("")

    # --- Per-source listing ---
    print_status("Per-source spin measurements:", "INFO")
    for r in records:
        print_status(
            f"    {r['source']:<18s} a_meas={r['spin_measured']:.2f} "
            f"+/- {r['spin_uncertainty']:.2f}  method={r['method']}",
            "INFO",
        )
        print_status(f"      reference: {r['reference']}", "DEBUG")
        print_status(f"      environment: {r['environment']}", "DEBUG")
        is_high = r['spin_measured'] > 0.9
        print_status(f"      {'[NEAR MAXIMAL]' if is_high else ''} spin={r['spin_measured']:.2f} ± {r['spin_uncertainty']:.2f}", "DEBUG")
    print_status("")
    print_status(f"  Over-maximal spin analysis: {n_over_maximal}/{std_stats['n']} sources have a/M > 0.95", "INFO")
    print_status(f"  Near-maximal clustering: {n_near_max}/{std_stats['n']} sources have a/M > 0.90", "INFO")
    print_status(f"  TEP spin correction: None (A ≈ 1 in exterior, B ≈ 0; no spectral modification)", "INFO")
    print_status(f"  TEP spin correction formula: a_corrected = a_measured (no correction needed)", "DEBUG")
    print_status(f"  Result: observed spin clustering is not explained by TEP", "DEBUG")

    # --- Save CSV ---
    csv_rows = [
        {
            "source": r["source"],
            "spin_measured": r["spin_measured"],
            "spin_uncertainty": r["spin_uncertainty"],
            "method": r["method"],
            "environment": r["environment"],
            "reference": r["reference"],
        }
        for r in records
    ]
    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, csv_rows)
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    tep_note = (
        "TEP does not modify exterior spectral measurements. The conformal "
        "factor A ~ 1 in the exterior and the disformal term B ~ 0, so spin "
        "measurements derived from X-ray reflection or continuum-fitting are "
        "unaffected by TEP. The observed spin clustering is not explained by TEP."
    )
    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "data_sources": DATA_SOURCES,
        "spin_measurements": records,
        "observed_distribution": {
            **std_stats,
            "n_near_maximal": n_near_max,
            "ks_test": std_ks,
        },
        "tep_assessment": {
            "modifies_exterior_spectra": False,
            "conformal_factor_exterior": "~1",
            "disformal_term_exterior": "~0",
            "spin_correction": None,
            "note": tep_note,
        },
        "note": tep_note,
    }
    json_path = step_json_path(STEP_ID)
    summary = finalize_result(
        STEP_ID, summary,
        description="Compile published BH spin measurements from X-ray reflection and continuum-fitting spectroscopy and analyze the observed clustering near maximal Kerr",
        key_result=f"Spin distribution (n={std_stats['n']}, mean={std_stats['mean']:.3f}) shows {n_near_max}/{std_stats['n']} near-maximal sources; KS test p={std_ks['p_value']:.4f} vs uniform; TEP does not modify exterior spin measurements",
        dependencies=["step_00_data_download"],
    )
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 09 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
