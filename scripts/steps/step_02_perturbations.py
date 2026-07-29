#!/usr/bin/env python3
"""Step 02: Perturbation Analysis.

Computes the supported nonspinning gravitational Regge-Wheeler and Zerilli
effective potentials and quasinormal-mode frequencies for canonical TEP.
Gravitational waves propagate on the geometric metric g, prescribed here as
Schwarzschild; the disformal matter metric is not used for GW predictions.

The perturbation equation for axial gravitational perturbations of the
geometric metric g_{mu nu} is the Regge-Wheeler equation:

    d^2 Psi / dr*^2 + [omega^2 - V_RW(r*)] Psi = 0

where r* is the Schwarzschild tortoise coordinate and V_RW is the effective
potential.  For a general spherically symmetric geometric metric

    ds^2 = -f dt^2 + f^{-1} dr^2 + R^2 dOmega^2

the Regge-Wheeler potential is:

    V_RW = f * [l(l+1)/R^2 - 6*M_eff/R^3 * dR/dr]

where M_eff = (R/2)(1 - (dR/dr)^2 * f - R * d^2R/dr^2 * f + R * (dR/dr)^2 * f'/2f)
is the effective mass function.

For this repository's geometric prescription:
    f = F = 1 - 2M/r, R = r, M_eff = M
    V_RW = F * [l(l+1)/r^2 - 6M/r^3]  (standard Schwarzschild RW)

Matter-metric conformal/disformal profiles remain background diagnostics only;
they do not rescale the gravitational potential or generate GW echoes.

Outputs (prefixed with step_02_perturbations):
  - results/step_02_perturbations.json   (summary diagnostics)
  - results/step_02_perturbations.csv     (potential profiles)
  - logs/step_02_perturbations.log        (verbose log)

Usage:
    python scripts/steps/step_02_perturbations.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

import numpy as np
from scipy.integrate import cumulative_trapezoid

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    step_json_path,
    step_csv_path,
    rel,
    write_json,
    finalize_result,
    TEPBHModel,
    solve_tep_bh,
)

STEP_ID = "step_02_perturbations"


# =============================================================================
# Tortoise coordinate for the geometric Schwarzschild metric
# =============================================================================

def compute_tortoise_coordinate(r, metric):
    """Compute r* for the prescribed geometric Schwarzschild metric.

    The supported gravitational perturbation problem uses
    dr*/dr = 1/|F| on the sampled grid, with F = 1 - 2M/r.  The absolute
    value merely permits the background grid to extend through r=2M;
    QNM calculations below use only the static exterior.
    """
    F = metric['F']
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
# Regge-Wheeler potential for gravitational waves on geometric g
# =============================================================================

def compute_regge_wheeler_potential(r, r_star, metric, l=2):
    """Compute the axial gravitational potential on geometric g.

    Canonical TEP places gravitational waves on g, and this repository
    prescribes g as Schwarzschild.  Therefore

        V_RW = F [l(l+1)/r^2 - 6M/r^3],

    with no matter-metric A^{-4} factor.  The second returned array is an
    explicit Schwarzschild reference and is identical by construction.
    """
    F = metric['F']
    M = 1.0  # geometric units
    F_safe = np.where(np.abs(F) > 1e-30, F, np.nan)
    V_RW_geometric = F_safe * (l * (l + 1) / r**2 - 6.0 * M / r**3)
    V_RW_geometric = np.where(np.isfinite(V_RW_geometric), V_RW_geometric, 0)

    return V_RW_geometric, V_RW_geometric.copy()


# =============================================================================
# Zerilli potential
# =============================================================================

def compute_zerilli_potential(r, r_star, metric, l=2):
    """Compute the polar gravitational potential on geometric g.

    The geometric metric is Schwarzschild, so the standard Zerilli
    potential is used without matter-metric conformal scaling.  The second
    returned array is an identical Schwarzschild reference.
    """
    F = metric['F']
    lam = (l - 1) * (l + 2) / 2.0
    M = 1.0

    F_safe = np.where(np.abs(F) > 1e-30, F, np.nan)
    numerator = 2.0 * F_safe * (
        lam**2 * (lam + 1) * r**3
        + 3.0 * lam**2 * M * r**2
        + 9.0 * lam * M**2 * r
        + 9.0 * M**3
    )
    denominator = r**3 * (lam * r + 3.0 * M)**2
    V_Z_geometric = numerator / np.where(denominator > 1e-50, denominator, np.nan)
    V_Z_geometric = np.where(np.isfinite(V_Z_geometric), V_Z_geometric, 0)

    return V_Z_geometric, V_Z_geometric.copy()


# =============================================================================
# QNM via WKB approximation
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
# Echo cavity check
# =============================================================================

def check_echo_cavity(r, V_RW, r_star=None, r_temporal_horizon=None, idx_temporal_horizon=None):
    """Check the supported exterior geometric potential for two barriers.

    A Schwarzschild RW potential has one exterior photon-sphere barrier.
    This routine deliberately makes no GW inference from the matter metric
    or from the interior structure (the TEP matter metric is globally
    Lorentzian with no boundary, but GWs propagate on the geometric metric).
    The r_temporal_horizon / idx_temporal_horizon arguments are retained for
    API compatibility but are expected to be None under the current
    Gaussian-bump model (no determinant-zero boundary).
    """
    idx_h = np.argmin(np.abs(r - 2.0))
    V_ext = V_RW[idx_h:]
    V_max = float(np.max(V_ext)) if len(V_ext) else 0.0
    peak_threshold = max(V_max * 0.01, 1e-10)
    n_peaks_ext = sum(
        V_ext[i] > V_ext[i - 1]
        and V_ext[i] > V_ext[i + 1]
        and V_ext[i] > peak_threshold
        for i in range(1, len(V_ext) - 1)
    )
    has_exterior_cavity = n_peaks_ext >= 2

    return {
        "domain": "geometric_schwarzschild_exterior",
        "n_peaks_exterior": n_peaks_ext,
        "has_exterior_echo_cavity": has_exterior_cavity,
        "gravitational_echo_prediction_supported": False,
        "note": "No matter-metric or black-hole-interior structure is used as a GW echo prediction.",
    }


# =============================================================================
# Main
# =============================================================================

def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status(f"STEP 02: Perturbation Analysis", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    model = TEPBHModel(
        beta_A=-1.0, B0=1.0, n_B=2.0,
        phi_0=2.0, delta=0.05, M=1.0,
    )
    print_status(f"Model parameters: {model.to_dict()}", "INFO")

    print_status("Solving TEP-BH background...", "PROCESS")
    solution = solve_tep_bh(model)
    if not solution["success"]:
        print_status("BACKGROUND SOLUTION FAILED", "ERROR")
        return {"step": STEP_ID, "status": "failed"}

    r = solution["r"]
    metric = solution["metric"]

    print_status(f"  Grid: {len(r)} points, r in [{r[0]:.6e}, {r[-1]:.6e}]", "DEBUG")
    print_status(f"  F at r=3M: {metric['F'][np.argmin(np.abs(r - 3.0))]:.6f} (Schwarzschild metric function)", "DEBUG")
    print_status(f"  A at r=3M: {metric['A'][np.argmin(np.abs(r - 3.0))]:.6f} (conformal factor, ~1 in exterior)", "DEBUG")

    # --- Tortoise coordinate ---
    print_status("Computing tortoise coordinate...", "PROCESS")
    r_star = compute_tortoise_coordinate(r, metric)
    print_status(f"  r* range: [{r_star[0]:.4e}, {r_star[-1]:.4e}]", "INFO")
    print_status(f"  r* at r=3M: {r_star[np.argmin(np.abs(r - 3.0))]:.4e}", "INFO")
    print_status(f"  dr*/dr = 1/|F| (Schwarzschild tortoise, A cancels)", "DEBUG")
    print_status(f"  r* at r=2M (horizon): {r_star[np.argmin(np.abs(r - 2.0))]:.4e} (-> -inf)", "DEBUG")
    print_status(f"  r* at r=10M: {r_star[np.argmin(np.abs(r - 10.0))]:.4e}", "DEBUG")

    # --- Regge-Wheeler potential ---
    print_status("Computing Regge-Wheeler potential (l=2)...", "PROCESS")
    V_RW, V_RW_schw = compute_regge_wheeler_potential(r, r_star, metric, l=2)

    idx_3M = np.argmin(np.abs(r - 3.0))
    print_status(f"  V_RW at r=3M (geometric TEP): {V_RW[idx_3M]:.6f}", "INFO")
    print_status(f"  V_RW at r=3M (Schw ref):      {V_RW_schw[idx_3M]:.6f}", "INFO")
    print_status(f"  Expected (Schw):     0.148148", "INFO")
    print_status(f"  V_RW formula: F * [l(l+1)/r^2 - 6M/r^3], l=2", "DEBUG")
    print_status(f"  V_RW at r=2M (horizon): {V_RW[np.argmin(np.abs(r - 2.0))]:.6f} (-> 0)", "DEBUG")
    print_status(f"  V_RW at r=10M: {V_RW[np.argmin(np.abs(r - 10.0))]:.6f}", "DEBUG")
    print_status(f"  V_RW at r=100M: {V_RW[np.argmin(np.abs(r - 100.0))]:.6f} (-> 0)", "DEBUG")

    # --- Zerilli potential ---
    print_status("Computing Zerilli potential (l=2)...", "PROCESS")
    V_Z, V_Z_schw = compute_zerilli_potential(r, r_star, metric, l=2)
    print_status(f"  V_Z at r=3M (geometric TEP): {V_Z[idx_3M]:.6f}", "DEBUG")
    print_status(f"  V_Z at r=3M (Schw ref):      {V_Z_schw[idx_3M]:.6f}", "DEBUG")
    print_status(f"  V_Z formula: 2F(lam^2(lam+1)r^3 + 3lam^2Mr^2 + 9lamM^2r + 9M^3) / [r^3(lam*r + 3M)^2]", "DEBUG")
    print_status(f"  lam = (l-1)(l+2)/2 = {(2-1)*(2+2)/2.0:.1f} for l=2", "DEBUG")

    # --- QNMs (exterior only) ---
    print_status("Computing quasinormal modes (WKB, n=0, exterior only)...", "PROCESS")
    qnm_rw = compute_qnm_wkb(r, r_star, V_RW, n=0)
    qnm_rw_schw = compute_qnm_wkb(r, r_star, V_RW_schw, n=0)
    qnm_z = compute_qnm_wkb(r, r_star, V_Z, n=0)
    qnm_z_schw = compute_qnm_wkb(r, r_star, V_Z_schw, n=0)
    print_status(f"  RW peak: V_max={qnm_rw['V_max']:.6f} at r={qnm_rw['r_max']:.4f}M, r*={qnm_rw.get('r_star_max', 0):.4f}", "DEBUG")
    print_status(f"  RW V''_0 = {qnm_rw.get('V_pp', 0):.6f} (curvature at peak)", "DEBUG")
    print_status(f"  Zerilli peak: V_max={qnm_z['V_max']:.6f} at r={qnm_z['r_max']:.4f}M", "DEBUG")
    print_status(f"  WKB formula: omega^2 = V_0 + i*sqrt(-2*V''_0)*(n+1/2)", "DEBUG")

    # --- Echo cavity check ---
    print_status("Checking for echo cavity...", "PROCESS")
    # r_temporal_horizon is None under the Gaussian-bump model (globally
    # Lorentzian, no determinant-zero boundary); retained for API compatibility.
    r_t = solution.get("r_temporal_horizon")
    idx_t = solution.get("idx_temporal_horizon")
    echo = check_echo_cavity(r, V_RW, r_star=r_star,
                             r_temporal_horizon=r_t,
                             idx_temporal_horizon=idx_t)

    # --- Diagnostics ---
    print_status("", "INFO")
    print_status("Regge-Wheeler Potential (l=2)", "TITLE")
    print_status(f"  V_max (geometric TEP): {qnm_rw['V_max']:.6f} at r={qnm_rw['r_max']:.4f}M", "INFO")
    print_status(f"  V_max (Schw ref):      {qnm_rw_schw['V_max']:.6f} at r={qnm_rw_schw['r_max']:.4f}M", "INFO")
    print_status(f"  Expected (Schw): V_max=0.148148 at r=3.0M", "INFO")

    print_status("", "INFO")
    print_status("QNM Frequencies (l=2, n=0)", "TITLE")
    if qnm_rw["omega_R"] is not None:
        print_status(f"  Geometric TEP/Schw: omega = {qnm_rw['omega_R']:.6f} + {qnm_rw['omega_I']:.6f}i", "INFO")
    else:
        print_status("  Geometric TEP/Schw: omega = None", "WARNING")
    print_status("  Expected accurate Schw: omega ~ 0.3737 - 0.0890i", "INFO")
    print_status("  Canonical TEP-Schwarzschild gravitational shift: 0", "INFO")

    print_status("", "INFO")
    print_status("Exterior Geometric Potential Check", "TITLE")
    print_status(f"  Peaks (exterior): {echo.get('n_peaks_exterior', 0)}", "INFO")
    print_status(f"  Exterior cavity:  {echo['has_exterior_echo_cavity']}", "INFO")
    print_status("  No matter-sector structure is reported as a GW echo prediction.", "INFO")

    # --- Save CSV ---
    csv_data = np.column_stack([
        r, r_star, V_RW, V_RW_schw, V_Z, V_Z_schw,
    ])
    csv_path = step_csv_path(STEP_ID)
    np.savetxt(
        csv_path, csv_data,
        header="r r_star V_RW_geometric V_RW_schwarzschild_reference V_Z_geometric V_Z_schwarzschild_reference",
        delimiter=" ", comments="# ",
    )
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    # Potential profile summary (V_max, r_peak) for both potentials
    potential_profiles = {
        "regge_wheeler": {
            "l": 2,
            "V_max": qnm_rw["V_max"],
            "r_peak_M": qnm_rw["r_max"],
            "r_star_peak": qnm_rw.get("r_star_max"),
            "V_pp_at_peak": qnm_rw.get("V_pp"),
            "V_at_3M": float(V_RW[idx_3M]),
            "V_at_2M": float(V_RW[np.argmin(np.abs(r - 2.0))]),
            "V_at_10M": float(V_RW[np.argmin(np.abs(r - 10.0))]),
            "formula": "V_RW = F * [l(l+1)/r^2 - 6M/r^3]",
        },
        "zerilli": {
            "l": 2,
            "V_max": qnm_z["V_max"],
            "r_peak_M": qnm_z["r_max"],
            "r_star_peak": qnm_z.get("r_star_max"),
            "V_pp_at_peak": qnm_z.get("V_pp"),
            "V_at_3M": float(V_Z[idx_3M]),
            "V_at_2M": float(V_Z[np.argmin(np.abs(r - 2.0))]),
            "V_at_10M": float(V_Z[np.argmin(np.abs(r - 10.0))]),
            "formula": "V_Z = 2F(lam^2(lam+1)r^3 + 3lam^2Mr^2 + 9lamM^2r + 9M^3) / [r^3(lam*r + 3M)^2]",
        },
    }

    # QNM frequency comparison table
    qnm_comparison_table = {
        "regge_wheeler_l2_n0": {
            "omega_R_tep": qnm_rw["omega_R"],
            "omega_I_tep": qnm_rw["omega_I"],
            "omega_R_schw": qnm_rw_schw["omega_R"],
            "omega_I_schw": qnm_rw_schw["omega_I"],
            "tep_schw_shift": 0.0,
            "expected_accurate_schw": "0.3737 - 0.0890i",
        },
        "zerilli_l2_n0": {
            "omega_R_tep": qnm_z["omega_R"],
            "omega_I_tep": qnm_z["omega_I"],
            "omega_R_schw": qnm_z_schw["omega_R"],
            "omega_I_schw": qnm_z_schw["omega_I"],
            "tep_schw_shift": 0.0,
        },
        "note": "Canonical TEP GWs propagate on geometric g (prescribed Schwarzschild); RW/Zerilli QNMs are identical to Schwarzschild.",
    }

    # Echo analysis result (explicit)
    echo_analysis = {
        "domain": echo.get("domain", "geometric_schwarzschild_exterior"),
        "n_peaks_exterior": echo.get("n_peaks_exterior", 0),
        "has_exterior_echo_cavity": echo.get("has_exterior_echo_cavity", False),
        "gravitational_echo_prediction_supported": echo.get("gravitational_echo_prediction_supported", False),
        "r_temporal_horizon": r_t,
        "note": echo.get("note", "No matter-metric or black-hole-interior structure is used as a GW echo prediction."),
    }

    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "model": model.to_dict(),
        "gravitational_metric": "geometric_g_schwarzschild",
        "matter_metric_used_for_gravitational_predictions": False,
        "n_grid_points": len(r),
        "r_range": [float(r[0]), float(r[-1])],
        "r_star_range": [float(r_star[0]), float(r_star[-1])],
        "potential_profiles": potential_profiles,
        "qnm_regge_wheeler_geometric": qnm_rw,
        "qnm_regge_wheeler_schwarzschild_reference": qnm_rw_schw,
        "qnm_zerilli_geometric": qnm_z,
        "qnm_zerilli_schwarzschild_reference": qnm_z_schw,
        "qnm_comparison_table": qnm_comparison_table,
        "canonical_tep_schwarzschild_qnm_shift": 0.0,
        "exterior_geometric_potential_diagnostic": echo,
        "echo_analysis": echo_analysis,
        "note": "Canonical TEP GWs propagate on geometric g; with g prescribed as Schwarzschild, RW/Zerilli potentials and nonspinning QNMs are Schwarzschild and are not A^-4 matter-metric scaled.",
    }
    summary = finalize_result(
        STEP_ID, summary,
        description="Compute Regge-Wheeler and Zerilli gravitational perturbation potentials and WKB quasinormal-mode frequencies on the prescribed Schwarzschild geometric metric",
        key_result=f"RW/Zerilli QNMs (l=2, n=0) are identical to Schwarzschild (omega_RW={qnm_rw['omega_R']:.4f}{qnm_rw['omega_I']:+.4f}i); no echo cavity (n_peaks={echo.get('n_peaks_exterior', 0)}, no GW echo prediction)",
        model=model.to_dict(),
        dependencies=["step_01_field_equations"],
    )
    json_path = step_json_path(STEP_ID)
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status(f"Step 02 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
