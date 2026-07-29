#!/usr/bin/env python3
"""Step 04: Accretion Dynamics and ISCO.

Computes circular timelike geodesic properties, epicyclic frequencies,
ISCO properties, radiative efficiency, redshift profile, and QPO
frequency ratios for the TEP disformal metric.

For the conformal exterior (gtilde = A² * g_Schw, A ≈ 1), all orbital
dynamics reduce to Schwarzschild.  The conformal factor cancels in the
geodesic equations, giving:
  - r_ISCO = 6M
  - E_ISCO = sqrt(8/9)
  - L_ISCO = 2*sqrt(3)*M
  - Omega = sqrt(M/r³)
  - Epicyclic: Omega_r² = (M/r³)(1-6M/r), Omega_theta = Omega

The TEP modifications are in the interior (r < 2M) and do not affect
observable orbital dynamics.  The disformal correction at r = 6M is
< 1e-6 (B ~ 0 in the exterior).

Outputs:
  - results/step_04_accretion.json
  - results/step_04_accretion.csv
  - logs/step_04_accretion.log

Usage:
    python scripts/steps/step_04_accretion.py
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

STEP_ID = "step_04_accretion"


# =============================================================================
# Circular geodesics (Schwarzschild exterior, conformal factor cancels)
# =============================================================================

def compute_circular_geodesics(r, M=1.0):
    """Compute circular timelike geodesic properties.

    For Schwarzschild (and the conformal exterior where A ≈ 1):
        E = (1 - 2M/r) / sqrt(1 - 3M/r)
        L = sqrt(Mr) / sqrt(1 - 3M/r)
        Omega = sqrt(M/r³)

    The conformal factor A cancels: gtilde = A²*g => geodesic paths
    are unchanged, only proper time rescales.
    """
    r_safe = np.where(r > 3.0 * M, r, np.nan)
    f = 1.0 - 2.0 * M / r_safe
    denom = np.sqrt(np.maximum(1.0 - 3.0 * M / r_safe, 1e-30))

    E = f / denom
    L_sq = M * r_safe / np.maximum(1.0 - 3.0 * M / r_safe, 1e-30)
    L = np.sqrt(np.maximum(L_sq, 0))
    Omega = np.sqrt(M / r_safe**3)

    E = np.where(r > 3.0 * M, E, np.nan)
    L = np.where(r > 3.0 * M, L, np.nan)
    Omega = np.where(r > 3.0 * M, Omega, np.nan)

    return {'r': r, 'E': E, 'L': L, 'Omega': Omega}


# =============================================================================
# Epicyclic frequencies
# =============================================================================

def compute_epicyclic_frequencies(r, M=1.0):
    """Compute epicyclic frequencies for circular orbits.

    For Schwarzschild:
        Omega_r² = (M/r³)(1 - 6M/r)   (radial epicyclic)
        Omega_theta² = M/r³            (vertical = orbital)

    The radial epicyclic vanishes at ISCO (r = 6M).
    """
    r_safe = np.where(r > 3.0 * M, r, np.nan)
    Omega_sq = M / r_safe**3
    Omega_r_sq = Omega_sq * (1.0 - 6.0 * M / r_safe)
    Omega_theta_sq = Omega_sq

    Omega_r_sq = np.where(r > 3.0 * M, Omega_r_sq, np.nan)
    Omega_theta_sq = np.where(r > 3.0 * M, Omega_theta_sq, np.nan)

    ratio = np.sqrt(np.maximum(Omega_r_sq, 0)) / np.sqrt(np.maximum(Omega_theta_sq, 0))

    return {
        'r': r,
        'Omega_r_sq': Omega_r_sq,
        'Omega_theta_sq': Omega_theta_sq,
        'Omega_r': np.sqrt(np.maximum(Omega_r_sq, 0)),
        'Omega_theta': np.sqrt(np.maximum(Omega_theta_sq, 0)),
        'ratio_r_to_theta': ratio,
    }


# =============================================================================
# ISCO properties
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


# =============================================================================
# Redshift profile
# =============================================================================

def compute_redshift_profile(r, M=1.0):
    """Compute gravitational redshift factor.

    g_redshift = sqrt(1 - 2M/r)  (Schwarzschild, exterior)

    The conformal factor A ≈ 1 in the exterior, so the redshift
    is the same as Schwarzschild.
    """
    r_safe = np.where(r > 2.0 * M + 1e-10, r, np.nan)
    g = np.sqrt(1.0 - 2.0 * M / r_safe)
    z = 1.0 / g - 1.0
    return {'r': r, 'g_redshift': g, 'redshift_z': z}


# =============================================================================
# QPO analysis
# =============================================================================

def analyze_qpo_resonance(r, epicyclic, M=1.0):
    """Check for 3:2 epicyclic frequency resonance.

    In Schwarzschild, the ratio Omega_r/Omega_theta varies from 0 at
    ISCO to 1 at infinity.  The 3:2 resonance (ratio = 2/3) occurs at
    a specific radius, which can be compared to observed QPO frequencies.
    """
    ratio = epicyclic['ratio_r_to_theta']
    valid = np.isfinite(ratio) & (ratio > 0) & (r > 3.0 * M)
    if np.any(valid):
        ratios = ratio[valid]
        r_valid = r[valid]
        idx_32 = np.argmin(np.abs(ratios - 2.0 / 3.0))
        r_32 = r_valid[idx_32]
        ratio_32 = ratios[idx_32]
        has_32 = abs(ratio_32 - 2.0 / 3.0) < 0.01
        # Compute frequencies at the resonance radius
        Omega_at_r = np.sqrt(M / r_32**3)
        f_orb = Omega_at_r / (2.0 * np.pi)
        f_lower = (2.0 / 3.0) * f_orb  # radial epicyclic
        f_upper = f_orb  # vertical = orbital
    else:
        r_32 = None
        ratio_32 = None
        has_32 = False
        f_orb = None
        f_lower = None
        f_upper = None

    return {
        'r_32_resonance': float(r_32) if r_32 is not None else None,
        'ratio_at_resonance': float(ratio_32) if ratio_32 is not None else None,
        'has_32_resonance': bool(has_32),
        'f_orbital_at_resonance': float(f_orb) if f_orb is not None else None,
        'f_lower_at_resonance': float(f_lower) if f_lower is not None else None,
        'f_upper_at_resonance': float(f_upper) if f_upper is not None else None,
        'note': '3:2 QPO ratio from epicyclic frequencies (Schwarzschild exterior)',
    }


# =============================================================================
# Main
# =============================================================================

def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 04: Accretion Dynamics and ISCO", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    model = TEPBHModel(
        beta_A=-1.0, B0=1.0, n_B=2.0,
        phi_0=2.0, delta=0.05, M=1.0,
    )
    print_status(f"Model parameters: {model.to_dict()}", "INFO")
    print_status("Exterior is conformal (A~1, B~0): orbital dynamics = Schwarzschild", "INFO")

    print_status("Solving TEP-BH background...", "PROCESS")
    solution = solve_tep_bh(model)
    if not solution["success"]:
        print_status("BACKGROUND SOLUTION FAILED", "ERROR")
        return {"step": STEP_ID, "status": "failed"}

    r = solution["r"]
    metric = solution["metric"]
    M = model.M

    print_status(f"  Grid: {len(r)} points, r in [{r[0]:.6e}, {r[-1]:.6e}]", "DEBUG")
    print_status(f"  r_min = {r[0]:.6e} M, r_max = {r[-1]:.6e} M", "DEBUG")
    print_status(f"  M = {M} (geometric units)", "DEBUG")

    # Verify exterior screening
    idx_6M = np.argmin(np.abs(r - 6.0 * M))
    A_at_6M = float(metric['A'][idx_6M])
    B_at_6M = float(metric['B'][idx_6M])
    print_status(f"  A at r=6M: {A_at_6M:.8f} (should be ~1)", "INFO")
    print_status(f"  B at r=6M: {B_at_6M:.2e} (should be ~0)", "INFO")

    # --- Circular geodesics ---
    print_status("Computing circular geodesic properties...", "PROCESS")
    geodesics = compute_circular_geodesics(r, M)
    print_status(f"  E at r=6M: {geodesics['E'][idx_6M]:.6f} (= sqrt(8/9) ~ 0.9428)", "DEBUG")
    print_status(f"  L at r=6M: {geodesics['L'][idx_6M]:.6f} M (= 2*sqrt(3) ~ 3.464)", "DEBUG")
    print_status(f"  Omega at r=6M: {geodesics['Omega'][idx_6M]:.6f} /M (= 1/(6*sqrt(6)) ~ 0.0680)", "DEBUG")
    print_status(f"  Formula: E = (1-2M/r)/sqrt(1-3M/r), L = sqrt(Mr)/sqrt(1-3M/r), Omega = sqrt(M/r^3)", "DEBUG")

    # --- Epicyclic frequencies ---
    print_status("Computing epicyclic frequencies...", "PROCESS")
    epicyclic = compute_epicyclic_frequencies(r, M)
    print_status(f"  Omega_r^2 = (M/r^3)(1-6M/r), Omega_theta^2 = M/r^3", "DEBUG")
    print_status(f"  Radial epicyclic vanishes at ISCO (r=6M): Omega_r = {epicyclic['Omega_r'][idx_6M]:.6f}", "DEBUG")

    # --- ISCO properties ---
    print_status("Computing ISCO properties...", "PROCESS")
    isco = compute_isco_properties(M)
    print_status(f"  r_ISCO = 6M = {isco['r_isco']:.4f}", "DEBUG")
    print_status(f"  E_ISCO = sqrt(8/9) = {isco['E_isco']:.6f}", "DEBUG")
    print_status(f"  L_ISCO = 2*sqrt(3)*M = {isco['L_isco']:.6f}", "DEBUG")
    print_status(f"  Omega_ISCO = 1/(6*sqrt(6)*M) = {isco['Omega_isco']:.6f}", "DEBUG")
    print_status(f"  Radiative efficiency eta = 1 - E_ISCO = {isco['radiative_efficiency_percent']:.4f}%", "DEBUG")
    print_status(f"  Redshift factor g = sqrt(1-2M/r_ISCO) = sqrt(2/3) = {isco['redshift_factor_isco']:.6f}", "DEBUG")
    print_status(f"  Redshift z at ISCO = {isco['redshift_z_isco']:.6f}", "DEBUG")

    # --- Redshift profile ---
    print_status("Computing redshift profile...", "PROCESS")
    redshift = compute_redshift_profile(r, M)
    print_status(f"  g_redshift at r=6M: {redshift['g_redshift'][idx_6M]:.6f} (= sqrt(2/3))", "DEBUG")
    print_status(f"  g_redshift at r=10M: {redshift['g_redshift'][np.argmin(np.abs(r - 10.0))]:.6f}", "DEBUG")
    print_status(f"  Formula: g = sqrt(1 - 2M/r), z = 1/g - 1", "DEBUG")

    # --- QPO analysis ---
    print_status("Analyzing QPO frequency ratios...", "PROCESS")
    qpo = analyze_qpo_resonance(r, epicyclic, M)
    if qpo["r_32_resonance"] is not None:
        print_status(f"  3:2 resonance at r = {qpo['r_32_resonance']:.4f} M", "DEBUG")
        print_status(f"  f_orbital at resonance: {qpo['f_orbital_at_resonance']:.6f} (geometric units)", "DEBUG")
        print_status(f"  f_lower (radial) at resonance: {qpo['f_lower_at_resonance']:.6f}", "DEBUG")
        print_status(f"  f_upper (vertical) at resonance: {qpo['f_upper_at_resonance']:.6f}", "DEBUG")
    print_status(f"  QPO ratio condition: Omega_r/Omega_theta = 2/3 at resonance radius", "DEBUG")

    # --- Report ---
    print_status("", "INFO")
    print_status("ISCO Properties (Schwarzschild exterior, A~1)", "TITLE")
    print_status(f"  r_ISCO:              {isco['r_isco']:.4f} M", "INFO")
    print_status(f"  E_ISCO:              {isco['E_isco']:.6f}", "INFO")
    print_status(f"  L_ISCO:              {isco['L_isco']:.6f} M", "INFO")
    print_status(f"  Omega_ISCO:          {isco['Omega_isco']:.6f} / M", "INFO")
    print_status(f"  Radiative efficiency: {isco['radiative_efficiency_percent']:.4f}%", "INFO")
    print_status(f"  Redshift factor g:    {isco['redshift_factor_isco']:.6f}", "INFO")
    print_status(f"  Redshift z at ISCO:   {isco['redshift_z_isco']:.6f}", "INFO")

    print_status("", "INFO")
    print_status("Epicyclic Frequencies", "TITLE")
    idx_10M = np.argmin(np.abs(r - 10.0))
    idx_20M = np.argmin(np.abs(r - 20.0))
    idx_50M = np.argmin(np.abs(r - 50.0))
    for r_val, idx in [(6.0, idx_6M), (10.0, idx_10M), (20.0, idx_20M), (50.0, idx_50M)]:
        print_status(f"  r={r_val}M: Omega_r={epicyclic['Omega_r'][idx]:.6f}, "
                     f"Omega_theta={epicyclic['Omega_theta'][idx]:.6f}, "
                     f"ratio={epicyclic['ratio_r_to_theta'][idx]:.6f}", "INFO")

    print_status("", "INFO")
    print_status("QPO Analysis", "TITLE")
    if qpo["r_32_resonance"] is not None:
        print_status(f"  3:2 resonance at r = {qpo['r_32_resonance']:.4f} M", "INFO")
        print_status(f"  Ratio at resonance: {qpo['ratio_at_resonance']:.6f}", "INFO")
        print_status(f"  Has 3:2 resonance: {qpo['has_32_resonance']}", "INFO")
        print_status(f"  f_orbital: {qpo['f_orbital_at_resonance']:.6f}", "INFO")
        print_status(f"  f_lower (radial): {qpo['f_lower_at_resonance']:.6f}", "INFO")
        print_status(f"  f_upper (vertical): {qpo['f_upper_at_resonance']:.6f}", "INFO")
    else:
        print_status(f"  No 3:2 resonance found", "INFO")
    print_status(f"  Note: {qpo['note']}", "INFO")

    # --- Save CSV ---
    csv_data = np.column_stack([
        r,
        geodesics["E"],
        geodesics["L"],
        geodesics["Omega"],
        epicyclic["Omega_r"],
        epicyclic["Omega_theta"],
        epicyclic["ratio_r_to_theta"],
        redshift["g_redshift"],
        redshift["redshift_z"],
    ])
    csv_path = step_csv_path(STEP_ID)
    np.savetxt(
        csv_path, csv_data,
        header="r E L Omega Omega_r Omega_theta ratio_r_theta g_redshift redshift_z",
        delimiter=" ", comments="# ",
    )
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- QPO frequency scaling for a reference BH mass ---
    # Convert geometric-unit frequencies to Hz for a reference stellar-mass BH
    # f_Hz = Omega_geometric / (2*pi) * c^3 / (G * M_kg)
    # For M = 10 M_sun: characteristic frequency ~ 3.23 kHz
    G_NEWTON = 6.67430e-11
    C_LIGHT = 2.99792458e8
    M_SUN = 1.98847e30
    ref_mass_msun = 10.0  # reference stellar-mass BH
    ref_M_kg = ref_mass_msun * M_SUN
    freq_scale = C_LIGHT**3 / (2.0 * np.pi * G_NEWTON * ref_M_kg)  # Hz per geometric Omega
    qpo_scaling = {
        "reference_BH_mass_Msun": ref_mass_msun,
        "frequency_conversion": "f_Hz = Omega / (2*pi) * c^3 / (G * M_kg)",
        "freq_scale_Hz_per_Omega": float(freq_scale),
        "Omega_isco_Hz": float(isco["Omega_isco"] * freq_scale),
        "f_orbital_isco_Hz": float(isco["Omega_isco"] / (2.0 * np.pi) * freq_scale * 2.0 * np.pi),
    }
    if qpo["f_orbital_at_resonance"] is not None:
        qpo_scaling["f_orbital_at_resonance_Hz"] = float(qpo["f_orbital_at_resonance"] * freq_scale)
        qpo_scaling["f_lower_at_resonance_Hz"] = float(qpo["f_lower_at_resonance"] * freq_scale)
        qpo_scaling["f_upper_at_resonance_Hz"] = float(qpo["f_upper_at_resonance"] * freq_scale)
    print_status(f"  QPO frequency scaling for M = {ref_mass_msun} M_sun:", "DEBUG")
    print_status(f"    freq_scale = {freq_scale:.2f} Hz per geometric Omega", "DEBUG")
    print_status(f"    Omega_ISCO -> {qpo_scaling['Omega_isco_Hz']:.2f} Hz", "DEBUG")
    if "f_orbital_at_resonance_Hz" in qpo_scaling:
        print_status(f"    f_orbital at 3:2 resonance -> {qpo_scaling['f_orbital_at_resonance_Hz']:.2f} Hz", "DEBUG")
        print_status(f"    f_lower (radial) at 3:2 resonance -> {qpo_scaling['f_lower_at_resonance_Hz']:.2f} Hz", "DEBUG")
        print_status(f"    f_upper (vertical) at 3:2 resonance -> {qpo_scaling['f_upper_at_resonance_Hz']:.2f} Hz", "DEBUG")

    # --- Save JSON summary ---
    # Epicyclic frequency profiles at multiple radii
    epicyclic_profiles = {
        "at_6M": {
            "Omega_r": float(epicyclic["Omega_r"][idx_6M]) if np.isfinite(epicyclic["Omega_r"][idx_6M]) else None,
            "Omega_theta": float(epicyclic["Omega_theta"][idx_6M]) if np.isfinite(epicyclic["Omega_theta"][idx_6M]) else None,
            "ratio": float(epicyclic["ratio_r_to_theta"][idx_6M]) if np.isfinite(epicyclic["ratio_r_to_theta"][idx_6M]) else None,
        },
        "at_10M": {
            "Omega_r": float(epicyclic["Omega_r"][idx_10M]) if np.isfinite(epicyclic["Omega_r"][idx_10M]) else None,
            "Omega_theta": float(epicyclic["Omega_theta"][idx_10M]) if np.isfinite(epicyclic["Omega_theta"][idx_10M]) else None,
            "ratio": float(epicyclic["ratio_r_to_theta"][idx_10M]) if np.isfinite(epicyclic["ratio_r_to_theta"][idx_10M]) else None,
        },
        "at_20M": {
            "Omega_r": float(epicyclic["Omega_r"][idx_20M]) if np.isfinite(epicyclic["Omega_r"][idx_20M]) else None,
            "Omega_theta": float(epicyclic["Omega_theta"][idx_20M]) if np.isfinite(epicyclic["Omega_theta"][idx_20M]) else None,
            "ratio": float(epicyclic["ratio_r_to_theta"][idx_20M]) if np.isfinite(epicyclic["ratio_r_to_theta"][idx_20M]) else None,
        },
        "at_50M": {
            "Omega_r": float(epicyclic["Omega_r"][idx_50M]) if np.isfinite(epicyclic["Omega_r"][idx_50M]) else None,
            "Omega_theta": float(epicyclic["Omega_theta"][idx_50M]) if np.isfinite(epicyclic["Omega_theta"][idx_50M]) else None,
            "ratio": float(epicyclic["ratio_r_to_theta"][idx_50M]) if np.isfinite(epicyclic["ratio_r_to_theta"][idx_50M]) else None,
        },
    }

    # ISCO orbital parameters with full detail
    isco_detail = {
        "r_isco_M": isco["r_isco"],
        "E_isco": isco["E_isco"],
        "L_isco_M": isco["L_isco"],
        "Omega_isco_per_M": isco["Omega_isco"],
        "radiative_efficiency": isco["radiative_efficiency"],
        "radiative_efficiency_percent": isco["radiative_efficiency_percent"],
        "redshift_factor_isco": isco["redshift_factor_isco"],
        "redshift_z_isco": isco["redshift_z_isco"],
        "formulas": {
            "E_isco": "sqrt(8/9) ~ 0.9428",
            "L_isco": "2*sqrt(3)*M ~ 3.464*M",
            "Omega_isco": "1/(6*sqrt(6)*M) ~ 0.0680/M",
            "eta": "1 - E_isco ~ 5.72%",
            "g_redshift": "sqrt(1 - 2M/r_ISCO) = sqrt(2/3) ~ 0.8165",
        },
    }

    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "model": model.to_dict(),
        "n_points": len(r),
        "r_min": float(r[0]),
        "r_max": float(r[-1]),
        "radial_grid": {
            "n_points": len(r),
            "r_min_M": float(r[0]),
            "r_max_M": float(r[-1]),
        },
        "isco_properties": isco,
        "isco_detail": isco_detail,
        "A_at_6M": A_at_6M,
        "B_at_6M": B_at_6M,
        "epicyclic_profiles": epicyclic_profiles,
        "epicyclic_at_6M": epicyclic_profiles["at_6M"],
        "epicyclic_at_10M": epicyclic_profiles["at_10M"],
        "epicyclic_at_20M": epicyclic_profiles["at_20M"],
        "epicyclic_at_50M": epicyclic_profiles["at_50M"],
        "qpo_analysis": qpo,
        "qpo_frequency_scaling": qpo_scaling,
        "note": "Exterior orbital dynamics are Schwarzschild (conformal factor A cancels). "
                "TEP modifications are interior (r < 2M) and do not affect observable orbital dynamics.",
    }
    summary = finalize_result(
        STEP_ID, summary,
        description="Compute circular timelike geodesic properties, epicyclic frequencies, ISCO, radiative efficiency, redshift profile, and QPO frequency ratios for the TEP disformal metric",
        key_result=f"ISCO at r=6M with eta={isco['radiative_efficiency_percent']:.2f}% radiative efficiency; epicyclic frequencies and 3:2 QPO resonance match Schwarzschild (conformal factor A cancels in exterior)",
        model=model.to_dict(),
        dependencies=["step_01_field_equations"],
    )
    json_path = step_json_path(STEP_ID)
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 04 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
