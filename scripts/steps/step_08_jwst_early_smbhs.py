#!/usr/bin/env python3
"""Step 08: JWST Early SMBHs — Eddington-limited Growth Analysis.

JWST and ground-based surveys have identified massive black holes at high
redshift, when the Universe was <1 Gyr old. This step compares their masses
with an idealized continuous-Eddington growth envelope from 100 M_sun seeds
using the standard 45 Myr Salpeter e-folding time.

TEP does not modify exterior mass measurements. The conformal factor
A(phi) ~ 1 in the exterior region where BH masses are measured
(virial / reverberation / dynamical), so observed masses are true masses.
The JWST early SMBH problem is not addressed by TEP through mass inflation.

Data sources:
  - JWST CEERS survey: Harikane et al. 2023, ApJ 958, 11
  - JWST JADES survey: Greene et al. 2023 (CEERS-AGN z~7)
  - Bogdan et al. 2024, Nature Astronomy (z=10.1 BH ~10^9 M_sun)
  - Wu et al. 2015 (high-z quasar masses)
  - Inayoshi, Haiman, Ostriker 2020 (Eddington-limited growth review)

Outputs:
  - data/processed/jwst_highz_bh_masses.json
  - results/step_08_jwst_early_smbhs.json
  - results/step_08_jwst_early_smbhs.csv
  - logs/step_08_jwst_early_smbhs.log
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

STEP_ID = "step_08_jwst_early_smbhs"

# Cosmological parameters (Planck 2018)
H0 = 67.4  # km/s/Mpc
OMEGA_M = 0.315
OMEGA_L = 0.685

# Eddington-limited growth parameters
ETA_RAD = 0.1  # Radiative efficiency
EPSILON = 0.1  # Accretion efficiency
T_EDD = 45e6  # Salpeter e-folding time in years (epsilon=0.1)


def age_at_redshift(z: float) -> float:
    """Age of universe at redshift z in years."""
    from math import sqrt, log
    # Flat LCDM: t(z) = (2/(3*H0*sqrt(Omega_L))) * asinh(sqrt(Omega_L/Omega_M) / (1+z)^(3/2))
    H0_s = H0 * 1000 / 3.08567758128e22  # H0 in s^-1 (Mpc in meters)
    t_H = 1.0 / H0_s  # Hubble time in seconds
    arg = sqrt(OMEGA_L / OMEGA_M) / (1 + z)**1.5
    t_sec = (2.0 / (3.0 * H0_s * sqrt(OMEGA_L))) * np.arcsinh(arg)
    return t_sec / (3.1557e7)  # Convert to years


def eddington_limited_growth(M_seed: float, t_available_yr: float,
                              f_edd: float = 1.0) -> float:
    """Eddington-limited growth from seed mass.

    M(t) = M_seed * exp(t / t_salpeter * f_edd)

    Parameters
    ----------
    M_seed : float
        Seed black hole mass in solar masses.
    t_available_yr : float
        Available time in years.
    f_edd : float
        Eddington ratio (1.0 = Eddington-limited).

    Returns
    -------
    float
        Final mass in solar masses.
    """
    return M_seed * np.exp(t_available_yr / T_EDD * f_edd)


def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 08: JWST Early SMBHs — Eddington-limited Growth Analysis", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status("")

    # --- Compiled JWST high-z BH measurements ---
    highz_bhs = {
        "CEERS_2782": {
            "redshift": 5.242,
            "mass_msun": {"value": 1.3e7, "uncertainty": 0.4e7},
            "citation": "Kocevski et al. 2023, ApJL 946, L14",
            "survey": "JWST CEERS",
        },
        "CEERS_746": {
            "redshift": 5.624,
            "mass_msun": {"value": 2.8e7, "uncertainty": 1.9e7},
            "citation": "Kocevski et al. 2023, ApJL 946, L14",
            "survey": "JWST CEERS",
        },
        "GN_z11": {
            "redshift": 10.6,
            "mass_msun": {"value": 1.6e6, "uncertainty": 0.8e6},
            "citation": "Maiolino et al. 2024, Nature 627, 59",
            "survey": "JWST",
        },
        "UHZ1": {
            "redshift": 10.1,
            "mass_msun": {"value": 4e7, "uncertainty": 1e7},
            "citation": "Bogdan et al. 2024, ApJL 960, L1",
            "survey": "Chandra+JWST",
        },
        "J1243_0100": {
            "redshift": 7.07,
            "mass_msun": {"value": 3.3e8, "uncertainty": 2.0e8},
            "citation": "Matsuoka et al. 2019, ApJL 872, L2",
            "survey": "Subaru HSC",
        },
        "J1342_0928": {
            "redshift": 7.54,
            "mass_msun": {"value": 8e8, "uncertainty": 2e8},
            "citation": "Banados et al. 2018, Nature 553, 473",
            "survey": "SDSS",
        },
    }

    # Save compiled data
    highz_path = PROCESSED_DIR / "jwst_highz_bh_masses.json"
    write_json(highz_path, highz_bhs)
    print_status(f"High-z BH masses saved to {rel(highz_path)}", "SUCCESS")

    # --- Standard Eddington-limited growth analysis ---
    print_status("", "INFO")
    print_status("Standard Eddington-Limited Growth", "TITLE")
    print_status(f"  Seed mass: 100 M_sun (stellar remnant)", "INFO")
    print_status(f"  Salpeter time: {T_EDD/1e6:.0f} Myr (epsilon=0.1)", "INFO")
    print_status(f"  Assumption: continuous Eddington-limited accretion", "INFO")
    print_status(f"  Radiative efficiency eta_rad = {ETA_RAD}", "DEBUG")
    print_status(f"  Accretion efficiency epsilon = {EPSILON}", "DEBUG")
    print_status(f"  Salpeter e-folding time: t_Salp = epsilon * M_Edd / L_Edd = {T_EDD/1e6:.0f} Myr", "DEBUG")
    print_status(f"  Growth formula: M(t) = M_seed * exp(t / t_Salp * f_edd)", "DEBUG")
    print_status(f"  Cosmology: H0={H0}, Omega_m={OMEGA_M}, Omega_L={OMEGA_L} (Planck 2018)", "DEBUG")
    print_status("")

    csv_rows = []
    M_seed = 100.0  # Solar masses

    print_status(f"  {'Source':<20} {'z':<8} {'Age(Gyr)':<10} {'M_obs(M_sun)':<15} {'M_Edd(M_sun)':<15} {'Ratio':<10} {'t_need(Myr)':<12}", "INFO")
    print_status(f"  {'-'*20} {'-'*8} {'-'*10} {'-'*15} {'-'*15} {'-'*10} {'-'*12}", "INFO")

    for name, bh in highz_bhs.items():
        z = bh["redshift"]
        M_obs = bh["mass_msun"]["value"]
        M_obs_err = bh["mass_msun"]["uncertainty"]
        age_yr = age_at_redshift(z)
        age_gyr = age_yr / 1e9

        # Eddington-limited growth (f_edd = 1)
        M_edd = eddington_limited_growth(M_seed, age_yr, f_edd=1.0)

        # How many e-folds needed?
        n_efolds = np.log(M_obs / M_seed)
        t_needed = n_efolds * T_EDD
        ratio = M_obs / M_edd

        print_status(f"  {name:<20} {z:<8.2f} {age_gyr:<10.3f} {M_obs:<15.2e} {M_edd:<15.2e} {ratio:<10.1f} {t_needed/1e6:<12.0f}", "INFO")
        print_status(f"    {name}: z={z:.3f}, age={age_gyr:.4f} Gyr, M_obs={M_obs:.2e} ± {M_obs_err:.2e} M_sun", "DEBUG")
        print_status(f"    {name}: M_Edd(100 M_sun seed, f_edd=1) = {M_edd:.2e} M_sun", "DEBUG")
        print_status(f"    {name}: e-folds needed = ln(M_obs/M_seed) = {n_efolds:.2f}, t_needed = {t_needed/1e6:.0f} Myr", "DEBUG")
        print_status(f"    {name}: M_obs/M_Edd = {ratio:.2f}x ({'EXCEEDS' if ratio > 1 else 'within'} Eddington envelope)", "DEBUG")
        print_status(f"    {name}: citation = {bh['citation']}", "DEBUG")

        csv_rows.append({
            "source": name,
            "redshift": z,
            "age_gyr": age_gyr,
            "mass_observed_msun": M_obs,
            "mass_observed_uncertainty": M_obs_err,
            "mass_eddington_limited_msun": M_edd,
            "ratio_obs_to_edd": ratio,
            "efolds_needed": n_efolds,
            "time_needed_myr": t_needed / 1e6,
            "citation": bh["citation"],
        })

    print_status("", "INFO")
    n_exceed = sum(row["ratio_obs_to_edd"] > 1.0 for row in csv_rows)
    print_status(
        f"Key finding: {n_exceed}/{len(csv_rows)} objects exceed the idealized ",
        "INFO",
    )
    print_status(
        "100-solar-mass, continuous-Eddington growth envelope when the standard ",
        "INFO",
    )
    print_status("45 Myr Salpeter e-folding time is used.", "INFO")

    # --- TEP assessment ---
    print_status("", "INFO")
    print_status("TEP Assessment", "TITLE")
    print_status("  TEP does not modify exterior mass measurements.", "INFO")
    print_status("  The conformal factor A(phi) ~ 1 in the exterior region where BH", "INFO")
    print_status("  masses are measured (virial / reverberation / dynamical methods),", "INFO")
    print_status("  so observed masses are true masses.", "INFO")
    print_status("")
    print_status("  The JWST early SMBH problem is not addressed by TEP through mass", "INFO")
    print_status("  inflation. The observed masses are real, and the challenge of", "INFO")
    print_status("  growing billion-solar-mass BHs in <1 Gyr remains an open", "INFO")
    print_status("  astrophysical problem that TEP does not resolve.", "INFO")

    # --- Compute statistics ---
    ratios = [r["ratio_obs_to_edd"] for r in csv_rows]
    print_status("", "INFO")
    print_status(f"  Mean M_obs/M_edd ratio: {np.mean(ratios):.1f}x", "INFO")
    print_status(f"  Min ratio: {np.min(ratios):.1f}x", "INFO")
    print_status(f"  Max ratio: {np.max(ratios):.1f}x", "INFO")
    print_status(f"  All ratios: {[f'{r:.2f}' for r in ratios]}", "DEBUG")
    print_status(f"  Sources exceeding Eddington envelope: {n_exceed}/{len(csv_rows)}", "DEBUG")
    print_status(f"  Salpeter timescale: {T_EDD/1e6:.0f} Myr per e-fold (epsilon=0.1, f_edd=1)", "DEBUG")
    print_status(f"  TEP assessment: A ≈ 1 in exterior → observed masses are true masses", "SUCCESS")

    # --- Save CSV ---
    csv_path = step_csv_path(STEP_ID)
    write_csv(csv_path, csv_rows)
    print_status(f"CSV saved to {rel(csv_path)}", "SUCCESS")

    # --- Save JSON summary ---
    summary = {
        "step": STEP_ID,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "data_sources": {name: bh["citation"] for name, bh in highz_bhs.items()},
        "eddington_parameters": {
            "seed_mass_msun": M_seed,
            "salpeter_time_myr": T_EDD / 1e6,
            "radiative_efficiency": ETA_RAD,
        },
        "cosmology": {
            "H0": H0,
            "Omega_m": OMEGA_M,
            "Omega_L": OMEGA_L,
        },
        "results": {
            "n_sources": len(highz_bhs),
            "n_exceed_continuous_eddington_envelope": n_exceed,
            "mean_ratio_obs_to_edd": float(np.mean(ratios)),
            "min_ratio_obs_to_edd": float(np.min(ratios)),
            "max_ratio_obs_to_edd": float(np.max(ratios)),
        },
        "tep_assessment": (
            "TEP does not modify exterior mass measurements. The conformal factor "
            "A(phi) ~ 1 in the exterior where BH masses are measured, so observed "
            "masses are true masses. This calculation uses the standard 45 Myr "
            "Salpeter e-folding time and reports which objects exceed an idealized "
            "continuous-Eddington envelope from a 100-solar-mass seed; it does not "
            "attribute any residual growth tension to TEP."
        ),
    }
    json_path = step_json_path(STEP_ID)
    summary = finalize_result(
        STEP_ID, summary,
        description="Compare JWST high-redshift SMBH masses against an idealized Eddington-limited growth envelope from 100 M_sun seeds using the 45 Myr Salpeter timescale",
        key_result=f"{n_exceed}/{len(csv_rows)} high-z SMBHs exceed the continuous-Eddington growth envelope (mean M_obs/M_Edd = {np.mean(ratios):.1f}x); TEP does not modify exterior mass measurements",
        dependencies=["step_00_data_download"],
    )
    write_json(json_path, summary)
    print_status(f"JSON summary saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 08 complete.", "SUCCESS")
    return summary


if __name__ == "__main__":
    main()
