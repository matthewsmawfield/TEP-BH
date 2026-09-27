#!/usr/bin/env python3
"""Step 00 (Inference): Download and compile S-star astrometry and spectroscopy.

THE PRIMARY FALSIFIABLE TEST of TEP at black-hole scales is a non-isochronous
refit of the raw S-star observations around Sgr A*. This step acquires the data.

DATA SOURCES:
  0. Gillessen et al. (2017): 25 years of S-star monitoring in the Galactic Center
     - Astrometry: 145 NACO/NTT/Keck/Gemini epochs for S2 relative to SgrA* (1992-2016)
     - Spectroscopy: 44 SINFONI/Keck/Gemini radial velocity epochs for S2
     - VizieR catalogue J/ApJ/837/30, table5.dat (machine-readable ASCII)
     - Bibcode: 2017ApJ...837...30G

  1. GRAVITY Collaboration (2018a): S2 pericentre passage, gravitational redshift
     - Astrometry: NACO (1992-2017) + GRAVITY (2017-2018)
     - Spectroscopy: SINFONI radial velocities
     - DOI: 10.1051/0004-6361/201833718

  2. GRAVITY Collaboration (2020): Schwarzschild precession detection
     - Astrometry: NACO (118 points) + GRAVITY (2017-2019) + flare positions (75)
     - Spectroscopy: SINFONI + Keck
     - DOI: 10.1051/0004-6361/202037813

  3. GRAVITY Collaboration (2022): Mass distribution, multiple stellar orbits
     - Astrometry: S2, S29, S38, S55 (NACO + GRAVITY)
     - Spectroscopy: SINFONI + ERIS
     - DOI: 10.1051/0004-6361/202142465

  4. GRAVITY Collaboration (2024): Extended mass constraints
     - Updated data through 2022.7
     - DOI: 10.1051/0004-6361/202452274

PUBLISHED ORBITAL PARAMETERS (S2):
  M_BH = 4.297 × 10^6 M_sun  (±0.013)
  R0 = 8277 pc  (±33)
  T_orb = 16.05 yr
  a = 0.1251 arcsec  (= 970 AU at R0)
  e = 0.8843
  i = 134.9°
  omega = 66.0°  (argument of perihelion)
  Omega = 226.0°  (longitude of ascending node)
  t_peri = 2018.37  (pericentre epoch)

DATA AVAILABILITY:
  The raw astrometric and spectroscopic data tables are published in the
  papers' appendices and electronic supplements. VizieR catalogues may
  contain machine-readable versions. This step attempts to download from
  public sources; if unavailable, it compiles the published orbital
  parameters and generates a faithful synthetic dataset for pipeline
  testing, clearly marked as such.

Outputs:
  data/raw/sstar/                    (downloaded raw data)
  data/processed/sstar_astrometry.json
  data/processed/sstar_spectroscopy.json
  data/processed/sstar_orbital_params.json
  results/step_43_sstar_data.json
  logs/step_43_sstar_data.log
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np

# Physical constants
G_SI = 6.674e-11  # m^3 kg^-1 s^-2
c_si = 2.998e8    # m/s
M_sun_kg = 1.989e30  # kg
AU_m = 1.496e11  # m

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    write_json,
    finalize_result,
    tep_s2_conformal_factor,
    RAW_DIR,
    PROCESSED_DIR,
    RESULTS_DIR,
)

STEP_ID = "step_43_sstar_data"
TIMEOUT_SECONDS = 60

# ============================================================
# Published S2 orbital parameters (GRAVITY Collaboration 2022)
# ============================================================
S2_ORBITAL_PARAMS = {
    "star": "S2",
    "M_BH_Msun": 4.297e6,
    "M_BH_err": 0.013e6,
    "R0_pc": 8277.0,
    "R0_err": 33.0,
    "T_orb_yr": 16.05,
    "T_orb_err": 0.18,
    "a_arcsec": 0.1251,
    "a_err": 0.0009,
    "eccentricity": 0.8843,
    "e_err": 0.0036,
    "inclination_deg": 134.9,
    "i_err": 0.4,
    "omega_deg": 66.0,  # argument of perihelion
    "omega_err": 0.7,
    "Omega_deg": 226.0,  # longitude of ascending node
    "Omega_err": 1.0,
    "t_peri_yr": 2018.37,
    "t_peri_err": 0.01,
    "references": [
        "GRAVITY Collaboration 2018a, A&A 615, L15 (DOI: 10.1051/0004-6361/201833718)",
        "GRAVITY Collaboration 2020, A&A 636, L5 (DOI: 10.1051/0004-6361/202037813)",
        "GRAVITY Collaboration 2022, A&A 657, L12 (DOI: 10.1051/0004-6361/202142465)",
    ],
}

# Other S-stars with published orbital parameters
OTHER_STARS = {
    "S29": {
        "T_orb_yr": 83.0,
        "a_arcsec": 0.293,
        "eccentricity": 0.75,
        "inclination_deg": 130.0,
        "t_peri_yr": 2021.4,
    },
    "S38": {
        "T_orb_yr": 19.2,
        "a_arcsec": 0.144,
        "eccentricity": 0.42,
        "inclination_deg": 168.0,
        "t_peri_yr": 2022.1,
    },
    "S55": {
        "T_orb_yr": 12.0,
        "a_arcsec": 0.107,
        "eccentricity": 0.73,
        "inclination_deg": 128.0,
        "t_peri_yr": 2021.7,
    },
}


# ============================================================
# Synthetic data generation (for pipeline testing)
# ============================================================
def _tep_conformal_factor(r_arcsec, M_BH_Msun, R0_pc, eta):
    """Compute the TEP conformal factor A(r) for the mass-inflation branch.

    The mass-inflation branch corresponds to a negative sGB scalar charge
    (alpha_GB < 0), giving A = e^{-phi} = exp(+eta/(3*r_rs)) > 1 for eta > 0.
    Spatial scales are magnified relative to GR, producing positive Phantom Mass.
    """
    R_s_m = 2 * G_SI * M_BH_Msun * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m
    R_s_arcsec = R_s_AU / R0_pc
    r_rs = np.maximum(r_arcsec / R_s_arcsec, 1.0)
    return np.exp(eta / (3.0 * r_rs))


def generate_s2_synthetic_astrometry(params, n_points=150, noise_mas=0.1, eta_inject=0.0):
    """Generate synthetic S2 astrometry from published orbital parameters.

    The astrometry consists of (t, x, y, sigma_x, sigma_y) where:
    - t is in years
    - x, y are angular offsets from Sgr A* in arcseconds
    - sigma is the measurement uncertainty in arcseconds

    The orbit is a Keplerian ellipse with GR Schwarzschild precession.

    If eta_inject > 0, the TEP conformal factor A(r) is applied to simulate
    what would be observed if TEP is true. This creates a dataset with a
    known TEP signal that the pipeline should recover.
    """
    M = params["M_BH_Msun"]  # solar masses
    R0 = params["R0_pc"]  # parsecs
    a = params["a_arcsec"]  # arcseconds
    e = params["eccentricity"]
    i = np.radians(params["inclination_deg"])
    omega = np.radians(params["omega_deg"])
    Omega = np.radians(params["Omega_deg"])
    T = params["T_orb_yr"]
    t_peri = params["t_peri_yr"]

    # Convert noise from mas to arcsec
    noise_arcsec = noise_mas / 1000.0

    # Time span: 1992 to 2024
    t_arr = np.linspace(1992.0, 2024.0, n_points)

    # Mean anomaly
    n = 2 * np.pi / T  # mean motion (rad/yr)
    M_anom = n * (t_arr - t_peri)

    # Solve Kepler's equation for eccentric anomaly E
    E = M_anom.copy()
    for _ in range(100):
        E = M_anom + e * np.sin(E)

    # True anomaly
    cos_nu = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_nu = np.sqrt(1 - e**2) * np.sin(E) / (1 - e * np.cos(E))
    nu = np.arctan2(sin_nu, cos_nu)

    # Orbital radius in arcseconds
    r = a * (1 - e**2) / (1 + e * np.cos(nu))

    # Schwarzschild precession rate (rad per year) — same formula as step_01
    a_AU = a * R0
    a_m = a_AU * AU_m
    d_omega_per_orbit = 6 * np.pi * G_SI * M * M_sun_kg / (c_si**2 * a_m * (1 - e**2))
    d_omega_per_yr = d_omega_per_orbit / T
    omega_precessed = omega + d_omega_per_yr * (t_arr - t_peri)

    # Position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Rotate by precessed argument of perihelion
    x_per_p = x_orb * np.cos(omega_precessed) - y_orb * np.sin(omega_precessed)
    y_per_p = x_orb * np.sin(omega_precessed) + y_orb * np.cos(omega_precessed)

    # Apply inclination (project onto sky)
    x_sky_p = x_per_p
    y_sky_p = y_per_p * np.cos(i)

    # Rotate by longitude of ascending node
    x_sky_rot_p = x_sky_p * np.cos(Omega) - y_sky_p * np.sin(Omega)
    y_sky_rot_p = x_sky_p * np.sin(Omega) + y_sky_p * np.cos(Omega)

    # Apply TEP conformal factor if eta_inject != 0
    # This simulates what would be observed if TEP is true:
    # the observed positions are scaled by A(r) relative to the local orbit
    if eta_inject != 0.0:
        r_obs = np.sqrt(x_sky_rot_p**2 + y_sky_rot_p**2)
        A = _tep_conformal_factor(r_obs, M, R0, eta_inject)
        x_sky_rot_p = A * x_sky_rot_p
        y_sky_rot_p = A * y_sky_rot_p

    # Add noise
    rng = np.random.default_rng(42)
    x_obs = x_sky_rot_p + rng.normal(0, noise_arcsec, n_points)
    y_obs = y_sky_rot_p + rng.normal(0, noise_arcsec, n_points)

    # Uncertainties (GRAVITY: ~50-100 μas, NACO: ~1-2 mas)
    # Early epochs (NACO) have larger errors
    sigma = np.where(t_arr < 2017.0, 0.002, noise_arcsec)  # 2 mas for NACO, noise for GRAVITY

    astrometry = []
    for k in range(n_points):
        astrometry.append({
            "t_yr": float(t_arr[k]),
            "x_arcsec": float(x_obs[k]),
            "y_arcsec": float(y_obs[k]),
            "sigma_x_arcsec": float(sigma[k]),
            "sigma_y_arcsec": float(sigma[k]),
            "instrument": "NACO" if t_arr[k] < 2017.0 else "GRAVITY",
        })

    return astrometry


def generate_s2_synthetic_spectroscopy(params, n_points=100, noise_kms=10, eta_inject=0.0):
    """Generate synthetic S2 radial velocities from published orbital parameters.

    The spectroscopy consists of (t, v_rad, sigma_v) where:
    - t is in years
    - v_rad is the radial velocity in km/s
    - sigma_v is the measurement uncertainty in km/s

    If eta_inject > 0, the TEP conformal factor is applied to simulate
    the additional gravitational redshift from the temporal well.
    """
    M = params["M_BH_Msun"]
    R0 = params["R0_pc"]
    a = params["a_arcsec"]
    e = params["eccentricity"]
    i = np.radians(params["inclination_deg"])
    omega = np.radians(params["omega_deg"])
    T = params["T_orb_yr"]
    t_peri = params["t_peri_yr"]

    # Time span
    t_arr = np.linspace(1992.0, 2024.0, n_points)

    # Mean anomaly
    n_motion = 2 * np.pi / T
    M_anom = n_motion * (t_arr - t_peri)

    # Solve Kepler's equation
    E = M_anom.copy()
    for _ in range(100):
        E = M_anom + e * np.sin(E)

    # True anomaly
    cos_nu = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_nu = np.sqrt(1 - e**2) * np.sin(E) / (1 - e * np.cos(E))
    nu = np.arctan2(sin_nu, cos_nu)

    # Schwarzschild precession (same as astrometry)
    a_AU = a * R0
    a_m = a_AU * AU_m
    d_omega_per_orbit = 6 * np.pi * G_SI * M * M_sun_kg / (c_si**2 * a_m * (1 - e**2))
    d_omega_per_yr = d_omega_per_orbit / T
    omega_precessed = omega + d_omega_per_yr * (t_arr - t_peri)

    # Radial velocity (simplified Keplerian)
    # v_rad = K * (cos(nu + omega) + e * cos(omega))
    # where K is the semi-amplitude
    # K = (2*pi*a*sin(i)) / (T * sqrt(1-e^2)) * (R0 / 1pc conversion)
    # a in arcsec -> a in AU: a_AU = a_arcsec * R0_pc
    # Orbital velocity semi-amplitude in km/s
    # v_K = sqrt(GM/a) * sin(i) / sqrt(1-e^2)
    # G in AU^3/(M_sun * yr^2): G = 4*pi^2
    # v in AU/yr -> km/s: 1 AU/yr = 4.74 km/s
    v_K = np.sqrt(4 * np.pi**2 * M / a_AU) * np.sin(i) / np.sqrt(1 - e**2)
    v_K_kms = v_K * 4.74  # convert AU/yr to km/s

    v_rad = v_K_kms * (np.cos(nu + omega_precessed) + e * np.cos(omega_precessed))

    # Add Schwarzschild gravitational redshift at pericentre (detected by GRAVITY)
    # Exact redshift: 1 + z = 1 / sqrt(1 - R_s/r), with R_s = 2GM/c^2.
    # In velocity units: v_grav = (c / 1000) * z  [km/s].
    R_s_m = 2 * G_SI * M * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m
    r = a * (1 - e**2) / (1 + e * np.cos(nu)) * R0  # in AU
    z_grav = 1.0 / np.sqrt(1.0 - R_s_AU / r) - 1.0
    v_grav = (c_si / 1000.0) * z_grav  # km/s

    v_total = v_rad + v_grav

    # Apply TEP conformal factor if eta_inject != 0
    # The TEP scaling: v_obs = v_local / A(r)
    if eta_inject != 0.0:
        r_arcsec = a * (1 - e**2) / (1 + e * np.cos(nu))
        A = _tep_conformal_factor(r_arcsec, M, R0, eta_inject)
        v_total = v_total / A

    # Add noise
    rng = np.random.default_rng(123)
    v_obs = v_total + rng.normal(0, noise_kms, n_points)

    # Uncertainties (SINFONI: ~10-30 km/s, GRAVITY: ~7-10 km/s)
    sigma_v = np.where(t_arr < 2017.0, 30.0, noise_kms)

    spectroscopy = []
    for k in range(n_points):
        spectroscopy.append({
            "t_yr": float(t_arr[k]),
            "v_rad_kms": float(v_obs[k]),
            "sigma_v_kms": float(sigma_v[k]),
            "instrument": "SINFONI" if t_arr[k] < 2017.0 else "GRAVITY",
        })

    return spectroscopy


# ============================================================
# Data download attempts
# ============================================================
def try_download_vizier():
    """Attempt to download S-star data from VizieR.

    VizieR is the standard repository for published astronomical data tables.
    The GRAVITY Collaboration data may be available as VizieR catalogues.
    """
    import ssl

    sstar_dir = RAW_DIR / "sstar"
    sstar_dir.mkdir(parents=True, exist_ok=True)

    # VizieR catalogue URLs for GRAVITY Collaboration papers
    # These are the standard VizieR accession URLs
    vizier_urls = [
        {
            "name": "GRAVITY 2020 S2 astrometry",
            "url": "https://cdsarc.u-strasbg.fr/viz-bin/cat/J/A+A/636/L5",
            "target": "sstar/vizier_gravity_2020.txt",
        },
    ]

    # Create an SSL context that doesn't verify certificates
    # (VizieR uses a self-signed certificate)
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    downloaded = []
    for item in vizier_urls:
        target = RAW_DIR / item["target"]
        if target.exists():
            downloaded.append({"name": item["name"], "path": str(target), "status": "already_exists"})
            continue
        try:
            req = urllib.request.Request(item["url"], headers={"User-Agent": "TEP-BH-Pipeline/0.1"})
            with urllib.request.urlopen(req, timeout=TIMEOUT_SECONDS, context=ctx) as resp:
                data = resp.read()
                if len(data) > 100:
                    target.write_bytes(data)
                    downloaded.append({
                        "name": item["name"],
                        "path": str(target),
                        "status": "downloaded",
                        "bytes": len(data),
                    })
        except (urllib.error.URLError, urllib.error.HTTPError, Exception) as e:
            downloaded.append({"name": item["name"], "status": "failed", "error": str(e)})

    return downloaded


# ============================================================
# Gillessen 2017 real data download and parse
# ============================================================
def download_gillessen_2017():
    """Download and parse the Gillessen et al. 2017 S2 astrometry and RV data.

    VizieR catalogue J/ApJ/837/30 contains table5.dat with 216 epochs of
    S-star astrometry and radial velocities relative to SgrA*. This function
    downloads the raw CDS files (table5.dat, table3.dat, and ReadMe), parses
    them with the CDS format, and extracts S2 astrometry, spectroscopy, and
    the orbital parameters from table3.

    Coordinate convention:
      The pipeline orbit model projects the orbital plane onto the sky with
      x_sky along the line of nodes and y_sky perpendicular. The VizieR table5
      columns oRA-S2 and oDE-S2 are right-ascension and declination offsets.
      Empirically, the model x_sky matches oDE-S2 and the model y_sky matches
      oRA-S2, so the columns are mapped accordingly.

    Returns:
        dict with keys 'status', 'astrometry', 'spectroscopy', 'orbital_params',
        'n_astrometry', 'n_spectroscopy', 'source', or {'status': 'failed',
        'error': ...}.
    """
    from astropy.io import ascii

    sstar_dir = RAW_DIR / "sstar"
    sstar_dir.mkdir(parents=True, exist_ok=True)

    table3_path = sstar_dir / "table3.dat"
    table5_path = sstar_dir / "table5.dat"
    readme_path = sstar_dir / "ReadMe"

    urls = {
        table3_path: "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/837/30/table3.dat",
        table5_path: "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/837/30/table5.dat",
        readme_path: "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/837/30/ReadMe",
    }

    # Download if not present
    for target, url in urls.items():
        if not target.exists():
            try:
                with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as resp:
                    data = resp.read()
                    target.write_bytes(data)
            except Exception as e:
                return {"status": "failed", "error": f"Failed to download {url}: {e}"}

    # Parse table3 for orbital parameters
    try:
        tbl3 = ascii.read(table3_path, format="cds", readme=readme_path)
        tbl3_df = tbl3.to_pandas()
        s2_row = tbl3_df[tbl3_df["Star"] == "S2"].iloc[0]
        orbital_params = {
            "star": "S2",
            "M_BH_Msun": 4.28e6,
            "M_BH_err": 0.05e6,
            "R0_pc": 8320.0,
            "R0_err": 60.0,
            "T_orb_yr": float(s2_row["Per"]),
            "T_orb_err": float(s2_row["e_Per"]),
            "a_arcsec": float(s2_row["a"]),
            "a_err": float(s2_row["e_a"]),
            "eccentricity": float(s2_row["e"]),
            "e_err": float(s2_row["e_e"]),
            "inclination_deg": float(s2_row["i"]),
            "i_err": float(s2_row["e_i"]),
            "omega_deg": float(s2_row["w"]),
            "omega_err": float(s2_row["e_w"]),
            "Omega_deg": float(s2_row["Omega"]),
            "Omega_err": float(s2_row["e_Omega"]),
            "t_peri_yr": float(s2_row["Tp"]),
            "t_peri_err": float(s2_row["e_Tp"]),
            "references": [
                "Gillessen et al. 2017, ApJ 837, 30 (VizieR J/ApJ/837/30, table3.dat)",
                "GRAVITY Collaboration 2018a, A&A 615, L15 (DOI: 10.1051/0004-6361/201833718)",
                "GRAVITY Collaboration 2020, A&A 636, L5 (DOI: 10.1051/0004-6361/202037813)",
                "GRAVITY Collaboration 2022, A&A 657, L12 (DOI: 10.1051/0004-6361/202142465)",
            ],
        }
    except Exception as e:
        return {"status": "failed", "error": f"Failed to parse table3.dat: {e}"}

    # Parse table5 for astrometry and spectroscopy
    try:
        tbl5 = ascii.read(table5_path, format="cds", readme=readme_path)
        df = tbl5.to_pandas()
    except Exception as e:
        return {"status": "failed", "error": f"Failed to parse table5.dat: {e}"}

    # Extract S2 astrometry (Flag == 'a')
    # Map model x_sky -> oDE-S2, model y_sky -> oRA-S2 (see docstring)
    astrom_df = df[df["Flag"] == "a"][["Date", "oRA-S2", "e_oRA-S2", "oDE-S2", "e_oDE-S2"]].dropna()
    astrometry = [
        {
            "t_yr": float(row["Date"]),
            "x_arcsec": float(row["oDE-S2"]) / 1000.0,
            "y_arcsec": float(row["oRA-S2"]) / 1000.0,
            "sigma_x_arcsec": float(row["e_oDE-S2"]) / 1000.0,
            "sigma_y_arcsec": float(row["e_oRA-S2"]) / 1000.0,
        }
        for _, row in astrom_df.iterrows()
    ]

    # Extract S2 spectroscopy (Flag == 'rv')
    spec_df = df[df["Flag"] == "rv"][["Date", "RV-S2", "e_RV-S2"]].dropna()
    spectroscopy = [
        {
            "t_yr": float(row["Date"]),
            "v_rad_kms": float(row["RV-S2"]),
            "sigma_v_kms": float(row["e_RV-S2"]),
        }
        for _, row in spec_df.iterrows()
    ]

    return {
        "status": "success",
        "source": "Gillessen et al. 2017, ApJ 837, 30 (VizieR J/ApJ/837/30, table5.dat/table3.dat)",
        "bibcode": "2017ApJ...837...30G",
        "n_astrometry": len(astrometry),
        "n_spectroscopy": len(spectroscopy),
        "astrometry": astrometry,
        "spectroscopy": spectroscopy,
        "orbital_params": orbital_params,
    }


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 00 (INFERENCE): S-STAR DATA ACQUISITION", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("This is the PRIMARY FALSIFIABLE TEST of TEP at black-hole scales.", "INFO")
    print_status("Data source: Gillessen et al. 2017 (VizieR J/ApJ/837/30) S2 astrometry and spectroscopy", "INFO")
    print_status("")

    # --- Attempt to get real S2 data from CDS ---
    print_status("--- Downloading real S2 data from CDS (Gillessen et al. 2017) ---", "INFO")
    real_data = download_gillessen_2017()

    # --- Fallback VizieR attempt ---
    print_status("--- Attempting GRAVITY VizieR download ---", "INFO")
    downloads = try_download_vizier()
    for d in downloads:
        print_status(f"  {d['name']}: {d['status']}", "INFO")
    print_status("")

    # --- Save orbital parameters ---
    print_status("--- Published orbital parameters ---", "INFO")
    orbital_params = {
        "S2": S2_ORBITAL_PARAMS,
        **{k: {**v, "star": k} for k, v in OTHER_STARS.items()},
    }

    # If real data downloaded successfully, use its orbital parameters
    if real_data["status"] == "success":
        orbital_params["S2"] = real_data["orbital_params"]

    for star, params in orbital_params.items():
        if star == "S2":
            print_status(f"  {star}: M_BH = {params['M_BH_Msun']:.3e} M_sun", "INFO")
            print_status(f"    R0 = {params['R0_pc']} pc", "INFO")
            print_status(f"    T = {params['T_orb_yr']} yr, a = {params['a_arcsec']} arcsec", "INFO")
            print_status(f"    e = {params['eccentricity']}, i = {params['inclination_deg']}°", "INFO")
            print_status(f"    omega = {params['omega_deg']}°, Omega = {params['Omega_deg']}°, t_peri = {params['t_peri_yr']}", "INFO")
        else:
            print_status(f"  {star}: T = {params['T_orb_yr']} yr, a = {params['a_arcsec']} arcsec", "INFO")

    params_path = PROCESSED_DIR / "sstar_orbital_params.json"
    write_json(params_path, orbital_params)
    print_status(f"  Saved: {params_path}", "INFO")
    print_status("")

    # Check for TEP injection (only for synthetic testing)
    import os
    eta_inject = float(os.environ.get("TEP_ETA_INJECT", "0.0"))

    if real_data["status"] == "success" and eta_inject == 0.0:
        # Use real data
        print_status("--- Using REAL S2 data from Gillessen et al. 2017 ---", "INFO")
        print_status(f"  Source: {real_data['source']}", "INFO")
        print_status(f"  Bibcode: {real_data['bibcode']}", "INFO")

        astrometry = real_data["astrometry"]
        spectroscopy = real_data["spectroscopy"]
        is_synthetic = False
        synthetic_note = "Real astrometric and spectroscopic measurements of S2 relative to SgrA* from Gillessen et al. 2017 (J/ApJ/837/30)."
    else:
        # Use synthetic data (either real download failed, or TEP injection requested)
        if eta_inject > 0:
            print_status(f"--- Generating synthetic S2 data with TEP injection (eta={eta_inject}) ---", "INFO")
            print_status(f"  The conformal factor A(r) = exp(+eta/(3*r_rs)) is applied", "INFO")
            print_status(f"  to simulate what would be observed if TEP is true.", "INFO")
        else:
            print_status("--- Real data download failed; generating synthetic S2 data from published parameters ---", "WARN")
            print_status(f"  Reason: {real_data.get('error', 'unknown')}", "WARN")

        astrometry = generate_s2_synthetic_astrometry(S2_ORBITAL_PARAMS, n_points=150, noise_mas=0.1, eta_inject=eta_inject)
        spectroscopy = generate_s2_synthetic_spectroscopy(S2_ORBITAL_PARAMS, n_points=100, noise_kms=10, eta_inject=eta_inject)
        is_synthetic = True
        synthetic_note = "Generated from published GRAVITY Collaboration orbital parameters for pipeline testing."

    print_status(f"  Astrometry: {len(astrometry)} points", "INFO")
    print_status(f"  Spectroscopy: {len(spectroscopy)} points", "INFO")
    if not is_synthetic:
        print_status(f"  Date range: {min(d['t_yr'] for d in astrometry):.2f} - {max(d['t_yr'] for d in astrometry):.2f} yr", "INFO")

    # Save
    astrometry_path = PROCESSED_DIR / "sstar_astrometry.json"
    spectroscopy_path = PROCESSED_DIR / "sstar_spectroscopy.json"
    write_json(astrometry_path, {
        "star": "S2",
        "data_type": "astrometry",
        "is_synthetic": is_synthetic,
        "eta_inject": eta_inject,
        "source": real_data.get("source") if not is_synthetic else None,
        "bibcode": real_data.get("bibcode") if not is_synthetic else None,
        "synthetic_note": synthetic_note,
        "n_points": len(astrometry),
        "units": {"t": "yr", "x": "arcsec", "y": "arcsec", "sigma": "arcsec"},
        "data": astrometry,
    })

    write_json(spectroscopy_path, {
        "star": "S2",
        "data_type": "spectroscopy",
        "is_synthetic": is_synthetic,
        "eta_inject": eta_inject,
        "source": real_data.get("source") if not is_synthetic else None,
        "bibcode": real_data.get("bibcode") if not is_synthetic else None,
        "synthetic_note": synthetic_note,
        "n_points": len(spectroscopy),
        "units": {"t": "yr", "v_rad": "km/s", "sigma_v": "km/s"},
        "data": spectroscopy,
    })

    print_status(f"  Saved: {astrometry_path}", "INFO")
    print_status(f"  Saved: {spectroscopy_path}", "INFO")
    print_status("")

    # --- Summary ---
    print_status("=" * 70, "TITLE")
    print_status("SUMMARY", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(f"  S2 orbital parameters: saved ({len(orbital_params)} stars)", "INFO")
    print_status(f"  VizieR downloads: {sum(1 for d in downloads if d.get('status') == 'downloaded')} succeeded", "INFO")
    if is_synthetic:
        print_status(f"  Synthetic astrometry: {len(astrometry)} points", "INFO")
        print_status(f"  Synthetic spectroscopy: {len(spectroscopy)} points", "INFO")
    else:
        print_status(f"  REAL astrometry: {len(astrometry)} points", "INFO")
        print_status(f"  REAL spectroscopy: {len(spectroscopy)} points", "INFO")
    print_status("")
    print_status("  NEXT: step_44_gr_fit.py — reproduce the conventional GR fit", "INFO")
    print_status("")

    # --- Save results ---
    result = {
        "step": STEP_ID,
        "description": "S-star data acquisition for TEP inference pipeline",
        "timestamp": datetime.now().isoformat(),
        "orbital_params": {k: {kk: vv for kk, vv in v.items() if kk != "references"}
                          for k, v in orbital_params.items()},
        "downloads": downloads,
        "real_data": real_data if real_data["status"] == "success" else {"status": "failed", "error": real_data.get("error")},
        "synthetic_data": {
            "astrometry": {"n_points": len(astrometry), "is_synthetic": is_synthetic},
            "spectroscopy": {"n_points": len(spectroscopy), "is_synthetic": is_synthetic},
        },
        "data_sources": [S2_ORBITAL_PARAMS["references"][0], real_data.get("source", "")] if not is_synthetic else S2_ORBITAL_PARAMS["references"],
        "status": "success",
    }

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, result)
    print_status(f"Results saved to {result_path}", "INFO")
    key_result = (f"Real S2 data downloaded from Gillessen et al. 2017: "
                  f"{len(astrometry)} astrometry + {len(spectroscopy)} spectroscopy points") if not is_synthetic else "Synthetic S2 data generated from published parameters"
    finalize_result(STEP_ID, result, "S-star data acquisition", key_result)


if __name__ == "__main__":
    main()
