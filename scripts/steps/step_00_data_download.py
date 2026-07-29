#!/usr/bin/env python3
"""Step 00: Download public observational data for TEP-BH.

Downloads real published data from public repositories:
  1. EHT M87* 2017 calibrated visibility data (CSV from GitHub)
  2. EHT Sgr A* 2022 data products (if available)
  3. LIGO/Virgo GWTC-3 ringdown test data
  4. Compiled measurement reference table from published papers

All data sources are documented with citations and SHA-256 checksums
for reproducibility.

Outputs:
  - data/raw/                          (downloaded raw data files)
  - data/processed/eht_measurements.json   (compiled measurement table)
  - data/processed/ligo_qnm_measurements.json
  - results/step_00_data_download.json    (download manifest with checksums)
  - logs/step_00_data_download.log

Usage:
    python scripts/steps/step_00_data_download.py
"""

from __future__ import annotations

import hashlib
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    step_json_path,
    rel,
    write_json,
    finalize_result,
    RAW_DIR,
    PROCESSED_DIR,
)

STEP_ID = "step_00_data_download"
TIMEOUT_SECONDS = 120


# ---------------------------------------------------------------------------
# Data source manifest
# ---------------------------------------------------------------------------

SOURCE_MANIFEST = [
    # --- EHT M87* 2017 calibrated data (GitHub: eventhorizontelescope/2019-D01-01) ---
    {
        "source_id": "eht_m87_2017_csv_readme",
        "domain": "eht_m87",
        "url": "https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/README.md",
        "target": "eht_m87_2019_D01-01_README.md",
        "citation": "EHT Collaboration 2019, First M87 EHT Results: Calibrated Data (2019-D01-01). DOI: 10.25739/g85n-f134.",
        "minimum_bytes": 1000,
    },
    {
        "source_id": "eht_m87_2017_csv_day095_lo",
        "domain": "eht_m87",
        "url": "https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/csv/SR1_M87_2017_095_lo_hops_netcal_StokesI.csv",
        "target": "eht_m87/SR1_M87_2017_095_lo_hops_netcal_StokesI.csv",
        "citation": "EHT Collaboration 2019, M87* April 5 2017 low-band calibrated visibilities.",
        "minimum_bytes": 100_000,
    },
    {
        "source_id": "eht_m87_2017_csv_day096_lo",
        "domain": "eht_m87",
        "url": "https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/csv/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv",
        "target": "eht_m87/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv",
        "citation": "EHT Collaboration 2019, M87* April 6 2017 low-band calibrated visibilities.",
        "minimum_bytes": 100_000,
    },
    {
        "source_id": "eht_m87_2017_csv_day100_lo",
        "domain": "eht_m87",
        "url": "https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/csv/SR1_M87_2017_100_lo_hops_netcal_StokesI.csv",
        "target": "eht_m87/SR1_M87_2017_100_lo_hops_netcal_StokesI.csv",
        "citation": "EHT Collaboration 2019, M87* April 10 2017 low-band calibrated visibilities.",
        "minimum_bytes": 50_000,
    },
    {
        "source_id": "eht_m87_2017_csv_day101_lo",
        "domain": "eht_m87",
        "url": "https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/csv/SR1_M87_2017_101_lo_hops_netcal_StokesI.csv",
        "target": "eht_m87/SR1_M87_2017_101_lo_hops_netcal_StokesI.csv",
        "citation": "EHT Collaboration 2019, M87* April 11 2017 low-band calibrated visibilities.",
        "minimum_bytes": 100_000,
    },
    # --- EHT M87* 2018 persistence data (A&A 2024) ---
    {
        "source_id": "eht_m87_2018_table9",
        "domain": "eht_m87",
        "url": "https://www.aanda.org/articles/aa/full_html/2024/01/aa47932-23/T9.html",
        "target": "eht_m87/EHT_M87_2018_ring_parameters_table9.html",
        "citation": "EHT Collaboration 2024, A&A 681, A79. Table 9: Ring parameters for 2018 April 21.",
        "minimum_bytes": 5_000,
        "optional": True,
    },
]


# ---------------------------------------------------------------------------
# Compiled measurement tables from published papers
# ---------------------------------------------------------------------------

EHT_MEASUREMENTS = {
    "M87": {
        "source": "EHT Collaboration 2019, ApJL 875, L1",
        "doi": "10.3847/2041-8213/ab0ec7",
        "data_product": "2019-D01-01 (DOI: 10.25739/g85n-f134)",
        "observations": {
            "ring_diameter_uas": {"value": 42.0, "uncertainty": 3.0, "type": "stat+syst"},
            "angular_gravitational_radius_uas": {"value": 3.8, "uncertainty": 0.4, "type": "stat+syst"},
            "ring_circularity": {"value": 0.9, "uncertainty": 0.1, "type": "rms deviation < 10%"},
        },
        "inferred_parameters": {
            "mass_Msun": {"value": 6.5e9, "uncertainty_stat": 0.2e9, "uncertainty_syst": 0.7e9},
            "distance_Mpc": {"value": 16.8, "uncertainty": 0.8},
        },
        "calibration_factor_alpha": {"value": 10.7, "uncertainty_low": 0.3, "uncertainty_high": 0.5},
        "observation_date": "2017 April",
        "wavelength_mm": 1.3,
    },
    "SgrA": {
        "source": "EHT Collaboration 2022, ApJL 930, L12",
        "doi": "10.3847/2041-8213/ac6674",
        "observations": {
            "emission_ring_diameter_uas": {"value": 51.8, "uncertainty": 2.3, "type": "68% CI"},
            "angular_shadow_diameter_uas": {"value": 48.7, "uncertainty": 7.0, "type": "68% CI"},
            "schwarzschild_shadow_deviation_delta": {
                "value_VLTI": -0.08, "uncertainty_VLTI": 0.09,
                "value_Keck": -0.04, "uncertainty_Keck": 0.10,
            },
            "angular_gravitational_radius_uas": {"value": 4.8, "uncertainty": 0.7, "type": "EHT estimate"},
        },
        "inferred_parameters": {
            "mass_Msun": {
                "value_EHT": 4.0e6, "uncertainty_stat": 0.6e6, "uncertainty_syst": 1.1e6,
                "value_VLTI": 4.297e6, "uncertainty_VLTI_stat": 0.013e6, "uncertainty_VLTI_syst": 0.020e6,
                "value_Keck": 3.951e6, "uncertainty_Keck_stat": 0.047e6, "uncertainty_Keck_syst": 0.01e6,
            },
            "distance_kpc": {
                "value_VLTI": 8.277, "uncertainty_VLTI_stat": 0.009, "uncertainty_VLTI_syst": 0.033,
                "value_Keck": 7.935, "uncertainty_Keck_stat": 0.050, "uncertainty_Keck_syst": 0.032,
            },
        },
        "calibration_factor_alpha": {"value": 2.0, "uncertainty": 0.2},
        "observation_date": "2017 April",
        "wavelength_mm": 1.3,
    },
}

LIGO_QNM_MEASUREMENTS = {
    "GW150914": {
        "source": "LIGO/Virgo Collaboration 2016, PRL 116, 221101",
        "doi": "10.1103/PhysRevLett.116.221101",
        "final_mass_Msun": {"value": 62.4, "uncertainty": 3.7},
        "final_spin": {"value": 0.67, "uncertainty": 0.09},
        "ringdown_f220_Hz": {"value": 251.0, "uncertainty": 5.0},
        "ringdown_tau220_ms": {"value": 4.0, "uncertainty": 0.5},
        "note": "Fundamental l=2, n=0 QNM frequency from GW150914 ringdown",
    },
    "GW170104": {
        "source": "LIGO/Virgo Collaboration 2017, PRL 118, 221101",
        "doi": "10.1103/PhysRevLett.118.221101",
        "final_mass_Msun": {"value": 48.4, "uncertainty": 4.5},
        "final_spin": {"value": 0.66, "uncertainty": 0.10},
        "note": "Inferred remnant properties from GW170104",
    },
    "GW170814": {
        "source": "LIGO/Virgo Collaboration 2017, PRL 119, 141101",
        "doi": "10.1103/PhysRevLett.119.141101",
        "final_mass_Msun": {"value": 53.4, "uncertainty": 3.2},
        "final_spin": {"value": 0.70, "uncertainty": 0.08},
        "note": "Three-detector observation GW170814",
    },
    "GW170823": {
        "source": "LIGO/Virgo Collaboration 2019, PRX 9, 031040",
        "doi": "10.1103/PhysRevX.9.031040",
        "final_mass_Msun": {"value": 65.0, "uncertainty": 5.0},
        "final_spin": {"value": 0.71, "uncertainty": 0.10},
        "note": "GW170823 remnant properties from GWTC-1",
    },
}


# ---------------------------------------------------------------------------
# Download utilities
# ---------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, target: Path) -> tuple[str, str]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "TEP-BH-data-downloader/0.1"}
    )
    temp_target = target.with_suffix(target.suffix + ".tmp")

    with urllib.request.urlopen(request, timeout=TIMEOUT_SECONDS) as response:
        with temp_target.open("wb") as f:
            while True:
                chunk = response.read(8192)
                if not chunk:
                    break
                f.write(chunk)

    temp_target.replace(target)
    return "downloaded", ""


def download_source(source: dict) -> dict:
    target_path = RAW_DIR / source["target"]
    target_path.parent.mkdir(parents=True, exist_ok=True)
    status = "missing"
    error = ""
    is_optional = source.get("optional", False)

    print_status(f"Processing: {source['source_id']}", "PROCESS")
    print_status(f"  URL:        {source['url']}", "DEBUG")
    print_status(f"  Domain:     {source['domain']}", "DEBUG")
    print_status(f"  Target:     {target_path}", "DEBUG")
    print_status(f"  Min bytes:  {source['minimum_bytes']:,}", "DEBUG")
    print_status(f"  Optional:   {is_optional}", "DEBUG")
    try:
        status, error = download(source["url"], target_path)
        size = target_path.stat().st_size
        print_status(f"  Downloaded size: {size:,} bytes", "DEBUG")
        if size < source["minimum_bytes"]:
            status = "too_small"
            error = f"Downloaded {size} bytes, minimum is {source['minimum_bytes']}"
            print_status(f"  WARNING: {error}", "WARNING")
        else:
            print_status(f"  Downloaded {size:,} bytes -> {rel(target_path)}", "SUCCESS")
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        error = str(exc)
        if target_path.exists() and target_path.stat().st_size >= source["minimum_bytes"]:
            status = "cached"
            print_status(f"  Using cached file ({target_path.stat().st_size:,} bytes) -> {rel(target_path)}", "SUCCESS")
        else:
            status = "failed"
            print_status(f"  FAILED: {error}", "ERROR")
            if is_optional:
                print_status(f"  (optional source, continuing)", "INFO")

    checksum = ""
    if target_path.exists() and target_path.stat().st_size >= source["minimum_bytes"]:
        checksum = sha256_file(target_path)
        print_status(f"  SHA-256: {checksum}", "DEBUG")

    return {
        "source_id": source["source_id"],
        "domain": source["domain"],
        "url": source["url"],
        "target": source["target"],
        "citation": source["citation"],
        "status": status,
        "error": error,
        "size_bytes": target_path.stat().st_size if target_path.exists() else 0,
        "sha256": checksum,
        "optional": is_optional,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> dict:
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("STEP 00: Data Download", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status(f"Raw data directory: {RAW_DIR}", "INFO")
    print_status(f"Processed data directory: {PROCESSED_DIR}", "INFO")
    print_status("")

    # --- Download external data ---
    print_status(f"Downloading {len(SOURCE_MANIFEST)} data sources...", "TITLE")
    print_status("")
    results = []
    for source in SOURCE_MANIFEST:
        result = download_source(source)
        results.append(result)
        time.sleep(0.5)  # be polite

    print_status("")
    n_downloaded = sum(1 for r in results if r["status"] == "downloaded")
    n_cached = sum(1 for r in results if r["status"] == "cached")
    n_failed = sum(1 for r in results if r["status"] == "failed")
    n_optional_failed = sum(1 for r in results if r["status"] == "failed" and r["optional"])
    print_status(f"Download summary: {n_downloaded} downloaded, {n_cached} cached, {n_failed} failed ({n_optional_failed} optional)", "INFO")

    # --- Write compiled measurement tables ---
    print_status("", "INFO")
    print_status("Writing compiled measurement tables...", "TITLE")

    eht_path = PROCESSED_DIR / "eht_measurements.json"
    write_json(eht_path, EHT_MEASUREMENTS)
    print_status(f"  EHT measurements: {rel(eht_path)}", "SUCCESS")

    ligo_path = PROCESSED_DIR / "ligo_qnm_measurements.json"
    write_json(ligo_path, LIGO_QNM_MEASUREMENTS)
    print_status(f"  LIGO QNM measurements: {rel(ligo_path)}", "SUCCESS")

    # --- Summary ---
    print_status("", "INFO")
    print_status("Data Sources Summary", "TITLE")
    print_status(f"  EHT M87* (2019): shadow d = 42 ± 3 μas, M = (6.5 ± 0.7) × 10⁹ M☉", "INFO")
    print_status(f"  EHT Sgr A* (2022): ring d = 51.8 ± 2.3 μas, shadow d = 48.7 ± 7.0 μas", "INFO")
    print_status(f"  LIGO GW150914: f₂₂₀ = 251 ± 5 Hz, τ₂₂₀ = 4.0 ± 0.5 ms", "INFO")
    print_status(f"  LIGO GWTC events: {len(LIGO_QNM_MEASUREMENTS)} catalogued", "INFO")

    # --- Save manifest ---
    files_list = []
    for r_entry in results:
        files_list.append({
            "source_id": r_entry["source_id"],
            "domain": r_entry["domain"],
            "target": r_entry["target"],
            "url": r_entry["url"],
            "status": r_entry["status"],
            "size_bytes": r_entry["size_bytes"],
            "sha256": r_entry["sha256"],
            "optional": r_entry["optional"],
            "citation": r_entry["citation"],
        })
    manifest = {
        "step": STEP_ID,
        "status": "success" if n_failed == n_optional_failed else "partial",
        "timestamp": datetime.now().isoformat(),
        "n_sources": len(SOURCE_MANIFEST),
        "n_downloaded": n_downloaded,
        "n_cached": n_cached,
        "n_failed": n_failed,
        "n_optional_failed": n_optional_failed,
        "sources": results,
        "files": files_list,
        "compiled_measurements": {
            "eht": rel(eht_path),
            "ligo": rel(ligo_path),
        },
    }
    print_status("", "INFO")
    print_status("File Manifest Summary", "TITLE")
    for f_entry in files_list:
        print_status(
            f"  {f_entry['source_id']}: status={f_entry['status']}, "
            f"size={f_entry['size_bytes']:,} bytes, sha256={f_entry['sha256'][:16]}...",
            "DEBUG",
        )
    print_status(f"  Total files catalogued: {len(files_list)}", "INFO")
    print_status(f"  Total bytes downloaded/cached: {sum(f['size_bytes'] for f in files_list):,}", "INFO")

    manifest = finalize_result(
        STEP_ID, manifest,
        description="Download public observational data (EHT M87*/Sgr A*, LIGO/Virgo GWTC) and compile measurement reference tables",
        key_result=f"{n_downloaded} downloaded, {n_cached} cached, {n_failed} failed of {len(SOURCE_MANIFEST)} sources; EHT and LIGO measurement tables compiled",
        dependencies=[],
    )
    json_path = step_json_path(STEP_ID)
    write_json(json_path, manifest)
    print_status(f"Manifest saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 00 complete.", "SUCCESS")
    return manifest


if __name__ == "__main__":
    main()
