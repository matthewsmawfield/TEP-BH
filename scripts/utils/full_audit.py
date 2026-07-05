#!/usr/bin/env python3
"""TEP-TH Full Audit: Integrity, Statistics, Manuscript Consistency."""

from __future__ import annotations

import json
import pathlib
import sys

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
MANUSCRIPT = PROJECT_ROOT / "manuscripts" / "27-TEP-TH-v0.2-Thika.md"

AUDIT_ISSUES: list[str] = []
AUDIT_WARNINGS: list[str] = []


def audit_result_files() -> None:
    print("=" * 70)
    print("1. RESULT FILE INTEGRITY")
    print("=" * 70)
    step_names = [
        "step_00_temporal_horizon_mapping",
        "step_01_matter_frame_curvature",
        "step_02_geodesic_completeness",
        "step_03_effective_stress_energy",
        "step_04_full_bbn_abundances",
        "step_05_recombination_visibility",
        "step_06_cmb_blackbody_origin",
        "step_07_entropy_arrow",
        "step_08_primordial_perturbation_boundary",
        "step_09b_native_tensor_integration",
        "step_10_cmb_lss_class",
    ]
    all_present = True
    for step in step_names:
        jpath = RESULTS_DIR / f"{step}.json"
        cpath = RESULTS_DIR / f"{step}.csv"
        if not jpath.exists():
            print(f"  [FAIL] {step}: JSON missing")
            AUDIT_ISSUES.append(f"{step}: missing JSON result")
            all_present = False
            continue
        try:
            data = json.loads(jpath.read_text())
        except json.JSONDecodeError as exc:
            print(f"  [FAIL] {step}: invalid JSON ({exc})")
            AUDIT_ISSUES.append(f"{step}: invalid JSON")
            all_present = False
            continue
        txt = json.dumps(data).lower()
        if "synthetic_fallback" in txt:
            print(f"  [FAIL] {step}: synthetic_fallback flag found")
            AUDIT_ISSUES.append(f"{step}: synthetic_fallback")
            all_present = False
        elif "fake_data" in txt:
            print(f"  [FAIL] {step}: fake_data flag found")
            AUDIT_ISSUES.append(f"{step}: fake_data")
            all_present = False
        else:
            print(f"  [PASS] {step}: valid JSON, no synthetic flags")
    if all_present:
        print(f"\n  All {len(step_names)} result files present and clean.")


def audit_codebase() -> None:
    print("\n" + "=" * 70)
    print("2. CODEBASE SCAN")
    print("=" * 70)
    issues = []
    for py_file in SCRIPTS_DIR.rglob("*.py"):
        if py_file.name in ("full_audit.py", "pipeline_audit.py", "check2_audit.py"):
            continue
        content = py_file.read_text()
        lower = content.lower()
        if "synthetic_fallback" in lower:
            issues.append(f"{py_file.relative_to(PROJECT_ROOT)}: synthetic_fallback")
        if "fake_data" in lower:
            issues.append(f"{py_file.relative_to(PROJECT_ROOT)}: fake_data")
        # Check for unguarded placeholder strings (not in comments)
        for i, line in enumerate(content.split("\n")):
            l = line.lower()
            if (
                "placeholder" in l
                and "notimplemented" not in l
                and "#" not in l
                and '"""' not in l
                and "'" not in l
            ):
                # This is a crude heuristic; skip obvious docstrings
                stripped = line.strip()
                if not stripped.startswith("#") and not stripped.startswith('"""'):
                    issues.append(f"{py_file.relative_to(PROJECT_ROOT)}:{i+1}: placeholder")

    if issues:
        for issue in issues[:10]:
            print(f"  [WARN] {issue}")
        if len(issues) > 10:
            print(f"  ... and {len(issues)-10} more")
        AUDIT_WARNINGS.extend(issues)
    else:
        print("  [PASS] No synthetic/placeholder flags in production code.")

    # Check for unused stub files that claim to be placeholders
    stub_files = [
        SCRIPTS_DIR / "utils" / "bbn_cross_validation.py",
        SCRIPTS_DIR / "utils" / "cmb_class_resolver.py",
    ]
    for stub in stub_files:
        if stub.exists():
            content = stub.read_text().lower()
            if "placeholder" in content or "stub" in content:
                print(f"  [WARN] {stub.relative_to(PROJECT_ROOT)}: marked as stub/placeholder")
                AUDIT_WARNINGS.append(f"{stub.name}: stub/placeholder module exists")


def audit_physical_numbers() -> None:
    print("\n" + "=" * 70)
    print("4. KEY PHYSICAL NUMBERS")
    print("=" * 70)

    # BBN
    bbn = json.loads((RESULTS_DIR / "step_04_full_bbn_abundances.json").read_text())
    tep = bbn["TEP_abundances"]
    print(f"  TEP Y_p:      {tep['Y_p']:.4f}   (Planck ~0.247)")
    print(f"  TEP D/H:      {tep['D_H']:.2e}   (Planck ~2.6e-5)")
    print(f"  TEP He3/H:    {tep['He3_H']:.2e}")
    print(f"  TEP Li7/H:    {tep['Li7_H']:.2e}")
    print(f"  TEP N_eff:    {tep['N_eff']:.3f}")

    # Recombination (top-level keys for screened TEP results)
    rec = json.loads((RESULTS_DIR / "step_05_recombination_visibility.json").read_text())
    z_star = rec.get("recombination_epoch", {}).get("z_star_TEP")
    r_s = rec.get("sound_horizon", {}).get("r_s_TEP")
    if z_star is not None:
        print(f"  TEP z_*:      {z_star:.1f}    (Planck 1089.92)")
    if r_s is not None:
        print(f"  TEP r_s:      {r_s:.1f} Mpc  (Planck ~147)")

    # CMB/LSS
    cmb = json.loads((RESULTS_DIR / "step_10_cmb_lss_class.json").read_text())
    lcdm = cmb["lcdm_class_output"]
    print(f"  CLASS z_*:    {lcdm['z_star']:.1f}")
    print(f"  CLASS r_s:    {lcdm['r_s_Mpc']:.1f} Mpc")
    print(f"  CLASS sigma8: {lcdm['sigma8']:.4f}")
    print(f"  CMB score:    {cmb['cmb_scorecard']['cmb_consistent']}")
    print(f"  LSS score:    {cmb['lss_scorecard']['lss_consistent']}")


def audit_manuscript() -> None:
    print("\n" + "=" * 70)
    print("5. MANUSCRIPT CONSISTENCY")
    print("=" * 70)
    if not MANUSCRIPT.exists():
        print(f"  [FAIL] Manuscript not found: {MANUSCRIPT}")
        AUDIT_ISSUES.append("Manuscript missing")
        return

    text = MANUSCRIPT.read_text().lower()

    # Pipeline step count
    if "nine-step" in text and "eleven-step" not in text:
        print(f"  [WARN] Manuscript says 'nine-step' pipeline; codebase has 11 steps (00-10)")
        AUDIT_WARNINGS.append("Manuscript step count outdated (says nine, should be eleven)")
    else:
        print(f"  [PASS] Manuscript step count appears updated")

    # Check if step_10 CMB/LSS is referenced
    if "step_10" in text or "cmb anisotropy" in text or "class boltzmann" in text:
        print(f"  [PASS] Manuscript references CMB/LSS checks")
    else:
        print(f"  [WARN] Manuscript does not reference step_10 CMB/LSS checks")
        AUDIT_WARNINGS.append("Manuscript missing step_10 CMB/LSS references")


def main() -> int:
    audit_result_files()
    audit_codebase()
    audit_physical_numbers()
    audit_manuscript()

    print("\n" + "=" * 70)
    print("AUDIT SUMMARY")
    print("=" * 70)
    print(f"Issues:   {len(AUDIT_ISSUES)}")
    print(f"Warnings: {len(AUDIT_WARNINGS)}")
    if AUDIT_ISSUES:
        print("\nIssues:")
        for issue in AUDIT_ISSUES:
            print(f"  - {issue}")
    if AUDIT_WARNINGS:
        print("\nWarnings:")
        for warning in AUDIT_WARNINGS:
            print(f"  - {warning}")
    if not AUDIT_ISSUES and not AUDIT_WARNINGS:
        print("\nSTATUS: ALL CHECKS PASS")
    elif not AUDIT_ISSUES:
        print("\nSTATUS: PASS WITH WARNINGS")
    else:
        print("\nSTATUS: FAIL — issues require attention")
    print("=" * 70)
    return 1 if AUDIT_ISSUES else 0


if __name__ == "__main__":
    sys.exit(main())
