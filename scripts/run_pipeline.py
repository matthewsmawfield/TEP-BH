#!/usr/bin/env python3
"""TEP-BH canonical pipeline runner — single entry point.

Orchestrates the full black-hole analysis pipeline for Paper 28 (Bahrain).
All scripts are numbered step_00 through step_49 in scripts/steps/:

  CORE      00–11  Data download, fixed-background metric, perturbations,
                   ray-tracing, accretion, observational constraints,
                   GW/QPO/JWST/TDE/EHT confrontations.
  COUPLED   12–15  Self-gravitating sGB, scalar perturbations, Kerr-TEP,
                   dynamical signatures.
  GEOM      16–19  Exact geometry, null expansions, observer redshift,
                   frequency transfer.
  OBS       20–24  Corrected observables, cT derivation, quadratic action,
                   characteristic matrices, conformal invariance.
  QNM       25–30  QNM solvers: perturbative, spectral, horizon branch,
                   validation, matrix Leaver, Taylor recurrence.
  INT       31–34,50  Interior: Frobenius, integration, analysis, solver,
                   darkness threshold.
  PHAN      35–39  Phantom mass, EHT visibility fit.
  KGW       40–41  Kerr-sGB true solution, GW250114 confrontation.
  FIG       42     Figure generation.
  INFER     43–49  S-star inference: data, GR fit, mass-bias, TEP transfer,
                   joint forward-model, likelihood, MCMC posterior.

Usage:
    python scripts/run_pipeline.py                       # Run everything
    python scripts/run_pipeline.py --start-step step_25   # Start from a step
    python scripts/run_pipeline.py --stop-step qnm        # Stop after QNM phase
    python scripts/run_pipeline.py --skip-steps step_02  # Skip specific steps
    python scripts/run_pipeline.py --no-derive           # Skip all derive phases
    python scripts/run_pipeline.py --no-inference        # Skip inference phase
    python scripts/run_pipeline.py --no-figures          # Skip figure generation
    python scripts/run_pipeline.py --list-steps          # List all steps

Author: Matthew Lukin Smawfield
Version: TEP-BH v0.2 (Bahrain)
License: CC-BY-4.0
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import runpy
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STEPS_DIR = PROJECT_ROOT / "scripts" / "steps"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(STEPS_DIR))

from scripts.utils.logger import TEPLogger, set_step_logger, print_status
from bh_common import ensure_dirs, write_json, rel, RESULTS_DIR, LOGS_DIR

# ---------------------------------------------------------------------------
# Phase registries — each entry is (step_id, description, deps)
# ---------------------------------------------------------------------------

# Phase 1: Core — Data & Fixed-Background (00–11)
CORE_STEPS = [
    ("step_00_data_download", "Data Download", []),
    ("step_01_field_equations", "Strong-Field Metric Construction", ["step_00_data_download"]),
    ("step_02_perturbations", "Perturbation Analysis", ["step_01_field_equations"]),
    ("step_03_raytracing", "Ray-Tracing & Observational Comparison", ["step_01_field_equations"]),
    ("step_04_accretion", "Accretion Dynamics & ISCO", ["step_01_field_equations"]),
    ("step_05_observational_constraints", "Observational Constraints", ["step_00_data_download", "step_01_field_equations", "step_02_perturbations", "step_03_raytracing"]),
    ("step_06_gw190521_mass_gap", "GW190521 Mass Gap — Exterior Consistency Check", ["step_00_data_download"]),
    ("step_07_qpo_frequency_lock", "QPO Frequency Lock — Epicyclic Resonance", ["step_04_accretion"]),
    ("step_08_jwst_early_smbhs", "JWST Early SMBHs — Eddington-limited Growth Analysis", []),
    ("step_09_spin_bias", "Spin Bias Analysis", ["step_03_raytracing"]),
    ("step_10_tde_missing_flares", "TDE Missing Flares — Sub-Eddington TDE Luminosities", []),
    ("step_11_eht_polarization", "EHT Polarization — Exterior-Null Result", ["step_00_data_download", "step_03_raytracing"]),
]

# Phase 2: Core — Coupled Solution (12–15)
COUPLED_STEPS = [
    ("step_12_self_gravitating", "Self-Gravitating Einstein-Scalar-Gauss-Bonnet Solution", ["step_01_field_equations"]),
    ("step_13_scalar_perturbations", "Scalar-Led & Mixed Scalar-Tensor Perturbations", ["step_12_self_gravitating"]),
    ("step_14_kerr_tep", "Kerr-TEP Generalization with Frame-Dragging", ["step_12_self_gravitating"]),
    ("step_15_dynamical_signatures", "Dynamical Signatures — PN, Inspiral, Cosmology", ["step_12_self_gravitating"]),
]

# Phase 3: Derive — Geometry & Invariants (16–19)
GEOMETRY_STEPS = [
    ("step_16_exact_geometry", "Exact Geometry, Tensors & Invariants", ["step_01_field_equations"]),
    ("step_17_null_expansions", "Null Expansions & Geodesic Structure", ["step_01_field_equations"]),
    ("step_18_observer_redshift", "Observer Redshift & Clock-Rate Profile", ["step_01_field_equations"]),
    ("step_19_observer_frequency_transfer", "Observer Frequency Transfer — Invariant Test", ["step_01_field_equations"]),
]

# Phase 4: Derive — Corrected Observables & Fundamentals (20–24)
OBSERVABLE_STEPS = [
    ("step_20_corrected_observables", "Corrected Exterior Observables (Fixed-ADM)", ["step_12_self_gravitating"]),
    ("step_21_cT_derivation", "Tensor Speed c_T Derivation", []),
    ("step_22_quadratic_action", "Quadratic Action Derivation", ["step_12_self_gravitating"]),
    ("step_23_characteristic_matrices", "Tensor Characteristic Matrices", []),
    ("step_24_conformal_invariance_check", "Conformal Invariance Check", ["step_34_solve_interior"]),
]

# Phase 5: Derive — QNM Solvers (25–30)
QNM_STEPS = [
    ("step_25_qnm_solver", "QNM Solver (Perturbative sGB)", ["step_20_corrected_observables"]),
    ("step_26_coupled_spectral_solver", "Coupled Spectral Solver", ["step_25_qnm_solver"]),
    ("step_27_qnm_horizon_branch", "QNM Horizon Branch", ["step_25_qnm_solver"]),
    ("step_28_qnm_validation", "QNM Validation", ["step_27_qnm_horizon_branch"]),
    ("step_29_matrix_leaver", "Matrix Leaver Solver", ["step_34_solve_interior"]),
    ("step_30_taylor_recurrence", "Taylor Recurrence Solver", ["step_34_solve_interior"]),
]

# Phase 6: Derive — Interior (31–34, 50)
INTERIOR_STEPS = [
    ("step_31_frobenius_analysis", "Frobenius Analysis", ["step_12_self_gravitating"]),
    ("step_32_interior_integration", "Interior Integration", ["step_31_frobenius_analysis"]),
    ("step_33_interior_analysis", "Interior Analysis", ["step_32_interior_integration"]),
    ("step_34_solve_interior", "TEP Interior Solver on Hayward Background", ["step_32_interior_integration"]),
    ("step_50_darkness_threshold", "Operational Darkness Threshold", ["step_34_solve_interior"]),
]

# Phase 7: Derive — Phantom Mass & EHT (35–39)
PHANTOM_STEPS = [
    ("step_35_mass_bias_sign", "Mass-Bias Sign Equation", ["step_01_field_equations"]),
    ("step_36_phantom_mass_critical", "Phantom Mass Critical Scale", ["step_20_corrected_observables", "step_34_solve_interior"]),
    ("step_37_phantom_mass_raytrace", "Phantom Mass Ray-Tracing", ["step_36_phantom_mass_critical"]),
    ("step_38_phantom_mass_scan", "Phantom Mass Scan", ["step_36_phantom_mass_critical"]),
    ("step_39_eht_visibility_fit", "EHT Visibility-Domain Fit", ["step_20_corrected_observables"]),
]

# Phase 8: Derive — Kerr & GW (40–41)
KERR_GW_STEPS = [
    ("step_40_kerr_sgb_true", "Kerr-sGB True Solution", ["step_14_kerr_tep"]),
    ("step_41_gw250114_confrontation", "GW250114 Corrected Confrontation", ["step_20_corrected_observables"]),
]

# Phase 9: Figures (42)
FIGURE_STEPS = [
    ("step_42_generate_figures", "Figure Generation", ["step_05_observational_constraints"]),
]

# Phase 10: S-star Inference (43–49)
INFERENCE_STEPS = [
    ("step_43_sstar_data", "S-star Data Acquisition", []),
    ("step_44_gr_fit", "Conventional GR Fit", ["step_43_sstar_data"]),
    ("step_45_mass_bias", "Mass-Bias Sign (Gate -1)", ["step_44_gr_fit"]),
    ("step_46_tep_transfer", "TEP Transfer-Function Fit", ["step_44_gr_fit"]),
    ("step_47_joint_forward", "Joint Forward-Model", ["step_46_tep_transfer"]),
    ("step_48_likelihood", "Likelihood Comparison", ["step_46_tep_transfer"]),
    ("step_49_posterior", "Phantom Mass Posterior (MCMC)", ["step_48_likelihood"]),
]

# Phase boundaries for --start-step / --stop-step by phase name
PHASE_BOUNDARIES = {
    "core": ("step_00_data_download", "step_11_eht_polarization"),
    "coupled": ("step_12_self_gravitating", "step_15_dynamical_signatures"),
    "geometry": ("step_16_exact_geometry", "step_19_observer_frequency_transfer"),
    "observables": ("step_20_corrected_observables", "step_24_conformal_invariance_check"),
    "qnm": ("step_25_qnm_solver", "step_30_taylor_recurrence"),
    "interior": ("step_31_frobenius_analysis", "step_50_darkness_threshold"),
    "phantom": ("step_35_mass_bias_sign", "step_39_eht_visibility_fit"),
    "kerr_gw": ("step_40_kerr_sgb_true", "step_41_gw250114_confrontation"),
    "figures": ("step_42_generate_figures", "step_42_generate_figures"),
    "inference": ("step_43_sstar_data", "step_49_posterior"),
}

ALL_PHASES = [
    ("CORE", "Core — Data & Fixed-Background", CORE_STEPS),
    ("COUPLED", "Core — Coupled Solution", COUPLED_STEPS),
    ("GEOM", "Derive — Geometry & Invariants", GEOMETRY_STEPS),
    ("OBS", "Derive — Corrected Observables", OBSERVABLE_STEPS),
    ("QNM", "Derive — QNM Solvers", QNM_STEPS),
    ("INT", "Derive — Interior", INTERIOR_STEPS),
    ("PHAN", "Derive — Phantom Mass & EHT", PHANTOM_STEPS),
    ("KGW", "Derive — Kerr & GW", KERR_GW_STEPS),
    ("FIG", "Figures", FIGURE_STEPS),
    ("INFER", "S-star Inference", INFERENCE_STEPS),
]


def _all_step_ids():
    ids = []
    for _, _, steps in ALL_PHASES:
        for sid, _, _ in steps:
            ids.append(sid)
    return ids


def _find_step_index(step_id: str):
    idx = 0
    for _, _, steps in ALL_PHASES:
        for sid, _, _ in steps:
            if sid == step_id:
                return idx
            idx += 1
    raise ValueError(f"Unknown step: {step_id}")


def _resolve_phase_boundary(name: str):
    start, stop = PHASE_BOUNDARIES[name]
    return _find_step_index(start), _find_step_index(stop)


def import_step(step_id: str):
    """Import a step module and return its main callable.

    All steps are now named step_NN_*.py in scripts/steps/.
    Most have main(); a few script-style modules are wrapped with runpy.
    """
    script_path = STEPS_DIR / f"{step_id}.py"
    if not script_path.exists():
        raise ValueError(f"Step script not found: {script_path}")
    try:
        mod = importlib.import_module(step_id)
        if hasattr(mod, "main"):
            return mod.main
    except Exception:
        pass
    # Script-style: wrap with runpy
    def _runpy_wrapper():
        runpy.run_path(str(script_path), run_name="__main__")
        return {"status": "success"}
    return _runpy_wrapper


# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------

def run_pipeline(args: argparse.Namespace) -> dict:
    ensure_dirs()

    pipeline_logger = TEPLogger(
        "tep_bh_pipeline",
        log_file_path=LOGS_DIR / "tep_bh_pipeline.log",
        reset_log=True,
        prefix="PIPELINE",
    )
    set_step_logger(pipeline_logger)

    all_ids = _all_step_ids()
    total_steps = len(all_ids)

    # Determine step range
    start_idx = 0
    stop_idx = total_steps - 1

    if args.start_step:
        if args.start_step in PHASE_BOUNDARIES:
            start_idx, _ = _resolve_phase_boundary(args.start_step)
        else:
            start_idx = _find_step_index(args.start_step)
    if args.stop_step:
        if args.stop_step in PHASE_BOUNDARIES:
            _, stop_idx = _resolve_phase_boundary(args.stop_step)
        else:
            stop_idx = _find_step_index(args.stop_step)

    skip_steps = set(args.skip_steps) if args.skip_steps else set()

    # Determine which phases to skip
    if args.no_derive:
        for _, _, steps in [GEOMETRY_STEPS, OBSERVABLE_STEPS, QNM_STEPS, INTERIOR_STEPS, PHANTOM_STEPS, KERR_GW_STEPS]:
            skip_steps.update(sid for sid, _, _ in steps)
    if args.no_interior:
        skip_steps.update(sid for sid, _, _ in INTERIOR_STEPS)
    if args.no_inference:
        skip_steps.update(sid for sid, _, _ in INFERENCE_STEPS)
    if args.no_figures:
        skip_steps.update(sid for sid, _, _ in FIGURE_STEPS)

    print_status("=" * 70, "TITLE")
    print_status("TEP-BH ANALYSIS PIPELINE", "TITLE")
    print_status("Paper 28 — Bahrain — v0.2", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(f"Project root: {PROJECT_ROOT}", "INFO")
    print_status(f"Started: {datetime.now().isoformat()}", "INFO")
    print_status(f"Python: {sys.version.split()[0]}", "INFO")
    print_status("")

    print_status("Pipeline Configuration:", "TITLE")
    print_status(f"  Start: {all_ids[start_idx]}", "INFO")
    print_status(f"  Stop:  {all_ids[stop_idx]}", "INFO")
    print_status(f"  Total steps: {stop_idx - start_idx + 1}", "INFO")
    if skip_steps:
        print_status(f"  Skip:  {sorted(skip_steps)}", "INFO")
    print_status("")

    results = {
        "pipeline_start": datetime.now().isoformat(),
        "steps_completed": [],
        "steps_failed": [],
        "steps_skipped": [],
        "execution_times": {},
        "status": "RUNNING",
    }

    total_start = time.time()
    global_idx = 0

    for phase_tag, phase_name, phase_steps in ALL_PHASES:
        for step_id, step_desc, deps in phase_steps:
            if global_idx < start_idx or global_idx > stop_idx:
                global_idx += 1
                continue
            global_idx += 1

            if step_id in skip_steps:
                print_status(f"[{phase_tag:6s}] {step_id:40s} SKIPPED", "WARNING")
                results["steps_skipped"].append(step_id)
                continue

            print_status(f"[{phase_tag:6s}] {step_id}: {step_desc.upper()}", "TITLE")
            print_status(f"  Dependencies: {deps if deps else 'none'}", "DEBUG")
            print_status("")

            # Set up a step-specific logger with prefix so all output is tagged
            step_logger = TEPLogger(
                name=f"tep_bh_{step_id}",
                log_file_path=LOGS_DIR / f"{step_id}.log",
                reset_log=True,
                prefix=step_id,
            )
            # Tee to pipeline log
            pipeline_log_path = LOGS_DIR / "tep_bh_pipeline.log"
            if pipeline_log_path.exists():
                import logging as _logging
                from scripts.utils.logger import TEPFileFormatter
                tee_fh = _logging.FileHandler(pipeline_log_path, mode='a', encoding='utf-8')
                tee_fh.setLevel(_logging.DEBUG)
                tee_fh.setFormatter(TEPFileFormatter(prefix=step_id))
                step_logger.logger.addHandler(tee_fh)
            set_step_logger(step_logger)

            step_start = time.time()
            try:
                step_main = import_step(step_id)
                step_result = step_main()
                results["steps_completed"].append(step_id)
                step_time = time.time() - step_start
                results["execution_times"][step_id] = round(step_time, 2)
                print_status(f"[{phase_tag:6s}] {step_id} completed in {step_time:.1f}s", "SUCCESS")
                if isinstance(step_result, dict):
                    kr = step_result.get("key_result", "")
                    desc = step_result.get("description", "")
                    if desc:
                        print_status(f"  Description: {desc}", "INFO")
                    if kr:
                        print_status(f"  Key result:  {kr}", "INFO")
                    json_keys = len(step_result)
                    print_status(f"  JSON output: {json_keys} keys", "DEBUG")
            except Exception as e:
                results["steps_failed"].append(step_id)
                step_time = time.time() - step_start
                results["execution_times"][step_id] = round(step_time, 2)
                print_status(f"[{phase_tag:6s}] {step_id} FAILED after {step_time:.1f}s", "ERROR")
                print_status(f"  Error: {e}", "ERROR")
                print_status(f"  Exception: {type(e).__name__}", "ERROR")
                traceback.print_exc()
                if not args.continue_on_error:
                    print_status("Pipeline halted. Use --continue-on-error to proceed.", "WARNING")
                    results["status"] = "FAILED"
                    break

            # Restore pipeline logger for pipeline-level messages
            set_step_logger(pipeline_logger)

            print_status("")
            print_status("-" * 70, "TITLE")
            print_status("")

        if results["status"] == "FAILED":
            break

    # Summary
    total_time = time.time() - total_start
    results["pipeline_end"] = datetime.now().isoformat()
    results["total_time_seconds"] = round(total_time, 2)

    print_status("=" * 70, "TITLE")
    print_status("PIPELINE EXECUTION SUMMARY", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status(f"  Steps completed: {len(results['steps_completed'])}", "INFO")
    print_status(f"  Steps skipped:   {len(results['steps_skipped'])}", "INFO")
    print_status(f"  Steps failed:    {len(results['steps_failed'])}", "INFO")
    print_status(f"  Total time:      {total_time:.1f}s", "INFO")
    print_status(f"  Avg per step:    {total_time/max(len(results['steps_completed']),1):.1f}s", "INFO")
    print_status("")

    if results["steps_completed"]:
        print_status("Completed Steps:", "TITLE")
        for sid in results["steps_completed"]:
            elapsed = results["execution_times"].get(sid, 0)
            print_status(f"  [OK]   {sid:40s} {elapsed:6.1f}s", "SUCCESS")
        print_status("")

    if results["steps_skipped"]:
        print_status("Skipped Steps:", "TITLE")
        for sid in results["steps_skipped"]:
            print_status(f"  [SKIP] {sid}", "WARNING")
        print_status("")

    if results["steps_failed"]:
        print_status("Failed Steps:", "TITLE")
        for sid in results["steps_failed"]:
            elapsed = results["execution_times"].get(sid, 0)
            print_status(f"  [FAIL] {sid:40s} {elapsed:6.1f}s", "ERROR")
        print_status("")

    if results["execution_times"]:
        print_status("Step Timing Breakdown (sorted by duration):", "TITLE")
        for sid, elapsed in sorted(results["execution_times"].items(), key=lambda x: -x[1]):
            bar = "#" * int(elapsed / max(results["execution_times"].values()) * 30)
            print_status(f"  {sid:40s} {elapsed:6.1f}s {bar}", "INFO")
    print_status("")

    if results["steps_failed"]:
        results["status"] = "COMPLETED_WITH_ERRORS" if args.continue_on_error else "FAILED"
        print_status(f"Status: {results['status']}", "WARNING")
    else:
        results["status"] = "SUCCESS"
        print_status("Status: ALL STEPS COMPLETED SUCCESSFULLY", "SUCCESS")

    print_status("")
    print_status(f"Results directory: {RESULTS_DIR}", "INFO")
    print_status(f"Logs directory:    {LOGS_DIR}", "INFO")
    print_status("")

    print_status("Output Files:", "TITLE")
    for path in sorted(RESULTS_DIR.glob("*.json")):
        size_kb = path.stat().st_size / 1024
        print_status(f"  {path.name:45s} {size_kb:8.1f} KB", "INFO")
    for path in sorted(RESULTS_DIR.glob("*.csv")):
        size_kb = path.stat().st_size / 1024
        print_status(f"  {path.name:45s} {size_kb:8.1f} KB", "INFO")
    print_status("")

    print_status("Log Files:", "TITLE")
    for path in sorted(LOGS_DIR.glob("*.log")):
        size_kb = path.stat().st_size / 1024
        n_lines = sum(1 for _ in open(path))
        print_status(f"  {path.name:45s} {size_kb:8.1f} KB ({n_lines} lines)", "INFO")
    print_status("")

    print_status("=" * 70, "TITLE")

    results_file = RESULTS_DIR / "pipeline_results.json"
    write_json(results_file, results)
    print_status(f"Pipeline results saved to: {rel(results_file)}", "INFO")

    checksums = {}
    for path in sorted(RESULTS_DIR.glob("*.json")):
        checksums[str(path.name)] = hashlib.sha256(path.read_bytes()).hexdigest()
    for path in sorted(RESULTS_DIR.glob("*.csv")):
        checksums[str(path.name)] = hashlib.sha256(path.read_bytes()).hexdigest()
    checksums_file = RESULTS_DIR / "checksums_sha256.json"
    write_json(checksums_file, checksums)
    print_status(
        f"SHA-256 checksums for {len(checksums)} result files "
        f"saved to: {rel(checksums_file)}",
        "INFO",
    )

    return results


def main():
    parser = argparse.ArgumentParser(
        description="TEP-BH Analysis Pipeline — Black-Hole Temporal Horizon",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=""\
"""Phases (step_NN):
  core       00-11  Data & fixed-background
  coupled    12-15  Self-grav sGB, Kerr-TEP, dynamical
  geometry   16-19  Exact geometry, null expansions, redshift
  observables 20-24 Corrected observables, cT, quadratic action
  qnm        25-30  QNM solvers
  interior   31-34,50  Interior analysis, solver & darkness threshold
  phantom    35-39  Phantom mass, EHT visibility
  kerr_gw    40-41  Kerr-sGB, GW250114
  figures    42     Figure generation
  inference  43-49  S-star inference

Examples:
  python scripts/run_pipeline.py                       # Run everything
  python scripts/run_pipeline.py --start-step step_25   # Start from QNM
  python scripts/run_pipeline.py --stop-step qnm        # Stop after QNM phase
  python scripts/run_pipeline.py --no-derive --no-inference  # Core only
  python scripts/run_pipeline.py --list-steps          # List all steps
""",
    )
    parser.add_argument("--start-step", type=str, default=None,
                        help="First step to execute (step_id or phase name)")
    parser.add_argument("--stop-step", type=str, default=None,
                        help="Last step to execute (step_id or phase name)")
    parser.add_argument("--skip-steps", type=lambda s: [x.strip() for x in s.split(",")],
                        default=None, help="Comma-separated step_ids to skip")
    parser.add_argument("--continue-on-error", action="store_true",
                        help="Continue pipeline even if a step fails")
    parser.add_argument("--no-derive", action="store_true",
                        help="Skip derive phase")
    parser.add_argument("--no-interior", action="store_true",
                        help="Skip interior solver phase")
    parser.add_argument("--no-inference", action="store_true",
                        help="Skip S-star inference phase")
    parser.add_argument("--no-figures", action="store_true",
                        help="Skip figure generation phase")
    parser.add_argument("--list-steps", action="store_true",
                        help="List all pipeline steps and exit")
    args = parser.parse_args()

    if args.list_steps:
        print("TEP-BH Pipeline Steps:")
        print("=" * 70)
        for phase_tag, phase_name, phase_steps in ALL_PHASES:
            print(f"\n  [{phase_tag}] {phase_name}:")
            for sid, desc, deps in phase_steps:
                dep_str = f" (requires: {', '.join(deps)})" if deps else ""
                print(f"    {sid:40s} - {desc}{dep_str}")
        return

    try:
        results = run_pipeline(args)
        if results["status"] == "SUCCESS":
            sys.exit(0)
        elif results["status"] == "COMPLETED_WITH_ERRORS":
            sys.exit(1)
        else:
            sys.exit(2)
    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user.")
        sys.exit(130)
    except Exception as e:
        print(f"\n\nPipeline crashed: {e}")
        traceback.print_exc()
        sys.exit(3)


if __name__ == "__main__":
    main()
