#!/usr/bin/env python3
"""TEP-TH temporal-horizon pipeline runner.

Runs the eleven-step pipeline that validates the temporal-horizon cosmology.
Each step writes JSON and CSV results to results/.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent / "steps"))

from th_common import step_json_path

PIPELINE_STEPS = [
    ("step_00_temporal_horizon_mapping", "Temporal-horizon mapping", []),
    ("step_01_matter_frame_curvature", "Matter-frame curvature", []),
    ("step_02_geodesic_completeness", "Geodesic completeness", []),
    ("step_03_effective_stress_energy", "Effective stress-energy", []),
    ("step_04_full_bbn_abundances", "BBN abundance validation", []),
    ("step_05_recombination_visibility", "Recombination visibility", []),
    ("step_06_cmb_blackbody_origin", "CMB blackbody origin", []),
    ("step_07_entropy_arrow", "Entropy and arrow of time", []),
    ("step_08_primordial_perturbation_boundary", "Primordial perturbation boundary", []),
    ("step_09_temporal_horizon_claim_gate", "Temporal-horizon claim gate", []),
    ("step_09b_native_tensor_integration", "Native tensor-mode integration", []),
    ("step_10_cmb_lss_class", "CMB and LSS consistency", []),
]


def run_step(step_module: str) -> dict:
    try:
        module = __import__(step_module)
        if not hasattr(module, "run"):
            return {"step": step_module, "error": "No run() function"}
        return module.run()
    except Exception as exc:
        return {"step": step_module, "error": str(exc), "traceback": traceback.format_exc()}


def selected_steps(specific_steps: list[str] | None) -> list[tuple[str, str, list[str]]]:
    if specific_steps is None:
        return PIPELINE_STEPS
    requested = set(specific_steps)
    step_deps = {name: deps for name, _, deps in PIPELINE_STEPS}
    known = set(step_deps)
    unknown = sorted(requested - known)
    if unknown:
        raise SystemExit(f"Unknown step(s): {', '.join(unknown)}")

    expanded = set(requested)
    changed = True
    while changed:
        changed = False
        for name in list(expanded):
            for dep in step_deps[name]:
                if dep not in expanded:
                    expanded.add(dep)
                    changed = True
    return [step for step in PIPELINE_STEPS if step[0] in expanded]


def run_pipeline(specific_steps: list[str] | None = None, resume: bool = False) -> dict:
    print("=" * 70)
    print("TEP-TH TEMPORAL-HORIZON PIPELINE")
    print("=" * 70)
    print()

    steps_to_run = selected_steps(specific_steps)
    results: dict[str, dict] = {}
    completed: set[str] = set()
    failed: list[str] = []
    skipped: list[str] = []

    for step_module, description, deps in steps_to_run:
        if resume and step_json_path(step_module).exists():
            print(f"[SKIP] {step_module}: existing result found")
            results[step_module] = {
                "step": step_module,
                "status": "skipped_existing_result",
                "result_path": str(step_json_path(step_module)),
            }
            completed.add(step_module)
            continue

        missing = [dep for dep in deps if dep not in completed]
        if missing:
            print(f"[SKIP] {step_module}: missing dependencies {missing}")
            results[step_module] = {
                "step": step_module,
                "error": "missing dependencies",
                "missing_dependencies": missing,
            }
            skipped.append(step_module)
            continue

        print(f"[RUN] {step_module}: {description}")
        start = time.time()
        result = run_step(step_module)
        elapsed = time.time() - start

        if result is None:
            result = {"step": step_module, "status": "failed", "error": "Step returned None (missing return value)"}
        results[step_module] = result

        if isinstance(result, dict) and "error" in result:
            print(f"  [FAIL] {elapsed:.1f}s: {result['error']}")
            failed.append(step_module)
        else:
            print(f"  [OK] {elapsed:.1f}s")
            completed.add(step_module)

    print()
    print("=" * 70)
    print(f"PIPELINE COMPLETE: {len(completed)}/{len(steps_to_run)} steps succeeded")
    if failed:
        print(f"FAILED: {', '.join(failed)}")
    if skipped:
        print(f"SKIPPED: {', '.join(skipped)}")
    print("=" * 70)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description="TEP-TH Pipeline")
    parser.add_argument("--steps", nargs="+", help="Specific step module names to run")
    parser.add_argument("--resume", action="store_true", help="Skip steps whose result JSON already exists")
    args = parser.parse_args()

    results = run_pipeline(args.steps, resume=args.resume)
    return 1 if any("error" in result for result in results.values()) else 0


if __name__ == "__main__":
    sys.exit(main())
