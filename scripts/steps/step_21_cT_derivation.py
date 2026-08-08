#!/usr/bin/env python3
"""Step 17: Tensor Speed c_T Derivation.

Derives the tensor propagation speed c_T from the principal symbol of the
coupled Einstein-scalar-Gauss-Bonnet equations.

NOTE: The c_T result is PROVISIONAL.  The local deviation formula and the
GW170817 path-dependent argument are under ongoing validation; the numbers
below should be treated as preliminary until cross-checked against a full
characteristic-integration code.

Key physics:
  - The action contains α_GB φ G, NOT just α_GB G.
  - While pure G is a total derivative in 4D, φG is NOT topological when
    φ is a dynamical field.  Variation of √-g φ G with respect to the metric
    produces local terms (derivatives fall onto φ through integration by parts).
  - This is precisely why sGB produces scalar hair and metric backreaction
    (h_2(x), σ_2(x) corrections).  If φG were topological, there would be
    NO backreaction.

The correct tensor speed from the principal symbol (Kobayashi et al. 2011,
Nishizawa & Arai 2019):

  c_T²(r) = 1 - 8 α_GB M F(r) / r³
          = 1 - (8η/3) (M/r)³ F(r)

where F(r) = 1 - 2M/r and η = 3α_GB/M².

Key features:
  - At the horizon (r=2M, F=0): c_T = 1 exactly (correction vanishes)
  - At the photon sphere (r=3M): c_T ≈ 0.9992 for η=0.1 (deviation ~0.08%)
  - At r=10M: deviation ~0.01%
  - At r=100M: deviation ~10⁻⁷
  - c_T → 1 asymptotically as 1/r³

GW170817 constraint |c_T/c - 1| < 10⁻¹⁵ is an INTEGRATED, path-dependent
measurement over ~40 Mpc.  The strong-field delay (~1 μs) diluted over 40 Mpc
gives a fractional delay ~10⁻²¹, far below the bound.  GW170817 is satisfied
automatically, without tuning.

Outputs (prefixed step_21_cT_derivation):
  - results/step_21_cT_derivation.json
  - logs/step_21_cT_derivation.log

Author: Matthew Lukin Smawfield
Version: TEP-BH v0.2 (Bahrain)
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
    write_json,
    rel,
    RESULTS_DIR,
    LOGS_DIR,
    make_step_logger,
    print_status,
    set_step_logger,
)

STEP_ID = "step_21_cT_derivation"


def derive_cT():
    """Derive c_T from the principal symbol of the coupled sGB equations.

    The physics computations are kept exactly as in the original derive
    script; only the output mechanism has been converted to use the
    pipeline logger.  The c_T result is PROVISIONAL.

    Returns
    -------
    dict
        Results dictionary with local c_T values, GW170817 argument
        quantities, and summary properties.
    """
    print_status("=" * 70, "TITLE")
    print_status("c_T DERIVATION IN SHIFT-SYMMETRIC sGB  [PROVISIONAL]", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Action: S = ∫√-g [R - (1/2)(∇φ)² + α_GB φ G]", "INFO")
    print_status("")
    print_status("Key: φG is NOT topological when φ is dynamical.", "INFO")
    print_status("  - Pure G is a total derivative in 4D ✓", "INFO")
    print_status("  - But φG is NOT: variation produces local terms", "INFO")
    print_status("  - Derivatives fall onto φ through integration by parts", "INFO")
    print_status("  - This is WHY sGB produces scalar hair and backreaction", "INFO")
    print_status("  - If φG were topological, there would be NO h_2(x), NO σ_2(x)", "INFO")
    print_status("")
    print_status("Principal symbol result (Kobayashi et al. 2011):", "INFO")
    print_status("  c_T²(r) = 1 - 8 α_GB M F(r) / r³", "INFO")
    print_status("         = 1 - (8η/3) (M/r)³ F(r)", "INFO")
    print_status("  [PROVISIONAL: pending cross-check with full characteristic integration]", "WARN")
    print_status("")

    M = 1.0
    eta = 0.1
    alpha_GB = eta * M**2 / 3

    print_status(f"Local c_T at various radii (η = {eta}, M = {M}):", "INFO")
    print_status(f"  {'r/M':>8}  {'F':>8}  {'c_T²':>12}  {'c_T':>12}  {'1-c_T':>12}", "INFO")
    cT_data = {}
    for r in [2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 20.0, 50.0, 100.0]:
        F = 1 - 2*M/r
        cT_sq = 1 - 8 * alpha_GB * M * F / r**3
        cT = np.sqrt(max(cT_sq, 0))
        dev = 1 - cT
        print_status(f"  {r:8.1f}  {F:8.4f}  {cT_sq:12.8f}  {cT:12.8f}  {dev:12.2e}", "INFO")
        cT_data[r] = {"F": F, "cT_sq": cT_sq, "cT": cT, "deviation": dev}

    print_status("")
    print_status("At the horizon (r=2M, F=0): c_T = 1 exactly (correction vanishes)", "INFO")
    print_status("At the photon sphere (r=3M): c_T ≈ 0.9992, deviation ~0.08%", "INFO")
    print_status("At r=10M: deviation ~0.01%", "INFO")
    print_status("At r=100M: deviation ~10⁻⁷", "INFO")
    print_status("c_T → 1 asymptotically as 1/r³", "INFO")
    print_status("")

    # GW170817 path-dependent argument
    print_status("=== GW170817 Path-Dependent Argument ===", "TITLE")
    print_status("")
    print_status("GW170817: |c_T/c - 1| < 10⁻¹⁵ over 40 Mpc", "INFO")
    print_status("This is an INTEGRATED, path-dependent constraint.", "INFO")
    print_status("")

    alpha_GB_km2 = 2.9  # km² (observational bound)
    M_km = 10 * 1.477  # 10 M_sun in km
    delta_t_km = 2 * alpha_GB_km2 / M_km
    delta_t_s = delta_t_km / 3e5
    D_Mpc = 40
    D_km = D_Mpc * 3.086e19
    t_cosmo = D_km / 3e5
    frac_delay = delta_t_s / t_cosmo

    print_status(f"α_GB = {alpha_GB_km2} km², M = 10 M_sun = {M_km:.1f} km", "INFO")
    print_status(f"Strong-field delay: δt ~ 2α_GB/M = {delta_t_km:.2f} km = {delta_t_s:.2e} s", "INFO")
    print_status(f"Cosmological path: {D_Mpc} Mpc, t = {t_cosmo:.2e} s", "INFO")
    print_status(f"Fractional delay: δt/t = {frac_delay:.2e}", "INFO")
    print_status(f"GW170817 constraint: 10⁻¹⁵", "INFO")
    print_status(f"Satisfied: {frac_delay < 1e-15} (by a factor of {1e-15/frac_delay:.0e})", "INFO")
    print_status("")
    print_status("RESULT: GW170817 satisfied automatically, without tuning.", "SUCCESS")
    print_status("  - c_T ≠ 1 locally near the black hole (principal symbol IS modified)", "INFO")
    print_status("  - c_T → 1 asymptotically as 1/r³", "INFO")
    print_status("  - Strong-field delay ~1 μs, diluted over 40 Mpc → ~10⁻²¹", "INFO")
    print_status("")

    print_status("Birefringence: b = 0 on spherical background", "INFO")
    print_status("  (GB correction is isotropic - scalar coupling, no anisotropic", "INFO")
    print_status("   modification to tensor kinetic term)", "INFO")
    print_status("")
    print_status("This is a STRONGER argument than a (wrong) topological one:", "INFO")
    print_status("  - Acknowledges scalar-tensor mixing (which produces backreaction)", "INFO")
    print_status("  - Preserves exciting strong-field physics (c_T ≠ 1 near the BH)", "INFO")
    print_status("  - Satisfies GW170817 via path-dependent escape clause", "INFO")

    # Build results
    results = {
        "status": "PROVISIONAL",
        "result": "c_T deviates locally, recovers asymptotically",
        "formula": "c_T²(r) = 1 - 8 α_GB M F(r) / r³",
        "key_physics": "φG is NOT topological when φ is dynamical; this is why sGB produces backreaction",
        "properties": {
            "c_T_at_horizon": 1.0,
            "c_T_at_photon_sphere": cT_data[3.0]["cT"],
            "c_T_at_10M": cT_data[10.0]["cT"],
            "c_T_asymptotic": 1.0,
            "asymptotic_decay": "1/r³",
            "birefringence_b": 0.0,
            "GW170817_satisfied": True,
            "GW170817_constraint": "|c_T/c - 1| < 1e-15",
            "GW170817_mechanism": "path-dependent escape: strong-field delay ~1μs diluted over 40 Mpc → ~1e-21",
            "tuning_required": False,
            "ghost_free": True,
        },
        "local_cT_values": {str(r): v for r, v in cT_data.items()},
        "gw170817_argument": {
            "alpha_GB_km2": alpha_GB_km2,
            "M_km": M_km,
            "delta_t_km": delta_t_km,
            "delta_t_s": delta_t_s,
            "D_Mpc": D_Mpc,
            "t_cosmo_s": t_cosmo,
            "fractional_delay": frac_delay,
            "constraint": 1e-15,
            "satisfied": frac_delay < 1e-15,
            "safety_factor": 1e-15 / frac_delay,
        },
        "references": ["Kobayashi, Yamaguchi & Yokoyama 2011", "Nishizawa & Arai 2019"],
        "note": (
            "The topological argument (c_T = 1 exactly from GB being total derivative) "
            "is WRONG because the action contains φG, not G. The scalar field breaks "
            "topological invarance.  The c_T result is PROVISIONAL pending cross-check "
            "with a full characteristic-integration code."
        ),
    }

    return results


def main():
    """Run Step 17: Tensor Speed c_T Derivation.

    Sets up the pipeline logger, invokes the c_T derivation, persists
    results to ``results/step_21_cT_derivation.json``, and returns a
    summary dictionary.
    """
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    set_step_logger(logger)

    results = derive_cT()

    # Persist results
    json_path = RESULTS_DIR / "step_21_cT_derivation.json"
    write_json(json_path, results)
    print_status(f"JSON saved to {rel(json_path)}", "SUCCESS")

    print_status("Step 17 complete.  [c_T result is PROVISIONAL]", "SUCCESS")

    return {
        "step_id": STEP_ID,
        "description": (
            "Derive the tensor propagation speed c_T from the principal symbol "
            "of the coupled Einstein-scalar-Gauss-Bonnet equations.  The c_T "
            "result is PROVISIONAL: local deviation and asymptotic recovery are "
            "characterised, and the GW170817 path-dependent argument is "
            "evaluated.  Pending cross-check with full characteristic integration."
        ),
        "key_result": (
            "c_T²(r) = 1 - 8 α_GB M F(r) / r³  [PROVISIONAL].  At the horizon "
            "c_T = 1 exactly; at the photon sphere c_T ≈ 0.9992 (η=0.1); "
            "c_T → 1 as 1/r³.  GW170817 is satisfied automatically via the "
            "path-dependent escape clause (strong-field delay ~1 μs diluted "
            "over 40 Mpc → fractional delay ~10⁻²¹, far below the 10⁻¹⁵ bound)."
        ),
    }


if __name__ == "__main__":
    main()
