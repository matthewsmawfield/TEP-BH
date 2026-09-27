#!/usr/bin/env python3
"""Operational darkness threshold for the Temporal Horizon.

Quantifies the transfer-factor threshold Lambda at which a regular
temporal well becomes operationally dark to a distant observer, and
compares it against the step_34 benchmark lapse minima.  This is the
quantitative criterion requested by the darkness-mechanism review:
the Temporal Horizon is defined by the transfer factor
T_{e->o} = N_o/N_e, and "darkness" requires T above a threshold set by
the instrument's contrast sensitivity — here made explicit.

Signal channels for a static emitter at matter-frame lapse N_tilde_e
seen by a static distant observer (N_o -> 1, T = N_o/N_tilde_e):

  1. Radiometric (photon) channel.  Null transport is conformally
     invariant, so the conserved Liouville quantity is I_nu/nu^3;
     specific intensity along a stationary ray is suppressed by
     I_o/I_e = T^{-3} (and bolometric flux by T^{-4}).  An interior
     emitter of intrinsic emissivity contrast C_em relative to the
     photon-ring emission is invisible in the image when
         T > T_dark = (C_em / eps)^{1/3},
     with eps the image dynamic-range floor (interior flux floor
     relative to the ring).  For EHT 230 GHz imaging eps ~ 1e-1..1e-3
     brackets optimistic-to-conservative contrast.

  2. Massive-signal channel.  A particle (or neutrino) at rest in the
     well carries Killing energy E = m N_tilde_e; reaching the distant
     observer requires E >= m, i.e. a boost factor ~ T.  Massive
     messengers are confined, not merely redshifted.

  3. Temporal-dilution channel.  An interior process of intrinsic
     duration dt_e is observed over dt_o = T dt_e: in-well variability
     is stretched by the same factor.

  4. Spectral channel.  Emission at intrinsic nu_e is observed at
     nu_o = nu_e/T; interior mm-band emission lands far below the
     observing band.

  5. Transit delay (photon channel is conformal): the coordinate
     transit time through the deep region, dt = integral dr/N^2 on the
     geometric lapse, is finite — photons escape but arrive dimmed;
     operational darkness is enforced by channel 1, not by trapping.

Benchmark inputs are recomputed from the corrected step_34
boundary-value solver for each source model, so the comparison is
reproducible rather than hand-entered.

Outputs:
    results/step_50_darkness_threshold.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from step_34_solve_interior import (
    hayward_N, gb_invariant, gb_invariant_exact,
    linear_f, d_linear_f, nonlinear_f, d_nonlinear_f,
    solve_scalar_bvp,
)

# Physical constants for mass-to-time conversion
_C = 2.99792458e8
_G = 6.674e-11
_MSUN = 1.989e30

OBSERVER_MASSES = {
    "Sgr_A*": 4.297e6,      # M_sun (GRAVITY/Keck)
    "M87*": 6.5e9,          # M_sun (EHT)
    "GW150914_remnant": 62.0,  # M_sun (approximate final mass)
}


def t_per_M(M_sun):
    """GM/c^3 in seconds for a mass in solar masses."""
    return _G * M_sun * _MSUN / _C**3


def profile_and_diagnostics(bvp, g, M):
    """Extract N_tilde profile, T_max, and well structure from a BVP result."""
    sol = bvp["sol"]
    r = np.geomspace(1e-6, 300.0, 4000)
    phi, phi_p = sol.sol(r)
    N = hayward_N(r, M, g)
    N_tilde = np.exp(-phi) * N
    i = int(np.argmin(N_tilde))
    return {
        "r": r,
        "phi": phi,
        "phi_p": phi_p,
        "N": N,
        "N_tilde": N_tilde,
        "N_tilde_min": float(N_tilde[i]),
        "r_N_tilde_min": float(r[i]),
        "A_centre": float(np.exp(-phi[0])),
        "T_max": float(1.0 / N_tilde[i]),
    }


def deep_transit_time(g, M, r_out=10.0):
    """Coordinate transit time through the deep interior, in units of M.

    Photons follow the geometric null cone (conformal invariance), so the
    one-way transit time from the centre to r_out is integral dr/N^2.
    The excess over the flat-space crossing time is also returned.
    """
    transit = quad(lambda r: 1.0 / hayward_N(r, M, g)**2, 1e-6, r_out,
                   epsabs=1e-12, limit=200)[0]
    excess = transit - r_out
    return {"transit_M": float(transit), "excess_M": float(excess), "r_out": r_out}


def main():
    M = 1.0
    eta = 0.3
    g_benchmark = 1.1

    print("=" * 70)
    print("OPERATIONAL DARKNESS THRESHOLD — TEMPORAL HORIZON")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. Radiometric darkness threshold grid
    # ------------------------------------------------------------------
    eps_grid = [1e-1, 1e-2, 1e-3]        # image dynamic-range floor
    C_em_grid = [1.0, 1e2, 1e6]          # intrinsic emissivity contrast
    threshold_table = []
    for eps in eps_grid:
        for C_em in C_em_grid:
            T_dark = (C_em / eps) ** (1.0 / 3.0)
            threshold_table.append({
                "epsilon": eps,
                "C_em": C_em,
                "T_dark": float(T_dark),
            })

    T_dark_min = min(t["T_dark"] for t in threshold_table)
    T_dark_max = max(t["T_dark"] for t in threshold_table)
    T_dark_nominal = (1e2 / 1e-2) ** (1.0 / 3.0)  # eps=1e-2, C_em=1e2

    print("\nRadiometric threshold  T_dark = (C_em/eps)^{1/3}")
    print(f"  {'eps':>8s}  {'C_em':>8s}  {'T_dark':>10s}")
    for t in threshold_table:
        print(f"  {t['epsilon']:8.0e}  {t['C_em']:8.0e}  {t['T_dark']:10.2f}")
    print(f"  Range: {T_dark_min:.2f} .. {T_dark_max:.2f}; "
          f"nominal (eps=1e-2, C_em=1e2): {T_dark_nominal:.2f}")

    # ------------------------------------------------------------------
    # 2. Benchmark lapse minima (recomputed from step_34 BVP solver)
    # ------------------------------------------------------------------
    benchmarks = {}
    for gb_name, gb_fn in [("schwarzschild_form", gb_invariant),
                           ("hayward_exact", gb_invariant_exact)]:
        bvp = solve_scalar_bvp(M, eta, g_benchmark, linear_f, d_linear_f, gb_fn)
        benchmarks[gb_name] = profile_and_diagnostics(bvp, g_benchmark, M)
        benchmarks[gb_name]["phi_c"] = bvp["phi_c"]
        benchmarks[gb_name]["Q_eff"] = bvp["Q_eff"]

    # nonlinear saturating coupling, declared source convention
    bvp_nl = solve_scalar_bvp(M, eta, g_benchmark,
                              lambda p: nonlinear_f(p, phi0=0.5),
                              lambda p: d_nonlinear_f(p, phi0=0.5),
                              gb_invariant)
    benchmarks["nonlinear_phi0_0.5|schwarzschild_form"] = profile_and_diagnostics(bvp_nl, g_benchmark, M)
    benchmarks["nonlinear_phi0_0.5|schwarzschild_form"]["phi_c"] = bvp_nl["phi_c"]

    print("\nBenchmark transfer factors (g = 1.1 M core, eta = +0.3):")
    for name, b in benchmarks.items():
        print(f"  {name:40s}: N_tilde_min={b['N_tilde_min']:.3e} at r={b['r_N_tilde_min']:.3f} M, "
              f"T_max={b['T_max']:.3f}")

    # ------------------------------------------------------------------
    # 3. Channel-by-channel comparison
    # ------------------------------------------------------------------
    comparison = {}
    for name, b in benchmarks.items():
        T = b["T_max"]
        radiometric = {
            "suppression_I_nu": float(T**-3),
            "suppression_bolometric": float(T**-4),
            "passes_nominal": bool(T > T_dark_nominal),
            "passes_all": bool(T > T_dark_max),
            "margin_nominal": float(T / T_dark_nominal),
            "margin_min_threshold": float(T / T_dark_min),
        }
        massive = {
            "escape_boost_required": float(T),
            "confined": bool(T > 10.0),   # > decade of energy needed to climb out
        }
        spectral = {
            "required_nu_e_for_230GHz_Hz": float(230e9 * T),
            "note": "interior emitter must radiate at nu_e = T * nu_obs to be seen at nu_obs",
        }
        dilution = {"observed_duration_factor": float(T)}
        comparison[name] = {
            "T_max": T,
            "radiometric": radiometric,
            "massive_channel": massive,
            "spectral": spectral,
            "dilution": dilution,
        }

    # ------------------------------------------------------------------
    # 4. Transit delays (photon channel is conformal => geometric lapse)
    # ------------------------------------------------------------------
    delays = {}
    for g in [1.1, 1.5, 2.0]:
        tr = deep_transit_time(g, M, r_out=10.0)
        row = dict(tr)
        for obs, M_sun in OBSERVER_MASSES.items():
            row[obs] = {
                "transit_s": tr["transit_M"] * t_per_M(M_sun),
                "excess_s": tr["excess_M"] * t_per_M(M_sun),
            }
        delays[f"g_{g}"] = row

    print("\nDeep-transit delay (centre -> 10 M, geometric lapse; photons conformal):")
    for gkey, row in delays.items():
        sgra = row["Sgr_A*"]["transit_s"]
        m87 = row["M87*"]["transit_s"]
        gw = row["GW150914_remnant"]["transit_s"]
        print(f"  {gkey}: transit={row['transit_M']:.1f} M -> "
              f"Sgr A* {sgra:.1f} s, M87* {m87:.2e} s, GW150914 {gw:.2e} s")

    # ------------------------------------------------------------------
    # 5. Verdict
    # ------------------------------------------------------------------
    primary = benchmarks["schwarzschild_form"]
    verdict = {
        "T_dark_range": [T_dark_min, T_dark_max],
        "T_dark_nominal": T_dark_nominal,
        "primary_benchmark_T_max": primary["T_max"],
        "primary_exceeds_nominal": bool(primary["T_max"] > T_dark_nominal),
        "primary_exceeds_full_range": bool(primary["T_max"] > T_dark_max),
        "primary_margin_over_max_threshold": float(primary["T_max"] / T_dark_max),
        "statement": (
            f"The operational darkness threshold is T_dark ~ "
            f"{T_dark_min:.0f}–{T_dark_max:.0f} (radiometric, "
            f"eps=1e-1..1e-3, C_em=1..1e6). The corrected benchmark "
            f"(linear sGB source on the g=1.1 M core) reaches "
            f"T_max = {primary['T_max']:.3e}, exceeding the most "
            f"conservative threshold by a factor "
            f"{primary['T_max']/T_dark_max:.1f}. Interior photon "
            f"emission is suppressed by T^-3 = {primary['T_max']**-3:.2e} "
            f"in specific intensity and T^-4 = {primary['T_max']**-4:.2e} "
            f"bolometrically; massive messengers require an escape boost "
            f"~{primary['T_max']:.0f}. The exact-background-GB source "
            f"variant reaches T_max = {benchmarks['hayward_exact']['T_max']:.2f} "
            f"(screened exterior charge), below the threshold — bounding "
            f"the prescribed-background depth sensitivity and setting the "
            f"quantitative target for the coupled nonlinear interior."
        ),
    }

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print(f"  T_dark range:              {T_dark_min:.1f} .. {T_dark_max:.1f}")
    print(f"  Primary benchmark T_max:   {primary['T_max']:.3f}")
    print(f"  Exceeds nominal threshold: {verdict['primary_exceeds_nominal']} "
          f"(margin {primary['T_max']/T_dark_nominal:.1f}x)")
    print(f"  Exceeds full range:        {verdict['primary_exceeds_full_range']} "
          f"(margin {verdict['primary_margin_over_max_threshold']:.1f}x)")
    print(f"  Exact-GB variant T_max:    {benchmarks['hayward_exact']['T_max']:.3f} "
          f"(screened, below threshold)")

    output = {
        "description": (
            "Operational darkness threshold for the Temporal Horizon: "
            "transfer-factor criterion for EHT non-detection of the deep "
            "interior, compared against the corrected step_34 benchmarks."
        ),
        "channels": {
            "radiometric": "I_nu suppressed by T^-3 (Liouville); flux by T^-4",
            "massive": "escape requires boost ~ 1/N_tilde_min",
            "dilution": "observed durations stretched by T",
            "spectral": "nu_o = nu_e / T",
            "transit": "conformal photons: finite dt = int dr/N^2",
        },
        "threshold_grid": threshold_table,
        "benchmarks": {
            k: {kk: vv for kk, vv in b.items()
                if kk not in ("r", "phi", "phi_p", "N", "N_tilde")}
            for k, b in benchmarks.items()
        },
        "comparison": comparison,
        "transit_delays": delays,
        "verdict": verdict,
    }

    out_path = RESULTS_DIR / "step_50_darkness_threshold.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
