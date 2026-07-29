#!/usr/bin/env python3
"""
GW250114 Confrontation and Universal α_GB Fit (CORRECTED)
==========================================================

Fixes from v1:
1. GW250114 final mass: 68 M_sun (detector frame), not 130 M_sun
2. β normalization: β = η/12, not η/3
3. α_GB unit conversion: GM_sun/c^2 = 1.477 km
4. Observational bound: α_EdGB < 2.9 km^2 (from Perkins et al. 2021)

The key question: is the frame-split signature observable within current constraints?
"""

import numpy as np
import json

# ============================================================
# Physical constants
# ============================================================
M_sun_kg = 1.989e30
G_Newton = 6.674e-11
c = 2.998e8
M_sun_sec = 4.9255e-6  # seconds
M_sun_km = G_Newton * M_sun_kg / c**2 / 1000  # km

# ============================================================
# GW250114 Parameters (from LIGO/Virgo posterior, arXiv:2509.08099)
# ============================================================
M_final_GW250114 = 68.1   # M_sun (detector frame, redshifted)
a_final_GW250114 = 0.68   # dimensionless spin
SNR_GW250114 = 76         # network SNR (approximately)

# Kerr QNM for a=0.68, l=2, n=0 (interpolated from Berti et al. 2006)
# omega*M ≈ 0.523 - 0.081i (approximate for a=0.68)
omega_kerr = 0.523 - 0.081j

# Physical frequency
f_GW250114_Hz = omega_kerr.real / (2 * np.pi * M_final_GW250114 * M_sun_sec)
sigma_f_GW250114 = f_GW250114_Hz / SNR_GW250114 * 2  # rough estimate

# ============================================================
# Observational bound on α_GB
# ============================================================
# From Perkins et al. (2021), arXiv:2104.11189:
# sqrt(α_EdGB) < 1.7 km → α_EdGB < 2.9 km²
# Note: EdGB ≠ shift-symmetric sGB, but this is the closest available bound.
alpha_GB_bound_km2 = 2.9  # km²

def alpha_GB_from_eta(eta, M_solar):
    """Convert η to α_GB in km².
    α_GB = η * M² / 3
    M in km = M_solar * GM_sun/c²
    """
    M_km = M_solar * M_sun_km
    return eta * M_km**2 / 3

def eta_from_alpha_GB(alpha_km2, M_solar):
    """Convert α_GB (km²) to η for a given source mass.
    η = 3 * α_GB / M²
    """
    M_km = M_solar * M_sun_km
    return 3 * alpha_km2 / M_km**2

def max_eta_for_mass(M_solar):
    """Maximum η within the observational bound α_GB < 2.9 km²."""
    return eta_from_alpha_GB(alpha_GB_bound_km2, M_solar)

# ============================================================
# TEP-sGB QNM shift (correct β normalization)
# ============================================================
def qnm_shift_fraction(eta):
    """QNM frequency shift = 19.6 * β² where β = η/12.
    This is a leading-order estimate from the horizon shift.
    NOT a full spectral calculation.
    """
    beta = eta / 12
    return 19.6 * beta**2

def confront_GW250114():
    """Confront TEP-sGB predictions with GW250114 data."""
    print("=" * 70)
    print("GW250114 CONFRONTATION (CORRECTED)")
    print("=" * 70)

    print(f"\nGW250114 parameters (from arXiv:2509.08099):")
    print(f"  Final mass (detector frame): M = {M_final_GW250114} M_sun")
    print(f"  Final spin: a = {a_final_GW250114}")
    print(f"  Network SNR ≈ {SNR_GW250114}")
    print(f"  Kerr QNM frequency (a≈0.68): f = {f_GW250114_Hz:.1f} Hz")
    print(f"  Measurement uncertainty: σ_f ≈ {sigma_f_GW250114:.2f} Hz")

    print(f"\n  Observational bound: α_GB < {alpha_GB_bound_km2} km²")
    print(f"  (from Perkins et al. 2021, EdGB bound)")
    print(f"  Maximum η for GW250114 mass: {max_eta_for_mass(M_final_GW250114):.6f}")

    print(f"\n--- TEP-sGB QNM shift for different α_GB ---")
    print(f"  {'α_GB (km²)':>12}  {'η':>10}  {'QNM shift':>12}  {'Δf (Hz)':>10}  {'σ':>8}")
    print(f"  {'----------':>12}  {'----------':>10}  {'----------':>12}  {'--------':>10}  {'------':>8}")

    for alpha_km in [0.0, 0.5, 1.0, 2.0, 2.9]:
        eta = eta_from_alpha_GB(alpha_km, M_final_GW250114)
        shift_frac = qnm_shift_fraction(eta)
        delta_f = f_GW250114_Hz * shift_frac
        n_sigma = abs(delta_f) / sigma_f_GW250114 if sigma_f_GW250114 > 0 else 0

        print(f"  {alpha_km:12.2f}  {eta:10.6f}  {100*shift_frac:11.4f}%  {delta_f:+10.4f}  {n_sigma:8.4f}")

    print(f"\n  Key finding: At the observational bound (α_GB = 2.9 km²),")
    print(f"  the QNM shift for GW250114 is {100*qnm_shift_fraction(max_eta_for_mass(M_final_GW250114)):.4f}%")
    print(f"  = {f_GW250114_Hz * qnm_shift_fraction(max_eta_for_mass(M_final_GW250114)):.4f} Hz")
    print(f"  = {abs(f_GW250114_Hz * qnm_shift_fraction(max_eta_for_mass(M_final_GW250114)))/sigma_f_GW250114:.4f}σ")
    print(f"  This is FAR below detectability for GW250114.")

# ============================================================
# Universal α_GB and observability of frame-split
# ============================================================
def universal_analysis():
    """Analyze the observability of the frame-split across source masses."""
    print("\n" + "=" * 70)
    print("UNIVERSAL α_GB AND FRAME-SPLIT OBSERVABILITY")
    print("=" * 70)

    print(f"\nObservational bound: α_GB < {alpha_GB_bound_km2} km²")
    print(f"GM_sun/c² = {M_sun_km:.4f} km")
    print()

    # For each source mass, compute max η and the resulting observables
    print(f"  {'M (M_sun)':>12}  {'η_max':>10}  {'Shadow %':>10}  {'ISCO tilde %':>14}  {'QNM shift %':>12}")
    print(f"  {'----------':>12}  {'------':>10}  {'--------':>10}  {'------------':>14}  {'-----------':>12}")

    sources = [
        ("10 M_sun (low-mass stellar)", 10.0),
        ("30 M_sun (stellar)", 30.0),
        ("60 M_sun (GW150914-like)", 60.0),
        ("68 M_sun (GW250114 remnant)", 68.0),
        ("100 M_sun (intermediate)", 100.0),
        ("10^6 M_sun (Sgr A*)", 1e6),
        ("6.5×10^9 M_sun (M87*)", 6.5e9),
    ]

    results = []
    for name, M_solar in sources:
        eta_max = max_eta_for_mass(M_solar)

        # Shadow shift: O(β²) = O((η/12)²)
        # From the v3 script at η = -0.1: shadow = -0.0440%
        # The sign is negative for both mass-inflation and mass-deflation
        # because the geometric metric receives the O(η²) inward shift.
        # Scale: shadow_shift = -0.0440 * (|eta|/0.1)²
        shadow_pct = -0.0440 * (abs(eta_max) / 0.1)**2

        # ISCO on tilde g: O(η) from conformal factor
        # From the v3 script at η = -0.1 (mass-inflation branch): ISCO = +1.95%
        # This scales as |η|, so: ISCO = +1.95 * (|eta_max|/0.1)
        isco_pct = +1.95 * (abs(eta_max) / 0.1)

        # QNM shift: O(β²) = 19.6 * (η/12)²
        qnm_pct = 100 * qnm_shift_fraction(eta_max)

        print(f"  {M_solar:12.2e}  {eta_max:10.6f}  {shadow_pct:+10.4f}  {isco_pct:+14.4f}  {qnm_pct:+12.4f}")
        results.append({
            "name": name, "M": M_solar, "eta_max": eta_max,
            "shadow_pct": shadow_pct, "isco_pct": isco_pct, "qnm_pct": qnm_pct
        })

    print()
    print("  KEY FINDINGS:")
    print()
    print("  1. The frame-split is REAL but SMALL within current constraints:")
    print(f"     - Best case (10 M_sun): shadow {results[0]['shadow_pct']:+.4f}%, ISCO {results[0]['isco_pct']:+.4f}%")
    print(f"     - GW250114 (68 M_sun): shadow {results[3]['shadow_pct']:+.6f}%, ISCO {results[3]['isco_pct']:+.6f}%")
    print(f"     - Sgr A* (10^6 M_sun): shadow {results[5]['shadow_pct']:+.2e}%, ISCO {results[5]['isco_pct']:+.2e}%")
    print()
    print("  2. The ISCO shift is O(η) and dominates over the O(η²) shadow shift.")
    print("     The frame-split (opposite signs, different orders) is a divergence:")
    print("     shadow contracts while matter-metric ISCO expands.")
    print()
    print("  3. The strongest signal is the ISCO on tilde_g for low-mass BHs:")
    print(f"     - 10 M_sun: ISCO shift = {results[0]['isco_pct']:.4f}% (potentially observable)")
    print(f"     - 30 M_sun: ISCO shift = {results[1]['isco_pct']:.6f}%")
    print(f"     - 60 M_sun: ISCO shift = {results[2]['isco_pct']:.6f}%")
    print()
    print("  4. The QNM shift is negligible for all sources within current bounds.")
    print()
    print("  5. The frame-split is NOT a detection by sign divergence alone.")
    print("     A quantitative consistency relation between shadow, ISCO, and QNM")
    print("     — derived from one frozen action — is needed as the discriminator.")
    print()
    print("  CONCLUSION: TEP-sGB is observationally viable but the frame-split")
    print("  signature is extremely small within current α_GB constraints.")
    print("  The strongest near-term probe is the ISCO on tilde_g for low-mass")
    print("  stellar black holes (M ~ 10 M_sun), where the O(η) conformal")
    print("  factor produces a ~0.75% ISCO shift at the observational bound.")
    print("  Next-generation detectors (ET, CE) and precision spectroscopy")
    print("  of low-mass mergers offer the best path forward.")

    return results

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    confront_GW250114()
    results = universal_analysis()

    output = {
        "GW250114": {
            "M_final_detector": M_final_GW250114,
            "spin": a_final_GW250114,
            "SNR": SNR_GW250114,
            "f_Kerr_Hz": f_GW250114_Hz,
            "sigma_f_Hz": sigma_f_GW250114,
            "source": "arXiv:2509.08099",
        },
        "observational_bound": {
            "alpha_GB_max_km2": alpha_GB_bound_km2,
            "source": "Perkins et al. 2021, arXiv:2104.11189 (EdGB)",
            "note": "EdGB bound applied as proxy for shift-symmetric sGB",
        },
        "unit_conversion": {
            "GM_sun_c2_km": M_sun_km,
            "alpha_GB_formula": "alpha_GB = eta * M^2 / 3, M in km",
        },
        "universal_analysis": [
            {"name": r["name"], "M_solar": r["M"], "eta_max": r["eta_max"],
             "shadow_pct": r["shadow_pct"], "isco_pct": r["isco_pct"],
             "qnm_pct": r["qnm_pct"]}
            for r in results
        ],
    }

    with open("results/step_41_gw250114_confrontation.json", "w") as f:
        json.dump(output, f, indent=2)

    print("\nResults saved to results/step_41_gw250114_confrontation.json")
