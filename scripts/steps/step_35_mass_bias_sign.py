#!/usr/bin/env python3
"""Derive the mass-bias sign equation for strong-field Phantom Mass.

THE CENTRAL QUESTION:
    Do slow deep clocks produce an INFLATED or DEFLATED inferred mass?

NAIVE ARGUMENT (gives the wrong sign):
    If deep clocks run slow, the observed orbital period P_o is longer
    than the local period P_s (P_o = T * P_s where T > 1).
    Kepler's third law: M ~ a^3 / P^2
    If we hold the spatial scale fixed: M_app = a^3 / P_o^2 = M_local / T^2 < M_local.
    So slow clocks alone DEFLATE the mass. This is the opposite of Phantom Mass.

THE RESOLUTION:
    The spatial scale is NOT held fixed. TEP modifies the inferred distance
    and spatial scale through the same temporal-transfer mechanism. The
    complete relation involves three factors:
    1. Temporal transfer: T_P = P_o / P_s (period stretching, T_P > 1)
    2. Spatial calibration: S_a = a_GR / a_local (apparent vs local scale)
    3. Dynamical law modification: D_dyn (modified equation of motion)

    The full relation is:
        M_app / M_local = (S_a^3 * D_dyn) / T_P^2

    Positive Phantom Mass requires:
        S_a^3 * D_dyn > T_P^2

This script derives this relation analytically and identifies the conditions
under which each sign is realised.

Outputs:
  results/step_35_mass_bias_sign.json
"""

from __future__ import annotations

import json
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")


def derive_mass_bias_sign():
    """Derive the mass-bias sign equation from first principles.

    SETUP:
        A source (star, hot spot, emitter) orbits in the deep region
        near a compact object. An observer at infinity measures:
        - Angular positions: theta(t_o)
        - Spectroscopic velocities: v_o(t_o)
        - Orbital periods: P_o (in observer time)

    The observer infers a mass using GR inversion:
        1. Convert observed angles to inferred distance: a_GR = D * theta
           where D is the (GR-inferred) distance to the source
        2. Apply Kepler's third law: M_app = (4 pi^2 a_GR^3) / (G P_o^2)
        3. Or use the velocity mass estimator: M_app = v_GR^2 * a_GR / G

    TEP replaces this with:
        1. The local orbital period P_s differs from P_o by the temporal
           transfer: P_o = T_P * P_s, where T_P = dtau_o / dtau_s > 1
        2. The local spatial scale a_local differs from a_GR by the
           spatial calibration: a_GR = S_a * a_local
        3. The local equation of motion may differ from Newtonian/GR
           by a dynamical factor D_dyn

    DERIVATION:
    """
    # Symbols
    P_s, P_o, T_P = sp.symbols('P_s P_o T_P', positive=True)
    a_local, a_GR, S_a = sp.symbols('a_local a_GR S_a', positive=True)
    v_local, v_GR = sp.symbols('v_local v_GR', real=True)
    M_local, M_app = sp.symbols('M_local M_app', positive=True)
    D_dyn = sp.symbols('D_dyn', positive=True)
    G_const = sp.Symbol('G', positive=True)

    # --- Step 1: Temporal transfer ---
    # The observer measures period P_o. The source has local period P_s.
    # P_o = T_P * P_s, where T_P > 1 (deep clocks run slow)
    P_o_expr = T_P * P_s

    # --- Step 2: Spatial calibration ---
    # The observer infers orbital radius a_GR from angular position theta
    # and distance D: a_GR = D * theta
    # The local radius a_local may differ: a_GR = S_a * a_local
    a_GR_expr = S_a * a_local

    # --- Step 3: GR mass inference (Kepler) ---
    # M_app = 4 pi^2 a_GR^3 / (G P_o^2)
    M_app_kepler = 4 * sp.pi**2 * a_GR_expr**3 / (G_const * P_o_expr**2)

    # --- Step 4: Local (TEP) mass ---
    # M_local = 4 pi^2 a_local^3 / (G P_s^2) * (1/D_dyn)
    # where D_dyn accounts for modified dynamics (D_dyn = 1 for Newtonian)
    M_local_expr = 4 * sp.pi**2 * a_local**3 / (G_const * P_s**2 * D_dyn)

    # --- Step 5: Ratio ---
    # M_app / M_local = (S_a^3 / T_P^2) * D_dyn
    ratio = sp.simplify(M_app_kepler / M_local_expr)
    print("M_app / M_local =", ratio)

    # --- Step 6: The sign condition ---
    # M_app > M_local  (positive Phantom Mass)  iff  S_a^3 * D_dyn > T_P^2
    # M_app = M_local  (zero Phantom Mass)      iff  S_a^3 * D_dyn = T_P^2
    # M_app < M_local  (negative Phantom Mass)  iff  S_a^3 * D_dyn < T_P^2

    # --- Step 7: What determines S_a? ---
    # The spatial calibration S_a = a_GR / a_local has two contributions:
    # (a) Distance calibration: D_GR / D_local
    #     The observer infers distance D_GR from standard candles, parallax, etc.
    #     If the temporal field affects photon propagation (lensing, Shapiro delay),
    #     the GR-inferred distance may differ from the true distance.
    # (b) Angular scale: theta_GR / theta_local
    #     The observed angular separation maps to a physical scale through
    #     the inferred distance. If the temporal field magnifies the apparent
    #     angular size (through lensing), a_GR > a_local.
    #
    # In the strong-field regime:
    #   - Gravitational lensing by the temporal field can magnify the apparent
    #     orbit: S_a > 1 (the orbit appears larger from outside)
    #   - The temporal transfer stretches the period: T_P > 1
    #   - The competition between S_a^3 and T_P^2 determines the sign

    # --- Step 8: Velocity-based estimator ---
    # M_app = v_GR^2 * a_GR / G
    # v_GR = v_observed (spectroscopic, already includes temporal transfer)
    # v_local = v_true (local velocity)
    # v_GR = S_v * v_local, where S_v is the velocity calibration
    # M_app = S_v^2 * v_local^2 * S_a * a_local / G
    # M_local = v_local^2 * a_local / (G * D_dyn)
    # M_app / M_local = S_v^2 * S_a * D_dyn
    S_v = sp.Symbol('S_v', positive=True)
    ratio_velocity = S_v**2 * S_a * D_dyn
    print("M_app / M_local (velocity) =", ratio_velocity)

    return {
        'kepler_ratio': str(ratio),
        'velocity_ratio': str(ratio_velocity),
        'sign_condition': {
            'positive_phantom': 'S_a^3 * D_dyn > T_P^2',
            'zero_phantom': 'S_a^3 * D_dyn = T_P^2',
            'negative_phantom': 'S_a^3 * D_dyn < T_P^2',
        },
        'naive_argument': 'Slow clocks alone (S_a=1, D_dyn=1) give M_app = M_local / T_P^2 < M_local (DEFLATION)',
        'resolution': 'Positive Phantom Mass requires spatial magnification S_a^3 to dominate temporal stretching T_P^2',
        'key_insight': 'The sign is NOT assumed; it is determined by the competition between spatial calibration and temporal transfer',
    }


def analyze_spatial_calibration():
    """Analyze what determines the spatial calibration S_a.

    The spatial calibration S_a = a_GR / a_local has three potential sources:

    1. LENSING MAGNIFICATION:
       The temporal field acts as a lens. If the photon geodesics are bent
       more strongly than in GR (or differently), the apparent angular size
       of the orbit is magnified: theta_GR > theta_true.
       This gives S_a > 1 (apparent orbit is larger).

       In the TEP framework, photons propagate on tilde_g (the matter metric).
       The conformal factor A(r) > 1 in the deep region means the photon
       paths are longer in the matter frame, potentially magnifying the
       apparent orbit.

    2. DISTANCE MiscalIBRATION:
       If the GR-inferred distance D_GR differs from the true distance D_local,
       then a_GR = D_GR * theta differs from a_local = D_local * theta_true.
       S_a = (D_GR / D_local) * (theta_GR / theta_true).

    3. ORBITAL DYNAMICS:
       The orbit itself may be different in TEP vs GR. If the temporal field
       modifies the effective potential, the local orbit may have a different
       radius than the GR-inferred orbit.

    For the strong-field case:
       The key question is whether the temporal field magnifies the apparent
       orbit enough to overcome the period stretching.

       Rough estimate: if the temporal field produces a lensing magnification
       of factor mu, then S_a ~ sqrt(mu) (angular magnification).
       The condition S_a^3 > T_P^2 becomes mu^(3/2) > T_P^2, i.e., mu > T_P^(4/3).

       For T_P ~ 10 (moderate temporal transfer), we need mu > 10^(4/3) ~ 21.5.
       This is a strong lensing requirement but not impossible in the
       strong-field regime near a compact object.
    """
    T_P, S_a, D_dyn, mu = sp.symbols('T_P S_a D_dyn mu', positive=True)

    # The critical threshold
    S_a_critical = T_P**(sp.Rational(2, 3))
    mu_critical = T_P**(sp.Rational(4, 3))

    return {
        'sources': [
            '1. Lensing magnification: temporal field bends photon paths (S_a > 1)',
            '2. Distance miscalibration: GR-inferred distance differs from true (S_a ≠ 1)',
            '3. Orbital dynamics: modified effective potential changes local orbit',
        ],
        'critical_threshold': {
            'S_a_critical': f'S_a > T_P^(2/3) = {S_a_critical}',
            'mu_critical': f'For lensing: mu > T_P^(4/3) = {mu_critical}',
            'example': 'For T_P = 10: S_a > 4.64, mu > 21.5 (strong but achievable near compact object)',
        },
        'key_point': 'The spatial calibration is NOT a free parameter — it is determined by the photon propagation on tilde_g and the orbital dynamics on the temporal-well geometry',
    }


def three_outcomes():
    """Define the three possible outcomes of the sign equation."""
    return {
        'outcome_1': {
            'name': 'Positive Phantom Mass',
            'condition': 'S_a^3 * D_dyn > T_P^2',
            'interpretation': 'Spatial magnification dominates temporal stretching. The conventional mass is inflated. M_matter < M_GR-fit.',
            'implication': 'The compact-object mass contains a temporal-inference component. Strong-field Phantom Mass is real.',
        },
        'outcome_2': {
            'name': 'Zero Phantom Mass',
            'condition': 'S_a^3 * D_dyn = T_P^2',
            'interpretation': 'Spatial and temporal effects exactly cancel. The conventional mass equals the local material mass.',
            'implication': 'TEP reduces to GR for mass inference. No Phantom Mass at black-hole scales.',
        },
        'outcome_3': {
            'name': 'Negative Phantom Mass',
            'condition': 'S_a^3 * D_dyn < T_P^2',
            'interpretation': 'Temporal stretching dominates. The conventional mass is deflated. M_matter > M_GR-fit.',
            'implication': 'Slow deep clocks make the inferred mass smaller. The object has MORE material mass than GR suggests.',
        },
        'key_principle': 'The pipeline must be allowed to return any of the three outcomes. A pipeline structured to only return outcome_1 is not testing the thesis — it is confirming it.',
    }


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 70)
    print("MASS-BIAS SIGN EQUATION (GATE -1)")
    print("=" * 70)
    print()
    print("THE CENTRAL QUESTION:")
    print("  Do slow deep clocks produce an INFLATED or DEFLATED inferred mass?")
    print()
    print("NAIVE ARGUMENT (wrong sign):")
    print("  Slow clocks → longer observed period P_o = T_P * P_s (T_P > 1)")
    print("  Kepler: M ~ a^3 / P^2")
    print("  If spatial scale fixed: M_app = M_local / T_P^2 < M_local (DEFLATION)")
    print("  This is the OPPOSITE of Phantom Mass!")
    print()

    print("--- Derivation ---")
    result = derive_mass_bias_sign()
    print()
    print(f"  Kepler ratio: M_app / M_local = {result['kepler_ratio']}")
    print(f"  Velocity ratio: M_app / M_local = {result['velocity_ratio']}")
    print()
    print("  SIGN CONDITIONS:")
    for k, v in result['sign_condition'].items():
        print(f"    {k}: {v}")
    print()
    print(f"  Naive: {result['naive_argument']}")
    print(f"  Resolution: {result['resolution']}")
    print(f"  Key insight: {result['key_insight']}")
    print()

    print("--- Spatial Calibration Analysis ---")
    spatial = analyze_spatial_calibration()
    for s in spatial['sources']:
        print(f"  {s}")
    print()
    for k, v in spatial['critical_threshold'].items():
        print(f"  {k}: {v}")
    print(f"  {spatial['key_point']}")
    print()

    print("--- Three Outcomes ---")
    outcomes = three_outcomes()
    for key in ['outcome_1', 'outcome_2', 'outcome_3']:
        o = outcomes[key]
        print(f"  {o['name']}:")
        print(f"    Condition: {o['condition']}")
        print(f"    Interpretation: {o['interpretation']}")
        print(f"    Implication: {o['implication']}")
        print()
    print(f"  KEY PRINCIPLE: {outcomes['key_principle']}")
    print()

    print("=" * 70)
    print("CONCLUSION")
    print("=" * 70)
    print()
    print("The sign of strong-field Phantom Mass is NOT assumed.")
    print("It is determined by the competition between:")
    print("  - Spatial magnification S_a^3 (from photon propagation on tilde_g)")
    print("  - Temporal stretching T_P^2 (from slow deep clocks)")
    print("  - Dynamical modification D_dyn (from modified equation of motion)")
    print()
    print("Positive Phantom Mass requires S_a^3 * D_dyn > T_P^2.")
    print("This is a testable condition, not an assumption.")
    print("The pipeline must be allowed to return any of the three outcomes.")
    print()

    results = {
        'description': 'Mass-bias sign equation for strong-field Phantom Mass',
        'central_question': 'Do slow deep clocks produce inflated or deflated inferred mass?',
        'naive_argument': result['naive_argument'],
        'derivation': {
            'kepler_ratio': result['kepler_ratio'],
            'velocity_ratio': result['velocity_ratio'],
            'formula': 'M_app / M_local = (S_a^3 * D_dyn) / T_P^2',
        },
        'sign_conditions': result['sign_condition'],
        'spatial_calibration': spatial,
        'three_outcomes': outcomes,
        'conclusion': 'The sign is determined by the competition between spatial magnification and temporal stretching. It is testable, not assumed.',
    }

    with open(os.path.join(RESULTS_DIR, 'step_35_mass_bias_sign.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {os.path.join(RESULTS_DIR, 'step_35_mass_bias_sign.json')}")


if __name__ == "__main__":
    main()
