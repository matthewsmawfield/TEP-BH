#!/usr/bin/env python3
"""Step 02: Null expansions and geodesic structure.

Computes the null expansions theta_+, theta_- for the TEP matter metric
on both Schwarzschild and Hayward backgrounds, verifying:

  - Temporal freeze: theta_+ -> 0 as r -> 0 (outgoing null geodesics
    have vanishing expansion, infinite affine parameter)
  - Regular centre: theta_- remains finite as r -> 0 (ingoing null
    geodesics reach the centre in finite affine parameter)

For a static spherically symmetric metric ds^2 = -F dt^2 + G dr^2 + H dOmega^2,
the null expansions for the two radial null congruences are:

  theta_+ = (1/sqrt(g_rr)) * (d/dr)(ln(sqrt(g_thth))) * 2
          = (1/sqrt(G)) * (H'/(H)) * ... 

More precisely, for a null congruence with tangent k^mu = (dt/dlambda, dr/dlambda),
the expansion is theta = nabla_mu k^mu. For the outgoing radial null geodesic
in a static spherical metric:

  ds^2 = 0 => dt/dr = +/- sqrt(G/F)

  k^mu_out = (sqrt(G/F), 1, 0, 0) / sqrt(G)  (affine-normalized)
  theta_+ = (1/sqrt(G)) * d/dr)(2*sqrt(H)) / sqrt(H) = 2/(sqrt(G)) * H'/(2H)
          = H'/(H * sqrt(G))

  theta_- = -H'/(H * sqrt(G))   (ingoing, opposite sign)

For the TEP matter metric g_tilde = A^2 g:
  F_tilde = A^2 F, G_tilde = A^2 G, H_tilde = A^2 r^2

  theta_+_tilde = H_tilde' / (H_tilde * sqrt(G_tilde))
                = (2A A' r^2 + 2A^2 r) / (A^2 r^2 * A sqrt(G))
                = (2 A'/A + 2/r) / (A sqrt(G))
                = 2(A' r + A) / (A^2 r^2 * A sqrt(G)) * r^2
                = 2(A' r + A) / (A^3 r * sqrt(G))

With A = (r_h/r)^phi0, A' = -phi0 * (r_h/r)^phi0 / r = -phi0 * A / r:
  A' r + A = -phi0 * A + A = A(1 - phi0)
  theta_+ = 2 * A(1-phi0) / (A^3 r sqrt(G)) = 2(1-phi0) / (A^2 r sqrt(G))

For phi0 = 1: theta_+ = 0 exactly! (The conformal factor exactly cancels
the geometric expansion.) This is the temporal freeze.

For phi0 < 1: theta_+ > 0 (expanding, normal)
For phi0 > 1: theta_+ < 0 (contracting)

The affine parameter integral:
  lambda = int sqrt(-det(g_2D)) dr / (g_vv * something) ...
  For the EF form: lambda = int sqrt(-det gtilde_2D) dr

Outputs:
  results/step_17_null_expansions.json
  results/step_17_null_expansions.csv
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

from scripts.utils.logger import TEPLogger  # noqa: E402

r, M, ell, phi0, r_h = sp.symbols("r M ell phi0 r_h", positive=True)
RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")


def null_expansions(F, G, H):
    """Compute null expansions theta_+, theta_- for ds^2 = -F dt^2 + G dr^2 + H dOmega^2.

    theta_+ = H' / (H * sqrt(G))    (outgoing)
    theta_- = -H' / (H * sqrt(G))   (ingoing)
    """
    Hp = sp.diff(H, r)
    theta_plus = Hp / (H * sp.sqrt(G))
    theta_minus = -Hp / (H * sp.sqrt(G))
    return sp.simplify(theta_plus), sp.simplify(theta_minus)


def affine_parameter_integrand(F, G, H):
    """The affine parameter integrand for radial null geodesics.

    For a radial null geodesic, the affine parameter is:
      lambda = int sqrt(G/F) dr   (outgoing, in standard coords)
      lambda = int sqrt(-det(g_2D)) / |g_vv| dr   (in EF coords)

    In standard coordinates with ds^2 = -F dt^2 + G dr^2 + H dOmega^2:
      Null condition: -F (dt/dlambda)^2 + G (dr/dlambda)^2 = 0
      => dt/dlambda = +/- sqrt(G/F) * dr/dlambda
      Affine parameter: lambda = int sqrt(G) dr / |k^r| ... 

    For the proper affine parameter of a radial null geodesic:
      lambda = int sqrt(G/F) * dr  (this is the coordinate time integral)
      But the actual affine parameter is:
      lambda = int sqrt(-g_{2D}) dr / (energy * g_{vv}) in EF coords.

    For standard coords, the affine parameter for outgoing null:
      lambda_out = int sqrt(G) dr  (with appropriate normalization)
      lambda_in  = int sqrt(G) dr  (same magnitude, opposite direction)

    The key quantity for geodesic completeness is whether this integral
    converges or diverges as r -> 0.
    """
    # The affine parameter integrand is sqrt(G) for radial null geodesics
    # in standard coordinates (with appropriate normalization).
    # For the conformal metric: sqrt(G_tilde) = A * sqrt(G)
    return sp.sqrt(G)


def compute_null_structure(name, F, G, H, A_factor=None, compute_limits=True):
    """Compute null expansions, affine parameter, and limits."""
    print(f"\n=== {name} ===")
    theta_p, theta_m = null_expansions(F, G, H)
    lambda_integrand = affine_parameter_integrand(F, G, H)

    limits = {}
    if compute_limits:
        for label, expr in [("theta_plus", theta_p), ("theta_minus", theta_m),
                            ("lambda_integrand", lambda_integrand)]:
            try:
                lim = sp.limit(expr, r, 0)
                limits[label + "_r0"] = str(lim)
            except (sp.PoleError, TypeError, ValueError):
                limits[label + "_r0"] = "undefined"

    print(f"  theta_+ = {sp.simplify(theta_p)}")
    print(f"  theta_- = {sp.simplify(theta_m)}")
    print(f"  lambda integrand = {sp.simplify(lambda_integrand)}")
    if limits:
        print(f"  limits: {limits}")

    return {
        "name": name,
        "theta_plus": str(sp.simplify(theta_p)),
        "theta_minus": str(sp.simplify(theta_m)),
        "lambda_integrand": str(sp.simplify(lambda_integrand)),
        "limits": limits,
        "A": str(A_factor) if A_factor else "1",
    }


# ---------------------------------------------------------------------------
# Metric definitions
# ---------------------------------------------------------------------------
def schwarzschild_F():
    return 1 - 2 * M / r


def schwarzschild_G():
    return 1 / (1 - 2 * M / r)


def hayward_F():
    return 1 - 2 * M * r**2 / (r**3 + 2 * M * ell**2)


def hayward_G():
    return 1 / hayward_F()


def conformal_A(phi0_val):
    return (r_h / r) ** phi0_val


# ---------------------------------------------------------------------------
# Regression tests
# ---------------------------------------------------------------------------
def regression_tests():
    """Verify known null expansion results."""
    print("\n" + "=" * 60)
    print("REGRESSION TESTS")
    print("=" * 60)
    results = {}

    # 1. Schwarzschild: theta_+ = 2*sqrt(F)/r (outside horizon, F>0)
    F_s = schwarzschild_F()
    G_s = schwarzschild_G()
    H_s = r**2
    tp_s, tm_s = null_expansions(F_s, G_s, H_s)
    tp_s_simplified = sp.simplify(tp_s)
    # theta_+ = H'/(H*sqrt(G)) = 2r/(r^2 * sqrt(1/F)) = 2*sqrt(F)/r
    expected_tp_s = 2 * sp.sqrt(F_s) / r
    # Compare squared (to avoid sign issues inside horizon)
    diff = sp.simplify(tp_s_simplified**2 - expected_tp_s**2)
    results["theta_plus_Schwarzschild"] = {
        "expected": "2*sqrt(F)/r",
        "computed": str(tp_s_simplified),
        "pass": diff == 0,
    }
    print(f"  theta+_Schwarzschild = {tp_s_simplified}  PASS={diff == 0}")

    # 2. TEP-Hayward phi0=1: theta_+ = 0 exactly (temporal freeze)
    F_h = hayward_F()
    G_h = hayward_G()
    A1 = conformal_A(sp.Integer(1))
    F_tilde = A1**2 * F_h
    G_tilde = A1**2 * G_h
    H_tilde = A1**2 * r**2
    tp_th, tm_th = null_expansions(F_tilde, G_tilde, H_tilde)
    tp_th_simplified = sp.simplify(tp_th)
    results["theta_plus_TEP_Hayward_phi0_1"] = {
        "expected": "0 (temporal freeze)",
        "computed": str(tp_th_simplified),
        "pass": tp_th_simplified == 0,
    }
    print(f"  theta+_TEP_Hayward(phi0=1) = {tp_th_simplified}  PASS={tp_th_simplified == 0}")

    # 3. TEP-Schwarzschild phi0=1: theta_+ = 0 exactly (temporal freeze)
    F_tilde_s = A1**2 * F_s
    G_tilde_s = A1**2 * G_s
    H_tilde_s = A1**2 * r**2
    tp_ts, tm_ts = null_expansions(F_tilde_s, G_tilde_s, H_tilde_s)
    tp_ts_simplified = sp.simplify(tp_ts)
    results["theta_plus_TEP_Schwarzschild_phi0_1"] = {
        "expected": "0 (temporal freeze)",
        "computed": str(tp_ts_simplified),
        "pass": tp_ts_simplified == 0,
    }
    print(f"  theta+_TEP_Schwarzschild(phi0=1) = {tp_ts_simplified}  PASS={tp_ts_simplified == 0}")

    # 4. TEP-Schwarzschild phi0=2: theta_+ != 0 (no temporal freeze)
    A2 = conformal_A(sp.Integer(2))
    F_tilde_s2 = A2**2 * F_s
    G_tilde_s2 = A2**2 * G_s
    H_tilde_s2 = A2**2 * r**2
    tp_ts2, tm_ts2 = null_expansions(F_tilde_s2, G_tilde_s2, H_tilde_s2)
    tp_ts2_simplified = sp.simplify(tp_ts2)
    results["theta_plus_TEP_Schwarzschild_phi0_2"] = {
        "expected": "nonzero (no freeze)",
        "computed": str(tp_ts2_simplified),
        "pass": tp_ts2_simplified != 0,
    }
    print(f"  theta+_TEP_Schwarzschild(phi0=2) = {tp_ts2_simplified}  PASS={tp_ts2_simplified != 0}")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    logger = TEPLogger("step_17_null_expansions")
    logger.info("Step 02: Null Expansions and Geodesic Structure")

    # --- Regression tests ---
    test_results = regression_tests()
    all_pass = all(t["pass"] for t in test_results.values())
    if not all_pass:
        print("\n*** REGRESSION TESTS FAILED ***")
        for name, t in test_results.items():
            if not t["pass"]:
                print(f"    FAIL: {name}: expected {t.get('expected')}, got {t.get('computed')}")

    # --- Per-metric null structure ---
    metrics = {}

    # Schwarzschild baseline
    F_s = schwarzschild_F()
    G_s = schwarzschild_G()
    metrics["schwarzschild"] = compute_null_structure(
        "Schwarzschild", F_s, G_s, r**2
    )

    # Hayward baseline
    F_h = hayward_F()
    G_h = hayward_G()
    metrics["hayward"] = compute_null_structure(
        "Hayward", F_h, G_h, r**2
    )

    # TEP-Schwarzschild phi0=1
    A1 = conformal_A(sp.Integer(1))
    metrics["tep_schwarzschild_phi0_1"] = compute_null_structure(
        "TEP-Schwarzschild (phi0=1)",
        A1**2 * F_s, A1**2 * G_s, A1**2 * r**2, A1
    )

    # TEP-Schwarzschild phi0=2
    A2 = conformal_A(sp.Integer(2))
    metrics["tep_schwarzschild_phi0_2"] = compute_null_structure(
        "TEP-Schwarzschild (phi0=2)",
        A2**2 * F_s, A2**2 * G_s, A2**2 * r**2, A2
    )

    # TEP-Hayward phi0=1 (regular benchmark) — skip limits (too complex for SymPy)
    metrics["tep_hayward_phi0_1"] = compute_null_structure(
        "TEP-Hayward (phi0=1)",
        A1**2 * F_h, A1**2 * G_h, A1**2 * r**2, A1, compute_limits=False
    )

    # --- Affine parameter analysis ---
    # For the TEP metrics, compute the affine parameter integral
    # lambda = int sqrt(G_tilde) dr from r_min to r_h
    # With G_tilde = A^2 * G = (r_h/r)^{2*phi0} * G
    # sqrt(G_tilde) = (r_h/r)^{phi0} * sqrt(G)
    # For phi0=1: sqrt(G_tilde) = (r_h/r) * sqrt(G)
    #   In deep interior (Schwarzschild): sqrt(G) ~ sqrt(r/(2M)), so
    #   sqrt(G_tilde) ~ (r_h/r) * sqrt(r/(2M)) = r_h / sqrt(2M r) ~ r^{-1/2}
    #   => lambda ~ int r^{-1/2} dr ~ r^{1/2} -> 0 (CONVERGES, finite affine parameter)
    # For phi0=2: sqrt(G_tilde) = (r_h/r)^2 * sqrt(G)
    #   sqrt(G_tilde) ~ (r_h/r)^2 * sqrt(r/(2M)) ~ r^{-3/2}
    #   => lambda ~ int r^{-3/2} dr ~ r^{-1/2} -> oo (DIVERGES, infinite affine parameter)

    print("\n" + "=" * 60)
    print("AFFINE PARAMETER ANALYSIS")
    print("=" * 60)

    # Deep interior: F ~ -2M/r, G ~ -r/(2M), sqrt(G) ~ sqrt(r/(2M))
    # sqrt(G_tilde) = A * sqrt(G) = (r_h/r)^phi0 * sqrt(r/(2M))
    # = r_h^phi0 * r^{-phi0} * r^{1/2} / sqrt(2M) = r_h^phi0 / sqrt(2M) * r^{1/2 - phi0}
    # lambda = int r^{1/2 - phi0} dr
    #   phi0 = 1: exponent = -1/2, lambda ~ r^{1/2} -> 0 (CONVERGES)
    #   phi0 = 2: exponent = -3/2, lambda ~ r^{-1/2} -> oo (DIVERGES)
    #   phi0 = 3/2: exponent = -1, lambda ~ ln(r) -> -oo (DIVERGES logarithmically)

    affine_analysis = {}
    for phi0_val in [sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2), sp.Integer(2)]:
        exponent = sp.Rational(1, 2) - phi0_val
        if exponent > -1:
            behavior = "converges (finite affine parameter)"
            limit = "finite"
        elif exponent == -1:
            behavior = "diverges logarithmically (infinite affine parameter)"
            limit = "ln(r) -> -oo"
        else:
            behavior = "diverges as power law (infinite affine parameter)"
            limit = f"r^{exponent + 1} -> oo"

        key = f"phi0_{phi0_val}"
        affine_analysis[key] = {
            "phi0": str(phi0_val),
            "exponent": str(exponent),
            "behavior": behavior,
            "limit": limit,
            "inward_finite": exponent > -1,
            "outward_infinite": exponent <= -1,  # same integral, both directions
        }
        print(f"  phi0={phi0_val}: sqrt(G_tilde) ~ r^{exponent}, lambda ~ r^{exponent+1}  -> {behavior}")

    # KEY FINDING: for phi0=1, BOTH inward and outward affine parameters
    # converge (finite). The "temporal freeze" (theta_+ = 0) means the
    # expansion vanishes, but the affine parameter is still finite.
    # The manuscript claims outward is infinite and inward is finite.
    # This needs careful analysis: theta_+ = 0 means the congruence has
    # zero expansion, but the affine parameter depends on the integral
    # of sqrt(G_tilde), which for phi0=1 converges from BOTH directions.
    # The asymmetry must come from the EF coordinate structure or the
    # disformal term, not the pure conformal metric.

    print("\n  KEY FINDING: For phi0=1, the pure conformal affine parameter")
    print("  converges from BOTH directions (inward and outward).")
    print("  The inward/outward asymmetry claimed in the manuscript requires")
    print("  the EF coordinate structure or the disformal term, not just the")
    print("  conformal factor. This is Finding #9 from the review.")

    # --- Write outputs ---
    output = {
        "step": "02_null_expansions",
        "description": "Null expansions and geodesic structure",
        "regression_tests": test_results,
        "all_tests_pass": all_pass,
        "metrics": metrics,
        "affine_parameter_analysis": affine_analysis,
        "key_finding": (
            "For phi0=1, theta_+ = 0 exactly (temporal freeze confirmed). "
            "However, the pure conformal affine parameter converges from BOTH "
            "directions (inward and outward). The inward-finite/outward-infinite "
            "asymmetry requires EF coordinates or the disformal term."
        ),
    }

    out_json = os.path.join(RESULTS_DIR, "step_17_null_expansions.json")
    with open(out_json, "w") as fh:
        json.dump(output, fh, indent=2, default=str)
    print(f"\nWrote {out_json}")

    logger.info("Step 02 complete")
    return output


if __name__ == "__main__":
    main()
