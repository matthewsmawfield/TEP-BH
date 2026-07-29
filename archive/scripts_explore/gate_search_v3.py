#!/usr/bin/env python3
"""Refined approach: combine mild conformal growth with disformal temporal effect.

The tension:
- A -> infinity: suppresses curvature but makes areal radius diverge (Gate 1 fail)
- A = 1: keeps areal radius finite but doesn't suppress curvature (Gate 5 fail)
- B large: creates temporal effect but breaks determinant (Gate 2 fail)

New idea: use A that grows MILDLY (e.g., A ~ r^{-alpha} with alpha < 1)
so that:
  - areal radius rho = A*r ~ r^{1-alpha} -> 0 (finite, doesn't diverge)
  - curvature K ~ A^{-12} * K_Schw ~ r^{12*alpha - 6} -> 0 if alpha > 0.5
  - determinant: A^2 dominates (as in conformal prototype)

So alpha in (0.5, 1.0) gives:
  - Finite areal radius (goes to 0, not infinity)
  - Vanishing curvature
  - Lorentzian determinant

The disformal term with q != 0 then adds the temporal effect on top.

But wait — if areal radius -> 0, that means the spatial size shrinks.
That's the opposite of "opening up" but it's also not "continuous regular space."
The user wants space to remain ordinary, not shrink to a point.

Actually, re-reading the user's requirements:
  "The physical areal radius must not diverge"
  "Preferably it should remain finite and monotonic, or approach a finite limiting value"

So areal radius -> 0 is acceptable (it's finite). The key is it doesn't DIVERGE.
And the user wants "continuous regular space" — which means no singularity,
not that the space has to stay large.

Let me try alpha = 0.75 (between 0.5 and 1.0):
  - A ~ r^{-0.75}, so A^2 ~ r^{-1.5}
  - rho = A*r ~ r^{0.25} -> 0 (finite, doesn't diverge)
  - K ~ r^{12*0.75 - 6} = r^{3} -> 0 (curvature vanishes!)
  - det_2d ~ -A^4 ~ -r^{-3} (negative, Lorentzian)

This could work! The areal radius goes to 0 but doesn't diverge.
The curvature vanishes. The determinant stays negative.
And with q != 0, we get temporal stretching from the disformal term.

To get A ~ r^{-0.75} instead of r^{-2}, we need phi_0 = 0.75 (since A ~ r^{phi_0}).
With beta_A = -1: A = exp(-phi) = exp(|phi|) ~ (r_h/r)^{phi_0} = r^{-phi_0}.
So phi_0 = 0.75 gives A ~ r^{-0.75}.

But we also need to check geodesic completeness:
  null integrand ~ sqrt(-det_2d) ~ A^2 ~ r^{-1.5}
  integral of r^{-1.5} dr ~ r^{-0.5} -> infinity as r -> 0. COMPLETE!

  timelike: tau ~ A * sqrt(r/(2M)) dr ~ r^{-0.75} * r^{0.5} dr = r^{-0.25} dr
  integral of r^{-0.25} dr ~ r^{0.75} -> 0. NOT COMPLETE for timelike!

Hmm, timelike geodesics would be incomplete. Need alpha > 1 for timelike
completeness (tau ~ r^{0.5-alpha} dr, integral diverges if alpha >= 1.5).

Wait, let me recalculate. For the conformal prototype with phi_0 = 2:
  null: integrand ~ A^2 ~ r^{-4}, integral ~ r^{-3} -> inf. Complete.
  timelike: tau ~ A * sqrt(r/(2M)) ~ r^{-2} * r^{0.5} = r^{-1.5}, integral ~ r^{-0.5} -> inf. Complete.

For phi_0 = 0.75:
  null: integrand ~ A^2 ~ r^{-1.5}, integral ~ r^{-0.5} -> inf. Complete.
  timelike: tau ~ A * r^{0.5} ~ r^{-0.75+0.5} = r^{-0.25}, integral ~ r^{0.75} -> 0. INCOMPLETE.

For timelike completeness: need integral of r^{0.5-alpha} dr to diverge
  => 0.5 - alpha <= -1 => alpha >= 1.5

But alpha >= 1.5 means areal radius rho = A*r ~ r^{1-alpha} -> r^{-0.5} -> infinity.
That's the conformal prototype again!

So there's a FUNDAMENTAL TENSION:
  - Areal radius finite: alpha < 1
  - Timelike completeness: alpha >= 1.5
  These are incompatible with pure conformal.

BUT: with the disformal term (q != 0), the geodesic integrands change!
The null integrand is sqrt(-det_2d), and with B*q^2 contributing:
  det_2d = -A^2 * [A^2 + B*(F*psi'^2 + 2*q*psi')]
  If B*q^2 is large, gtilde_vv is large, and the determinant can be
  dominated by the disformal term rather than A^4.

Let me compute the actual geodesic integrands with the disformal term.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid

def compute_metric_v3(r, params):
    """Metric with mild conformal growth + disformal temporal effect."""
    r = np.asarray(r, dtype=float)
    M = params['M']
    q = params['q']
    phi_0 = params['phi_0']  # controls A ~ r^{-phi_0}
    delta = params['delta']
    beta_A = params.get('beta_A', -1.0)
    B0 = params['B0']
    n_B = params.get('n_B', 2.0)
    sigma_B = params.get('sigma_B', 1.5)

    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    # Scalar: psi(r) = phi_0 * ln(r/r_h) * S(r)  (logarithmic, same as prototype)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)
    psi = phi_0 * np.log(r_safe / r_h) * S
    dpsi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    # Conformal factor: A = exp(beta_A * psi) = exp(|psi|) for beta_A = -1
    # With phi_0 < 1, A grows mildly: A ~ r^{-phi_0}
    A = np.exp(beta_A * psi)
    A2 = A ** 2

    # Disformal: quartic Gaussian damped (same as prototype)
    # B activates near horizon, vanishes in deep interior
    abs_psi = np.abs(psi)
    B = B0 * abs_psi ** n_B / (1.0 + abs_psi ** n_B) * np.exp(-(psi ** 4) / (2.0 * sigma_B ** 4))

    # Metric components with phi = q*v + psi(r)
    gtilde_vv = -A2 * F + B * q ** 2
    gtilde_vr = A2 + B * q * dpsi
    gtilde_rr = B * dpsi ** 2
    gtilde_thth = A2 * r ** 2

    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2
    bracket = A2 + B * (F * dpsi ** 2 + 2 * q * dpsi)

    areal = A * r
    lorentzian = det_2d < 0

    return {
        'r': r, 'F': F, 'psi': psi, 'dpsi': dpsi,
        'A': A, 'B': B, 'A2': A2,
        'gtilde_vv': gtilde_vv, 'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr, 'gtilde_thth': gtilde_thth,
        'det_2d': det_2d, 'bracket': bracket,
        'areal_radius': areal, 'lorentzian': lorentzian,
    }


def analyze_geodesics_with_disformal(r, metric, M=1.0):
    """Compute actual geodesic integrands with the full disformal metric."""
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    det_2d = metric['det_2d']
    A = metric['A']

    disc = -det_2d
    disc_safe = np.where(disc > 0, disc, np.nan)
    sqrt_disc = np.sqrt(disc_safe)

    # Null affine parameter: lambda = int sqrt(-det_2d) dr
    null_integrand = sqrt_disc

    idx_h = np.argmin(np.abs(r - 2.0 * M))
    r_int = r[:idx_h + 1]
    null_int = null_integrand[:idx_h + 1]
    r_rev = r_int[::-1]
    int_rev = np.where(np.isfinite(null_int[::-1]), null_int[::-1], 0)
    lambda_null = cumulative_trapezoid(int_rev, r_rev, initial=0)
    lambda_total = float(lambda_null[-1]) if len(lambda_null) > 0 else 0.0

    # Power law
    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    int_inner = null_integrand[:n_check]
    mask = (int_inner > 0) & (r_inner > 0) & np.isfinite(int_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_int = np.log(int_inner[mask])
        alpha_null = -np.polyfit(log_r, log_int, 1)[0]
    else:
        alpha_null = np.nan

    # Timelike: for the full metric, the proper time integrand is more complex.
    # For radial infall with E=1 (from rest at infinity):
    #   The conserved energy: E = -(gvv * v' + gvr * r') = 1
    #   The normalization: gvv*v'^2 + 2*gvr*v'*r' + grr*r'^2 = -1
    #   Solving: v' = -(1 + gvr*r') / gvv
    #   Substituting into normalization:
    #   (1 + gvr*r')^2 / gvv - 2*(1 + gvr*r')*r' + grr*r'^2 = -1
    #   This is quadratic in r'. Let me solve it properly.

    # From E = -(gvv*v' + gvr*r') = 1:
    #   v' = -(1 + gvr*r') / gvv
    # Null condition: gvv*v'^2 + 2*gvr*v'*r' + grr*r'^2 = -1
    # Substituting v':
    #   gvv * (1 + gvr*r')^2 / gvv^2 - 2*gvr*(1+gvr*r')*r'/gvv + grr*r'^2 = -1
    #   (1 + gvr*r')^2 / gvv - 2*gvr*(1+gvr*r')*r'/gvv + grr*r'^2 = -1
    #   [1 + 2*gvr*r' + gvr^2*r'^2 - 2*gvr*r' - 2*gvr^2*r'^2] / gvv + grr*r'^2 = -1
    #   [1 - gvr^2*r'^2] / gvv + grr*r'^2 = -1
    #   1/gvv - gvr^2*r'^2/gvv + grr*r'^2 = -1
    #   r'^2 * (grr - gvr^2/gvv) = -1 - 1/gvv
    #   r'^2 = (-1 - 1/gvv) / (grr - gvr^2/gvv)
    #   r'^2 = (-(gvv + 1)/gvv) / ((grr*gvv - gvr^2)/gvv)
    #   r'^2 = -(gvv + 1) / det_2d
    #   Since det_2d < 0 (Lorentzian) and gvv < 0 outside, gvv > 0 inside:
    #   Outside: gvv+1 < 0 (gvv ~ -1), det_2d < 0, so r'^2 = negative/negative = positive ✓
    #   Inside: gvv > 0, gvv+1 > 0, det_2d < 0, so r'^2 = positive/negative = NEGATIVE!
    #
    # Wait, that can't be right. Inside the horizon, the energy E=1 from
    # rest at infinity doesn't make sense in the same way because r is
    # timelike. Let me reconsider.
    #
    # Actually, for the EF metric, the conserved energy E = -g_{va} x'^a is
    # always valid. Inside the horizon, the particle must move inward (r decreasing).
    # The issue is that E=1 corresponds to a particle at rest at infinity,
    # and inside the horizon, such a particle has dr/dtau < 0 (moving inward).
    #
    # Let me just compute r'^2 = -(gvv + 1) / det_2d numerically.
    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)
    rdot_sq = -(gvv + 1.0) / np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    rdot_sq = np.where(rdot_sq > 0, rdot_sq, np.nan)
    rdot = np.sqrt(rdot_sq)

    # dtau = dr / rdot
    # Integrate from horizon inward
    rdot_int = rdot[:idx_h + 1]
    rdot_rev = np.where(np.isfinite(rdot_int[::-1]), rdot_int[::-1], 1e-30)
    # dtau = dr / rdot, so tau = int dr/rdot
    dtau_dr = 1.0 / np.where(rdot_rev > 1e-30, rdot_rev, 1e-30)
    tau_inward = cumulative_trapezoid(dtau_dr, r_int[::-1], initial=0)
    tau_total = float(tau_inward[-1]) if len(tau_inward) > 0 else 0.0

    # Power law for timelike
    rdot_inner = rdot[:n_check]
    mask_t = (rdot_inner > 0) & (r_inner > 0) & np.isfinite(rdot_inner)
    if np.sum(mask_t) > 2:
        log_r = np.log(r_inner[mask_t])
        log_rdot = np.log(rdot_inner[mask_t])
        alpha_rdot = np.polyfit(log_r, log_rdot, 1)[0]  # rdot ~ r^{alpha_rdot}
        # dtau/dr ~ r^{-alpha_rdot}, integral diverges if alpha_rdot >= 1
        timelike_diverges = alpha_rdot >= 1.0
    else:
        alpha_rdot = np.nan
        timelike_diverges = False

    return {
        'null_affine_total': lambda_total,
        'null_power_law': float(alpha_null) if np.isfinite(alpha_null) else None,
        'null_diverges': bool(alpha_null >= 1.0) if np.isfinite(alpha_null) else False,
        'timelike_proper_time_total': tau_total,
        'timelike_power_law': float(alpha_rdot) if np.isfinite(alpha_rdot) else None,
        'timelike_diverges': bool(timelike_diverges),
        'null_integrand': null_integrand,
        'rdot': rdot,
    }


def test_mild_conformal():
    r = np.logspace(np.log10(1e-8), np.log10(50.0), 20000)

    print("=" * 80)
    print("MILD CONFORMAL + DISFORMAL TEMPORAL: phi_0 sweep")
    print("A ~ r^{-phi_0}, need phi_0 in (0.5, 1.5) for curvature vanishing")
    print("Areal radius rho = A*r ~ r^{1-phi_0}")
    print("=" * 80)

    for phi_0 in [0.5, 0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5]:
        params = {
            'M': 1.0, 'q': 0.1, 'phi_0': phi_0, 'delta': 0.05,
            'beta_A': -1.0, 'B0': 1.0, 'n_B': 2.0, 'sigma_B': 1.5,
        }
        metric = compute_metric_v3(r, params)
        geodesics = analyze_geodesics_with_disformal(r, metric, M=1.0)

        areal = metric['areal_radius']
        det = metric['det_2d']
        idx_h = np.argmin(np.abs(r - 2.0))

        frac_lor = np.sum(det < 0) / len(det)
        sign_changes = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]

        # Areal radius behavior
        areal_inner = areal[0]
        areal_max = np.max(areal)
        areal_diverges = areal_max > 1e6

        # Curvature (conformal estimate: K ~ A^{-12} * K_Schw ~ r^{12*phi_0 - 6})
        K_estimate_power = 12 * phi_0 - 6
        K_vanishes = K_estimate_power > 0  # i.e., phi_0 > 0.5

        print(f"\n  phi_0 = {phi_0:.2f}:")
        print(f"    Areal: innermost={areal_inner:.2e}, max={areal_max:.2e}, "
              f"diverges={areal_diverges}")
        print(f"    Lorentzian: frac={frac_lor:.6f}, sign_changes={len(sign_changes)}")
        print(f"    K estimate: ~r^{{{K_estimate_power:.1f}}}, vanishes={K_vanishes}")
        print(f"    Null: power={geodesics['null_power_law']}, "
              f"diverges={geodesics['null_diverges']}, "
              f"total={geodesics['null_affine_total']:.2e}")
        print(f"    Timelike: power={geodesics['timelike_power_law']}, "
              f"diverges={geodesics['timelike_diverges']}, "
              f"total={geodesics['timelike_proper_time_total']:.2e}")

        if len(sign_changes) > 0:
            print(f"    *** SIGN CHANGE at r={r[sign_changes[0]]:.4e} ***")

    # Now test the sweet spot: phi_0 = 0.75 with various q
    print("\n" + "=" * 80)
    print("SWEET SPOT: phi_0 = 0.75, varying q and B0")
    print("  A ~ r^{-0.75}, rho ~ r^{0.25} -> 0 (finite)")
    print("  K ~ r^{3} -> 0 (vanishes)")
    print("  Null: integrand ~ A^2 ~ r^{-1.5}, integral ~ r^{-0.5} -> inf (complete)")
    print("  Timelike: needs disformal help")
    print("=" * 80)

    for q in [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]:
        for B0 in [1.0, 10.0, 50.0]:
            params = {
                'M': 1.0, 'q': q, 'phi_0': 0.75, 'delta': 0.05,
                'beta_A': -1.0, 'B0': B0, 'n_B': 2.0, 'sigma_B': 1.5,
            }
            metric = compute_metric_v3(r, params)
            geodesics = analyze_geodesics_with_disformal(r, metric, M=1.0)

            det = metric['det_2d']
            frac_lor = np.sum(det < 0) / len(det)
            sign_changes = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]

            areal = metric['areal_radius']
            areal_max = np.max(areal)
            areal_div = areal_max > 1e6

            print(f"\n  q={q:5.1f}, B0={B0:5.1f}: "
                  f"lor_frac={frac_lor:.6f}, sc={len(sign_changes)}, "
                  f"areal_max={areal_max:.2e}, areal_div={areal_div}, "
                  f"null_div={geodesics['null_diverges']}, "
                  f"tau_div={geodesics['timelike_diverges']}, "
                  f"null_pwr={geodesics['null_power_law']}, "
                  f"tau_pwr={geodesics['timelike_power_law']}")

            if len(sign_changes) > 0:
                print(f"    *** SIGN CHANGE at r={r[sign_changes[0]]:.4e} ***")


if __name__ == '__main__':
    test_mild_conformal()
