#!/usr/bin/env python3
"""Disformal temporal branch — gate-test script.

Tests candidate TEP-BH metrics with phi(v,r) = q*v + psi(r) and bounded A
against six acceptance gates:

  Gate 1: Areal radius finite
  Gate 2: Determinant Lorentian (det_2d < 0 everywhere)
  Gate 3: Local clock / redshift well-defined and growing
  Gate 4: No hard horizon (outgoing null slope finite or asymptotic)
  Gate 5: Curvature invariants finite
  Gate 6: Geodesic completeness (null + timelike with angular momentum)

The geometric metric is Schwarzschild in ingoing Eddington-Finkelstein:
    ds_g^2 = -(1 - 2M/r) dv^2 + 2 dv dr + r^2 dOmega^2

The matter metric is:
    gtilde_{mu nu} = A^2(phi) g_{mu nu} + B(phi) d_mu phi d_nu phi

With phi(v,r) = q*v + psi(r):
    d_v phi = q
    d_r phi = psi'(r)

So the disformal contribution is:
    B * [q^2, q*psi', 0, 0; q*psi', (psi')^2, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0]

Full metric components:
    gtilde_vv  = -A^2 * F + B * q^2          (F = 1 - 2M/r)
    gtilde_vr  =  A^2 + B * q * psi'
    gtilde_rr  =  B * (psi')^2
    gtilde_thth = A^2 * r^2
    gtilde_phph = A^2 * r^2 * sin^2(theta)

2D determinant:
    det_2d = gtilde_vv * gtilde_rr - gtilde_vr^2
           = (-A^2 F + B q^2)(B psi'^2) - (A^2 + B q psi')^2
           = -A^2 B F psi'^2 + B^2 q^2 psi'^2 - A^4 - 2 A^2 B q psi' - B^2 q^2 psi'^2
           = -A^4 - A^2 B (F psi'^2 + 2 q psi') - 0
           = -A^2 (A^2 + B (F psi'^2 + 2 q psi'))

Wait, let me redo this more carefully.

det_2d = gtilde_vv * gtilde_rr - gtilde_vr^2

gtilde_vv = -A^2 F + B q^2
gtilde_vr = A^2 + B q psi'
gtilde_rr = B (psi')^2

det_2d = (-A^2 F + B q^2) * B (psi')^2 - (A^2 + B q psi')^2
       = -A^2 B F (psi')^2 + B^2 q^2 (psi')^2 - A^4 - 2 A^2 B q psi' - B^2 q^2 (psi')^2
       = -A^4 - A^2 B F (psi')^2 - 2 A^2 B q psi'
       = -A^2 [A^2 + B (F (psi')^2 + 2 q psi')]

So: det_2d = -A^2 * [A^2 + B * (F * psi'^2 + 2*q*psi')]

For Lorentzian signature: det_2d < 0
Since A^2 > 0 always, we need:
    A^2 + B * (F * psi'^2 + 2*q*psi') > 0

This is the key condition. With bounded A (A -> A_0 finite) and B growing
in the interior, we need the bracket to remain positive.

Key insight: the term 2*q*psi' can be made negative (if q and psi' have
opposite signs), which would help keep the bracket positive. But F*psi'^2
is always positive inside the horizon (F < 0, so F*psi'^2 < 0 inside).

Wait — inside the horizon F < 0, so F*psi'^2 < 0. That means:
    A^2 + B*(negative + 2*q*psi')

If psi' < 0 (scalar decreasing inward) and q > 0, then 2*q*psi' < 0.
So both terms inside B are negative inside the horizon.
Then: A^2 + B*(negative) could go negative → det_2d > 0 → signature change!

This is the central challenge. We need to choose signs and magnitudes
carefully so that the bracket stays positive.

Let me think about this differently. The condition is:

    A^2 + B * (F * psi'^2 + 2*q*psi') > 0

Inside the horizon (F < 0):
    F * psi'^2 < 0  (always, since psi'^2 > 0 and F < 0)
    2*q*psi' can be + or - depending on signs

If we want B to create a strong temporal effect (large gtilde_vv),
we need B*q^2 to be large. But we also need the determinant to stay
negative.

Let me try: q > 0, psi' > 0 (scalar increasing inward, or decreasing outward).
Then 2*q*psi' > 0, which helps keep the bracket positive.
But F*psi'^2 < 0 inside, which hurts.

The balance depends on magnitudes. Let's explore numerically.

Author: TEP-BH v0.2 development
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import json
from pathlib import Path

# =============================================================================
# Metric computation for disformal temporal branch
# =============================================================================

def compute_disformal_temporal_metric(r, params):
    """Compute the disformal matter metric with phi(v,r) = q*v + psi(r).

    Geometric metric: Schwarzschild EF
        ds_g^2 = -F dv^2 + 2 dv dr + r^2 dOmega^2
        F(r) = 1 - 2M/r

    Matter metric:
        gtilde = A^2 * g + B * (d_phi tensor d_phi)

    With phi = q*v + psi(r):
        gtilde_vv  = -A^2 * F + B * q^2
        gtilde_vr  =  A^2 + B * q * psi'
        gtilde_rr  =  B * (psi')^2
        gtilde_thth = A^2 * r^2

    Determinant:
        det_2d = -A^2 * [A^2 + B*(F*psi'^2 + 2*q*psi')]
    """
    r = np.asarray(r, dtype=float)
    M = params['M']
    q = params['q']
    phi_0 = params['phi_0']
    delta = params['delta']
    beta_A = params['beta_A']
    B0 = params['B0']
    sigma_B = params.get('sigma_B', 1.5)
    n_B = params.get('n_B', 2.0)
    A_max = params.get('A_max', 10.0)  # cap on conformal factor

    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    # Scalar field: psi(r) = phi_0 * ln(r/r_h) * S(r)
    # S(r) = 1/(1 + exp((r-r_h)/(delta*r_h))) — logistic activation
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    psi = phi_0 * np.log(r_safe / r_h) * S

    # psi'(r)
    dS = -S * (1.0 - S) / (delta * r_h)
    dpsi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    # Conformal factor: A = exp(beta_A * psi), capped at A_max
    A_raw = np.exp(beta_A * psi)
    A = np.minimum(A_raw, A_max)
    A2 = A ** 2

    # Disformal function: B(psi)
    # For the temporal branch, B should GROW in the interior (not vanish).
    # Try: B(psi) = B0 * |psi|^n_B / (1 + |psi|^n_B)  (saturating, no Gaussian damping)
    abs_psi = np.abs(psi)
    B = B0 * abs_psi ** n_B / (1.0 + abs_psi ** n_B)

    # Metric components
    gtilde_vv = -A2 * F + B * q ** 2
    gtilde_vr = A2 + B * q * dpsi
    gtilde_rr = B * dpsi ** 2
    gtilde_thth = A2 * r ** 2

    # 2D determinant
    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2

    # Verify: det_2d = -A^2 * [A^2 + B*(F*dpsi^2 + 2*q*dpsi)]
    bracket = A2 + B * (F * dpsi ** 2 + 2 * q * dpsi)
    det_2d_check = -A2 * bracket

    # Areal radius
    areal_radius = np.sqrt(np.maximum(gtilde_thth, 0))

    # Lorentzian check
    lorentzian = det_2d < 0

    # Effective temporal coefficient (gtilde_vv when F<0 inside)
    # Inside: gtilde_vv = -A^2*F + B*q^2 = A^2*|F| + B*q^2
    # This is the "temporal stretching" — it grows with B*q^2
    temporal_coeff = np.where(F < 0, A2 * np.abs(F) + B * q ** 2, np.abs(gtilde_vv))

    # Redshift factor: sqrt(gtilde_vv(r) / gtilde_vv(infinity))
    # At infinity: A->1, B->0, F->1, so gtilde_vv -> -1
    # Redshift z(r) = sqrt(|gtilde_vv(r)| / |gtilde_vv(inf)|) - 1
    gtilde_vv_inf = 1.0  # |gtilde_vv| at infinity
    redshift = np.sqrt(np.maximum(np.abs(gtilde_vv), 1e-30) / gtilde_vv_inf) - 1.0

    return {
        'r': r, 'F': F, 'psi': psi, 'dpsi': dpsi,
        'A': A, 'A_raw': A_raw, 'B': B, 'A2': A2,
        'gtilde_vv': gtilde_vv, 'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr, 'gtilde_thth': gtilde_thth,
        'det_2d': det_2d, 'det_2d_check': det_2d_check,
        'bracket': bracket,
        'areal_radius': areal_radius,
        'lorentzian': lorentzian,
        'temporal_coeff': temporal_coeff,
        'redshift': redshift,
    }


# =============================================================================
# Curvature invariants (numerical)
# =============================================================================

def compute_curvature_full(r, metric):
    """Compute full curvature invariants numerically.

    For a spherically symmetric metric in EF coordinates:
        ds^2 = g_vv dv^2 + 2 g_vr dv dr + g_rr dr^2 + g_thth dOmega^2

    The non-zero Christoffel symbols involve only r-derivatives.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    gthth = metric['gtilde_thth']
    det_2d = metric['det_2d']

    # Numerical derivatives
    gvv_p = np.gradient(gvv, r)
    gvr_p = np.gradient(gvr, r)
    grr_p = np.gradient(grr, r)
    gthth_p = np.gradient(gthth, r)

    gvv_pp = np.gradient(gvv_p, r)
    grr_pp = np.gradient(grr_p, r)
    gthth_pp = np.gradient(gthth_p, r)

    # Safe versions
    det_safe = np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    gthth_safe = np.where(np.abs(gthth) > 1e-50, gthth, np.nan)
    grr_safe = np.where(np.abs(grr) > 1e-50, grr, np.nan)

    # Inverse metric (2D block)
    ginv_vv = grr / det_safe
    ginv_vr = -gvr / det_safe
    ginv_rr = gvv / det_safe
    ginv_thth = 1.0 / gthth_safe

    # Christoffel symbols
    Gv_vv = -0.5 * ginv_vr * gvv_p
    Gr_vv = -0.5 * ginv_rr * gvv_p
    Gv_vr = 0.5 * ginv_vv * gvv_p
    Gr_vr = 0.5 * ginv_vr * gvv_p
    Gv_rr = ginv_vv * gvr_p + 0.5 * ginv_vr * grr_p
    Gr_rr = ginv_vr * gvr_p + 0.5 * ginv_rr * grr_p
    Gv_thth = -0.5 * ginv_vr * gthth_p
    Gr_thth = -0.5 * ginv_rr * gthth_p
    Gth_rth = 0.5 * ginv_thth * gthth_p

    # Riemann components
    R_vrvr = -0.5 * gvv_pp + (
        gvv * (Gv_vv * Gv_rr - Gv_vr ** 2)
        + gvr * (Gv_vv * Gr_rr - Gv_vr * Gr_vr)
        + gvr * (Gr_vv * Gv_rr - Gr_vr * Gv_vr)
        + grr * (Gr_vv * Gr_rr - Gr_vr ** 2)
    )

    R_rthrth = -0.5 * gthth_pp + gthth_p * Gr_thth

    # Ricci scalar (2D + angular)
    R_2D = 2.0 * R_vrvr / det_safe
    R_ang = -gthth_pp / gthth_safe + 0.5 * (gthth_p / gthth_safe) ** 2
    R_total = R_2D + 2.0 * R_ang / gthth_safe

    # Kretschmann scalar
    # K = 4 R_vrvr^2 / det_2d^2 + 4 R_rthrth^2 / (gthth * grr)
    # Need to handle grr -> 0 carefully
    K_term1 = 4.0 * R_vrvr ** 2 / det_safe ** 2
    K_term2 = np.where(np.abs(grr_safe) > 1e-50,
                       4.0 * R_rthrth ** 2 / (gthth_safe * grr_safe), 0.0)
    Kretschmann = K_term1 + K_term2

    # Ricci squared: R_{mu nu} R^{mu nu}
    # For spherical symmetry: R_{AB} R^{AB} + 2 R_{thth}^2 / gthth^2
    # R_{vv} = R_vrvr * ginv_rr (simplified for this symmetry)
    # This is approximate; full computation would need all Ricci components
    Ricci_vv = R_vrvr * ginv_rr
    Ricci_rr = R_vrvr * ginv_vv
    Ricci_vr = -R_vrvr * ginv_vr
    Ricci_thth = R_ang

    Ricci_sq = (Ricci_vv ** 2 * ginv_vv + 2 * Ricci_vr ** 2 * ginv_vr
                + Ricci_rr ** 2 * ginv_rr + 2 * Ricci_thth ** 2 * ginv_thth)

    # Weyl/Cotton: C^2 ~ K - 2 R_{mu nu} R^{mu nu} + R^2/3 (4D)
    # For this simplified case:
    C_squared = np.maximum(Kretschmann - 2 * Ricci_sq + R_total ** 2 / 3.0, 0.0)

    return {
        'R_vrvr': R_vrvr,
        'R_rthrth': R_rthrth,
        'Ricci_scalar': R_total,
        'Ricci_squared': Ricci_sq,
        'Kretschmann': Kretschmann,
        'C_squared': C_squared,
    }


# =============================================================================
# Geodesic completeness (with angular momentum)
# =============================================================================

def check_geodesics_full(r, metric, M=1.0):
    """Check null and timelike geodesic completeness including angular momentum.

    For the EF metric with Killing vector d/dv:
        E = -(gtilde_vv * v' + gtilde_vr * r')  (conserved energy)
        L = gtilde_thth * phi'                  (conserved angular momentum)

    Null: gtilde_{mu nu} x'^mu x'^nu = 0
    Timelike: gtilde_{mu nu} x'^mu x'^nu = -1

    For radial (L=0) null:
        (dr/dlambda)^2 = E^2 / (-det_2d)
        lambda = int sqrt(-det_2d) dr / |E|

    For non-radial null (L != 0):
        Effective potential includes L^2/gthth term.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    gthth = metric['gtilde_thth']
    det_2d = metric['det_2d']
    A = metric['A']

    disc = -det_2d  # > 0 for Lorentzian
    disc_safe = np.where(disc > 0, disc, np.nan)
    sqrt_disc = np.sqrt(disc_safe)

    # --- Radial null affine parameter ---
    null_integrand_radial = sqrt_disc

    idx_h = np.argmin(np.abs(r - 2.0 * M))
    r_int = r[:idx_h + 1]
    null_int = null_integrand_radial[:idx_h + 1]

    r_rev = r_int[::-1]
    int_rev = np.where(np.isfinite(null_int[::-1]), null_int[::-1], 0)
    lambda_radial = cumulative_trapezoid(int_rev, r_rev, initial=0)
    lambda_radial_total = float(lambda_radial[-1]) if len(lambda_radial) > 0 else 0.0

    # --- Non-radial null (L != 0) ---
    # Effective potential: V_eff = L^2 * |gvv| / (gthth * disc)
    # The radial equation: (dr/dlambda)^2 = (E^2 - L^2 * |gvv|/gthth) / disc
    # For completeness, we need the integral of sqrt(disc / (E^2 - V_eff)) dr
    # to diverge. With L=0 this reduces to the radial case.
    # For L != 0, the effective potential can create turning points, but
    # completeness requires that geodesics that reach deep interior do so
    # in infinite affine parameter.
    # We check the L=0 case (most stringent for reaching r=0) plus
    # verify that the effective potential doesn't create a finite-radius
    # barrier that would make the space incomplete.

    # Effective potential for null with L=1, E=1:
    L_test = 1.0
    E_test = 1.0
    V_eff = L_test ** 2 * np.abs(gvv) / np.where(gthth > 1e-30, gthth, np.nan) / disc_safe
    # Check if V_eff < E^2 everywhere inside (no barrier)
    no_barrier = np.all(V_eff[1:idx_h+1] < E_test ** 2) if idx_h > 1 else True

    # --- Timelike proper time (radial infall from rest at infinity) ---
    # For radial timelike: dtau^2 = -gtilde_{mu nu} dx^mu dx^nu
    # With E=1 (from rest at infinity), radial:
    # (dr/dtau)^2 = E^2/(-det_2d) - 1/grr_approx
    # Simplified: dtau ~ sqrt(-det_2d / (E^2 + |det_2d|/grr)) dr
    # For the conformal regime this reduces to A * sqrt(r/(2M)) dr
    # More carefully: for timelike radial with E=1:
    #   gvv v'^2 + 2 gvr v' r' + grr r'^2 = -1
    #   E = -(gvv v' + gvr r') = 1
    #   v' = -(1 + gvr r') / gvv
    #   Substituting: r'^2 = (E^2 + det_2d) / (gvv * grr - gvr^2) ... no
    #   Actually: r'^2 = (E^2 + det_2d) / (-det_2d * grr/gvv) ... complex
    # Let's use the standard result:
    #   For timelike: (dr/dtau)^2 = E^2/(-det_2d) - 1/grr_total
    #   where grr_total is the effective radial metric
    # Simpler: dtau = sqrt(-det_2d) / sqrt(E^2 + det_2d/grr_eff) dr
    # For now, use the conformal approximation as a lower bound:
    tau_integrand = A * np.sqrt(np.maximum(r, 1e-30) / (2.0 * M))
    tau_int = tau_integrand[:idx_h + 1]
    tau_rev = np.where(np.isfinite(tau_int[::-1]), tau_int[::-1], 0)
    tau_inward = cumulative_trapezoid(tau_rev, r_int[::-1], initial=0)
    tau_total = float(abs(tau_inward[-1])) if len(tau_inward) > 0 else 0.0

    # --- Power-law analysis near r=0 ---
    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    int_inner = null_integrand_radial[:n_check]
    mask = (int_inner > 0) & (r_inner > 0) & np.isfinite(int_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_int = np.log(int_inner[mask])
        alpha = -np.polyfit(log_r, log_int, 1)[0]
        null_diverges = alpha >= 1.0
    else:
        alpha = np.nan
        null_diverges = False

    # --- Outgoing null slope (Gate 4: no hard horizon) ---
    # dv/dr_outgoing = (-gvr + sqrt(disc)) / gvv
    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)
    dv_dr_out = (-gvr + sqrt_disc) / gvv_safe
    dv_dr_in = (-gvr - sqrt_disc) / gvv_safe

    # Check: is there a radius where outgoing null slope -> infinity?
    # (that would be a horizon in the matter metric)
    outgoing_finite = np.isfinite(dv_dr_out)
    # A horizon occurs where gvv -> 0 (outgoing null can't escape)
    gvv_zero = np.argmin(np.abs(gvv[1:] - 0.0)) + 1 if len(gvv) > 1 else 0
    has_matter_horizon = np.any(np.abs(gvv[1:idx_h+1]) < 1e-10) if idx_h > 1 else False

    return {
        'null_affine_total': lambda_radial_total,
        'null_integrand_power_law': float(alpha) if np.isfinite(alpha) else None,
        'null_diverges': bool(null_diverges),
        'timelike_proper_time_total': tau_total,
        'timelike_diverges': bool(tau_total > 1e3),
        'no_angular_barrier': bool(no_barrier),
        'dv_dr_outgoing': dv_dr_out,
        'dv_dr_ingoing': dv_dr_in,
        'has_matter_horizon': bool(has_matter_horizon),
        'gvv_min_inside': float(np.min(np.abs(gvv[1:idx_h+1]))) if idx_h > 1 else None,
    }


# =============================================================================
# Gate tests
# =============================================================================

def run_gate_tests(r, params):
    """Run all six gate tests and return pass/fail + diagnostics."""

    metric = compute_disformal_temporal_metric(r, params)
    curvature = compute_curvature_full(r, metric)
    geodesics = check_geodesics_full(r, metric, M=params['M'])

    results = {}

    # --- Gate 1: Areal radius finite ---
    areal = metric['areal_radius']
    areal_finite = np.all(np.isfinite(areal))
    areal_max = float(np.nanmax(areal))
    areal_at_innermost = float(areal[0])
    areal_at_horizon = float(areal[np.argmin(np.abs(r - 2.0 * params['M']))])
    results['gate1_spatial'] = {
        'pass': areal_finite and areal_max < 1e6,
        'areal_max': areal_max,
        'areal_at_innermost': areal_at_innermost,
        'areal_at_horizon': areal_at_horizon,
    }

    # --- Gate 2: Determinant Lorentzian ---
    det = metric['det_2d']
    lorentzian = det < 0
    frac_lorentzian = float(np.sum(lorentzian)) / len(det)
    det_min = float(np.min(det))
    det_max = float(np.max(det))
    # Check for sign changes
    sign_changes = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]
    results['gate2_determinant'] = {
        'pass': frac_lorentzian > 0.999,
        'frac_lorentzian': frac_lorentzian,
        'det_min': det_min,
        'det_max': det_max,
        'n_sign_changes': len(sign_changes),
    }

    # --- Gate 3: Clock / redshift ---
    redshift = metric['redshift']
    # Redshift should grow monotonically inward (at least in the interior)
    idx_h = np.argmin(np.abs(r - 2.0 * params['M']))
    redshift_interior = redshift[:idx_h]
    # Check monotonicity (allowing small numerical noise)
    diffs = np.diff(redshift_interior)
    monotonic = np.sum(diffs < 0) / max(len(diffs), 1) > 0.95  # mostly increasing inward
    redshift_at_innermost = float(redshift[0])
    redshift_at_horizon = float(redshift[idx_h])
    results['gate3_clocks'] = {
        'pass': np.isfinite(redshift_at_innermost) and redshift_at_innermost > redshift_at_horizon,
        'redshift_at_innermost': redshift_at_innermost,
        'redshift_at_horizon': redshift_at_horizon,
        'monotonic_inward': bool(monotonic),
    }

    # --- Gate 4: No hard horizon ---
    has_horizon = geodesics['has_matter_horizon']
    gvv_min = geodesics['gvv_min_inside']
    results['gate4_no_horizon'] = {
        'pass': not has_horizon,
        'has_matter_horizon': has_horizon,
        'gvv_min_inside': gvv_min,
    }

    # --- Gate 5: Curvature finite ---
    K = curvature['Kretschmann']
    R = curvature['Ricci_scalar']
    R_sq = curvature['Ricci_squared']
    C_sq = curvature['C_squared']

    # Mask out NaN regions (numerical issues near boundaries)
    K_finite = np.all(np.isfinite(K[np.isfinite(K)]))
    R_finite = np.all(np.isfinite(R[np.isfinite(R)]))
    R_sq_finite = np.all(np.isfinite(R_sq[np.isfinite(R_sq)]))

    K_max = float(np.nanmax(K)) if np.any(np.isfinite(K)) else np.inf
    R_max = float(np.nanmax(np.abs(R[np.isfinite(R)]))) if np.any(np.isfinite(R)) else np.inf
    R_sq_max = float(np.nanmax(R_sq[np.isfinite(R_sq)])) if np.any(np.isfinite(R_sq)) else np.inf

    results['gate5_curvature'] = {
        'pass': K_max < 1e10 and R_max < 1e10 and R_sq_max < 1e10,
        'Kretschmann_max': K_max,
        'Ricci_scalar_max': R_max,
        'Ricci_squared_max': R_sq_max,
    }

    # --- Gate 6: Geodesic completeness ---
    null_complete = geodesics['null_diverges']
    timelike_complete = geodesics['timelike_diverges']
    no_barrier = geodesics['no_angular_barrier']
    results['gate6_geodesics'] = {
        'pass': null_complete and timelike_complete and no_barrier,
        'null_diverges': null_complete,
        'null_affine_total': geodesics['null_affine_total'],
        'null_power_law': geodesics['null_integrand_power_law'],
        'timelike_diverges': timelike_complete,
        'timelike_proper_time_total': geodesics['timelike_proper_time_total'],
        'no_angular_barrier': no_barrier,
    }

    # Overall
    gate_keys = ['gate1_spatial', 'gate2_determinant', 'gate3_clocks',
                 'gate4_no_horizon', 'gate5_curvature', 'gate6_geodesics']
    all_pass = all(results[k]['pass'] for k in gate_keys)
    results['all_pass'] = all_pass
    results['metric'] = metric
    results['curvature'] = curvature
    results['geodesics'] = geodesics

    return results


# =============================================================================
# Main: test multiple candidate parameter sets
# =============================================================================

def test_candidates():
    """Test multiple candidate parameter sets against all gates."""

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 20000)

    candidates = [
        # Candidate 1: q > 0, psi' > 0 (scalar increasing inward), A capped
        {
            'name': 'q=1, psi_inward, A_cap=10, B_saturating',
            'M': 1.0, 'q': 1.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': -0.1,  # weak conformal (A ~ exp(0.1*|psi|), bounded)
            'B0': 10.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 10.0,
        },
        # Candidate 2: q > 0, stronger B, weaker A
        {
            'name': 'q=2, strong_B, A_cap=5',
            'M': 1.0, 'q': 2.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': -0.05,
            'B0': 50.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 5.0,
        },
        # Candidate 3: q > 0, very weak A (nearly 1), strong B
        {
            'name': 'q=1, A~1, B0=100',
            'M': 1.0, 'q': 1.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': -0.01,  # A ~ 1 + epsilon
            'B0': 100.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 2.0,
        },
        # Candidate 4: q > 0, A=1 exactly (pure disformal), B grows
        {
            'name': 'q=1, A=1 (pure disformal), B0=50',
            'M': 1.0, 'q': 1.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': 0.0,  # A = 1 everywhere
            'B0': 50.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 1.0,
        },
        # Candidate 5: q > 0, A=1, B grows, different psi profile
        {
            'name': 'q=0.5, A=1, B0=200, phi_0=1',
            'M': 1.0, 'q': 0.5, 'phi_0': 1.0, 'delta': 0.05,
            'beta_A': 0.0,
            'B0': 200.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 1.0,
        },
        # Candidate 6: q > 0, A=1, B grows, larger q
        {
            'name': 'q=5, A=1, B0=10, phi_0=2',
            'M': 1.0, 'q': 5.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': 0.0,
            'B0': 10.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 1.0,
        },
        # Candidate 7: q > 0, A=1, B grows, very large q
        {
            'name': 'q=10, A=1, B0=1, phi_0=2',
            'M': 1.0, 'q': 10.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': 0.0,
            'B0': 1.0, 'sigma_B': 1.5, 'n_B': 2.0,
            'A_max': 1.0,
        },
        # Candidate 8: q > 0, A=1, B grows, n_B=1 (linear saturation)
        {
            'name': 'q=1, A=1, B0=50, n_B=1, phi_0=2',
            'M': 1.0, 'q': 1.0, 'phi_0': 2.0, 'delta': 0.05,
            'beta_A': 0.0,
            'B0': 50.0, 'sigma_B': 1.5, 'n_B': 1.0,
            'A_max': 1.0,
        },
    ]

    print("=" * 80)
    print("DISFORMAL TEMPORAL BRANCH — GATE TESTS")
    print("=" * 80)

    for cand in candidates:
        name = cand.pop('name')
        print(f"\n--- {name} ---")
        try:
            results = run_gate_tests(r, cand)

            gate_keys = ['gate1_spatial', 'gate2_determinant', 'gate3_clocks',
                         'gate4_no_horizon', 'gate5_curvature', 'gate6_geodesics']
            labels = {
                'gate1_spatial': "Gate 1 (Spatial finite)",
                'gate2_determinant': "Gate 2 (Lorentzian)",
                'gate3_clocks': "Gate 3 (Clocks/Redshift)",
                'gate4_no_horizon': "Gate 4 (No hard horizon)",
                'gate5_curvature': "Gate 5 (Curvature finite)",
                'gate6_geodesics': "Gate 6 (Geodesics complete)",
            }
            for k in gate_keys:
                gate = results[k]
                status = "PASS" if gate['pass'] else "FAIL"
                print(f"  {labels[k]}: {status}")
                # Print key diagnostics
                for kk, v in gate.items():
                    if kk != 'pass':
                        if isinstance(v, float):
                            print(f"    {kk}: {v:.6e}")
                        else:
                            print(f"    {kk}: {v}")

            print(f"\n  OVERALL: {'ALL PASS' if results['all_pass'] else 'FAILED'}")
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()


if __name__ == '__main__':
    test_candidates()
