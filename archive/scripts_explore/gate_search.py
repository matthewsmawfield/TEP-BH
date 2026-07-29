#!/usr/bin/env python3
"""Refined search: saturating psi with smooth activation and controlled psi'.

The problem: near the horizon, psi' is large (activation zone), and
F*psi'^2 + 2*q*psi' can go negative, making the determinant flip sign.

Solution approaches:
1. Use a smoother activation (larger delta) to reduce psi' near horizon
2. Use a psi profile where psi' is always small
3. Make B activate more gradually
4. Use a B that depends on r directly (not just psi) to control where it's active

Key: we need B large in the DEEP interior (for temporal effect) but
small in the transition zone near the horizon (to keep determinant negative).
This is the OPPOSITE of the conformal prototype, where B was large near
the horizon and small in the deep interior.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import json

def compute_metric_v2(r, params):
    """Compute metric with saturating psi and controlled B activation."""

    r = np.asarray(r, dtype=float)
    M = params['M']
    q = params['q']
    phi_0 = params['phi_0']
    delta = params['delta']
    A0 = params.get('A0', 1.0)  # constant conformal factor
    B0 = params['B0']
    n_B = params.get('n_B', 2.0)

    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    # Saturating psi: tanh activation
    # psi(r) = phi_0 * tanh((r_h - r)/(delta*r_h)) * S(r)
    # S(r) = 1/(1 + exp((r-r_h)/(delta*r_h)))
    x = np.maximum((r_h - r) / (delta * r_h), 0)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)

    tanh_x = np.tanh(x)
    sech2_x = 1.0 / np.cosh(x) ** 2

    psi = phi_0 * tanh_x * S
    dpsi = phi_0 * (-sech2_x / (delta * r_h)) * S + phi_0 * tanh_x * dS

    # Conformal factor: constant A0
    A = np.full_like(r, A0)
    A2 = A ** 2

    # Disformal function: B grows in the deep interior
    # Use B(psi) = B0 * |psi|^n_B / (1 + |psi|^n_B)  (saturating)
    # With saturating psi, |psi| -> phi_0 as r -> 0, so B -> B0 * phi_0^n/(1+phi_0^n)
    # This is large in the deep interior and small near the horizon.
    abs_psi = np.abs(psi)
    B = B0 * abs_psi ** n_B / (1.0 + abs_psi ** n_B)

    # Metric components
    gtilde_vv = -A2 * F + B * q ** 2
    gtilde_vr = A2 + B * q * dpsi
    gtilde_rr = B * dpsi ** 2
    gtilde_thth = A2 * r ** 2

    # 2D determinant
    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2
    bracket = A2 + B * (F * dpsi ** 2 + 2 * q * dpsi)

    # Areal radius
    areal_radius = A * r  # = A0 * r (finite!)

    # Lorentzian
    lorentzian = det_2d < 0

    # Redshift
    gtilde_vv_inf = A0 ** 2  # |gtilde_vv| at infinity
    redshift = np.sqrt(np.maximum(np.abs(gtilde_vv), 1e-30) / gtilde_vv_inf) - 1.0

    return {
        'r': r, 'F': F, 'psi': psi, 'dpsi': dpsi,
        'A': A, 'B': B, 'A2': A2,
        'gtilde_vv': gtilde_vv, 'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr, 'gtilde_thth': gtilde_thth,
        'det_2d': det_2d, 'bracket': bracket,
        'areal_radius': areal_radius,
        'lorentzian': lorentzian,
        'redshift': redshift,
    }


def compute_curvature(r, metric):
    """Compute curvature invariants numerically."""
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    gthth = metric['gtilde_thth']
    det_2d = metric['det_2d']

    gvv_p = np.gradient(gvv, r)
    gvr_p = np.gradient(gvr, r)
    grr_p = np.gradient(grr, r)
    gthth_p = np.gradient(gthth, r)
    gvv_pp = np.gradient(gvv_p, r)
    grr_pp = np.gradient(grr_p, r)
    gthth_pp = np.gradient(gthth_p, r)

    det_safe = np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    gthth_safe = np.where(np.abs(gthth) > 1e-50, gthth, np.nan)
    grr_safe = np.where(np.abs(grr) > 1e-50, grr, np.nan)

    ginv_vv = grr / det_safe
    ginv_vr = -gvr / det_safe
    ginv_rr = gvv / det_safe
    ginv_thth = 1.0 / gthth_safe

    Gv_vv = -0.5 * ginv_vr * gvv_p
    Gr_vv = -0.5 * ginv_rr * gvv_p
    Gv_vr = 0.5 * ginv_vv * gvv_p
    Gr_vr = 0.5 * ginv_vr * gvv_p
    Gv_rr = ginv_vv * gvr_p + 0.5 * ginv_vr * grr_p
    Gr_rr = ginv_vr * gvr_p + 0.5 * ginv_rr * grr_p
    Gv_thth = -0.5 * ginv_vr * gthth_p
    Gr_thth = -0.5 * ginv_rr * gthth_p
    Gth_rth = 0.5 * ginv_thth * gthth_p

    R_vrvr = -0.5 * gvv_pp + (
        gvv * (Gv_vv * Gv_rr - Gv_vr ** 2)
        + gvr * (Gv_vv * Gr_rr - Gv_vr * Gr_vr)
        + gvr * (Gr_vv * Gv_rr - Gr_vr * Gv_vr)
        + grr * (Gr_vv * Gr_rr - Gr_vr ** 2)
    )
    R_rthrth = -0.5 * gthth_pp + gthth_p * Gr_thth

    R_2D = 2.0 * R_vrvr / det_safe
    R_ang = -gthth_pp / gthth_safe + 0.5 * (gthth_p / gthth_safe) ** 2
    R_total = R_2D + 2.0 * R_ang / gthth_safe

    K_term1 = 4.0 * R_vrvr ** 2 / det_safe ** 2
    K_term2 = np.where(np.abs(grr_safe) > 1e-50,
                       4.0 * R_rthrth ** 2 / (gthth_safe * grr_safe), 0.0)
    Kretschmann = K_term1 + K_term2

    return {
        'Ricci_scalar': R_total,
        'Kretschmann': Kretschmann,
        'R_vrvr': R_vrvr,
        'R_rthrth': R_rthrth,
    }


def check_geodesics(r, metric, M=1.0):
    """Check geodesic completeness."""
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    det_2d = metric['det_2d']
    A = metric['A']

    disc = -det_2d
    disc_safe = np.where(disc > 0, disc, np.nan)
    sqrt_disc = np.sqrt(disc_safe)

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
        alpha = -np.polyfit(log_r, log_int, 1)[0]
        null_diverges = alpha >= 1.0
    else:
        alpha = np.nan
        null_diverges = False

    # Timelike (conformal approx with A0)
    A0 = A[0] if A[0] == A[-1] else 1.0
    tau_integrand = A0 * np.sqrt(np.maximum(r, 1e-30) / (2.0 * M))
    tau_int = tau_integrand[:idx_h + 1]
    tau_rev = np.where(np.isfinite(tau_int[::-1]), tau_int[::-1], 0)
    tau_inward = cumulative_trapezoid(tau_rev, r_int[::-1], initial=0)
    tau_total = float(abs(tau_inward[-1])) if len(tau_inward) > 0 else 0.0

    # Matter horizon check
    gvv_inside = gvv[1:idx_h+1]
    has_matter_horizon = np.any(np.abs(gvv_inside) < 1e-10)

    return {
        'null_affine_total': lambda_total,
        'null_power_law': float(alpha) if np.isfinite(alpha) else None,
        'null_diverges': bool(null_diverges),
        'timelike_proper_time_total': tau_total,
        'timelike_diverges': bool(tau_total > 1e3),
        'has_matter_horizon': bool(has_matter_horizon),
    }


def run_gates(r, params):
    """Run all gates."""
    metric = compute_metric_v2(r, params)
    curvature = compute_curvature(r, metric)
    geodesics = check_geodesics(r, metric, M=params['M'])

    idx_h = np.argmin(np.abs(r - 2.0 * params['M']))

    # Gate 1: Areal radius finite
    areal = metric['areal_radius']
    g1 = np.all(np.isfinite(areal)) and np.max(areal) < 1e6

    # Gate 2: Lorentzian
    det = metric['det_2d']
    frac_lor = np.sum(det < 0) / len(det)
    sign_changes = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]
    g2 = frac_lor > 0.999 and len(sign_changes) == 0

    # Gate 3: Redshift
    redshift = metric['redshift']
    g3 = np.isfinite(redshift[0]) and redshift[0] > redshift[idx_h]

    # Gate 4: No hard horizon
    g4 = not geodesics['has_matter_horizon']

    # Gate 5: Curvature finite
    K = curvature['Kretschmann']
    R = curvature['Ricci_scalar']
    K_max = float(np.nanmax(K[np.isfinite(K)])) if np.any(np.isfinite(K)) else np.inf
    R_max = float(np.nanmax(np.abs(R[np.isfinite(R)]))) if np.any(np.isfinite(R)) else np.inf
    g5 = K_max < 1e10 and R_max < 1e10

    # Gate 6: Geodesics
    g6 = geodesics['null_diverges'] and geodesics['timelike_diverges']

    all_pass = g1 and g2 and g3 and g4 and g5 and g6

    return {
        'g1': g1, 'g2': g2, 'g3': g3, 'g4': g4, 'g5': g5, 'g6': g6,
        'all_pass': all_pass,
        'frac_lorentzian': frac_lor,
        'n_sign_changes': len(sign_changes),
        'areal_max': float(np.max(areal)),
        'areal_innermost': float(areal[0]),
        'redshift_innermost': float(redshift[0]),
        'redshift_horizon': float(redshift[idx_h]),
        'K_max': K_max,
        'R_max': R_max,
        'null_total': geodesics['null_affine_total'],
        'null_power': geodesics['null_power_law'],
        'null_diverges': geodesics['null_diverges'],
        'tau_total': geodesics['timelike_proper_time_total'],
        'has_horizon': geodesics['has_matter_horizon'],
        'det_innermost': float(det[0]),
        'det_horizon': float(det[idx_h]),
        'bracket_innermost': float(metric['bracket'][0]),
        'bracket_horizon': float(metric['bracket'][idx_h]),
        'B_innermost': float(metric['B'][0]),
        'B_horizon': float(metric['B'][idx_h]),
        'psi_innermost': float(metric['psi'][0]),
        'psi_horizon': float(metric['psi'][idx_h]),
        'dpsi_innermost': float(metric['dpsi'][0]),
        'dpsi_horizon': float(metric['dpsi'][idx_h]),
    }


def search():
    r = np.logspace(np.log10(1e-8), np.log10(50.0), 20000)

    # Strategy: vary delta (smoothness of activation), q, B0
    # The key problem is psi' near the horizon. Larger delta = smoother = smaller psi'.
    candidates = []

    for delta in [0.05, 0.1, 0.2, 0.5, 1.0]:
        for q in [0.1, 0.5, 1.0, 5.0]:
            for B0 in [1.0, 10.0, 50.0, 100.0]:
                for phi_0 in [0.5, 1.0, 2.0]:
                    candidates.append({
                        'M': 1.0, 'q': q, 'phi_0': phi_0, 'delta': delta,
                        'A0': 1.0, 'B0': B0, 'n_B': 2.0,
                    })

    print(f"Testing {len(candidates)} candidates...")
    print()

    passes = []
    best = None
    best_frac = 0

    for i, cand in enumerate(candidates):
        try:
            res = run_gates(r, cand)
            if res['frac_lorentzian'] > best_frac:
                best_frac = res['frac_lorentzian']
                best = (cand, res)
            if res['all_pass']:
                passes.append((cand, res))
                print(f"  PASS #{len(passes)}: q={cand['q']}, delta={cand['delta']}, "
                      f"B0={cand['B0']}, phi_0={cand['phi_0']}")
            elif res['g2'] and res['g1'] and res['g3']:
                # At least passes spatial + Lorentzian + clocks
                print(f"  NEAR-MISS: q={cand['q']}, delta={cand['delta']}, "
                      f"B0={cand['B0']}, phi_0={cand['phi_0']}")
                print(f"    frac_lor={res['frac_lorentzian']:.6f}, "
                      f"K_max={res['K_max']:.2e}, null_div={res['null_diverges']}, "
                      f"tau={res['tau_total']:.2e}")
        except Exception as e:
            pass

    print(f"\n{'='*80}")
    print(f"RESULTS: {len(passes)} full passes out of {len(candidates)} candidates")
    print(f"{'='*80}")

    if passes:
        for cand, res in passes[:5]:
            print(f"\n  PASS: q={cand['q']}, delta={cand['delta']}, B0={cand['B0']}, phi_0={cand['phi_0']}")
            for k, v in res.items():
                if k not in ['g1','g2','g3','g4','g5','g6','all_pass']:
                    if isinstance(v, float):
                        print(f"    {k}: {v:.6e}")
                    else:
                        print(f"    {k}: {v}")
    else:
        print(f"\n  No full passes. Best Lorentzian fraction: {best_frac:.6f}")
        if best:
            cand, res = best
            print(f"  Best: q={cand['q']}, delta={cand['delta']}, B0={cand['B0']}, phi_0={cand['phi_0']}")
            for k, v in res.items():
                if k not in ['g1','g2','g3','g4','g5','g6','all_pass']:
                    if isinstance(v, float):
                        print(f"    {k}: {v:.6e}")
                    else:
                        print(f"    {k}: {v}")


if __name__ == '__main__':
    search()
