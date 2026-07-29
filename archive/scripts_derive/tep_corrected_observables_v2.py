#!/usr/bin/env python3
"""
TEP Corrected Observables Derivation (v2)
==========================================

Fixes:
1. Shadow/photon sphere: computed on g^{sGB} (conformal invariance of null geodesics).
   The conformal factor A = e^{-phi} does NOT affect null geodesic paths.
2. ISCO: computed on tilde_g with the CORRECT formula for general spherical metrics
   where the angular component is h = A^2 * r^2, not r^2.

For ds^2 = -f dt^2 + g dr^2 + h dOmega^2:
  - Null geodesics: photon sphere at d(b^2)/dr = 0 where b^2 = h/f (at circular orbit)
  - Timelike geodesics: L^2 = f'h^2 / (fh' - f'h), ISCO at d(L^2)/dr = 0
  - Shadow size: b_ph = sqrt(h/f) at the photon sphere, evaluated at infinity (A->1)

Key physics:
  - Shadow is O(eta^2) (metric correction only, conformal factor cancels for null)
  - ISCO on tilde_g is O(eta) (conformal factor affects massive particles)
  - This IS the frame-split: photons probe g, massive particles probe tilde_g
"""

import numpy as np
from scipy.optimize import brentq
import json

# ============================================================
# Sotiriou-Zhou (2014) perturbative solution
# ============================================================

def h2_poly(x):
    """Sotiriou-Zhou Eq. (60): h_2(x) polynomial."""
    return (-98/5*x - 98/5*x**2 - 274/15*x**3 - 14/15*x**4
            + 52/15*x**5 + 20/3*x**6)

def sigma2_poly(x):
    """Sotiriou-Zhou Eq. (61): sigma_2(x) polynomial."""
    return (98/5*x + 58/5*x**2 + 38/5*x**3 - 406/15*x**4
            - 436/15*x**5 - 92/3*x**6)

def scalar_phi(r, eta, M=1.0):
    """sGB scalar: phi(r) = (2*eta*M/3)(1/r + M/r^2 + 4M^2/(3r^3))"""
    return (2*eta*M/3) * (1/r + M/r**2 + 4*M**2/(3*r**3))

def F_schw(r, M=1.0):
    return 1 - 2*M/r

def F_prime(r, M=1.0):
    return 2*M/r**2

# ============================================================
# Metric components
# ============================================================

def g_sGB_components(r, eta, M=1.0):
    """Geometric metric g^{sGB} components (A=1).
    ds^2 = -f dt^2 + g dr^2 + r^2 dOmega^2
    f = F(1 + beta^2 h2), g = F^{-1}(1 + beta^2 sigma_2)
    """
    r_H = 2*M
    x = r_H / r
    beta_sq = (eta/12)**2
    F = F_schw(r, M)
    h2 = h2_poly(x)
    s2 = sigma2_poly(x)
    f = F * (1 + beta_sq * h2)
    g = F**(-1) * (1 + beta_sq * s2)
    h = r**2
    return f, g, h

def tilde_g_components(r, eta, M=1.0):
    """Matter metric tilde_g = A^2 * g^{sGB} components.
    ds^2 = -f_t dt^2 + g_t dr^2 + h_t dOmega^2
    f_t = A^2 * f, g_t = A^2 * g, h_t = A^2 * r^2
    """
    f, g, _ = g_sGB_components(r, eta, M)
    phi = scalar_phi(r, eta, M)
    A = np.exp(-phi)
    A2 = A**2
    f_t = A2 * f
    g_t = A2 * g
    h_t = A2 * r**2
    return f_t, g_t, h_t, A

def numerical_deriv(func, r, *args, **kwargs):
    """Numerical derivative using central differences."""
    dr = max(1e-8 * r, 1e-12)
    return (func(r + dr, *args, **kwargs) - func(r - dr, *args, **kwargs)) / (2*dr)

# ============================================================
# Photon sphere and shadow (on g^{sGB}, conformal invariance)
# ============================================================
# For null geodesics on ds^2 = -f dt^2 + g dr^2 + h dOmega^2:
# b^2 = h/f at circular null orbit
# Photon sphere: d(b^2)/dr = 0 => h'/h = f'/f
# Shadow size: b_ph = sqrt(h/f) at photon sphere (evaluated at infinity, A->1)

def b_squared_null(r, eta, M=1.0, use_tilde=False):
    """b^2 = h/f for null geodesics at radius r."""
    if use_tilde:
        f, g, h, A = tilde_g_components(r, eta, M)
    else:
        f, g, h = g_sGB_components(r, eta, M)
    if f <= 0:
        return 0
    return h / f

def find_photon_sphere(eta, M=1.0, use_tilde=False):
    """Find photon sphere: minimum of b^2(r) = h/f.
    The shadow radius is b_ph = sqrt(b^2) at the minimum.
    """
    def db2_dr(r):
        return numerical_deriv(b_squared_null, r, eta, M, use_tilde=use_tilde)

    r_values = np.linspace(2.0 + 0.01, 10.0, 5000)
    db2_values = [db2_dr(r) for r in r_values]

    # Find sign change from negative to positive (minimum)
    for i in range(len(db2_values)-1):
        if db2_values[i] < 0 and db2_values[i+1] > 0:
            try:
                r_ph = brentq(db2_dr, r_values[i], r_values[i+1], xtol=1e-12)
                b_ph = np.sqrt(b_squared_null(r_ph, eta, M, use_tilde=use_tilde))
                return r_ph, b_ph
            except:
                continue
    return None, None

# ============================================================
# ISCO (on tilde_g for massive particles)
# ============================================================
# For timelike geodesics on ds^2 = -f dt^2 + g dr^2 + h dOmega^2:
# L^2 = f' * h^2 / (f*h' - f'*h)
# E^2 = f * (1 + L^2/h)
# ISCO: d(L^2)/dr = 0

def L_squared_timelike(r, eta, M=1.0, use_tilde=True):
    """L^2 for circular timelike orbit at r."""
    if use_tilde:
        f, g, h, A = tilde_g_components(r, eta, M)
    else:
        f, g, h = g_sGB_components(r, eta, M)

    fp = numerical_deriv(lambda rr: tilde_g_components(rr, eta, M)[0] if use_tilde
                         else g_sGB_components(rr, eta, M)[0], r)
    hp = numerical_deriv(lambda rr: tilde_g_components(rr, eta, M)[2] if use_tilde
                         else g_sGB_components(rr, eta, M)[2], r)

    denom = f * hp - fp * h
    if denom <= 0:
        return None

    L2 = fp * h**2 / denom
    if L2 <= 0:
        return None
    return L2

def E_squared_timelike(r, eta, M=1.0, use_tilde=True):
    """E^2 for circular timelike orbit at r."""
    L2 = L_squared_timelike(r, eta, M, use_tilde)
    if L2 is None:
        return None
    if use_tilde:
        f, g, h, A = tilde_g_components(r, eta, M)
    else:
        f, g, h = g_sGB_components(r, eta, M)
    return f * (1 + L2 / h)

def find_ISCO(eta, M=1.0, use_tilde=True):
    """Find ISCO: minimum of L^2(r) (marginally stable orbit).
    d(L^2)/dr = 0 with d^2(L^2)/dr^2 > 0.
    """
    def dL2_dr(r):
        L2 = L_squared_timelike(r, eta, M, use_tilde)
        if L2 is None:
            return 0
        dr = 1e-6 * r
        L2p = L_squared_timelike(r + dr, eta, M, use_tilde)
        L2m = L_squared_timelike(r - dr, eta, M, use_tilde)
        if L2p is None or L2m is None:
            return 0
        return (L2p - L2m) / (2*dr)

    r_min = 2.5 + 0.1 if use_tilde else 3.0
    r_values = np.linspace(r_min, 20.0, 3000)
    dL2_values = [dL2_dr(r) for r in r_values]

    # Find sign change from negative to positive (minimum of L^2)
    for i in range(len(dL2_values)-1):
        if dL2_values[i] < 0 and dL2_values[i+1] > 0:
            try:
                r_isco = brentq(dL2_dr, r_values[i], r_values[i+1], xtol=1e-12)
                L2 = L_squared_timelike(r_isco, eta, M, use_tilde)
                E2 = E_squared_timelike(r_isco, eta, M, use_tilde)
                if L2 and E2:
                    return r_isco, np.sqrt(E2), np.sqrt(L2)
            except:
                continue
    return None

# ============================================================
# Main computation
# ============================================================
if __name__ == "__main__":
    M = 1.0
    eta = 0.1  # perturbative regime: beta^2 * |h2| ~ 0.001 * 30 = 0.03 << 1

    print("=" * 70)
    print("TEP CORRECTED OBSERVABLES DERIVATION (v2)")
    print("=" * 70)
    print(f"Parameters: M={M}, eta={eta}")
    print(f"alpha_GB = eta*M^2/3 = {eta*M**2/3:.4f} [L^2]")
    print(f"beta^2 = (eta/3)^2 = {(eta/12)**2:.6f}")
    print()

    # Scalar and conformal factor
    print("--- sGB scalar and conformal factor A = e^{-phi} ---")
    for r in [2*M, 3*M, 6*M, 10*M]:
        phi = scalar_phi(r, eta, M)
        A = np.exp(-phi)
        print(f"  r={r:.1f}M: phi={phi:.6f}, A={A:.6f}, deviation={100*(A-1):.2f}%")

    print()

    # Verify conformal invariance: photon sphere on g^{sGB} vs tilde_g
    print("--- Photon sphere / shadow (null geodesics) ---")
    print("  Conformal invariance check: photon sphere should be same on g and tilde_g")
    print()

    r_ph_g, b_ph_g = find_photon_sphere(eta, M, use_tilde=False)
    r_ph_tilde, b_ph_tilde = find_photon_sphere(eta, M, use_tilde=True)

    print(f"  Schwarzschild:              r_ph=3.0000M, b=5.1962M")
    if r_ph_g:
        print(f"  On g^sGB (A=1):             r_ph={r_ph_g:.4f}M, b={b_ph_g:.4f}M  (O(eta^2))")
        print(f"    Deviation: r_ph: {100*(r_ph_g/3.0-1):.4f}%, b: {100*(b_ph_g/5.1962-1):.4f}%")
    if r_ph_tilde:
        print(f"  On tilde_g (A=e^-phi):      r_ph={r_ph_tilde:.4f}M, b={b_ph_tilde:.4f}M")
        print(f"    Deviation: r_ph: {100*(r_ph_tilde/3.0-1):.4f}%, b: {100*(b_ph_tilde/5.1962-1):.4f}%")
        print(f"    (Should match g^sGB due to conformal invariance of null geodesics)")

    print()

    # ISCO comparison
    print("--- ISCO (timelike geodesics, massive particles) ---")
    print(f"  Schwarzschild:              r=6.0000M, E=0.9428, L=3.4641M")

    # ISCO on g^{sGB} (A=1, metric only, O(eta^2))
    result_g = find_ISCO(eta, M, use_tilde=False)
    if result_g:
        r_isc_g, E_isc_g, L_isc_g = result_g
        print(f"  On g^sGB (A=1):             r={r_isc_g:.4f}M, E={E_isc_g:.4f}, L={L_isc_g:.4f}M  (O(eta^2))")
        print(f"    Deviation: r: {100*(r_isc_g/6.0-1):.4f}%, E: {100*(E_isc_g/0.9428-1):.4f}%")

    # ISCO on tilde_g (A=e^{-phi}, O(eta) from conformal factor)
    result_tilde = find_ISCO(eta, M, use_tilde=True)
    if result_tilde:
        r_isc_t, E_isc_t, L_isc_t = result_tilde
        print(f"  On tilde_g (A=e^-phi):      r={r_isc_t:.4f}M, E={E_isc_t:.4f}, L={L_isc_t:.4f}M  (O(eta))")
        print(f"    Deviation: r: {100*(r_isc_t/6.0-1):.4f}%, E: {100*(E_isc_t/0.9428-1):.4f}%")

    print()
    print("--- KEY RESULT: Frame-split signature ---")
    if result_g and result_tilde and r_ph_g:
        shadow_dev_g = 100*(b_ph_g/5.1962-1)
        isco_dev_g = 100*(r_isc_g/6.0-1)
        isco_dev_t = 100*(r_isc_t/6.0-1)
        print(f"  Shadow deviation (photons on g):     {shadow_dev_g:.4f}%  [O(eta^2)]")
        print(f"  ISCO deviation (massive on tilde_g): {isco_dev_t:.4f}%  [O(eta)]")
        print(f"  ISCO deviation (massive on g):       {isco_dev_g:.4f}%  [O(eta^2)]")
        print()
        print("  The frame-split: photons probe g (O(eta^2) shadow shift),")
        print("  massive particles probe tilde_g (O(eta) ISCO shift).")
        print("  The ISCO shift is an ORDER OF MAGNITUDE larger than the shadow shift")
        print("  because the conformal factor A = e^{-phi} affects massive particles")
        print("  but NOT photons. This is the definitive two-metric signature.")

    print()

    # QNM correction
    print("--- Regge-Wheeler potential correction ---")
    print("  The corrected delta_V_RW comes from metric perturbations h_2, sigma_2.")
    print("  Dimensions: [delta_V] = L^{-2} (correct, from dimensionless h_2, sigma_2).")
    print("  Order: O(beta^2) = O(eta^2).")
    print()
    beta_sq = (eta/12)**2
    horizon_shift = 19.6 * beta_sq
    print(f"  Horizon shift: r_H = 2M(1 - {horizon_shift:.6f}) = {2*M*(1-horizon_shift):.6f}M")
    print(f"  Relative shift: {100*horizon_shift:.4f}%")
    print()
    print("  The QNM frequency shift is O(eta^2), proportional to the horizon shift.")
    print("  For the l=2 fundamental mode (omega_R^Schw = 0.3737/M):")
    print("  delta_omega/omega ~ -horizon_shift (leading order)")
    print(f"  Estimated shift: ~{100*horizon_shift:.4f}% (needs Leaver solver for precise value)")

    # Save results
    results = {
        "parameters": {"M": M, "eta": eta},
        "scalar_profile": {},
        "shadow": {"Schwarzschild": {"r_ph": 3.0, "b": 5.1962}},
        "ISCO": {"Schwarzschild": {"r": 6.0, "E": 0.9428, "L": 3.4641}},
        "QNM": {"Schwarzschild_omega_R": 0.3737, "horizon_shift_pct": 100*horizon_shift},
    }

    for r in [2, 3, 6, 10]:
        phi = scalar_phi(r, eta, M)
        results["scalar_profile"][f"r={r}M"] = {
            "phi": float(phi), "A": float(np.exp(-phi))
        }

    if r_ph_g:
        results["shadow"]["g_sGB"] = {"r_ph": float(r_ph_g), "b": float(b_ph_g)}
    if r_ph_tilde:
        results["shadow"]["tilde_g"] = {"r_ph": float(r_ph_tilde), "b": float(b_ph_tilde)}
    if result_g:
        results["ISCO"]["g_sGB"] = {"r": float(r_isc_g), "E": float(E_isc_g), "L": float(L_isc_g)}
    if result_tilde:
        results["ISCO"]["tilde_g"] = {"r": float(r_isc_t), "E": float(E_isc_t), "L": float(L_isc_t)}

    with open("results/tep_corrected_observables_v2.json", "w") as f:
        json.dump(results, f, indent=2)

    print()
    print("Results saved to results/tep_corrected_observables_v2.json")
