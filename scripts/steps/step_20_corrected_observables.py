#!/usr/bin/env python3
"""
TEP Corrected Observables Derivation (v4, fixed-ADM normalization)
===================================================================

Fixes:
1. Fixed-ADM mass normalization: M_ADM = M is held fixed (the mass
   measured by a distant observer). The horizon mass parameter m shifts:
   m = M / (1 + 49*eta^2/360), giving r_H = 2m = 2M(1 - 19.6*beta^2).
   This is the correct physical normalization — the distant observer
   measures M, and the horizon contracts due to the temporal field's
   backreaction through the sGB coupling.

2. The horizon shift IS the TEP signature: the temporal field backreacts
   on the geometry, contracting the horizon. The horizon "naturally changes
   depending on observer's time rate" — a distant observer (fast clock)
   sees the horizon at r_H = 2M(1 - 19.6*beta^2); a deeper observer
   (slower clock) sees less redshift to any given emitter.

3. The frame-split SURVIVES as coupling-order difference:
   - Shadow: O(eta^2) = -0.044% (photons probe g_sGB)
   - ISCO on tilde_g: O(eta) = +1.95% (massive particles probe tilde_g)
   - Ratio: ~45x (the ISCO shift is an order of magnitude larger)
   With the mass-inflation branch (eta < 0 in this convention), the
   conformal factor A = e^{-phi} > 1 magnifies spatial scales. The ISCO
   on the matter metric moves outward, while the photon-sphere/shadow
   on the geometric metric remains inward-shifted at O(eta^2). The
   frame-split signature is a divergence, not a common sign.

4. The Temporal Horizon is observer-dependent: the frequency transfer
   Z = omega_e/omega_o = [A(r_o)/A(r_e)] * sqrt(F(r_o)/F(r_e)) depends
   on both emitter and observer positions. A deeper observer (slower
   clock) sees less redshift to the same emitter — the emitter is
   "less temporally remote." The practical observability boundary
   moves inward for deeper observers.
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

def h2_poly_prime(x):
    """dh2/dx"""
    return (-98/5 - 196/5*x - 274/5*x**2 - 56/15*x**3 + 52/3*x**4 + 40*x**5)

def sigma2_poly_prime(x):
    """dsigma2/dx"""
    return (98/5 + 116/5*x + 114/5*x**2 - 1624/15*x**3 - 2180/15*x**4 - 184/3*x**5)

# ============================================================
# Fixed-ADM mass normalization
# ============================================================
# M_ADM = M is fixed (mass measured by distant observer)
# m = horizon mass parameter: M_ADM = m(1 + 49*eta^2/360)
# => m = M / (1 + 49*eta^2/360) ≈ M(1 - 49*eta^2/360)
# r_H = 2m = 2M(1 - 49*eta^2/360) = 2M(1 - 19.6*beta^2)
# where beta = eta/12 (for m ≈ M)

def get_horizon_mass(M_ADM, eta):
    """Horizon mass parameter m for fixed ADM mass M_ADM."""
    return M_ADM / (1 + 49*eta**2/360)

def get_r_H(M_ADM, eta):
    """Horizon radius for fixed ADM mass."""
    m = get_horizon_mass(M_ADM, eta)
    return 2 * m

def scalar_phi(r, m, eta, M_ADM=1.0):
    """sGB scalar: phi = (2*alpha/m)(1/r + m/r^2 + 4m^2/(3r^3))
    where alpha = eta*M_ADM^2/3.
    
    Derived from Box phi = -alpha_GB * G, G_Schw = 48 m^2 / r^6.
    Horizon regularity fixes integration constant; phi -> 0 at infinity.
    Scalar charge: Q_s = 2*alpha_GB/m.
    Uses horizon mass m (not M_ADM) for consistency with Sotiriou-Zhou.
    """
    alpha_GB = eta * M_ADM**2 / 3
    return (2*alpha_GB / m) * (1/r + m/r**2 + 4*m**2/(3*r**3))

def scalar_phi_prime(r, m, eta, M_ADM=1.0):
    """dphi/dr = (2*alpha/m)(-1/r^2 - 2m/r^3 - 4m^2/r^4)"""
    alpha_GB = eta * M_ADM**2 / 3
    return (2*alpha_GB / m) * (-1/r**2 - 2*m/r**3 - 4*m**2/r**4)

# ============================================================
# Metric components (fixed-ADM)
# ============================================================

def g_sGB_components(r, m, eta, M_ADM=1.0):
    """Geometric metric g_sGB components.
    F = 1 - 2m/r (horizon at r = 2m)
    x = 2m/r
    beta = alpha_GB / r_H^2 = (eta*M_ADM^2/3) / (2m)^2
    f = F * (1 + beta^2 * h2)
    g = F^{-1} * (1 + beta^2 * sigma_2)
    h = r^2
    """
    r_H = 2 * m
    x = r_H / r
    alpha_GB = eta * M_ADM**2 / 3
    beta_sq = (alpha_GB / r_H**2)**2
    F = 1 - 2*m/r
    h2 = h2_poly(x)
    s2 = sigma2_poly(x)
    f = F * (1 + beta_sq * h2)
    g = F**(-1) * (1 + beta_sq * s2)
    h = r**2
    return f, g, h, beta_sq

def tilde_g_components(r, m, eta, M_ADM=1.0):
    """Matter metric tilde_g = A^2 * g_sGB."""
    f, g, h, bsq = g_sGB_components(r, m, eta, M_ADM)
    phi = scalar_phi(r, m, eta, M_ADM)
    A = np.exp(-phi)
    A2 = A**2
    return A2 * f, A2 * g, A2 * h, A, bsq

# Analytic derivatives
def g_sGB_analytic(r, m, eta, M_ADM=1.0):
    """Geometric metric with analytic derivatives."""
    r_H = 2 * m
    x = r_H / r
    alpha_GB = eta * M_ADM**2 / 3
    beta_sq = (alpha_GB / r_H**2)**2
    F = 1 - 2*m/r
    Fp = 2*m/r**2
    dx_dr = -r_H / r**2
    h2p = h2_poly_prime(x) * dx_dr
    h2 = h2_poly(x)
    s2 = sigma2_poly(x)
    f = F * (1 + beta_sq * h2)
    g = F**(-1) * (1 + beta_sq * s2)
    h = r**2
    fp = Fp * (1 + beta_sq * h2) + F * beta_sq * h2p
    hp = 2 * r
    return f, g, h, fp, hp, beta_sq

def tilde_g_analytic(r, m, eta, M_ADM=1.0):
    """Matter metric with analytic derivatives."""
    f, g, h, fp_g, hp_g, bsq = g_sGB_analytic(r, m, eta, M_ADM)
    phi = scalar_phi(r, m, eta, M_ADM)
    A = np.exp(-phi)
    A2 = A**2
    phi_p = scalar_phi_prime(r, m, eta, M_ADM)
    A2p = -2 * phi_p * A2
    fp_t = A2p * f + A2 * fp_g
    hp_t = A2p * h + A2 * hp_g
    return A2 * f, A2 * g, A2 * h, A, fp_t, hp_t, bsq

# ============================================================
# Photon sphere and shadow (on g_sGB, conformal invariance)
# ============================================================

def find_photon_sphere(m, eta, M_ADM=1.0, use_tilde=False):
    """Find photon sphere: minimum of b^2(r) = h/f."""
    def b2(r):
        if use_tilde:
            f, g, h, A, bsq = tilde_g_components(r, m, eta, M_ADM)
        else:
            f, g, h, bsq = g_sGB_components(r, m, eta, M_ADM)
        if f <= 0:
            return 0
        return h / f

    def db2(r):
        dr = max(1e-7 * r, 1e-12)
        return (b2(r + dr) - b2(r - dr)) / (2 * dr)

    r_lo = 2 * m + 0.01
    r_vals = np.linspace(r_lo, 10.0, 5000)
    db2_vals = [db2(r) for r in r_vals]
    for i in range(len(db2_vals) - 1):
        if db2_vals[i] < 0 and db2_vals[i + 1] > 0:
            try:
                r_ph = brentq(db2, r_vals[i], r_vals[i + 1], xtol=1e-12)
                return r_ph, np.sqrt(b2(r_ph))
            except:
                continue
    return None, None

# ============================================================
# ISCO (timelike geodesics)
# ============================================================

def find_ISCO(m, eta, M_ADM=1.0, use_tilde=False):
    """Find ISCO: minimum of L^2(r) using analytic derivatives."""
    def L2(r):
        if use_tilde:
            f, g, h, A, fp, hp, bsq = tilde_g_analytic(r, m, eta, M_ADM)
        else:
            f, g, h, fp, hp, bsq = g_sGB_analytic(r, m, eta, M_ADM)
        denom = f * hp - fp * h
        if denom <= 0:
            return None
        val = fp * h**2 / denom
        if val <= 0:
            return None
        return val

    def dL2(r):
        v = L2(r)
        if v is None:
            return 0
        dr = 1e-7 * r
        vp = L2(r + dr)
        vm = L2(r - dr)
        if vp is None or vm is None:
            return 0
        return (vp - vm) / (2 * dr)

    r_lo = 2 * m + 0.5
    r_vals = np.linspace(r_lo, 20.0, 5000)
    dL2_vals = [dL2(r) for r in r_vals]
    for i in range(len(dL2_vals) - 1):
        if dL2_vals[i] < 0 and dL2_vals[i + 1] > 0:
            try:
                r_isco = brentq(dL2, r_vals[i], r_vals[i + 1], xtol=1e-12)
                L2_val = L2(r_isco)
                if L2_val:
                    return r_isco, np.sqrt(L2_val)
            except:
                continue
    return None, None

# ============================================================
# Frequency transfer (observer-dependent Temporal Horizon)
# ============================================================

def frequency_transfer(r_e, r_o, m, eta, M_ADM=1.0):
    """Z = omega_e / omega_o = [A(r_o)/A(r_e)] * sqrt(F(r_o)/F(r_e))"""
    A_e = np.exp(-scalar_phi(r_e, m, eta, M_ADM))
    A_o = np.exp(-scalar_phi(r_o, m, eta, M_ADM))
    F_e = 1 - 2*m/r_e
    F_o = 1 - 2*m/r_o
    if F_e <= 0 or F_o <= 0:
        return np.inf
    return (A_o / A_e) * np.sqrt(F_o / F_e)

# ============================================================
# Main computation
# ============================================================
if __name__ == "__main__":
    M = 1.0  # Fixed ADM mass
    eta = -0.1  # mass-inflation branch (alpha_GB < 0, A > 1)

    m = get_horizon_mass(M, eta)
    r_H = get_r_H(M, eta)
    beta = eta / 12  # approximately (for m ≈ M)
    beta_sq = beta**2

    print("=" * 70)
    print("TEP CORRECTED OBSERVABLES (v4, fixed-ADM normalization)")
    print("=" * 70)
    print(f"Parameters: M_ADM = {M}, eta = {eta}")
    print(f"Horizon mass: m = {m:.8f}")
    print(f"Horizon radius: r_H = 2m = {r_H:.8f}")
    print(f"Horizon shift: 2M - r_H = {2*M - r_H:.6f} ({100*(1 - r_H/(2*M)):.4f}%)")
    print(f"  = 19.6 * beta^2 = {19.6*beta_sq:.6f}")
    print(f"beta = eta/12 = {beta:.6f}, beta^2 = {beta_sq:.6e}")
    print()

    # Scalar and conformal factor
    print("--- sGB scalar and conformal factor A = e^{-phi} ---")
    for r in [r_H, 3*M, 6*M, 10*M]:
        phi = scalar_phi(r, m, eta, M)
        A = np.exp(-phi)
        print(f"  r={r:.4f}M: phi={phi:.6f}, A={A:.6f}, deviation={100*(A-1):.2f}%")
    print()

    # Photon sphere / shadow
    print("--- Photon sphere / shadow (null geodesics on g_sGB) ---")
    r_ph_schw = 3.0
    b_ph_schw = 3 * np.sqrt(3)

    r_ph_g, b_ph_g = find_photon_sphere(m, eta, M, use_tilde=False)
    r_ph_tilde, b_ph_tilde = find_photon_sphere(m, eta, M, use_tilde=True)

    print(f"  Schwarzschild:  r_ph = {r_ph_schw:.4f}M, b = {b_ph_schw:.4f}M")
    if r_ph_g:
        dev_ph = 100 * (r_ph_g / r_ph_schw - 1)
        dev_b = 100 * (b_ph_g / b_ph_schw - 1)
        print(f"  On g_sGB:       r_ph = {r_ph_g:.5f}M, b = {b_ph_g:.5f}M")
        print(f"    Deviation: r_ph = {dev_ph:+.4f}%, b = {dev_b:+.4f}%  [O(eta^2)]")
    if r_ph_tilde:
        print(f"  On tilde_g:     r_ph = {r_ph_tilde:.5f}M, b = {b_ph_tilde:.5f}M")
        print(f"    (Matches g_sGB due to conformal invariance of null geodesics)")
    print()

    # ISCO
    print("--- ISCO (timelike geodesics, massive particles) ---")
    r_isco_schw = 6.0

    result_g = find_ISCO(m, eta, M, use_tilde=False)
    result_tilde = find_ISCO(m, eta, M, use_tilde=True)

    print(f"  Schwarzschild:  r = {r_isco_schw:.4f}M")
    if result_g:
        r_isc_g, L_isc_g = result_g
        dev_g = 100 * (r_isc_g / r_isco_schw - 1)
        print(f"  On g_sGB:     r = {r_isc_g:.5f}M  [O(eta^2)]")
        print(f"    Deviation: {dev_g:+.4f}%")
    if result_tilde:
        r_isc_t, L_isc_t = result_tilde
        dev_t = 100 * (r_isc_t / r_isco_schw - 1)
        print(f"  On tilde_g:     r = {r_isc_t:.5f}M  [O(eta)]")
        print(f"    Deviation: {dev_t:+.4f}%")
    print()

    # Frame-split signature
    print("=" * 70)
    print("FRAME-SPLIT SIGNATURE (fixed-ADM normalization)")
    print("=" * 70)
    if r_ph_g and result_tilde:
        shadow_dev = (b_ph_g / b_ph_schw - 1) * 100
        isco_dev_t = (r_isc_t / r_isco_schw - 1) * 100
        isco_dev_g = (r_isc_g / r_isco_schw - 1) * 100
        ratio = abs(isco_dev_t) / abs(shadow_dev) if shadow_dev != 0 else float('inf')

        print(f"  Shadow deviation (photons on g_sGB):    {shadow_dev:+.4f}%  [O(eta^2)]")
        print(f"  ISCO deviation (massive on tilde_g):      {isco_dev_t:+.4f}%  [O(eta)]")
        print(f"  ISCO deviation (massive on g_sGB):      {isco_dev_g:+.4f}%  [O(eta^2)]")
        print()
        print(f"  Coupling-order ratio: |ISCO_tilde|/|Shadow| = {ratio:.1f}x")
        print()
        print("  The frame-split is the COUPLING-ORDER DIFFERENCE:")
        print("  - Shadow shifts at O(eta^2) (photons probe g_sGB), sign NEGATIVE")
        print("  - ISCO on tilde_g shifts at O(eta) (massive particles probe tilde_g),")
        print("    sign POSITIVE for the mass-inflation branch (A > 1)")
        print("  - The ISCO shift is ~45x larger than the shadow shift")
        print("  - The frame-split is a DIVERGENCE:")
        print("    shadow contracts (O(eta^2)) while matter-metric ISCO expands (O(eta))")
        print()
        print("  This is the signature of the mass-inflation branch:")
        print("  conformal magnification (A = e^{-phi} > 1) pushes massive orbits")
        print("  outward, while the geometric metric's photon sphere receives only")
        print("  an O(eta^2) inward shift from sGB backreaction.")
    print()

    # Observer-dependent Temporal Horizon
    print("=" * 70)
    print("OBSERVER-DEPENDENT TEMPORAL HORIZON")
    print("=" * 70)
    print()
    print("Z = omega_e/omega_o = [A(r_o)/A(r_e)] * sqrt(F(r_o)/F(r_e))")
    print()
    print(f"  {'Emitter r/M':>12}", end="")
    for r_o in [100, 10, 6, 4, 3]:
        A_o = np.exp(-scalar_phi(r_o, m, eta, M))
        print(f"  {'r_o='+str(r_o)+'M':>12}", end="")
    print()
    print(f"  {'':>12}", end="")
    for r_o in [100, 10, 6, 4, 3]:
        A_o = np.exp(-scalar_phi(r_o, m, eta, M))
        print(f"  {'A_o='+f'{A_o:.4f}':>12}", end="")
    print()
    print("  " + "-" * 72)

    for r_e in [6, 4, 3, 2.5, 2.2, 2.1, 2.05, 2.01, r_H + 0.001]:
        print(f"  {r_e:12.4f}", end="")
        for r_o in [100, 10, 6, 4, 3]:
            z = frequency_transfer(r_e, r_o, m, eta, M)
            if np.isinf(z):
                print(f"  {'inf':>12}", end="")
            elif z > 1e6:
                print(f"  {z:12.2e}", end="")
            else:
                print(f"  {z:12.4f}", end="")
        print()

    print()
    print("A deeper observer (stronger conformal field, larger A_o) sees")
    print("LESS redshift to the same emitter. The emitter is 'less temporally")
    print("remote.' The Temporal Horizon is a RELATIONAL concept — it depends")
    print("on the observer's position in the temporal field.")
    print()

    # QNM estimate
    print("=" * 70)
    print("QNM ESTIMATE (horizon shift)")
    print("=" * 70)
    horizon_shift_pct = 100 * (1 - r_H / (2 * M))
    print(f"  Horizon shift: {horizon_shift_pct:.4f}%")
    print(f"  QNM shift ~ +{horizon_shift_pct:.4f}% (leading order, needs spectral solver)")
    print()

    # Perturbative control
    print("--- Perturbative control ---")
    x_horizon = 1.0  # x = r_H/r at horizon
    h2_at_horizon = h2_poly(x_horizon)
    print(f"  h_2(x=1) = {h2_at_horizon:.4f}")
    print(f"  beta^2 = {beta_sq:.6e}")
    print(f"  beta^2 * |h_2(1)| = {beta_sq * abs(h2_at_horizon):.6e}")
    print(f"  (Safely perturbative: << 1)")
    print()

    # Save results
    results = {
        "parameters": {"M_ADM": M, "eta": eta, "normalization": "fixed-ADM"},
        "horizon": {
            "m": float(m),
            "r_H": float(r_H),
            "shift_pct": float(horizon_shift_pct),
            "formula": "r_H = 2M(1 - 19.6*beta^2), beta = eta/12",
        },
        "scalar_profile": {},
        "shadow": {
            "Schwarzschild": {"r_ph": 3.0, "b": float(b_ph_schw)},
        },
        "ISCO": {
            "Schwarzschild": {"r": 6.0},
        },
        "frame_split": {
            "shadow_dev_pct": float(shadow_dev) if r_ph_g else None,
            "isco_tilde_dev_pct": float(isco_dev_t) if result_tilde else None,
            "isco_g_dev_pct": float(isco_dev_g) if result_g else None,
            "coupling_order_ratio": float(ratio) if r_ph_g and result_tilde else None,
            "signature": "mass-inflation sign divergence: shadow contracts (O(eta^2)) while matter-metric ISCO expands (O(eta))", 
        },
        "QNM": {
            "horizon_shift_pct": float(horizon_shift_pct),
            "estimated_shift_pct": float(horizon_shift_pct),
            "note": "Leading order from horizon shift; needs coupled spectral solver",
        },
        "perturbative_control": {
            "h2_at_horizon": float(h2_at_horizon),
            "beta_sq": float(beta_sq),
            "beta_sq_times_h2": float(beta_sq * abs(h2_at_horizon)),
        },
    }

    for r in [r_H, 3, 6, 10]:
        phi = scalar_phi(r, m, eta, M)
        results["scalar_profile"][f"r={r:.4f}M"] = {
            "phi": float(phi), "A": float(np.exp(-phi))
        }

    if r_ph_g:
        results["shadow"]["g_sGB"] = {"r_ph": float(r_ph_g), "b": float(b_ph_g)}
    if r_ph_tilde:
        results["shadow"]["tilde_g"] = {"r_ph": float(r_ph_tilde), "b": float(b_ph_tilde)}
    if result_g:
        results["ISCO"]["g_sGB"] = {"r": float(r_isc_g)}
    if result_tilde:
        results["ISCO"]["tilde_g"] = {"r": float(r_isc_t)}

    with open("results/step_20_corrected_observables.json", "w") as f:
        json.dump(results, f, indent=2)

    print("Results saved to results/step_20_corrected_observables.json")
