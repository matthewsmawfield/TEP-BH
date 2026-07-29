#!/usr/bin/env python3
"""Track A: QNM solver on the horizon-bearing shift-symmetric sGB branch.

Uses first-order perturbation theory with the exact Schwarzschild QNM
from the `qnm` package (Leaver's continued fraction method) and the
sGB correction to the potential from Sotiriou & Zhou (2014).

The first-order QNM perturbation formula is:
  δω = ∫ Ψ_0² δV dr* / (2ω_0 ∫ Ψ_0² dr*)

where Ψ_0 is the unperturbed QNM wave function, ω_0 is the unperturbed
QNM frequency, and δV is the sGB correction to the potential.

The wave function Ψ_0 is computed by integrating the Schwarzschild wave
equation from the peak of the potential (r ≈ 3M) outward in both
directions, which is numerically stable near the peak.

For the coupled polar-scalar system, the correction also includes the
mixing potential V_mix, which gives an additional O(α) correction to
the polar mode (the scalar-led mode gets an O(α) correction).

Outputs:
  results/step_27_qnm_horizon_branch.json
"""

from __future__ import annotations
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
import qnm

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _PROJECT_ROOT)
sys.path.insert(0, os.path.join(_PROJECT_ROOT, "scripts"))

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")

M = 1.0
L = 2
ETA_VALUES = [-0.05, -0.1, -0.15, -0.2]


# ---------------------------------------------------------------------------
# Sotiriou-Zhou (2014) perturbative sGB metric
# ---------------------------------------------------------------------------
def sgb_F(r, eta, M=1.0):
    r_h0 = 2.0 * M
    beta = eta / (3.0 * r_h0**2)
    beta2 = beta**2
    r_h = r_h0 * (1.0 - 19.6 * beta2)
    r_safe = np.maximum(r, r_h * 1.0001)
    F_schw = 1.0 - 2.0 * M / r_safe
    x = 2.0 * M / r_safe  # use Schwarzschild r_H = 2M for the polynomial
    h_2 = (-98.0/5.0 * x - 98.0/5.0 * x**2 - 274.0/15.0 * x**3
           - 14.0/15.0 * x**4 + 52.0/15.0 * x**5 + 20.0/3.0 * x**6)
    # F_corrected = (1 - r_H_sGB/r) * (1 + beta^2 * h_2)
    F = (1.0 - r_h / r_safe) * (1.0 + beta2 * h_2)
    return r_h, F


def sgb_F_prime(r, eta, M=1.0):
    r = np.asarray(r, dtype=float)
    dr = 1e-7 * np.maximum(np.abs(r), 0.1)
    _, F_plus = sgb_F(r + dr, eta, M)
    _, F_minus = sgb_F(r - dr, eta, M)
    return (F_plus - F_minus) / (2 * dr)


def gauss_bonnet_schw(r, M=1.0):
    return 48.0 * M**2 / np.maximum(r, 1e-30)**6


# ---------------------------------------------------------------------------
# Potentials
# ---------------------------------------------------------------------------
def V_RW_schw(r, M=1.0, l=2):
    """Regge-Wheeler potential for Schwarzschild."""
    F = 1.0 - 2.0 * M / r
    return F * (l*(l+1)/r**2 - 6.0*M/r**3)


def V_scalar_schw(r, M=1.0, l=2):
    """Scalar potential for Schwarzschild: V = F(l(l+1)/r^2 + F'/r)."""
    F = 1.0 - 2.0 * M / r
    Fp = 2.0 * M / r**2
    return F * (l*(l+1)/r**2 + Fp/r)


def V_RW_sgb(r, eta, M=1.0, l=2):
    """Regge-Wheeler potential with sGB corrections."""
    r_h, F = sgb_F(r, eta, M)
    m_eff = r_h / 2.0
    return F * (l*(l+1)/r**2 - 6.0*m_eff/r**3)


def V_scalar_sgb(r, eta, M=1.0, l=2):
    """Scalar potential with sGB corrections."""
    r_h, F = sgb_F(r, eta, M)
    Fp = sgb_F_prime(r, eta, M)
    return F * (l*(l+1)/r**2 + Fp/r)


def V_mix_sgb(r, eta, M=1.0, l=2):
    """Polar-scalar mixing from Gauss-Bonnet coupling."""
    r_h, F = sgb_F(r, eta, M)
    alpha_gb = eta * M**2 / 3.0
    ang = l * (l+1) * (l-1) * (l+2)
    G = gauss_bonnet_schw(r, M)
    return alpha_gb * ang * G * F * r**2 / (48.0 * M**2)


def delta_V_grav(r, eta, M=1.0, l=2):
    """sGB correction to the gravitational potential."""
    return V_RW_sgb(r, eta, M, l) - V_RW_schw(r, M, l)


def delta_V_scalar(r, eta, M=1.0, l=2):
    """sGB correction to the scalar potential."""
    return V_scalar_sgb(r, eta, M, l) - V_scalar_schw(r, M, l)


# ---------------------------------------------------------------------------
# Tortoise coordinate for Schwarzschild
# ---------------------------------------------------------------------------
def rstar_schw(r, M=1.0):
    """Schwarzschild tortoise coordinate r* = r + 2M*ln(r/(2M) - 1)."""
    return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0)


# ---------------------------------------------------------------------------
# Compute QNM wave function by integrating from the peak
# ---------------------------------------------------------------------------
def compute_qnm_wavefunction(omega, V_func, M=1.0, l=2, r_peak=None):
    """Compute the QNM wave function by integrating from the potential peak.

    The wave equation in r* is: d²Ψ/dr*² + (ω² - V)Ψ = 0

    We integrate from the peak of the potential outward in both directions.
    At the peak, we set Ψ = 1 and dΨ/dr* = 0 (arbitrary normalization).

    The wave function diverges at both boundaries, but the divergence is
    exponential and can be truncated for the perturbation integrals.
    """
    if r_peak is None:
        r_peak = 3.0 * M  # approximate peak of RW potential

    rs_peak = rstar_schw(r_peak, M)

    # Build r* grid
    r_inner = 2.01 * M
    r_outer = 50.0 * M
    rs_inner = rstar_schw(r_inner, M)
    rs_outer = rstar_schw(r_outer, M)

    # Grid from peak to inner (toward horizon)
    rs_grid_inner = np.linspace(rs_peak, rs_inner, 2000)
    # Grid from peak to outer (toward infinity)
    rs_grid_outer = np.linspace(rs_peak, rs_outer, 2000)

    # Interpolator r(r*)
    r_grid_fine = np.linspace(r_inner, r_outer, 5000)
    rs_grid_fine = rstar_schw(r_grid_fine, M)
    def r_of_rs(rs):
        return float(np.interp(rs, rs_grid_fine, r_grid_fine))

    # ODE: d²Ψ/dr*² = -(ω² - V) * Ψ
    def ode(rs, y):
        Psi, dPsi = y
        r = r_of_rs(rs)
        V = V_func(r)
        d2Psi = -(omega**2 - V) * Psi
        return [dPsi, d2Psi]

    # Integrate inward (from peak to horizon)
    sol_in = solve_ivp(ode, [rs_peak, rs_inner], [1.0, 0.0],
                       method='DOP853', rtol=1e-12, atol=1e-14,
                       t_eval=rs_grid_inner, max_step=0.1)
    # Integrate outward (from peak to infinity)
    sol_out = solve_ivp(ode, [rs_peak, rs_outer], [1.0, 0.0],
                        method='DOP853', rtol=1e-12, atol=1e-14,
                        t_eval=rs_grid_outer, max_step=0.1)

    if not sol_in.success or not sol_out.success:
        return None

    # Combine: inner part (reversed) + outer part
    rs_all = np.concatenate([sol_in.t[::-1], sol_out.t[1:]])
    Psi_all = np.concatenate([sol_in.y[0][::-1], sol_out.y[0][1:]])

    return rs_all, Psi_all


# ---------------------------------------------------------------------------
# Perturbative QNM correction
# ---------------------------------------------------------------------------
def compute_qnm_correction(omega_0, V_0_func, dV_func, M=1.0, l=2):
    """Compute first-order perturbative correction to QNM frequency.

    δω = ∫ Ψ_0² δV dr* / (2ω_0 ∫ Ψ_0² dr*)

    The integrals are truncated at finite distance to avoid the exponential
    divergence of the QNM wave function at the boundaries.
    """
    # Compute wave function
    result = compute_qnm_wavefunction(omega_0, V_0_func, M, l)
    if result is None:
        return None

    rs, Psi = result

    # Compute δV at the grid points
    r_grid_fine = np.linspace(2.01, 50.0, 5000)
    rs_grid_fine = rstar_schw(r_grid_fine, M)
    def r_of_rs(rs):
        return float(np.interp(rs, rs_grid_fine, r_grid_fine))

    r_vals = np.array([r_of_rs(rs_val) for rs_val in rs])
    dV_vals = np.array([dV_func(r) for r in r_vals])

    # Truncate where |Ψ| becomes too large (divergence)
    # Keep only where |Ψ| < 1e6 * max(|Ψ| near peak)
    Psi_near_peak = np.abs(Psi[len(Psi)//2 - 100:len(Psi)//2 + 100])
    threshold = 1e3 * np.max(Psi_near_peak)
    mask = np.abs(Psi) < threshold

    if np.sum(mask) < 10:
        return None

    rs_trunc = rs[mask]
    Psi_trunc = Psi[mask]
    dV_trunc = dV_vals[mask]

    # Compute integrals
    integrand_num = Psi_trunc**2 * dV_trunc
    integrand_den = Psi_trunc**2

    num = np.trapz(integrand_num, rs_trunc)
    den = np.trapz(integrand_den, rs_trunc)

    delta_omega = num / (2.0 * omega_0 * den)
    return delta_omega


# ---------------------------------------------------------------------------
# Coupled system correction (polar-scalar mixing)
# ---------------------------------------------------------------------------
def compute_coupled_correction(omega_0_grav, omega_0_scalar, eta, M=1.0, l=2):
    """Compute the coupled polar-scalar QNM corrections.

    For the polar gravitational-led mode, the correction has two parts:
    1. The potential correction δV_grav (same as axial)
    2. The mixing with the scalar mode: δω_mix = |<Ψ_grav|V_mix|Ψ_scalar>|² / (ω_grav - ω_scalar)

    For the scalar-led mode:
    1. The potential correction δV_scalar
    2. The mixing with the gravitational mode: δω_mix = |<Ψ_scalar|V_mix|Ψ_grav>|² / (ω_scalar - ω_grav)
    """
    # Compute wave functions
    result_g = compute_qnm_wavefunction(omega_0_grav, lambda r: V_RW_schw(r, M, l), M, l)
    result_s = compute_qnm_wavefunction(omega_0_scalar, lambda r: V_scalar_schw(r, M, l), M, l)

    if result_g is None or result_s is None:
        return None

    rs_g, Psi_g = result_g
    rs_s, Psi_s = result_s

    # Interpolate onto common grid
    rs_common = np.linspace(max(rs_g[0], rs_s[0]), min(rs_g[-1], rs_s[-1]), 3000)
    Psi_g_interp = np.interp(rs_common, rs_g, Psi_g)
    Psi_s_interp = np.interp(rs_common, rs_s, Psi_s)

    # Compute V_mix at grid points
    r_grid_fine = np.linspace(2.01, 50.0, 5000)
    rs_grid_fine = rstar_schw(r_grid_fine, M)
    def r_of_rs(rs):
        return float(np.interp(rs, rs_grid_fine, r_grid_fine))

    r_vals = np.array([r_of_rs(rs_val) for rs_val in rs_common])
    Vm_vals = np.array([V_mix_sgb(r, eta, M, l) for r in r_vals])

    # Truncate
    threshold = 1e3
    mask = (np.abs(Psi_g_interp) < threshold) & (np.abs(Psi_s_interp) < threshold)

    if np.sum(mask) < 10:
        return None

    rs_t = rs_common[mask]
    Pg = Psi_g_interp[mask]
    Ps = Psi_s_interp[mask]
    Vm = Vm_vals[mask]

    # Mixing matrix element: <grav|V_mix|scalar>
    mix_gs = np.trapz(Pg * Vm * Ps, rs_t)
    # Normalization integrals
    norm_g = np.trapz(Pg**2, rs_t)
    norm_s = np.trapz(Ps**2, rs_t)

    # Second-order correction from mixing:
    # δω_grav_mix = |mix_gs|² / (norm_g * (ω_grav - ω_scalar))
    # δω_scalar_mix = |mix_gs|² / (norm_s * (ω_scalar - ω_grav))

    delta_omega = omega_0_grav - omega_0_scalar
    if abs(delta_omega) < 1e-10:
        return None

    dw_grav_mix = abs(mix_gs)**2 / (norm_g * delta_omega)
    dw_scalar_mix = abs(mix_gs)**2 / (norm_s * (-delta_omega))

    return {
        "grav_mix_correction": dw_grav_mix,
        "scalar_mix_correction": dw_scalar_mix,
        "mix_matrix_element": mix_gs,
    }


# ---------------------------------------------------------------------------
# Eikonal formula from Bryant et al. 2021
# ---------------------------------------------------------------------------
def eikonal_qnm(l, alpha_gb, f0_prime=1.0):
    eps = l + 0.5
    sqrt3 = np.sqrt(3.0)
    omega_R_0 = l / (3.0 * sqrt3)
    delta_R = (4.0 / 27.0) * (alpha_gb**2 * l**2 * f0_prime**2 / eps**2) * (1.0 + 3.0 * eps / (2.0 * l))
    omega_R_plus = omega_R_0 * (1.0 + eps / (2.0 * l) + delta_R)
    omega_R_minus = omega_R_0 * (1.0 + eps / (2.0 * l) - delta_R)
    S_m_pp = (l / 27.0) * (1.0 - (560.0 / 2187.0) * (alpha_gb**2 * l * f0_prime**2 / eps))
    delta_I = (44.0 / 729.0) * (alpha_gb**2 * l**2 * f0_prime**2 / eps**2)
    omega_I_0 = -(3.0 * sqrt3 * eps) / (2.0 * l)
    omega_I_plus = omega_I_0 * (1.0 - delta_I) * S_m_pp
    omega_I_minus = omega_I_0 * (1.0 + delta_I) * S_m_pp
    return {
        "omega_R_plus": float(omega_R_plus),
        "omega_R_minus": float(omega_R_minus),
        "omega_I_plus": float(omega_I_plus),
        "omega_I_minus": float(omega_I_minus),
        "splitting_percent": float(abs(omega_R_plus - omega_R_minus) / omega_R_0 * 100),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("TRACK A: HORIZON-BEARING sGB QNM (perturbative)")
    print("Shift-symmetric sGB (Sotiriou-Zhou 2014)")
    print("First-order perturbation theory with exact Schwarzschild QNM base")
    print(f"M={M}, l={L}")
    print("=" * 70)

    # Get exact Schwarzschild QNMs from qnm package
    print("\n--- Exact Schwarzschild QNMs (Leaver) ---")
    mode_grav = qnm.modes_cache(s=-2, l=L, m=L, n=0)
    omega_0_grav, _, _ = mode_grav(a=0.0)
    print(f"  Gravitational l={L}: {omega_0_grav:.8f}")

    # Scalar QNM (s=0, l=2) - test scalar field on Schwarzschild
    mode_scalar = qnm.modes_cache(s=0, l=L, m=0, n=0)
    omega_0_scalar, _, _ = mode_scalar(a=0.0)
    print(f"  Scalar l={L}:       {omega_0_scalar:.8f}")

    all_results = []

    for eta in [0.0] + ETA_VALUES:
        print(f"\n{'='*50}")
        print(f"eta = {eta}")
        print(f"{'='*50}")

        r_h, _ = sgb_F(2.0*M, eta, M)
        alpha_gb = eta / 3.0
        beta = eta / (3.0 * (2.0*M)**2)
        print(f"  alpha_GB = {alpha_gb:.6f}, beta = {beta:.6f}")
        print(f"  r_H = {r_h:.6f}  (shift: {(r_h-2.0*M)/(2.0*M)*100:.4f}%)")

        # 1. Axial (gravitational) correction
        print(f"\n  Axial (gravitational) correction:")
        dw_axial = compute_qnm_correction(
            omega_0_grav,
            lambda r: V_RW_schw(r, M, L),
            lambda r: delta_V_grav(r, eta, M, L),
            M, L
        )
        if dw_axial is not None:
            omega_axial = omega_0_grav + dw_axial
            dR = (omega_axial.real - omega_0_grav.real) / omega_0_grav.real * 100
            dI = (omega_axial.imag - omega_0_grav.imag) / omega_0_grav.imag * 100
            print(f"    δω = {dw_axial:.6e}")
            print(f"    ω_axial = {omega_axial:.6f}")
            print(f"    Shift: Re={dR:+.4f}%  Im={dI:+.4f}%")
        else:
            print(f"    FAILED")
            omega_axial = omega_0_grav
            dR = dI = 0.0

        # 2. Scalar correction
        print(f"\n  Scalar correction:")
        dw_scalar = compute_qnm_correction(
            omega_0_scalar,
            lambda r: V_scalar_schw(r, M, L),
            lambda r: delta_V_scalar(r, eta, M, L),
            M, L
        )
        if dw_scalar is not None:
            omega_scalar = omega_0_scalar + dw_scalar
            dR_s = (omega_scalar.real - omega_0_scalar.real) / omega_0_scalar.real * 100
            dI_s = (omega_scalar.imag - omega_0_scalar.imag) / omega_0_scalar.imag * 100
            print(f"    δω = {dw_scalar:.6e}")
            print(f"    ω_scalar = {omega_scalar:.6f}")
            print(f"    Shift: Re={dR_s:+.4f}%  Im={dI_s:+.4f}%")
        else:
            print(f"    FAILED")
            omega_scalar = omega_0_scalar
            dR_s = dI_s = 0.0

        # 3. Polar-scalar mixing (second-order correction)
        print(f"\n  Polar-scalar mixing:")
        mix_result = compute_coupled_correction(omega_0_grav, omega_0_scalar, eta, M, L)
        if mix_result is not None:
            dw_polar_mix = mix_result["grav_mix_correction"]
            dw_scalar_mix = mix_result["scalar_mix_correction"]
            print(f"    Mix matrix element: {mix_result['mix_matrix_element']:.6e}")
            print(f"    Grav mix correction: {dw_polar_mix:.6e}")
            print(f"    Scalar mix correction: {dw_scalar_mix:.6e}")
        else:
            dw_polar_mix = 0.0
            dw_scalar_mix = 0.0
            print(f"    FAILED")

        # Total polar-led correction = potential correction + mixing correction
        omega_polar = omega_0_grav + dw_axial + dw_polar_mix
        omega_scalar_led = omega_0_scalar + dw_scalar + dw_scalar_mix

        dR_polar = (omega_polar.real - omega_0_grav.real) / omega_0_grav.real * 100
        dI_polar = (omega_polar.imag - omega_0_grav.imag) / omega_0_grav.imag * 100
        dR_scalar_led = (omega_scalar_led.real - omega_0_scalar.real) / omega_0_scalar.real * 100
        dI_scalar_led = (omega_scalar_led.imag - omega_0_scalar.imag) / omega_0_scalar.imag * 100

        print(f"\n  Total corrections:")
        print(f"    Axial:       {omega_axial:.6f}  (Re={dR:+.4f}%, Im={dI:+.4f}%)")
        print(f"    Polar-led:   {omega_polar:.6f}  (Re={dR_polar:+.4f}%, Im={dI_polar:+.4f}%)")
        print(f"    Scalar-led:  {omega_scalar_led:.6f}  (Re={dR_scalar_led:+.4f}%, Im={dI_scalar_led:+.4f}%)")
        print(f"    Iso breaking: Re={dR_polar-dR:+.4f}%  Im={dI_polar-dI:+.4f}%")

        result = {
            "eta": float(eta),
            "r_h": float(r_h),
            "alpha_GB": float(alpha_gb),
            "axial_qnm": {"real": float(omega_axial.real), "imag": float(omega_axial.imag)},
            "scalar_qnm": {"real": float(omega_scalar.real), "imag": float(omega_scalar.imag)},
            "polar_led_qnm": {"real": float(omega_polar.real), "imag": float(omega_polar.imag)},
            "scalar_led_qnm": {"real": float(omega_scalar_led.real), "imag": float(omega_scalar_led.imag)},
            "axial_shift": {"real_percent": float(dR), "imag_percent": float(dI)},
            "polar_led_shift": {"real_percent": float(dR_polar), "imag_percent": float(dI_polar)},
            "scalar_led_shift": {"real_percent": float(dR_scalar_led), "imag_percent": float(dI_scalar_led)},
            "isospectrality_breaking": {
                "real_percent": float(dR_polar - dR),
                "damping_percent": float(dI_polar - dI),
            },
            "corrections": {
                "axial_potential": {"real": float(dw_axial.real), "imag": float(dw_axial.imag)} if dw_axial else None,
                "scalar_potential": {"real": float(dw_scalar.real), "imag": float(dw_scalar.imag)} if dw_scalar else None,
                "polar_mixing": {"real": float(dw_polar_mix.real), "imag": float(dw_polar_mix.imag)} if isinstance(dw_polar_mix, complex) else None,
                "scalar_mixing": {"real": float(dw_scalar_mix.real), "imag": float(dw_scalar_mix.imag)} if isinstance(dw_scalar_mix, complex) else None,
            },
            "eikonal_bryant_2021": eikonal_qnm(L, alpha_gb),
        }
        all_results.append(result)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Horizon-bearing sGB QNM shifts (perturbative)")
    print("=" * 70)
    print(f"\n{'eta':>6s}  {'alpha':>8s}  {'Axial Re%':>10s}  {'Polar Re%':>10s}  {'Scalar Re%':>10s}  {'Iso Re%':>9s}  {'Eik%':>8s}")
    print("-" * 75)
    for r in all_results:
        eta = r["eta"]
        alpha = r["alpha_GB"]
        ar = r["axial_shift"]["real_percent"]
        pr = r["polar_led_shift"]["real_percent"]
        sr = r["scalar_led_shift"]["real_percent"]
        ir = r["isospectrality_breaking"]["real_percent"]
        es = r["eikonal_bryant_2021"]["splitting_percent"]
        print(f"{eta:6.3f}  {alpha:8.5f}  {ar:+10.4f}  {pr:+10.4f}  {sr:+10.4f}  {ir:+9.4f}  {es:8.4f}")

    # O(alpha^2) scaling
    print("\nO(alpha^2) scaling:")
    for r in all_results[1:]:
        alpha2 = r["alpha_GB"]**2
        ar = r["axial_shift"]["real_percent"]
        pr = r["polar_led_shift"]["real_percent"]
        if alpha2 > 0:
            print(f"  eta={r['eta']:.3f}: axial/alpha²={ar/alpha2:.2f}, polar/alpha²={pr/alpha2:.2f}")

    # Comparison to published data
    print("\n--- Comparison to published results ---")
    print("  Bryant et al. 2021 (eikonal, theory-universal):")
    print("    Axial O(alpha^2) shift, symmetric Zeeman splitting")
    print("  Blazquez-Salcedo et al. 2016 (EdGB, dilatonic coupling):")
    print("    Axial R_2 = +1.002e-3, Polar grav R_2 = -3.135e-2 (opposite signs)")
    print("    Scalar-led R_1 = -1.408e-2 (LINEAR in coupling)")
    print()
    print("  Our shift-symmetric sGB results should show:")
    print("    - O(alpha^2) for axial (no O(alpha) term)")
    print("    - O(alpha^2) for polar gravitational-led (from potential correction)")
    print("    - O(alpha) for scalar-led mixing correction (from V_mix ~ alpha)")
    print("    - Isospectrality breaking between axial and polar")

    output = {
        "method": "First-order perturbation theory with exact Schwarzschild QNM (qnm/Leaver) + sGB potential correction + second-order polar-scalar mixing",
        "parameters": {"M": M, "l": L},
        "schwarzschild_reference": {
            "gravitational": {"real": float(omega_0_grav.real), "imag": float(omega_0_grav.imag)},
            "scalar": {"real": float(omega_0_scalar.real), "imag": float(omega_0_scalar.imag)},
        },
        "eta_values": [0.0] + ETA_VALUES,
        "results": all_results,
        "references": {
            "sotiriou_zhou_2014": "Sotiriou & Zhou, PRD 90, 124063 (2014)",
            "bryant_2021": "Bryant et al., PRD 104, 044051 (2021)",
            "chung_yunes_2024": "Chung & Yunes, PRD 110, 064019 (2024)",
            "blazquez_2016": "Blazquez-Salcedo et al., PRD 94, 104024 (2016)",
            "leaver_1985": "Leaver, Proc. R. Soc. A 402, 285 (1985)",
        },
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "step_27_qnm_horizon_branch.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
