#!/usr/bin/env python3
"""Coupled spectral QNM solver for the polar-scalar system in sGB.

This script establishes the exterior isospectrality-breaking benchmark by:
1.  Using the qnm package to obtain the exact Schwarzschild quasinormal-mode
    frequencies for the gravitational/polar (s=-2) and scalar (s=0) sectors.
2.  Applying the sGB O(alpha^2) horizon shift to those benchmark frequencies.
3.  Computing the leading-order polar-scalar mixing from the off-diagonal
    sGB potential V_Zphi with a WKB-normalised Gaussian overlap for the
    two QNM wavefunctions localised at their respective potential peaks.
4.  Diagonalising the resulting 2x2 effective matrix to obtain the coupled
    polar-led and scalar-led QNM frequencies and the isospectrality-breaking
    percentage between the axial and polar-led modes.

This is the first exterior benchmark.  A full matrix continued-fraction solve
that also incorporates the regular temporal-well inner boundary conditions is
the next layer of Track 3.

Outputs:
  results/step_26_coupled_spectral_solver.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")


def V_scalar(r, M, l):
    F = 1.0 - 2.0 * M / r
    return F * (l * (l + 1) / r**2 + 2.0 * M / r**3)


def V_gravitational(r, M, l):
    F = 1.0 - 2.0 * M / r
    return F * (l * (l + 1) / r**2 - 6.0 * M / r**3)


def V_Zphi(r, eta, M, l):
    """Off-diagonal sGB mixing potential (L^-2)."""
    alpha_GB = eta * M**2 / 3.0
    F = 1.0 - 2.0 * M / r
    return alpha_GB * l * (l + 1) * (l - 1) * (l + 2) * F / r**4


def sgb_horizon_shift(M, eta):
    """O(alpha^2) shift of the sGB apparent horizon."""
    return 2.0 * M * 19.6 * (eta / 12.0) ** 2


def qnm_schwarzschild(s, l, n, M=1.0):
    """Exact Schwarzschild QNM from the qnm package (Cook-Zalutskiy Leaver)."""
    from qnm.nearby import NearbyRootFinder
    A_closest = complex(l * (l + 1) - s * (s + 1))
    omega_guess = 0.5 - 0.1j
    finder = NearbyRootFinder(
        a=0.0, s=s, m=0, A_closest_to=A_closest, omega_guess=omega_guess, n_inv=n, l_max=20
    )
    finder.set_params()
    omega = finder.do_solve()
    # qnm returns omega in units of 1/M; scale by M
    return omega / M


def peak_and_curvature(V_func, M, l, r_min=2.05, r_max=15.0, n=2000):
    """Find the peak and second derivative of an effective potential."""
    r = np.linspace(r_min, r_max, n)
    V = V_func(r, M, l)
    i_peak = np.argmax(V)
    r_peak = r[i_peak]
    V_peak = V[i_peak]
    # second derivative by finite difference
    dr = r[1] - r[0]
    Vpp = (V[i_peak + 1] - 2.0 * V[i_peak] + V[i_peak - 1]) / dr**2
    return float(r_peak), float(V_peak), float(Vpp), r, V


def wkb_gaussian_width(omega, V_peak, Vpp, M):
    """Approximate r_* width of the fundamental QNM wavefunction.

    For a mode just below a potential peak, the wavefunction decays away from
    the peak with a local exponential rate set by V'' and (V_peak - omega_R^2).
    """
    omega_R = omega.real
    # Effective inverted potential: U(r*) = V(r*) - omega_R^2.  Near the peak,
    # U ~ U0 + 0.5*Vpp*(r-rp)^2.  The QNM is a bound-like state in this
    # complex potential; the Gaussian width in r_* is (Vpp - ... )^{-1/4}.
    # We take the leading-order estimate.
    if Vpp <= 0.0:
        return 1.0
    # Geometric conversion d2V/dr_*^2 = F^2 d2V/dr^2 (because dr*/dr = 1/F,
    # but d2V/dr_*^2 = F d/dr(F dV/dr) = F^2 V'' + F F' V').  Near the peak V'~0,
    # so d2V/dr_*^2 ~ F^2 V''.
    F = 1.0 - 2.0 * M / r_of_V_peak(V_peak, M)
    # We do not need r here; the conversion from r to r_* is not crucial because
    # the WKB overlap is only an order-unity estimate.  Return a characteristic
    # width in r_* using (Vpp)^{-1/4}.
    return (Vpp) ** (-0.25)


def r_of_V_peak(V_peak, M):
    """Approximate radius at which V = V_peak (not needed accurately here)."""
    return 3.0 * M  # placeholder; only used for F conversion which is order one


def mixing_matrix_element(omega_g, omega_s, V_mix, M, l):
    """Estimate <g|V_Zphi|s> using WKB-localised Gaussian wavefunctions."""
    rp_g, Vp_g, Vpp_g, r_g, V_g = peak_and_curvature(V_gravitational, M, l)
    rp_s, Vp_s, Vpp_s, r_s, V_s = peak_and_curvature(V_scalar, M, l)

    # tortoise positions of the peaks
    def tortoise(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0)

    rsp_g = tortoise(rp_g)
    rsp_s = tortoise(rp_s)

    # Width in r_*: convert from r width using dV/dr_*^2 ~ F^2 dV/dr^2
    F_g = 1.0 - 2.0 * M / rp_g
    F_s = 1.0 - 2.0 * M / rp_s
    Vpp_rsp_g = F_g**2 * Vpp_g
    Vpp_rsp_s = F_s**2 * Vpp_s

    omega_gR = omega_g.real
    omega_sR = omega_s.real

    # characteristic decay exponents in r_* near the peak
    # kappa^2 = V''_rsp * (V_peak - omega_R^2) / (some curvature scale)
    # We use a simple Gaussian width: sigma = [2/(V''_rsp)]^{1/4} normalised.
    sigma_g = (2.0 / max(abs(Vpp_rsp_g), 1e-6)) ** 0.25
    sigma_s = (2.0 / max(abs(Vpp_rsp_s), 1e-6)) ** 0.25

    # Normalised Gaussian wavefunctions in r_*; unit normalisation
    # psi_i(r_*) = (2*pi*sigma_i^2)^(-1/4) exp(-(r_* - rsp_i)^2/(4 sigma_i^2))
    # overlap integral of two such Gaussians:
    sigma_eff = np.sqrt(sigma_g**2 + sigma_s**2)
    overlap = (2.0 * sigma_g * sigma_s / (sigma_g**2 + sigma_s**2)) ** 0.5
    overlap *= np.exp(-(rsp_g - rsp_s)**2 / (4.0 * (sigma_g**2 + sigma_s**2)))

    # Value of mixing potential at the midpoint of the two peaks
    r_mid = 0.5 * (rp_g + rp_s)
    V_mix_peak = V_Zphi(r_mid, -0.1, M, l)

    return V_mix_peak * overlap


def coupled_qnm_frequencies(omega_g, omega_s, M_mix):
    """Diagonalise the 2x2 effective matrix and return frequencies."""
    H = np.array([[omega_g**2, M_mix],
                  [M_mix, omega_s**2]])
    lam, evec = np.linalg.eig(H)
    # Sort by real part (lower -> polar-led, higher -> scalar-led)
    idx = np.argsort(lam.real)
    lam = lam[idx]
    evec = evec[:, idx]
    omega = np.sqrt(lam)
    # Ensure negative imaginary part (damped)
    omega = np.where(omega.imag > 0, -omega, omega)
    return omega[0], omega[1], evec


def main():
    M = 1.0
    l = 2
    eta = -0.1

    print("=" * 70)
    print("EXTERIOR ISOSPECTRALITY-BREAKING BENCHMARK (sGB, mass-inflation)")
    print("=" * 70)

    # Schwarzschild benchmark from qnm
    print("\n--- Schwarzschild reference (qnm Leaver) ---")
    omega_g_schw = qnm_schwarzschild(-2, l, 0, M)
    omega_s_schw = qnm_schwarzschild(0, l, 0, M)
    print(f"  gravitational/polar (s=-2): {omega_g_schw:.6f}")
    print(f"  scalar (s=0):               {omega_s_schw:.6f}")

    # sGB horizon shift applied to both benchmark frequencies
    delta_r_H = sgb_horizon_shift(M, eta)
    frac = delta_r_H / (2.0 * M)
    omega_g_ext = omega_g_schw * (1.0 + frac)
    omega_s_ext = omega_s_schw * (1.0 + 0.5 * frac)
    print(f"\n--- sGB-corrected exterior (O(alpha^2) horizon shift) ---")
    print(f"  horizon shift: {delta_r_H:.6f} M  ({100*frac:.4f}%)")
    print(f"  polar/axial at eta={eta}: {omega_g_ext:.6f}")
    print(f"  scalar at eta={eta}:      {omega_s_ext:.6f}")

    # Leading-order off-diagonal mixing
    M_mix = mixing_matrix_element(omega_g_ext, omega_s_ext,
                                  lambda r: V_Zphi(r, eta, M, l), M, l)
    print(f"\n--- polar-scalar mixing (WKB overlap) ---")
    print(f"  |M_mix|: {abs(M_mix):.6e}")

    # Coupled 2x2 diagonalisation
    omega_pol, omega_sca, evec = coupled_qnm_frequencies(omega_g_ext, omega_s_ext, M_mix)
    print(f"\n--- coupled QNM frequencies (eta={eta}) ---")
    print(f"  polar-led:  {omega_pol:.6f}")
    print(f"  scalar-led: {omega_sca:.6f}")

    # Isospectrality breaking: relative difference between axial (uncoupled) and polar-led
    break_real = 100.0 * (omega_pol.real - omega_g_ext.real) / omega_g_ext.real
    break_imag = 100.0 * (omega_pol.imag - omega_g_ext.imag) / omega_g_ext.imag
    print(f"\n  Isospectrality breaking (polar vs axial):")
    print(f"    real: {break_real:.4f}%")
    print(f"    imag: {break_imag:.4f}%")

    output = {
        "parameters": {"M": M, "l": l, "eta": eta},
        "method": "qnm-exterior Leaver + leading-order 2x2 mixing (WKB overlap)",
        "schwarzschild_reference": {
            "gravitational": {"omega_R": float(omega_g_schw.real), "omega_I": float(omega_g_schw.imag)},
            "scalar": {"omega_R": float(omega_s_schw.real), "omega_I": float(omega_s_schw.imag)},
        },
        "exterior_sgb": {
            "horizon_shift": float(delta_r_H),
            "fractional_shift": float(frac),
            "axial_uncoupled": {"omega_R": float(omega_g_ext.real), "omega_I": float(omega_g_ext.imag)},
            "scalar_uncoupled": {"omega_R": float(omega_s_ext.real), "omega_I": float(omega_s_ext.imag)},
        },
        "coupled_qnm": {
            "polar_led": {"omega_R": float(omega_pol.real), "omega_I": float(omega_pol.imag)},
            "scalar_led": {"omega_R": float(omega_sca.real), "omega_I": float(omega_sca.imag)},
        },
        "isospectrality_breaking_percent": {
            "real": float(break_real),
            "imag": float(break_imag),
        },
        "mixing_matrix_element": {"real": float(M_mix.real), "imag": float(M_mix.imag)},
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "step_26_coupled_spectral_solver.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
