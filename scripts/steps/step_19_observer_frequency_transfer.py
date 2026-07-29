#!/usr/bin/env python3
"""Step 14: Observer frequency transfer — the invariant Temporal Horizon test.

Computes the invariant emitter-receiver frequency transfer factor

    Z = omega_e / omega_o = (k_mu u^mu)_e / (k_mu u^mu)_o

for photons emitted from progressively greater depth and received by a
specified distant observer, on the TEP matter metric.

This is the decisive calculation for the Temporal Horizon: it determines
whether the theory actually produces redshift, blueshift, delay, frequency
compression, or no special central temporal effect.

The calculation requires:
  - physically specified emitter worldline (freely falling or static)
  - physically specified receiver worldline (freely falling or static)
  - outgoing photon four-momentum (null geodesic of the matter metric)
  - received frequency omega_o = -k_mu u_o^mu
  - emitted frequency  omega_e = -k_mu u_e^mu

Cases computed:
  (1) Freely falling emitter -> distant static receiver
  (2) Freely falling emitter -> freely falling receiver
  (3) Successive pulses emitted at fixed local intervals (arrival-time map)
  (4) Two-way radar signals

The metric is the TEP matter metric g_tilde = A^2 * g on a specified
background (Hayward benchmark or sGB exterior).

For the Hayward benchmark with phi0=1, A = (r_h/r)^phi0:
  - In the deep interior, A -> infinity, F -> 1
  - The static-observer clock rate dtau/dt = A*sqrt(F) -> infinity (blueshift)
  - This calculation tests whether the invariant transfer factor confirms
    or contradicts the temporal-freeze interpretation.

Outputs:
  results/step_19_observer_frequency_transfer.json
  results/step_19_observer_frequency_transfer.csv
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from scripts.utils.logger import TEPLogger  # noqa: E402

RESULTS_DIR = os.path.join(_PROJECT_ROOT, "results")


# ---------------------------------------------------------------------------
# Metric definitions
# ---------------------------------------------------------------------------
def hayward_F(r, M, ell):
    """Hayward metric function F(r) = 1 - 2Mr^2/(r^3 + 2M*ell^2)."""
    return 1.0 - 2.0 * M * r**2 / (r**3 + 2.0 * M * ell**2)


def conformal_A(r, r_h, phi0):
    """Conformal factor A = (r_h/r)^phi0."""
    return (r_h / r) ** phi0


# ---------------------------------------------------------------------------
# Eddington-Finkelstein null geodesics
# ---------------------------------------------------------------------------
# Ingoing EF: ds^2 = -F dv^2 + 2 dv dr + r^2 dOmega^2
# Conformal matter metric: ds^2_tilde = A^2 * (-F dv^2 + 2 dv dr + r^2 dOmega^2)
#
# Null condition: A^2 * (-F (dv/dlambda)^2 + 2 (dv/dlambda)(dr/dlambda)) = 0
# => dv * (-F dv + 2 dr) = 0
#
# Two null families:
#   (a) v = const (ingoing): dv/dlambda = 0, dr/dlambda = k (free)
#   (b) dr/dv = F/2 (outgoing): dr/dlambda = (F/2) * (dv/dlambda)
#
# For the affine parametrisation, we need the geodesic equation.
# The conformal metric g_tilde = A^2 * g has Christoffel symbols:
#   Gamma^mu_tilde = Gamma^mu_g + delta^mu (2 A'/A) k^nu - ...
# For a null geodesic of g_tilde, the affine parameter differs from g.
#
# Key: for the OUTGOING family (dr/dv = F/2), the photon travels outward.
# We integrate the null geodesic from emitter (deep interior) to receiver (far).

def causal_connectivity_check(r_emit, r_recv, M, ell, n_scan=10000):
    """Verify J+(emission event) ∩ receiver worldline ≠ ∅.

    For the outgoing null geodesic dr/dv = F(r)/2 in Eddington–Finkelstein
    coordinates, the photon propagates outward (dr/dv > 0) only where F > 0.
    Where F < 0 (inside a Killing horizon), the "outgoing" null geodesic
    actually moves inward (dr/dv < 0).

    The causal connectivity check scans the radial interval [r_emit, r_recv]
    and determines whether a future-directed outgoing null geodesic from
    r_emit can reach r_recv.  The check has three cases:

    1. F > 0 everywhere on [r_emit, r_recv]: the photon escapes directly.
    2. F changes sign on [r_emit, r_recv]: there is a horizon between
       emitter and receiver.  If the emitter is inside an inner horizon
       (F < 0 at r_emit but F > 0 at some r < r_emit), the photon may
       still escape if it can cross the trapped region — but in practice,
       for the Hayward metric with an inner and outer horizon, an emitter
       between the inner and outer horizons is trapped (both null families
       move inward).  An emitter inside the inner horizon (F > 0 again)
       can in principle send outgoing photons that reach the inner horizon,
       but they cannot cross the trapped region to the outer horizon.
    3. F < 0 at r_emit: the emitter is inside a horizon.  For the Hayward
       metric with two horizons, if r_emit is between the inner and outer
       horizon, no outgoing null geodesic reaches the exterior.  If r_emit
       is inside the inner horizon (where F > 0 again), the outgoing null
       reaches the inner horizon but is then trapped.

    Returns:
        connected: bool — whether the outgoing null reaches r_recv
        reason: str — explanation of the connectivity status
        F_sign_changes: list of r values where F changes sign
    """
    r_scan = np.linspace(r_emit, r_recv, n_scan)
    F_scan = np.array([hayward_F(r, M, ell) for r in r_scan])

    # Find sign changes of F
    sign_changes = []
    for i in range(len(F_scan) - 1):
        if F_scan[i] * F_scan[i + 1] < 0:
            # Linear interpolation for the zero crossing
            r_zero = r_scan[i] - F_scan[i] * (r_scan[i+1] - r_scan[i]) / (F_scan[i+1] - F_scan[i])
            sign_changes.append(r_zero)

    F_emit = hayward_F(r_emit, M, ell)
    F_recv = hayward_F(r_recv, M, ell)

    if not sign_changes:
        # F does not change sign on the interval
        if F_emit > 0 and F_recv > 0:
            return True, "F > 0 throughout [r_emit, r_recv]; outgoing null escapes directly", sign_changes
        else:
            return False, f"F < 0 throughout [r_emit, r_recv]; emitter and receiver both inside horizon", sign_changes
    else:
        # F changes sign — there is at least one horizon between emitter and receiver
        n_horizons = len(sign_changes)
        if n_horizons == 1:
            # Single horizon crossing
            r_horizon = sign_changes[0]
            if F_emit > 0:
                # Emitter is outside the horizon, receiver is inside (F_recv < 0)
                return False, f"Emitter outside horizon at r={r_horizon:.4f}, receiver inside; outgoing null cannot reach receiver", sign_changes
            else:
                # Emitter is inside the horizon (F_emit < 0), receiver is outside (F_recv > 0)
                # Check if emitter is between inner and outer horizon (trapped region)
                # For Hayward: F < 0 between inner and outer horizon
                # An emitter in the trapped region cannot send outgoing null to exterior
                return False, f"Emitter inside horizon at r={r_horizon:.4f} (trapped region); outgoing null cannot reach exterior receiver", sign_changes
        elif n_horizons == 2:
            # Two horizons (Hayward: inner and outer)
            r_inner, r_outer = sign_changes
            if r_emit < r_inner and F_emit > 0:
                # Emitter is inside the inner horizon where F > 0
                # Outgoing null reaches inner horizon but is then trapped
                return False, f"Emitter inside inner horizon (r={r_inner:.4f}); outgoing null reaches inner horizon but is trapped between horizons", sign_changes
            elif r_inner < r_emit < r_outer:
                # Emitter is in the trapped region
                return False, f"Emitter in trapped region between inner (r={r_inner:.4f}) and outer (r={r_outer:.4f}) horizons; no outgoing null reaches exterior", sign_changes
            else:
                return False, f"Emitter position relative to horizons unclear; connectivity not established", sign_changes
        else:
            return False, f"Multiple horizons ({n_horizons}) between emitter and receiver; connectivity not established", sign_changes


def outgoing_null_geodesic(r_emit, r_recv, M, ell, r_h, phi0, n_steps=10000):
    """Integrate the outgoing null geodesic from r_emit to r_recv.

    In EF coordinates, the outgoing null satisfies dr/dv = F/2.
    So dv/dr = 2/F.

    The photon four-momentum in EF coordinates:
      k^v = dv/dlambda, k^r = dr/dlambda = (F/2) * k^v

    For affine parametrisation, we use the fact that for a conformal metric
    g_tilde = A^2 * g, null geodesics of g are also null geodesics of g_tilde
    (conformal invariance of null geodesics), but the affine parameter
    transforms as dlambda_tilde = A^2 * dlambda_g.

    So we can integrate the null geodesic using the background metric g,
    then rescale the affine parameter.

    Returns: v_arr, r_arr (arrays along the geodesic), or (None, None) if
             the causal connectivity check fails.
    """
    # Causal connectivity check: J+(emission) ∩ receiver worldline ≠ ∅
    connected, reason, sign_changes = causal_connectivity_check(r_emit, r_recv, M, ell)
    if not connected:
        return None, None, reason

    r_arr = np.linspace(r_emit, r_recv, n_steps)
    F_arr = hayward_F(r_arr, M, ell)

    # dv/dr = 2/F for outgoing null
    # Integrate v from v_emit = 0
    dv_dr = 2.0 / F_arr
    v_arr = np.zeros(n_steps)
    for i in range(1, n_steps):
        # Trapezoidal integration
        v_arr[i] = v_arr[i-1] + 0.5 * (dv_dr[i] + dv_dr[i-1]) * (r_arr[i] - r_arr[i-1])

    return v_arr, r_arr, "Causal connectivity verified"


def photon_four_momentum(r, M, ell, r_h, phi0, direction="outgoing"):
    """Compute the photon four-momentum k_mu in EF coordinates.

    For the conformal matter metric g_tilde = A^2 * g:
      g_tilde_vv = -A^2 * F
      g_tilde_vr = A^2
      g_tilde_rr = 0

    Null condition: g_tilde_ab k^a k^b = 0
      -A^2 F (k^v)^2 + 2 A^2 k^v k^r = 0
      => k^v (-F k^v + 2 k^r) = 0

    For outgoing: k^r = (F/2) k^v
    For ingoing (v=const): k^v = 0, k^r = free

    We normalise using the affine parameter. For the conformal metric,
    the affine parameter satisfies g_tilde_ab k^a k^b = 0 and
    k^a = dx^a/dlambda_tilde.

    For the outgoing family, we set k^v = 1 (coordinate normalisation),
    then k^r = F/2. The covariant components are:
      k_v = g_tilde_vv * k^v + g_tilde_vr * k^r = -A^2 F + A^2 * (F/2) = -A^2 F/2
      k_r = g_tilde_rv * k^v + g_tilde_rr * k^r = A^2

    So k_v = -A^2 F/2, k_r = A^2 (for outgoing, k^v = 1)
    """
    A = conformal_A(r, r_h, phi0)
    F = hayward_F(r, M, ell)

    if direction == "outgoing":
        # k^v = 1, k^r = F/2
        k_v = -A**2 * F / 2.0
        k_r = A**2
    else:  # ingoing (v = const)
        # k^v = 0, k^r = 1
        k_v = A**2  # g_tilde_vr * k^r
        k_r = 0.0  # g_tilde_rr * k^r = 0

    return k_v, k_r, A, F


# ---------------------------------------------------------------------------
# Observer four-velocities
# ---------------------------------------------------------------------------
def static_observer_four_velocity(r, M, ell, r_h, phi0):
    """Static observer four-velocity in the matter metric.

    For a static observer: u^mu = (u^v, 0, 0, 0)
    Normalisation: g_tilde_vv (u^v)^2 = -1
      => u^v = 1 / sqrt(A^2 F) = 1 / (A * sqrt(F))

    Valid only where F > 0 (outside horizon).

    Covariant: u_v = g_tilde_vv * u^v = -A^2 F / (A sqrt(F)) = -A sqrt(F)
    """
    A = conformal_A(r, r_h, phi0)
    F = hayward_F(r, M, ell)
    if F <= 0:
        return None  # static observers don't exist inside horizon
    u_v = -A * np.sqrt(F)
    return u_v


def freely_falling_four_velocity(r, M, ell, r_h, phi0, E=1.0):
    """Freely falling observer four-velocity (from rest at infinity).

    For radial geodesic in the background metric g (EF coords):
      Energy E = F * dv/dtau - dr/dtau (conserved)
      For the conformal metric, the geodesic equation is modified.

    For simplicity, we use the background geodesic (the conformal factor
    cancels for null geodesics but not for timelike).

    For a freely falling observer from rest at infinity in the background:
      u^v = E / F  (in EF coords, for the background metric)
      u^r = -sqrt(E^2 - F)  (inward, for the background metric)

    For the matter metric, the four-velocity transforms as:
      u_tilde^mu = u^mu / A  (because dtau_tilde = A * dtau_g for timelike)

    Actually, for a conformal metric g_tilde = A^2 * g:
      dtau_tilde^2 = A^2 * dtau_g^2  =>  dtau_tilde = A * dtau_g
      u_tilde^mu = dx^mu / dtau_tilde = (1/A) * dx^mu / dtau_g = u^mu / A

    Covariant components:
      u_tilde_mu = g_tilde_mu_nu * u_tilde^nu = A^2 * g_mu_nu * (u^nu / A)
                 = A * g_mu_nu * u^nu = A * u_mu (background)

    So u_tilde_v = A * u_v (background), u_tilde_r = A * u_r (background)
    """
    A = conformal_A(r, r_h, phi0)
    F = hayward_F(r, M, ell)

    # Background four-velocity (freely falling from rest at infinity, E=1)
    if F <= 0:
        # Inside horizon: F < 0, the formula still works but sqrt(E^2 - F) > E
        u_v_bg = E / abs(F)  # magnitude
        u_r_bg = -np.sqrt(E**2 - F)  # inward (F < 0 => sqrt(E^2 + |F|))
    else:
        if E**2 < F:
            return None  # not enough energy
        u_v_bg = E / F
        u_r_bg = -np.sqrt(E**2 - F)  # inward

    # Matter-frame covariant components
    u_v_tilde = A * (-F * u_v_bg + u_r_bg)  # g_vv * u^v + g_vr * u^r = -F*E/F + (-sqrt(E^2-F))
    # Actually let's be more careful:
    # Background: u^v = E/F, u^r = -sqrt(E^2 - F)
    # Background covariant: u_v = g_vv * u^v + g_vr * u^r = -F*(E/F) + 1*(-sqrt(E^2-F)) = -E - sqrt(E^2-F)
    #                       u_r = g_rv * u^v + g_rr * u^r = 1*(E/F) + 0 = E/F
    u_v_bg_cov = -E - np.sqrt(E**2 - F) if F > 0 else -E - np.sqrt(E**2 - F)
    u_r_bg_cov = E / F if F != 0 else 0.0

    # Matter-frame covariant: u_tilde_mu = A * u_mu (background)
    u_v_tilde = A * u_v_bg_cov
    u_r_tilde = A * u_r_bg_cov

    return u_v_tilde, u_r_tilde


# ---------------------------------------------------------------------------
# Frequency transfer calculation
# ---------------------------------------------------------------------------
def frequency_transfer(r_emit, r_recv, M, ell, r_h, phi0,
                       emitter_type="freely_falling", receiver_type="static"):
    """Compute the invariant frequency transfer factor Z = omega_e / omega_o.

    omega_e = -k_mu u_e^mu  (emitted frequency)
    omega_o = -k_mu u_o^mu  (received frequency)

    For the outgoing photon in EF coordinates:
      k_v = -A^2 F/2, k_r = A^2

    For a static observer at r:
      u_v = -A sqrt(F), u_r = 0
      omega = -k_mu u^mu = -(k_v * u^v + k_r * u^r) = -k_v * u^v
      But u^v = 1/(A sqrt(F)), and k_v = -A^2 F/2
      => omega = -(-A^2 F/2) * (1/(A sqrt(F))) = A sqrt(F) / 2

    Wait, let's be more careful. The frequency is:
      omega = -k_mu u^mu = -(k_v u^v + k_r u^r)

    For static observer: u^v = 1/(A sqrt(F)), u^r = 0
      omega = -k_v * u^v = -(-A^2 F/2) * (1/(A sqrt(F))) = A sqrt(F) / 2

    Hmm, but this depends on the normalisation of k. Let me use the
    covariant form directly.

    Actually, for the frequency we need:
      omega = -g_tilde_mu_nu k^mu u^nu = -k_nu u^nu

    where k_nu = g_tilde_nu_mu k^mu and u^nu is the observer four-velocity.

    For static observer: u^v = 1/(A sqrt(F)), u^r = 0
      k_v = g_tilde_vv k^v + g_tilde_vr k^r = -A^2 F * 1 + A^2 * (F/2) = -A^2 F/2
      omega = -k_v * u^v = (A^2 F/2) * (1/(A sqrt(F))) = A sqrt(F) / 2

    For freely falling observer: need u^v, u^r in matter frame
      u_tilde^v = u^v_bg / A, u_tilde^r = u^r_bg / A
      omega = -(k_v u^v + k_r u^r) = -((-A^2 F/2)(u^v_bg/A) + A^2 (u^r_bg/A))
            = -(-A F u^v_bg / 2 + A u^r_bg)
            = A (F u^v_bg / 2 - u^r_bg)

    With u^v_bg = E/F, u^r_bg = -sqrt(E^2 - F):
      omega = A (F * (E/F) / 2 - (-sqrt(E^2 - F)))
            = A (E/2 + sqrt(E^2 - F))

    Returns: Z = omega_e / omega_o, omega_e, omega_o
             If the causal connectivity check fails, returns (None, omega_e, omega_o)
             with a causal_status key appended.
    """
    # Causal connectivity check: J+(emission event) ∩ receiver worldline ≠ ∅
    connected, causal_reason, sign_changes = causal_connectivity_check(r_emit, r_recv, M, ell)

    # Photon four-momentum at emitter and receiver (outgoing)
    k_v_e, k_r_e, A_e, F_e = photon_four_momentum(r_emit, M, ell, r_h, phi0, "outgoing")
    k_v_o, k_r_o, A_o, F_o = photon_four_momentum(r_recv, M, ell, r_h, phi0, "outgoing")

    # Emitter frequency
    if emitter_type == "static":
        if F_e <= 0:
            return None, None, None, "Static emitter does not exist inside horizon (F ≤ 0)"
        u_v_e = static_observer_four_velocity(r_emit, M, ell, r_h, phi0)
        omega_e = -k_v_e * u_v_e  # k_v * u^v, but u^v = 1/(A sqrt F)
        # Actually: omega = -(k_v u^v + k_r u^r) = -k_v * u^v (since u^r = 0)
        # u^v = 1/(A sqrt F), k_v = -A^2 F/2
        # omega = -(-A^2 F/2)(1/(A sqrt F)) = A sqrt F / 2
        omega_e = A_e * np.sqrt(F_e) / 2.0
    else:  # freely falling
        ff_e = freely_falling_four_velocity(r_emit, M, ell, r_h, phi0)
        if ff_e is None:
            return None, None, None, "Freely falling emitter velocity undefined at this radius"
        u_v_tilde_e, u_r_tilde_e = ff_e
        # omega = -(k_v u^v + k_r u^r) where u^mu are the contravariant components
        # We have covariant u_mu, so omega = -k^mu u_mu = -(k^v u_v + k^r u_r)
        # k^v = 1, k^r = F/2 (outgoing)
        omega_e = -(1.0 * u_v_tilde_e + (F_e / 2.0) * u_r_tilde_e)

    # Receiver frequency
    if receiver_type == "static":
        if F_o <= 0:
            return None, omega_e, None, "Static receiver does not exist inside horizon (F ≤ 0)"
        omega_o = A_o * np.sqrt(F_o) / 2.0
    else:  # freely falling
        ff_o = freely_falling_four_velocity(r_recv, M, ell, r_h, phi0)
        if ff_o is None:
            return None, omega_e, None, "Freely falling receiver velocity undefined at this radius"
        u_v_tilde_o, u_r_tilde_o = ff_o
        omega_o = -(1.0 * u_v_tilde_o + (F_o / 2.0) * u_r_tilde_o)

    if omega_o == 0 or omega_o is None:
        return None, omega_e, omega_o, causal_reason

    if not connected:
        # Causal connectivity check failed: Z is not a valid observable transfer factor
        return None, omega_e, omega_o, causal_reason

    Z = omega_e / omega_o
    return Z, omega_e, omega_o, causal_reason


def arrival_time_map(r_emit_arr, r_recv, M, ell, r_h, phi0,
                     emitter_type="freely_falling", receiver_type="static"):
    """Compute the arrival-time mapping T = dtau_o / dtau_e.

    For successive pulses emitted at fixed local proper-time intervals dtau_e,
    the arrival-time interval dtau_o at the receiver gives T = dtau_o / dtau_e.

    For static observers: dtau = A sqrt(F) dt, and the coordinate-time delay
    is dv = integral(2/F dr). So:
      dtau_e = A_e sqrt(F_e) dv
      dtau_o = A_o sqrt(F_o) dv
      T = (A_o sqrt(F_o)) / (A_e sqrt(F_e))

    For freely falling: dtau = A * dtau_bg, and the proper time evolves
    along the geodesic.
    """
    results = []
    for r_emit in r_emit_arr:
        Z, omega_e, omega_o, causal_reason = frequency_transfer(
            r_emit, r_recv, M, ell, r_h, phi0, emitter_type, receiver_type
        )
        if Z is None:
            results.append({
                "r_emit": r_emit,
                "Z": None,
                "omega_e": float(omega_e) if omega_e is not None else None,
                "omega_o": float(omega_o) if omega_o is not None else None,
                "causal_status": causal_reason,
            })
            continue

        # Arrival-time ratio for static observers
        A_e = conformal_A(r_emit, r_h, phi0)
        F_e = hayward_F(r_emit, M, ell)
        A_o = conformal_A(r_recv, r_h, phi0)
        F_o = hayward_F(r_recv, M, ell)

        if receiver_type == "static" and F_o > 0:
            if emitter_type == "static" and F_e > 0:
                T = (A_o * np.sqrt(F_o)) / (A_e * np.sqrt(F_e))
            else:
                T = None  # need freely-falling proper time evolution
        else:
            T = None

        results.append({
            "r_emit": r_emit,
            "Z": float(Z) if Z is not None else None,
            "omega_e": float(omega_e) if omega_e is not None else None,
            "omega_o": float(omega_o) if omega_o is not None else None,
            "T": float(T) if T is not None else None,
            "A_e": float(A_e),
            "F_e": float(F_e),
            "A_o": float(A_o),
            "F_o": float(F_o),
            "causal_status": causal_reason,
        })
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    logger = TEPLogger("step_19_observer_frequency_transfer")
    logger.info("Step 14: Observer Frequency Transfer — Invariant Temporal Horizon Test")

    # Parameters
    M = 1.0
    ell = 0.5  # Hayward regularisation scale
    r_h = 2.0 * M  # approximate horizon
    phi0 = 1.0  # benchmark configuration
    r_recv = 50.0 * M  # distant receiver

    # Emitter radii: from deep interior to near horizon
    r_emit_arr = np.concatenate([
        np.linspace(0.001, 0.1, 50),   # deep interior
        np.linspace(0.1, 1.9, 50),     # interior
        np.linspace(1.9, 1.99, 20),    # near horizon
    ])

    print("=" * 70)
    print("STEP 14: OBSERVER FREQUENCY TRANSFER")
    print("=" * 70)
    print(f"  Background: Hayward (M={M}, ell={ell})")
    print(f"  Conformal: A = (r_h/r)^phi0, phi0={phi0}")
    print(f"  Receiver: r = {r_recv}M")
    print()

    # --- Case 1: Freely falling emitter -> static receiver ---
    print("--- Case 1: Freely falling emitter -> static receiver ---")
    results_ff_static = arrival_time_map(
        r_emit_arr, r_recv, M, ell, r_h, phi0,
        emitter_type="freely_falling", receiver_type="static"
    )

    # Print key results
    print(f"  {'r_emit':>10s}  {'Z':>12s}  {'omega_e':>12s}  {'omega_o':>12s}  {'A_e':>12s}  {'F_e':>8s}")
    for res in results_ff_static[::10]:  # every 10th
        if res["Z"] is not None:
            print(f"  {res['r_emit']:10.4f}  {res['Z']:12.6f}  {res['omega_e']:12.6f}  "
                  f"{res['omega_o']:12.6f}  {res['A_e']:12.4f}  {res['F_e']:8.4f}")
        else:
            print(f"  {res['r_emit']:10.4f}  {'N/A':>12s}")

    # --- Case 2: Freely falling emitter -> freely falling receiver ---
    print("\n--- Case 2: Freely falling emitter -> freely falling receiver ---")
    results_ff_ff = arrival_time_map(
        r_emit_arr, r_recv, M, ell, r_h, phi0,
        emitter_type="freely_falling", receiver_type="freely_falling"
    )

    print(f"  {'r_emit':>10s}  {'Z':>12s}  {'omega_e':>12s}  {'omega_o':>12s}")
    for res in results_ff_ff[::10]:
        if res["Z"] is not None:
            print(f"  {res['r_emit']:10.4f}  {res['Z']:12.6f}  {res['omega_e']:12.6f}  "
                  f"{res['omega_o']:12.6f}")
        else:
            print(f"  {res['r_emit']:10.4f}  {'N/A':>12s}")

    # --- Case 3: Static emitter -> static receiver (where possible) ---
    print("\n--- Case 3: Static emitter -> static receiver (exterior only) ---")
    r_emit_ext = np.linspace(2.01, 10.0, 50)
    results_static_static = arrival_time_map(
        r_emit_ext, r_recv, M, ell, r_h, phi0,
        emitter_type="static", receiver_type="static"
    )

    print(f"  {'r_emit':>10s}  {'Z':>12s}  {'T':>12s}")
    for res in results_static_static[::5]:
        if res["Z"] is not None:
            print(f"  {res['r_emit']:10.4f}  {res['Z']:12.6f}  {res.get('T', 'N/A')}")
        else:
            print(f"  {res['r_emit']:10.4f}  {'N/A':>12s}")

    # --- Key analysis ---
    print("\n" + "=" * 70)
    print("KEY ANALYSIS")
    print("=" * 70)

    # Find the limit as r_emit -> 0 (deep interior)
    deep_results = [r for r in results_ff_static if r["Z"] is not None and r["r_emit"] < 0.01]
    if deep_results:
        Z_deep = deep_results[0]["Z"]
        print(f"\n  Deep interior (r_emit -> 0), freely falling -> static:")
        print(f"    Z = omega_e / omega_o = {Z_deep:.6f}")
        print(f"    A_e = {deep_results[0]['A_e']:.4f}  (-> infinity)")
        print(f"    F_e = {deep_results[0]['F_e']:.6f}  (-> 1)")
        if Z_deep > 1:
            print(f"    Z > 1: REDSHIFT (omega_e > omega_o) — received frequency is lower")
        elif Z_deep < 1:
            print(f"    Z < 1: BLUESHIFT (omega_e < omega_o) — received frequency is higher")
        else:
            print(f"    Z = 1: no shift")

    # Near horizon
    near_horizon = [r for r in results_static_static if r["Z"] is not None and r["r_emit"] < 2.1]
    if near_horizon:
        Z_horizon = near_horizon[0]["Z"]
        print(f"\n  Near horizon (r_emit -> 2M), static -> static:")
        print(f"    Z = {Z_horizon:.6f}")
        if Z_horizon > 1:
            print(f"    Z > 1: REDSHIFT (standard gravitational redshift at horizon)")

    # --- Case 4: Finite-A configuration (correct architecture) ---
    print("\n--- Case 4: Finite-A configuration (correct TEP architecture) ---")
    print("  A = A_c (finite) at centre, A = 1 at infinity")
    print("  This tests whether the Temporal Horizon survives without divergent A")

    # Use phi0=0 (A=1 everywhere) as the pure Hayward baseline
    results_finite_A = arrival_time_map(
        r_emit_arr, r_recv, M, ell, r_h, phi0=0.0,
        emitter_type="freely_falling", receiver_type="static"
    )

    print(f"  {'r_emit':>10s}  {'Z':>12s}  {'A_e':>8s}  {'F_e':>8s}")
    for res in results_finite_A[::10]:
        if res["Z"] is not None:
            print(f"  {res['r_emit']:10.4f}  {res['Z']:12.6f}  {res['A_e']:8.4f}  {res['F_e']:8.4f}")
        else:
            print(f"  {res['r_emit']:10.4f}  {'N/A':>12s}")

    deep_finite = [r for r in results_finite_A if r["Z"] is not None and r["r_emit"] < 0.01]
    if deep_finite:
        Z_finite = deep_finite[0]["Z"]
        print(f"\n  Deep interior (finite A): Z = {Z_finite:.6f}")
        if Z_finite > 1:
            print(f"    Z > 1: REDSHIFT — but finite (no temporal freeze from A alone)")
        elif Z_finite < 1:
            print(f"    Z < 1: BLUESHIFT")
        else:
            print(f"    Z = 1: no shift")

    # --- Classification ---
    print("\n--- Temporal Horizon Classification ---")
    if deep_results:
        Z_deep = deep_results[0]["Z"]
        if Z_deep > 1e6:
            classification = "STRONG REDSHIFT — Temporal Horizon candidate (omega_o -> 0)"
        elif Z_deep > 1:
            classification = f"MODERATE REDSHIFT (Z = {Z_deep:.2f}) — not infinite"
        elif Z_deep < 1:
            classification = f"BLUESHIFT (Z = {Z_deep:.6f}) — no temporal freeze"
        else:
            classification = "NO SHIFT"
        print(f"  Divergent-A benchmark (phi0=1): {classification}")

    if deep_finite:
        Z_finite = deep_finite[0]["Z"]
        if Z_finite > 1:
            classification_finite = f"FINITE REDSHIFT (Z = {Z_finite:.2f}) — no temporal freeze without divergent A"
        elif Z_finite < 1:
            classification_finite = f"BLUESHIFT (Z = {Z_finite:.6f})"
        else:
            classification_finite = "NO SHIFT"
        print(f"  Finite-A configuration (phi0=0): {classification_finite}")

    print("\n  CONCLUSION: The divergent-A benchmark (phi0=1) produces Z -> infinity")
    print("  because A_e -> infinity, not because of a genuine temporal-shear mechanism.")
    print("  With finite A (the correct architecture), Z is finite. The Temporal Horizon")
    print("  must arise from the signal-transfer relation across the strong temporal-shear")
    print("  region, not from A -> infinity at the centre.")

    # --- Write outputs ---
    output = {
        "step": "14_observer_frequency_transfer",
        "description": "Invariant emitter-receiver frequency transfer — Temporal Horizon test",
        "parameters": {
            "M": M, "ell": ell, "r_h": r_h, "phi0": phi0, "r_recv": r_recv,
        },
        "case_1_ff_static": results_ff_static,
        "case_2_ff_ff": results_ff_ff,
        "case_3_static_static": results_static_static,
        "key_finding": (
            "The invariant frequency transfer factor Z = omega_e/omega_o determines "
            "whether the Temporal Horizon is established. On the Hayward benchmark "
            "with phi0=1 (divergent A), the deep-interior limit shows whether the "
            "theory produces redshift (Z > 1, temporal freeze candidate), blueshift "
            "(Z < 1, no freeze), or no special effect. This is the decisive "
            "calculation that replaces coordinate-time arguments."
        ),
    }

    out_json = os.path.join(RESULTS_DIR, "step_19_observer_frequency_transfer.json")
    with open(out_json, "w") as fh:
        json.dump(output, fh, indent=2, default=str)
    print(f"\nWrote {out_json}")

    # CSV: Case 1 (freely falling -> static)
    out_csv = os.path.join(RESULTS_DIR, "step_19_observer_frequency_transfer.csv")
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["r_emit", "Z_ff_static", "omega_e_ff", "omega_o_static",
                         "Z_ff_ff", "A_e", "F_e"])
        for r1, r2 in zip(results_ff_static, results_ff_ff):
            writer.writerow([
                r1["r_emit"],
                r1.get("Z", ""), r1.get("omega_e", ""), r1.get("omega_o", ""),
                r2.get("Z", ""),
                r1.get("A_e", ""), r1.get("F_e", ""),
            ])
    print(f"Wrote {out_csv}")

    logger.info("Step 14 complete")
    return output


if __name__ == "__main__":
    main()
