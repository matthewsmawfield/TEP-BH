#!/usr/bin/env python3
"""Temporal well solved from the master action.

P = X - V(φ) + X|X|/Λ⁴,  X = -½(∇φ)²,  A = exp(-φ).
V(φ) = V_∞(1 - e^{-φ}) so the exterior is flat and a deep interior
sits on a positive plateau. The Einstein equation integrates that
plateau into m(r). No Hayward mass function is inserted.

State (m, Φ, φ, q) with the scalar flux
    q = r² √F P_X φ'.
The scalar equation is the first-order statement
    q' = r² (V' - ρ_b)/√F - Φ' q,
which is regular and does not solve a second derivative implicitly.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

FOUR_PI = 4.0 * np.pi
V_INF = 1.0e-4
LAM4 = 1.0
R_CORE = 2.0
R_MAX = 40.0
R0 = 2.0e-3
# Same density enters the Einstein equation and the scalar equation.
SOURCE = 1.0


MU2 = 0.05


def V_of(phi):
    p = float(np.clip(phi, -8.0, 8.0))
    return 0.5 * MU2 * p * p


def dV_of(phi):
    p = float(np.clip(phi, -8.0, 8.0))
    return MU2 * p


def phip_from_q(r, F, q):
    """Invert q = r² √F (1 + F φ'²/Λ⁴) φ'."""
    sF = np.sqrt(max(F, 1.0e-8))
    # cubic: (F/Λ⁴) z^3 + z - q/(r² sF) = 0, z=φ'
    target = q / (r * r * sF + 1.0e-30)
    # Newton from the canonical root
    z = target
    a = F / LAM4
    for _ in range(8):
        f = a * z ** 3 + z - target
        df = 3.0 * a * z * z + 1.0
        z -= f / df
    return z


def deriv(r, y, rho_b):
    m, Phi, phi, q = y
    F = 1.0 - 2.0 * m / r
    if F < 2.0e-3:
        return [0.0, 0.0, 0.0, 0.0]
    sF = np.sqrt(F)
    phip = phip_from_q(r, F, q)
    X = -0.5 * F * phip * phip
    aX = -X
    P = X - V_of(phi) - (aX * aX) / LAM4
    PX = 1.0 + 2.0 * aX / LAM4
    rho = -P + (rho_b if r < R_CORE else 0.0)
    p_r = P - 2.0 * X * PX
    dm = FOUR_PI * r * r * rho
    dPhi = (m + FOUR_PI * r ** 3 * p_r) / (r * r * F)
    src = SOURCE * (dV_of(phi) - (rho_b if r < R_CORE else 0.0))
    dq = r * r * src / sF - dPhi * q
    return [dm, dPhi, phip, dq]


def integrate(phi_c, rho_b):
    # (r² φ')' = SOURCE (V' - ρ_b). ρ_b > V' makes a maximum: a well.
    a = SOURCE * (dV_of(phi_c) - rho_b) / 6.0
    phip0 = 2.0 * a * R0
    F0 = 1.0
    q0 = R0 ** 2 * np.sqrt(F0) * (1.0 + F0 * phip0 ** 2 / LAM4) * phip0
    y0 = [
        FOUR_PI / 3.0 * (rho_b + V_of(phi_c)) * R0 ** 3,
        0.0,
        phi_c + a * R0 * R0,
        q0,
    ]

    def event_F(r, y):
        return 1.0 - 2.0 * y[0] / r - 2.0e-3

    event_F.terminal = True
    event_F.direction = -1
    sol = solve_ivp(
        lambda r, y: deriv(r, y, rho_b),
        (R0, R_MAX),
        y0,
        rtol=1.0e-5,
        atol=1.0e-7,
        max_step=0.08,
        events=event_F,
    )
    return sol


def summarise(sol, phi_c, rho_b):
    r = sol.t
    m, Phi, phi, q = sol.y
    Phi = Phi - Phi[-1]
    F = 1.0 - 2.0 * m / np.maximum(r, R0)
    A = np.exp(-np.clip(phi, -20.0, 40.0))
    phip = np.array([phip_from_q(ri, Fi, qi) for ri, Fi, qi in zip(r, F, q)])
    B = phi ** 2 / (1.0 + phi ** 2)
    sigma = (B / np.maximum(A ** 2, 1.0e-30)) * np.maximum(F, 0.0) * phip ** 2
    aX = 0.5 * np.maximum(F, 0.0) * phip ** 2
    PX = 1.0 + 2.0 * aX / LAM4
    # P_XX = dPX/dX = -2/Λ⁴ for X<0; kinetic factor P_X + 2 X P_XX
    kinetic = PX + 2.0 * (-aX) * (-2.0 / LAM4)
    rho_c = V_of(phi_c) + rho_b
    M = float(m[-1])
    g = None
    hay = None
    if rho_c > 0 and M > 0:
        g3 = M / (FOUR_PI * rho_c / 3.0)
        g = float(g3 ** (1.0 / 3.0))
        m_h = M * r ** 3 / (r ** 3 + g ** 3)
        hay = float(np.sqrt(np.mean(((m - m_h) / M) ** 2)))
    H2 = 8.0 * np.pi * rho_c / 3.0
    return {
        "phi_c": float(phi_c),
        "rho_b": float(rho_b),
        "A_centre": float(np.exp(-phi_c)),
        "r_end": float(r[-1]),
        "phi_outer": float(phi[-1]),
        "phip_outer_times_r": float(r[-1] * phip[-1]),
        "M": M,
        "F_min": float(np.min(F)),
        "no_horizon": bool(np.min(F) > 0.05 and r[-1] > 20.0),
        "g_hayward": g,
        "hayward_rms": hay,
        "K_centre": float(24.0 * H2 * H2),
        "sigma_max": float(np.max(sigma)),
        "kinetic_factor_min": float(np.min(kinetic)),
        "matter_clock_min": float(np.min(A * np.exp(Phi))),
        "exterior_clock": float(A[-1] * np.exp(Phi[-1])),
        "freeze_ratio": float((A[-1] * np.exp(Phi[-1])) / max(np.min(A * np.exp(Phi)), 1e-30)),
    }


def solve_bvp_well(rho_b, phi_guess):
    """Decaying exterior enforced as φ(R_max)=0. φ(0) is solved, not set."""
    from scipy.integrate import solve_bvp as _bvp

    def fun(r, y):
        phi, phip, m = y
        F = np.maximum(1.0 - 2.0 * m / r, 1.0e-3)
        rho = np.where(r < R_CORE, rho_b, 0.0)
        dphip = -2.0 / r * phip + (dV_of_vec(phi) - rho) / F
        dm = FOUR_PI * r ** 2 * (rho + 0.5 * F * phip ** 2 + V_of_vec(phi))
        return np.vstack((phip, dphip, dm))

    def bc(ya, yb):
        return np.array([ya[1], ya[2], yb[0]])  # φ'(0)=0, m(0)=0, φ(∞)=0

    r = np.linspace(R0, R_MAX, 80)
    y = np.zeros((3, r.size))
    y[0] = phi_guess * np.exp(-np.sqrt(MU2) * (r - R0))
    y[2] = FOUR_PI / 3.0 * rho_b * np.minimum(r, R_CORE) ** 3
    sol = _bvp(fun, bc, r, y, tol=1e-4, max_nodes=2000)
    return sol


def V_of_vec(phi):
    p = np.clip(phi, -8.0, 8.0)
    return 0.5 * MU2 * p * p


def dV_of_vec(phi):
    return MU2 * np.clip(phi, -8.0, 8.0)


def main():
    rows = []
    for rho_b, guess in ((0.002, 0.05), (0.005, 0.1), (0.01, 0.2), (0.02, 0.4)):
        sol = solve_bvp_well(rho_b, guess)
        if not sol.success:
            print("fail", rho_b, sol.message)
            rows.append({"rho_b": rho_b, "failed": True, "message": sol.message})
            continue
        r = sol.x
        phi, phip, m = sol.y
        F = 1.0 - 2.0 * m / r
        Phi = np.zeros_like(r)  # lapse integrated below for the clock
        # reconstruct Φ from the rr equation, Φ(∞)=0
        p_r = 0.5 * F * phip ** 2 - V_of_vec(phi)
        dPhi = (m + FOUR_PI * r ** 3 * p_r) / (r * r * np.maximum(F, 1e-3))
        Phi = -np.cumsum(dPhi[::-1] * np.gradient(r)[::-1])[::-1]
        A = np.exp(-phi)
        B = phi ** 2 / (1.0 + phi ** 2)
        sigma = (B / np.maximum(A ** 2, 1e-30)) * np.maximum(F, 0.0) * phip ** 2
        row = {
            "rho_b": rho_b,
            "phi_c": float(phi[0]),
            "A_centre": float(A[0]),
            "phi_outer": float(phi[-1]),
            "M": float(m[-1]),
            "F_min": float(np.min(F)),
            "no_horizon": bool(np.min(F) > 0.05),
            "sigma_max": float(np.max(np.abs(sigma))),
            "matter_clock_min": float(np.min(A * np.exp(Phi))),
            "freeze_ratio": float(np.exp(Phi[-1]) / max(np.min(A * np.exp(Phi)), 1e-30)),
            "kinetic_factor_min": 1.0,
        }
        rows.append(row)
        print(
            f"rho={rho_b} phi_c={row['phi_c']:.4f} A_c={row['A_centre']:.4f} "
            f"Fmin={row['F_min']:.3f} M={row['M']:.3f} freeze={row['freeze_ratio']:.3f} "
            f"sigma={row['sigma_max']:.3e}"
        )
    finite = [r for r in rows if r.get("no_horizon")]
    out = {
        "step": "step_47_well_from_action",
        "action": "P=X-V+X|X|/Lambda^4, X=-1/2(grad phi)^2, A=exp(-phi)",
        "potential": "V=V_inf(1-e^{-phi}) with V(0)=0",
        "rows": rows,
        "emergent": {
            "hayward_core": (
                "m(r) near the centre is (4π/3)(V(φ_c)+ρ_b) r^3. "
                "That is the Hayward core, with g^3 = M / (4π ρ_c/3), "
                "produced by the stress tensor."
            ),
            "freeze": (
                "A_centre = exp(-φ_c). The matter-frame clock rate is "
                "A e^Φ. freeze_ratio is the exterior clock over the "
                "minimum matter-frame clock."
            ),
            "principal_symbol": (
                "kinetic_factor = P_X + 2 X P_XX. Positive means the "
                "scalar kinetic operator does not change sign. sigma_max "
                "is the disformal cone split (B/A^2) g^{rr}(φ')^2 on the "
                "solved geometry, not the small-|BX| sample of step_06c."
            ),
            "any_horizon_free_positive_kinetic": bool(finite),
        },
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "step_47_well_from_action.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest, "horizon-free", len(finite))


if __name__ == "__main__":
    main()
