#!/usr/bin/env python3
"""Analyze the determinant condition for the disformal temporal branch.

det_2d = -A^2 * [A^2 + B*(F*psi'^2 + 2*q*psi')]

For Lorentzian: det_2d < 0, so we need:
    A^2 + B*(F*psi'^2 + 2*q*psi') > 0

Inside the horizon (F < 0):
    F*psi'^2 < 0  (always negative)
    2*q*psi' sign depends on q and psi'

The problem: if B grows large in the interior and both terms inside B
are negative, the bracket goes negative → signature change.

Key insight: we need 2*q*psi' to be POSITIVE and LARGER than |F*psi'^2|.
That means q and psi' must have the SAME SIGN.

Current psi(r) = phi_0 * ln(r/r_h) * S(r)
Inside horizon: ln(r/r_h) < 0, S ~ 1, so psi < 0 and psi' < 0 (becoming
more negative as r decreases).

So we need q < 0 to make 2*q*psi' > 0.

But wait — if q < 0, then B*q^2 is still positive (good for temporal
stretching of gtilde_vv), and 2*q*psi' > 0 (good for determinant).

Let me also check: what if we use a DIFFERENT psi(r) that has psi' > 0?
For example, psi(r) = -phi_0 * ln(r/r_h) * S(r), which is positive inside
and has psi' > 0. Then q > 0 gives 2*q*psi' > 0.

Let me test both approaches.
"""

import numpy as np

def analyze_determinant_condition():
    """Systematically analyze what parameter choices keep det_2d < 0."""

    r = np.logspace(np.log10(1e-8), np.log10(50.0), 20000)
    M = 1.0
    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    # Current psi profile
    phi_0 = 2.0
    delta = 0.05
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)

    # Option A: psi = phi_0 * ln(r/r_h) * S  (current: psi < 0, psi' < 0 inside)
    psi_A = phi_0 * np.log(r_safe / r_h) * S
    dpsi_A = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    # Option B: psi = -phi_0 * ln(r/r_h) * S  (flipped: psi > 0, psi' > 0 inside)
    psi_B = -phi_0 * np.log(r_safe / r_h) * S
    dpsi_B = -phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    idx_h = np.argmin(np.abs(r - r_h))
    idx_inner = 0  # innermost

    print("=" * 80)
    print("DETERMINANT CONDITION ANALYSIS")
    print("det_2d = -A^2 * [A^2 + B*(F*psi'^2 + 2*q*psi')]")
    print("Need: A^2 + B*(F*psi'^2 + 2*q*psi') > 0")
    print("=" * 80)

    for label, psi, dpsi in [("Option A (psi<0, psi'<0)", psi_A, dpsi_A),
                              ("Option B (psi>0, psi'>0)", psi_B, dpsi_B)]:
        print(f"\n--- {label} ---")
        print(f"  At innermost (r={r[idx_inner]:.2e}):")
        print(f"    psi = {psi[idx_inner]:.6e}")
        print(f"    psi' = {dpsi[idx_inner]:.6e}")
        print(f"    F = {F[idx_inner]:.6e}")
        print(f"    F*psi'^2 = {F[idx_inner] * dpsi[idx_inner]**2:.6e}")
        print(f"  At horizon (r={r[idx_h]:.2e}):")
        print(f"    psi = {psi[idx_h]:.6e}")
        print(f"    psi' = {dpsi[idx_h]:.6e}")
        print(f"    F = {F[idx_h]:.6e}")
        print(f"    F*psi'^2 = {F[idx_h] * dpsi[idx_h]**2:.6e}")

        # For each option, what sign of q makes 2*q*psi' > 0?
        if dpsi[idx_inner] < 0:
            print(f"  psi' < 0 inside → need q < 0 for 2*q*psi' > 0")
            q_sign = "negative"
        else:
            print(f"  psi' > 0 inside → need q > 0 for 2*q*psi' > 0")
            q_sign = "positive"

        # Check: can 2*q*psi' dominate F*psi'^2?
        # Need: 2*q*psi' > |F*psi'^2| = |F|*psi'^2
        # i.e., 2*q > |F|*psi'  (dividing by psi' with appropriate sign)
        # i.e., |q| > |F|*|psi'| / 2

        # At innermost:
        F_inner = F[idx_inner]
        dpsi_inner = dpsi[idx_inner]
        q_min = np.abs(F_inner) * np.abs(dpsi_inner) / 2.0
        print(f"  Minimum |q| to dominate F*psi'^2 at innermost: {q_min:.6e}")

        # At several interior radii
        for r_test in [0.1, 0.5, 1.0, 1.5, 1.9]:
            idx_t = np.argmin(np.abs(r - r_test))
            F_t = F[idx_t]
            dpsi_t = dpsi[idx_t]
            q_min_t = np.abs(F_t) * np.abs(dpsi_t) / 2.0
            print(f"  r={r_test:.1f}: F={F_t:.4f}, psi'={dpsi_t:.6e}, "
                  f"|q|_min={q_min_t:.6e}")

    # Now test: with the right sign of q, does the determinant stay negative?
    print("\n" + "=" * 80)
    print("TESTING WITH CORRECT q SIGN")
    print("=" * 80)

    # Use Option A (psi < 0, psi' < 0) with q < 0
    # Or Option B (psi > 0, psi' > 0) with q > 0
    # Both give 2*q*psi' > 0

    for label, psi, dpsi, q_val in [
        ("Option A, q=-1", psi_A, dpsi_A, -1.0),
        ("Option A, q=-5", psi_A, dpsi_A, -5.0),
        ("Option A, q=-10", psi_A, dpsi_A, -10.0),
        ("Option A, q=-50", psi_A, dpsi_A, -50.0),
        ("Option B, q=+1", psi_B, dpsi_B, 1.0),
        ("Option B, q=+5", psi_B, dpsi_B, 5.0),
        ("Option B, q=+10", psi_B, dpsi_B, 10.0),
        ("Option B, q=+50", psi_B, dpsi_B, 50.0),
    ]:
        # A = 1 (pure disformal), B = B0 * |psi|^2 / (1 + |psi|^2)
        B0 = 50.0
        A = 1.0
        A2 = 1.0
        abs_psi = np.abs(psi)
        B = B0 * abs_psi ** 2 / (1.0 + abs_psi ** 2)

        bracket = A2 + B * (F * dpsi ** 2 + 2 * q_val * dpsi)
        det_2d = -A2 * bracket

        frac_neg = np.sum(det_2d < 0) / len(det_2d)
        sign_changes = np.where(np.signbit(det_2d[:-1]) != np.signbit(det_2d[1:]))[0]

        # Check inside horizon specifically
        det_inside = det_2d[:idx_h]
        frac_neg_inside = np.sum(det_inside < 0) / len(det_inside)

        # Check the temporal coefficient
        gtilde_vv = -A2 * F + B * q_val ** 2
        gtilde_vv_innermost = gtilde_vv[idx_inner]

        print(f"\n  {label}:")
        print(f"    frac Lorentzian (all): {frac_neg:.4f}")
        print(f"    frac Lorentzian (inside): {frac_neg_inside:.4f}")
        print(f"    sign changes: {len(sign_changes)}")
        print(f"    gtilde_vv at innermost: {gtilde_vv_innermost:.6e}")
        print(f"    bracket at innermost: {bracket[idx_inner]:.6e}")
        print(f"    bracket at horizon: {bracket[idx_h]:.6e}")

        if len(sign_changes) > 0:
            r_sc = r[sign_changes[0]]
            print(f"    first sign change at r = {r_sc:.6e}")

    # The issue: F*psi'^2 grows as 1/r^2 near r=0 (since psi' ~ 1/r)
    # while 2*q*psi' also grows as 1/r. So F*psi'^2 ~ 1/r^2 dominates
    # 2*q*psi' ~ 1/r for small r. This means the bracket ALWAYS goes
    # negative near r=0 if B is nonzero there.
    #
    # SOLUTION: make B -> 0 as r -> 0 (just like the conformal prototype!),
    # but in a way that the TEMPORAL effect (B*q^2 in gtilde_vv) persists.
    #
    # Wait — if B -> 0, then B*q^2 -> 0 too, so no temporal effect.
    #
    # The real solution: we need B to be large enough to create temporal
    # stretching but small enough that the determinant stays negative.
    # OR: we need a different psi(r) where psi' doesn't blow up as 1/r.
    #
    # Key idea: use psi(r) that SATURATES (psi' -> 0) as r -> 0.
    # Then F*psi'^2 -> 0 and 2*q*psi' -> 0, so the bracket -> A^2 > 0.
    # But B can still be large (saturating), so B*q^2 gives temporal stretch.
    #
    # Example: psi(r) = phi_0 * (1 - r/r_h) * S(r)  — linear, saturates
    # Or: psi(r) = phi_0 * tanh((r_h - r)/delta) * S(r) — smooth saturation

    print("\n" + "=" * 80)
    print("TESTING WITH SATURATING psi(r)")
    print("=" * 80)

    # Saturating psi: psi(r) = phi_0 * tanh((r_h - r)/(delta*r_h)) * S(r)
    # This gives psi -> phi_0 (finite) as r -> 0, and psi' -> 0 as r -> 0.
    psi_sat = phi_0 * np.tanh(np.maximum((r_h - r) / (delta * r_h), 0)) * S
    # psi' for saturating profile
    sech2 = 1.0 / np.cosh(np.maximum((r_h - r) / (delta * r_h), 0)) ** 2
    dpsi_sat = phi_0 * (-sech2 / (delta * r_h)) * S + phi_0 * np.tanh(
        np.maximum((r_h - r) / (delta * r_h), 0)) * dS

    print(f"\n  Saturating psi at innermost: {psi_sat[idx_inner]:.6e}")
    print(f"  Saturating psi' at innermost: {dpsi_sat[idx_inner]:.6e}")
    print(f"  Saturating psi at horizon: {psi_sat[idx_h]:.6e}")
    print(f"  Saturating psi' at horizon: {dpsi_sat[idx_h]:.6e}")

    for q_val in [1.0, 5.0, 10.0, 50.0]:
        B0 = 50.0
        A = 1.0
        A2 = 1.0
        abs_psi = np.abs(psi_sat)
        B = B0 * abs_psi ** 2 / (1.0 + abs_psi ** 2)

        bracket = A2 + B * (F * dpsi_sat ** 2 + 2 * q_val * dpsi_sat)
        det_2d = -A2 * bracket

        frac_neg = np.sum(det_2d < 0) / len(det_2d)
        sign_changes = np.where(np.signbit(det_2d[:-1]) != np.signbit(det_2d[1:]))[0]
        det_inside = det_2d[:idx_h]
        frac_neg_inside = np.sum(det_inside < 0) / len(det_inside)

        gtilde_vv = -A2 * F + B * q_val ** 2
        gtilde_vv_innermost = gtilde_vv[idx_inner]

        # Areal radius
        areal = A * r  # since A=1, areal = r (finite!)

        print(f"\n  Saturating psi, q={q_val}:")
        print(f"    frac Lorentzian (all): {frac_neg:.4f}")
        print(f"    frac Lorentzian (inside): {frac_neg_inside:.4f}")
        print(f"    sign changes: {len(sign_changes)}")
        print(f"    gtilde_vv at innermost: {gtilde_vv_innermost:.6e}")
        print(f"    areal at innermost: {areal[idx_inner]:.6e}")
        print(f"    areal max: {np.max(areal):.6e}")

        if len(sign_changes) > 0:
            r_sc = r[sign_changes[0]]
            print(f"    first sign change at r = {r_sc:.6e}")
        else:
            print(f"    NO SIGN CHANGES — determinant stays negative!")

        # Check bracket at key radii
        for r_test in [1e-8, 1e-4, 0.1, 1.0, 1.9]:
            idx_t = np.argmin(np.abs(r - r_test))
            print(f"    r={r_test:.1e}: bracket={bracket[idx_t]:.6e}, "
                  f"B={B[idx_t]:.6e}, psi'={dpsi_sat[idx_t]:.6e}")


if __name__ == '__main__':
    analyze_determinant_condition()
