#!/usr/bin/env python3
"""Focus on the timelike completeness problem.

With phi_0 = 0.9:
  - Areal radius: finite (rho ~ r^{0.1} -> 0) ✓
  - Lorentzian: yes (0 sign changes with q=10, B0=1) ✓
  - Curvature: K ~ r^{4.8} -> 0 ✓
  - Null: complete (integrand ~ r^{-1.8}, integral ~ r^{-0.8} -> inf) ✓
  - Timelike: INCOMPLETE (power law 0.4, need >= 1.0) ✗

The timelike integrand in the conformal regime is:
  dtau/dr ~ 1/(A * sqrt(r/(2M))) ~ r^{phi_0 - 0.5}

For phi_0 = 0.9: dtau/dr ~ r^{0.4}, integral ~ r^{1.4} -> 0. INCOMPLETE.

But with the disformal term (q != 0), the metric is NOT purely conformal.
The exact timelike equation is:
  rdot^2 = -(gvv + 1) / det_2d

Inside the horizon (F < 0):
  gvv = -A^2*F + B*q^2 = A^2*|F| + B*q^2  (both positive inside)
  det_2d = -A^2 * [A^2 + B*(F*psi'^2 + 2*q*psi')]

In the deep interior (B -> 0 due to Gaussian damping):
  gvv -> A^2 * |F| ~ A^2 * 2M/r
  det_2d -> -A^4
  rdot^2 = -(A^2*|F| + 1) / (-A^4) = (A^2*|F| + 1) / A^4
         ~ A^2 * 2M/r / A^4 = 2M / (r * A^2)
         ~ 2M / (r * r^{-2*phi_0}) = 2M * r^{2*phi_0 - 1}

  rdot ~ r^{phi_0 - 0.5}
  dtau/dr ~ r^{0.5 - phi_0}
  integral ~ r^{1.5 - phi_0} -> 0 if phi_0 < 1.5. INCOMPLETE.

So in the deep interior where B -> 0, the timelike geodesic is the same
as the conformal case, and it's incomplete for phi_0 < 1.5.

The ONLY way to fix timelike completeness with bounded areal radius is
to have B NOT vanish in the deep interior, so the disformal term
contributes to the timelike dynamics.

If B stays finite (say B -> B_inf > 0) as r -> 0, then:
  gvv = A^2*|F| + B_inf*q^2 ~ A^2*2M/r + B_inf*q^2
  For the B_inf*q^2 term to dominate: B_inf*q^2 >> A^2*2M/r
  i.e., B_inf*q^2 >> r^{-2*phi_0} * r^{-1} = r^{-2*phi_0-1}
  This requires B_inf*q^2 to grow faster than r^{-2*phi_0-1}, which
  means B_inf can't be constant — it needs to grow.

But if B grows, the determinant condition becomes:
  A^2 + B*(F*psi'^2 + 2*q*psi') > 0

With psi' -> 0 (saturating psi) and F -> -inf:
  F*psi'^2 -> 0 (if psi' -> 0 fast enough)
  2*q*psi' -> 0
  So bracket -> A^2 > 0. OK!

So the plan: use SATURATING psi (psi' -> 0) + B that GROWS in the interior
+ mild conformal A (phi_0 ~ 0.9) + q != 0.

With saturating psi, psi' -> 0 in the deep interior, so:
  - The determinant bracket -> A^2 > 0 (always Lorentzian) ✓
  - B can grow without breaking the determinant ✓
  - gvv = A^2*|F| + B*q^2, and if B*q^2 >> A^2*|F|, then:
    gvv ~ B*q^2
    det_2d ~ -A^2 * (A^2 + B*0) = -A^4 (since psi' -> 0, the B terms vanish)
    Wait, that's not right. Let me recalculate.

    det_2d = gvv*grr - gvr^2
    gvv = -A^2*F + B*q^2
    gvr = A^2 + B*q*psi'
    grr = B*psi'^2

    As psi' -> 0:
    grr -> 0
    gvr -> A^2
    gvv -> A^2*|F| + B*q^2

    det_2d = (A^2*|F| + B*q^2) * 0 - (A^2)^2 = -A^4

    So det_2d -> -A^4 regardless of B. Good, always Lorentzian.

    But grr -> 0 means the radial metric component vanishes!
    The inverse metric ginv_rr = gvv/det_2d = (A^2*|F| + B*q^2) / (-A^4)
    = -(A^2*|F| + B*q^2) / A^4

    For timelike: rdot^2 = -(gvv + 1) / det_2d = (gvv + 1) / A^4
    gvv = A^2*|F| + B*q^2

    If B*q^2 >> A^2*|F| (B growing, q nonzero):
    gvv ~ B*q^2
    rdot^2 ~ B*q^2 / A^4
    rdot ~ sqrt(B) * q / A^2

    dtau/dr ~ A^2 / (q * sqrt(B))

    For timelike completeness: need integral of A^2/(q*sqrt(B)) dr to diverge.
    A^2 ~ r^{-2*phi_0}, so need:
    integral of r^{-2*phi_0} / sqrt(B(r)) dr to diverge.

    If B ~ r^{-beta} (growing), then:
    integrand ~ r^{-2*phi_0 + beta/2}
    integral diverges if -2*phi_0 + beta/2 <= -1
    i.e., beta >= 2*(2*phi_0 - 1) = 4*phi_0 - 2

    For phi_0 = 0.9: beta >= 4*0.9 - 2 = 1.6
    So B needs to grow at least as r^{-1.6} for timelike completeness.

    But we also need B to not break the determinant in the transition zone.
    With saturating psi, psi' -> 0 in the deep interior, so B can grow
    freely there. The constraint is only in the transition zone near the
    horizon where psi' is large.

    PLAN:
    1. Use saturating psi: psi = phi_0 * tanh((r_h-r)/(delta*r_h)) * S(r)
    2. Use mild conformal: A = exp(-psi) with phi_0 = 0.9
       (A ~ r^{-0.9} in the transition, saturating to A ~ exp(-phi_0) in deep interior)
       Wait — with saturating psi, psi -> phi_0 as r -> 0, so A -> exp(-phi_0) = const!
       That means A is BOUNDED, not growing. Good for areal radius.
       But then A^2*|F| ~ const * 2M/r -> inf, and we need B*q^2 to dominate.

    Actually, with saturating psi:
      psi -> phi_0 (finite), A -> exp(beta_A * phi_0) = const
      psi' -> 0
      F -> -2M/r -> -inf
      |F| -> 2M/r -> inf

    So gvv = A^2*|F| + B*q^2 ~ A_0^2 * 2M/r + B*q^2
    The A^2*|F| term grows as 1/r. For B*q^2 to dominate, B must grow
    faster than 1/r.

    det_2d = -A^4 (since psi' -> 0, grr -> 0, gvr -> A^2)
    So det_2d -> -A_0^4 = const. Lorentzian. ✓

    rdot^2 = (gvv + 1) / A_0^4 ~ (A_0^2 * 2M/r + B*q^2) / A_0^4

    If B ~ r^{-beta} with beta > 1:
    rdot^2 ~ B*q^2 / A_0^4 ~ r^{-beta} * q^2 / A_0^4
    rdot ~ r^{-beta/2}
    dtau/dr ~ r^{beta/2}
    integral ~ r^{beta/2 + 1} -> 0. STILL INCOMPLETE!

    Wait, that's wrong. If rdot ~ r^{-beta/2}, then dr/dtau ~ r^{-beta/2},
    so dtau = dr / rdot = r^{beta/2} dr, and integral of r^{beta/2} dr
    = r^{beta/2 + 1} / (beta/2 + 1) -> 0 as r -> 0. INCOMPLETE.

    The issue: faster rdot means LESS proper time, not more!
    For completeness, we need rdot -> 0 (particle slows down),
    so dtau = dr/rdot -> inf.

    rdot^2 = (gvv + 1) / A_0^4
    For rdot -> 0: need gvv + 1 -> 0, i.e., gvv -> -1.
    But inside the horizon, gvv = A^2*|F| + B*q^2 > 0 (both terms positive).
    So gvv + 1 > 1 > 0, and rdot^2 > 1/A_0^4 > 0.
    rdot is bounded below by 1/A_0^2. The particle never slows down!

    This means: with bounded A and B growing, the timelike geodesic
    reaches r=0 in FINITE proper time. The singularity is NOT avoided
    for timelike curves.

    FUNDAMENTAL PROBLEM: Inside the horizon, gvv > 0 (spacelike direction).
    The timelike direction is r. For geodesic completeness, we need the
    r-direction to "stretch" so that it takes infinite proper time to
    reach r=0. This requires the radial metric to blow up, which means
    either A or B*psi'^2 must diverge.

    But if psi' -> 0 (saturating psi), then B*psi'^2 -> 0, and the
    radial metric is controlled by A alone. With bounded A, the radial
    metric is bounded, and timelike geodesics are incomplete.

    CONCLUSION: We CANNOT have all of:
    1. Bounded A (finite areal radius)
    2. Saturating psi (psi' -> 0, for determinant safety)
    3. Timelike geodesic completeness

    At least one must be relaxed. The question is which.

    Option 1: Let A grow (but mildly, phi_0 < 1) — areal radius -> 0 but finite.
              Need phi_0 >= 1.5 for timelike completeness. But then areal radius
              diverges for phi_0 > 1. DEAD END with pure conformal.

    Option 2: Let psi' NOT saturate — keep psi' ~ 1/r (logarithmic psi).
              Then B*psi'^2 can diverge, providing radial stretching.
              But then F*psi'^2 ~ 1/r^3 inside, which breaks the determinant
              unless A^2 dominates. Need A^2 >> B*|F|*psi'^2, which means
              A must grow — back to the conformal prototype.

    Option 3: Use a DIFFERENT mechanism for timelike completeness.
              Instead of stretching the radial metric, make the effective
              potential create a barrier. If there's a turning point at
              r > 0, the geodesic never reaches r=0.
              But this would mean the singularity is avoided by a potential
              barrier, not by geodesic completeness. Different physics.

    Option 4: Accept that timelike geodesics reach r=0 in finite proper
              time, but argue that the curvature is finite there (no
              singularity). This is the "regular black hole" approach:
              geodesically incomplete but curvature-regular.
              But the user requires geodesic completeness (Gate 6).

    Option 5: Use a non-saturating psi with controlled psi' that keeps
              the determinant safe. The key is to have psi' grow, but
              not too fast. If psi' ~ r^{-gamma} with gamma < 1, then:
              - B*psi'^2 ~ r^{-2*gamma} (grows, provides radial stretching)
              - F*psi'^2 ~ r^{-1-2*gamma} (also grows)
              - Need A^2 + B*(F*psi'^2 + 2*q*psi') > 0
              - If B is small enough, A^2 dominates.

              With A ~ r^{-phi_0} and B*F*psi'^2 ~ r^{-1-2*gamma}:
              Need r^{-2*phi_0} >> r^{-1-2*gamma}
              i.e., 2*phi_0 < 1 + 2*gamma
              i.e., gamma > phi_0 - 0.5

              For phi_0 = 0.9: gamma > 0.4
              For timelike completeness via B*psi'^2:
                grr = B*psi'^2 ~ r^{-2*gamma}
                rdot^2 ~ gvv / (-det_2d) ~ A^2*|F| / A^4 = |F|/A^2 ~ r^{-1+2*phi_0}
                Hmm, this doesn't depend on B*psi'^2 in the conformal regime.

              Actually, let me think about this differently.
              The timelike proper time is:
                dtau = sqrt(-g_{mu nu} dx^mu dx^nu)
              For radial infall: dtau^2 = -gvv*dv^2 - 2*gvr*dv*dr - grr*dr^2
              With E = -(gvv*v' + gvr*r') = 1:
                The proper time per dr is:
                dtau = sqrt(-det_2d / (gvv + 1)) * dr / |E|
                     = sqrt(A^4 / (gvv + 1)) * dr  (in conformal regime)
                     = A^2 / sqrt(gvv + 1) * dr

              In the deep interior with bounded A and B -> 0:
                gvv ~ A_0^2 * 2M/r -> inf
                dtau ~ A_0^2 / sqrt(A_0^2 * 2M/r) * dr = A_0 * sqrt(r/(2M)) * dr
                integral ~ r^{3/2} -> 0. INCOMPLETE.

              With B NOT vanishing (B -> B_inf > 0) and psi' NOT vanishing:
                gvv ~ A_0^2 * 2M/r + B_inf * q^2
                grr ~ B_inf * psi'^2
                det_2d = gvv*grr - gvr^2

                If psi' ~ r^{-gamma}:
                grr ~ B_inf * r^{-2*gamma}
                gvr ~ A_0^2 + B_inf * q * r^{-gamma}
                gvv ~ A_0^2 * 2M/r + B_inf * q^2

                For the determinant:
                det_2d = (A_0^2*2M/r + B_inf*q^2)(B_inf*r^{-2*gamma}) - (A_0^2 + B_inf*q*r^{-gamma})^2

                This is complex. Let me just test numerically.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid

def compute_metric_v4(r, params):
    """Metric with logarithmic psi (non-saturating) + mild A + damped B + q.

    psi(r) = phi_0 * ln(r/r_h) * S(r)  (logarithmic, psi' ~ 1/r)
    A = exp(beta_A * psi)  (mild growth with phi_0 < 1)
    B = B0 * |psi|^n / (1+|psi|^n) * exp(-psi^4/(2*sigma^4))  (quartic damped)
    phi = q*v + psi(r)
    """
    r = np.asarray(r, dtype=float)
    M = params['M']
    q = params['q']
    phi_0 = params['phi_0']
    delta = params['delta']
    beta_A = params.get('beta_A', -1.0)
    B0 = params['B0']
    n_B = params.get('n_B', 2.0)
    sigma_B = params.get('sigma_B', 1.5)

    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    F = 1.0 - 2.0 * M / r_safe

    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)
    psi = phi_0 * np.log(r_safe / r_h) * S
    dpsi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)

    A = np.exp(beta_A * psi)
    A2 = A ** 2

    abs_psi = np.abs(psi)
    B = B0 * abs_psi ** n_B / (1.0 + abs_psi ** n_B) * np.exp(-(psi ** 4) / (2.0 * sigma_B ** 4))

    gtilde_vv = -A2 * F + B * q ** 2
    gtilde_vr = A2 + B * q * dpsi
    gtilde_rr = B * dpsi ** 2
    gtilde_thth = A2 * r ** 2

    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2

    areal = A * r
    lorentzian = det_2d < 0

    return {
        'r': r, 'F': F, 'psi': psi, 'dpsi': dpsi,
        'A': A, 'B': B, 'A2': A2,
        'gtilde_vv': gtilde_vv, 'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr, 'gtilde_thth': gtilde_thth,
        'det_2d': det_2d, 'areal_radius': areal, 'lorentzian': lorentzian,
    }


def exact_timelike(r, metric, M=1.0):
    """Compute exact timelike geodesic using full metric."""
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    det_2d = metric['det_2d']

    # rdot^2 = -(gvv + 1) / det_2d  (for E=1, radial)
    # Inside: gvv > 0, det_2d < 0, so rdot^2 = (gvv+1)/|det_2d| > 0
    det_safe = np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    rdot_sq = -(gvv + 1.0) / det_safe
    rdot_sq = np.where(rdot_sq > 0, rdot_sq, np.nan)
    rdot = np.sqrt(rdot_sq)

    # dtau = dr / rdot
    idx_h = np.argmin(np.abs(r - 2.0 * M))
    r_int = r[:idx_h + 1]
    rdot_int = rdot[:idx_h + 1]
    rdot_rev = np.where(np.isfinite(rdot_int[::-1]) & (rdot_int[::-1] > 0),
                        rdot_int[::-1], 1e-30)
    dtau_dr = 1.0 / rdot_rev
    tau = cumulative_trapezoid(dtau_dr, r_int[::-1], initial=0)
    tau_total = float(tau[-1]) if len(tau) > 0 else 0.0

    # Power law
    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    rdot_inner = rdot[:n_check]
    mask = (rdot_inner > 0) & (r_inner > 0) & np.isfinite(rdot_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_rdot = np.log(rdot_inner[mask])
        alpha = np.polyfit(log_r, log_rdot, 1)[0]
        # dtau/dr ~ r^{-alpha}, integral diverges if alpha >= 1
        diverges = alpha >= 1.0
    else:
        alpha = np.nan
        diverges = False

    return {
        'tau_total': tau_total,
        'rdot_power_law': float(alpha) if np.isfinite(alpha) else None,
        'timelike_diverges': bool(diverges),
        'rdot': rdot,
    }


def exact_null(r, metric, M=1.0):
    """Compute exact null geodesic."""
    det_2d = metric['det_2d']
    disc = -det_2d
    sqrt_disc = np.sqrt(np.where(disc > 0, disc, np.nan))

    idx_h = np.argmin(np.abs(r - 2.0 * M))
    r_int = r[:idx_h + 1]
    null_int = sqrt_disc[:idx_h + 1]
    r_rev = r_int[::-1]
    int_rev = np.where(np.isfinite(null_int[::-1]), null_int[::-1], 0)
    lam = cumulative_trapezoid(int_rev, r_rev, initial=0)
    lam_total = float(lam[-1]) if len(lam) > 0 else 0.0

    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    int_inner = sqrt_disc[:n_check]
    mask = (int_inner > 0) & (r_inner > 0) & np.isfinite(int_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_int = np.log(int_inner[mask])
        alpha = -np.polyfit(log_r, log_int, 1)[0]
        diverges = alpha >= 1.0
    else:
        alpha = np.nan
        diverges = False

    return {
        'lambda_total': lam_total,
        'null_power_law': float(alpha) if np.isfinite(alpha) else None,
        'null_diverges': bool(diverges),
    }


def curvature_estimate(r, metric, M=1.0):
    """Estimate curvature using conformal formula (valid where B->0)."""
    A = metric['A']
    K_schw = 48.0 * M ** 2 / r ** 6
    K_conf = np.where(A > 1e-30, A ** (-12) * K_schw, np.inf)

    # Power law
    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    K_inner = K_conf[:n_check]
    mask = (K_inner > 0) & (r_inner > 0) & np.isfinite(K_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_K = np.log(K_inner[mask])
        alpha = np.polyfit(log_r, log_K, 1)[0]
    else:
        alpha = np.nan

    return {
        'K_max': float(np.nanmax(K_conf[np.isfinite(K_conf)])) if np.any(np.isfinite(K_conf)) else np.inf,
        'K_power_law': float(alpha) if np.isfinite(alpha) else None,
        'K_vanishes': bool(alpha > 0) if np.isfinite(alpha) else False,
    }


def test_v4():
    r = np.logspace(np.log10(1e-8), np.log10(50.0), 20000)

    print("=" * 80)
    print("LOGARITHMIC psi + MILD A + QUARTIC-DAMPED B + q != 0")
    print("Testing phi_0 = 0.9 (A ~ r^{-0.9}, areal ~ r^{0.1} -> 0)")
    print("=" * 80)

    for q in [0.0, 0.1, 1.0, 10.0, 100.0]:
        for B0 in [0.01, 0.1, 1.0]:
            params = {
                'M': 1.0, 'q': q, 'phi_0': 0.9, 'delta': 0.05,
                'beta_A': -1.0, 'B0': B0, 'n_B': 2.0, 'sigma_B': 1.5,
            }
            metric = compute_metric_v4(r, params)
            null = exact_null(r, metric, M=1.0)
            tl = exact_timelike(r, metric, M=1.0)
            curv = curvature_estimate(r, metric, M=1.0)

            det = metric['det_2d']
            frac_lor = np.sum(det < 0) / len(det)
            sc = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]
            areal = metric['areal_radius']
            areal_max = np.max(areal)
            areal_div = areal_max > 1e6

            print(f"\n  q={q:6.1f}, B0={B0:6.2f}: "
                  f"lor={frac_lor:.6f}, sc={len(sc)}, "
                  f"areal_max={areal_max:.2e}, areal_div={areal_div}, "
                  f"K_vanish={curv['K_vanishes']}, K_pwr={curv['K_power_law']}, "
                  f"null_div={null['null_diverges']}, null_pwr={null['null_power_law']}, "
                  f"tau_div={tl['timelike_diverges']}, tau_pwr={tl['rdot_power_law']}")

            if len(sc) > 0:
                print(f"    SIGN CHANGE at r={r[sc[0]]:.4e}")

    # Also test phi_0 = 1.0 (areal -> const, not diverging)
    print("\n" + "=" * 80)
    print("phi_0 = 1.0 (A ~ r^{-1}, areal = A*r -> const)")
    print("=" * 80)

    for q in [0.0, 0.1, 1.0, 10.0, 100.0]:
        for B0 in [0.01, 0.1, 1.0]:
            params = {
                'M': 1.0, 'q': q, 'phi_0': 1.0, 'delta': 0.05,
                'beta_A': -1.0, 'B0': B0, 'n_B': 2.0, 'sigma_B': 1.5,
            }
            metric = compute_metric_v4(r, params)
            null = exact_null(r, metric, M=1.0)
            tl = exact_timelike(r, metric, M=1.0)
            curv = curvature_estimate(r, metric, M=1.0)

            det = metric['det_2d']
            frac_lor = np.sum(det < 0) / len(det)
            sc = np.where(np.signbit(det[:-1]) != np.signbit(det[1:]))[0]
            areal = metric['areal_radius']
            areal_max = np.max(areal)
            areal_div = areal_max > 1e6

            print(f"\n  q={q:6.1f}, B0={B0:6.2f}: "
                  f"lor={frac_lor:.6f}, sc={len(sc)}, "
                  f"areal_max={areal_max:.2e}, areal_div={areal_div}, "
                  f"K_vanish={curv['K_vanishes']}, K_pwr={curv['K_power_law']}, "
                  f"null_div={null['null_diverges']}, null_pwr={null['null_power_law']}, "
                  f"tau_div={tl['timelike_diverges']}, tau_pwr={tl['rdot_power_law']}")

            if len(sc) > 0:
                print(f"    SIGN CHANGE at r={r[sc[0]]:.4e}")


if __name__ == '__main__':
    test_v4()
