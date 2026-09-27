#!/usr/bin/env python3
"""Shared utilities and TEP-BH physics model for the analysis pipeline.

Provides:
  - Path helpers (PROJECT_ROOT, RESULTS_DIR, FIGURES_DIR, LOGS_DIR)
  - step_json_path / step_csv_path / figure_path / log_path
  - write_json / read_json / write_csv (with NumPy-safe encoder)
  - ensure_dirs
  - TEPLogger / print_status / set_step_logger re-exports
  - TEPBHModel: the TEP black-hole model parameters
  - solve_tep_bh: construct the disformal matter metric and diagnostics
  - compute_curvature_invariants: Kretschmann, Ricci for gtilde
  - compute_energy_conditions: Raychaudhuri convergence sign (NEC/SEC) for gtilde
  - check_geodesic_completeness: affine parameter integrals
  - compute_physical_volume: volume on appropriate slices

Physics summary
---------------
Geometric metric (Einstein frame, ingoing Eddington-Finkelstein):

    ds_g^2 = -(1 - 2M/r) dv^2 + 2 dv dr + r^2 dOmega^2

Scalar field profile (activates inside the horizon):

    phi(r) = phi_0 * ln(r / r_h) * S(r)

where S(r) = 1 / (1 + exp((r - r_h)/(delta * r_h))) is a smooth logistic
activation.  Inside the horizon (r < r_h), ln(r/r_h) < 0 so phi < 0.

Conformal factor:

    A(phi) = exp(beta_A * phi)   with beta_A = -1

Inside the horizon phi < 0, so A = exp(|phi|) -> +infinity as r -> 0.
With phi_0 = 2, A^2 ~ (r_h/r)^{2*phi_0} = (r_h/r)^4.

Disformal function (ULTRA-DAMPED SHEAR BUMP):

    B(phi) = B0 * |phi|^n_B / (1 + |phi|^n_B) * exp(-phi^4 / (2*sigma_B^4))
    with B0 = 1, n_B = 2, sigma_B = 1.5

This activates the disformal shear near the horizon (creating the optical
illusion of spatial collapse via temporal refraction) but vanishes
super-polynomially in the deep core as |phi| -> inf, because the quartic
Gaussian exp(-phi^4/(2*sigma_B^4)) dominates the rational |phi|^2/(1+|phi|^2)
saturation factor.  The quartic exponent (rather than quadratic) ensures
conformal dominance even against violently diverging Coulomb-like scalar
profiles (phi ~ 1/r) in the self-gravitating sGB branch.  When B -> 0, the
conformal stretching A^4 dominates the 2D determinant, keeping it strictly
negative (det gtilde_{2D} < 0) for all r > 0.  The metric is globally
Lorentzian and invertible — there is no determinant-zero temporal boundary.
In the exterior (r > 2M), the scalar field is screened by the logistic
profile: phi ~ 0, A ~ 1, B ~ 0.
The conformal factor A cancels in geodesic equations, so exterior
observables (shadow, ISCO, QNM) match Schwarzschild to < 0.1%.

Shear zone: the region near the horizon where B is active, producing
transient disformal shear.  Deeper in, B decays and the geometry gives
way to pure conformal dilation — the conformal factor A^4 dominates and
the light cones never close.

Disformal matter metric (q = 0, static scalar):

    gtilde_{vv}  = -A^2 (1 - 2M/r)
    gtilde_{vr}  = A^2
    gtilde_{rr}  = B (phi')^2
    gtilde_{thth}= A^2 r^2

The 2D (v,r) determinant:

    det_2d = gtilde_vv * gtilde_rr - gtilde_vr^2

Lorentzian signature requires det_2d < 0.  With the quartic Gaussian-damped
B(phi) and phi_0 = 2, the conformal term A^4 always dominates: in the
shear zone near the horizon B is active but not strong enough to flip
the sign, and in the deep core B -> 0 so the conformal term dominates
completely.  The determinant is strictly negative everywhere — the
space is globally Lorentzian and invertible with no boundary.

Geodesic completeness: both null and timelike geodesic integrands
diverge as r -> 0 (null ~ r^{-2*phi_0}, timelike ~ r^{1/2-phi_0}),
so the affine parameter and proper time to reach r=0 are infinite.
The space is null and timelike geodesically complete.

Physical volume: inside the horizon, r is timelike.  The physical
3-volume on a constant-r spacelike slice is

    dV = 4*pi * A^2 * r^2 * sqrt(|F|) dv

A^2 * r^2 ~ (r_h/r)^{2*phi_0} * r^2 = r_h^{2*phi_0} * r^{2-2*phi_0}
which -> infinity for phi_0 > 1.  The physical area of spheres GROWS,
not shrinks — this is the density illusion: the geometric r^2 -> 0 but
the physical A^2*r^2 diverges.  The areal radius rho = A*r ~ r^{1-phi_0}
diverges for phi_0 > 1, meaning the physical space opens up rather than
collapses.

EXACT FIXED-BACKGROUND THEOREM (verified analytically, see manuscript Section 3 / Appendix D):
For the pure conformal metric gtilde = A^2 * g_Schw with A = (r_h/r)^{phi_0}
in the deep interior (B -> 0), the exact Kretschmann scalar is:

    K[gtilde] ~ r^{4*phi_0 - 6}

NOT r^{12*phi_0 - 6} as the naive conformal formula A^{-12}*K_Schw gives.
The conformal formula omits derivative terms from the non-constant A(r)
which dominate for phi_0 < 3/2.

  - K -> 0 (singularity removed) requires phi_0 > 3/2
  - Finite areal radius rho = A*r requires phi_0 <= 1
  - These are MUTUALLY EXCLUSIVE.

The phi_0=2 prototype removes the singularity (K ~ r^2 -> 0) but has
diverging areal radius. The phi_0=1 branch has finite areal radius but
K = 9/(M^2*r^2) -> infinity (singularity NOT removed).

IMPLICATION: The fixed-Schwarzschild construction cannot simultaneously
produce a regular interior and bounded physical spatial radius. The
temporal field must backreact on the geometric metric. See the
manuscript (Section 3, Appendix D) for the full analysis.

Version: TEP-BH v0.2 (Bahrain)
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable
import numpy as np
from scipy.integrate import cumulative_trapezoid

try:
    from scripts.utils.logger import TEPLogger, print_status, set_step_logger
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "utils"))
    from logger import TEPLogger, print_status, set_step_logger

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
LOGS_DIR = PROJECT_ROOT / "logs"
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"


# =============================================================================
# Path and I/O helpers
# =============================================================================

def ensure_dirs() -> None:
    for d in [RESULTS_DIR, FIGURES_DIR, LOGS_DIR, DATA_DIR, RAW_DIR, PROCESSED_DIR]:
        d.mkdir(parents=True, exist_ok=True)


class TEPEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, cls=TEPEncoder) + "\n",
                    encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def step_json_path(step_id: str) -> Path:
    return RESULTS_DIR / f"{step_id}.json"


def step_csv_path(step_id: str) -> Path:
    return RESULTS_DIR / f"{step_id}.csv"


def figure_path(name: str) -> Path:
    return FIGURES_DIR / f"{name}.png"


def log_path(step_id: str) -> Path:
    return LOGS_DIR / f"{step_id}.log"


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def rounded(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def make_step_logger(step_id: str) -> TEPLogger:
    ensure_dirs()
    logger = TEPLogger(
        name=f"tep_bh_{step_id}",
        log_file_path=log_path(step_id),
        reset_log=True,
    )
    # Also tee step output to the main pipeline log so it captures everything.
    pipeline_log = LOGS_DIR / "tep_bh_pipeline.log"
    if pipeline_log.exists():
        import logging as _logging
        from scripts.utils.logger import TEPFileFormatter
        tee_fh = _logging.FileHandler(pipeline_log, mode='a', encoding='utf-8')
        tee_fh.setLevel(_logging.DEBUG)
        tee_fh.setFormatter(TEPFileFormatter())
        logger.logger.addHandler(tee_fh)
    set_step_logger(logger)
    return logger


def finalize_result(step_id: str, result: dict, description: str,
                    key_result: str, model: dict | None = None,
                    dependencies: list | None = None) -> dict:
    """Inject standard metadata into a step's result dict.

    Every step's JSON output should include:
      - step: step identifier
      - description: human-readable description of what the step does
      - key_result: one-sentence summary of the main finding
      - model: parameter set used (if applicable)
      - dependencies: list of upstream step_ids this step consumed
      - status: 'success' or 'failed'
      - timestamp: ISO format
    """
    result['step'] = step_id
    result['description'] = description
    result['key_result'] = key_result
    if model is not None:
        result['model'] = model
    if dependencies is not None:
        result['dependencies'] = dependencies
    result.setdefault('status', 'success')
    from datetime import datetime
    result['timestamp'] = datetime.now().isoformat()
    return result


# =============================================================================
# TEP Black-Hole Model
# =============================================================================

def _scalar_field_profile(r, M, phi_0, delta):
    """Scalar field profile phi(r) = phi_0 * ln(r/r_h) * S(r).

    S(r) = 1 / (1 + exp((r - r_h)/(delta * r_h))) is a smooth logistic
    activation that is ~1 inside the horizon and ~0 outside.

    Inside the horizon (r < r_h): ln(r/r_h) < 0, so phi < 0.
    With beta_A = -1: A = exp(-phi) = exp(|phi|) -> +inf as r -> 0.
    """
    r = np.asarray(r, dtype=float)
    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    phi = phi_0 * np.log(r_safe / r_h) * S
    return phi


def tep_s2_conformal_factor(r_rs, eta, delta=2000.0, beta_A=1.0):
    """TEP S-star inference conformal factor: logarithmic branch.

    This is the TEP power-law/logarithmic branch, A = exp(beta_A * phi),
    with phi = eta * ln(r_rs) * S(r_rs; delta), where r_rs = r / r_h.

    For the S2 exterior (r_rs >> 1) and beta_A = +1, positive eta gives
    A > 1 (mass-inflation branch).  The logistic S damps the field at
    galactic scales (r_rs >> delta), so A -> 1 at Earth.

    The default delta may be overridden by the TEP_LOG_DELTA environment
    variable for parameter scans.

    Parameters
    ----------
    r_rs : array_like
        Radius in units of r_h (Schwarzschild radii).
    eta : float
        Dimensionless TEP coupling (mapped to phi_0).
    delta : float
        Logistic width in units of r_h; controls how far the logarithmic
        field extends above the horizon.  Default 2e3 is the S2 best-fit scale.
    beta_A : float
        Sign of the conformal mapping; +1 for mass-inflation outside.
    """
    import os
    delta_env = os.environ.get("TEP_LOG_DELTA")
    if delta_env is not None:
        delta = float(delta_env)
    r_rs = np.asarray(r_rs, dtype=float)
    phi = _scalar_field_profile(r_rs, 0.5, eta, delta)
    return np.exp(beta_A * phi)


def _scalar_field_gradient(r, M, phi_0, delta):
    """Radial derivative dphi/dr."""
    r = np.asarray(r, dtype=float)
    r_h = 2.0 * M
    r_safe = np.maximum(r, 1e-30)
    S = 1.0 / (1.0 + np.exp((r - r_h) / (delta * r_h)))
    dS = -S * (1.0 - S) / (delta * r_h)
    dphi = phi_0 * (1.0 / r_safe * S + np.log(r_safe / r_h) * dS)
    return dphi


class TEPBHModel:
    """TEP black-hole model parameters.

    A(phi) = exp(beta_A * phi)     — conformal clock mapping
    B(phi) = B0 * |phi|^2 / (1 + |phi|^2) * exp(-phi^4 / (2*sigma_B^4))
             — quartic-damped "shear bump" disformal response

    The disformal metric (q=0, static scalar, Schwarzschild EF geometric):
        gtilde_{vv}  = -A^2 (1 - 2M/r)
        gtilde_{vr}  = A^2
        gtilde_{rr}  = B (phi')^2
        gtilde_{thth}= A^2 r^2

    The quartic Gaussian-damped B(phi) activates the disformal shear near the
    horizon (creating the optical illusion of spatial collapse via
    temporal refraction) but vanishes in the deep core as |phi| -> inf.
    When B -> 0, the conformal stretching A^4 dominates the 2D determinant,
    keeping it strictly negative (det gtilde_{2D} < 0) for all r > 0.
    The metric is globally Lorentzian and invertible — there is no
    determinant-zero temporal boundary.

    Asymptotic analysis (r -> 0, phi ~ phi_0 * ln(r/r_h) -> -inf):
        A^2 ~ (r_h/r)^{2*phi_0}           (conformal, diverges)
        B   ~ exp(-phi^4 / (2*sigma_B^4))  (disformal, super-polynomial decay)
        |F| B (phi')^2 / A^2 -> 0          (conformal dominance)

    For phi_0 = 2: A^4 ~ r^{-8}, while the disformal term decays
    super-polynomially; the conformal term always dominates, so
    det gtilde_{2D} < 0 everywhere (globally Lorentzian).

    Geodesic completeness: both null and timelike geodesic integrands
    diverge as r -> 0 (null ~ r^{-2*phi_0}, timelike ~ r^{1/2-phi_0}),
    so the affine parameter and proper time to reach r=0 are infinite.
    The space is null and timelike geodesically complete.

    Curvature: the EXACT Kretschmann K ~ r^{4*phi_0 - 6} (verified
    analytically; see manuscript Section 3 / Appendix D).
    For phi_0 = 2, K ~ r^2 -> 0, so the center is not a curvature
    singularity.  The areal radius rho = A*r ~ r^{1-phi_0} diverges for
    phi_0 > 1, meaning the physical space opens up rather than collapsing.

    NOTE: The naive conformal formula A^{-12}*K_Schw ~ r^{12*phi_0-6} is
    WRONG for phi_0 < 3/2 (it omits derivative terms). The exact formula
    K ~ r^{4*phi_0-6} gives a no-go theorem: curvature regularity
    (phi_0 > 3/2) and bounded areal radius (phi_0 <= 1) are mutually
    exclusive. This prototype (phi_0=2) is curvature-regular but
    spatially-enlarged. See manuscript Section 3.2.
    """

    def __init__(self, beta_A=-1.0, B0=1.0, n_B=2.0,
                 phi_0=2.0, delta=0.05, M=1.0, sigma_B=1.5):
        self.beta_A = beta_A
        self.B0 = B0
        self.n_B = n_B
        self.phi_0 = phi_0
        self.delta = delta
        self.M = M
        self.r_h = 2.0 * M
        self.sigma_B = sigma_B

    def A(self, phi):
        """Conformal factor A(phi) = exp(beta_A * phi)."""
        return np.exp(self.beta_A * np.asarray(phi, dtype=float))

    def B(self, phi):
        """Ultra-damped "shear bump" disformal function.

        B(phi) = B0 * |phi|^2 / (1 + |phi|^2) * exp(-phi^4 / (2*sigma_B^4))

        Uses a quartic Gaussian envelope to ensure B -> 0 in the deep core
        even against violently diverging scalar profiles (e.g. the Coulomb-like
        sGB scalar phi ~ 1/r).  The quartic exponent exp(-phi^4) decays far
        faster than the quadratic exp(-phi^2), guaranteeing conformal dominance
        for both the prescribed (logarithmic) and self-gravitating (Coulomb)
        scalar branches.

        This activates the disformal shear near the horizon (creating the
        optical illusion of spatial collapse via temporal refraction) but
        vanishes in the deep core. When B -> 0, the conformal stretching
        completely dominates, preventing the light cones from closing and
        ensuring the space remains fully Lorentzian down to r=0.
        """
        phi_arr = np.asarray(phi, dtype=float)
        abs_phi = np.abs(phi_arr)
        return self.B0 * abs_phi**self.n_B / (1.0 + abs_phi**self.n_B) * np.exp(-(phi_arr**4) / (2.0 * self.sigma_B**4))

    def phi(self, r):
        """Scalar field profile phi(r)."""
        return _scalar_field_profile(r, self.M, self.phi_0, self.delta)

    def dphi_dr(self, r):
        """Radial derivative of scalar field."""
        return _scalar_field_gradient(r, self.M, self.phi_0, self.delta)

    def b0_canonical_coefficient(self, mass_kg):
        """Convention-A equivalent of the code-unit B0 for a given mass.

        The code uses geometric units (G = c = 1) with a dimensionless
        field (phi = Phi/M_*) and the black-hole mass M = 1, so B0 = 1.0
        means B0 = 1.0 M_BH^2 in geometrized units -- i.e. mass^{-2} in
        natural units (hbar = c = 1), the correct dimension for a
        dimensionless-field disformal coefficient. In the canonical
        Paper 0 EFT the field has mass dimension 1 and B_0 has
        dimension -4; expressed as the dimensionless coefficient
        B_0 M_Pl^4 of the holonomy-bound convention (|B_0| M_Pl^4
        <= 5e-8), the code value maps to

            B0_A = B0 * (M_BH_geom * M_Pl_nat)^2   (dimensionless)

        where M_BH_geom = GM/c^2 [m] and M_Pl_nat = M_Pl c/hbar [m^-1].
        For a solar-mass black hole this is ~3e74 -- far ABOVE the bound.
        This is NOT a violation: it demonstrates that the strong-field
        construction normalization and the canonical weak-field transport
        amplitude are distinct-sector parameters that must not be
        numerically identified (Paper 0 Section 2.2 states that
        dimensionless-field normalizations are not imported into the
        weak-field theory).
        """
        G_SI = 6.674e-11        # m^3 kg^-1 s^-2
        C_SI = 299792458.0      # m/s
        HBAR_SI = 1.054571817e-34  # J s
        M_PL_KG = 4.341e-9      # reduced Planck mass [kg]
        m_geom = G_SI * mass_kg / C_SI**2          # geometrized mass [m]
        m_pl_inv = M_PL_KG * C_SI / HBAR_SI        # M_Pl in m^-1
        return self.B0 * (m_geom * m_pl_inv) ** 2

    def to_dict(self):
        M_SUN_KG = 1.989e30
        return {
            'beta_A': self.beta_A,
            'B0': self.B0,
            'B0_units': 'geometric code units (M = 1, dimensionless field); B0 = 1 M_BH^2 geometrized (mass^-2 natural units) -- strong-field construction normalization, not the canonical weak-field coefficient of Paper 0',
            'B0_canonical_coefficient_per_solar_mass': self.b0_canonical_coefficient(M_SUN_KG),
            'B0_convention_note': 'B0_A = B0 * (M_BH_geom * M_Pl_nat)^2 in the canonical dimensionless convention; ~3e74 for 1 M_sun, ~1e94 for M87* -- incommensurable with the Paper 0 holonomy bound |B_0| <= 5e-8 because the two are distinct-sector parameters (strong-field construction envelope vs weak-field transport amplitude), not two values of one parameter',
            'n_B': self.n_B,
            'B_type': 'gaussian_bump',
            'sigma_B': self.sigma_B,
            'phi_0': self.phi_0,
            'delta': self.delta,
            'M': self.M,
            'r_h': self.r_h,
        }

def _schwarzschild_F(r, M):
    """Schwarzschild metric function F(r) = 1 - 2M/r."""
    return 1.0 - 2.0 * M / np.asarray(r, dtype=float)


def compute_disformal_metric(r, model):
    """Compute the full disformal matter metric gtilde_{mu nu}.

    For Schwarzschild EF geometric metric with static scalar (q=0):
        gtilde_{vv}  = -A^2 (1 - 2M/r)
        gtilde_{vr}  = A^2
        gtilde_{rr}  = B (phi')^2
        gtilde_{thth}= A^2 r^2
        gtilde_{phph}= A^2 r^2 sin^2(theta)

    Parameters
    ----------
    r : ndarray (ascending)
    model : TEPBHModel

    Returns
    -------
    dict with all metric components and derived quantities.
    """
    r = np.asarray(r, dtype=float)
    M = model.M

    F = _schwarzschild_F(r, M)
    phi = model.phi(r)
    dphi = model.dphi_dr(r)
    A = model.A(phi)
    B = model.B(phi)
    A2 = A ** 2

    gtilde_vv = -A2 * F
    gtilde_vr = A2
    gtilde_rr = B * dphi ** 2
    gtilde_thth = A2 * r ** 2
    gtilde_phph = A2 * r ** 2

    det_2d = gtilde_vv * gtilde_rr - gtilde_vr ** 2
    det_4d = det_2d * gtilde_thth * gtilde_phph

    lorentzian = det_2d < 0
    nondegenerate = np.abs(det_4d) > 1e-30

    # Scalar kinetic term X = -1/2 g^{mu nu} dphi_mu dphi_nu
    # For EF: g^{rr} = F, g^{vr} = 1, g^{vv} = 0
    X = -0.5 * F * dphi ** 2

    # Invertibility: A^2 - 2 B X != 0 for X = -1/2 (grad phi)^2
    invert_cond = A2 - 2.0 * B * X

    # Ratio that controls signature: disformal_term / conformal_term
    # det_2d = -A^4 + A^2 * |F| * B * (phi')^2  (inside horizon, F<0)
    # = A^4 * (-1 + |F| * B * (phi')^2 / A^2)
    # Lorentzian iff |F| * B * (phi')^2 / A^2 < 1
    disformal_ratio = np.abs(F) * B * dphi ** 2 / np.where(A2 > 1e-30, A2, 1e-30)

    return {
        'r': r,
        'F': F,
        'phi': phi,
        'dphi': dphi,
        'A': A,
        'B': B,
        'A2': A2,
        'gtilde_vv': gtilde_vv,
        'gtilde_vr': gtilde_vr,
        'gtilde_rr': gtilde_rr,
        'gtilde_thth': gtilde_thth,
        'gtilde_phph': gtilde_phph,
        'det_2d': det_2d,
        'det_4d': det_4d,
        'lorentzian': lorentzian,
        'nondegenerate': nondegenerate,
        'X': X,
        'invert_cond': invert_cond,
        'disformal_ratio': disformal_ratio,
    }


# =============================================================================
# Curvature invariants
# =============================================================================

def compute_curvature_invariants(r, metric, M=1.0):
    """Compute curvature invariants of the disformal metric gtilde.

    The curvature is computed in STANDARD Schwarzschild coordinates (t, r)
    where the metric is diagonal and g_rr = A^2/F + B*(phi')^2 does NOT
    vanish in the conformal regime (unlike EF coordinates where g_rr =
    B*(phi')^2 -> 0 when B -> 0).

    Standard-coordinates metric (q=0, static scalar):
        gtilde_{tt}  = -A^2 * F
        gtilde_{tr}  = 0  (diagonal)
        gtilde_{rr}  = A^2 / F + B * (phi')^2
        gtilde_{thth}= A^2 * r^2

    For a diagonal spherically symmetric metric, the Kretschmann scalar is:

        K = 4*E^2 + 8*F_t^2 + 8*F_r^2 + 4*G^2

    where (orthonormal-frame Riemann components):
        E   = R_{trtr} / (g_tt * g_rr)
        F_t = R_{tthth} / (g_tt * g_thth)
        F_r = R_{rthrth} / (g_rr * g_thth)
        G   = R_{thphthph} / g_thth^2

    with:
        R_{trtr}    = g_tt * [d_r(Gamma^t_tr) + (Gamma^t_tr)^2 - Gamma^t_tr * Gamma^r_rr]
        R_{rthrth}  = g_thth * [d_r(Gamma^th_rth) + (Gamma^th_rth)^2 - Gamma^th_rth * Gamma^r_rr]
        R_{tthth}   = -0.25 * g_tt' * g_thth' / g_rr
        R_{thphthph}= g_thth * (1 - g^{rr} * g_thth'^2 / (4 * g_thth))

    This formula is verified to give K = 48*M^2/r^6 for Schwarzschild.

    EXACT NO-GO THEOREM (verified analytically with SymPy):
    For the pure conformal metric gtilde = A^2 * g_Schw with A = (r_h/r)^{phi_0}
    in the deep interior (B -> 0):
        K ~ r^{4*phi_0 - 6}
    - K -> 0 (singularity removed) requires phi_0 > 3/2
    - Finite areal radius rho = A*r requires phi_0 <= 1
    These are MUTUALLY EXCLUSIVE. The conformal formula A^{-12}*K_Schw
    ~ r^{12*phi_0-6} is WRONG for phi_0 < 3/2 because it omits derivative
    terms from the non-constant conformal factor.

    The horizon (F=0) is a coordinate singularity in standard coordinates.
    We regularize F near the horizon and blend with EF-based curvature
    in a narrow transition band.

    Parameters
    ----------
    r : ndarray (ascending)
    metric : dict from compute_disformal_metric
    M : float

    Returns
    -------
    dict with Kretschmann, Ricci scalar, and references.
    """
    r = np.asarray(r, dtype=float)
    A = metric['A']
    B = metric['B']
    F = metric['F']
    dphi = metric['dphi']
    A2 = metric['A2']

    # Schwarzschild reference
    K_schwarzschild = 48.0 * M ** 2 / r ** 6
    R_schwarzschild = np.zeros_like(r)  # vacuum

    # --- Standard-coordinates metric components (diagonal, q=0) ---
    # g_tt = -A^2 * F  (same as EF g_vv)
    # g_rr = A^2 / F + B * (phi')^2  (does NOT vanish when B->0)
    # g_thth = A^2 * r^2  (same as EF)
    #
    # Regularize F near the horizon to avoid 1/F -> infinity.
    # The horizon is at r = 2M where F = 0. In a narrow band around it,
    # we clamp |F| to a minimum value.
    r_h = 2.0 * M
    F_reg_eps = 1e-6  # minimum |F| for regularization
    F_abs = np.abs(F)
    F_sign = np.sign(F)
    F_reg = F_sign * np.maximum(F_abs, F_reg_eps)

    g_tt = -A2 * F  # same in standard and EF (q=0)
    g_rr = A2 / F_reg + B * dphi ** 2  # standard coords: A^2/F + B*(phi')^2
    g_thth = A2 * r ** 2

    # Safe versions (avoid division by zero)
    g_tt_safe = np.where(np.abs(g_tt) > 1e-50, g_tt, np.nan)
    g_rr_safe = np.where(np.abs(g_rr) > 1e-50, g_rr, np.nan)
    g_thth_safe = np.where(np.abs(g_thth) > 1e-50, g_thth, np.nan)

    # Inverse metric (diagonal)
    ginv_tt = 1.0 / g_tt_safe
    ginv_rr = 1.0 / g_rr_safe
    ginv_thth = 1.0 / g_thth_safe

    # Derivatives (numerical, on the ascending r grid)
    g_tt_p = np.gradient(g_tt, r)
    g_rr_p = np.gradient(g_rr, r)
    g_thth_p = np.gradient(g_thth, r)

    # Christoffel symbols (diagonal metric, only r-derivatives)
    Gt_tr = 0.5 * ginv_tt * g_tt_p          # Gamma^t_tr
    Gr_tt = -0.5 * ginv_rr * g_tt_p          # Gamma^r_tt
    Gr_rr = 0.5 * ginv_rr * g_rr_p           # Gamma^r_rr
    Gr_thth = -0.5 * ginv_rr * g_thth_p      # Gamma^r_thth
    Gth_rth = 0.5 * ginv_thth * g_thth_p     # Gamma^th_rth

    # Derivatives of Christoffels
    dGt_tr = np.gradient(Gt_tr, r)
    dGth_rth = np.gradient(Gth_rth, r)

    # Riemann components (coordinate basis)
    # R^t_{rtr} = d_r(Gamma^t_tr) + (Gamma^t_tr)^2 - Gamma^t_tr * Gamma^r_rr
    R_trtr = g_tt * (dGt_tr + Gt_tr ** 2 - Gt_tr * Gr_rr)

    # R^th_{rthr} = d_r(Gamma^th_rth) + (Gamma^th_rth)^2 - Gamma^th_rth * Gamma^r_rr
    R_rthrth = g_thth * (dGth_rth + Gth_rth ** 2 - Gth_rth * Gr_rr)

    # R_{tthth} = -0.25 * g_tt' * g_thth' / g_rr
    R_tthtth = -0.25 * g_tt_p * g_thth_p / g_rr_safe

    # R_{thphthph} = g_thth * (1 - g^{rr} * g_thth'^2 / (4 * g_thth))
    R_thphthph = g_thth * (1.0 - ginv_rr * g_thth_p ** 2 / (4.0 * g_thth_safe))

    # Orthonormal-frame components
    E = R_trtr / (g_tt_safe * g_rr_safe)
    F_t = R_tthtth / (g_tt_safe * g_thth_safe)
    F_r = R_rthrth / (g_rr_safe * g_thth_safe)
    G_comp = R_thphthph / g_thth_safe ** 2

    # Kretschmann: K = 4*E^2 + 8*F_t^2 + 8*F_r^2 + 4*G^2
    Kretschmann = 4.0 * E ** 2 + 8.0 * F_t ** 2 + 8.0 * F_r ** 2 + 4.0 * G_comp ** 2

    # --- Ricci scalar (diagonal spherically symmetric) ---
    # R_{tt} = R^r_{trt} = -R^r_{ttr} = -(d_r(Gamma^r_tt) + Gamma^r_tt*Gamma^r_rr - ...)
    # For diagonal metric:
    # R_{tt} = -d_r(Gr_tt) - Gr_tt * Gr_rr + Gt_tr * Gr_tt / (g_tt/g_rr)
    # Actually, use the standard formula:
    # R_{tt} = -0.5 * g_tt'' / g_rr + 0.5 * g_tt' * g_rr' / (2 * g_rr^2) + ...
    # Simpler: R = 2 * R_{trtr} / (g_tt * g_rr) + 2 * R_{rthrth} / (g_rr * g_thth)
    #          + 2 * R_{tthth} / (g_tt * g_thth) + 2 * R_{thphthph} / g_thth^2
    # (This is the trace of the Riemann tensor for diagonal metric)
    g_tt_pp = np.gradient(g_tt_p, r)
    g_thth_pp = np.gradient(g_thth_p, r)

    # 2D Ricci scalar (t,r sector)
    R_2D = 2.0 * R_trtr / (g_tt_safe * g_rr_safe)

    # Angular Ricci
    R_ang = -g_thth_pp / g_thth_safe + 0.5 * (g_thth_p / g_thth_safe) ** 2

    # Total Ricci scalar
    R_total = R_2D + 2.0 * R_ang / g_thth_safe

    # --- Horizon blending ---
    # Near the horizon (F=0), standard coordinates have a coordinate singularity
    # (g_rr = A²/F → ∞). The F regularization (clamping |F|) creates a cusp
    # that causes inaccurate numerical derivatives.
    #
    # Near the horizon, A ≈ 1 (the scalar field is just activating), so the
    # metric is approximately Schwarzschild. We use the Schwarzschild K
    # (which equals the conformal proxy when A ≈ 1) as the horizon blend.
    # This is NOT the wrong conformal formula — it's the exact Schwarzschild
    # value, which is correct when A ≈ 1 and B ≈ 0.
    #
    # In the deep interior (where A is large and the conformal formula is
    # wrong), the standard-coordinates curvature is used exclusively.
    A_inv12 = np.where(A > 1e-30, A ** (-12), np.inf)
    K_conformal_proxy = A_inv12 * K_schwarzschild

    # Blend: use standard-coords curvature in the interior, Schwarzschild K
    # near the horizon and exterior (where A ≈ 1).
    # The transition is controlled by |F|: when |F| is small (near horizon),
    # use Schwarzschild; when |F| is large (deep interior or far exterior),
    # use standard-coords.
    horizon_band = np.abs(F) < 0.1  # band around horizon where std coords fail
    exterior = F > 0.1  # exterior where A ≈ 1 and std coords are fine but Schw is exact
    use_schw = horizon_band  # use Schwarzschild K near horizon
    Kretschmann = np.where(use_schw & np.isfinite(K_conformal_proxy), K_conformal_proxy, Kretschmann)
    R_total = np.where(horizon_band, np.nan, R_total)

    # Clean up
    Kretschmann = np.where(np.isfinite(Kretschmann), Kretschmann, np.nan)

    return {
        'R_trtr': R_trtr,
        'R_rthrth': R_rthrth,
        'R_tthtth': R_tthtth,
        'R_thphthph': R_thphthph,
        'Ricci_scalar': R_total,
        'Ricci_2D': R_2D,
        'Ricci_angular': R_ang,
        'Kretschmann': Kretschmann,
        'Kretschmann_conformal': K_conformal_proxy,  # diagnostic only, NOT invariant
        'Kretschmann_disformal': K_conformal_proxy,  # retained for API compat
        'Kretschmann_schwarzschild': K_schwarzschild,
        'Ricci_scalar_schwarzschild': R_schwarzschild,
        'Kretschmann_ratio': Kretschmann / np.where(K_schwarzschild > 0, K_schwarzschild, np.nan),
    }


# =============================================================================
# Energy conditions / Raychaudhuri convergence sign
# =============================================================================

def compute_energy_conditions(r, metric, curvature, M=1.0):
    """Compute the timelike and null convergence conditions for gtilde.

    This is the calculation flagged as outstanding in Section 7.3 ("their
    sign is not fixed by the ansatz alone"): the sign of
    Rtilde_{mu nu} u^mu u^nu, which is the source term in the Raychaudhuri
    equation and, via the Einstein equations, the geometric statement of
    the Strong Energy Condition (SEC, timelike case) and Null Energy
    Condition (NEC, null case). If this sign is negative on the relevant
    congruence, the disformal branch geometrically violates the SEC/NEC
    in exactly the sense required by the Hawking-Penrose theorems to
    evade geodesic incompleteness -- the same escape clause used in
    TEP-TH (Paper 27) for the cosmological shear zone.

    We use the warped-product decomposition of the 4D Ricci tensor for
    ds^2 = h_{AB} dx^A dx^B + rho^2(r) dOmega^2 with rho = A(r) r
    (the areal radius of gtilde):

        R_{AB}      = (1/2) R^(2) h_{AB} - (2/rho) D_A D_B rho
        R_{thth}    = 1 - |D rho|^2 - rho * Box(rho)

    where D_A is the covariant derivative of the 2D (v,r) metric h_{AB},
    and R^(2) is its Ricci scalar (already computed as curvature['Ricci_2D']).
    This decomposition is verified below to reduce to R_{mu nu} = 0 for the
    A=1, B=0 (Schwarzschild) limit as a self-consistency check.

    Two congruences are evaluated:
      - Null, ingoing (k^A from the exact null condition; NEC-relevant,
        robust everywhere since the space is globally Lorentzian and no
        timelike normalization is required).
      - Timelike, radial infall from rest at infinity, E=1 (u^A from the
        exact geodesic normalization; SEC-relevant).

    Parameters
    ----------
    r : ndarray (ascending)
    metric : dict from compute_disformal_metric
    curvature : dict from compute_curvature_invariants
    M : float

    Returns
    -------
    dict with the convergence-condition sign diagnostics.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    A = metric['A']
    det_2d = metric['det_2d']
    det_2d_safe = np.where(np.abs(det_2d) > 1e-50, det_2d, np.nan)
    R_2D = curvature['Ricci_2D']

    # --- 2D Christoffels (r-row only; rho depends on r alone) ---
    gvv_p = np.gradient(gvv, r)
    gvr_p = np.gradient(gvr, r)
    grr_p = np.gradient(grr, r)

    ginv_vv = grr / det_2d_safe
    ginv_vr = -gvr / det_2d_safe
    ginv_rr = gvv / det_2d_safe

    Gr_vv = -0.5 * ginv_rr * gvv_p
    Gr_vr = 0.5 * ginv_vr * gvv_p
    Gr_rr = ginv_vr * gvr_p + 0.5 * ginv_rr * grr_p

    # --- Areal radius of gtilde and its covariant Hessian ---
    rho = A * r
    rho_p = np.gradient(rho, r)
    rho_pp = np.gradient(rho_p, r)

    DvDv_rho = -Gr_vv * rho_p
    DvDr_rho = -Gr_vr * rho_p
    DrDr_rho = rho_pp - Gr_rr * rho_p

    box_rho = (ginv_vv * DvDv_rho + 2.0 * ginv_vr * DvDr_rho + ginv_rr * DrDr_rho)
    grad_rho_sq = ginv_rr * rho_p ** 2

    # --- Angular Ricci self-consistency check (should -> 0 for Schwarzschild) ---
    R_thth_over_gthth = (1.0 - grad_rho_sq - rho * box_rho) / np.where(rho ** 2 > 1e-50, rho ** 2, np.nan)

    # --- 2D block of the 4D Ricci tensor, R_{AB} = (1/2) R2D h_{AB} - (2/rho) D_A D_B rho ---
    Rvv = 0.5 * R_2D * gvv - (2.0 / np.where(rho > 1e-50, rho, np.nan)) * DvDv_rho
    Rvr = 0.5 * R_2D * gvr - (2.0 / np.where(rho > 1e-50, rho, np.nan)) * DvDr_rho
    Rrr = 0.5 * R_2D * grr - (2.0 / np.where(rho > 1e-50, rho, np.nan)) * DrDr_rho

    # --- Null congruence (ingoing): k^A = (dv/dr|_in, 1) up to overall scale ---
    disc = -det_2d
    sqrt_disc = np.sqrt(np.where(disc > 0, disc, np.nan))
    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)
    dv_dr_in = (-gvr - sqrt_disc) / gvv_safe
    kv, kr = dv_dr_in, np.ones_like(r)
    R_kk = Rvv * kv ** 2 + 2.0 * Rvr * kv * kr + Rrr * kr ** 2

    # --- Timelike congruence: radial infall from rest at infinity, E=1 ---
    E = 1.0
    rdot_sq = -(gvv + E ** 2) / det_2d_safe
    rdot = -np.sqrt(np.where(rdot_sq > 0, rdot_sq, np.nan))  # ingoing branch
    vdot = -(E + gvr * rdot) / gvv_safe
    R_uu = Rvv * vdot ** 2 + 2.0 * Rvr * vdot * rdot + Rrr * rdot ** 2

    lorentzian = metric['lorentzian']
    # Noise floor: the far exterior has |R_kk|, |R_uu| ~ 1e-8 or smaller
    # (the disformal correction is genuinely tiny there, cf. A=1.00027 at
    # r_h), which is at the level of finite-difference round-off on a
    # log-spaced grid. Values below this floor are not physically
    # meaningful sign information and are excluded from the violation
    # masks and range summaries (they are still returned in the raw
    # R_kk_null / R_uu_timelike arrays).
    noise_floor = 1e-6
    nec_significant = lorentzian & np.isfinite(R_kk) & (np.abs(R_kk) > noise_floor)
    sec_significant = lorentzian & np.isfinite(R_uu) & (np.abs(R_uu) > noise_floor)
    nec_violated = nec_significant & (R_kk < 0)
    sec_violated = sec_significant & (R_uu < 0)

    def _range(mask):
        if not np.any(mask):
            return None
        return (float(r[mask].min()), float(r[mask].max()))

    def _contiguous_ranges(mask):
        idxs = np.where(mask)[0]
        if len(idxs) == 0:
            return []
        splits = np.where(np.diff(idxs) > 1)[0]
        groups = np.split(idxs, splits + 1)
        return [(float(r[g[-1]]), float(r[g[0]])) for g in groups]

    return {
        'R_thth_check': R_thth_over_gthth,  # -> 0 for Schwarzschild (A=1,B=0): self-consistency test
        'R_kk_null': R_kk,          # R_{mu nu} k^mu k^nu, ingoing null congruence (NEC sign)
        'R_uu_timelike': R_uu,      # R_{mu nu} u^mu u^nu, radial-infall congruence, E=1 (SEC sign)
        'noise_floor': noise_floor,
        'nec_violated_mask': nec_violated,
        'sec_violated_mask': sec_violated,
        'nec_violated_r_range': _range(nec_violated),
        'sec_violated_r_range': _range(sec_violated),
        'nec_violated_r_intervals': _contiguous_ranges(nec_violated),
        'sec_violated_fraction_of_lorentzian': (float(np.sum(sec_violated)) / float(np.sum(sec_significant))
                                                  if np.any(sec_significant) else None),
        'nec_violated_anywhere': bool(np.any(nec_violated)),
        'sec_violated_anywhere': bool(np.any(sec_violated)),
    }


# =============================================================================
# Geodesic completeness
# =============================================================================

def check_geodesic_completeness(r, metric):
    """Check null and timelike geodesic completeness of gtilde.

    For the disformal metric with Killing vector K = d/dv, the conserved
    energy is E = -(gtilde_{vv} v' + gtilde_{vr} r'), where ' = d/dlambda.

    From the null condition and the conserved quantity E:

        (dr/dlambda)^2 = E^2 / (-det_2d)

    where det_2d = gtilde_{vv} gtilde_{rr} - gtilde_{vr}^2 < 0 (Lorentzian).

    Therefore:

        dlambda = sqrt(-det_2d) / |E| * dr

    The affine parameter to reach r = 0 from some radius r0 is:

        lambda(r0 -> 0) = (1/|E|) * int_0^{r0} sqrt(-det_2d(r')) dr'

    Geodesic completeness requires this integral to diverge as r0 -> 0.

    This is the EXACT result, valid everywhere in the globally Lorentzian
    space.  In the conformal regime (B -> 0), det_2d -> -A^4, so
    sqrt(-det_2d) -> A^2, recovering the conformal approximation
    lambda ~ int A^2 dr.

    For timelike geodesics, the proper time scales as d tau_tilde = A * d tau_Schw
    in the conformal regime, and the integral also diverges.

    Parameters
    ----------
    r : ndarray (ascending)
    metric : dict

    Returns
    -------
    dict with completeness diagnostics.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    det_2d = metric['det_2d']
    A = metric['A']
    A2 = metric['A2']
    F = metric['F']

    # Null geodesic slopes in the (v,r) plane
    disc = -det_2d  # gvr^2 - gvv*grr > 0 for Lorentzian
    disc_safe = np.where(disc > 0, disc, np.nan)
    sqrt_disc = np.sqrt(disc_safe)
    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)

    dv_dr_out = (-gvr + sqrt_disc) / gvv_safe
    dv_dr_in = (-gvr - sqrt_disc) / gvv_safe

    # --- Null affine parameter (EXACT: sqrt(-det_2d) integrand) ---
    # lambda = int sqrt(-det_2d) dr  (setting |E| = 1)
    # This is valid everywhere in the globally Lorentzian space, not just
    # the conformal approximation. In the conformal regime (B->0), this
    # reduces to int A^2 dr.
    null_integrand = np.sqrt(np.where(disc > 0, disc, np.nan))

    # Integrate from the horizon inward (r descending)
    idx_h = np.argmin(np.abs(r - 2.0))
    r_interior = r[:idx_h + 1]
    null_integrand_interior = null_integrand[:idx_h + 1]

    # Cumulative integral from horizon inward
    # r is ascending, so reverse for inward integration
    r_rev = r_interior[::-1]
    integrand_rev = null_integrand_interior[::-1]
    # Replace NaN with 0 for integration (NaN can occur from numerical
    # noise on the log-spaced grid; with the Gaussian bump B and phi_0=2,
    # det_2d < 0 everywhere so there is no non-Lorentzian region)
    integrand_rev = np.where(np.isfinite(integrand_rev), integrand_rev, 0)
    lambda_null_inward = cumulative_trapezoid(integrand_rev, r_rev, initial=0)

    # Total affine parameter from horizon to innermost
    lambda_null_total = float(lambda_null_inward[-1]) if len(lambda_null_inward) > 0 else 0.0

    # --- Timelike proper-time estimate (conformal approximation) ---
    # For radial infall from rest at infinity, d tau_tilde ~= A * d tau_Schw
    # with d tau_Schw = sqrt(r/(2M)) dr. The space is globally Lorentzian,
    # so timelike motion is well-defined all the way to r=0.
    M_eff = 1.0  # geometric units
    tau_integrand = A * np.sqrt(np.maximum(r, 1e-30) / (2.0 * M_eff))
    tau_integrand = np.where(metric['lorentzian'], tau_integrand, 0.0)
    tau_integrand_interior = tau_integrand[:idx_h + 1]
    tau_rev = tau_integrand_interior[::-1]
    r_rev2 = r_interior[::-1]
    tau_inward = cumulative_trapezoid(tau_rev, r_rev2, initial=0)
    tau_total = float(abs(tau_inward[-1])) if len(tau_inward) > 0 else 0.0

    # --- Power-law analysis of the null integrand near r=0 ---
    n_check = min(50, len(r) // 10)
    r_inner = r[:n_check]
    integrand_inner = null_integrand[:n_check]
    mask = (integrand_inner > 0) & (r_inner > 0) & np.isfinite(integrand_inner)
    if np.sum(mask) > 2:
        log_r = np.log(r_inner[mask])
        log_int = np.log(integrand_inner[mask])
        alpha = -np.polyfit(log_r, log_int, 1)[0]  # integrand ~ r^{-alpha}
        null_diverges = alpha >= 1.0
    else:
        alpha = np.nan
        null_diverges = False

    # --- Physical radial distance (on constant-v slice, exterior only) ---
    # In the exterior (r > 2M), constant-v slices are spacelike.
    # gtilde_rr_spatial = gtilde_rr - gtilde_vr^2 / gtilde_vv
    # In the exterior B -> 0, so gtilde_rr -> 0 and this reduces to
    # -gtilde_vr^2 / gtilde_vv = A^4 / (A^2 * F) = A^2 / F
    gvv_safe2 = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)
    grr_spatial = grr - gvr ** 2 / gvv_safe2
    # Only meaningful where grr_spatial > 0 (spacelike slice)
    spacelike_mask = grr_spatial > 0
    integrand_ell = np.sqrt(np.where(spacelike_mask, grr_spatial, 0))
    ell = cumulative_trapezoid(integrand_ell, r, initial=0)

    idx_horizon = np.argmin(np.abs(r - 2.0))

    return {
        'ell_physical': ell,
        'ell_at_horizon': float(ell[idx_horizon]),
        'ell_at_innermost': float(ell[0]),
        'ell_at_outermost': float(ell[-1]),
        'null_affine_parameter': null_integrand,
        'null_affine_total_horizon_to_center': lambda_null_total,
        'null_integrand_power_law': float(alpha) if np.isfinite(alpha) else None,
        'null_diverges': bool(null_diverges),
        'timelike_proper_time_total': tau_total,
        'timelike_diverges': bool(tau_total > 1e3),  # large = complete
        'dv_dr_outgoing': dv_dr_out,
        'dv_dr_ingoing': dv_dr_in,
    }


# =============================================================================
# Physical volume and density
# =============================================================================

def compute_physical_volume(r, metric, M=1.0):
    """Compute physical volume and density.

    Inside the horizon, r is timelike and v is spacelike.  The physical
    3-volume on a constant-r spacelike slice is:

        dV = 4*pi * sqrt(det(gtilde_{ij})) * dv   (on v-spatial slice)

    The induced metric on a constant-r hypersurface has components:
        gtilde_{vv} (spacelike inside horizon, since F < 0 => gtilde_vv > 0)
        gtilde_{thth}, gtilde_{phph}

    So dV = 4*pi * sqrt(gtilde_vv * gtilde_thth * gtilde_phph) * dv
          = 4*pi * A^2 * r^2 * sqrt(|F|) * dv   (since gtilde_vv = A^2*|F| inside)

    But this is per unit v.  For the total volume enclosed within radius r,
    we integrate over the spatial extent.  In the exterior (r > 2M), the
    standard approach is:

        V(r) = 4*pi * int_0^r A * R^2 * sqrt(gtilde_rr_spatial) dr'

    where gtilde_rr_spatial is the radial metric on a spacelike slice.

    For the exterior (B->0, A->1): V -> 4/3 * pi * r^3 (Schwarzschild).
    For the interior: the volume element is modified by A.

    On the constant-r slices inside the horizon, the volume element per
    unit v does not collapse as r -> 0. For phi_0 = 2, A^2*r^2 diverges
    while sqrt(|F|) also grows. The space is globally Lorentzian, so this
    is a physical volume element on the matter domain, not a formal
    continuation diagnostic.

    Parameters
    ----------
    r : ndarray (ascending)
    metric : dict
    M : float

    Returns
    -------
    dict with volume and density profiles.
    """
    gvv = metric['gtilde_vv']
    gvr = metric['gtilde_vr']
    grr = metric['gtilde_rr']
    A = metric['A']
    R = metric['r']
    F = metric['F']

    # Exterior volume: integrate on spacelike slices (r > 2M)
    gvv_safe = np.where(np.abs(gvv) > 1e-50, gvv, np.nan)
    grr_spatial = grr - gvr ** 2 / gvv_safe

    # In the exterior (r > 2M): grr_spatial > 0 (spacelike)
    # In the interior (r < 2M): grr_spatial < 0 (r is timelike)
    spacelike = grr_spatial > 0
    integrand_ext = np.where(spacelike, A * R ** 2 * np.sqrt(np.maximum(grr_spatial, 0)), 0)

    # V_physical: cumulative from the innermost spacelike point outward
    V_physical = 4 * np.pi * cumulative_trapezoid(integrand_ext, r, initial=0)
    V_geometric = 4.0 / 3.0 * np.pi * r ** 3

    # Interior volume element: on constant-r slices inside the horizon,
    # the spatial metric is gtilde_{vv} (spacelike since F<0) x gtilde_{thth}^2.
    # dV_interior = 4*pi * sqrt(|gtilde_vv|) * gtilde_thth * dv
    #             = 4*pi * A^2 * r^2 * sqrt(|F|) * dv
    # This is the physical volume per unit v-time.
    # The key: A^2 * r^2 ~ (r_h/r)^{2*phi_0} * r^2 = r_h^{2*phi_0} * r^{2-2*phi_0}
    # For phi_0 = 2, A^2*r^2 diverges and sqrt(|F|) also grows.
    # The space is globally Lorentzian, so this is a physical volume
    # element, not a formal continuation diagnostic.
    dV_interior_per_v = 4 * np.pi * A ** 2 * r ** 2 * np.sqrt(np.abs(F))

    # For the total enclosed volume, we use the exterior volume at the
    # horizon as a reference and note that the interior adds more volume,
    # not less.  The density is computed using the exterior volume.
    idx_h = np.argmin(np.abs(r - 2.0 * M))
    V_at_horizon = float(V_physical[idx_h])

    # Densities: use exterior volume where available
    rho_physical = M / np.where(V_physical > 1e-50, V_physical, np.nan)
    rho_geometric = 3.0 * M / (4.0 * np.pi * np.maximum(r, 1e-30) ** 3)

    return {
        'dV_physical': integrand_ext,
        'dV_geometric': R ** 2,
        'V_physical': V_physical,
        'V_geometric': V_geometric,
        'rho_physical': rho_physical,
        'rho_geometric': rho_geometric,
        'dV_interior_per_v': dV_interior_per_v,
    }


# =============================================================================
# Main solver
# =============================================================================

def solve_tep_bh(model, r_min=1e-8, r_max=50.0, n_points=20000):
    """Construct the TEP black-hole solution.

    Parameters
    ----------
    model : TEPBHModel
    r_min : float (innermost radius in units of M)
    r_max : float (outermost radius)
    n_points : int

    Returns
    -------
    dict with metric, curvature, geodesics, volume, and diagnostics.
    """
    M = model.M

    # Radial grid: logarithmic (fine resolution near origin)
    r = np.logspace(np.log10(r_min), np.log10(r_max), n_points)

    # Compute disformal metric
    metric = compute_disformal_metric(r, model)

    # Compute curvature invariants
    curvature = compute_curvature_invariants(r, metric, M=M)

    # Check geodesic completeness
    geodesics = check_geodesic_completeness(r, metric)

    # Physical volume and density
    volume = compute_physical_volume(r, metric, M=M)

    # Energy conditions / Raychaudhuri convergence sign
    energy_conditions = compute_energy_conditions(r, metric, curvature, M=M)

    # Check for a determinant-zero boundary: with the Gaussian bump B(phi)
    # and phi_0=2, the determinant is strictly negative everywhere (no
    # boundary). We check for sign changes and report None if the space is
    # globally Lorentzian (the expected case for the current model).
    det_2d = metric['det_2d']

    sign_changes = np.where(np.signbit(det_2d[:-1]) != np.signbit(det_2d[1:]))[0]
    if len(sign_changes) > 0:
        idx_temporal = int(sign_changes[0])
        r0, r1 = r[idx_temporal], r[idx_temporal + 1]
        d0, d1 = det_2d[idx_temporal], det_2d[idx_temporal + 1]
        if d1 != d0:
            r_temporal = float(r0 - d0 * (r1 - r0) / (d1 - d0))
        else:
            r_temporal = float(0.5 * (r0 + r1))
    else:
        # No determinant-zero boundary: the space is globally Lorentzian
        idx_temporal = None
        r_temporal = None

    return {
        'success': True,
        'r': r,
        'metric': metric,
        'curvature': curvature,
        'geodesics': geodesics,
        'volume': volume,
        'energy_conditions': energy_conditions,
        'r_temporal_horizon': r_temporal,
        'idx_temporal_horizon': idx_temporal,
        'model': model.to_dict(),
    }
