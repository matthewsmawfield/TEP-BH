#!/usr/bin/env python3
"""Step 06 (Inference): Phantom Mass posterior via MCMC.

Sample the posterior distribution of the TEP parameters using
Markov Chain Monte Carlo (MCMC), with emcee if available, otherwise
a simple Metropolis-Hastings sampler.

The key output is the posterior distribution of M_phantom^T,
which quantifies the uncertainty on the Phantom Mass residual.

POSTERIOR:
  p(theta | data) ∝ L(data | theta) * p(theta)

  where theta = [M_BH, R0, a, e, i, omega, Omega, T, t_peri, x0, y0, eta_TEP]

  Priors:
    M_BH:   Uniform(1e5, 1e7) M_sun
    R0:     Uniform(7000, 10000) pc
    a:      Uniform(0.05, 0.5) arcsec
    e:      Uniform(0, 0.99)
    i:      Uniform(0, 180) deg
    omega:  Uniform(0, 360) deg
    Omega:  Uniform(0, 360) deg
    T:      Uniform(5, 30) yr
    t_peri: Uniform(1990, 2030) yr
    x0:     Uniform(-1, 1) arcsec
    y0:     Uniform(-1, 1) arcsec
    eta_TEP: Uniform(0, 1)  (mass-inflation branch, alpha_GB < 0)

  Derived quantity:
    M_phantom^T = M_GR - M_matter(eta_TEP, M_BH, ...)

Outputs:
  data/processed/sstar_phantom_mass_posterior.json
  results/step_49_posterior.json
  logs/step_49_posterior.log
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps" / "inference"))

from step_46_tep_transfer import tep_residuals
from step_45_mass_bias import (
    temporal_transfer_factor,
    spatial_calibration_factor,
    dynamical_modification_factor,
    G_SI, c_si, M_sun_kg, AU_m,
)

from bh_common import (
    ensure_dirs,
    make_step_logger,
    print_status,
    write_json,
    finalize_result,
    PROCESSED_DIR,
    RESULTS_DIR,
)

STEP_ID = "step_49_posterior"


# ============================================================
# Prior and posterior
# ============================================================
PRIOR_BOUNDS = {
    "M_BH": (1e5, 1e7),
    "R0": (7000, 10000),
    "a": (0.05, 0.5),
    "e": (0.0, 0.99),
    "i": (0, 180),
    "omega": (0, 360),
    "Omega": (0, 360),
    "T": (5, 30),
    "t_peri": (1990, 2030),
    "x0": (-1, 1),
    "y0": (-1, 1),
    "eta_TEP": (0.0, 1.0),
}

PARAM_NAMES = ["M_BH", "R0", "a", "e", "i", "omega", "Omega", "T", "t_peri", "x0", "y0", "eta_TEP"]

# R0 prior (same soft prior used in steps 01 and 03 to break the M_BH--R0 degeneracy)
R0_PRIOR = 8320.0
SIGMA_R0 = 60.0


def log_prior(theta):
    """Log prior (uniform within bounds)."""
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PRIOR_BOUNDS[name]
        if theta[i] < lo or theta[i] > hi:
            return -np.inf
    return 0.0


def log_likelihood(theta, astrometry, spectroscopy):
    """Log likelihood = -chi^2 / 2."""
    res = tep_residuals(theta, astrometry, spectroscopy)
    chi2 = np.sum(res**2)
    return -chi2 / 2


def log_posterior(theta, astrometry, spectroscopy):
    """Log posterior = log prior + log likelihood + soft R0 prior."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    R0 = theta[1]
    r0_prior_term = -0.5 * ((R0 - R0_PRIOR) / SIGMA_R0) ** 2
    return lp + log_likelihood(theta, astrometry, spectroscopy) + r0_prior_term


# ============================================================
# MCMC sampler
# ============================================================
def run_metropolis_hastings(log_post, theta0, n_samples, n_burn, astrometry, spectroscopy):
    """Simple Metropolis-Hastings MCMC sampler."""
    ndim = len(theta0)
    samples = np.zeros((n_samples, ndim))
    current = theta0.copy()
    current_lp = log_post(current, astrometry, spectroscopy)

    # Proposal scale (relative to parameter range)
    scales = np.array([
        0.01 * (PRIOR_BOUNDS[name][1] - PRIOR_BOUNDS[name][0])
        for name in PARAM_NAMES
    ])

    n_accept = 0
    for i in range(n_samples + n_burn):
        proposal = current + np.random.normal(0, scales)
        proposal_lp = log_post(proposal, astrometry, spectroscopy)

        if proposal_lp > current_lp:
            accept = True
        else:
            accept = np.random.random() < np.exp(proposal_lp - current_lp)

        if accept:
            current = proposal
            current_lp = proposal_lp
            n_accept += 1

        if i >= n_burn:
            samples[i - n_burn] = current

    accept_rate = n_accept / (n_samples + n_burn)
    return samples, accept_rate


def run_emcee(log_post, theta0, n_walkers, n_samples, n_burn, astrometry, spectroscopy):
    """Run emcee ensemble sampler if available."""
    import emcee

    ndim = len(theta0)

    # Initialize walkers around theta0
    perturbations = np.array([
        0.01 * (PRIOR_BOUNDS[name][1] - PRIOR_BOUNDS[name][0])
        for name in PARAM_NAMES
    ])
    p0 = theta0 + np.random.normal(0, perturbations, (n_walkers, ndim))

    # Ensure within bounds
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PRIOR_BOUNDS[name]
        p0[:, i] = np.clip(p0[:, i], lo, hi)

    sampler = emcee.EnsembleSampler(
        n_walkers, ndim, log_post,
        args=(astrometry, spectroscopy),
    )

    # Burn-in
    state = sampler.run_mcmc(p0, n_burn, progress=False)
    sampler.reset()

    # Production
    sampler.run_mcmc(state, n_samples, progress=False)

    samples = sampler.get_chain(discard=0, flat=True)
    accept_rate = np.mean(sampler.acceptance_fraction)

    return samples, accept_rate


# ============================================================
# Derived quantities
# ============================================================
def compute_phantom_mass(theta, M_GR):
    """Compute M_phantom from a parameter sample.

    The mass-bias ratio is evaluated relative to GR (eta=0) so that
    eta=0 gives ratio=1 and M_phantom=0 identically.
    """
    M_BH = theta[0]  # M_g for this sample
    eta = theta[11]
    a = theta[2]
    e = theta[3]
    R0 = theta[1]

    # Compute transfer factors at pericentre
    R_s_m = 2 * G_SI * M_GR * M_sun_kg / c_si**2
    R_s_AU = R_s_m / AU_m
    r_peri_AU = a * (1 - e) * R0
    r_peri_rs = r_peri_AU / R_s_AU

    T_P = temporal_transfer_factor(r_peri_rs, np.inf, eta=eta)
    S_a = spatial_calibration_factor(r_peri_rs, M_GR, R0, eta=eta)
    D_dyn = dynamical_modification_factor(r_peri_rs, eta=eta)

    T_P_GR = temporal_transfer_factor(r_peri_rs, np.inf, eta=0.0)
    S_a_GR = spatial_calibration_factor(r_peri_rs, M_GR, R0, eta=0.0)
    D_dyn_GR = dynamical_modification_factor(r_peri_rs, eta=0.0)

    ratio = ((S_a / S_a_GR)**3 * (D_dyn / D_dyn_GR)) / (T_P / T_P_GR)**2
    M_matter = M_BH / ratio
    M_phantom = M_BH - M_matter

    return M_phantom, M_matter, ratio


# ============================================================
# Main
# ============================================================
def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)
    print_status("=" * 70, "TITLE")
    print_status("STEP 06 (INFERENCE): PHANTOM MASS POSTERIOR (MCMC)", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("")
    print_status("Sampling the posterior distribution of TEP parameters.", "INFO")
    print_status("Key output: posterior of M_phantom^T", "INFO")
    print_status("")

    # --- Load data and fits ---
    astrometry_path = PROCESSED_DIR / "sstar_astrometry.json"
    spectroscopy_path = PROCESSED_DIR / "sstar_spectroscopy.json"
    tep_fit_path = PROCESSED_DIR / "sstar_tep_fit.json"
    gr_fit_path = PROCESSED_DIR / "sstar_gr_fit.json"

    if not all(p.exists() for p in [astrometry_path, spectroscopy_path, tep_fit_path, gr_fit_path]):
        print_status("ERROR: Previous fits not found. Run steps 00-05 first.", "ERROR")
        return

    astrometry = json.loads(astrometry_path.read_text())["data"]
    spectroscopy = json.loads(spectroscopy_path.read_text())["data"]
    tep_fit = json.loads(tep_fit_path.read_text())
    gr_fit = json.loads(gr_fit_path.read_text())

    M_GR = gr_fit["M_BH_GR_Msun"]

    # Set R0 prior from GR fit to break the M_BH--R0 degeneracy
    global R0_PRIOR, SIGMA_R0
    R0_PRIOR = gr_fit["best_fit"]["R0 (pc)"]["value"]
    SIGMA_R0 = gr_fit["best_fit"]["R0 (pc)"]["error"]
    if SIGMA_R0 <= 0:
        SIGMA_R0 = 60.0

    # Initial point from TEP fit
    theta0 = np.array([
        tep_fit["best_fit"][f"{name} (M_sun)" if name == "M_BH" else
                            f"{name} (pc)" if name == "R0" else
                            f"{name} (arcsec)" if name in ["a", "x0", "y0"] else
                            f"{name} (deg)" if name in ["i", "omega", "Omega"] else
                            f"{name} (yr)" if name in ["T", "t_peri"] else
                            name]["value"]
        for name in PARAM_NAMES
    ])

    is_synthetic = json.loads(astrometry_path.read_text()).get("is_synthetic", False)
    print_status(f"  Data: {'SYNTHETIC' if is_synthetic else 'REAL'}", "INFO")
    print_status(f"  M_GR-fit = {M_GR:.4e} M_sun", "INFO")
    print_status("")

    # --- Run MCMC ---
    # Fix random seed for reproducible posterior numbers
    np.random.seed(42)

    n_samples = 10000
    n_burn = 1000

    # Try emcee first, fall back to Metropolis-Hastings
    use_emcee = False
    try:
        import emcee
        use_emcee = True
    except ImportError:
        pass

    if use_emcee:
        n_walkers = 24
        print_status(f"--- Running emcee MCMC ---", "INFO")
        print_status(f"  Walkers: {n_walkers}", "INFO")
        print_status(f"  Samples: {n_samples}", "INFO")
        print_status(f"  Burn-in: {n_burn}", "INFO")
        print_status("")

        samples, accept_rate = run_emcee(
            log_posterior, theta0, n_walkers, n_samples, n_burn,
            astrometry, spectroscopy,
        )
    else:
        print_status(f"--- Running Metropolis-Hastings MCMC ---", "INFO")
        print_status(f"  (emcee not available; using simple MH sampler)", "WARN")
        print_status(f"  Samples: {n_samples}", "INFO")
        print_status(f"  Burn-in: {n_burn}", "INFO")
        print_status("")

        samples, accept_rate = run_metropolis_hastings(
            log_posterior, theta0, n_samples, n_burn,
            astrometry, spectroscopy,
        )

    print_status(f"  Acceptance rate: {accept_rate:.3f}", "INFO")
    print_status(f"  Total samples: {len(samples)}", "INFO")
    print_status("")

    # --- Compute posterior statistics ---
    print_status("--- Posterior statistics ---", "TITLE")

    posterior_stats = {}
    for i, name in enumerate(PARAM_NAMES):
        vals = samples[:, i]
        median = np.median(vals)
        p16, p84 = np.percentile(vals, [16, 84])
        mean = np.mean(vals)
        std = np.std(vals)
        posterior_stats[name] = {
            "median": float(median),
            "mean": float(mean),
            "std": float(std),
            "p16": float(p16),
            "p84": float(p84),
        }
        unit = {"M_BH": "M_sun", "R0": "pc", "a": "arcsec", "e": "", "i": "deg",
                "omega": "deg", "Omega": "deg", "T": "yr", "t_peri": "yr",
                "x0": "arcsec", "y0": "arcsec", "eta_TEP": ""}[name]
        print_status(f"  {name:10s} = {median:12.4f} +{p84-median:8.4f} / -{median-p16:8.4f}  {unit}", "INFO")

    print_status("")

    # --- Compute Phantom Mass posterior ---
    print_status("--- Phantom Mass posterior ---", "TITLE")

    M_phantom_samples = np.zeros(len(samples))
    M_matter_samples = np.zeros(len(samples))
    ratio_samples = np.zeros(len(samples))

    for i in range(len(samples)):
        M_ph, M_mat, ratio = compute_phantom_mass(samples[i], M_GR)
        M_phantom_samples[i] = M_ph
        M_matter_samples[i] = M_mat
        ratio_samples[i] = ratio

    M_ph_median = np.median(M_phantom_samples)
    M_ph_p16, M_ph_p84 = np.percentile(M_phantom_samples, [16, 84])
    M_ph_mean = np.mean(M_phantom_samples)
    M_ph_std = np.std(M_phantom_samples)

    M_mat_median = np.median(M_matter_samples)
    M_mat_p16, M_mat_p84 = np.percentile(M_matter_samples, [16, 84])

    ratio_median = np.median(ratio_samples)
    ratio_p16, ratio_p84 = np.percentile(ratio_samples, [16, 84])

    # Sign determination
    n_positive = np.sum(M_phantom_samples > 0)
    n_negative = np.sum(M_phantom_samples < 0)
    n_zero = np.sum(np.abs(M_phantom_samples) < 1.0)  # within 1 M_sun of zero
    frac_positive = n_positive / len(samples)
    frac_negative = n_negative / len(samples)

    print_status(f"  M_phantom^T:", "INFO")
    print_status(f"    median = {M_ph_median:.4e} M_sun", "INFO")
    print_status(f"    mean   = {M_ph_mean:.4e} M_sun", "INFO")
    print_status(f"    std    = {M_ph_std:.4e} M_sun", "INFO")
    print_status(f"    68% CI = [{M_ph_p16:.4e}, {M_ph_p84:.4e}] M_sun", "INFO")
    print_status(f"    fraction positive: {frac_positive:.3f}", "INFO")
    print_status(f"    fraction negative: {frac_negative:.3f}", "INFO")
    print_status("")

    print_status(f"  M_matter^TEP:", "INFO")
    print_status(f"    median = {M_mat_median:.4e} M_sun", "INFO")
    print_status(f"    68% CI = [{M_mat_p16:.4e}, {M_mat_p84:.4e}] M_sun", "INFO")
    print_status("")

    print_status(f"  Ratio M_g/M_matter:", "INFO")
    print_status(f"    median = {ratio_median:.8f}", "INFO")
    print_status(f"    68% CI = [{ratio_p16:.8f}, {ratio_p84:.8f}]", "INFO")
    print_status("")

    # --- Sign determination ---
    print_status("--- SIGN DETERMINATION (Gate -1) ---", "TITLE")
    if frac_positive > 0.95:
        sign = "positive"
        print_status(f"  M_phantom^T > 0 (fraction positive: {frac_positive:.3f})", "INFO")
        print_status("  -> Spatial magnification dominates", "INFO")
    elif frac_negative > 0.95:
        sign = "negative"
        print_status(f"  M_phantom^T < 0 (fraction negative: {frac_negative:.3f})", "INFO")
        print_status("  -> Temporal stretching dominates", "INFO")
    else:
        sign = "indeterminate"
        print_status(f"  M_phantom^T sign is INDETERMINATE", "WARN")
        print_status(f"  (positive: {frac_positive:.3f}, negative: {frac_negative:.3f})", "WARN")
    print_status("")

    # --- Save ---
    # Save a subsample of the chain
    n_save = min(500, len(samples))
    indices = np.random.choice(len(samples), n_save, replace=False)
    chain_subset = samples[indices].tolist()

    result = {
        "step": STEP_ID,
        "description": "Phantom Mass posterior via MCMC",
        "is_synthetic_data": is_synthetic,
        "sampler": "emcee" if use_emcee else "metropolis-hastings",
        "n_samples": len(samples),
        "acceptance_rate": float(accept_rate),
        "posterior_stats": posterior_stats,
        "phantom_mass_posterior": {
            "median": float(M_ph_median),
            "mean": float(M_ph_mean),
            "std": float(M_ph_std),
            "p16": float(M_ph_p16),
            "p84": float(M_ph_p84),
            "fraction_positive": float(frac_positive),
            "fraction_negative": float(frac_negative),
            "sign": sign,
        },
        "matter_mass_posterior": {
            "median": float(M_mat_median),
            "p16": float(M_mat_p16),
            "p84": float(M_mat_p84),
        },
        "ratio_posterior": {
            "median": float(ratio_median),
            "p16": float(ratio_p16),
            "p84": float(ratio_p84),
        },
        "M_GR_Msun": float(M_GR),
        "chain_subset": chain_subset,
        "chain_subset_n": n_save,
        "timestamp": datetime.now().isoformat(),
        "status": "success",
    }

    posterior_path = PROCESSED_DIR / "sstar_phantom_mass_posterior.json"
    write_json(posterior_path, result)
    print_status(f"  Saved: {posterior_path}", "INFO")

    result_path = RESULTS_DIR / f"{STEP_ID}.json"
    write_json(result_path, result)
    finalize_result(STEP_ID, result, "Phantom Mass posterior via MCMC", "Posterior distribution of M_phantom^T sampled")

    print_status("")
    print_status("  INFERENCE PIPELINE COMPLETE.", "TITLE")
    print_status("  All steps 00-06 executed.", "INFO")


if __name__ == "__main__":
    main()
