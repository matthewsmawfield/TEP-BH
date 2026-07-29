#!/usr/bin/env python3
"""TEP-BH Figure Generation.

Generates the 8 figures specified in the TEP-BH paper plan.

Outputs (in results/figures/):
  fig1_two_metric.png/pdf, fig2_radial_distance.png/pdf, ...

Usage:
    python scripts/figures/generate_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bh_common import (
    ensure_dirs, make_step_logger, print_status, rel, FIGURES_DIR,
    TEPBHModel, solve_tep_bh,
)
from scipy.integrate import cumulative_trapezoid


# --- Inline accretion functions (Schwarzschild exterior, A cancels) ---

def _compute_circular_geodesics(r, M=1.0):
    r_safe = np.where(r > 3.0 * M, r, np.nan)
    f = 1.0 - 2.0 * M / r_safe
    denom = np.sqrt(np.maximum(1.0 - 3.0 * M / r_safe, 1e-30))
    E = f / denom
    L = np.sqrt(np.maximum(M * r_safe / np.maximum(1.0 - 3.0 * M / r_safe, 1e-30), 0))
    Omega = np.sqrt(M / r_safe**3)
    return {'r': r, 'E': np.where(r > 3.0 * M, E, np.nan),
            'L': np.where(r > 3.0 * M, L, np.nan),
            'Omega': np.where(r > 3.0 * M, Omega, np.nan)}


def _compute_epicyclic_frequencies(r, M=1.0):
    r_safe = np.where(r > 3.0 * M, r, np.nan)
    Omega_sq = M / r_safe**3
    Omega_r_sq = Omega_sq * (1.0 - 6.0 * M / r_safe)
    return {
        'r': r,
        'Omega_r': np.sqrt(np.maximum(np.where(r > 3.0 * M, Omega_r_sq, np.nan), 0)),
        'Omega_theta': np.sqrt(np.maximum(np.where(r > 3.0 * M, Omega_sq, np.nan), 0)),
        'ratio_r_to_theta': np.sqrt(np.maximum(np.where(r > 3.0 * M, Omega_r_sq, np.nan), 0)) /
                            np.sqrt(np.maximum(np.where(r > 3.0 * M, Omega_sq, np.nan), 0)),
    }


def _compute_isco_properties(M=1.0):
    r_isco = 6.0 * M
    return {
        'r_isco': float(r_isco),
        'E_isco': float(np.sqrt(8.0 / 9.0)),
        'L_isco': float(2.0 * np.sqrt(3.0) * M),
        'Omega_isco': float(1.0 / (6.0 * np.sqrt(6.0) * M)),
        'radiative_efficiency': float(1.0 - np.sqrt(8.0 / 9.0)),
        'radiative_efficiency_percent': float((1.0 - np.sqrt(8.0 / 9.0)) * 100),
        'redshift_factor_isco': float(np.sqrt(1.0 - 2.0 * M / r_isco)),
        'redshift_z_isco': float(1.0 / np.sqrt(1.0 - 2.0 * M / r_isco) - 1.0),
    }


def _compute_redshift_profile(r, M=1.0):
    r_safe = np.where(r > 2.0 * M + 1e-10, r, np.nan)
    g = np.sqrt(1.0 - 2.0 * M / r_safe)
    return {'r': r, 'g_redshift': g, 'redshift_z': 1.0 / g - 1.0}


# --- Inline perturbation functions ---

def _compute_tortoise_coordinate(r, metric):
    F = metric['F']
    F_abs = np.abs(F)
    F_safe = np.where(F_abs > 1e-30, F_abs, 1e-30)
    drstar_dr = 1.0 / F_safe
    idx_ref = np.argmin(np.abs(r - 3.0))
    r_star = np.zeros_like(r)
    r_star[idx_ref:] = cumulative_trapezoid(drstar_dr[idx_ref:], r[idx_ref:], initial=0)
    r_star[:idx_ref+1] = -cumulative_trapezoid(
        drstar_dr[idx_ref::-1], r[idx_ref::-1], initial=0)[::-1]
    return r_star


def _compute_regge_wheeler_potential(r, r_star, metric, l=2):
    F = metric['F']
    F_safe = np.where(np.abs(F) > 1e-30, F, np.nan)
    V_RW_schw = F_safe * (l * (l + 1) / r**2 - 6.0 / r**3)
    # GWs propagate on geometric g (Schwarzschild); no A^{-4} matter-metric scaling.
    V_RW = V_RW_schw
    V_RW = np.where(np.isfinite(V_RW), V_RW, 0)
    V_RW_schw = np.where(np.isfinite(V_RW_schw), V_RW_schw, 0)
    return V_RW, V_RW_schw

STEP_ID = "generate_figures"

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.linewidth': 0.8,
    'figure.figsize': (8, 6),
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

FIG_DIR = FIGURES_DIR


def get_solution():
    """Get the standard TEP-BH solution for plotting."""
    model = TEPBHModel(
        beta_A=-1.0, B0=1.0, n_B=2.0,
        phi_0=2.0, delta=0.05, M=1.0,
        sigma_B=1.5,
    )
    return model, solve_tep_bh(model)


def save_fig(fig, name):
    """Save figure as PNG and PDF."""
    fig.savefig(FIG_DIR / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)
    print_status(f"  Saved {name}.png and {name}.pdf", "SUCCESS")


def fig1_two_metric(model, sol):
    """Figure 1 — Two-metric causal architecture."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    r = sol['r']
    m = sol['metric']

    # Geometric metric (Schwarzschild)
    ax1.plot(r, -m['F'], 'b-', label=r'$-g_{vv} = F = 1 - 2M/r$')
    ax1.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$ (horizon)')
    ax1.set_xscale('log')
    ax1.set_yscale('symlog', linthresh=1e-6)
    ax1.set_xlabel(r'$r / M$')
    ax1.set_ylabel(r'metric component')
    ax1.set_title(r'Geometric metric $g_{\mu\nu}$ (Schwarzschild)')
    ax1.legend(fontsize=10)
    ax1.set_xlim(1e-4, 50)

    # Physical metric (disformal)
    ax2.plot(r, np.abs(m['gtilde_vv']), 'r-', label=r'$|\tilde{g}_{vv}|$')
    ax2.plot(r, m['gtilde_rr'], 'g-', label=r'$\tilde{g}_{rr} = B(\phi\')^2$')
    ax2.plot(r, m['gtilde_thth'], 'm-', label=r'$\tilde{g}_{\theta\theta} = A^2 r^2$')
    r_t = sol.get('r_temporal_horizon')
    if r_t:
        ax2.axvline(r_t, color='orange', ls=':', linewidth=2, label=r'$r_t$ (temporal horizon)')
    ax2.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$r / M$')
    ax2.set_ylabel(r'metric component')
    ax2.set_title(r'Physical metric $\tilde{g}_{\mu\nu}$ (disformal)')
    ax2.legend(fontsize=10)
    ax2.set_xlim(1e-4, 50)

    plt.tight_layout()
    save_fig(fig, 'fig1_two_metric')


def fig2_radial_distance(model, sol):
    """Figure 2 — Geometric versus physical radial distance."""
    fig, ax = plt.subplots(figsize=(8, 6))
    r = sol['r']
    ell = sol['geodesics']['ell_physical']

    ax.plot(r, r, 'b--', label=r'Geometric $r$ (Schwarzschild)')
    ax.plot(r, ell, 'r-', label=r'Physical $\tilde{\ell}(r)$ (disformal)')

    ax.axvline(2.0, color='gray', ls='--', alpha=0.5)
    r_t = sol.get('r_temporal_horizon')
    if r_t:
        ax.axvline(r_t, color='orange', ls=':', linewidth=2, label=r'$r_t$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'distance $/ M$')
    ax.set_title('Geometric vs physical radial distance')
    ax.legend(loc='upper left', fontsize=10)
    ax.set_xlim(1e-4, 50)

    plt.tight_layout()
    save_fig(fig, 'fig2_radial_distance')


def fig3_clock_projection(model, sol):
    """Figure 3 — Clock projection."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    r = sol['r']
    m = sol['metric']

    ax1.plot(r, m['phi'], 'b-', label=r'$\phi(r)$')
    ax1.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    ax1.set_xscale('log')
    ax1.set_xlabel(r'$r / M$')
    ax1.set_ylabel(r'$\phi$')
    ax1.set_title('Scalar field profile')
    ax1.legend(fontsize=10)
    ax1.set_xlim(1e-4, 50)

    ax2.plot(r, m['A'], 'r-', label=r'$A(\phi) = e^{\beta_A \phi}$')
    ax2.axhline(1.0, color='gray', ls='--', alpha=0.5)
    ax2.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    r_t = sol.get('r_temporal_horizon')
    if r_t:
        ax2.axvline(r_t, color='orange', ls=':', linewidth=2, label=r'$r_t$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$r / M$')
    ax2.set_ylabel(r'$A(\phi)$')
    ax2.set_title('Conformal clock factor')
    ax2.legend(fontsize=10)
    ax2.set_xlim(1e-4, 50)

    plt.tight_layout()
    save_fig(fig, 'fig3_clock_projection')


def fig4_curvature(model, sol):
    """Figure 4 — Curvature invariants."""
    fig, ax = plt.subplots(figsize=(8, 6))
    r = sol['r']
    K = sol['curvature']['Kretschmann']
    K_schw = sol['curvature']['Kretschmann_schwarzschild']

    ax.plot(r, K, 'r-', linewidth=2, label=r'$\tilde{K}$ (TEP physical metric)')
    ax.plot(r, K_schw, 'b--', linewidth=2, label=r'$K$ (Schwarzschild)')
    ax.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    r_t = sol.get('r_temporal_horizon')
    if r_t:
        ax.axvline(r_t, color='orange', ls=':', linewidth=2, label=r'$r_t$ (temporal horizon)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'Kretschmann scalar $K$')
    ax.set_title('Curvature invariants: TEP vs Schwarzschild')
    ax.legend(fontsize=10)
    ax.set_xlim(1e-4, 50)
    ax.set_ylim(1e-10, 1e55)

    plt.tight_layout()
    save_fig(fig, 'fig4_curvature')


def fig5_volume_density(model, sol):
    """Figure 5 — Physical volume and density."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    r = sol['r']
    vol = sol['volume']

    ax1.plot(r, vol['V_physical'], 'r-', label=r'$\tilde{V}$ (physical)')
    ax1.plot(r, vol['V_geometric'], 'b--', label=r'$V$ (geometric)')
    ax1.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$r / M$')
    ax1.set_ylabel(r'Volume')
    ax1.set_title('Physical vs geometric volume')
    ax1.legend(fontsize=10)
    ax1.set_xlim(1e-4, 50)

    rho_phys = vol['rho_physical']
    rho_geom = vol['rho_geometric']
    # Mask non-finite values for clean log-scale plotting
    rho_phys_plot = np.where(np.isfinite(rho_phys) & (rho_phys > 0), rho_phys, np.nan)
    rho_geom_plot = np.where(np.isfinite(rho_geom) & (rho_geom > 0), rho_geom, np.nan)
    ax2.plot(r, rho_phys_plot, 'r-', label=r'$\tilde{\rho}$ (physical)')
    ax2.plot(r, rho_geom_plot, 'b--', label=r'$\rho$ (geometric)')
    ax2.axvline(2.0, color='gray', ls='--', alpha=0.5, label=r'$r = 2M$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$r / M$')
    ax2.set_ylabel(r'Density')
    ax2.set_title('Physical vs geometric density')
    ax2.legend(fontsize=10)
    ax2.set_xlim(1e-4, 50)

    plt.tight_layout()
    save_fig(fig, 'fig5_volume_density')


def fig6_penrose(model, sol):
    """Figure 6 — Schematic Penrose diagrams."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    labels = ['Schwarzschild', 'Gravastar', 'Regular BH', 'TEP temporal horizon']
    descriptions = [
        r'Singularity\nat $r=0$\n(curvature $\to \infty$)',
        'Thin shell\nat $r_0$\n(material surface)',
        r'De Sitter\ncore\n($r \to 0$ regular)',
        r'Temporal\nhorizon at $r_t$\n($\tilde{\ell} \to \infty$)',
    ]

    for ax, label, desc in zip(axes, labels, descriptions):
        # Draw diamond-shaped Penrose diagram
        diamond = plt.Polygon([(0.5, 1), (1, 0.5), (0.5, 0), (0, 0.5)],
                              fill=False, edgecolor='black', linewidth=1.5)
        ax.add_patch(diamond)

        if 'Schwarzschild' in label:
            # Singularity: jagged line at top
            ax.plot([0.3, 0.7], [0.95, 0.95], 'r-', linewidth=2)
            ax.plot([0.35, 0.65], [0.98, 0.98], 'r-', linewidth=2)
        elif 'Gravastar' in label:
            # Thin shell
            ax.plot([0.35, 0.65], [0.65, 0.65], 'b-', linewidth=3)
        elif 'Regular' in label:
            # Smooth core
            ax.plot([0.3, 0.7], [0.1, 0.1], 'g-', linewidth=2)
        elif 'TEP' in label:
            # Temporal horizon
            ax.plot([0.35, 0.65], [0.3, 0.3], color='m', linewidth=2, linestyle=':')

        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.set_aspect('equal')
        ax.set_title(label, fontsize=11)
        ax.text(0.5, -0.05, desc, ha='center', va='top', fontsize=9,
                transform=ax.transData)
        ax.set_xticks([])
        ax.set_yticks([])

    plt.suptitle('Causal structure comparison', fontsize=14, y=1.02)
    plt.tight_layout()
    save_fig(fig, 'fig6_penrose')


def fig7_perturbations(model, sol):
    """Figure 7 — Perturbation potentials."""
    fig, ax = plt.subplots(figsize=(8, 6))
    r = sol['r']
    m = sol['metric']
    r_star = _compute_tortoise_coordinate(r, m)
    V_RW, V_RW_schw = _compute_regge_wheeler_potential(r, r_star, m, l=2)

    # Only plot the exterior region where the potential is physically meaningful
    # (r > 2M, where the potential has the standard barrier shape)
    valid = np.isfinite(V_RW) & np.isfinite(V_RW_schw) & np.isfinite(r_star)
    valid &= (r > 2.0) & (r < 20.0)
    valid &= (V_RW > -1e10) & (V_RW < 1e10)

    ax.plot(r_star[valid], V_RW[valid], 'r-', linewidth=2,
            label=r'$V_{\rm RW}$ (TEP)')
    ax.plot(r_star[valid], V_RW_schw[valid], 'b--', linewidth=2,
            label=r'$V_{\rm RW}$ (Schwarzschild)')
    ax.set_xlabel(r'$r_*$ (tortoise coordinate)')
    ax.set_ylabel(r'$V_{\rm RW}$')
    ax.set_title('Regge-Wheeler potential (no reflective cavity)')
    ax.legend(fontsize=10)
    ax.set_ylim(-0.05, 0.25)

    plt.tight_layout()
    save_fig(fig, 'fig7_perturbations')


def fig8_observable_shifts(model, sol):
    """Figure 8 — Observable shifts across TEP parameter domain."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Vary phi_0 and compute key diagnostics
    phi_0_values = np.linspace(1.5, 4.0, 20)
    rt_ratios = []
    K_values = []
    ell_values = []

    for phi_0 in phi_0_values:
        m_test = TEPBHModel(
            beta_A=-1.0, B0=1.0, n_B=2.0,
            phi_0=phi_0, delta=0.05, M=1.0,
            sigma_B=1.5,
        )
        sol_test = solve_tep_bh(m_test, r_min=1e-6, r_max=50.0, n_points=5000)
        r_t = sol_test.get('r_temporal_horizon')
        if r_t is not None:
            rt_ratios.append(r_t / m_test.r_h)
        else:
            rt_ratios.append(np.nan)

        r_test = sol_test['r']
        K_test = sol_test['curvature']['Kretschmann']
        idx_01 = np.argmin(np.abs(r_test - 0.1))
        K_values.append(K_test[idx_01])
        ell_values.append(sol_test['geodesics']['ell_at_horizon'])

    # Plot 1: temporal horizon location vs phi_0
    ax = axes[0]
    ax.plot(phi_0_values, rt_ratios, 'ro-')
    ax.set_xlabel(r'$\phi_0$ (scalar amplitude)')
    ax.set_ylabel(r'$r_t / r_{\rm Schw}$')
    ax.set_title('Temporal horizon location')
    ax.set_ylim(0, 1)

    # Plot 2: Kretschmann at r=0.1M vs phi_0
    ax = axes[1]
    ax.plot(phi_0_values, K_values, 'bs-')
    ax.set_xlabel(r'$\phi_0$')
    ax.set_ylabel(r'$\tilde{K}$ at $r = 0.1M$')
    ax.set_title('Curvature regularization')
    ax.set_yscale('log')

    # Plot 3: Physical distance at horizon vs phi_0
    ax = axes[2]
    ax.plot(phi_0_values, ell_values, 'g^-')
    ax.set_xlabel(r'$\phi_0$')
    ax.set_ylabel(r'$\tilde{\ell}$ at $r = 2M$')
    ax.set_title('Physical distance at horizon')

    plt.tight_layout()
    save_fig(fig, 'fig8_observable_shifts')


def fig9_accretion(model, sol):
    """Figure 9 — Accretion dynamics: circular geodesics and epicyclic frequencies."""
    r = sol['r']
    metric = sol['metric']
    M = model.M

    geodesics = _compute_circular_geodesics(r, M)
    epicyclic = _compute_epicyclic_frequencies(r, M)
    isco = _compute_isco_properties(M)
    redshift = _compute_redshift_profile(r, M)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Specific energy and angular momentum
    ax = axes[0, 0]
    r_plot = r[r > 3.0 * M]
    E = geodesics['E'][r > 3.0 * M]
    L = geodesics['L'][r > 3.0 * M]
    ax.plot(r_plot, E, 'b-', linewidth=2, label=r'$E$ (specific energy)')
    ax.axhline(isco['E_isco'], color='b', linestyle='--', alpha=0.5,
               label=r'$E_{\rm ISCO} = \sqrt{8/9}$')
    ax2 = ax.twinx()
    ax2.plot(r_plot, L, 'r-', linewidth=2, label=r'$L$ (specific ang. mom.)')
    ax2.axhline(isco['L_isco'], color='r', linestyle='--', alpha=0.5,
                label=r'$L_{\rm ISCO} = 2\sqrt{3}\,M$')
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'$E$', color='b')
    ax2.set_ylabel(r'$L / M$', color='r')
    ax.set_title('Circular orbit energy and angular momentum')
    ax.set_xlim(3, 20)
    ax.set_ylim(0.85, 1.05)
    ax2.set_ylim(3.0, 5.0)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='center right')
    ax.axvline(6.0, color='gray', linestyle=':', alpha=0.5, label='ISCO')

    # Plot 2: Epicyclic frequencies
    ax = axes[0, 1]
    valid = (r > 3.0 * M) & (r < 30.0)
    ax.plot(r[valid], epicyclic['Omega_r'][valid], 'b-', linewidth=2,
            label=r'$\Omega_r$ (radial)')
    ax.plot(r[valid], epicyclic['Omega_theta'][valid], 'r--', linewidth=2,
            label=r'$\Omega_\theta$ (vertical)')
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(6.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'$\Omega$')
    ax.set_title('Epicyclic frequencies')
    ax.legend(fontsize=10)
    ax.set_xlim(3, 30)
    ax.set_ylim(-0.01, 0.08)
    ax.text(7, 0.06, 'ISCO', fontsize=10, color='gray')

    # Plot 3: Redshift profile
    ax = axes[1, 0]
    valid = (r > 2.0 * M) & (r < 30.0)
    ax.plot(r[valid], redshift['g_redshift'][valid], 'b-', linewidth=2,
            label=r'$g = \sqrt{1 - 2M/r}$')
    ax.axvline(6.0, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(isco['redshift_factor_isco'], color='b', linestyle='--', alpha=0.5,
               label=r'$g_{\rm ISCO} = \sqrt{2/3}$')
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'$g$ (redshift factor)')
    ax.set_title('Gravitational redshift from circular orbit')
    ax.legend(fontsize=10)
    ax.set_xlim(2, 30)
    ax.set_ylim(0, 1.1)

    # Plot 4: QPO frequency ratio
    ax = axes[1, 1]
    valid = (r > 6.0 * M) & (r < 30.0) & np.isfinite(epicyclic['ratio_r_to_theta'])
    ax.plot(r[valid], epicyclic['ratio_r_to_theta'][valid], 'g-', linewidth=2,
            label=r'$\Omega_r / \Omega_\theta$')
    ax.axhline(2.0/3.0, color='r', linestyle='--', alpha=0.7,
               label=r'$2/3$ (3:2 resonance)')
    # Find the 3:2 resonance point
    ratios = epicyclic['ratio_r_to_theta'][valid]
    r_vals = r[valid]
    idx_32 = np.argmin(np.abs(ratios - 2.0/3.0))
    r_32 = r_vals[idx_32]
    ax.axvline(r_32, color='r', linestyle=':', alpha=0.5)
    ax.plot(r_32, ratios[idx_32], 'ro', markersize=8)
    ax.text(r_32 + 0.5, 0.3, f'r = {r_32:.1f}M', fontsize=10, color='r')
    ax.set_xlabel(r'$r / M$')
    ax.set_ylabel(r'$\Omega_r / \Omega_\theta$')
    ax.set_title('QPO frequency ratio (3:2 resonance)')
    ax.legend(fontsize=10)
    ax.set_xlim(6, 30)
    ax.set_ylim(0, 1.1)

    plt.tight_layout()
    save_fig(fig, 'fig9_accretion')


def fig10_observational_constraints(model, sol):
    """Figure 10 — Observational constraints: EHT shadows and LIGO QNMs."""
    import json
    import os

    PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
    eht_path = PROCESSED_DIR / "eht_measurements.json"
    ligo_path = PROCESSED_DIR / "ligo_qnm_measurements.json"

    if not eht_path.exists() or not ligo_path.exists():
        print_status("  Skipping fig10: measurement data not found (run step_00)", "WARNING")
        return

    with open(eht_path) as f:
        eht_data = json.load(f)
    with open(ligo_path) as f:
        ligo_data = json.load(f)

    # Load step_05 results
    step05_path = PROJECT_ROOT / "results" / "step_05_observational_constraints.json"
    if step05_path.exists():
        with open(step05_path) as f:
            step05 = json.load(f)
    else:
        print_status("  Skipping fig10: step_05 results not found", "WARNING")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Plot 1: EHT shadow comparison
    ax = axes[0]
    sources = ['M87*', 'Sgr A*']
    measured = [step05['M87']['shadow_measured_uas'],
                step05['SgrA']['shadow_measured_uas']]
    measured_err = [step05['M87']['shadow_uncertainty_uas'],
                    step05['SgrA']['shadow_uncertainty_uas']]
    tep_pred = [step05['M87']['shadow_tep_uas'],
                step05['SgrA']['shadow_tep_uas']]
    schw_pred = [step05['M87']['shadow_schw_uas'],
                 step05['SgrA']['shadow_schw_uas']]

    x = np.arange(len(sources))
    width = 0.25
    ax.bar(x - width, measured, width, yerr=measured_err, label='Measured',
           color='gray', capsize=5, edgecolor='black')
    ax.bar(x, tep_pred, width, label='TEP prediction', color='red', alpha=0.7)
    ax.bar(x + width, schw_pred, width, label='Schwarzschild', color='blue', alpha=0.7)
    ax.set_ylabel(r'Shadow diameter ($\mu$as)')
    ax.set_title('EHT shadow diameter comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(sources)
    ax.legend(fontsize=9)

    # Plot 2: LIGO QNM frequency comparison
    ax = axes[1]
    gw = step05['GW150914']
    f_tep = gw['f220_tep_Hz']
    f_schw = gw['f220_schw_Hz']
    f_kerr = gw['f220_kerr_Hz']
    f_meas = gw['f220_measured_Hz']
    labels = [f'Measured\n({f_meas:.0f}±{gw["f220_uncertainty_Hz"]:.0f} Hz)',
              f'TEP\n({f_tep:.1f} Hz)',
              f'Schwarzschild\n({f_schw:.1f} Hz)',
              f'Kerr\n({f_kerr:.0f} Hz)']
    values = [gw['f220_measured_Hz'], gw['f220_tep_Hz'],
              gw['f220_schw_Hz'], gw['f220_kerr_Hz']]
    errors = [gw['f220_uncertainty_Hz'], 0, 0, 0]
    colors = ['gray', 'red', 'blue', 'green']
    x = np.arange(len(labels))
    ax.bar(x, values, yerr=errors, color=colors, alpha=0.7, capsize=5, edgecolor='black')
    ax.set_ylabel(r'$f_{220}$ (Hz)')
    ax.set_title('GW150914 ringdown frequency')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)

    # Plot 3: Chi-squared comparison
    ax = axes[2]
    tests = ['M87\nshadow', 'Sgr A*\nshadow', 'GW150914\nf220', 'GW150914\ntau220']
    chi2_tep = [step05['M87']['chi2_tep'], step05['SgrA']['chi2_tep'],
                gw['f220_chi2_tep'], gw['tau220_chi2_tep']]
    chi2_schw = [step05['M87']['chi2_schw'], step05['SgrA']['chi2_schw'],
                 gw['f220_chi2_schw'], gw['tau220_chi2_schw']]
    chi2_kerr = [step05['M87']['chi2_schw'], step05['SgrA']['chi2_schw'],
                 gw['f220_chi2_kerr'], gw['tau220_chi2_kerr']]

    x = np.arange(len(tests))
    width = 0.25
    ax.bar(x - width, chi2_tep, width, label='TEP', color='red', alpha=0.7)
    ax.bar(x, chi2_schw, width, label='Schwarzschild', color='blue', alpha=0.7)
    ax.bar(x + width, chi2_kerr, width, label='Kerr', color='green', alpha=0.7)
    ax.axhline(4.0, color='gray', linestyle='--', alpha=0.5, label='2σ threshold')
    ax.set_ylabel(r'$\chi^2$')
    ax.set_title('Constraint comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(tests, fontsize=9)
    ax.legend(fontsize=9)
    ax.set_yscale('log')

    plt.tight_layout()
    save_fig(fig, 'fig10_observational_constraints')


def main():
    ensure_dirs()
    logger = make_step_logger(STEP_ID)

    print_status("TEP-BH Figure Generation", "TITLE")
    print_status(f"Step ID: {STEP_ID}", "INFO")
    print_status(f"Timestamp: {datetime.now().isoformat()}", "INFO")
    print_status(f"Output directory: {FIG_DIR}", "INFO")
    print_status("")

    model, sol = get_solution()

    print_status("Figure 1: Two-metric causal architecture...", "PROCESS")
    fig1_two_metric(model, sol)

    print_status("Figure 2: Geometric vs physical radial distance...", "PROCESS")
    fig2_radial_distance(model, sol)

    print_status("Figure 3: Clock projection...", "PROCESS")
    fig3_clock_projection(model, sol)

    print_status("Figure 4: Curvature invariants...", "PROCESS")
    fig4_curvature(model, sol)

    print_status("Figure 5: Physical volume and density...", "PROCESS")
    fig5_volume_density(model, sol)

    print_status("Figure 6: Penrose diagrams...", "PROCESS")
    fig6_penrose(model, sol)

    print_status("Figure 7: Perturbation potentials...", "PROCESS")
    fig7_perturbations(model, sol)

    print_status("Figure 8: Observable shifts...", "PROCESS")
    fig8_observable_shifts(model, sol)

    print_status("Figure 9: Accretion dynamics...", "PROCESS")
    fig9_accretion(model, sol)

    print_status("Figure 10: Observational constraints...", "PROCESS")
    fig10_observational_constraints(model, sol)

    print_status(f"All 10 figures saved to {rel(FIG_DIR)}", "SUCCESS")


if __name__ == '__main__':
    main()
