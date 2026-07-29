# Temporal Equivalence Principle: Black Holes and the Temporal Horizon

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Status:** Preprint — Paper 28 (Bahrain), v0.1

## Abstract

The Temporal Equivalence Principle (TEP) asks: *can the phenomena attributed to black holes arise from an extreme relative-time gradient without physical collapse into an ultra-dense object?* The framework distinguishes the Einstein-frame gravitational metric $g_{\mu\nu}$ from the universal causal matter metric
$$\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}+B(\phi)\nabla_\mu\phi\nabla_\nu\phi,$$
so that matter, photons, and ideal clocks propagate on $\tilde g_{\mu\nu}$ while tensor and scalar characteristics are derived from the principal symbol of the coupled gravitational–scalar system.

The mass conventionally assigned to an astrophysical black-hole candidate is not directly weighed. It is reconstructed from observed angular positions, spectral shifts, signal periods and image scales under the Isochrony Axiom — the assumption that source clocks, photon propagation and observer clocks can be mapped onto a single universal time coordinate. TEP defines a strong-field Phantom Mass residual $M_{\rm phantom}^{T} \equiv M_{\rm fit}^{\rm GR} - M_{\rm matter}^{\rm TEP}$ and proposes to determine its sign and magnitude through a non-isochronous refit. The sign is not assumed: slow deep clocks alone deflate the inferred mass; positive Phantom Mass requires spatial magnification to dominate temporal stretching. This is the strong-field counterpart of Phantom Mass at galactic scales. A "black hole" is a *temporal well* — a regular spatial region in which the rate of proper time differs radically from the exterior, with physical lapse $N(r) \geq N_{\min} > 0$ everywhere, finite but extremely small at the centre. *Cosmological expansion and black-hole collapse are dual misinterpretations of dynamical proper time.*

The fixed-background conformal theorem proves that the temporal field must backreact on the geometric metric: holding $g_{\mu\nu}$ fixed to Schwarzschild makes finite curvature ($\phi_0 \geq 3/2$) and bounded areal radius ($\phi_0 \leq 1$) mutually exclusive. The inverse reconstruction identifies $w_r = -1$ as a sufficient regularity target, demonstrated by the Hayward benchmark and now realised by the regularised sGB-coupled interior solver. The sGB coupling serves as a low-energy effective field theory for the exterior, providing conditional observable shifts at fixed GR-calibrated mass: shadow $-0.044\%$ (photons, $\mathcal{O}(\eta^2)$), ISCO on $\tilde g$ $+1.95\%$ (massive particles, $\mathcal{O}(\eta)$) at $\eta = -0.1$. The coupling-order difference between the photon and massive-particle sectors reflects conformal invariance of null geodesics in the sGB benchmark (Section 6.2): both sectors are governed by the same unified temporal field, and the difference arises because null geodesics are conformally invariant while timelike geodesics feel the conformal factor $A = e^{-\phi} > 1$. The tensor speed $c_T = 1$ on the Schwarzschild background; the full tensor characteristic metric on the TEP solution requires the complete coupled perturbation system.

The decisive analysis is a *TEP-native refit of the raw observations*, not a perturbation of a pre-assumed Schwarzschild mass. We demonstrate the pipeline on the 25-year astrometric and spectroscopic record of the S2 star orbiting Sgr A*, using GRAVITY/VLT and Keck data through a 7-step non-isochronous inference. At S2 scales the mass-inflation branch is selected with $100\%$ positive Phantom Mass posterior, while the coupling is consistent with zero at current precision, showing TEP reduces to GR in the weak field. At horizon scale, a direct visibility-domain fit to $46{,}846$ binned EHT M87$^\ast$ and Sgr A$^\ast$ long-baseline amplitudes excludes the extreme horizonless threshold at $>2\sigma$, bounding the temporal well's regularisation scale to the sub-critical regime ($g < g_{\rm crit}$). The full observational replacement requires the nonlinear TEP solution and the raw-observable refit.

## Overview

TEP-BH (Paper 28, Bahrain) derives the causal completion of gravitational collapse in the Temporal Equivalence Principle framework. The geometric Schwarzschild endpoint is not the physical endpoint of matter propagation. The complete conformal-disformal matter metric is globally Lorentzian and nondegenerate, with the conformal factor diverging as $A \to \infty$ and the disformal shear vanishing as $B \to 0$ in the deep interior. The center is continuous, curvature-free regular space with diverging areal radius, not a spatial singularity.

The paper provides:

- **Spherical ansatz and field equations** for the strong-field scalar configuration
- **Causal regularity and completeness** of the physical metric $\tilde g_{\mu\nu}$ at $r\to 0$
- **The density illusion**: physical volume element does not collapse as $r^3$
- **Singularity theorem reinterpretation** in a two-metric framework
- **Horizon thermodynamics and information** without a physical singularity
- **Comparative anatomy** distinguishing TEP from gravastars, temporal wells, fuzzballs and firewalls
- **Linear stability and ringdown** signatures
- **Electromagnetic propagation and black-hole imaging** signatures
- **Accretion, ISCO and redshift transfer** in the disformal causal metric
- **Minimal empirical test programme**

## Manuscript Sections

1. Introduction: Which Conclusions Follow from a Two-Metric Construction?
2. TEP Foundations and Frame Dictionary
3. Spherical Geometry and the Fixed-Background Theorem
4. Coupled Solution: Einstein–Scalar–Gauss–Bonnet Exterior
5. Global Geometry: Benchmark and Validation
6. Causal Structure and the Temporal Horizon
7. Perturbation Analysis: Three-Channel Ringdown and Tensor Speed
8. Observable Predictions
9. Interpretation
10. Conclusion

Appendices A–L: Conventions and Disformal Identities, Reduced Field Equations, Strong-Field Asymptotic Expansion, Curvature Invariants and the Fixed-Background Theorem, Geodesic Structure, Physical Volume and Density Diagnostics, Scope of Gravitational Perturbations, Ray Tracing and Orbital Mechanics, Comparison Table, Reproducibility, Hayward Benchmark, Explicit Derivation of Corrected Exterior Observables.

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-BH.git
cd TEP-BH

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Site Build

```bash
cd site
npm run build
```

This generates:
- `site/dist/index.html` — static manuscript
- `28-TEP-BH-v0.1-Bahrain.md` — root markdown manuscript

## PDF Generation

```bash
python scripts/utils/generate_site_pdf.py --quality high --wait-time 8
```

Generates `28-TEP-BH-v0.1-Bahrain.pdf` in both the root and `site/public/docs/`.

## Project Structure

```
TEP-BH/
├── core/                 # Physics modules
├── scripts/              # Pipeline and utility scripts
│   ├── run_pipeline.py
│   ├── steps/            # Individual pipeline steps
│   └── utils/            # Utilities (PDF generation, PDF processing, logging)
├── results/              # Pipeline outputs
├── site/                 # Manuscript site
│   ├── components/       # HTML components (edit these)
│   ├── dist/             # Built static site (auto-generated)
│   ├── public/           # Static assets, PDF, sitemap
│   ├── build.js          # Site build script
│   └── manifest.json     # Site manifest
├── manuscripts/          # Markdown manuscripts
├── 28-TEP-BH-v0.1-Bahrain.md  # Auto-generated root manuscript
└── 28-TEP-BH-v0.1-Bahrain.pdf # Generated PDF
```

## Citation

If using this work, please cite:

```bibtex
@software{tep_bh_2026,
  author       = {Matthew Lukin Smawfield},
  title        = {Temporal Equivalence Principle: Black Holes and the Temporal Horizon},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1 (Bahrain)},
  doi          = {10.5281/zenodo.xxxxxxxx},
  url = {https://mlsmawfield.com/tep/bh}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Related Papers

- **TEP-C0 (Paper 26, Athens)**: Distance-redshift and supernova evidence
- **TEP-HC (Paper 18, Cambridge)**: Acoustic-sector perturbations via hi_class
- **TEP-BH (Paper 28, Bahrain)**: Black-hole causal completion (this work)

## Contact

For questions or issues, please open a GitHub issue at https://github.com/matthewsmawfield/TEP-BH
