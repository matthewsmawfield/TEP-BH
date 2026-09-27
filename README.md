# Temporal Equivalence Principle: Black Holes and the Temporal Horizon

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21677826.svg)](https://doi.org/10.5281/zenodo.21677826)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](LICENSE)

![TEP-BH: Black Holes and the Temporal Horizon](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.3 (Bahrain)  
**First published:** 29 July 2026 · **Last updated:** 23 September 2026
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.21677826](https://doi.org/10.5281/zenodo.21677826)  
**Website:** [https://mlsmawfield.com/tep/bh](https://mlsmawfield.com/tep/bh)  
**Paper Series:** TEP Series: Paper 28 (Black Holes and the Temporal Horizon)

## Abstract

The Temporal Equivalence Principle treats proper time as a dynamical field. Matter, light, and ideal clocks couple to the universal causal metric $\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}+B(\phi)\nabla_\mu\phi\nabla_\nu\phi$. The Einstein-frame metric carries the gravitational dynamics; tensor propagation is determined by the principal symbol of the coupled temporal–geometric equations. This paper develops the strong-field consequence: Under TEP, a black hole is modeled as a temporal well—a target regular spatial region in which matter-frame clock rates become strongly suppressed. A continuous, extreme but finite gradient in proper-time rate provides a unified mechanism for the principal observational signatures attributed to a black hole.

Darkness, apparent compactness, large inferred mass, and practical inaccessibility are exterior reconstructions of temporal decoupling. The operational boundary is the Temporal Horizon — the observer-relative threshold beyond which the clock-transfer factor renders signals practically undetectable — not the event horizon. Inferred gravitational mass and local material mass need not coincide; a strong-field phantom mass residual $M_{\rm phantom}^{T} \equiv M_{\rm fit}^{\rm GR} - M_{\rm matter}^{\rm TEP}$ measures this discrepancy, its sign determined by the data. Standard black-hole ontology assumes isochrony: source clocks, photon propagation, and observer clocks mapped onto a single general-relativistic time coordinate. TEP drops that closure. Under TEP, the conventional reconstruction — compact mass, event horizon, singular collapse — is no longer the unique reading of the same observations.

The temporal field cannot sit passively on fixed Schwarzschild geometry: finite curvature and bounded areal radius are mutually exclusive when $g_{\mu\nu}$ is held fixed. Strong temporal structure forces gravitational backreaction. The canonical TEP matter coupling fixes the temporal sector; the leading curvature operator supplies backreaction; the regularising nonlinear coefficients are selected by global regularity and observation. The construction programme is solving the global field equations of the fixed EFT architecture — not selecting among competing theories.

Four observational consequences follow from one temporal field: time transfer, mass inference, photon accessibility, and ringdown. Weak-field data recover GR; horizon-scale images constrain the photon-region geometry but do not directly establish an event horizon; gravitational-wave ringdown provides the sharpest test, because TEP replaces the purely ingoing event-horizon condition with propagation through a regular temporal domain, changing the late-time spectral problem. Cosmological expansion and black-hole collapse are dual misreadings of dynamical proper time.

Keywords: Temporal Equivalence Principle, temporal well, temporal horizon, dynamical proper time, phantom mass, isochrony axiom, conformal-disformal metric, Schwarzschild incompatibility, apparent compactness, black holes, modified gravity, temporal shear, black-hole observations

# 1. Introduction: Black Holes under the Temporal Equivalence Principle

## Overview

TEP-BH (Paper 28, Bahrain) develops the strong-field consequence of the Temporal Equivalence Principle. Under TEP, a black hole is a temporal well: a regular spatial region with an extreme gradient in relative proper-time accumulation. The standard strong-field template — event horizon, singularity, ultradense core — is a reconstruction from observational data under an implicit isochronous transfer model. Once proper time is treated as dynamical, the conventional reconstruction is no longer the unique reading of the same observations.

The paper provides:

- **TEP theory**: two-metric action, conformal-disformal matter metric, frame dictionary, isochrony axiom identified as the closure TEP replaces
- **Schwarzschild incompatibility**: finite curvature and bounded areal radius are mutually exclusive when the geometric metric is held fixed — strong temporal structure forces gravitational backreaction
- **Existence demonstrations**: regular geometry (Hayward class) and curvature-coupled scalar (sGB) as attainability proofs
- **Minimal temporal-well criteria**: $0 < N(r)$, $N_{\min} \ll N_o$, no observer-independent one-way boundary
- **The Temporal Horizon**: observer-relative accessibility threshold, not a null boundary
- **Unified phenomenology**: redshift, darkness, mass inference, photon region, and ringdown from one temporal field
- **Data confrontation**: S2 as weak-field consistency, EHT as photon-region constraint, gravitational waves as sharpest test
- **The construction programme**: TEP-selected global solution, coupled characteristics, raw non-isochronous refit — next stage, not precondition

## Manuscript Sections

1. Introduction: Black Holes under the Temporal Equivalence Principle
2. TEP Theory: Action, Matter Metric, and Frame Dictionary
3. Why Schwarzschild Is Not Enough
4. Existence Demonstrations
5. Global Geometry and Curvature Regularity
6. Causal Structure
7. Perturbation Analysis: Ringdown Structure
8. Four Consequences of One Principle
9. Interpretation and Scope
10. Conclusion

Appendices A–L: Conventions and Disformal Identities, Reduced Field Equations, Strong-Field Asymptotic Expansion, Curvature Invariants and the Schwarzschild Incompatibility, Geodesic Structure, Physical Volume and Density Diagnostics, Fixed-Schwarzschild Perturbation Baseline, Ray Tracing and Orbital Mechanics, Comparison Table, Reproducibility, Regular-Geometry Validation Benchmark, Explicit Derivation of Corrected Exterior Observables.

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
- `28-TEP-BH-v0.3-Bahrain.md` — root markdown manuscript

## PDF Generation

```bash
python scripts/utils/generate_site_pdf.py --quality high --wait-time 8
```

Generates `28-TEP-BH-v0.3-Bahrain.pdf` in both the root and `site/public/docs/`.

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
├── 28-TEP-BH-v0.3-Bahrain.md  # Auto-generated root manuscript
└── 28-TEP-BH-v0.3-Bahrain.pdf # Generated PDF
```

## License

Creative Commons Attribution 4.0 International License (CC BY 4.0) — see [LICENSE](LICENSE) file.

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 14** | [TEP-GNSS-MGEX](https://github.com/matthewsmawfield/TEP-GNSS-MGEX) | MGEX Multi-GNSS Clock Replication | — |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | [TEP-J0437](https://github.com/matthewsmawfield/TEP-J0437) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |
| **Paper 18** | [TEP-HC](https://github.com/matthewsmawfield/TEP-HC) | EFT Mapping and Acoustic Peak Constraints via hi_class | [10.5281/zenodo.20572722](https://doi.org/10.5281/zenodo.20572722) |
| **Paper 19** | [TEP-LENS](https://github.com/matthewsmawfield/TEP-LENS) | Blind-Prediction Residual Test in Multiply-Imaged Supernovae | — |
| **Paper 26** | [TEP-C0](https://github.com/matthewsmawfield/TEP-C0) | A Covariant Alternative to Cosmic Expansion | [10.5281/zenodo.20370143](https://doi.org/10.5281/zenodo.20370143) |
| **Paper 28** | **TEP-BH** (This repo) | Black Holes and the Temporal Horizon | [10.5281/zenodo.21677826](https://doi.org/10.5281/zenodo.21677826) |

## Citation

```bibtex
@article{smawfield2026bh,
  title={Temporal Equivalence Principle: Black Holes and the Temporal Horizon},
  author={Smawfield, Matthew Lukin},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.21677826},
  note={Preprint v0.3 (Bahrain)},
  license={CC-BY-4.0}
}
```

---

## Open Science Statement

These are working preprints shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under a Creative Commons Attribution 4.0 International license (CC BY 4.0) to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)

