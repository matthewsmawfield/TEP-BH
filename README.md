# Temporal Equivalence Principle: Temporal Horizon Cosmology and the Absence of a Physical Big Bang Singularity

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Status:** Complete — Paper 27 (Thika), v0.2

## Abstract

Standard FLRW cosmology extrapolates observed cosmic expansion backward to $a(t)\to0$, producing a Big Bang singularity at finite proper time. This paper demonstrates that this singularity is a reconstruction artifact of imposing a globally isochronous expanding-frame description on a conformal temporal geometry. In the Temporal Equivalence Principle (TEP), the observational role of FLRW expansion is reconstructed through conformal temporal transport: the effective scale factor $a_{\rm eff}$ arises from accumulated open-path conformal temporal shear along cosmological lines of sight rather than from physical expansion of space. TEP-C0 (Paper 26) established the distance-redshift and supernova evidence and deferred full nonsingular matter-frame closure to a dedicated temporal-horizon analysis; here that closure is delivered.

The Temporal Horizon Cosmology framework is developed here, proving that the apparent $a_{\rm eff}\to0$ limit is not a physical curvature singularity but a temporal horizon. Two distinct projections of the temporal field are required: $A_{\rm clock}(z)=(1+z)^{-1}$ is the exact observational clock/redshift mapping that drives $a_{\rm eff}\to0$ as $z\to\infty$, while $A_{\rm dyn}(z)=\left(1+z/z_{t}\right)^{-\epsilon_{\rm eff}(z)}$ is the dynamically screened shear response that modifies expansion, BBN, recombination, and perturbations only at late times. Proposition 1 establishes curvature regularity of the temporal conformal boundary: for $A_{\rm clock}(\eta)=C\eta^{-p}$ with $0 < p\le\tfrac12$, all polynomial curvature invariants vanish at the boundary, timelike proper time diverges, and null geodesics have divergent affine parameter. Figure 1 (Section 4.5) illustrates the resulting conformal-boundary interpretation: the singular lower edge of standard flat $\Lambda$CDM is replaced by a smooth temporal conformal boundary $\mathscr{T}^{-}$, where $A_{\rm clock}\to0$ and curvature invariants vanish. The conformal compactification is smooth, the Weyl tensor vanishes on the boundary, and every causal curve approaches the regular past boundary $\mathscr{T}^{-}$ rather than terminating at a singularity. The temporal horizon is therefore simultaneously curvature-empty, timelike-complete, and null-complete in this branch.

The effective stress-energy tensor of the temporal field violates the Strong Energy Condition, an explicit prerequisite of the Hawking-Penrose singularity theorems. The thermal screening scale is derived from the transition redshift: $T_{\rm lock}=T_{0}(1+z_{t})\sim 0.03\,\mathrm{eV}$ for $z_{t}\sim 130$ and $T_{0}=2.725\,\mathrm{K}$, giving strong epoch-by-epoch screening ($S_{\rm epoch}\sim 10^{-12}$ at BBN, $\sim 10^{-2}$ at recombination); $T_{\rm lock}$ is not an independent parameter but is fixed once $z_{t}$ is fixed. Screened TEP reproduces the standard BBN successful sector and inherits the standard lithium anomaly. Recombination is computed with the full non-equilibrium Peebles/RECFAST treatment. The temporal-horizon thermal mapping preserves a FIRAS-compatible blackbody with no spectral distortions.

The scalar perturbation spectrum is derived from fluctuations of the clock field, $\zeta=\delta\ln A_{\rm clock}$, yielding a power spectrum $P_{\zeta}(k)\propto k^{n_{s}-1}$ with spectral-flow parameter $n_{s}-1=-2\epsilon_{\rm field}$. The observed Planck value $n_{s}=0.965$ constrains $\epsilon_{\rm field}=0.0175$. Tensor modes are derived directly from the temporal-conformal metric: for $A_{\rm clock}(\eta)\sim\eta^{-p}$ the tensor source term $A_{\rm clock}''/A_{\rm clock}=p(p+1)/\eta^{2}\to 0$ at the horizon, so the tensor equation approaches the Minkowski vacuum. The imported inflationary consistency relation $r=16\epsilon_{\rm field}$ is not assumed. Numerical integration of the native tensor equation across the finite transition profile (Step 09b) yields $r(k_{\rm pivot})=9\times 10^{-6}$ and $r_{\rm max}=6.26\times 10^{-4}$, both well below the BICEP/Keck 2021 bound $r<0.036$; tensor power is controlled only by the finite transition region. CMB anisotropy and LSS observables are reproduced in the screened-limit reduction, inheriting agreement with Planck 2018 and BOSS DR12 by construction rather than as independent empirical confirmation.

The causal matter-frame universe is curvature-regular at the temporal conformal boundary. The apparent Big Bang is a temporal horizon, not a physical curvature singularity. All background and thermal observational pillars are preserved in the screened-limit reduction; the scalar perturbation shape is reproduced, and the tensor-to-scalar ratio is computed from the native temporal-conformal wave equation, yielding values well below observational bounds.

## Overview

TEP-TH (Paper 27, Thika) delivers the full temporal-horizon closure of the Temporal Equivalence Principle framework, proving that the apparent Big Bang singularity is a reconstruction artifact. Building on TEP-C0 (Paper 26, Athens) which established the distance-redshift and supernova evidence, and TEP-HC (Paper 18, Cambridge) which validated the acoustic-sector perturbations via hi_class, TEP-TH provides:

- **Temporal-horizon curvature analysis** proving regularity at $A_{\rm clock}\to 0$ (Proposition 1)
- **Geodesic completeness** demonstrating infinite affine/proper time at the temporal boundary
- **Penrose diagram** (Figure 1) placing $\mathscr{T}^{-}$ on the same rigorous footing as $\mathscr{I}^{+}$
- **BBN nucleosynthesis** with AlterBBN network and epoch-screened TEP expansion
- **Recombination physics** with full non-equilibrium Peebles/RECFAST treatment
- **Temporal-horizon thermal mapping** preserving FIRAS-compatible blackbody
- **Scalar perturbations** derived from clock-field fluctuations with spectral-flow parameter $\epsilon_{\rm field}=0.0175$
- **Tensor perturbations** computed from native temporal-conformal wave equation yielding $r(k_{\rm pivot})=9\times 10^{-6}$
- **CMB anisotropy** (TT, TE, EE) and **LSS observables** validated in screened-limit reduction

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-TH.git
cd TEP-TH

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Initialize external submodules (AlterBBN, CLASS)
git submodule update --init --recursive
```

## Quick Start

```bash
# Run full pipeline
cd TEP-TH
python scripts/run_pipeline.py

# Individual steps (see scripts/README.md)
python scripts/steps/step_00_temporal_horizon_mapping.py
python scripts/steps/step_01_matter_frame_curvature.py
python scripts/steps/step_02_geodesic_completeness.py
```

## Pipeline Execution

### Full Pipeline
```bash
cd TEP-TH
python scripts/run_pipeline.py
```

### Individual Steps
```bash
# Temporal horizon mapping
python scripts/steps/step_00_temporal_horizon_mapping.py

# Matter-frame curvature
python scripts/steps/step_01_matter_frame_curvature.py

# Geodesic completeness
python scripts/steps/step_02_geodesic_completeness.py

# BBN with AlterBBN
python scripts/steps/step_04_full_bbn_abundances.py

# Recombination
python scripts/steps/step_05_recombination_visibility.py

# Scalar perturbations
python scripts/steps/step_08_primordial_perturbation_boundary.py

# Tensor perturbations
python scripts/steps/step_09b_native_tensor_integration.py
```

## Data Sources

All data is downloaded from public repositories:

- **BBN data**: [BBN Frontiers 2020 Review](https://arxiv.org/abs/2005.14167)
  - Light-element abundance constraints
  - Nuclear reaction rates

- **FIRAS CMB**: [NASA LAMBDA](https://lambda.gsfc.nasa.gov)
  - COBE/FIRAS monopole spectrum
  - Perfect blackbody validation

- **CMB/LSS**: Planck 2018 and BOSS DR12 public releases
  - Used for screened-limit validation

## Physics Validation

| Observable | Computed | Target | Status |
|------------|----------|--------|--------|
| $\tilde{\mathcal{K}}$ at horizon | 0 | 0 (regular) | ✅ |
| Timelike proper time | $\infty$ | $\infty$ (complete) | ✅ |
| Null affine parameter | $\infty$ | $\infty$ (complete) | ✅ |
| $Y_p$ | 0.247 | Planck 0.247 | ✅ |
| D/H | $2.51\times 10^{-5}$ | PDG $2.6\times 10^{-5}$ | ✅ |
| $n_s$ | 0.965 | Planck 0.965 | ✅ |
| $r(k_{\rm pivot})$ | $9\times 10^{-6}$ | BICEP/Keck $<0.036$ | ✅ |
| $r_{\rm max}$ | $6.26\times 10^{-4}$ | BICEP/Keck $<0.036$ | ✅ |

## Project Structure

```
TEP-TH/
├── core/                 # Physics modules
│   ├── conformal_scaling.py  # Conformal clock field
│   ├── constants.py         # Physical constants
│   └── cosmology.py          # Background cosmology
├── scripts/              # Pipeline scripts
│   ├── steps/               # Individual pipeline steps
│   │   ├── step_00_temporal_horizon_mapping.py
│   │   ├── step_01_matter_frame_curvature.py
│   │   ├── step_02_geodesic_completeness.py
│   │   └── ...
│   └── utils/               # Utilities
├── external/              # External codes
│   ├── AlterBBN/            # BBN network
│   └── class/               # CLASS Boltzmann code
├── data/                 # Data directory
│   ├── raw/                 # Downloaded datasets
│   └── processed/          # Intermediate outputs
├── results/              # Pipeline outputs
│   ├── figures/             # Generated plots
│   └── outputs/             # JSON results
├── site/                 # Manuscript site
│   ├── components/          # HTML components
│   └── public/              # Static assets
└── manuscripts/          # Markdown manuscripts
```

## Testing

```bash
# Run individual step with verbose logging
python scripts/steps/step_01_matter_frame_curvature.py --verbose

# Check results
cat results/step_01_matter_frame_curvature.json
```

## Methodology

### Temporal-Horizon Mapping

The TEP framework distinguishes two projections of the temporal field:

```
A_clock(z) = (1 + z)^(-1)  # Exact observational clock/redshift mapping
A_dyn(z) = (1 + z/z_t)^(-epsilon_eff(z))  # Dynamically screened shear response
epsilon_eff(z) = epsilon_dyn * S_epoch(T)
S_epoch(T) = 1 / (1 + (T_eV/T_lock)^n_epoch)
```

Where:
- `A_clock(z)`: drives a_eff → 0 as z → ∞ (temporal horizon)
- `A_dyn(z)`: modifies expansion, BBN, recombination, perturbations at late times
- `epsilon_dyn`: baseline shear amplitude
- `z_t`: transition redshift (=100)
- `T_lock = 0.03` eV: thermal screening threshold
- `S_epoch(T)`: epoch-screening function (S → 0 in hot early universe, S → 1 at late times)

### Curvature Regularity

For the temporal-horizon profile $A_{\rm clock}(\eta)=C\eta^{-p}$ with $0<p\le\tfrac12$:
- All polynomial curvature invariants vanish at the boundary
- Timelike proper time diverges for $0<p\le 1$
- Null affine parameter diverges for $0<p\le\tfrac12$
- The boundary is a regular conformal-temporal endpoint

### BBN Computation

- **Network**: AlterBBN (full nuclear network)
- **Physics**: TEP-modified expansion rate with epoch screening
- **Abundances**: $Y_p$, D/H, $^3$He/H, $^7$Li/H, $N_{\rm eff}$

### Recombination

- **Treatment**: Full non-equilibrium Peebles/RECFAST
- **Thermal mapping**: Temporal-horizon conformal transformation preserves blackbody
- **Validation**: $x_e(z)$, $g(z)$, $r_s$, $\theta_s$

### Perturbations

- **Scalar**: Derived from clock-field fluctuations $\zeta=\delta\ln A_{\rm clock}$
- **Spectral flow**: $n_s-1=-2\epsilon_{\rm field}$ with $\epsilon_{\rm field}=0.0175$
- **Tensor**: Native temporal-conformal wave equation, source term $\to 0$ at horizon
- **Results**: $r(k_{\rm pivot})=9\times 10^{-6}$, $r_{\rm max}=6.26\times 10^{-4}$

## Citation

If using this work, please cite:

```bibtex
@software{tep_th_2026,
  author       = {Matthew Lukin Smawfield},
  title        = {Temporal Equivalence Principle: Temporal Horizon Cosmology and the Absence of a Physical Big Bang Singularity},
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1 (Thika)},
  doi          = {10.5281/zenodo.20723059},
  url = {https://mlsmawfield.com/tep/th}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Status

- ✅ Temporal-horizon curvature analysis: Complete
- ✅ Geodesic completeness: Complete
- ✅ Penrose diagram: Complete
- ✅ BBN: Complete (AlterBBN network)
- ✅ Recombination: Complete (Peebles/RECFAST)
- ✅ Scalar perturbations: Complete
- ✅ Tensor perturbations: Complete (native wave equation)
- ✅ CMB/LSS validation: Complete (screened-limit reduction)

## Related Papers

- **TEP-C0 (Paper 26, Athens)**: Distance-redshift and supernova evidence
- **TEP-HC (Paper 18, Cambridge)**: Acoustic-sector perturbations via hi_class
- **TEP-TH (Paper 27, Thika)**: Temporal-horizon closure (this work)

## Contact

For questions or issues, please open a GitHub issue at https://github.com/matthewsmawfield/TEP-TH