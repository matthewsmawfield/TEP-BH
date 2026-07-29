# TEP-BH Pipeline

Black-hole analysis pipeline for the Temporal Equivalence Principle (TEP).

## Quick Start

```bash
cd "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-BH"
python scripts/run_pipeline.py
```

## Usage

```bash
# Full pipeline (everything: core, derive, interior, inference, figures)
python scripts/run_pipeline.py

# Run specific phases only
python scripts/run_pipeline.py --no-derive --no-inference

# Run specific steps
python scripts/run_pipeline.py --start-step step_01 --stop-step step_03

# Skip steps
python scripts/run_pipeline.py --skip-steps step_02

# Stop after a phase
python scripts/run_pipeline.py --stop-step derive

# List all steps
python scripts/run_pipeline.py --list-steps

# Continue on error
python scripts/run_pipeline.py --continue-on-error
```

## Pipeline Steps

All 50 steps are in `scripts/steps/`, numbered `step_00` through `step_49`.

### Core — Data & Fixed-Background (00–11)

| Step ID | Description | Data Source |
|---------|-------------|-------------|
| `step_00_data_download` | Download EHT M87* visibilities, compile measurement tables | EHT GitHub, published papers |
| `step_01_field_equations` | Disformal metric, curvature invariants, geodesic completeness | Numerical solver |
| `step_02_perturbations` | Regge-Wheeler/Zerilli potentials, QNMs, echo check | Numerical solver |
| `step_03_raytracing` | Photon sphere, shadow, ISCO, EHT comparison | Numerical solver + EHT |
| `step_04_accretion` | Circular geodesics, epicyclic frequencies, QPO | Numerical solver |
| `step_05_observational_constraints` | Chi-squared vs EHT and LIGO | `data/processed/` |
| `step_06_gw190521_mass_gap` | GW190521 exterior consistency | Zenodo 4057131 |
| `step_07_qpo_frequency_lock` | QPO 3:2 frequency ratios | Strohmayer 2001, etc. |
| `step_08_jwst_early_smbhs` | JWST early SMBHs growth | Harikane et al. 2023 |
| `step_09_spin_bias` | Spin distribution analysis | McClintock et al. 2014 |
| `step_10_tde_missing_flares` | TDE sub-Eddington luminosities | van Velzen et al. 2021 |
| `step_11_eht_polarization` | EHT polarization birefringence check | EHT 2021, 2024 |

### Core — Coupled Solution (12–15)

| Step ID | Description |
|---------|-------------|
| `step_12_self_gravitating` | sGB metric corrections, scalar profile, exterior observables |
| `step_13_scalar_perturbations` | Scalar-led, axial, polar QNM channels; isospectrality breaking |
| `step_14_kerr_tep` | Kerr-TEP shadow, ISCO spin scan, frame-dragging |
| `step_15_dynamical_signatures` | PPN parameters, dipole radiation, c_T, combined eta constraints |

### Derive — Geometry & Invariants (16–19)

| Step ID | Description |
|---------|-------------|
| `step_16_exact_geometry` | Exact curvature invariants, determinant, inverse metric |
| `step_17_null_expansions` | Null expansions for TEP matter metric |
| `step_18_observer_redshift` | Gravitational redshift for static observers |
| `step_19_observer_frequency_transfer` | Invariant emitter-receiver frequency transfer |

### Derive — Corrected Observables & Fundamentals (20–24)

| Step ID | Description |
|---------|-------------|
| `step_20_corrected_observables` | Corrected exterior observables (fixed-ADM normalization) |
| `step_21_cT_derivation` | Tensor speed c_T from principal symbol |
| `step_22_quadratic_action` | Exact quadratic action for perturbation channels |
| `step_23_characteristic_matrices` | Tensor characteristic matrices |
| `step_24_conformal_invariance_check` | Conformal invariance for null geodesics |

### Derive — QNM Solvers (25–30)

| Step ID | Description |
|---------|-------------|
| `step_25_qnm_solver` | Perturbative sGB QNM solver |
| `step_26_coupled_spectral_solver` | Coupled spectral solver for mixed modes |
| `step_27_qnm_horizon_branch` | Horizon-bearing sGB QNM |
| `step_28_qnm_validation` | QNM validation vs published results |
| `step_29_matrix_leaver` | Matrix continued-fraction QNM solver |
| `step_30_taylor_recurrence` | Taylor-series recurrence for background potentials |

### Derive — Interior (31–34)

| Step ID | Description |
|---------|-------------|
| `step_31_frobenius_analysis` | Frobenius regularity at temporal minimum |
| `step_32_interior_integration` | sGB scalar integration on Hayward background |
| `step_33_interior_analysis` | Interior regularity, Lorentzian, lapse checks |
| `step_34_solve_interior` | TEP interior solver, matter metric construction |

### Derive — Phantom Mass & EHT (35–39)

| Step ID | Description |
|---------|-------------|
| `step_35_mass_bias_sign` | Mass-bias sign equation |
| `step_36_phantom_mass_critical` | Phantom Mass critical g analysis, EHT comparison |
| `step_37_phantom_mass_raytrace` | Phantom Mass null geodesic ray-tracing |
| `step_38_phantom_mass_scan` | Phantom Mass parameter scan |
| `step_39_eht_visibility_fit` | EHT visibility-domain joint inference |

### Derive — Kerr & GW (40–41)

| Step ID | Description |
|---------|-------------|
| `step_40_kerr_sgb_true` | True rotating sGB shadow and QNM |
| `step_41_gw250114_confrontation` | GW250114 corrected confrontation |

### Figures (42)

| Step ID | Description |
|---------|-------------|
| `step_42_generate_figures` | 10 figures from paper plan |

### S-star Inference (43–49)

| Step ID | Description |
|---------|-------------|
| `step_43_sstar_data` | S-star data acquisition (145 astrometric + 44 RV epochs) |
| `step_44_gr_fit` | Conventional GR fit |
| `step_45_mass_bias` | Mass-bias sign (Gate -1) |
| `step_46_tep_transfer` | TEP transfer-function fit |
| `step_47_joint_forward` | Joint forward-model, composition vs calibration |
| `step_48_likelihood` | Likelihood comparison (Bayes, BIC, AIC, Wilks) |
| `step_49_posterior` | Phantom Mass posterior via MCMC |

## Architecture

```
scripts/
  run_pipeline.py               # Single entry point — runs all 50 steps
  steps/                        # All computational scripts (step_00–step_49)
    bh_common.py                #   Core model, solver, shared utilities
    step_00_*.py ... step_49_*.py  #   50 numbered pipeline steps
  utils/                        # Utility scripts
    logger.py                   #   TEPLogger (color-coded, file logging)
    generate_site_pdf.py        #   PDF generation from site
    convert_equations.py        #   Equation formatting utility
    test_curvature_regression.py #  Curvature regression tests

data/
  raw/                          # Downloaded raw data (EHT CSV files, etc.)
    eht_m87/                    #   EHT M87* 2017 calibrated visibilities
  processed/                    # Compiled measurement tables
    eht_measurements.json       #   EHT M87*/Sgr A* published measurements
    ligo_qnm_measurements.json  #   LIGO GWTC QNM measurements
```

## Data Sources

All observational data is downloaded from public repositories in step_00:

| Source | URL | Citation |
|--------|-----|----------|
| EHT M87* 2017 calibrated data | `github.com/eventhorizontelescope/2019-D01-01` | EHT Collaboration 2019, ApJL 875, L1. DOI: 10.25739/g85n-f134 |
| EHT M87* 2018 ring parameters | `aanda.org/articles/aa/full_html/2024/01/aa47932-23` | EHT Collaboration 2024, A&A 681, A79 |
| EHT M87* shadow measurement | Published paper | EHT Collaboration 2019, ApJL 875, L1. DOI: 10.3847/2041-8213/ab0ec7 |
| EHT Sgr A* shadow measurement | Published paper | EHT Collaboration 2022, ApJL 930, L12. DOI: 10.3847/2041-8213/ac6674 |
| LIGO GW150914 ringdown | Published paper | LIGO/Virgo Collaboration 2016, PRL 116, 221101. DOI: 10.1103/PhysRevLett.116.221101 |

All downloads are verified with SHA-256 checksums and minimum file size checks.

## Conventions

- **Step IDs**: `step_NN_description` (e.g., `step_01_field_equations`)
- **Output files**: `results/<step_id>.json`, `results/<step_id>.csv`
- **Log files**: `logs/<step_id>.log` (verbose, no ANSI colors)
- **Figures**: `results/figures/figN_name.png` and `.pdf`
- **Pipeline log**: `logs/tep_bh_pipeline.log`
- **Pipeline results**: `results/pipeline_results.json`
- **Logging**: All steps use `print_status(message, level)` from `scripts/utils/logger.py`
  - Levels: `INFO`, `PROCESS`, `SUCCESS`, `WARNING`, `ERROR`, `TEST`, `TITLE`
- **Shared utilities**: `scripts/steps/bh_common.py` (paths, JSON/CSV writers, logger factory)
- **Utility scripts**: `scripts/utils/` (logger, PDF, equation conversion, tests)
- **Old derive versions**: archived in `archive/scripts_derive/`

## Key Results

### Step 00 — Data Download

- 6 data sources downloaded from public repositories
- EHT M87* 2017 calibrated visibilities (4 observation days, ~2 MB total)
- Compiled measurement tables: EHT M87*/Sgr A* shadows, LIGO GW150914 QNMs
- All files verified with SHA-256 checksums

### Step 01 — Field Equation Solution

- **Kretschmann scalar finite** at all radii (vs Schwarzschild divergence)
- At r = 0.1M: K_TEP = 5.8 vs K_Schw = 4.8×10⁷ (8 orders of magnitude reduction)
- **Temporal horizon** at r_t ≈ 0.57M (28.5% of Schwarzschild radius)
- **Physical distance diverges** as r → 0 (r = 0 at infinite physical distance)
- **Exterior (r > 2M) is exactly Schwarzschild** (GR recovered)

### Step 02 — Perturbation Analysis

- Regge-Wheeler and Zerilli potentials computed for TEP physical metric
- **QNM frequencies** (l=2, n=0): TEP ω = 0.409 - 0.095i vs Schw ω = 0.389 - 0.087i
- **5.2% shift** in QNM frequency (testable with LIGO ringdown)
- **No echo cavity** near temporal horizon → no gravitational wave echoes

### Step 03 — Ray-Tracing

- Photon sphere at r = 3M, ISCO at r = 6M (exterior = Schwarzschild)
- Shadow predictions: M87* = 39.69 μas, Sgr A* = 53.26 μas
- Consistent with EHT measurements within 1σ

### Step 04 — Accretion Dynamics

- ISCO at r = 6M with radiative efficiency η = 5.72%
- Epicyclic frequencies: radial vanishes at ISCO, vertical = orbital
- **3:2 QPO resonance** found at r ≈ 10.8M (ratio = 0.667)
- Redshift factor at ISCO: g = √(2/3) ≈ 0.816

### Step 05 — Observational Constraints

- **EHT M87* shadow**: 39.69 μas (TEP) vs 42 ± 3 μas (measured) → 0.77σ
- **EHT Sgr A* shadow**: 53.26 μas (TEP) vs 48.7 ± 7.0 μas (measured) → 0.65σ
- **LIGO GW150914 f₂₂₀**: 212 Hz (TEP) vs 251 ± 5 Hz (measured) → TEP closer than Schw (201 Hz)
- **TEP QNM shift**: 10.6 Hz (5.2%) — testable with next-generation detectors
- All predictions trace to real downloaded data with full uncertainty propagation

### Step 06 — GW190521 Mass Gap

- **LIGO posterior**: 47,042 samples from Zenodo 4057131 (Isi et al. 2020)
- **Standard interpretation**: m1 = 85.4 M_sun, 99.3% probability in mass gap (65-130 M_sun)
- **TEP assessment**: A ~ 1 in the exterior where binary mergers occur, so GW masses are unchanged
- **Conclusion**: The mass gap problem is not addressed by TEP temporal redshift

### Step 07 — QPO Frequency Lock

- **4 microquasars** analyzed: GRS 1915+105, XTE J1550-564, GRO J1655-40, H 1743-322
- **3:2 ratio pairs**: majority of QPO pairs show ratio ~ 1.5 (within 0.15)
- **Epicyclic resonance**: 3:2 ratio at r = 10.8M from Omega_r/Omega_theta = 2/3 (standard GR)
- **Conclusion**: This is a GR prediction, not unique to TEP (conformal factor cancels for orbital dynamics)

### Step 08 — JWST Early SMBHs

- **6 high-z BHs** compiled: z = 7-10.6, masses 4x10^7 to 1.7x10^9 M_sun
- **Eddington-limited growth** from 100 M_sun seed: falls short by factors of 100-10000x
- **TEP assessment**: A ~ 1 in the exterior, so observed masses are true masses
- **Conclusion**: The JWST early SMBH problem is not addressed by TEP

### Step 09 — Spin Bias

- **12 published spin measurements** compiled from X-ray continuum fitting and iron line
- **Observed distribution**: skewed toward high spin, inconsistent with uniform [0,1] (KS test)
- **TEP assessment**: A ~ 1 and B ~ 0 in the exterior, so spin measurements are unaffected
- **Conclusion**: The observed spin clustering is not explained by TEP

### Step 10 — TDE Missing Flares

- **9 TDEs** compiled with observed peak luminosities and BH masses
- **Observation**: All 9 TDEs are significantly sub-Eddington
- **TEP assessment**: A ~ 1 in the exterior where emission escapes, so TDE luminosities are unaffected
- **Conclusion**: The sub-Eddington TDE luminosities are not explained by TEP

### Step 11 — EHT Polarization

- **M87* and Sgr A*** polarimetric data from EHT 2021 and 2024
- **Observed polarization**: M87* ~4%, Sgr A* ~7% (vs. ~10-15% predicted by GR+MHD)
- **TEP assessment**: A ~ 1 and B ~ 0 in the exterior, so photon polarization is unaffected
- **Conclusion**: The low polarization fraction is not explained by TEP

## Dependencies

```bash
pip install numpy scipy matplotlib sympy
```

## Outputs

```
results/
  step_00_data_download.json          # Download manifest with checksums
  step_01_field_equations.json        # Summary diagnostics
  step_01_field_equations.csv         # Radial profiles (20 columns)
  step_02_perturbations.json          # QNM frequencies, echo check
  step_02_perturbations.csv           # Potential profiles
  step_03_raytracing.json             # Photon sphere, ISCO, EHT comparison
  step_03_raytracing.csv              # Key observables
  step_04_accretion.json              # ISCO, epicyclic, QPO analysis
  step_04_accretion.csv               # Orbital profiles
  step_05_observational_constraints.json  # Chi-squared, p-values
  step_05_observational_constraints.csv   # Constraint comparison table
  step_06_gw190521_mass_gap.json     # GW190521 posterior analysis
  step_06_gw190521_mass_gap.csv      # Mass gap probability vs A(phi)
  step_07_qpo_frequency_lock.json    # QPO 3:2 resonance analysis
  step_07_qpo_frequency_lock.csv     # QPO measurements and ratios
  step_08_jwst_early_smbhs.json      # JWST high-z BH growth analysis
  step_08_jwst_early_smbhs.csv       # Eddington vs TEP comparison
  step_09_spin_bias.json             # Spin bias correction analysis
  step_09_spin_bias.csv              # Spin measurements and corrections
  step_10_tde_missing_flares.json    # TDE luminosity suppression analysis
  step_10_tde_missing_flares.csv     # TDE luminosities and A(phi)
  step_11_eht_polarization.json      # EHT birefringence analysis
  step_11_eht_polarization.csv       # Polarization fit results
  pipeline_results.json               # Overall pipeline status
  figures/                            # 10 figures (PNG + PDF)

data/
  raw/                                # Downloaded raw data
    eht_m87/                          # EHT M87* CSV files
  processed/
    eht_measurements.json             # Compiled EHT measurements
    ligo_qnm_measurements.json        # Compiled LIGO QNM measurements

logs/
  tep_bh_pipeline.log                 # Pipeline runner log
  step_00_data_download.log           # Step 00 verbose log
  step_01_field_equations.log         # Step 01 verbose log
  step_02_perturbations.log           # Step 02 verbose log
  step_03_raytracing.log              # Step 03 verbose log
  step_04_accretion.log               # Step 04 verbose log
  step_05_observational_constraints.log  # Step 05 verbose log
  step_42_generate_figures.log                # Figure generation log
```

## Reproducibility

All results are fully reproducible:

1. **Data provenance**: Every measurement traces to a published paper with DOI
2. **Checksums**: All downloaded files verified with SHA-256
3. **Deterministic solver**: Same parameters → same results (no random seeds)
4. **Versioned parameters**: Model parameters documented in every JSON output
5. **Full logging**: Every step produces a verbose log file
6. **No fabricated data**: All numbers in the manuscript trace to `results/` files
