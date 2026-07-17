# Detectability of Sudden Dark-Energy Transitions — Code & Data

Code, data, and manuscript for the paper on whether expansion-history probes
(cosmic chronometers, Type Ia supernovae, anisotropic BAO) can detect and
characterize a sudden step in the dark-energy density, plus the horizon-entropy
(holographic) dark-energy extension.

All detection rates are computed end-to-end from the **real** covariances of each
probe. Every pipeline script is self-validating: it checks its null distribution
against the analytic chi^2_1 expectation before reporting rates.

## Layout

```
realdata/     Real observational data
  desi_dr2.py         DESI DR2 BAO (6 anisotropic tracers, DM/rd + DH/rd, correlated)
  PantheonSH0ES.dat   Pantheon+ SH0ES supernova compilation (1701 SNe)
  PPcov.cov           Full 1701x1701 Pantheon+ stat+sys covariance

pipeline/     Core detection-hierarchy pipeline (main result, Table 1)
  cc_data.py          Cosmic-chronometer compilation + Moresco covariance (validated
                      against Moresco 2023 CC-only H0=66.7+/-5.3)
  REGEN_lib.py        Shared library: step model, grid LRT, null validation
  REGEN_lib2.py       Extended library (signed amplitudes)
  REGEN_cc.py         CC detection rates
  REGEN_sne.py        Pantheon+ detection rates
  REGEN_bao.py        DESI DR2 BAO detection rates
  REGEN_fig.py        Detectability-hierarchy figure

analyses/     Supporting analyses (one script per paper result)
  mechanism_test.py   DM-only / DH-only / both — proves the 2-observable mechanism
  combine2.py         SNe+BAO vs BAO-alone (SNe add nothing for step detection)
  exclusion_test.py   What step amplitudes real DESI DR2 already excludes (A in [-0.22,0.23])
  exclusion_verify.py   Coverage check for the exclusion interval
  exclusion_fig.py    Exclusion figure
  location_scan.py    Detection rate vs transition redshift z_t
  location_fig.py     Location-scan figure
  width_test.py       Can BAO recover transition WIDTH? (detection yes, discrimination no)
  sign_cc/sne/bao.py  Signed-amplitude (thinning vs thickening) rates
  sne_residual.py     SNe residual figure — visual proof Om absorbs the step
  observability3.py   Unifying Fisher subspace-overlap criterion (orthogonal fraction)

models/       Horizon-entropy (holographic) dark-energy models
  barrow_hde3.py        Barrow HDE (validated: Delta=0 -> w0=-0.891)
  barrow_exclusion2.py  Barrow Delta vs DESI DR2 BAO
  barrow_cmb.py         + CMB shift parameter R (Delta<0.05)
  barrow_cmb_full2.py   Full 3-param CMB distance priors {R, l_A, omega_b} (Delta<0.028)
  multimodel.py         Unified power-law framework (Barrow & Tsallis share one ODE)
  multimodel_fit.py     Multi-model fit vs DESI BAO + CMB
  multimodel_final.py   Final comparison incl. Renyi + H0 scan
  renyi_coupled.py      Renyi HDE (coupled 2D system; validated low-z, CMB unreliable)
  renyi_fit.py          Renyi fit driver

data/         Result arrays (.npz) produced by the scripts above
paper/        Manuscript
  step_detectability.tex / .pdf   The paper (15 pp)
  references.bib                   Bibliography
  figures/                         All figures
```

## Reproducing

Requires: python3, numpy, scipy, matplotlib.

```
cd pipeline
python3 REGEN_cc.py      # CC detection rates
python3 REGEN_sne.py     # Pantheon+ rates
python3 REGEN_bao.py     # DESI DR2 rates
```

Each prints its null-median check (should match chi^2_1 median ~0.455) before
the rates. The models/ and analyses/ scripts run standalone from their own
directory and import realdata/ via a relative path (run from repo root, or
adjust sys.path).

## Key results

- Detection hierarchy (3sigma, 37.5% shed): CC ~0%, Pantheon+ ~1%, DESI DR2 BAO ~83%.
- Mechanism: BAO's DM/rd + DH/rd pair breaks the step-vs-Omega_m degeneracy that
  blinds single-distance probes. DM-only collapses to SNe-level blindness.
- Exclusion: real DESI DR2 already excludes a z=0.7 step shedding >19% (95% CL).
- Width: BAO detects transitions and prefers sharp ones, but cannot measure the
  width (sharp vs gradual near-degenerate at current precision).
- Barrow HDE: BAO-only Delta<0.88; +CMB Delta<0.028 (consistent with literature).
- Multi-model: data discriminate between entropy variants — steeper rho(L)
  (Tsallis) helps, flatter (Barrow) hurts, none beats LCDM.
- Hubble tension: a z=0.7 late-time step is empirically excluded as a resolution
  (shifts inferred H0 by ~2% of the gap).

## Caveats (see paper Sec. 12)

Reference details in references.bib were assembled from arXiv IDs and should be
verified against the published records before submission. The multi-model CMB
uses distance priors (not the full Planck likelihood) and fixes omega_b; the
Renyi CMB extension needs a more robust integrator (only its BAO-only fit is
validated).
