# Dark Energy Observability Program — Four Linked Papers

Eric Schultz (ORCID 0009-0006-6283-1696)

One question at four scales: **what can cosmology actually learn about dark energy?**
Each paper feeds the next; together they form a coherent program running from the
detectability of a single background feature up to an epistemic bound on whether
the cosmological-constant hypothesis can ever be confirmed.

---

## The four papers

### Paper A — Step paper (detectability)
*"Detectability and exclusion of localized dark-energy features in the expansion history."*
**Keeper result:** for a localized *step* in the dark-energy density, anisotropic
BAO vastly outperforms supernovae, because BAO's two observables (D_M/r_d and
D_H/r_d) break the step-vs-Ω_m degeneracy a single SN magnitude cannot. Verified
on real Pantheon+ and DESI DR2, cross-checked analytically.
**Status:** essentially done; needs two trivial bibliography fixes (Landim→2017,
de-Souza-not-Abdalla). *(Not included in this zip — lives in the step-paper
working directory.)*

### Paper B — Companion (`hde_companion.pdf`)
*"Which generalized entropy? Discrimination among horizon-entropy dark-energy
models with DESI DR2 and Planck, and the closure of the structure-growth channel."*
**Two results:** (1) a uniform discrimination of Barrow/Tsallis/Rényi/Kaniadakis
horizon-entropy models (steeper-helps / flatter-hurts / none beats ΛCDM); (2) the
distinctive result — the **closed structure-growth channel**: for event-horizon
holographic dark energy, δ_DE/δ_m ~ 1/k³ via light-cone phase-averaging, so fσ8
carries no information beyond the background H(z).
**Status:** complete in shape (6pp, compiles clean), bug-checked. One mechanical
gate remains before submission — see "Pending" below.

### Paper C — Null-horizon classification (`null_horizon.pdf`)
*"Which dark energy can structure growth see? A response-kernel classification of
geometric dark energy."*
**Class-level theorem:** the growth-observability of any geometric dark energy is
set by the UV scaling of the kernel through which its density responds to the
gravitational potential. Proven as a **locality dichotomy via Riemann–Lebesgue**:
local kernels cluster (1/k²), genuinely nonlocal (null-line) kernels are
growth-empty (1/k³). The event horizon is the unique cutoff with an
*unconditionally* closed channel. The particle-horizon case closes both
numerically (high-precision 1/k³) and structurally (resonance would require the
null ray to travel at the acoustic sound speed — impossible, since light is
always faster than sound).
**Status:** complete, referee-grade draft (5pp, compiles clean).

### Paper D — Capstone (`observability.pdf`)
*"Can cosmology confirm that dark energy is a cosmological constant?
An observable-subspace bound."*
**Unifying object:** observability of a DE feature = ‖Π_probe R‖, the response
projected onto the probe's sensitivity subspace. Papers A–C are one construction
under three projections. **New theorem:** finitely many probes span a
finite-dimensional observable subspace, so most of the infinite-dimensional DE
response is invisible in principle.
**The vacuum-energy bound (the headline):** a cosmological constant *is*
falsifiable in the directions the data cover — clustering (obs 0.87) and low-z
background evolution (obs 0.92, where DESI's evolving-w hint lives). But a
structured family of near-Λ deviations is blind — high-z/early DE (0.21),
super-horizon (0.45), the mid-z desert (0.22). So the strongest honest statement
is **"dark energy is a cosmological constant *within the observable subspace*,"**
never without the qualifier. Confirming DE = vacuum energy is not about measuring
w to more decimals; it needs a probe with support in the blind zones.
**Status:** 3pp draft, compiles clean.

---

## How they connect

- A establishes the observability functional for **one background feature**.
- B applies it to a **model family's growth signatures** and proves the closed channel.
- C generalizes the closed channel to a **whole class of geometric models**.
- D abstracts all three into an **epistemic bound**, with the vacuum-energy
  hypothesis as the sharp stake.

Internal cross-consistency verified: Paper B ("event-horizon HDE growth channel
closed") and Paper D ("distinguishing signal under-observed") tell the same story
from two angles — an event-horizon model is the smooth branch, so it has no
distinguishing signal for D's probes to see.

---

## Pending items (all minor / mechanical unless noted)

**Paper B:**
1. Run the om_b marginalization on the real pipeline: wire the two hooks in
   `omb_marginalize_RECIPE.py` into `barrow_cmb_full2.py`, replace the fixed-om_b
   χ² with `profiled_chi2` per model, regenerate the discrimination table with
   profiled om_b (turns the red provisional numbers black).
2. Confirm the steeper-helps/flatter-hurts ordering survives the non-uniform
   per-model om_b penalties. *(Bug-check showed the ordering **could** shift, so
   this must be checked, not assumed.)*
3. Pin the Rényi z_match ≥ 50 with a convergence check (z_match=15 carries a
   ~0.2 l_A ≈ 2σ systematic).
4. Swap the article preamble → revtex for journal submission.
5. Confirm Luciano 2026 (arXiv:2603.21218) exists and the exact author lists on
   both Luciano papers (post-dates the knowledge cutoff used when drafting).

**Papers C and D:**
6. Swap article → revtex preamble for submission.
7. **The one substantive upgrade (shared by C's optional check and D):** replace
   the schematic (physically-motivated) probe kernels with real survey Fisher
   kernels (DESI fσ8, Euclid/LSST lensing, Planck/SO/CMB-S4). The *structural*
   results do not depend on this; only the quantitative fractions (0.74, 0.47,
   0.21, etc.) would shift.

**Submission order:** A ships first (worked example the others cite) → B → C
(cites B as its worked event-horizon example) → D (cites all three).

---

## Bug-check record (this program)

Adversarial checks found and fixed:
- **om_b module bounds bug** (Paper B): default bounds (0.0205, 0.0240) clipped
  models with extreme best-fit om_b, inflating χ² by ~30000 and potentially
  corrupting the ordering. Fixed: widened to (0.016, 0.028) + bound-hit guard.
- **Growth-index derivation bug** (Paper B §4.6): the γ≈4/7 conclusion was right
  but the reasoning bridging it via the clustered formula was wrong (that formula
  gives 0 at c_eff²=1). Rewrote as two separate branches; citation corrected to
  Mehrabi et al. PRD 92 123513 (2015).
- **Rényi z_match systematic** (Paper B): characterized, not a blow-up; fix is
  z_match ≥ 50.
- **r_s calibration** (Paper B): confirmed the "144 vs 147 Mpc" was a false alarm
  (r⋆ at photon decoupling = 144.7 Mpc for l_A vs r_drag = 147.1 Mpc for BAO),
  verified against Planck.
- **Cross-paper consistency** (B ↔ D): verified consistent.
- **Cross-citation title spelling**: fixed localized/localised drift.

All three included papers compile clean with zero errors and zero undefined
references.

---

## File manifest

- `hde_companion.pdf` / `.tex` — Paper B
- `null_horizon.pdf` / `.tex` — Paper C
- `observability.pdf` / `.tex` — Paper D
- `omb_marginalize_RECIPE.py` — drop-in om_b marginalization module for Paper B
  (self-tested, bounds-bug fixed)
- `references.bib` — shared bibliography (companion)
- `perturbation_channel_closed.md` — standing note on the closed-channel result
