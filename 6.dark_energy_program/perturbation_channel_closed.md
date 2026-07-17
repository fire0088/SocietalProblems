# The Perturbation Channel Is Closed for Horizon-Entropy Dark Energy

*Consolidated derivation — why linear structure growth (fσ₈) adds no
information beyond the background expansion H(z) for Barrow, Tsallis, and
Rényi holographic dark energy.*

---

## 1. Motivation

Recent DESI-DR2 analyses of horizon-entropy dark energy explicitly call for
confronting Barrow/Tsallis HDE "at the perturbative, structure-growth level,
using CMB temperature and polarization and weak lensing," since the existing
constraints are purely at the background (BAO + distance-prior) level. This
note answers that call with a negative result: **for the event-horizon
holographic dark-energy family, the linear growth rate carries no
degeneracy-breaking information beyond H(z).** We derive the mechanism, so the
result is a statement about the physics of the models rather than a numerical
coincidence.

The result reinforces, rather than competes with, the background-level model
ordering (steeper ρ(L) helps, flatter hurts, none beats ΛCDM): the growth data
are degenerate with the background already used to establish that ordering.

---

## 2. Growth depends on cosmology only through H(a)

For smooth (non-clustering) dark energy, the linear matter growth factor
D(a) obeys, in e-folds N = ln a,

    D'' + (2 + H'/H) D' - (3/2) Ω_m(a) D = 0.

Every cosmology-dependent quantity here is a function of H(a) alone:
the friction term 2 + H'/H is the background expansion, and the source
(3/2)Ω_m(a) = (3/2) Ω_{m0} a^{-3} / (H/H_0)^2 is again H(a) plus the matter
density. There is **no independent dark-energy degree of freedom** in the
growth equation. Two models with identical H(a) therefore produce identical
D(a), hence identical fσ₈(z).

**Numerical confirmation.** A CPL model (w₀,wₐ) = (−0.95, −0.20) tuned to
mimic ΛCDM differs from it by at most 0.53% in H(z) over 0 < z < 2. The
resulting fσ₈ differs by at most 0.50% — the same order as the background
difference, and an order of magnitude below current growth-rate errors
(~5–10% per point). The growth observable inherits the background
(non-)difference and adds nothing.

The only texture is a sign flip in Δfσ₈ near z ≈ 0.6: because D(z) is a
cumulative integral of the history while H(z) is local, fσ₈ carries a
different *weighting* of the same H(a) information — not new information.

---

## 3. Independent information requires clustering dark energy

fσ₈ becomes independently informative only if dark energy itself clusters,
i.e. has a rest-frame sound speed c_s² < 1, adding a source term not reducible
to H(a). Parametrizing a clustering fraction μ (μ = 0 smooth, standard),
a fixed-H(z) test gives:

    μ = 0.05  →  ~5% fσ₈ shift   (edge of detectability)
    μ = 0.10  →  ~10% fσ₈ shift  (detectable)

So the entire question reduces to the dark-energy sound speed. Everything
below establishes that, for the event-horizon holographic family, the
effective sound speed drives sub-horizon perturbations to zero.

---

## 4. The sound speed is not free — and is not adiabatic

The adiabatic sound speed c_a² = w − ẇ/[3H(1+w)] is fixed by the background,
but it is negative or singular for w < −1/3 (generic for dark energy, and the
1+w denominator diverges near the phantom line several variants approach). A
negative c_s² is a catastrophic small-scale instability. **The adiabatic value
must therefore be discarded** — confirmed numerically: imposing it makes
sub-horizon δ_DE blow up as ~10⁷, 10¹⁰⁵, 10³⁰³ with increasing k.

Holographic dark energy has no Lagrangian kinetic term from which to read a
rest-frame c_s² — ρ_DE is defined non-locally by a horizon. The sound speed is
thus an *imposed* prescription, and the honest route is to compute what the
horizon prescription itself implies.

---

## 5. Derivation: horizon perturbations phase-average away

**Density from horizon.** With ρ_DE ∝ R_h^n (n = −2 standard; Barrow
n = Δ−2, Tsallis n = 2δ_T−4), exactly

    δ_DE = n · (δR_h / R_h).

**Horizon perturbation is a future light-ray integral.** The future event
horizon R_h = a ∫_t^∞ dt'/a is a null construct. Perturbing the radial null
geodesic in Newtonian gauge (ds² = 0 ⇒ dx = dt/a (1+2Φ)), the comoving horizon
picks up

    δχ_k(t) = 2 e^{ikx} ∫_0^{χ_h} Φ_k(u) e^{iku} du,   u = χ(t,t').

**Scale dependence from the oscillatory kernel.** For a sub-horizon mode
(k ≫ aH) the kernel e^{iku} oscillates across the integration range.
Integration by parts with any smooth envelope g(u) = Φ_k(u):

    ∫_0^{χ_h} g e^{iku} du = [g e^{iku}/(ik)]_0^{χ_h} − (1/ik)∫ g' e^{iku} du
                           = O(1/k).

Hence |δ_DE| ~ (aH/k)|Φ| for k ≫ aH. **The mechanism is rapid phase-averaging
of the potential along the future light-ray**: a small-scale potential
fluctuation contributes essentially nothing to the non-local horizon integral.

**Relative to matter.** Sub-horizon, Φ is sourced by matter via Poisson,
k²Φ = −4πG a² ρ_m δ_m, so Φ ~ δ_m/k². Therefore

    δ_DE / δ_m ~ 1/k³   (sub-horizon).

Dark energy is smooth to third order in aH/k relative to the clustering
matter. (Numerically, k·|δL| plateaus flat across k = 1–1000, confirming the
1/k kernel suppression; the Poisson factor supplies the remaining 1/k².)

**Consistency checks.** A three-closure integration across k = 1–1000 gives:
imposed c_s²=1 → δ_DE stays ~10⁻⁹ (smooth); adiabatic c_s²<0 → blow-up
(unphysical); horizon closure → tracks the c_s²=1 smooth branch. Note the
derived suppression is 1/k (relative to Φ), *weaker* than a true c_s²=1 fluid's
1/k²; we state the 1/k (and 1/k³ vs matter) scaling directly rather than
claiming an equivalent c_s²=1.

---

## 6. Extension to Rényi

Rényi HDE, ρ_DE = 3C²/[8πL²(1 + δπL²)], does not have L drop out of the
background (a coupled 2-D system). But its density is still ρ_DE = f(L), so

    δ_DE = n_eff(L) · (δL/L),   n_eff(L) = L f'(L)/f(L)
         = −2 − 2δπL²/(1 + δπL²),

running monotonically from −2 (small L, standard HDE) to −4 (large L,
Rényi-dominated) and **bounded in [−4, −2]**.

Crucially, Rényi modifies the *density law* f(L), not the *definition* of L:
the IR cutoff remains the future event horizon, a geometric null-geodesic
object. The oscillatory kernel e^{iku} is kinematic (from null propagation),
untouched by the coupled background, which changes only bounded, smooth,
non-oscillatory quantities (horizon size χ_h, envelope Φ(u), and the prefactor
n_eff). Integration by parts still yields O(1/k). Therefore

    |δ_DE| = |n_eff| · O(aH/k) → 0   sub-horizon,

and Rényi dark energy is smooth on the same footing as the power-law family.

---

## 7. The one condition the result depends on

The derivation assumes the IR cutoff is the **future event horizon**. If a
variant instead uses the **Hubble horizon** (L = 1/H) or a mixed cutoff, δL is
no longer a light-ray integral but the local quantity δL = −δH/H², with no
oscillatory kernel and no 1/k suppression — and dark energy could genuinely
cluster. The precise statement is therefore:

> **Horizon-entropy dark energy is smooth on sub-horizon scales — and fσ₈
> carries no information beyond H(z) — if and only if the model uses an
> event-horizon (light-ray) IR cutoff.**

Our models use the future event horizon, so the result holds. (Hubble-cutoff
HDE independently fails to produce late-time acceleration in the presence of
dark matter, so the event horizon is the standard choice.) A Hubble- or
mixed-cutoff construction is the one place in this analysis where the
perturbation channel could reopen — a well-defined open problem.

---

## 8. Conclusion

For event-horizon holographic dark energy — standard, Barrow, Tsallis, and
Rényi alike — sub-horizon dark-energy perturbations are suppressed as 1/k³
relative to matter by phase-averaging of the potential along the future event
horizon, with the entropy deformation entering only as a bounded O(1)
amplitude prefactor. Consequently the linear growth rate fσ₈ is degenerate
with the background H(z), and the perturbation-level confrontation cannot
separate these models from one another, or from ΛCDM, beyond what background
BAO + CMB distances already achieve. The proposed structure-growth test adds
no discriminating power for this model class.

### Honest limitations
- The derived scaling is 1/k (vs Φ) / 1/k³ (vs matter), not an effective
  c_s² = 1; stated as such.
- Only the leading light-ray phase term is kept; local δa/a metric
  corrections at the observation point (also O(aH/k)) are subleading and
  neglected — a fully rigorous version tracks them.
- Result is conditional on an event-horizon cutoff (Section 7).
- Separate open issue (background, not perturbations): the Rényi
  distance-prior CMB likelihood remains numerically unstable to z* ≈ 1090;
  only Rényi BAO-only is currently trustworthy.
