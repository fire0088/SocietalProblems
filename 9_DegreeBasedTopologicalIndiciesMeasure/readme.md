# Degree-based topological indices measure molecular size: a null-baseline protocol for QSPR

**Eric Schultz** (ORCID 0009-0006-6283-1696), Independent Researcher

*Draft — prepared for submission to Journal of Cheminformatics*

---

## Abstract

In quantitative structure–property relationship (QSPR) studies using degree-based topological indices, the reported coefficient of determination and the descriptor's actual contribution move in **opposite** directions: as R² rises from 0.34 to 0.94, the improvement over simply counting heavy atoms collapses from +0.148 to +0.005. The most impressive published results are those in which the descriptor is doing the least work.

We explain this. Writing a degree-based index as `TI_f(G) = m·⟨p(G), f⟩`, where `m` is the edge count and `p(G)` the edge-degree-pair profile, we bound the correlation between any such index and `m` in terms of three quantities: the variation of the weighting function along directions in which real molecular profiles differ, the mean edge weight, and the coefficient of variation of molecular size in the sample. A corollary quantified over *all* properties follows: no degree-based index can differ from edge count in its correlation with **any** property by more than `√(2(1−ρ))`, where `ρ` is its correlation with edge count. On 8,583 drug-like molecules the bound is valid for all 14 indices tested and tight to within 0.013 for ten; the corollary holds without exception across 70 index–property pairs.

Empirically, 14 standard indices span one dimension (PC1 = 97.5% on 6,000 ZINC molecules) and each correlates with edge count at r ≥ 0.947. Against molecular weight — exactly additive, and among the most commonly used response variables in this literature — no index improves on heavy-atom count on any dataset tested (ΔR² between −0.030 and +0.014). Against measured melting points (N = 8,583) the best single index reaches R² = 0.157, below a two-number baseline of atom count plus ring count (0.218); the 14 indices jointly reach 0.336, and adding atom and ring counts to them changes nothing.

At the sample sizes typical of this literature (N ≈ 15) the reported statistic is barely informative: all 14 indices win at some frequency, and the expected best-of-14 R² exceeds the full-sample value by roughly 0.10. A permutation null places the 95th percentile at R² ≈ 0.35 for linear fits — and, because many studies fit quadratic or cubic models to fifteen molecules, at ≈ 0.50 and ≈ 0.63 respectively. Raising the fit order raises the bar a claim must clear rather than lowering it.

We recommend a referee-facing protocol whose central element is a single diagnostic: report ΔR², the improvement over a heavy-atom-count null. This reconciles the present literature with Gutman and Tošović's 2013 size-controlled benchmark, which reached the opposite conclusion on isomer sets and is widely cited by the studies it undercuts.

**Keywords:** topological index, QSPR, molecular descriptor, degree-based index, null baseline, chemical graph theory

---

## 1. Introduction

A molecular graph `G` represents a hydrogen-suppressed molecule with atoms as vertices and bonds as edges. A *degree-based topological index* assigns a real number to `G` by summing a symmetric function of the endpoint degrees over all edges:

```
TI_f(G) = Σ_{uv ∈ E(G)} f(d_u, d_v)
```

Instances include the Zagreb indices (Gutman & Trinajstić 1972), the Randić connectivity index (Randić 1975), the atom–bond connectivity (ABC) and geometric–arithmetic (GA) indices, the forgotten index (Furtula & Gutman 2015), and the Sombor index (Gutman 2021). Several thousand such invariants have now been proposed.

### 1.1 The template under examination

A recognisable publication template has become common in pharmacy, pharmaceutical-chemistry and applied-mathematics journals. Its steps are:

1. select a therapeutic class and 10–25 drug molecules from it;
2. compute a battery of 6–20 degree-based indices;
3. obtain 5–13 physicochemical properties, typically from an online compilation;
4. regress each index against each property;
5. report, per property, the index with the largest R², described as the best predictor.

Representative recent instances span *PLOS ONE*, *Scientific Reports*, *ACS Omega*, *Frontiers in Chemistry*, *BioNanoScience* and *Journal of Mathematical Chemistry*. The properties most often used are molecular weight, molar refractivity, molar volume, polarizability, surface tension, boiling point and flash point.

### 1.2 The difficulty

Four of the most frequently used properties — molecular weight, molar refractivity, molar volume, polarizability — are *additive*: to good approximation each is a sum of atomic or group contributions, and hence a near-linear function of molecular size. Degree-based indices are also, as we show, near-linear functions of size. A high R² between the two is therefore close to guaranteed and carries little information about the index.

This is not a new suspicion. Gutman and Tošović (2013) tested the correlation ability of 20 vertex-degree-based indices against heats of formation and normal boiling points of **octane isomers**, and found that for many of the indices the correlation ability is weak or nil, with the augmented Zagreb and ABC indices performing best. Their design fixes `n` and `m`, which removes size variation entirely — and, we note, forces the choice of non-additive properties, since molecular weight and molar refractivity are literally constant within a constitutional-isomer class.

That paper is cited in a large fraction of the studies described in §1.1, often as introductory boilerplate, without its implication being drawn.

Randić and co-workers had earlier addressed descriptor intercorrelation directly, arguing that apparent parallelism between descriptors need not be a concern unless they are strictly collinear, and that orthogonalization reveals genuine structural differences even between nearly collinear descriptors, so that intercorrelated descriptors retain interest. We accept this argument as stated. We observe only that it concerns *orthogonalized combinations* of descriptors. It provides no licence for the comparative claim — "index X is the best predictor of property Y" — which is what the template produces.

### 1.3 Contribution

1. **A theorem** (§2) bounding the correlation between any degree-based index and edge count, in terms of three interpretable quantities, and identifying exactly what must be broken for such an index to carry non-size information.
2. **Dimensionality measurement** (§4.1): 14 standard indices span one effective dimension across two unrelated chemical spaces.
3. **A null-baseline analysis** (§4.3–4.4): against additive properties no index improves on heavy-atom count; against measured non-additive properties the whole battery is weak.
4. **A mechanism for the published literature** (§4.5): the reported R² is a readout of sample size dispersion, the index's marginal contribution falls as R² rises, and at N ≈ 15 the reported value is largely sampling noise.
5. **A protocol** (§5.3): report ΔR² over a heavy-atom-count null.

---

## 2. Theory

### 2.1 Profile decomposition

Let `P = {(i,j) : 1 ≤ i ≤ j ≤ Δ}` be the set of unordered degree pairs available at maximum degree `Δ` (for hydrogen-suppressed organic molecules, `Δ = 4`, so `|P| = 10`). For a molecular graph `G` with `m` edges define the **edge-degree-pair profile**

```
p(G) ∈ Δ(P),    p_π(G) = #{uv ∈ E(G) : {d_u, d_v} = π} / m
```

a point in the probability simplex over `P`. Regarding `f` as a vector in `ℝ^P`, every degree-based index factorises exactly:

```
TI_f(G) = m(G) · ⟨p(G), f⟩                                    (1)
```

The index is the product of an **extensive** factor (edge count) and an **intensive** factor (mean edge weight). All structural information beyond size is carried by `⟨p(G), f⟩`.

### 2.2 The bound

Consider a dataset `G_1, …, G_N` with edge counts `m_i` and profiles `p_i`. Let

```
p̄ = (1/N) Σ p_i,   u_i = p_i − p̄,   ḡ = ⟨p̄, f⟩,   CV = sd(m)/mean(m)
```

and assume `ḡ > 0`. Write `ṽ` for the mean-centred version of a vector `v ∈ ℝ^N`.

> **Proposition.** Let `E_i = m_i⟨u_i, f⟩` and `τ_f = ‖Ẽ‖ / (ḡ‖m̃‖)`. If `τ_f ≤ 1`, then
>
> ```
> r²(TI_f, m) ≥ 1 − τ_f²                                      (2)
> ```
>
> Moreover, letting `Σ` be the `m²`-weighted covariance of the profiles, `Σ = Σ_i w_i u_i u_iᵀ` with `w_i = m_i² / Σ_j m_j²`, and `f_c = f − c**1**` for any scalar `c`,
>
> ```
> τ_f ≤ √(f_cᵀ Σ f_c) · rms(m) / (ḡ · sd(m))                   (3)
> ```

**Proof.** From (1), `TI_f(G_i) = ḡ m_i + E_i`. Set `S = ḡ m`, so `TI_f = S + E` and, since centring is linear, `T̃I_f = S̃ + Ẽ` with `S̃ = ḡ m̃` parallel to `m̃`.

For the correlation bound, note `r(TI_f, m) = cos ∠(T̃I_f, m̃) = cos ∠(T̃I_f, S̃)`. The set `{S̃ + e : ‖e‖ ≤ ρ}` is a Euclidean ball of radius `ρ` centred at `S̃`; the maximum angle it subtends at the origin from `S̃` is `arcsin(ρ/‖S̃‖)` when `ρ < ‖S̃‖`. Taking `ρ = ‖Ẽ‖` gives `∠(T̃I_f, S̃) ≤ arcsin(τ_f)` and hence `r ≥ √(1 − τ_f²)`, which is (2).

For (3): each `u_i` satisfies `⟨u_i, **1**⟩ = 0`, so `⟨u_i, f⟩ = ⟨u_i, f_c⟩` for every `c` — the weighting function is determined only up to an additive constant. Then `‖Ẽ‖ ≤ ‖E‖` (centring is an orthogonal projection) and

```
‖E‖² = Σ_i m_i² ⟨u_i, f_c⟩² = (Σ_j m_j²) · f_cᵀ Σ f_c
```

so `‖E‖ = √(Σ_j m_j²) · √(f_cᵀ Σ f_c) = √N · rms(m) · √(f_cᵀ Σ f_c)`. Dividing by `ḡ‖m̃‖ = ḡ√N · sd(m)` gives (3). ∎

### 2.3 A corollary over all properties

The Proposition constrains the relationship between `TI_f` and edge count. The claim of interest concerns properties. The gap closes in one step.

> **Lemma.** For any random variables `X, Y, Z`, `|r(X,Z) − r(Y,Z)| ≤ √(2(1 − r(X,Y)))`.
>
> *Proof.* Let `x, y, z` be the centred, unit-normalised data vectors. Then `r(X,Z) − r(Y,Z) = ⟨x − y, z⟩ ≤ ‖x − y‖‖z‖ = ‖x − y‖`, and `‖x − y‖² = 2 − 2⟨x,y⟩ = 2(1 − r(X,Y))`. Symmetry gives the absolute value. ∎

> **Corollary.** Let `ρ = r(TI_f, m) ≥ √(1 − τ_f²)` by (2). Then for **every** property `P` whatsoever,
>
> ```
> |r(TI_f, P) − r(m, P)| ≤ √(2(1 − ρ)) ≤ √(2(1 − √(1 − τ_f²)))        (4)
> ```

This is the operative statement. It is not a claim about the five properties we happened to test: it says that for any response variable at all — measured, computed, physicochemical, biological — a degree-based index and a bond count are interchangeable to within a margin fixed by `τ_f`, which is a property of the index and the molecule set alone.

For the atom–bond connectivity index on our melting-point set, `ρ = 0.9983` and the margin is 0.058. Across 14 indices and 5 properties (measured melting point, molecular weight, molar refractivity, topological polar surface area, cLogP) there are no violations, and the largest observed deviation for ABC is 0.015 — well inside the guarantee (§4.2).

### 2.3 Interpretation

The bound factorises into three parts, each with a chemical reading.

- **`√(f_cᵀ Σ f_c)` — index spread along realised variation.** Not how much `f` varies over all conceivable degree pairs, but how much it varies along the directions in which real molecular profiles actually differ. This is why the naive sandwich `f_min·m ≤ TI_f ≤ f_max·m`, which uses the worst-case range of `f`, is far too weak to be useful: it ignores the covariance structure entirely.
- **`1/ḡ` — normalisation.** Only the *relative* spread of `f` matters; rescaling `f` leaves the index's correlation structure unchanged, as it must.
- **`rms(m)/sd(m) ≈ √(1 + CV²)/CV` — size dispersion.** This diverges as `CV → 0`.

The three factors are the complete list of things one could break in order to construct a degree-based index that is not a size proxy. Two are properties of the index; the third is a property of the *study design*, and is under the analyst's control.

The `CV` term recovers Gutman and Tošović's result immediately. On an isomer set, `CV = 0` exactly, the bound is vacuous, and the indices are free to differ arbitrarily from edge count and from one another — which is precisely the regime in which a comparison among them is meaningful. Their study and the studies of §1.1 were never in conflict; they sit at opposite ends of a continuum that (2)–(3) parameterises.

---

## 3. Methods

### 3.1 Index battery

Fourteen degree-based indices with weighting functions `f(a,b)`: first Zagreb `a+b`; second Zagreb `ab`; forgotten `a²+b²`; hyper-Zagreb `(a+b)²`; Randić `(ab)^(−1/2)`; sum-connectivity `(a+b)^(−1/2)`; ABC `√((a+b−2)/(ab))`; geometric–arithmetic `2√(ab)/(a+b)`; harmonic `2/(a+b)`; Sombor `√(a²+b²)`; augmented Zagreb `(ab/(a+b−2))³`; inverse sum indeg `ab/(a+b)`; symmetric division deg `a/b + b/a`; redefined third Zagreb `ab(a+b)`.

### 3.2 Datasets

| set | N | source | properties |
|---|---|---|---|
| Drug set | 38 | curated from the therapeutic classes recurring in the literature of §1.1 | computed (RDKit) |
| ZINC drug-like | 6,000 | random subsample of a 250k drug-like collection | — (dimensionality only) |
| ESOL | 1,105 | Delaney aqueous solubility | measured log S |
| Melting point | 8,583 | Open Notebook Science melting-point compilation, joined to structures via PubChem CAS lookup | measured T_m |
| Boiling point | 6,400 | Yaws boiling-point compilation, same join | T_b |

Filters for the T_m/T_b sets: 10–60 heavy atoms, maximum degree ≤ 4, must contain carbon, salts and mixtures (`.` in SMILES) excluded.

### 3.3 Null baseline and statistics

The null descriptor is the **heavy-atom count** `n`. For each property we report `R²(n)`, `R²` of the best of the 14 indices, and `ΔR² = R²(best) − R²(n)`. Out-of-sample performance uses leave-one-out or repeated 70/30 splits. Small-sample behaviour is characterised by repeated random draws of 15 molecules, matching the modal published study size.

---

## 4. Results

### 4.1 The indices span one dimension

| | drug set (N = 38) | ZINC (N = 6,000) |
|---|---|---|
| PC1 variance explained | 98.64% | 97.49% |
| PC2 | 1.13% | 1.76% |
| PC3 | 0.21% | 0.69% |
| median pairwise \|r\| | 0.989 | 0.979 |
| minimum pairwise \|r\| | 0.934 | 0.908 |
| min r(TI, m) | 0.961 | 0.947 |

The result replicates across two unrelated chemical spaces, and is not an artefact of size range: the ZINC sample has roughly half the size dispersion of the drug set (CV(m) = 0.213 vs 0.408) yet shows the same collapse.

A rank statement makes the point exactly. On ESOL, cross-validated R² for measured log S:

```
n alone              0.3199
best single index    0.3348
all 14 indices       0.3803
n + all 14 indices   0.3803      ← identical to four decimal places
```

Adding heavy-atom count to the fourteen indices changes out-of-sample performance not at all: `n` already lies in their span.

The observed edge-degree-pair profile is highly concentrated. Across the drug set, four of the ten available pairs account for 91% of all edges: (2,3) 41.3%, (2,2) 17.7%, (3,3) 17.5%, (1,3) 14.2%.

### 4.2 The bound is valid and tight

Evaluated on the 8,583-molecule melting-point set (CV(m) = 0.417):

| index | τ (exact) | τ (spectral bound) | 1 − τ²_spec | actual r²(TI, m) | slack |
|---|---|---|---|---|---|
| ABC | 0.0605 | 0.0608 | 0.9963 | 0.99663 | +0.0003 |
| GA | 0.0613 | 0.0614 | 0.9962 | 0.99648 | +0.0003 |
| sum-conn | 0.0924 | 0.0929 | 0.9914 | 0.99231 | +0.0010 |
| SDD | 0.2680 | 0.2686 | 0.9279 | 0.92885 | +0.0010 |
| harmonic | 0.1817 | 0.1828 | 0.9666 | 0.96868 | +0.0021 |
| Randić | 0.1567 | 0.1585 | 0.9749 | 0.97974 | +0.0049 |
| M1 | 0.1975 | 0.1986 | 0.9606 | 0.97271 | +0.0121 |
| Sombor | 0.2235 | 0.2242 | 0.9497 | 0.96221 | +0.0125 |
| ISI | 0.1739 | 0.1761 | 0.9690 | 0.98370 | +0.0147 |
| AZI | 0.2289 | 0.2316 | 0.9464 | 0.97101 | +0.0246 |
| M2 | 0.3836 | 0.3867 | 0.8504 | 0.92675 | +0.0763 |
| HM | 0.4183 | 0.4204 | 0.8232 | 0.90313 | +0.0799 |
| F | 0.4669 | 0.4683 | 0.7807 | 0.87174 | +0.0910 |
| ReZG3 | 0.6405 | 0.6447 | 0.5844 | 0.83514 | +0.2507 |

The spectral bound (3) is valid in 14/14 cases and non-vacuous in 14/14. For comparison, a Hölder bound using `‖u_i‖₁` and the range of `f` is valid but vacuous (τ > 1) for 6 of the 14 indices; the covariance structure is essential.

### 4.3 Against additive properties, no index beats counting atoms

We anchor on **molecular weight**, which is exactly additive and is not the output of any fitted model. Molar refractivity (Crippen) and Labute ASA are themselves fitted additive schemes, so correlations with them are partly circular; we report them for comparability with the literature, which uses such values routinely, but they carry no independent weight here.

| dataset | property | R²(n) | R²(best) | ΔR² |
|---|---|---|---|---|
| Drug set (38) | molecular weight | 0.9850 | 0.9826 | **−0.0025** |
| | Labute ASA | 0.9911 | 0.9930 | +0.0019 |
| | molar refractivity | 0.9564 | 0.9696 | +0.0131 |
| ESOL (1,105) | molecular weight | 0.8510 | 0.8645 | +0.0136 |
| | molar refractivity | 0.9414 | 0.9432 | +0.0018 |
| | Labute ASA | 0.9618 | 0.9521 | **−0.0096** |
| T_m set (8,583) | molecular weight | 0.7977 | 0.7941 | **−0.0036** |
| | molar refractivity | 0.8856 | 0.9128 | +0.0272 |
| T_b set (6,400) | molecular weight | 0.8430 | 0.8136 | **−0.0295** |
| | molar refractivity | 0.8233 | 0.8824 | +0.0591 |

In four of ten cases the null beats every one of the 14 indices outright. The largest gain anywhere is +0.059.

### 4.4 Against measured non-additive properties, the battery is weak

| property | N | R²(n) | R²(best of 14) | best index |
|---|---|---|---|---|
| melting point | 8,583 | 0.0728 | 0.1569 | ReZG3 |
| boiling point | 6,400 | 0.1514 | 0.2116 | harmonic |
| log S (ESOL) | 1,105 | 0.3272 | 0.3415 | GA |

On the same molecules for which additive properties reach R² = 0.79–0.91, **no single index** exceeds R² = 0.212 on a measured non-additive property.

We state the countervailing result plainly. Used *jointly*, the 14 indices do better than any of them alone. Cross-validated R² for measured melting point (N = 8,583):

| model | CV R² |
|---|---|
| heavy-atom count `n` | 0.0713 |
| molecular weight | 0.0626 |
| best single index (ReZG3) | 0.1545 |
| **`n` + ring count** | **0.2179** |
| all 14 indices | 0.3364 |
| all 14 indices + `n` + ring count | 0.3359 |

Two things follow. First, the best single index is beaten by a two-number baseline that requires no chemical graph theory at all. Second, the 14 indices jointly *do* carry information beyond simple counts (0.336 vs 0.218) — a real effect we do not dispute. But adding `n` and ring count to the battery changes cross-validated performance by −0.0005, confirming that trivial size descriptors already lie within the span of the indices. The joint model is therefore a legitimate descriptor set; the *comparative single-index ranking* that the template produces is not a legitimate inference.

### 4.5 Reported R² is a readout of sample size dispersion

Drawing subsamples of 15 molecules from ESOL at controlled size dispersion (Figure 1):

| property | CV(n) | R²(best) | R²(n) | ΔR² |
|---|---|---|---|---|
| molecular weight | 0.035 | 0.330 | 0.212 | +0.118 |
| | 0.180 | 0.530 | 0.462 | +0.068 |
| | 0.475 | 0.856 | 0.835 | +0.021 |
| molar refractivity | 0.035 | 0.336 | 0.189 | +0.148 |
| | 0.180 | 0.685 | 0.668 | +0.017 |
| | 0.475 | 0.945 | 0.940 | **+0.005** |
| log S (measured) | 0.035 | 0.253 | 0.061 | +0.191 |
| | 0.475 | 0.444 | 0.398 | +0.046 |

Nothing changes across a row except how much size variation the sample contains.

**The inversion** (Figure 2a). As the headline R² rises, ΔR² falls. For molar refractivity the reported R² climbs from 0.34 to 0.94 while the index's contribution over atom count collapses from +0.148 to +0.005. *The most impressive-looking results are those in which the descriptor is doing the least work.*

**At published sample sizes the winner is noise.** Drawing 15 molecules and reporting best-of-14:

| property | distinct winners | modal winner share | mean R² | 5th–95th percentile |
|---|---|---|---|---|
| melting point | 14 of 14 | 48% (ReZG3) | 0.250 | 0.020 – 0.583 |
| boiling point | 14 of 14 | 42% (harmonic) | 0.311 | 0.026 – 0.738 |
| log S, CV < 0.05 | 14 of 14 | 14% (harmonic) | — | — |

Figure 2b contrasts this sampling distribution with the true permutation null, obtained by pairing each 15-molecule draw with responses from unrelated molecules:

| | real 15-draws | permutation null |
|---|---|---|
| melting point, mean | 0.250 | 0.102 |
| melting point, 95th percentile | 0.576 | 0.322 |
| boiling point, mean | 0.306 | 0.116 |
| boiling point, 95th percentile | 0.730 | 0.352 |

Two conclusions, and we are careful to separate them. **(i)** Any best-of-14 R² below roughly 0.35 at N = 15 falls inside the permutation null and is not evidence of a structure–property relationship. A substantial share of published values lies in this range. **(ii)** Values well above that — such as a recently reported r = 0.836 (R² = 0.698) for boiling point — are *not* explicable as noise; they sit at the 100th percentile of the null. What they are explicable by is upward bias: at N = 15 the expected best-of-14 R² is 0.306 against a full-sample value of 0.212, and the 5th–95th range spans 0.026–0.730. A single reported value from a 15-molecule study is therefore compatible with almost any underlying effect size, and the publication filter selects the upper tail. We do not claim published results are noise; we claim they are uninformative as point estimates.

---

### 4.6 Functional form: the threshold is not fixed

Studies in the target literature frequently fit quadratic or cubic rather than linear models. We therefore repeated the analysis at polynomial degrees 1–3, fitting the *same* order to the index and to the null so the comparison stays like-for-like.

**At full sample size the verdict is unchanged, with one honest exception.** For molecular weight — our anchor additive property — ΔR² is −0.0036, −0.0036 and −0.0041 at degrees 1, 2 and 3: the null beats every index at every order, and curvature is irrelevant. For boiling point ΔR² is likewise flat (+0.060, +0.057, +0.057). For melting point, however, the indices benefit from curvature more than atom count does, and ΔR² *rises* from +0.084 to +0.112. Our linear-only analysis was mildly conservative there. In absolute terms both remain poor (R² 0.085 for the null, 0.197 for the best index at degree 3), and cross-validated values track in-sample values to within 0.001 at N = 8,583, so this is a real effect rather than overfitting.

**At N = 15 the effect is dramatic, and it runs our way.** Fitting a cubic means estimating four parameters from fifteen points, and the permutation null inflates accordingly:

| property | fit degree | null mean | **null 95th pct** | real-draw mean | real-draw 95th |
|---|---|---|---|---|---|
| melting point | 1 | 0.104 | **0.323** | 0.248 | 0.582 |
| | 2 | 0.207 | **0.475** | 0.351 | 0.664 |
| | 3 | 0.316 | **0.605** | 0.442 | 0.729 |
| boiling point | 1 | 0.117 | **0.358** | 0.311 | 0.734 |
| | 2 | 0.237 | **0.533** | 0.435 | 0.782 |
| | 3 | 0.349 | **0.654** | 0.543 | 0.848 |

The evidential threshold is therefore **not a fixed number but a function of the fit order the study used**. A paper reporting R² = 0.60 from a cubic fit to 15 molecules has produced a value that lies inside its own permutation null; the same value from a linear fit would be well outside it. Higher-order fits do not strengthen a small-sample QSPR claim — they raise the bar it must clear, and papers reporting improved R² after moving from linear to cubic at N ≈ 15 are reporting overfitting.

The recently published r = 0.836 boiling-point claim survives this at every order (null 95th percentile 0.654 at degree 3 against an observed 0.698), though at degree 3 its margin is slight.

---

## 5. Discussion

### 5.1 What is and is not claimed

We do **not** claim degree-based indices carry no information. The selection-bias null (shuffled response, max over 14 indices) gives mean R² = 0.0011 and 95th percentile 0.0040 on ESOL, against an observed ΔR² of +0.0143 — a small but real margin. Out of sample the indices genuinely help on some additive properties (molar refractivity, leave-one-out RMSE improved 16.6% over the null on the drug set).

What we claim is narrower and, we think, harder to dispute: **the comparative claims produced by the template of §1.1 are not interpretable**, because the objects being compared span one dimension, the properties most often used are size proxies, and the sample sizes used cannot distinguish among 14 collinear predictors.

### 5.2 Reconciliation with prior work

Gutman and Tošović (2013) are not contradicted; they are explained. Their fixed-`n` design sets CV = 0, makes the bound vacuous, and thereby creates the only conditions under which the indices can be told apart — which is why they were able to find real differences and report that many indices correlate weakly or not at all. The subsequent literature restored size variation, adopted additive properties, and reported the resulting size correlation as predictive success.

Randić's orthogonalization argument stands, and licenses continued interest in intercorrelated descriptors *used in combination*. It does not license single-index comparative rankings.

Roy and Gramatica's warning about R² concerns external validation and applies here with additional force, but is orthogonal to the size-proxy mechanism we identify.

### 5.3 A checklist for editors and referees

The following can be applied without specialist knowledge of chemical graph theory, and instantiates rather than extends existing QSAR validation practice. For any submission using degree-based topological indices:

1. **Report ΔR² over a heavy-atom-count null** for every index–property pair. If ΔR² ≈ 0, the result concerns molecular size and should be described that way.
2. **Report CV(n)** for the molecule set. Without it, R² is uninterpretable.
3. **State whether the property is additive.** Molecular weight, molar refractivity, molar volume and polarizability are; correlations with them are not evidence about an index.
4. **Do not crown a best index at N < 30** without a resampling analysis of winner stability.
5. **Compare R² against a permutation null computed at the study's own sample size and fit order.** Indicative 95th percentiles at N ≈ 15 (§4.6): linear ≈ 0.35, quadratic ≈ 0.50, cubic ≈ 0.63. A reported value below the applicable threshold is not evidence of a structure–property relationship. Note that raising the fit order raises this bar rather than lowering it.
6. **Report the provenance of property values** — measured, or estimated by an additive group-contribution scheme. Regressing an additive estimate on an additive descriptor is circular.

### 5.4 Limitations

The boiling-point compilation is a handbook aggregation that may include estimated as well as measured values; the melting-point set is curated experimental data and should be treated as primary. The CAS-to-structure join through PubChem has an unmeasured mismatch rate; a manual audit of a random sample would bound it. The maximum-degree ≤ 4 filter, which we had expected to bias the sample, excludes exactly one molecule of 8,584 and changes no reported quantity in the third decimal place. Our index battery covers degree-based indices only — distance-based, eccentricity-based and neighbourhood-degree-based families are not tested here, though the profile decomposition (1) suggests analogous arguments. Finally, the drug set of 38 molecules was assembled by hand from classes recurring in the literature and is not a random sample of drug space; all headline claims are replicated on the larger unbiased sets.

---

## Declaration of generative AI use

Portions of this work were carried out with the assistance of Claude, a large language model developed by Anthropic. Its contribution was substantial and is described here in full so that readers and referees can weight it appropriately.

The model was used to: survey the prior literature and identify the relevant precedents, including Gutman & Tošović (2013) and Randić et al. (1994); propose and formalise the profile decomposition (1) and the bound (2)–(3); identify and prove the corollary (4); write all analysis and figure-generation code; execute the computational experiments reported in Section 4; and draft substantial portions of the manuscript text.

The research direction, the choice of which lines of enquiry to pursue and abandon, and all decisions about scope, framing and claim strength were made by the author. Two of the model's initial hypotheses — a worst-case sandwich bound, and a claim that heavy-atom count would outperform every index out of sample — were falsified during the work and retracted; the results reported here reflect the surviving claims. An adversarial review pass, also model-assisted, identified and corrected a mischaracterised null distribution and one factual error in an earlier draft.

The author has verified the reported results, takes full responsibility for the content of this paper, and affirms that the generative AI tool is not and cannot be an author. All code and data required to reproduce every number and figure are provided as supplementary material, and independent reproduction is encouraged.

---

## 6. Reproducibility

All datasets are public. Analysis code, molecule lists and figure scripts are provided as supplementary material; the full pipeline runs in a few minutes on a laptop using RDKit and NumPy. No proprietary software or data is required to reproduce any number in this paper.

---

## Figures

**Figure 1.** Reported R² in degree-based-index QSPR as a function of the size dispersion of the sample, for two additive properties and one measured property. Subsamples of 15 molecules from ESOL. Red: best of 14 indices; black dashed: heavy-atom-count null. Shaded band: interquartile range. Orange region marks the dispersion typical of published drug sets. `size_dispersion.png`

**Figure 2.** (a) The inversion: as the headline R² rises with size dispersion, the index's improvement over the atom-count null falls toward zero. (b) Null distribution of best-of-14 R² from random 15-molecule draws, for measured melting and boiling points; dashed line marks a recently published boiling-point claim. `figure2.png`

---

## References

*(to be completed in journal format)*

- Delaney JS. ESOL: estimating aqueous solubility directly from molecular structure. *J Chem Inf Comput Sci* 2004;44(3):1000–1005. doi:10.1021/ci034243x
- Furtula B, Gutman I. A forgotten topological index. *J Math Chem* 2015;53:1184–1190.
- Gutman I. Geometric approach to degree-based topological indices: Sombor indices. *MATCH Commun Math Comput Chem* 2021;86(1):11–16.
- Gutman I, Trinajstić N. Graph theory and molecular orbitals. Total π-electron energy of alternant hydrocarbons. *Chem Phys Lett* 1972;17(4):535–538.
- Gutman I, Tošović J. Testing the quality of molecular structure descriptors. Vertex-degree-based topological indices. *J Serb Chem Soc* 2013;78(6):805–810. doi:10.2298/JSC121002134G
- Hayat S, Suhaili N, Jamil H. Statistical significance of valency-based topological descriptors for correlating thermodynamic properties of benzenoid hydrocarbons with applications. *Comput Theor Chem* 2023;1227:114259.
- Randić M. On characterization of molecular branching. *J Am Chem Soc* 1975;97(23):6609–6615. **[verify title against ACS record]**
- Randić M et al. Graphical bond orders: novel structural descriptors. *J Chem Inf Comput Sci* 1994;34(2). **[verify author list and page range]**
- Roy K, Gramatica P. Beware of R²: simple, unambiguous assessment of the prediction accuracy of QSAR and QSPR models. *J Chem Inf Model* 2014;54:1696–1701.
- Wiener H. Structural determination of paraffin boiling points. *J Am Chem Soc* 1947;69(1):17–20.
