# Two Regimes of Metastatic Organotropism: Corroborating and Extending Biophysical and Network Accounts with a Venous-Dominance Criterion

**Eric Schultz** — Independent Researcher — ORCID: 0009-0006-6283-1696

> **AI-assistance disclosure.** An AI system (Anthropic's Claude) played a significant role in this work under the author's direction: it implemented the compartment model and analysis code, executed the numerical experiments, retrieved and collated the literature and clinical-incidence data, and drafted the manuscript text. The author conceived the research direction, made all substantive judgments, verified the sources, and takes full responsibility for the content. The AI is not an author and is not credited as a source, consistent with COPE and arXiv authorship guidance.

*Exploratory note. This work does not claim a new framework; it provides independent corroboration of two recent lines of work — the biophysical/mechanical account of organotropism [6] and network flow models of metastatic spread [7] — and adds three modest increments: (i) convergent support, from an independent algebraic model, for the inverse flow–retention relationship reported by [7]; (ii) a decomposition of "mechanical routing" into specific valveless venous first-pass conduits (Batson's, portal, caval) that resolves the axial-bone gradient generic flow models do not address; and (iii) an anatomically-computable dominance ratio D that operationalizes, as an explicit switch, the qualitative claim that mechanical and molecular determinants are complementary [6]. Not a validated predictive model; scope and limitations at the end.*

## Introduction and biological context

Metastatic organotropism — the reproducible tendency of each primary cancer to colonize particular distant organs — has been debated since Paget's 1889 "seed and soil" hypothesis (cells grow only where the microenvironment is compatible). The competing view (Ewing) held that organotropism reflects mechanical routing by the circulation. A century of evidence supports both; the field treats them as complementary. A third, underweighted fact is metastatic inefficiency: most tumor cells reaching an organ never colonize it [1], so **arrival, retention, and colonization are distinct steps** that need not share a dominant driver — indeed, the canonical statement of the field is that mechanical factors govern delivery while molecular interactions govern the probability of growth [1].

**Relation to prior work.** Azubuike & Tanner [6] revisit the seed/soil and mechanical hypotheses jointly, arguing they are complementary and that endothelial mechanobiology regulates organ selection at extravasation — leaving *which* determinant dominates largely qualitative. Singh & Jacobs [7] present a PDE-on-network flow model reproducing gut/liver tropism and an inverse velocity–concentration relationship, using general cardiac-output flow without distinguishing valveless venous first-pass conduits or the molecular-subtype axis. This note corroborates and extends both: independent (algebraic) corroboration of the flow–retention law of [7]; decomposition of mechanical routing into named venous conduits addressing the axial-bone gradient [7] does not model; and the dominance ratio D as an explicit, computable form of the complementarity [6] describe qualitatively.

What has been missing is a principled account of *when* each mechanism dominates. Some tropisms follow vascular anatomy with striking fidelity (prostate → axial skeleton via Batson's valveless vertebral venous plexus [2]; colorectal → liver via portal vein). Others track molecular biology (breast bone-vs-liver depends on ER/HER2 subtype, via colonization programs such as ER-regulated SCUBE2 signaling that primes the osteoblast niche) [3].

**Contribution and its limits.** One conceptual proposal, one tested prediction. The proposal: these behaviors are two *regimes* of a single process — separable along the arrival/retention/colonization axis — distinguished by whether the primary's venous drainage is dominated by a single first-pass conduit. A minimal compartment calculation makes the distinction explicit and generates predictions; it is not itself a quantitative predictor. Stated once here; not re-hedged throughout.

**Abbreviations.** BBB, blood–brain barrier; CTC, circulating tumor cell; ER, estrogen receptor; HER2, human epidermal growth factor receptor 2; HR, hormone receptor; IVC, inferior vena cava; RCC, renal cell carcinoma; SEER, Surveillance, Epidemiology, and End Results; TNBC, triple-negative breast cancer.

## Core claim and the dominance criterion

- **Routing-dominated regime.** A single venous first-pass conduit carries most disseminating cells to one bed; venous topology sets destination, soil effects secondary. Examples: prostate → axial bone (Batson's), gut → liver (portal), kidney/sarcoma → lung (caval).
- **Colonization-dominated regime.** Venous inflow distributed across comparable routes; subtype-specific molecular programs drive destination. Example: breast (bone-vs-liver tracks ER/HER2 subtype).

**Making the criterion non-circular.** The framework is meaningful only if regime membership is assignable *independently* of the outcome it explains — otherwise classifying prostate as routing-dominated *because* it seeds bone is circular. Define a dominance ratio from anatomy alone:

    D = (venous outflow fraction through the single largest first-pass bed)
        / (summed outflow fraction through all other beds)

computed from published venous-drainage anatomy, with no reference to metastatic sites. Provisionally, D ≳ 2 denotes the routing regime and D ≲ 1 the colonization regime, with an explicit intermediate band. **This threshold is a stipulation, not a fitted value**; fixing it against an outcome-independent anatomical dataset is required future work, and near-boundary assignments are provisional. Stating the criterion explicitly converts the framework from post-hoc relabeling into a falsifiable one.

## Model structure (a framing device, stated honestly)

The compartment calculation instantiates the regime idea and generates predictions; it is not a fitted predictor, and outputs are used for ranking only, never magnitude. Seven compartments (lung, liver, axial bone, appendicular bone, brain, adrenal, lymph node); burden M_o satisfies

    dM_o/dt = f_o (1 − M_o/K_o) − c_o M_o

**The compartments do not interact**, so each equation is independent with closed-form steady state M_o* = f_o K_o / (f_o + c_o K_o). No integration required; the dynamical notation is expository. A coupled variant (finite CTC pool) changed nothing material (breast bone:liver 1.65 → 1.52). The model captures relative delivery and retention, deliberately not inter-site competition, sequential spread, dormancy, or CTC survival.

**Inflow — three anatomical edge types**, all tropism-independent: (1) arterial perfusion (background); (2) venous first-pass routes (Batson's → axial bone, portal → liver, caval → lung; Batson's valveless, pressure-gated); (3) barrier-crossing competence (per-tumor gate, e.g. BBB protease competence — not an attraction term). An earlier hand-coded-affinity version was discarded as circular.

**Acknowledged gap: lymphatic routing.** The framework is venous-centric, yet the compartment set includes lymph node and clinical data show substantial nodal involvement (e.g. breast). Lymphatic drainage is a third routing channel not represented in the mechanism; the lymph-node compartment is a passive sink. Extending D to lymphatic conduits is a needed generalization, not done here.

## Evidence, and what it does and does not show

**Ablation (within-model).** Zeroing venous first-pass terms drops top-site accuracy 6/7 → 1/7, all-correct 76% → 0% across 3000 draws. This establishes the venous terms carry the discriminating signal *within this model* — nearly forced, since they are the only per-cancer inputs (arterial perfusion, clearance, capacity are shared). It is an internal consistency check, **not independent evidence that venous routing dominates in vivo**; that stronger claim is tested for no cancer here (see P5).

**Mechanism separation.** Under ablation, prostate→bone and gut→liver collapse while melanoma→brain is unchanged (93%→93%) — brain tropism depends on barrier-crossing competence, not routing.

**Calibration.** Top-site correct 6/7; cosine 0.94 (prostate, autopsy), 0.91 (gastric, SEER), on the representable subset and thus optimistic. Ranking/ablation conclusions do not depend on these.

## The discriminating prediction and the marrow-volume confounder

The routing regime predicts that in Batson-routed prostate carcinoma, vertebral involvement declines monotonically from the lumbar entry upward. Bubendorf (n=1,589): ~97% (lumbar) → ~38% (cervical). Multiple decay rates fit comparably; the claim is directional.

**Countering the red-marrow-volume explanation.** Vertebral metastasis might track red-marrow volume (greater lumbar) — a seed-and-soil alternative needing no retrograde flow. Three arguments against marrow volume as *sole* driver: (1) red-marrow fraction doesn't fall as steeply/monotonically as the metastasis gradient; (2) most-cited-as-decisive, the gradient is *pressure-dependent* — Coman–DeLong [4] raised lumbar mets to ~70% with abdominal pressure vs ~none without, and marrow volume is pressure-invariant. **We flag this keystone rests on mid-20th-century animal work whose modern replication we have not established; discount accordingly until reproduced.** (3) Bubendorf report an inverse spine/lung relationship and earlier spine seeding. Not exclusive; the claim is only that the pressure route is needed for the full pattern.

**A deeper confound: drainage and soil may be correlated.** Venous first-pass beds and compatible niches need not be independent — a tumor chronically seeding an organ may co-evolve receptors for it, so "routing" and "soil" could be two readings of one coupled process. If so, the two-regime distinction is a useful descriptive axis rather than a separation of independent causes. D is defined from anatomy precisely to keep the routing axis measurable even if the two are entangled, but regime independence is an assumption, not established.

## Regime boundary: breast cancer, and an honest reading of P4

Breast is the model's one top-site failure, locating the regime boundary: bone-vs-liver is subtype-governed, and all subtypes share venous drainage, so routing can't explain a subtype-dependent split.

**What P4 does and does not establish.** One SEER cohort [5] (Xiao et al. 2010–2013; population incidence, one denominator, all four subtypes, no estimates): bone:liver varies 3.5-fold (CV 0.50), tracking ER→bone / HER2→liver (ER+ mean 2.85 vs ER− 1.37; HER2-enriched most liver-shifted). **Explicit logical status:** this confirms breast bone-vs-liver is subtype-driven — the *mainstream* seed-and-soil expectation — and refutes only a routing-*only* account of breast that few defend. **P4 is consistency evidence locating the colonization regime, not a test of the framework's novel content.** The novel claim — routing *dominates* in prostate/gut/kidney — has no comparably rigorous outcome-independent test here and is the principal open item. The paper's rigor currently sits where the stakes are lowest.

| Breast subtype | Bone (%) | Liver (%) | Bone:Liver |
|---|---|---|---|
| HR+/HER2− (luminal A-like) | 3.1 | 0.8 | 3.88 |
| HR+/HER2+ (luminal B-like) | 5.1 | 2.8 | 1.82 |
| HR−/HER2+ (HER2-enriched) | 4.6 | 4.2 | 1.10 |
| HR−/HER2− (TNBC) | 2.8 | 1.7 | 1.65 |

*Table 1. Subtype-specific breast metastasis and bone:liver ratio, single SEER cohort (Xiao et al., Oncotarget 2017; population incidence). Identical venous drainage across subtypes means the 3.5-fold variation cannot be a routing effect; it tracks receptor subtype, locating the colonization regime. All values same cohort/denominator; none estimated.*

## Testable predictions

- **P1.** Axial-bone concentration scales with venous-pressure exposure (Coman–DeLong, pending modern replication).
- **P2.** Liver-vs-lung first-pass dominance follows portal-vs-caval drainage fraction.
- **P3.** Brain tropism correlates with tumor protease/BBB-disruption capacity, not any venous measure.
- **P4 (tested; consistency evidence).** Breast bone-vs-liver tracks ER/HER2 subtype, not anatomy (Table 1) — locating the colonization regime, not validating routing.
- **P5 (the key untested claim).** Across cancers, the anatomically-defined dominance ratio D predicts top metastatic site in the routing regime *independently of molecular markers*. This is the discriminating test the framework most needs and does not yet have.

## Significance, generalization, and limitations

**Transferable idea.** Beyond metastasis, the reusable content is decomposing tissue localization into separable *arrival*, *retention*, and *colonization* operators, with different diseases dominated by different operators — applicable to infection tropism, embolic disease, drug biodistribution. D is *computable today* from existing venous-drainage anatomy with no new data, making regime assignment an immediately testable, low-cost hypothesis.

**Translational hook.** If the regime concept holds even directionally, it argues for stratifying metastasis-*prevention* trials by regime: anti-niche agents (bisphosphonates, denosumab, SCUBE2/Hedgehog inhibitors) predicted to benefit colonization-regime tumors more than routing-regime ones, reframing some regime-mismatched negative results as design artifacts. A concrete, fundable design, offered as hypothesis.

**Limitations.** (i) D-threshold stipulated, not fitted; near-boundary assignments provisional. (ii) Model is rank-only, uncoupled, seven-node; tails (esp. bone) inflated; aggregate breast over-predicts bone. (iii) Within-model ablation is near-forced, not in vivo evidence. (iv) P4 confirms mainstream subtype view, not the novel routing claim (untested, P5). (v) Lymphatic routing unmodeled. (vi) Routing and soil may be confounded. (vii) Cross-dataset comparisons mix autopsy (end-stage) and SEER (clinically detected), used for rank only; P4 avoids this within one SEER cohort. (viii) Small sample (seven cancers). (ix) The reference list is limited to the primary anatomical and epidemiological sources directly supporting the argument [1–5]; a full submission would expand the molecular-mechanism citations (e.g. SCUBE2, CXCR4) beyond the single subtype-pattern reference used here. Framing is a structured hypothesis, supported in part, not confirmed; no clinical use implied.

## References

1. Chambers AF, Groom AC, MacDonald IC. Dissemination and growth of cancer cells in metastatic sites. *Nat Rev Cancer*. 2002;2(8):563–572. doi:10.1038/nrc865.
2. Batson OV. The function of the vertebral veins and their role in the spread of metastases. *Ann Surg*. 1940;112(1):138–149. doi:10.1097/00000658-194007000-00016.
3. Kennecke H, Yerushalmi R, Woods R, et al. Metastatic behavior of breast cancer subtypes. *J Clin Oncol*. 2010;28(20):3271–3277. doi:10.1200/JCO.2009.25.9820.
4. Coman DR, DeLong RP. The role of the vertebral venous system in the metastasis of cancer to the spinal column: experiments with tumor-cell suspensions in rats and rabbits. *Cancer*. 1951;4(3):610–618. doi:10.1002/1097-0142(195105)4:3<610::AID-CNCR2820040312>3.0.CO;2-Q.
5. Wu Q, Li J, Zhu S, et al. Breast cancer subtypes predict the preferential site of distant metastases: a SEER based study. *Oncotarget*. 2017;8(17):27990–27996. doi:10.18632/oncotarget.15856.
6. Azubuike UF, Tanner K. Biophysical determinants of cancer organotropism. *Trends Cancer*. 2023;9(3):188–197. doi:10.1016/j.trecan.2022.11.002.
7. Singh K, Jacobs BA. A network based model for predicting spatial progression of metastasis. *Bull Math Biol*. 2025;87:60. doi:10.1007/s11538-025-01441-1.
