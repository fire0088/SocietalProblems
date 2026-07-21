DARK-ENERGY INVESTIGATION — DOCUMENT BUNDLE
============================================

Two documents, both compile-ready LaTeX with matching PDFs.

paper/
  de_ledger.tex / .pdf
    "A Standing Ledger of Dark-Energy Positions"
    The research-result document: the tiered ledger (Supported / Open /
    Conjectural) of positions on the DESI DR2 phantom-crossing signal,
    with the pre-registered DR3 test, AP-basis analysis, executive
    summary, figures, and a verified bibliography.

methodology/
  methodology.tex / .pdf
    "Anatomy of an AI-Assisted Investigation"
    The meta / methodology document: how the investigation was conducted
    (chronological trajectory + reversals), AI-vs-human comparative
    analysis, limitations & reproducibility statement, and an adversarial
    red-team / threat-model section.

figures/
  fig_stack.pdf        — the ~4-sigma signal decomposed into weak pieces
  fig_leaveoneout.pdf  — per-tracer leave-one-out significances
  fig_apbasis.pdf      — LRG residuals in the Alcock-Paczynski basis
  fig_dr3test.pdf      — pre-registered DR3 test with systematic floor
  (Both .tex files expect the figure PDFs in the same directory as the
   .tex when compiling — see note below.)

TO COMPILE
----------
Each .tex expects the four figure PDFs alongside it. Either copy the
figures into paper/ and methodology/ before compiling, or adjust the
\includegraphics paths. Then:
    pdflatex de_ledger.tex      (run twice for refs)
    pdflatex methodology.tex    (run three times for TOC + refs)

CAVEATS (see the documents' own caveat boxes for detail)
--------
- Absolute significances come from a simplified pipeline (compressed CMB
  priors, magnitude proxies) and are method-conditional; robust content
  is in directions, relative shifts, and reproduced structure.
- Bibliographies were verified this session against primary sources, but
  a few 2026 volume/page numbers and one arXiv ID should be re-checked
  against ADS before any submission.
- methodology.tex contains intentional [PLACEHOLDER] fields (repo URL,
  environment, seeds, model version) to fill on release.
