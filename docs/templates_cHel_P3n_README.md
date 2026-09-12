# cHel_P3n (D̃) reweighting templates

First look at the boosted semileptonic **D̃** templates built with the spin-correlation
**reweighting method** (our `docs/reweighting_method.md`, following your dilepton analysis
[arXiv:2406.03976](https://arxiv.org/abs/2406.03976)). Produced by
`make_templates_cHel_P3n.py` from the per-event skim + the inclusive UHH2 strength map.

## What it does

For every event in the boosted-central **reco** region it looks up
`α(m_tt, cosθ*)` from the strength map using the event's **gen** production kinematics,
reads its **gen** observable `O = gen_cHel_P3n`, forms the per-event factor

```
w(f) = (1 − α·O·(1−f)) / (1 − α·O)
```

and fills the **reco** observable with `eventweight · w(f)`. Three templates:

| template | f | meaning | ROOT name |
|---|---|---|---|
| nominal | 0 | SM (factor = 1) | `cHel_P3n` |
| up | +1 | spin correlation **removed** (NoSC) | `cHel_P3n_NoSCUp` |
| down | −1 | spin correlation **enhanced** | `cHel_P3n_NoSCDown` |

The lookup and O are pure gen quantities (inclusive gen collection, so migrated events pull the
weight from their true gen bin); the fill is at reco. → `templates_cHel_P3n_boostedCentral.root`
feeds Combine directly (morph parameter f; calibration maps fitted f → physical D̃).

## Choices made (please sanity-check)

- **Region:** boosted + central, `m_tt > 800 GeV` **and** `|cosθ*| < 0.4`. Encoded by the
  pre-windowed reco branch `reco_cHel_P3n_Mtt800_Inf_cosThetaLT0p4`.
- **Fit variable / binning:** `reco_cHel_P3n`, **6 bins on [−1, 1]** (Fig. 6 style).
- **f = ±1** as the first up/down pair (within the |f|≤1 validity range).

## The prefactor — two distinct roles (resolved)

From the physical PDF of the opening-angle cosine `O`,
`(1/σ) dσ/dO = ½(1 + α_a·α_b·D·O)`, with spin-analyzing powers `α_a·α_b = (1)(0.4) = 0.4` for our
lepton + hadronic-b pair (the dilepton pair is `(1)(−1)`, |·| = 1 — which is exactly why the
dilepton team never has to separate the two roles below; we do).

- **Map (extraction):** `α_map = 5·A = D`, where `A = (N₊−N₋)/(N₊+N₋) = α_a·α_b·D/2`. The `5 =
  2/(α_a·α_b)` inverts the analyzing-power dilution to recover the *true* D̃. ✔ keep for the map.
- **Weight (density):** the reweighting density must be the *actual PDF* of `O`, whose slope is the
  **physical** `λ = α_a·α_b·D = 0.4·D = 0.4·α_map = 2·A` — prefactor **2**, not 5.

The script therefore builds `d = 1 + λ·O` with `λ = 0.4·α_map`. This is the only choice that (a)
flattens the distribution at `f=+1` and (b) stays positive: `|λ| ≤ 0.4 < 1 ⇒ d ≥ 0.6`. Prefactor 5
would put the weight-hyperbola pole inside `|O|≤1` (your red Desmos curves) → divergent/negative
weights.

**Gen-closure proof (`gen_closure_cHel_P3n.py`).** Directly measures whether `w|_{f=+1}` flattens
the gen `cHel_P3n` slope, and scans the prefactor. On a known modulation the scan is unambiguous:

```
 p(on A)   density        A_up(f=1)     verdict
   1       0.20·α_map      +0.052       under-flattened
   2       0.40·α_map      +0.001       slope FLATTENED   <== correct (= α_a·α_b)
   3       0.60·α_map      −0.051       over-flattened
   5       1.00·α_map      −0.162       inverted; density pole enters |O|≤1
```

Run it on the real skim to get the data version of this table for the dilepton group.

## How to run (EL7 container, bare pyROOT — no UHH2 dicts)

```bash
cd .../ZprimeSemiLeptonic/src/entanglementScripts/reweighting
# 1) the templates for Combine
python make_templates_cHel_P3n.py \
       spincorr_skim_TTToSemiLeptonic_UL18.root \
       strengthmaps_uhh2/NoSC_strengthmap_lb.root
# -> templates_cHel_P3n_boostedCentral.root  (3 TH1F for Combine)
# -> templates_cHel_P3n_boostedCentral.pdf   (overlay + ratio)

# 2) the prefactor / slope-flattening proof (inclusive gen, ~2M sampled events)
python gen_closure_cHel_P3n.py \
       spincorr_skim_TTToSemiLeptonic_UL18.root \
       strengthmaps_uhh2/NoSC_strengthmap_lb.root
# -> gen_closure_cHel_P3n.pdf  + the printed prefactor scan (A_up -> 0 at p=2)
```

Template diagnostics: yields, per-template D̃ estimator with the ordering check, min density,
guard/negative-weight counts, low-stat map-cell count. Gen-closure diagnostics: nominal gen A/D̃ and
the `A_up` prefactor scan.

## Open items for after this look

- Extend to the full `f` grid Combine needs (e.g. −1…+1 in steps) once the ±1 shapes look right.
- Systematic variations: the skim carries all 146 correction weights — each template can be
  reproduced under any of them for the datacard.
- Same machinery works for `cHel` (D) and the full `C_ij`/`B_i` set — only the map object changes.
