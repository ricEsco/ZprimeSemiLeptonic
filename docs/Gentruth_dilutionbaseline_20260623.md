# Gen-Truth Spin-Correlation Reference & Reconstruction-Dilution Baseline

## **Goal context:** Standard-Model measurement of top-quark entanglement in the **boosted** semileptonic regime for the UL18 era

**Date:** 2026-06-23  
This document records what was established in the 2026-06-23 session and locks in the reconstruction-dilution **baseline** that the optimization (P1) will be judged against. The priority observable is **D̃** (the coefficient from `cHel_P3n`); **D** (from `cHel`) is kept as a sanity check.

---

## 1. Necessary Infrastructure for a gen v.s. reco Comparison Present in AnalysisTree

### 1.1 Generator Information

UHH2 provides a `TTbarGen` class that finds the generator-level ttbar system from `GenParticles` (`vector<GenParticle>`).

#### Issue

During processing, the analysis uses the `TTbarGenProducer` object to point to the `ttbargen` object via transient `get_handle` calls to use this gen-info during the analysis however, this does **not** save the `ttbargen` info to memory like `declare_event_output` would have.

#### Solution

`GenParticles` **is saved** in the AnalysisTree (confirmed via `InspectAnalysisTree.py`) and the gen-level ttbar system is rebuilt from `GenParticles`. Used `testingTTbarGen_v2.py` to confirm the available fields (properties) and indices (generator history) for each gen-particle:

- `m_pdgId`
- `m_charge`
- `m_pt/eta/phi/energy`
- `m_index`
- `m_mother1`
- `m_mother2`
- `m_daughter1`
- `m_daughter2`

Using these I created a pyROOT script to find the gen-level ttbar system post-analysis:

**`ttbargen_py.py`** — a faithful pyROOT port of UHH2 `TTbarGen.cxx`, navigating the decay chain via `GenParticle.daughter(gps,i)` / `mother(gps,i)`.

#### Validation

Running `validate_ttbargen.py` over the two muon files for UL18:

- 100 % of events classified semileptonic (muhad 92.1 %, tauhad 7.9 %, zero `notfound`)
- Checked event 0 to confirm the top quark assignment gets the charge correct (`TopLep` = antitop for a µ⁻)

>&nbsp;&nbsp;&nbsp;&nbsp; Tau handling  
>&nbsp;&nbsp;&nbsp;&nbsp; Inspected τ-decays using `tauDecayProbe.py` and found they are pruned from the gen-record, meaning 0 % of τ→ℓ chains are reachable.  
>&nbsp;&nbsp;&nbsp;&nbsp; **Chose to include tauhad events and use the τ direction as the lepton analyzer**  
>&nbsp;&nbsp;&nbsp;&nbsp; This will induce a small, ~8 % physical smear with respect to the reco µ; can be split out later.

Finally, created `gen_spincorr.py` which recreates the gen-level spin correlation variables from the gen-level ttbar system and stores the cHel & cHel_P3n histograms in `gen_spincorr_TTToSemiLeptonic_all.root`.

### 1.2 Reconstructed Information

The `SpinCorrelations` module (`ZprimeSemiLeptonicModules.cxx:1323+`) is ran during analysis from `ZprimeAnalysisModule.cxx:1248`.

This module uses `declare_event_output` and fills per-event reco branches with the **same recipe as `Hists.cxx:1896-1897`** (`evt.set(h_cHel,…)` at 1814–1815) saving the event-level information for the following variables:

- `chi2`
- `M_tt`
- `beta`
- `dyreco`
- `cosTheta1/2{k,r,n}_(anti)Lep`, `cosTheta1/2{k,r,n}`
- `C{ij}`, `C{ij}_plus/minus`
- `cHel`, `cHel_P3n`
- `Sigma/Delta_phi`
- `eventweight`
- need to add Scattering angle: `cos_PosTop_beam`, found in  `ZprimeSemiLeptonicModules.cxx:1656`
  
>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Default value `−10` for events that fail reconstruction/b-tagging.*

#### Since we have access to the event weight branch via `eventweight`, this is what makes the per-event reco-vs-gen comparison possible entirely offline, no reprocessing

---

## 2. Validation of the gen reference

Every check says the gen computation is sound:

- Using `recogen_diag.py` the values of the entanglement witnesses were calculated  **on the same exact events for both reco and gen and were found to be positively correlated** (on 141 892 matched reconstructed events):
  - **+0.42** for D
  - **+0.29** for D̃
- D(m_tt)_gen is smooth and physical (monotonic −0.53 → +0.08) and **agrees with reco near threshold** (`[0,350]`: −0.53 vs −0.44)
  - Since D is basis-independent, this validates the gen tops and rest-frame boosts.
- D(m_tt)_reco is not-smooth, this indicates reconstruction artifacts affecting the physics
- D̃(m_tt)_gen is generally smaller at low m_tt (save the `[0,350]` bin) and increases in boosted region
- D̃(m_tt)_reco is inflated in all bins (save the `[1000,1500]` bin) wrt gen-values
  - this suggests the smearing is coming from the Bernreuther basis
- Using `component_diag.py` the values of the spin-analyzer projections onto the Bernreuther basis were computed, again **on the same exact events for both reco and gen and were found to be positively correlated** (`cosTheta1/2 {k,r,n}`) within the range [0.33–0.45]:
  - The mean values of `cosTheta1/2_r` and `cosTheta1/2_n` for gen and reco match indicating they're well modeled
  - The `cosTheta1/2_` projection is the most diluted by reco:
    - ⟨cosθ1k⟩_gen = +0.092 while ⟨cosθ1k⟩_reco = +0.039 — indicates helicity-polarization washout

### → All-in-all the gen reference is trustworthy and used as truth and shows we have some optimization to do

---

## 3. Pipeline validation (offline reco branches reproduce the analysis)

The per-event reco spin-correlation values are stored as AnalysisTree branches (see §1.2), so reco can be read directly and compared to the analysis' own `extractCoeff.py` output:

| Quantity | offline (this work) | analysis `extractCoeff.py` | status |
|---|---|---|---|
| D̃_reco inclusive (`cHel_P3n`) | 0.3519 | 0.3519 | exact! ✓ |
| D̃_reco, reco `M_tt>800` | 0.1537 | 0.1539 | 0.13% ✓ |

### → Reading reco from the branches faithfully reproduces the analysis. Combined with the gen validation (§2), the **entire reco-vs-gen comparison chain is trustworthy**

---

## 4. Baseline Bias

Using the offline reco and gen variables in `dilution_map.py`, we look in the **gen-level boosted** region, `m_tt^gen > 800 GeV`, and compare the values of our observables computed from reconstruction/generator variables. The 2D dilution map and gen/reco D̃ maps in the boosted regime are stored in `dilution_map_boosted.root`.

Since the values are computed on **identical, event-weighted events** defined by the event's gen-level m_tt value, any difference in these values indicates a bias induced by our reconstruction:

| Observable | reco | gen (truth) | reco / gen |
|---|---|---|---|
| **D̃** (`cHel_P3n`) | **0.2771** | **0.2195** | **≈ 1.26** |
| D (`cHel`, sanity) | 0.0414 | 0.0444 | ~0.93 |

### The reconstruction overestimates the boosted D̃ by ~26 %

Because the bias is an *overestimate*, the current reconstruction would bias a measured entanglement value **high** — a real systematic, not benign dilution-toward-zero. Thus, we define the **P1 success criterion** be to drive this ratio toward 1.0 (and flatten the per-cell 2D ratio map toward unity).  

> The 1.26 reco/gen ratio is from a full-pass over the `TTToSemileptonic` events that pass our baseline selection and chi2 cut to contribute to our UL18 signal yield

---

## 5. Effects of Mass Migration

Using `dilution_map.py` it was also found that the reconstructed `m_tt` is substantially smeared relative to truth (dominated by the neutrino-pₙ two-fold ambiguity-- *need to clarify exactly what this means and if can be mitigated by a physically motivated solution*).  
This shows up directly in the different D̃ values computed in the equivalent reco vs gen bins:

- D̃_reco over **reco-binned** `M_tt>800` = **0.1537**
- D̃_reco over **gen-binned** `m_tt>800` = **0.2771**

Same reco-level `cHel_P3n` values — **only the event population differs** (reco-classified-boosted vs truly-boosted).  

### The large gap indicates significant `m_tt` migration across the 800 GeV boundary

 ---

## Summary

### The reco/gen ratio at fixed truth binning (0.2771/0.2195 = 1.26) isolates the coefficient bias and is the number P1 must drive to 1.0. Separately, the reco-binned→gen-binned shift (0.1537→0.2771, ×1.80) is the migration impact — same reco values, different populations — which a final measurement will absorb via a response/unfolding treatment

---
---

## A. Physics takeaways (reconstruction distortion)

- The reconstruction **heavily smears the b-side spin analyzer** — per-event reco-vs-gen correlation only ~0.4 (the lepton is well-measured; the b is the problem). This is the b-assignment issue P1 targets.
- **Low mass / resolved regime:** reco **manufactures a spurious `cHel_P3n` correlation** far above truth (reco D̃ 0.34–0.49 vs gen 0.04–0.13), from combinatorics + χ² selection sculpting + wrong-b picking. Strong argument for the boosted focus and for excluding low mass.
- **Boosted regime (>800):** the distortion is milder (~26 % overestimate) — the signal region is on much firmer ground, but still biased high and worth optimizing.
- reco D itself is **erratic vs mass** (a +0.36/+0.39 bump at 400–500 GeV) while gen D is smooth — another signature of reconstruction artifacts in the resolved regime.

---

## B. Reference numbers

**Analysis reco values** (full UL18, both flavors, `extractCoeff`, multiplier = 5):

| Histogram | coefficient |
|---|---|
| `cHel` (D) | 0.1217 |
| `cHel_Mtt300_400` | −0.0439 |
| `cHel_Mtt300_400_betaLT0p9` | −0.0477 |
| `cHel_P3n` (D̃) | 0.3519 |
| `cHel_P3n_Mtt800_Inf` | 0.1539 |
| `cHel_P3n_Mtt800_Inf_cosThetaLT0p4` | 0.3935 |

**D̃ and D vs reco `M_tt`** (matched events, reco vs gen) — reco D̃ inflated at low mass, converging in the boosted bins:

| m_tt bin [GeV] | D_reco | D_gen | D̃_reco | D̃_gen |
|---|---|---|---|---|
| 0–350   | −0.440 | −0.534 | 0.344 | 0.130 |
| 350–400 | +0.056 | −0.405 | 0.430 | 0.078 |
| 400–450 | +0.360 | −0.287 | 0.493 | 0.043 |
| 450–500 | +0.394 | −0.249 | 0.417 | 0.074 |
| 500–600 | +0.112 | −0.205 | 0.363 | 0.012 |
| 600–700 | −0.132 | −0.134 | 0.221 | 0.104 |
| 700–800 | −0.111 | −0.104 | 0.175 | 0.102 |
| 800–1000| −0.113 | −0.015 | 0.144 | 0.096 |
| 1000–1500| −0.180 | +0.082 | 0.090 | 0.150 |
| 1500+   | +0.300 | +0.077 | 0.300 | 0.186 |

(Coefficient = 5 × forward-backward asymmetry of the distribution split at 0; equivalently 5 × ⟨sign⟩.)

---

## C. Tooling produced (pyROOT, CMSSW_10_6_28 / Python 2.7)

| Script | Purpose |
|---|---|
| `ttbargen_py.py` | `TTbarGenPy` — pyROOT port of `TTbarGen.cxx` (build per event from `t.GenParticles`) |
| `gen_spincorr.py` | gen-level spin-correlation (mirror of `Hists.cxx:1582-2000`) + gen coefficient maps |
| `recogen_diag.py` | per-event reco-vs-gen comparison (correlation, sign agreement, mass-binned D/D̃) |
| `component_diag.py` | per-axis `cosTheta1/2{k,r,n}` reco-branch vs gen comparison |
| `dilution_map.py` | the baseline: validations + boosted yardstick + per-(cosΘ,m_tt) dilution map (`dilution_map_boosted.root`) |

Notes: build `TTbarGenPy(t.GenParticles)` and use it **within the same entry** (proxies point into the current event); read only `GenParticles*` + needed scalar branches via `SetBranchStatus` for speed; the `TList::Clear`/`THashList::Delete` lines at teardown are harmless.

---

## D. C++ changes from this session (await reprocessing)

These are committed but **not yet in the processed files** used for the baseline (which relied only on pre-existing branches):

- [x] Coefficient `TProfile2D` maps in `Hists.cxx`: `cHel(_P3n)_coeff_Mtt_vs_cosThetaStar` and `…_vs_beta` (filled with `5·sign(cHel)`; per-cell mean = coefficient — validated to machine precision against `extractCoeff`). M_tt binning extended to 200–2000 GeV (22 var bins; includes a 300 edge).
- [x] `deltaR_hadTop_bGen` fill moved out of the merged-only block → now fills **both** resolved and merged topologies.
- [ ] 1D `cos_ThetaStar` (Bernreuther production angle) histogram.

---

## E. Open items / next steps

1. [x] **Full pass** of `dilution_map.py` (`NMAX=-1`) to tighten the boosted reco/gen ratio and fill the per-cell `Dt_ratio_boosted` map with usable statistics.
2. [ ] **P1 — b-tag-aware χ² candidate selection** (the optimization). Success = boosted **D̃ reco/gen → 1** and the ratio map flattening toward unity.
3. [x] **Reprocess** with the latest code so the new maps, the resolved `deltaR_hadTop_bGen`, and the β maps are available on disk.
4. [ ] Optional: investigate the low-mass `cHel_P3n` inflation via the per-term ⟨c1ᵢ·c2ᵢ⟩ reco-vs-gen comparison
