# Control and systematics plots

Three plotting scripts with official CMS styling (`tdrstyle.py` and `CMS_lumi.py`):

| Script | What it shows | Data |
|---|---|---|
| `plotPreselectionKinematics.py` | Stacked MC of the Preselection kinematics at one selection step | never |
| `plotAnalysisControlPlots.py` | Stacked MC (+ data) of the Analysis kinematics at one selection step | `--data` |
| `plotSystematicVariations.py` | Every up/down variation against the nominal, for the SC observables | MC only |

Run them inside an environment with PyROOT, e.g. after `cmsenv` in the Combine release.

## Common features

- **Input:** the per-job SFrame outputs are summed directly, so no `hadd` is needed.
  - Files are opened one at a time, and only the requested histograms are copied out.
  - The size of the `AnalysisTree` in the Preselection files therefore doesn't matter.
- **Cache:** the summed histograms go into a small `cache*.root` file in the output folder.
  - A rerun, for example after a style or binning change, reads only the cache.
  - The scripts re-read the job outputs on their own when a run needs histograms that aren't cached yet; `--rebuild` forces a fresh read anyway.
- **Output layout:** `plots/<year>/<script>/...`, relative to where you run the script (`--outdir` changes it). Every plot is saved as `.pdf` and `.png`. One multi-page PDF per run is written at the end, so that only one PDF file is ever open.
- **Labels:**
  - "CMS" and "Private Work" sit above the frame on the left (or "Simulation Private Work" when no data is drawn).
  - The luminosity sits above the frame on the right.
  - `--year` (default `UL18`) sets the input paths and the luminosity label. The luminosity values are the ones in `CMS_lumi.py`.
- **Legend:** spread over 3 columns across the top of the frame, with the channel and step written below it.
- **Uncertainty band:** MC stat only.
- **Y axis:** linear by default; `--logy` switches to a log scale.
- **Stack groups:** tt̄ (ℓ+jets), tt̄ (2ℓ), tt̄ (had.), single t, W+jets, Drell–Yan, diboson, QCD. They use the CMS-recommended 10-colour palette.
- **Skipped histograms:** a histogram whose axis range differs between job files is skipped with a warning. This happens when it was booked with `xlow == xup`, which tells ROOT to choose the range from the data, separately in each job; currently that's `phi_mu1` and `phi_mu2`.

### Binning (`plotConfig.binning(name, source)`)

`binning()` returns `(nbins, lo, hi)`, or `None` to keep the booked binning. `source` is one of:
- `"preselection"`: the Preselection plots;
- `"analysis"`: the Analysis control plots;
- `"sc"`: the systematics plots.

The same histogram name can be booked differently in the Preselection and the Analysis. For example, `pt_mu` is booked with 50 bins on 0–1500 in one and 90 bins on 0–900 in the other.

Because the histograms are already binned, a new binning can only:
- merge booked bins, and
- narrow the range.

Every new bin edge has to fall exactly on a booked edge. If one doesn't, that plot is skipped with a message that gives the booked binning, and the run carries on. Example:

```python
def binning(name, source):
    if source == "analysis" and name.startswith("pt_mu"): return (45, 0., 900.)   # booked (90, 0, 900)
    ...
    return None
```

The `"sc"` default (`None`) merges bins by `SC_REBIN = 4` (24 → 6), like the Combine input.

### Under/overflow (`plotConfig.fold_overflow(name)`)

- **Open-ended quantities** (pt, masses, MET, ST, multiplicities): under- and overflow are added to the first and last bin.
- **Bounded quantities** (cos*, C*, cHel*, η, φ, Σφ/Δφ, tagger scores, β, τ): never added. Their under- and overflow only hold the sentinel values (−10 / 99) of events where the observable isn't defined.

To change the defaults:
- `BOUNDED_PATTERNS` holds the name patterns of the bounded quantities.
- `OVERFLOW_OVERRIDE` sets the behaviour for one histogram name.
- `--no-overflow` turns the folding off everywhere, in all three scripts.

## 1. Preselection kinematics

```bash
python3 plotPreselectionKinematics.py --list-steps                     # available <step>_General folders
python3 plotPreselectionKinematics.py [--step MET_General] [--logy] [--no-overflow] [--plots pt_jet MET]
```

- Plots the distributions in `PRESELECTION_KINEMATICS`.
- Simulation only: the Preselection isn't split by lepton channel, and the corrections applied in the Analysis aren't applied yet.
- The ratio panel shows the relative MC statistical uncertainty.

## 2. Analysis control plots

```bash
python3 plotAnalysisControlPlots.py --list-steps
python3 plotAnalysisControlPlots.py [--channel muon electron] [--step EndOfBaselineSelection_General] [--data] [--logy]
```

- **Distributions:** those in `ANALYSIS_KINEMATICS`, from `<step>_General`, one per plot.
- **Data:** drawn only with `--data`. It is never drawn for generator-level quantities (`GEN_ONLY`). For the fit observables (`BLINDED`, at any step), data also needs `--unblind`.
- **Ratio panel:** data/prediction, or the relative MC stat uncertainty when data isn't drawn.

## 3. Systematic variations vs nominal

```bash
python3 plotSystematicVariations.py [--channel muon] [--region Merged] [--observables 'cos_theta1*' Ckk] \
                                    [--nuisances 'btag_*' pdf] [--process TTbar] [--no-overflow]
```

- **What is plotted:** one plot per region × SC observable × Combine nuisance. Up is blue and down is red, with the change in the integral in the legend. The lower panel shows variation / nominal, with the nominal MC stat uncertainty in grey.
- **Scope:** by default, the sum of all MC. A nuisance that acts on some processes only (`muR_ttbar`, `NoSC_f`, ...) shifts just those, as in the fit. `--process` restricts the plots to one datacard process.
- **Histograms:** the raw histograms from the `<region>_SC_*` folders. PDF is shown as the Hessian ±1σ of the 100 replicas.
- **Output:** files are named after the Combine nuisances (e.g. `btag_lfstats1_UL18.pdf`). There is one folder and one multi-page PDF per region and observable.
- **Size of a full run:** about 36 observables × ~40 nuisances per channel. `--observables` and `--nuisances` (wildcards allowed) narrow it down.

## Notes

- **JEC/JER:** these come from separate runs and aren't included in the systematics plots yet.
- **`no dictionary for class ...` warnings:** harmless. The input files were written with a different CMSSW version than the one loaded.
- **`CMS_lumi.py` changes:** two small changes for Python 3 and newer ROOT versions: `iPosX//10` and explicit ROOT imports.
- **No THStack:** the stack is drawn as cumulative histograms instead of a `THStack`, because a `THStack` drawn over an axis frame crashes in some recent ROOT builds. The result looks the same.
