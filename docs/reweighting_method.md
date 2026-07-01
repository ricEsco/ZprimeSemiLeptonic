# The Spin-Correlation Reweighting Method

This note explains the reweighting used to produce the up/down template variations
for the spin-correlation measurement, with the relevant code. It works on
**generator-level** quantities to build the weights, which are then applied to the
reco-level templates.

---

## 1. What it does, end to end

The measurement extracts a spin-correlation observable by comparing data to
templates in which the spin-correlation content of the signal is varied up and
down. Reweighting is how those variations are made: it changes the **weight** of
each simulated event so the sample follows a scaled spin-correlation hypothesis,
**without regenerating any events**. The kinematics stay frozen; only the
spin-correlation modulation moves if and only if the f-values remain within the range of validity, namely between -1 and 1. If larger |f| values are used one must inspect the effect of these values on the kinematics.

The observable being varied is **not fixed to one choice**. It can be any
spin-correlation component we want a measurement of: the diagonal correlations
`c_kk`, `c_rr`, `c_nn`, the helicity opening angle `cos φ` (`gen_ll_cHel`), and so
on.

The full chain for one observable $O$:

1. **Build the strength map** $\alpha(m_{t\bar t}, \cos\theta^*)$ from the
   generator-level distribution of $O$ — Section 3.
2. **Generate the weight vector** on the generator-level events: for every event,
   combine $\alpha$ (looked up from its gen production kinematics) with the gen
   value of $O$ and a variation fraction $f$ — Section 4.
3. **Apply the same weight vector at reco**: fill the reco-level template histogram
   with these weights, giving the up/down reco templates — Section 5.
4. **Fit in Combine**: nominal + up + down templates define a morphing parameter;
   the fit extracts its value $f$.
5. **Interpret via calibration**: a calibration equation maps the fitted $f$ to the
   physical value of the observable.

The essential point: the weight is **generated from the generator-level
observable**, then that identical per-event weight is **applied to the reco fill**,
so the gen-level spin variation is carried into the reco template that Combine sees.

---

## 2. The quantity being varied

In the helicity basis $\{\hat k, \hat r, \hat n\}$ the spin-dependent production
density of the $t\bar t$ system is

$$
\frac{1}{\sigma}\frac{d\sigma}{d\Omega_1\, d\Omega_2}
\;\propto\;
1 + B^{+}_i\,\cos{\theta}^{(1)}_i + B^{-}_i\,\cos{\theta}^{(2)}_i
- C_{ij}\,\cos{\theta}^{(1)}_i\,\cos{\theta}^{(2)}_j ,
\qquad i,j \in \{k,r,n\},
$$

with $B^{\pm}$ the polarization coefficients and $C_{ij}$ the $3\times3$
spin-correlation matrix, and $\cos{\theta}^{(1,2)}$ the direction cosines of the top/antitop
decay products along the basis axes. The measured observables are components of
this structure — e.g. `c_kk`, `c_rr`, `c_nn` are the diagonal correlations, and
`cos φ` (`gen_ll_cHel`) is the combination $-(c_{kk}+c_{rr}-c_{nn})$ up to
convention.

The coefficients depend on the production kinematics $(m_{t\bar t},\cos\theta^*)$,
where $\cos\theta^*$ is the top scattering angle in the $t\bar t$ frame. That
dependence is exactly what the $\alpha$ map captures, so the reweighting acts
correctly bin-by-bin in production phase space.

Helicity basis: $\hat k$ = top direction in the $t\bar t$ CM frame; $\hat n$ =
normal to the production plane; $\hat r$ completes the right-handed set.

---

## 3. Building the strength map $\alpha(m_{t\bar t}, \cos\theta^*)$

For an observable $O$, the spin-correlation strength as a function of production
kinematics is a 2-D map $\alpha(m_{t\bar t}, \cos\theta^*)$, extracted from the
gen-level distribution of the **sign** of $O$. The observable axis is two bins split
at zero (`[-1, 0, 1]`), so each $(m_{t\bar t}, \cos\theta^*)$ cell holds $N_+$ and
$N_-$, and

$$
A = \frac{N_+ - N_-}{N_+ + N_-}, \qquad \alpha = (\text{prefactor})\times A .
$$

The prefactor depends only on the observable family, so the **same machinery builds a
map for any component** — `c_ij`, `cos φ`, the polarization terms, etc.:

The 2-bin observable axis is then collapsed and the map is stored as a 2-D
histogram of $\alpha$ (value and variance) over $(m_{t\bar t}, \cos\theta^*)$.

The only inputs are the gen distributions of $m_{t\bar t}$, $\cos\theta^*$, and the
observable.

---

## 4. Generating the weight vector (generator level)

With the map in hand, the weight vector is one per-event factor. For each event:
look up $\alpha$ from its gen production kinematics, read its gen observable value
$O$, form the nominal density $d = 1 - \alpha O$ and the varied density
$d_f = 1 - \alpha O (1-f)$, and take the ratio:

$$
w \;\longrightarrow\; w \cdot
\frac{1 - \alpha(m_{t\bar t},\cos\theta^*)\, O \,(1-f)}
     {1 - \alpha(m_{t\bar t},\cos\theta^*)\, O}.
$$

The lookup uses the gen production variables and the weight uses the gen observable
value, so the **entire weight vector is a generator-level quantity**. The variation
fraction sets which template you get: $f = 0$ leaves the weight unchanged (nominal),
$f = +x$ is the up variation, $f = -x$ the down.

---

## 5. From weights to reco templates and the fit

The same per-event weight vector then multiplies the nominal event weight when the
**reco-level** observable histogram is filled.

This gives three reco templates:

- $f = 0$ → **nominal**
- `NoSC_percent = +x` → **up**
- `NoSC_percent = -x` → **down**

They go to Combine as the morphing inputs. The fit floats the associated parameter
and returns its best-fit value $f$.

A **calibration equation** — the known (linear)
relation between the injected variation and the resulting template — maps $f$ to the
physical value of the observable.

---

## 6. Binning and config

The production axes use a fixed binning; the observable axis is always the 2-bin
sign split `[-1, 0, 1]`.

- $m_{t\bar t}$: 52 uniform bins, **300 → 2000 GeV** (53 edges).
- $\cos\theta^*$: 48 uniform bins, **−1 → 1** (49 edges).

The map is defined by a gen-only config of this form (here the observable is
`gen_ll_cHel`, but it can be any component):

```yaml
variables:

gen_variables:
  gen_ttbar_mass:
    - [300.0, 332.69, 365.38, ... , 1967.31, 2000.0]   # 53 edges
  gen_top_scatteringangle_ttbarframe:
    - [-1.0, -0.958, -0.917, ... , 0.958, 1.0]         # 49 edges
  gen_ll_cHel:
    - [-1, 0, 1]                                        # sign split
```

Units / ranges: $m_{t\bar t}$ in GeV; $\cos\theta^* \in [-1,1]$; all component
observables (`c_*`, `cos φ`, polarization projections) dimensionless and bounded in
$[-1,1]$.

---

## 7. The map-building script, step by step

`generate_nosc_reweight.py` produces `NoSC_reweight.root`. It is observable-agnostic
— it loops over whatever list of observables you give it and writes one $\alpha$ map
per observable:

```python
# Any set of spin-correlation observables to build up/down variations for
variables = [gen_c_kk, gen_c_rr, gen_c_nn, gen_ll_cHel]   # extend as needed
```

The sequence per observable:

**Step 1 — define the axes.** Write the three-axis config: $m_{t\bar t}$ (53 edges),
$\cos\theta^*$ (49 edges), and the observable as a 2-bin sign split. That sign axis
is the trick — it records only the sign of $O$ per $(m_{t\bar t},\cos\theta^*)$ cell.

**Step 2 — fill the 3-D gen histogram.** From the gen-level signal events, fill

$$
H[\,m_{t\bar t},\ \cos\theta^*,\ \text{sign}(O)\,],
$$

weighted by the generator weight. Each $(m_{t\bar t},\cos\theta^*)$ cell ends up with
the pair $(N_+, N_-)$.

**Step 3 — counts → $\alpha$.** Load the histogram, read $N_+$ (`...,1`) and $N_-$
(`...,0`), form the asymmetry, apply the observable prefactor (the two code blocks in
Section 3).

**Step 4 — collapse and store.** Sum away the sign axis, pack value + variance, and
write each map keyed by process and observable. The file is appended to if it exists,
created fresh otherwise:

The maps written here are exactly what `NoSC_reweight` (Section 4) loads at
application time. To add observables, you only extend the `variables` list.

---

## 9. Generator-level inputs

| Branch | Meaning |
|---|---|
| `gen_ttbar_mass` | $m_{t\bar t}$ (GeV), $t\bar t$ invariant mass |
| `gen_top_scatteringangle_ttbarframe` | $\cos\theta^*$, top scattering angle in the $t\bar t$ frame, $[-1,1]$ |
| `trueLevelWeight` | generator-level event weight |
| `gen_ll_cHel` | $\cos\varphi$, helicity opening angle |
| `gen_c_kk` … `gen_c_rr` | spin-correlation components $C_{ij}$ (9 entries) |
| `gen_b1k`, `gen_b1n`, `gen_b1r` / `gen_b2k`, `gen_b2n`, `gen_b2r` | polarization projections (full-density-matrix case) |

The production variables and the gen weight are always used. Beyond those, the
observables in this list are whichever spin-correlation components you want up/down
variations of — a single one (e.g. `cos φ`) for one measurement, or the full set of
$C_{ij}$ for the steerability.

---

### Note on gen-particle collection

The gen-particles used to construct the strength map and weight vector **need to be an inclusive gen-collection**.

The reason this is necessary is that the $(m_{t\bar t},\cos\theta^*)$ values for a given event can change drastically from gen-level to reco-level (i.e. event migration) due to hadronization and detector effects. This means if we are looking at a particular event with a given reco-level $(m_{t\bar t},\cos\theta^*)$ value which was not originally produced in that $(m_{t\bar t},\cos\theta^*)$ bin at gen-level, then it will require the weight from a different gen-level $(m_{t\bar t},\cos\theta^*)$ bin that would have been lost due to reco-level cuts.

Hence, the inclusive gen-level ttbar system of gen-particles needs to be saved and accessible and cannot be touched while the pre-selection and selection applies cuts on reco-events. One alternative I can think of is to simply create the gen-level strength map and weight vectors for every event and save those instead of the ttbar gen-particles themselves, this of course leads to less flexibility for potential future changes.
