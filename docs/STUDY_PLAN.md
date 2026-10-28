# Study plan: from TMVA tuning to report figures

Report structure: **Part I** TMVA (how it works, validation, parameter explorations),
**Part II** application to data and comparison with other methods.

## 0. Ground rules (what makes the results defensible)

1. **One figure of merit, defined up front** (section 1), the same in every study.
2. **Every number has an uncertainty.** Use K-fold cross-validation (`scripts/run_cv.sh`);
   differences smaller than the fold-to-fold spread are not differences.
3. **Tune on folds, report on the holdout, once.** All tuning runs use `--holdout=1`, which
   sets aside the events with `event % 7 == 6`. The final configuration is evaluated on those
   events a single time (`applyTMVA --events=holdout`). Re-tuning after looking at the holdout
   invalidates it.
4. **Change one thing at a time** (`scripts/run_scan.sh`), everything else fixed in `TRAIN_OPTS`.
5. **Every run records itself**: `TMVA-output/models/<dataset>/train_config.txt` and one row in
   `TMVA-output/results/tmva_summary.csv`. Keep the code under git and note the commit.
6. **Report negative results** (a parameter that does not matter is a finding).

## 1. Figures of merit

All are computed on test events and written to `tmva_summary.csv`.

- **ROC integral** (AUC): probability that a random signal candidate scores above a random
  background one. Independent of the mass window and of the luminosity.
- **Background rejection at fixed signal efficiency**, $1-\varepsilon_B$ at $\varepsilon_S = 30, 50, 70\,\%$:
  what an analyst sees when choosing a working point.
- **Significance in the signal window**
  $$\frac{S}{\sqrt{S+B}},\qquad S,B = \text{expected counts in } M_{K\pi}\in M_0\pm 2\sigma_M \text{ at } L$$
  with its Monte Carlo error. The best-cut value is optimistic (the cut is chosen on the same
  events); quote it only with that caveat or use the holdout.

For the report text: AUC and rejection are the headline numbers for method comparisons;
$S/\sqrt{S+B}$ links to the physics goal but depends on the window, the luminosity and the
smearing.

## 2. Study map

Commands run from the project root; build with `make` first. `python3 scripts/make_report.py`
turns all results into `report/figures/*.png|pdf` and `report/tables/*.tex|md`.

### Part I: the method

| ID | Question | Run | Produces | Notes |
|---|---|---|---|---|
| I.1 | What does the classifier see? | `./features data/candidates_inc_5_100.root` | signal vs background feature distributions | explains each feature (`docs/TMVA_REFERENCE.md`) |
| V1 | Baseline performance and its uncertainty | `./scripts/run_cv.sh cv_baseline` | `cv_summary` table, BDT response, significance, mass-shape plots | 5 folds, mean $\pm$ std |
| V2 | Overtraining? | `./applyTMVA <cand> --dataset=ds_cv_baseline_f0 --nfolds=5 --fold=0 --holdout=1 --events=train --tag=cv_baseline_f0_train` | train vs test AUC | gap well below the fold spread = no overtraining |
| V3 | Does the BDT sculpt the background mass shape? | `make importance TAG=<tag>` and the `bkgmass` plot | `corr_mass_bkg` per variable, background fraction in window before/after | required before any sideband method |
| V4 | Is there enough training data? | `./scripts/run_scan.sh bkg-subsample "0.05 0.1 0.2 0.5 1"` | learning curve (AUC vs training size) | flat curve at 1 = more data will not help |
| V5 | Which features matter? | `make importance TAG=<tag>`, `make ablation TAG=<tag>` | importance and ablation figures | read groups, not single variables |
| V6 | Do the training weights matter? | `./scripts/run_scan.sh train-weights "flat xs"` | scan table | `xs` = physical pThat mixture |
| V7 | Robustness to track smearing | `make study` (or `LEVELS="0 0.01 0.03 0.05 0.10"`) | `smearing` figure and table | the BDT is unaffected; the mass window is |
| V8 | Do sidebands represent the background? | closure test in `docs/TMVA_REFERENCE.md` section 6 | AUC and significance vs full-mass training | rehearsal for real data |

### Part I: parameter explorations

Run each as a scan, with cross-validation. Suggested ranges (the default is in bold):

| Parameter | Command | Values | Question |
|---|---|---|---|
| depth | `run_scan.sh depth "2 3 4 5 6 8"` | **3** | interaction order the problem needs |
| trees | `run_scan.sh ntrees "100 200 400 800 1600"` | **400** | convergence; pair with learning rate |
| leaf size | `run_scan.sh min-node "0.25 0.5 1 2 4"` | **2** | regularization |
| split grid | `run_scan.sh ncuts "10 20 40 80"` | **20** | resolution of the cut search |
| learning rate | `run_scan.sh ada-beta "0.2 0.5 1"` | **0.5** | AdaBoost step size |
| boost type | `run_scan.sh boost-type "AdaBoost RealAdaBoost Grad"` | **AdaBoost** | algorithm |
| Grad shrinkage | `TRAIN_OPTS="--boost-type=Grad" run_scan.sh shrinkage "0.03 0.1 0.3"` | **0.1** | (Grad only) |
| bagging | `run_scan.sh bag-fraction "0.3 0.5 0.8"` | **0.5** | stochastic boosting |
| separation | `run_scan.sh separation "GiniIndex CrossEntropy"` | **GiniIndex** | split criterion |

Practical:

- Time one training first (`time ./trainTMVA ...`), then budget (values $\times$ folds) trainings.
- For a first pass use `FOLDS=3 TRAIN_OPTS="--bkg-subsample=0.2"`, then repeat the
  interesting region at full statistics.
- Depth and trees interact; after the one-at-a-time scans, run a small 2-D grid (depth $\times$
  trees) around the best region.

**Choosing the final configuration:** take the *simplest* setting whose mean is within one
fold-to-fold standard deviation of the best (fewer trees, shallower, larger leaves). Then train it
on all non-holdout events and evaluate on the holdout once.

### Part II: data and comparison

| ID | Question | How | Status |
|---|---|---|---|
| A1 | Train with embedding signal and data-sideband background | `trainTMVA --signal-file=... --background-file=... --bkg-cut="mass<1.80\|\|mass>1.93"` | implemented; closure-tested on simulation (V8) |
| A2 | Do embedding and data agree for every input variable? | compare sideband-subtracted data with embedding signal, variable by variable | to be written |
| A3 | Yield with and without the BDT | fit the mass spectrum before and after the BDT cut; compare signal yield, statistical error and S/B | needs the data candidate maker and a fit |
| A4 | Comparison with other methods | `TRAIN_OPTS="--also=Fisher,LD,Likelihood,MLP,BDTG" ./scripts/run_cv.sh cv_methods` | implemented; plus a cut-based baseline (to be written) |
| A5 | Systematic uncertainty from the BDT | vary the cut, the training sample, the hyperparameters; take the spread of the corrected yield | after A3 |

For A4, compare all methods on the same folds at the same figure of merit. Include a
rectangular-cut baseline tuned on the same variables, since that is what the BDT replaces.
Comparisons with published results should quote the published selection, acceptance and
luminosity; the numbers here are generator-level and will not match directly.

## 3. Suggested order

1. V1 baseline (also gives the timing). 2. V2, V3. 3. V5 (importance and ablation).
4. Parameter scans with `--bkg-subsample`, then refine. 5. V4, V6. 6. Final configuration
and a single holdout evaluation. 7. V7 and V8, then Part II.

## 4. What I would add next (in order of value)

1. **Signal efficiency and background rejection versus $p_T$** at the chosen cut: the standard
   plot, and it exposes pT sculpting.
2. **Paired differences in scans** (each configuration minus the baseline on the same folds),
   which are much tighter than comparing independent means.
3. **Rectangular-cut baseline** for A4.
4. **Topological variables** (decay length, DCA, pointing angle) from generator vertices. This
   will change the physics result more than any hyperparameter, so it should come before the
   final scans if it is going to happen at all.
5. **A random or Bayesian hyperparameter search** if the one-at-a-time scans show strong
   interactions.

## 5. Caveats to state in the report

- Truth particle IDs (perfect PID) and no detector response: the background is smaller than in a
  real detector. The simulation establishes *relative* behaviour (which parameters matter, how the
  BDT compares with cuts), not absolute performance.
- Jets, isolation and event-level features use every final-state particle with true momenta.
- Yields scale with the luminosity and depend on the signal window; quote ratios and trends.
- The 5-10 GeV/$c$ pThat bin carries about 97% of the physical weight, so physical-mixture
  numbers rest on few weighted rows; the Monte Carlo error is printed with every significance.
