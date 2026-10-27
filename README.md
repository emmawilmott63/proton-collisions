# proton-collisions: D0 → Kπ with Pythia8, ROOT, FastJet and TMVA

Generates pp collisions at √s = 200 GeV (Pythia8), builds K–π candidates, and
trains a BDT (ROOT TMVA) to separate true D0 → Kπ pairs from combinatorial
background. Optionally smears the K and π momenta to study detector effects.

## Layout

```
macros/      all C++ sources and ROOT macros (plotD0.C)
scripts/     run_generation.sh, run_smear_study.sh, run_ablation.sh
docs/        TMVA_REFERENCE.md: features, importance, all TMVA parameters, STAR preparation
data/        generated particle files and candidate files (created)
plots/       other plots, e.g. from mass and sanity       (created)
logs/        generation and candidate-building logs       (created)
TMVA-output/ everything the TMVA programs produce         (created)
  models/<dataset>/   trained model: weights/*.xml and train_config.txt
  training/           TMVA training files (open with TMVAGui)
  results/            response histograms, importance CSVs, tmva_summary.csv
  plots/              BDT response, significance, background-mass PNGs
  logs/               train and apply logs
```

## 1. Build

Edit `PYTHIA` and `FASTJET` at the top of the `Makefile` if your install paths
differ, make sure ROOT is set up in the shell, then:

```
make            # everything
make buildtree  # or a single program
```

Run everything below from the project root.

## 2. Generate data

```
./scripts/run_generation.sh inc 1000000 2     # sample, events per bin, parallel jobs
```

Generates five pTHat bins (5–10, 10–20, 20–40, 40–60, 60–100 GeV/c), then merges
them into `data/d0_inc_5_100.root` with `hadd`. Finished bins are skipped on a
rerun. Optional environment variables: `EDGES`, `EVENTS_PER_BIN`, `SEED_BASE`,
`FORCE=1` (see the script header). The 5–10 bin carries about 97% of the cross
section, so give it the most events.

Single run: `./jets-inc 100000 54321 data/d0_inclusive.root --pthat-min=5 --pthat-max=100`
(`./jets-charm` is the charm-enriched version).

## 3. Quick looks (optional)

```
./sanity data/d0_inc_5_100.root                  # jet pT, particle IDs
./mass   data/d0_inc_5_100.root --min-daughter-pt=0.6 --max-eta=1   # K-pi mass, OS - SS
./reconstruct-d0 data/d0_inc_5_100.root data/d0_reco_5_100.root
root -l -b -q 'macros/plotD0.C("data/d0_reco_5_100.root")'
```

`mass` and `reconstruct-d0` make mass plots only; TMVA needs the candidate tree
from step 4, not `d0_reco_*.root`.

## 4. Candidates → train → apply

```
./buildtree data/d0_inc_5_100.root                       # -> data/candidates_inc_5_100.root
./features  data/candidates_inc_5_100.root               # optional: signal vs background features
./trainTMVA data/candidates_inc_5_100.root --dataset=ds_inc
./applyTMVA data/candidates_inc_5_100.root --dataset=ds_inc
```

Outputs: `TMVA-output/plots/tmva_BDT_inc_5_100_{response,significance,bkgmass}.png`,
`TMVA-output/results/bdt_response_inc_5_100.root` and one row in `TMVA-output/results/tmva_summary.csv`.
`--dataset` must be the same for training and applying, and so must `--method`,
`--vars` and `--drop`.

### What the programs do

| Program | Job |
|---|---|
| `buildtree` | One row per opposite-sign K–π pair (daughter pT > 0.5, \|η\| < 1, 1.6 < M < 2.1). Writes the label `isSignal`, 27 features, and weights. |
| `trainTMVA` | BDT on even event numbers; odd event numbers are the test sample. `mass` is a spectator only, never a training variable. |
| `importanceTMVA` | Permutation importance of each variable and feature group for a trained BDT, plus each variable's correlation with the mass. |
| `applyTMVA` | Scores the test events, scans the BDT cut, and reports S, B, S/B and $S/\sqrt{S+B}$ (with Monte Carlo error) in a signal window. Prints a mass-sculpting check. |

### Weights

- `weight`: 1 for signal, 1/`bkg-keep` for background.
- `xsWeight`: σ_gen / N_generated of the pTHat bin (mb per event), read from `RunInfo`.
- Yields in `applyTMVA` = weight × xsWeight × L, with `--lumi-pb` (default 1 pb⁻¹).
- Training uses `weight` only by default; `--train-weights=xs` uses the physical mixture.

### Common options (written `--key=value`; `buildtree`, `trainTMVA` and `applyTMVA` reject unknown options)

| Program | Options |
|---|---|
| `buildtree` | `--smear-pt --smear-eta --smear-phi`, `--min-daughter-pt --max-eta --mass-min --mass-max`, `--bkg-keep`, `--tag`, `--out`, `--seed` |
| `trainTMVA` | BDT: `--method --boost-type --ntrees --depth --min-node --ncuts --ada-beta --shrinkage --bagging --bag-fraction --separation --prune-method --bdt-options`. Rows: `--presel --sig-cut --bkg-cut --signal-file --background-file`. Other: `--vars --drop --drop-group --train-weights --norm-mode --dataset --tag` (all described in `docs/TMVA_REFERENCE.md`) |
| `applyTMVA` | `--events=test\|train\|all --lumi-pb --sig-window --window-nsigma --min-bkg-rows --cut --presel --dataset --tag` (`--cut` is the BDT cut; `--presel` must match training) |

Output names are tagged from the input file name (`d0_inc_5_100.root` →
`_inc_5_100`); `--tag=name` overrides that.

## 5. TMVA diagnostics (no smearing needed)

```
make tmva                         # candidates -> train -> apply, no smearing (tag smear0)
make importance TAG=smear0        # permutation importance + mass correlation -> TMVA-output/results/importance_smear0.csv
make ablation   TAG=smear0        # retrain without each feature group; compare auc in TMVA-output/results/tmva_summary.csv
grep -A 32 "Ranking result" TMVA-output/logs/smear0_train.log      # TMVA's own variable ranking
```

Every training writes its model straight into `TMVA-output/models/<dataset>/`: `weights/TMVAClassification_<method>.weights.xml`
(read by `applyTMVA`) and `train_config.txt` (command, selections, TMVA options, rows, variables). Nothing is moved
afterwards: TMVA is run from inside `TMVA-output/models/`, so its output lands there. Training again with the same
`--dataset` overwrites that model. Results from before the `TMVA-output/` folder existed (`ds_*` folders, `weight/`,
TMVA files in `data/`, `plots/` and `logs/`) are moved there once with `make migrate-output`.
Feature definitions, importance interpretation, every BDT parameter, and the plan for STAR embedding
(two-file training with `--signal-file` / `--background-file` / `--bkg-cut`, and a sideband closure
test you can run on simulation) are in `docs/TMVA_REFERENCE.md`.

## 6. Smearing study (automated)

```
make study                                   # levels 0, 1, 3, 5, 10% relative pT smearing
make study LEVELS="0 0.02 0.05"              # other levels
make study TRAIN_OPTS="--depth=5 --ntrees=600 --min-node=0.5"
make study METHOD=BDTG
make reapply                                 # rerun only applyTMVA (keeps trained BDTs)
make clean-study                             # delete this study's outputs
make help
```

For each level, `make` runs build, train and apply, and prints `TMVA-output/results/tmva_summary.csv`
at the end. It only redoes steps whose inputs changed, so you can interrupt and
resume. Rebuilding a program does not redo finished results: run
`make clean-study` after changing the code. Write levels consistently
(`0.10` and `0.1` produce different tags). `./scripts/run_smear_study.sh` is a
wrapper around `make study`.

Useful extras: `BUILD_OPTS="--smear-eta=0.002 --smear-phi=0.002"`,
`APPLY_OPTS="--lumi-pb=10"`.

### Reading the plots

- **response**: BDT score for true D0 (blue) vs fake pairs (red), each normalized to unit area.
- **significance**: $S/\sqrt{S+B}$ vs BDT cut. The flat left part is "no cut"; the peak is the best cut; the grey band is the Monte Carlo statistical error.
- **bkgmass**: fake-pair mass before and after the cut. A bump at 1.865 GeV/c² would mean the BDT learned the mass.

### Summary CSV columns

`tag, method, events, lumi_pb, mass_peak, window, sigma_eff, S_nocut, B_nocut, signif_nocut, best_cut, S_best, B_best, signif_best, signif_best_err, effS_best, bkgrej_best, auc`.
`auc` (flat weights) measures separation power independently of the mass window.
The signal window is set automatically from the signal mass distribution
(median ± 2 robust σ, minimum ±0.005 GeV/c²), so it follows the smearing.

Overtraining check: `./applyTMVA <candidates> --dataset=ds_X --events=train --tag=X_train`
and compare its `auc` with the test row.

## Known limitations

- K and π are identified by true Pythia `pid` (perfect PID), so the background is lower than in a real detector.
- Only the K and π momenta are smeared. Jets, isolation and event-level features use unsmeared particles.
- Absolute yields depend on `--lumi-pb` and the window; compare ratios and trends, not single numbers. Differences of a few percent in $S/\sqrt{S+B}$ are Monte Carlo noise.
- `isoRel*` divides by the pair pT and has very long tails.

## Troubleshooting

| Problem | Fix |
|---|---|
| `permission denied: ./scripts/...sh` | `chmod +x scripts/*.sh` |
| `Unknown option --...` | Check spelling; options must be `--key=value`. |
| `Cannot find TMVA-output/models/<dataset>/weights/...xml` | Run `trainTMVA` first, with the same `--dataset` and `--method`. Results from before `TMVA-output/` existed: `make migrate-output`. |
| `No xsWeight branch` | Rebuild the candidates with the current `buildtree`. |
| `Too few signal or background rows` | Generate more events. Training needs at least 100 of each class in both halves. |
| Linker error about TMVA | The Makefile links `-lTMVA` for `trainTMVA` and `applyTMVA`. Make sure ROOT was built with TMVA. |
