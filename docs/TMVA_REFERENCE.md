# TMVA reference: features, importance, parameters, and the road to STAR embedding

## 1. What one row is (the candidate)

One row = one opposite-sign K–π pair from **one event**, with the K and π
assigned by their particle IDs, passing the preselection (daughter pT > 0.5
GeV/c, |η| < 1, 1.6 < M < 2.1 GeV/c²). The BDT scores each row independently.

- **Label** `isSignal = 1`: both tracks are the daughters of the same generated D0 → Kπ.
- **Hard conditions** (opposite sign, same event, acceptance, mass window) are part of
  how rows are built (`buildtree`) or of the row selection (`--presel`), not of the BDT.
  See section 5.

### Input contract (what any candidate maker must write)

For simulation or real STAR data, `trainTMVA` and `applyTMVA` only need a TTree
`Candidates` with these float/int branches. Anything that writes them can feed the pipeline.

| Branch | Type | Meaning |
|---|---|---|
| the feature branches of section 2 | float | one per training variable |
| `mass` | float | K–π invariant mass; spectator only, never trained on |
| `event` | int | unique per event within a file; parity (even/odd) decides train/test |
| `weight` | float | per-row weight (1 for data; MC reweighting factors for embedding) |
| `isSignal` | int | 1 for signal rows; 0 otherwise (may be absent in a background-only file) |
| `xsWeight` | double | optional: cross section per generated event [mb], for the pThat-binned simulation |

## 2. Features

27 are used by default (`tmvaCommon.h`). Group = the unit used by `--drop-group`
and by group importance.

| Feature | Definition | Group | Default | Same definition possible in real data? |
|---|---|---|---|---|
| `ptK`, `ptPi` | transverse momentum of K, π | daughters | yes | yes (tracks) |
| `etaK`, `etaPi` | pseudorapidity of K, π | daughters | yes | yes |
| `ptHard`, `ptSoft` | harder / softer daughter pT | daughters | no (redundant) | yes |
| `ptSum` | scalar pT sum of the two | daughters | no (redundant) | yes |
| `ptBalance` | \|pT_K − pT_π\| / (pT_K + pT_π) | daughters | yes | yes |
| `dEta`, `dPhi` | η_K − η_π, φ_K − φ_π (signed, wrapped) | pair | yes | yes |
| `dR` | √(dEta² + dPhi²) | pair | yes | yes |
| `openAngle` | 3D opening angle between K and π momenta | pair | yes | yes |
| `pairPt`, `pairEta` | pT, η of the K+π four-vector | pair | yes | yes |
| `cosThetaStar` | cos of the K angle in the pair rest frame, relative to the pair direction (K mass hypothesis) | pair | yes | yes |
| `sameJet` | 1 if K and π are constituents of the same anti-kT (R = 0.4, jet pT > 2 GeV/c) jet | jet | yes | only if jets are built the same way |
| `jetPt`, `jetEta`, `jetMass`, `jetNConst` | reference jet: the shared jet if `sameJet`, else the jet nearest the pair | jet | all but `jetEta` | as above |
| `zJet` | pairPt / jetPt | jet | yes | as above |
| `dRJet` | ΔR from pair direction to the reference-jet axis (9 if no jet) | jet | yes | as above |
| `isoRel02`, `isoRel03`, `isoRel`, `isoRel06` | charged pT in a cone R = 0.2 / 0.3 / 0.4 / 0.6 around the pair axis, divided by pairPt (K, π excluded) | isolation | yes | tracks only, with real track cuts |
| `isoNeutralRel` | same for neutral particles, R = 0.4 | isolation | yes | needs calorimeter towers |
| `nCharged` | charged multiplicity, \|η\| < 1, pT > 0.2 | event | yes | yes, but depends on track cuts and efficiency |
| `nJets` | number of jets in \|η\| < 1 | event | yes | as for jets |
| `sphericityT`, `thrustT` | transverse sphericity / thrust of those charged particles and jets | event | yes | as above |

**Simulation-only shortcuts to be aware of** (they make the simulated
problem easier or different from data):

- K and π are identified by true particle ID (perfect PID).
- Jets, isolation and event-level features use every final-state particle with true momenta, with no detector acceptance or efficiency.
- There are no topological variables (decay length, DCA of the daughters, pointing angle). In standard D0 analyses these are usually the strongest discriminators; they are missing here only because the generator output stores momenta, not vertices.
- `isoRel*` divides by `pairPt`, which can be tiny, giving values up to about 1400.

## 3. Importance

Three complementary measures. Run them on the same trained BDT.

| Measure | How | What it tells you |
|---|---|---|
| TMVA's own ranking | `grep -A 32 "Ranking result" logs/<tag>_train.log` | how often and how usefully each variable is used in the trees (gain-based; favours variables with many distinct values) |
| Permutation importance | `make importance TAG=smear0` (or `./importanceTMVA <candidates> --dataset=ds_<tag>`) | AUC drop when a variable, or a whole group, is shuffled; also prints each variable's correlation with the mass in background rows. Writes `data/importance<tag>.csv` |
| Ablation | `make ablation TAG=smear0` (or `./scripts/run_ablation.sh <candidates>`) | retrains without each group and compares AUC and best significance in `data/tmva_summary.csv`. The most reliable, because the BDT cannot fall back on a correlated partner |

How to read them:

- Single-variable numbers underestimate importance when variables overlap (for example `ptK`, `ptPi`, `pairPt`). Use the group numbers or the ablation for decisions.
- A variable that is **both important and strongly mass-correlated in the background** (`corr_mass_bkg` large) is the one that can sculpt the background mass shape.
- The AUC differences between configurations are small numbers; trust a difference only if it exceeds the `err` column or the scatter between reruns.

Importance results from your trained BDT go here once you have run the commands above:

| Feature or group | TMVA ranking | Permutation drop | Ablation ΔAUC |
|---|---|---|---|
| (fill in) | | | |

## 4. TMVA parameters

All are command-line options of `trainTMVA` (defaults in brackets); the exact string passed to TMVA is saved in `weight/<dataset>/train_config.txt`.

### BDT

| Option | TMVA name | Default | Effect | Tuning notes |
|---|---|---|---|---|
| `--ntrees` | NTrees | 400 | number of boosted trees | more trees: lower bias; overtrains eventually, so watch train vs test AUC |
| `--depth` | MaxDepth | 3 | levels per tree, so the order of feature interactions | 3–5 is typical; your train/test matched at 3, so try 4–5 |
| `--min-node` | MinNodeSize | 2.0 % | smallest allowed leaf, as a percent of training events | regularizer; lower it (0.5) when going deeper |
| `--ncuts` | nCuts | 20 | grid points per variable when searching splits | 20–50; higher is slower and slightly finer |
| `--boost-type` | BoostType | AdaBoost | AdaBoost, RealAdaBoost, Grad, Bagging | `--method=BDTG` selects Grad |
| `--ada-beta` | AdaBoostBeta | 0.5 | AdaBoost learning rate | smaller is slower but more robust against noise |
| `--shrinkage` | Shrinkage | 0.1 | Grad learning rate | smaller needs more trees |
| `--bagging` | UseBaggedBoost | 1 | each tree sees a random subsample | reduces overtraining |
| `--bag-fraction` | BaggedSampleFraction | 0.5 | fraction used per tree | |
| `--separation` | SeparationType | GiniIndex | split criterion (CrossEntropy, MisClassificationError, SDivSqrtSPlusB) | rarely matters; ignored by Grad |
| `--prune-method` | PruneMethod | NoPruning | ExpectedError, CostComplexity | with shallow trees pruning is unnecessary |
| `--prune-strength` | PruneStrength | 0 | pruning strength | only with pruning |
| `--bdt-options` | any | | raw string appended last | for anything not listed |

### Samples, weights, variables

| Option | Default | Effect |
|---|---|---|
| train / test split | even / odd `event` | split by event, so event-level features cannot leak; fixed in the code |
| `--norm-mode` | EqualNumEvents | training sees equal total weight for signal and background; test yields are unaffected |
| `--train-weights` | flat | `weight` only; `xs` multiplies by `xsWeight` (physical pThat mixture, noisier) |
| `--presel`, `--sig-cut`, `--bkg-cut` | none | row selections (section 5) |
| `--signal-file`, `--background-file` | none | two-file mode (section 6) |
| `--vars`, `--drop`, `--drop-group` | 27 features | choose the input variables; must match in `applyTMVA` |
| `--dataset`, `--method`, `--tag`, `--out` | | names of the weights directory, method title, and output files |

### Fixed choices

- Classification with two classes; no variable transformations (BDT splits are unchanged by monotonic rescaling of one variable).
- The test sample is never used for training or tuning in `trainTMVA`; the best-cut choice in `applyTMVA` is made on the test sample, so quote it with its Monte Carlo error.

## 5. Can conditions like "same event" and "opposite sign" be hard-coded?

The BDT is a function of one row's feature vector, so it cannot enforce logical or structural
conditions by itself. There are four ways to build such conditions in, and the first two are the
right tools for yours:

1. **Define the row** (candidate builder). "Same event" is automatic because pairs are only formed
   inside one event; "opposite sign", the K/π assignment, acceptance and mass window are enforced
   in `buildtree`. A STAR candidate maker would enforce the same. Mixed-event or like-sign pairs
   are *different background constructions*, so they would be separate candidate sets, not a
   switch inside TMVA.
2. **Select rows** (hard cuts at training and application): `--presel` (both classes),
   `--sig-cut`, `--bkg-cut`. TMVA then trains and tests only on rows passing them, and
   `applyTMVA --presel` scores only those rows, so the cut is built into the whole chain.
3. **Give the BDT the physics as a feature**: `cosThetaStar`, `dR`, `zJet` and so on. For D0 the
   missing and most important ones are the topological variables (decay length, DCA, pointing angle).
4. **Constrain the model itself** (monotonicity, forced structure, embedded conservation laws):
   as far as I know, TMVA's BDT does not offer this. Other libraries do (for example XGBoost
   through TMVA's PyMVA wrapper).

Put the charge sign in only if you want the BDT to use it. Training a BDT on like-sign background
in the same sample does not require it as a feature.

## 6. Preparing for STAR embedding data

A typical D0 setup: **signal from embedding, background from real data**.

```
./trainTMVA --signal-file=embedding_candidates.root \
            --background-file=data_candidates.root \
            --bkg-cut="mass<1.80||mass>1.93" \
            --dataset=ds_star --tag=star
```

- Signal rows are `isSignal == 1` in the embedding file; background rows are every data candidate outside the signal mass window. If the data file has no `isSignal` branch it is not required.
- Weights: put the embedding reweighting (to the true D0 pT spectrum) and any data weights in `weight`.

**Closure test you can run now on simulation**: train the background from the sidebands only,
then check it still works on the true background under the peak.

```
./trainTMVA --signal-file=data/candidates_inc_5_100.root --background-file=data/candidates_inc_5_100.root \
            --bkg-cut="mass<1.80||mass>1.93" --dataset=ds_sideband --tag=sideband
./applyTMVA data/candidates_inc_5_100.root --dataset=ds_sideband --tag=sideband
```

Compare its AUC and best significance with the full-mass training. A clear loss means the
sidebands do not represent the background under the peak (typically because features depend on
mass), and the same would happen with real data.

### Checklist before the first real-data training

| Item | Status in this code |
|---|---|
| Candidate format and branch names | defined (section 1) |
| Two-file training with a sideband selection | implemented (untested on real data) |
| Same row selection in training and application | `--presel` |
| Importance, ablation, mass-correlation and sculpting checks | implemented |
| Reproducible training record | `weight/<dataset>/train_config.txt` |
| Features computable identically in embedding and data | open: jet, isolation, event-level and neutral features need detector-level definitions; track kinematics are fine |
| Topological variables | open: need vertex information from the generator (regeneration) or come directly from STAR |
| Data vs embedding agreement for each input variable (sideband-subtracted data against embedding signal) | open |
| Realistic PID instead of truth IDs | open |
| Cross-validation and a hyperparameter scan | open (the test half is currently used once) |
