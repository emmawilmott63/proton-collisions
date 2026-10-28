#!/usr/bin/env python3
"""Figures and tables for the report, from the TMVA result CSVs.

Reads   TMVA-output/results/tmva_summary.csv     (one row per applyTMVA run)
        TMVA-output/results/importance*.csv      (importanceTMVA)
Writes  report/figures/*.png and *.pdf
        report/tables/*.tex (booktabs) and *.md

Which study each tag belongs to is decided by its name:

  <name>_f<k>               one fold of a cross-validation run
  scan_<param>_<value>_f<k> one fold of a parameter scan point
  smear<x>                  smearing study (smear0p03 = 3 %)
  ablation_<base>_all / _no_<group>   feature-group ablation
  any row with events == holdout      final holdout evaluation

Usage (from the project root):
  python3 scripts/make_report.py
  python3 scripts/make_report.py --results TMVA-output/results --out report
Only matplotlib, pandas and numpy are needed.
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

METRICS = {
    "auc": ("ROC integral", 4),
    "rej_effS50": (r"bkg. rejection at $\varepsilon_S = 50\%$", 4),
    "signif_best": (r"best $S/\sqrt{S+B}$", 1),
}
FOLD_RE = re.compile(r"^(.*)_f(\d+)$")
SCAN_RE = re.compile(r"^scan_(.+)_([^_]+)$")
SMEAR_RE = re.compile(r"^smear(\d+)(?:p(\d+))?$")

plt.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 100, "savefig.bbox": "tight",
})


# ------------------------------------------------------------------ helpers
def save_fig(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(out / "figures" / f"{name}.{ext}", dpi=200)
    plt.close(fig)
    print(f"  figure  {name}")


def pm(mean, std, nd):
    if pd.isna(mean):
        return ""
    if pd.isna(std) or std == 0:
        return f"{mean:.{nd}f}"
    return f"{mean:.{nd}f} $\\pm$ {std:.{nd}f}"


def pm_md(text):
    return text.replace("$\\pm$", "±")


def write_table(rows, columns, out, name, caption=""):
    """rows: list of dicts; columns: list of (key, header). LaTeX + Markdown."""
    if not rows:
        return
    heads = [h for _, h in columns]
    body = [[str(r.get(k, "")) for k, _ in columns] for r in rows]

    def esc(s):
        return s.replace("_", r"\_").replace("%", r"\%")

    tex = ["\\begin{table}[ht]", "\\centering",
           "\\begin{tabular}{l" + "r" * (len(columns) - 1) + "}", "\\toprule",
           " & ".join(esc(h) if "$" not in h else h for h in heads) + r" \\",
           "\\midrule"]
    for row in body:
        cells = [esc(c) if "$" not in c else c.replace("_", r"\_") for c in row]
        tex.append(" & ".join(cells) + r" \\")
    tex += ["\\bottomrule", "\\end{tabular}"]
    if caption:
        tex.append(f"\\caption{{{caption}}}")
    tex.append("\\end{table}")
    (out / "tables" / f"{name}.tex").write_text("\n".join(tex) + "\n")

    md = ["| " + " | ".join(heads) + " |",
          "|" + "|".join("---" for _ in heads) + "|"]
    md += ["| " + " | ".join(pm_md(c) for c in row) + " |" for row in body]
    (out / "tables" / f"{name}.md").write_text("\n".join(md) + "\n")
    print(f"  table   {name}")


def parse_tags(df):
    base, fold = [], []
    for t in df["tag"].astype(str):
        m = FOLD_RE.match(t)
        base.append(m.group(1) if m else t)
        fold.append(int(m.group(2)) if m else -1)
    df = df.copy()
    df["base"], df["foldid"] = base, fold
    return df


def agg(group, metric):
    v = group[metric].dropna().to_numpy(dtype=float)
    if len(v) == 0:
        return np.nan, np.nan, 0
    return v.mean(), (v.std(ddof=1) if len(v) > 1 else 0.0), len(v)


# ----------------------------------------------------------------- studies
def cv_summary(df, out):
    """Mean +- std over folds for every cross-validated configuration."""
    d = df[(df["events"] == "test") & (df["foldid"] >= 0)]
    rows = []
    for (base, method), g in d.groupby(["base", "method"], sort=True):
        row = {"config": base, "method": method, "folds": len(g)}
        for m, (_, nd) in METRICS.items():
            mean, std, _ = agg(g, m)
            row[m] = pm(mean, std, nd)
        rows.append(row)
    cols = [("config", "Configuration"), ("method", "Method"),
            ("folds", "Folds")] + [(m, lab) for m, (lab, _) in METRICS.items()]
    write_table(rows, cols, out, "cv_summary",
                "Cross-validation results (mean $\\pm$ standard deviation over folds).")


def scans(df, out):
    d = df[(df["events"] == "test") & (df["foldid"] >= 0)]
    params = {}
    for base in d["base"].unique():
        m = SCAN_RE.match(base)
        if m:
            params.setdefault(m.group(1), []).append((m.group(2), base))
    for param, items in sorted(params.items()):
        def num(v):
            try:
                return float(v)
            except ValueError:
                return None
        numeric = all(num(v) is not None for v, _ in items)
        items.sort(key=lambda it: num(it[0]) if numeric else it[0])
        fig, axes = plt.subplots(len(METRICS), 1, sharex=True,
                                 figsize=(6.2, 2.7 * len(METRICS)))
        rows = []
        for ax, (metric, (label, nd)) in zip(np.atleast_1d(axes), METRICS.items()):
            xs, means, stds = [], [], []
            for k, (val, base) in enumerate(items):
                g = d[(d["base"] == base) & (d["method"] == "BDT")]
                mean, std, n = agg(g, metric)
                x = num(val) if numeric else k
                xs.append(x); means.append(mean); stds.append(std)
                ax.plot([x] * n, g[metric].dropna(), ".", color="0.6", ms=4)
                if metric == "auc":
                    rows.append({"value": val, "folds": n})
                rows[[r["value"] for r in rows].index(val)][metric] = pm(mean, std, nd)
            ax.errorbar(xs, means, yerr=stds, fmt="o-", color="C0", capsize=3)
            ax.set_ylabel(label)
        ax = np.atleast_1d(axes)[-1]
        if numeric:
            ax.set_xlabel(param)
            pos = [x for x in xs if x > 0]
            if pos and len(pos) == len(xs) and max(pos) / min(pos) > 20:
                ax.set_xscale("log")
        else:
            ax.set_xticks(range(len(items)))
            ax.set_xticklabels([v for v, _ in items], rotation=20)
            ax.set_xlabel(param)
        np.atleast_1d(axes)[0].set_title(f"Scan of {param} (points: folds, bars: mean $\\pm$ std)")
        save_fig(fig, out, f"scan_{param}")
        write_table(rows, [("value", param), ("folds", "Folds")] +
                    [(m, lab) for m, (lab, _) in METRICS.items()],
                    out, f"scan_{param}", f"Scan of {param}.")


def smearing(df, out):
    d = df[(df["events"] == "test") & (df["foldid"] < 0) & (df["method"] == "BDT")]
    pts = []
    for _, r in d.iterrows():
        m = SMEAR_RE.match(str(r["tag"]))
        if m:
            x = float(f"{m.group(1)}.{m.group(2)}") if m.group(2) else float(m.group(1))
            pts.append((x, r))
    if len(pts) < 2:
        return
    pts.sort(key=lambda p: p[0])
    x = np.array([p[0] for p in pts]) * 100
    col = lambda c: np.array([p[1][c] for p in pts], dtype=float)  # noqa: E731
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    axes[0].plot(x, col("sigma_eff") * 1000, "o-")
    axes[0].set_ylabel(r"signal mass width $\sigma$ [MeV/$c^2$]")
    axes[1].errorbar(x, col("signif_nocut"), fmt="s--", label="no BDT cut")
    axes[1].errorbar(x, col("signif_best"), yerr=col("signif_best_err"),
                     fmt="o-", capsize=3, label="best BDT cut")
    axes[1].set_ylabel(r"$S/\sqrt{S+B}$"); axes[1].set_yscale("log"); axes[1].legend()
    axes[2].plot(x, col("auc"), "o-"); axes[2].set_ylabel("ROC integral")
    for ax in axes:
        ax.set_xlabel(r"relative $p_T$ smearing [%]")
    fig.suptitle("Effect of track-momentum smearing")
    save_fig(fig, out, "smearing")
    rows = [{"smear": f"{p[0]*100:g}", "sigma": f"{p[1]['sigma_eff']*1000:.1f}",
             "nocut": f"{p[1]['signif_nocut']:.1f}",
             "best": pm(p[1]["signif_best"], p[1]["signif_best_err"], 1),
             "gain": f"{p[1]['signif_best'] / p[1]['signif_nocut']:.2f}",
             "auc": f"{p[1]['auc']:.4f}"} for p in pts]
    write_table(rows, [("smear", r"$p_T$ smearing [\%]"),
                       ("sigma", r"$\sigma_M$ [MeV/$c^2$]"),
                       ("nocut", r"$S/\sqrt{S+B}$ no cut"),
                       ("best", r"$S/\sqrt{S+B}$ best cut"),
                       ("gain", "Gain"), ("auc", "ROC integral")],
                out, "smearing", "Smearing study.")


def ablation(df, out):
    d = df[(df["events"] == "test") & (df["method"] == "BDT")
           & df["tag"].astype(str).str.startswith("ablation_")]
    for base in sorted({re.sub(r"_(all|no_[a-z]+)$", "", t) for t in d["tag"]}):
        ref = d[d["tag"] == f"{base}_all"]
        if ref.empty:
            continue
        ref = ref.iloc[0]
        rows, labels, dauc, dsig = [], [], [], []
        for _, r in d[d["tag"].str.startswith(base + "_no_")].iterrows():
            g = r["tag"].split("_no_")[-1]
            labels.append(g)
            dauc.append(r["auc"] - ref["auc"])
            dsig.append(r["signif_best"] / ref["signif_best"] - 1)
            rows.append({"group": g, "auc": f"{r['auc']:.4f}",
                         "dauc": f"{dauc[-1]:+.4f}",
                         "dsig": f"{100 * dsig[-1]:+.1f}\\%"})
        if not rows:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
        axes[0].barh(labels, dauc, color="C3"); axes[0].set_xlabel(r"$\Delta$ ROC integral")
        axes[1].barh(labels, np.array(dsig) * 100, color="C3")
        axes[1].set_xlabel(r"$\Delta$ best $S/\sqrt{S+B}$ [%]")
        fig.suptitle("Retraining without one feature group (vs. all features)")
        name = "ablation_" + base[len("ablation_"):]
        save_fig(fig, out, name)
        write_table(rows, [("group", "Group removed"), ("auc", "ROC integral"),
                           ("dauc", r"$\Delta$ ROC"), ("dsig", r"$\Delta S/\sqrt{S+B}$")],
                    out, name,
                    f"Ablation (reference: all features, ROC integral {ref['auc']:.4f}).")


def importance(results, out):
    for path in sorted(results.glob("importance*.csv")):
        d = pd.read_csv(path)
        groups = d[d["type"] == "group"].sort_values("auc_drop")
        var = d[d["type"] == "variable"].sort_values("auc_drop").tail(15)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        axes[0].barh(groups["name"], groups["auc_drop"], xerr=groups["err"], color="C0")
        axes[0].set_title("Feature groups"); axes[0].set_xlabel("AUC drop when shuffled")
        axes[1].barh(var["name"], var["auc_drop"], xerr=var["err"], color="C1")
        axes[1].set_title("Top variables"); axes[1].set_xlabel("AUC drop when shuffled")
        base = d["baseline_auc"].iloc[0] if "baseline_auc" in d else float("nan")
        fig.suptitle(f"Permutation importance ({path.stem}, baseline AUC {base:.4f})")
        save_fig(fig, out, path.stem)
        rows = [{"name": r["name"], "group": r["group"],
                 "drop": pm(r["auc_drop"], r["err"], 4),
                 "corr": "" if pd.isna(r["corr_mass_bkg"]) else f"{r['corr_mass_bkg']:+.2f}"}
                for _, r in d.sort_values("auc_drop", ascending=False).iterrows()]
        write_table(rows, [("name", "Feature / group"), ("group", "Group"),
                           ("drop", "AUC drop"), ("corr", r"corr. with $M$ (bkg)")],
                    out, path.stem, "Permutation importance.")


def holdout(df, out):
    d = df[df["events"] == "holdout"]
    rows = []
    for _, r in d.iterrows():
        rows.append({"tag": r["tag"], "method": r["method"],
                     "auc": f"{r['auc']:.4f}",
                     "rej": f"{r['rej_effS50']:.4f}",
                     "sig": pm(r["signif_best"], r["signif_best_err"], 1)})
    write_table(rows, [("tag", "Configuration"), ("method", "Method"),
                       ("auc", "ROC integral"),
                       ("rej", r"bkg. rejection at $\varepsilon_S = 50\%$"),
                       ("sig", r"best $S/\sqrt{S+B}$")],
                out, "holdout",
                "Final evaluation on the holdout events (never used for tuning).")


def method_comparison(df, out):
    d = df[(df["events"] == "test") & (df["foldid"] >= 0)]
    for base, g in d.groupby("base"):
        methods = sorted(g["method"].unique())
        if len(methods) < 2:
            continue
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
        for ax, metric in zip(axes, ["auc", "rej_effS50"]):
            means, stds = [], []
            for m in methods:
                mean, std, _ = agg(g[g["method"] == m], metric)
                means.append(mean); stds.append(std)
            ax.bar(methods, means, yerr=stds, capsize=3, color="C2")
            ax.set_ylabel(METRICS[metric][0])
            lo = np.nanmin(means) - 3 * np.nanmax(stds + [0.005])
            ax.set_ylim(max(0, lo), 1.0)
        fig.suptitle(f"Method comparison ({base})")
        save_fig(fig, out, f"methods_{base}")


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="TMVA-output/results")
    ap.add_argument("--out", default="report")
    args = ap.parse_args()

    results, out = Path(args.results), Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    summary = results / "tmva_summary.csv"
    if not summary.exists():
        sys.exit(f"{summary} not found: run trainTMVA/applyTMVA first")
    df = parse_tags(pd.read_csv(summary))
    needed = {"events", "auc", "rej_effS50", "nfolds"}
    if not needed <= set(df.columns):
        sys.exit("tmva_summary.csv has the old column layout: rerun applyTMVA "
                 "(make reapply) so it writes the new columns "
                 f"({', '.join(sorted(needed - set(df.columns)))} missing)")
    print(f"{len(df)} rows in {summary}")

    cv_summary(df, out)
    scans(df, out)
    method_comparison(df, out)
    smearing(df, out)
    ablation(df, out)
    importance(results, out)
    holdout(df, out)
    print(f"Done: {out}/figures and {out}/tables")


if __name__ == "__main__":
    main()
