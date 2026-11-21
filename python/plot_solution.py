#
# Plots of a HALO solution, from the `solution_plot_data_*.json` file that the
# solver writes with `"generate_solution_plots": true` (it also runs this script):
#
#  * epoch_drift_<case>.pdf/.png: how far the solution's epochs at the constraint
#    points (where segments 2k-1 and 2k end) moved from the initial guess's.
#  * defects_<case>.pdf/.png: the position and velocity defects at the constraint
#    points, for the initial guess and the solution.
#
# Usage: python python/plot_solution.py solution_plot_data_<case>.json [--outdir DIR]
#
# Or, to compare the epoch drift of several solutions (e.g., with different options)
# in one plot (epoch_drift_comparison.pdf/.png, or --name):
#
#   python python/plot_solution.py --compare a.json b.json --labels "default" "offsets" [--name NAME]
#
# Needs only numpy and matplotlib (in the pixi environment).
#

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")  # (no display needed: the solver runs this)
import matplotlib.pyplot as plt
import numpy as np

SEC2DAY = 86400.0
FLOOR = 1.0e-16  # zero on the log scale

GUESS_COLOR = "0.6"
SOLUTION_COLOR = "tab:blue"
VELOCITY_COLOR = "tab:red"


def load(filename):
    with open(filename) as f:
        d = json.load(f)
    for k in ["et_guess", "dr_guess", "dv_guess", "et_solution", "dr_solution", "dv_solution"]:
        d[k] = np.atleast_1d(np.asarray(d[k], dtype=float))
    return d


def save(fig, outdir, name):
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(path, dpi=200 if ext == "png" else None, bbox_inches="tight")
        print(f"   {path}")
    plt.close(fig)


def plot_epoch_drift(d, outdir):
    """The solution's epochs at the constraint points minus the initial guess's."""
    t = (d["et_guess"] - d["et_ref"]) / SEC2DAY
    drift = (d["et_solution"] - d["et_guess"]) / 3600.0
    rms = np.sqrt(np.mean(drift**2))
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    ax.axhline(0.0, color="0.75", linewidth=0.8)
    ax.plot(t, drift, color=SOLUTION_COLOR, linewidth=0.8, marker=".", markersize=2.0,
            label=f"solution $-$ initial guess (rms {rms:.2f} h, max {np.max(np.abs(drift)):.2f} h)")
    ax.set_xlabel("Constraint point epoch (days from the reference epoch)")
    ax.set_ylabel("Epoch drift (hours)")
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    save(fig, outdir, f"epoch_drift_{d['case']}")


def plot_defects(d, outdir):
    """|dr| and |dv| at the constraint points, initial guess and solution (log scale)."""
    tg = (d["et_guess"] - d["et_ref"]) / SEC2DAY
    ts = (d["et_solution"] - d["et_ref"]) / SEC2DAY
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 4.2), sharex=True)
    for ax, key, unit in [(ax1, "dr", "km"), (ax2, "dv", "km/s")]:
        g = np.maximum(d[f"{key}_guess"], FLOOR)
        s = np.maximum(d[f"{key}_solution"], FLOOR)
        ax.semilogy(tg, g, color=GUESS_COLOR, linestyle="none", marker=".", markersize=2.0,
                    label=f"initial guess (max {np.max(g):.1e} {unit})")
        ax.semilogy(ts, s, color=SOLUTION_COLOR if key == "dr" else VELOCITY_COLOR,
                    linestyle="none", marker=".", markersize=2.0,
                    label=f"solution (max {np.max(s):.1e} {unit})")
        ax.set_ylabel(f"|d{key[1]}| ({unit})")
        ax.grid(True, which="major", linewidth=0.4, alpha=0.5)
        ax.legend(loc="center right", fontsize=8, frameon=True, framealpha=0.9)
    ax2.set_xlabel("Constraint point epoch (days from the reference epoch)")
    ax1.set_xlim(min(tg[0], ts[0]), max(tg[-1], ts[-1]))
    fig.tight_layout()
    save(fig, outdir, f"defects_{d['case']}")


def running_mean(t, y, window):
    """Mean of y over a centered window (same units as t), at each t."""
    c = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(t, t - window / 2.0, side="left")
    hi = np.searchsorted(t, t + window / 2.0, side="right")
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)


def plot_epoch_drift_comparison(datasets, labels, outdir, name, window=60.0):
    """The epoch drift of several solutions: each as a thin line, with its running
    mean over `window` days (bold)."""
    colors = ["tab:orange", "tab:blue", "tab:green", "tab:red", "tab:purple"]
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.axhline(0.0, color="0.75", linewidth=0.8)
    span = max((d["et_guess"][-1] - d["et_guess"][0]) / SEC2DAY for d in datasets)
    years = span > 1000.0
    tmax = 0.0
    for i, (d, label) in enumerate(zip(datasets, labels)):
        c = colors[i % len(colors)]
        t = (d["et_guess"] - d["et_ref"]) / SEC2DAY
        drift = (d["et_solution"] - d["et_guess"]) / 3600.0
        rms = np.sqrt(np.mean(drift**2))
        tp = t / 365.25 if years else t
        tmax = max(tmax, tp[-1])
        ax.plot(tp, drift, color=c, linewidth=0.4, alpha=0.45)
        ax.plot(tp, running_mean(t, drift, window), color=c, linewidth=1.6,
                label=f"{label} (rms {rms:.2f} h, max {np.max(np.abs(drift)):.2f} h)")
    ax.set_xlabel("Constraint point epoch (" + ("years" if years else "days") + " from the reference epoch)")
    ax.set_ylabel("Epoch drift (hours)")
    ax.set_xlim(0.0, tmax)
    ax.grid(True, linewidth=0.4, alpha=0.5)
    ax.legend(loc="best", fontsize=8, frameon=True, framealpha=0.9,
              title=f"thin: each constraint point; bold: {window:g}-day mean", title_fontsize=7)
    fig.tight_layout()
    save(fig, outdir, name)


def main():
    parser = argparse.ArgumentParser(description="Plots of a HALO solution (epoch drift and defects).")
    parser.add_argument("datafile", nargs="?", help="solution_plot_data_*.json file written by the solver")
    parser.add_argument("--outdir", default=None, help="where to write the plots (default: the data file's folder)")
    parser.add_argument("--compare", nargs="+", metavar="FILE",
                        help="compare the epoch drift of these solution_plot_data_*.json files in one plot")
    parser.add_argument("--labels", nargs="+", help="the labels of the --compare files")
    parser.add_argument("--name", default="epoch_drift_comparison", help="file name of the --compare plot")
    args = parser.parse_args()
    matplotlib.rcParams["font.family"] = "serif"
    matplotlib.rcParams["font.size"] = 9
    if args.compare:
        labels = args.labels if args.labels else [os.path.basename(f) for f in args.compare]
        if len(labels) != len(args.compare):
            parser.error("--labels must have one label for each --compare file")
        outdir = args.outdir if args.outdir else "."
        os.makedirs(outdir, exist_ok=True)
        plot_epoch_drift_comparison([load(f) for f in args.compare], labels, outdir, args.name)
        return
    if not args.datafile:
        parser.error("a data file (or --compare) is required")
    d = load(args.datafile)
    outdir = args.outdir if args.outdir else (os.path.dirname(args.datafile) or ".")
    os.makedirs(outdir, exist_ok=True)
    plot_epoch_drift(d, outdir)
    plot_defects(d, outdir)


if __name__ == "__main__":
    main()
