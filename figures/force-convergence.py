#!/usr/bin/env python3
"""fig:force-convergence -- analytic force against a central difference.

Six panels, one per all-terms Cartesian sweep named by the existing script
paper/figures/force-convergence.py (its PANELS list).  Each point is a row of
paper/evidence/forces/<stem>.all.tsv: the largest |numerical - analytic| over
every atom and component (Ha/Bohr) at that finite-difference step (Bohr).  A
central difference is wrong at order step^2, so a correct analytic force gives
slope 2 on the falling branch.  The exponent is fitted by least squares over
the fixed window 1e-2 to 3e-4 Bohr, the same window forces.txt declares.
"""
import csv
import math
import pathlib
import sys

import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _style as S

FDIR = pathlib.Path("/Users/bright/Documents/MOSAICS/paper/evidence/forces")
OUT = HERE / "force-convergence"

PANELS = [
    ("1efs_bsc1",         "1efs, parmbsc1"),
    ("1efs_ol21_ol3",     "1efs, OL21/OL3"),
    ("dna_acgt_ol21_ol3", "d(ACGT), OL21/OL3"),
    ("rna_acgu_ol21_ol3", "r(ACGU), OL21/OL3"),
    ("peptide_ff14sb",    "ACE-AGSV-NME, ff14SB"),
    ("peptide_ff19sb",    "ACE-AGSV-NME, ff19SB"),
]
WINDOW = (3.0e-4, 1.0e-2)


def read(stem):
    path = FDIR / f"{stem}.all.tsv"
    x, y = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            x.append(float(row["step_bohr"]))
            y.append(float(row["max_abs_diff_ha_per_bohr"]))
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [y[i] for i in order], path


def slope(x, y):
    pts = [(math.log10(a), math.log10(b)) for a, b in zip(x, y)
           if WINDOW[0] * 0.999 <= a <= WINDOW[1] * 1.001]
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    return sxy / sxx, n


def main():
    S.apply()
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 3.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.92, bottom=0.15,
                        wspace=0.12, hspace=0.32)
    XLIM = (1.5e-7, 2.5e-2)
    YLIM = (1.0e-11, 2.0e-2)
    files, n_pts, report = [], 0, []
    for i, (stem, label) in enumerate(PANELS):
        ax = axes[i // 3][i % 3]
        x, y, path = read(stem)
        files.append(path)
        n_pts += len(x)
        # slope-2 reference through the first (coarsest) point
        x0, y0 = x[-1], y[-1]
        gx = [XLIM[0] * 4, x0 * 1.6]
        ax.plot(gx, [y0 * (g / x0) ** 2 for g in gx], color=S.LIGHT, lw=0.8,
                ls="-", zorder=1)
        ax.plot(x, y, color=S.BLUE, lw=0.8, marker="o", ms=3, mfc=S.BLUE,
                mec=S.BLUE, mew=0.6, zorder=3, clip_on=False)
        m, n = slope(x, y)
        report.append((stem, m, n, min(y), max(y)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        ax.set_xticks([1e-6, 1e-4, 1e-2])
        ax.set_yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
        ax.tick_params(axis="both", which="minor", length=0)
        S.ygrid(ax)
        ax.text(0.04, 0.95, label, transform=ax.transAxes, fontsize=7,
                color=S.TEXT, ha="left", va="top")
        ax.text(0.04, 0.80, f"slope {m:.3f}", transform=ax.transAxes, fontsize=7,
                color="#222222", ha="left", va="top")
        S.letter(ax, f"({'abcdef'[i]})", dx=-0.02, dy=0.03)
    fig.text(0.535, 0.02, "finite-difference step (Bohr)", ha="center", va="bottom")
    fig.text(0.012, 0.535, "max |numerical − analytic| force  (Ha/Bohr)",
             rotation=90, ha="left", va="center")
    pdf, png = S.save(fig, OUT)
    plt.close(fig)

    print(f"wrote {pdf}\nwrote {png}")
    print("files read:")
    for f in files:
        print(f"   {f}")
    print(f"panels {len(PANELS)}   markers {n_pts}")
    for stem, m, n, lo, hi in report:
        print(f"   {stem:<20} slope {m:.3f} over {n} points in "
              f"[{WINDOW[0]:.0e}, {WINDOW[1]:.0e}]   residual range {lo:.2e} .. {hi:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
