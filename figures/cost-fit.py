#!/usr/bin/env python3
"""fig:cost-fit -- wall clock against MC step count, 1efs heteroduplex, 818 atoms.

Reads work/preprint-I/runs/2026-08-31-cost-818/result.json, written by
measure.py beside it.  Every timing is plotted: the five repeats at each N as
open markers and their median as a filled marker.  The line is a least-squares
fit t(N) = a + b N on the medians, recomputed here and checked against the
slope, intercept and R^2 the run recorded.
"""
import json
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _style as S

RUN = pathlib.Path("/Users/bright/Documents/MOSAICS/work/preprint-I/runs/"
                   "2026-08-31-cost-818/result.json")
OUT = HERE / "cost-fit"


def fit(N, t):
    n = len(N)
    mx, my = sum(N) / n, sum(t) / n
    sxx = sum((x - mx) ** 2 for x in N)
    sxy = sum((x - mx) * (y - my) for x, y in zip(N, t))
    b = sxy / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in t)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(N, t))
    return a, b, 1 - ss_res / ss_tot


def main():
    d = json.loads(RUN.read_text())
    rows = d["rows"]
    N = [r["steps"] for r in rows]
    med = [r["median"] for r in rows]
    a, b, r2 = fit(N, med)
    for name, got, want in (("per_step_seconds", b, d["per_step_seconds"]),
                            ("setup_seconds", a, d["setup_seconds"]),
                            ("r_squared", r2, d["r_squared"])):
        if abs(got - want) > abs(want) * 1e-9:
            sys.exit(f"{name}: recomputed {got!r}, result.json holds {want!r}")

    S.apply()
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.17)
    xs = [-2, max(N) + 2]
    ax.plot(xs, [a + b * x for x in xs], "-", color=S.BLUE, lw=0.9, zorder=2)
    n_raw = 0
    for r in rows:
        ax.plot([r["steps"]] * len(r["times"]), r["times"], marker="o", ls="none",
                ms=3, mfc="none", mec=S.BLUE, mew=0.6, zorder=3, alpha=0.8)
        n_raw += len(r["times"])
    ax.plot(N, med, marker="o", ls="none", ms=5, mfc=S.BLUE, mec=S.BLUE, mew=0.6,
            zorder=4)
    ax.set_xlim(-4, max(N) + 4)
    ax.set_xlabel("Monte Carlo steps $N$")
    ax.set_ylabel("wall clock (s)")
    S.ygrid(ax)
    ax.text(0.97, 0.05,
            f"$t(N) = a + bN$\n"
            f"$b$ = {b * 1000:.3f} ms per step\n"
            f"$a$ = {a * 1000:.1f} ms\n"
            f"$R^2$ = {r2:.5f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
            color="#222222", linespacing=1.35)
    handles = [
        Line2D([], [], marker="o", ls="none", ms=5, mfc=S.BLUE, mec=S.BLUE, mew=0.6,
               label=f"median of {d['repeats']} repeats"),
        Line2D([], [], marker="o", ls="none", ms=3, mfc="none", mec=S.BLUE, mew=0.6,
               label="single repeat"),
        Line2D([], [], color=S.BLUE, lw=0.9, label="least squares on medians"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=7,
              handletextpad=0.5, borderaxespad=0.2, labelspacing=0.3,
              bbox_to_anchor=(0.2, 1.0))
    pdf, png = S.save(fig, OUT)
    plt.close(fig)

    print(f"wrote {pdf}\nwrote {png}")
    print(f"files read:\n   {RUN}")
    print(f"markers: {len(N)} medians + {n_raw} repeats = {len(N) + n_raw}")
    print(f"per step {b * 1000:.3f} ms   setup {a * 1000:.2f} ms   R^2 {r2:.6f}   "
          f"atoms {d['atoms']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
