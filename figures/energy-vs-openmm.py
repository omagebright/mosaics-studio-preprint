#!/usr/bin/env python3
"""fig:energy-vs-openmm -- the three programs measured against OpenMM.

Peter Minary asked (8 September 2026) for the comparison to be shown as three
sets, one per program, so that the MOSAICS spread can be read against the
spread the two established engines show among themselves.

OpenMM is the reference rather than sander for three reasons, all of which
this figure makes visible:

  * sander prints its per-term energies to four decimals and its total is the
    sum of those rounded terms, so a sander difference carries up to 2.0e-4
    kcal/mol of print quantisation -- the size of the differences being
    compared.  OpenMM's writer resolves to 5e-13.
  * sander was never run on the CHARMM36 cells, so any sander-referenced
    comparison silently drops them, and they are the closest agreement here.
  * OpenMM and GROMACS return the same Coulomb constant, so referencing
    OpenMM leaves only sander needing the rescale of Section 'constant'.

Values come from data/MATRIX.json, whose electrostatics are already on the
CODATA constant 332.0637133.  Nothing numeric is hard-coded except the
physical constants, the print floors, and the F-12 angle offset, each of which
is cited in the caption.
"""
import json
import pathlib
import statistics
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _style as S

DATA = HERE.parent / "data"
OUT = HERE / "energy-vs-openmm"

KT300 = 0.0019872041 * 300.0
BAND = 1.0e-3
# Print resolution of each format, kcal/mol, as figures/energy-agreement.py.
FLOOR = {"sander": 2.0e-4, "gromacs": 2.4e-7, "mosaics": 0.0}
# F-12: the AMBER topology builder writes equilibrium angles through a
# seven-digit pi/180; worth this much on 1efs and the largest remaining
# difference on every passing AMBER row.  evidence/findings.txt.
F12 = 1.8e-4

LABEL = {"1efs": "1efs", "1BNA": "1BNA", "2KOC": "2KOC", "dna_acgt": "d(ACGT)",
         "rna_acgu": "r(ACGU)", "1UBQ": "1UBQ", "5PTI": "5PTI",
         "peptide": "AGSV", "tax9": "Tax9", "1HHK": "1HHK", "1F7Y": "1F7Y"}
GROUP = {"DNA:RNA hybrid": "DNA and hybrid", "B-DNA duplex": "DNA and hybrid",
         "DNA 4-mer, 5′/3′ termini": "DNA and hybrid",
         "RNA hairpin": "RNA", "RNA 4-mer, 5′/3′ termini": "RNA",
         "protein": "proteins and peptides", "peptide": "proteins and peptides",
         "protein-peptide": "complexes", "protein-RNA": "complexes"}
GROUP_ORDER = ["DNA and hybrid", "RNA", "proteins and peptides", "complexes"]

MOSAICS_COLOUR = "#6f42a1"
PANELS = [("mosaics", "MOSAICS", MOSAICS_COLOUR, "o"),
          ("sander",  "sander",  S.SANDER,       "D"),
          ("gromacs", "GROMACS", S.GROMACS,      "^")]


def cells():
    rows = json.loads((DATA / "MATRIX.json").read_text())
    for r in rows:
        r["group"] = GROUP[r["klass"]]
        r["cell"] = f"{LABEL.get(r['system'], r['system'])} {r['force_field']}"
    rows.sort(key=lambda r: (GROUP_ORDER.index(r["group"]), r["system"],
                             r["force_field"]))
    return rows


def main():
    S.apply()
    rows = cells()
    x = list(range(len(rows)))

    fig, axes = plt.subplots(3, 1, figsize=(7.4, 6.8), sharex=True,
                             gridspec_kw={"hspace": 0.30})

    for ax, (key, name, colour, marker) in zip(axes, PANELS):
        S.ygrid(ax)
        ax.axhline(KT300, color=S.LIGHT, lw=0.7, ls="--", zorder=1)
        ax.axhline(BAND, color=S.LIGHT, lw=0.6, ls=":", zorder=1)

        floor = FLOOR[key]
        if floor > 0:
            ax.axhspan(1e-9, floor, color=S.GRID, alpha=0.55, zorder=0)
            txt = ("2×10$^{-4}$" if key == "sander" else "2×10$^{-7}$")
            ax.text(0.15, floor * 1.5,
                    f"print resolution of the {name} output file, {txt} kcal/mol",
                    ha="left", va="bottom", fontsize=6.2, color=S.TEXT)
        if key == "mosaics":
            ax.axhline(F12, color="#8a6d3b", lw=0.7, ls="-.", zorder=1)
            ax.text(0.1, F12 / 1.5,
                    "F-12, 1.8×10$^{-4}$: equilibrium angles as the AMBER\n"
                    "topology builder writes them",
                    ha="left", va="top", fontsize=6.2, color="#8a6d3b",
                    linespacing=1.35)

        xs, ys, filled = [], [], []
        for i, r in enumerate(rows):
            v, o = r.get(key), r.get("openmm")
            if not isinstance(v, (int, float)) or not isinstance(o, (int, float)):
                continue
            d = abs(v - o)
            xs.append(i); ys.append(max(d, 1e-9)); filled.append(d >= floor)
        ax.scatter([xs[i] for i in range(len(xs)) if filled[i]],
                   [ys[i] for i in range(len(xs)) if filled[i]],
                   s=26, marker=marker, facecolor=colour, edgecolor="white",
                   linewidth=0.5, zorder=3)
        ax.scatter([xs[i] for i in range(len(xs)) if not filled[i]],
                   [ys[i] for i in range(len(xs)) if not filled[i]],
                   s=26, marker=marker, facecolor="none", edgecolor=colour,
                   linewidth=0.8, zorder=3)

        med = statistics.median([abs(r[key] - r["openmm"]) for r in rows
                                 if isinstance(r.get(key), (int, float))
                                 and isinstance(r.get("openmm"), (int, float))])
        ax.set_yscale("log"); ax.set_ylim(1e-7, 2.0)
        ax.set_xlim(-0.7, len(rows) - 0.3)
        ax.set_yticks([1e-6, 1e-4, 1e-2, 1e0])
        ax.set_title(f"({'abc'[PANELS.index((key, name, colour, marker))]})   "
                     f"{name} − OpenMM      n = {len(xs)},  "
                     f"median {med:.1e} kcal/mol".replace("e-04", "×10$^{-4}$")
                                                 .replace("e-05", "×10$^{-5}$")
                                                 .replace("e-06", "×10$^{-6}$"),
                     loc="left", fontsize=8, color=S.TEXT, pad=4)

        for g in GROUP_ORDER[1:]:
            first = min(i for i, r in enumerate(rows) if r["group"] == g)
            ax.axvline(first - 0.5, color=S.GRID, lw=0.6, zorder=0)

    for ax in axes:
        ax.text(len(rows) - 0.3, KT300 * 1.35, "kT at 300 K, 0.596",
                ha="right", va="bottom", fontsize=6.2, color=S.LIGHT)
    for g in GROUP_ORDER:
        idx = [i for i, r in enumerate(rows) if r["group"] == g]
        lo, hi = min(idx) - 0.45, max(idx) + 0.45
        axes[0].annotate(g, xy=((lo + hi) / 2, 1.34), xycoords=("data", "axes fraction"),
                         ha="center", va="bottom", fontsize=7, color=S.TEXT,
                         annotation_clip=False)
        axes[0].annotate("", xy=(lo, 1.30), xytext=(hi, 1.30),
                         xycoords=("data", "axes fraction"),
                         textcoords=("data", "axes fraction"),
                         arrowprops=dict(arrowstyle="-", color=S.GRID, lw=0.8),
                         annotation_clip=False)

    fig.supylabel("|  total energy − OpenMM  |   (kcal/mol)", fontsize=8,
                  color=S.TEXT, x=0.045)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([r["cell"] for r in rows], rotation=60,
                             ha="right", fontsize=6.2)

    fig.legend(handles=[Line2D([], [], marker="o", ls="none", mfc="none",
                               mec=S.TEXT, ms=5,
                               label="hollow: below the print resolution of "
                                     "that program's own output file")],
               loc="lower left", bbox_to_anchor=(0.115, -0.055),
               frameon=False, fontsize=6.5)

    pdf, png = S.save(fig, OUT)
    print(f"wrote {pdf.name} and {png.name}")
    for key, name, _, _ in PANELS:
        v = sorted(abs(r[key] - r["openmm"]) for r in rows
                   if isinstance(r.get(key), (int, float))
                   and isinstance(r.get("openmm"), (int, float)))
        print(f"  {name:8s} n={len(v):2d} median={statistics.median(v):.3e} "
              f"max={max(v):.3e}")


if __name__ == "__main__":
    main()
