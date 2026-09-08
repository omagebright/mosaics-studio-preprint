#!/usr/bin/env python3
"""fig:term-agreement -- MOSAICS against three engines, term by term.

The totals moved to figures/energy-vs-openmm.py on 8 September 2026, which
shows the three programs against a common reference the way the comparison was
asked for. What stays here is the part that figure cannot carry:

    the comparisons resolved into the terms the engines share: bond,
    angle, dihedral (proper + improper), Lennard-Jones, electrostatics, and
    their sum.  Per-term values are re-derived from each engine's own output
    file, located by the cell-resolution logic of the existing script
    work/preprint-I/figures/energy-agreement.py (imported, not copied).

Nothing numeric is hard-coded except the physical constants and the print
resolution of each engine's output format.  Every path is read from the
evidence tree; nothing under it is written.
"""
import contextlib
import importlib.util
import io
import json
import pathlib
import re
import statistics
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _style as S

REPO = pathlib.Path("/Users/bright/Documents/MOSAICS")
MATRIX = REPO / "work/preprint-I/energies/MATRIX.json"
EXISTING = REPO / "work/preprint-I/figures/energy-agreement.py"
OUT = HERE / "term-agreement"

# Reuse the existing script's cell-resolution and constants.
spec = importlib.util.spec_from_file_location("existing_ea", EXISTING)
ea = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ea)
K_CODATA, K_MOSAICS, K_SANDER, KJ = ea.K_CODATA, ea.K_MOSAICS, ea.K_SANDER, ea.KJ_PER_KCAL
KT300 = 0.0019872041 * 300.0
BAND = 1.0e-3

# The two OL24 cells: the existing script reads a derived TSV that carries
# MOSAICS on its raw constant and no Coulomb split.  The engine logs below are
# the runs those TSVs were derived from (FINDINGS-B1-two-OL24-cells-2026-09-05.md,
# normalize_matrix_constant.py FILE_OVERRIDE); both print a Total Coulomb term,
# so the CODATA rescale is reproducible.  OpenMM for these cells is the
# per-term summary beside the TSV, which has no LJ/Coulomb split.
OL24 = {
    ("1efs", "OL24/OL3"): {
        "mosaics": "record/08-mosaics-source/MOSAICS/examples/1efs_heteroduplex/"
                   "mosaics_ol24_ol3_vac.log",
        "openmm": "record/TODO/01-force-field-resolved-heteroduplex-landscape/"
                  "validation/reference_engines/ol24_ol3/openmm_amber_compare/"
                  "openmm_energy_summary.tsv",
    },
    ("dna_acgt", "OL24"): {
        "mosaics": "engine/examples/pure_terminal_controls/"
                   "mcmc_dna_acgt_ol24_ol3_terminal_vac.log",
        "openmm": "record/TODO/01-force-field-resolved-heteroduplex-landscape/"
                  "validation/reference_engines/ol24_ol3_terminal/pure_controls/"
                  "dna_acgt/openmm_amber_compare/openmm_energy_summary.tsv",
    },
}

TERMS = ["bond", "angle", "dihedral", "lj", "elec", "nonbonded"]
TERM_LABEL = {"bond": "bond", "angle": "angle", "dihedral": "dihedral",
              "lj": "Lennard-Jones", "elec": "electrostatics",
              "nonbonded": "LJ + elec"}

SHORT = {"1efs": "1efs", "1BNA": "1BNA", "dna_acgt": "d(ACGT)", "2KOC": "2KOC",
         "rna_acgu": "r(ACGU)", "1UBQ": "1UBQ", "5PTI": "5PTI",
         "peptide": "AGSV", "tax9": "Tax9", "1HHK": "1HHK", "1F7Y": "1F7Y"}
GROUP = {"DNA:RNA hybrid": "DNA and hybrid", "B-DNA duplex": "DNA and hybrid",
         "DNA 4-mer, 5′/3′ termini": "DNA and hybrid",
         "RNA hairpin": "RNA", "RNA 4-mer, 5′/3′ termini": "RNA",
         "protein": "proteins and peptides", "peptide": "proteins and peptides",
         "protein-peptide": "complexes", "protein-RNA": "complexes"}
GROUP_ORDER = ["DNA and hybrid", "RNA", "proteins and peptides", "complexes"]

# Print resolution of each reference format, kcal/mol.  sander prints four
# decimals; GROMACS six decimals in kJ/mol; OpenMM is half a unit in the last
# printed place, which is 5e-13 for the twelve-decimal writer.
FLOOR_SANDER = 2.0e-4
FLOOR_GROMACS = 2.4e-7


def read(rel):
    return (REPO / rel).read_text(errors="replace")


def num(s):
    return float(s)


# --------------------------------------------------------------- parsers
def mosaics_terms(rel):
    text = read(rel)

    def grab(key):
        hits = re.findall(r"^\s*-\s*" + re.escape(key) + r"\s+(-?[\d.]+(?:[eE][-+]?\d+)?)\s*$",
                          text, re.M)
        return num(hits[-1]) if hits else None

    bond, angle, tors = grab("Bond energy"), grab("Bend energy"), grab("Torsion energy")
    onefour, inter, coul = grab("Onefour energy"), grab("Inter energy"), grab("Total Coulomb energy")
    total = grab("Total energy")
    if None in (bond, angle, tors, onefour, inter, coul, total):
        raise ValueError(f"incomplete MOSAICS energy block in {rel}")
    elec = coul * (K_CODATA / K_MOSAICS)
    lj = onefour + inter - coul
    out = {"bond": bond, "angle": angle, "dihedral": tors, "lj": lj, "elec": elec,
           "nonbonded": lj + elec, "total": total - coul + elec}
    cmap = grab("Cmap energy")
    if cmap is not None:
        out["cmap"] = cmap
    return out


def sander_terms(rel):
    text = read(rel)
    block = text[text.rfind("FINAL RESULTS"):]

    def col(name):
        hit = re.search(re.escape(name) + r"\s*=\s*(-?\d+\.\d+)", block)
        return num(hit.group(1)) if hit else None

    lj = col("VDWAALS") + col("1-4 VDW")
    elec = (col("EEL") + col("1-4 EEL")) * (K_CODATA / K_SANDER)
    out = {"bond": col("BOND"), "angle": col("ANGLE"), "dihedral": col("DIHED"),
           "lj": lj, "elec": elec, "nonbonded": lj + elec}
    if col("CMAP") is not None:
        out["cmap"] = col("CMAP")
    out["total"] = sum(v for k, v in out.items() if k != "nonbonded")
    return out, FLOOR_SANDER


def gromacs_terms(rel):
    text = read(rel)
    legend = {int(m.group(1)): m.group(2)
              for m in re.finditer(r'@ s(\d+) legend "([^"]+)"', text)}
    rows = [ln for ln in text.splitlines() if ln.strip() and ln.strip()[0] not in "#@"]
    fields = rows[-1].split()
    v = {legend[i]: num(fields[i + 1]) for i in legend if i + 1 < len(fields)}
    lj = (v.get("LJ-14", 0.0) + v.get("LJ (SR)", 0.0)) / KJ
    elec = (v.get("Coulomb-14", 0.0) + v.get("Coulomb (SR)", 0.0)) / KJ
    out = {"bond": v["Bond"] / KJ, "angle": v["Angle"] / KJ,
           "dihedral": (v.get("Proper Dih.", 0.0) + v.get("Per. Imp. Dih.", 0.0)
                        + v.get("Improper Dih.", 0.0)) / KJ,
           "lj": lj, "elec": elec, "nonbonded": lj + elec}
    out["total"] = (v["Potential"] / KJ if "Potential" in v
                    else sum(out[k] for k in TERMS if k != "nonbonded"))
    return out, FLOOR_GROMACS


def openmm_floor(text):
    dp = max(len(x.split(".")[1]) for x in re.findall(r"-?\d+\.\d+", text))
    return 0.5 * 10.0 ** -dp


def openmm_terms(rel):
    text = read(rel)
    if rel.endswith(".tsv"):                       # OL24 per-term summary
        v = {}
        for line in text.splitlines()[1:]:
            k, x = line.split("\t")[:2]
            v[k] = num(x)
        out = {"bond": v["bond"], "angle": v["angle"], "dihedral": v["dihedral"],
               "nonbonded": v["nonbonded_total"], "total": v["total"]}
        return out, openmm_floor(text)
    if "openmm force" in text:                     # CHARMM36 writer
        v = {}
        for line in text.splitlines():
            hit = re.match(r"\s{2}(.+?)\s{2,}(-?\d+\.\d+)\s", line + " ")
            if hit:
                v[hit.group(1).strip()] = num(hit.group(2))
        lj = v["Lennard-Jones"] + v.get("1-4 Lennard-Jones", 0.0)
        elec = v["electrostatics"]
        out = {"bond": v["bond"],
               "angle": v["angle"] + v.get("Urey-Bradley (1-3)", 0.0),
               "dihedral": v["torsion"] + v.get("improper (harmonic)", 0.0),
               "lj": lj, "elec": elec, "nonbonded": lj + elec, "total": v["TOTAL"]}
        if v.get("correction map"):
            out["cmap"] = v["correction map"]
        return out, openmm_floor(text)
    force = {"HarmonicBondForce": "bond", "HarmonicAngleForce": "angle",
             "PeriodicTorsionForce": "dihedral", "CMAPTorsionForce": "cmap",
             "NonbondedForce": "nonbonded", "TOTAL": "total",
             "Lennard-Jones (charges zeroed)": "lj",
             "Coulomb (epsilons zeroed)": "elec"}
    out = {}
    for key, name in force.items():
        hit = re.search(r"^" + re.escape(key) + r"\s+(-?\d+\.\d+)\s*$", text, re.M)
        if hit:
            out[name] = num(hit.group(1))
    if "total" not in out:
        out["total"] = sum(out[k] for k in ("bond", "angle", "dihedral", "cmap", "nonbonded")
                           if k in out)
    return out, openmm_floor(text)


PARSER = {"sander": sander_terms, "openmm": openmm_terms, "gromacs": gromacs_terms}


# ------------------------------------------------------------- resolution
def load():
    matrix = json.loads(MATRIX.read_text())
    # build() notes that the two OL24 TSVs are on the raw constant; those two
    # cells are re-read from their engine logs below, so the note is expected.
    with contextlib.redirect_stderr(io.StringIO()):
        resolved = {(e["cell"]["system"], e["cell"]["force_field"]): e for e in ea.build()}
    cells = []
    files_read = {str(MATRIX), str(EXISTING)}
    for cell in matrix:
        key = (cell["system"], cell["force_field"])
        entry = resolved[key]
        if key in OL24:
            mfile = OL24[key]["mosaics"]
            ref_files = {"openmm": OL24[key]["openmm"]}
        else:
            mfile = entry["mosaics_file"]
            ref_files = {eng: path for eng, (_, path, _) in entry["refs"].items()}
        mos = mosaics_terms(mfile)
        files_read.add(str(REPO / mfile))
        refs = {}
        for eng, path in ref_files.items():
            if cell.get(eng) is None:
                continue
            terms, floor = PARSER[eng](path)
            refs[eng] = (terms, floor, path)
            files_read.add(str(REPO / path))
        cells.append({"cell": cell, "mosaics": mos, "mfile": mfile, "refs": refs})
    return cells, sorted(files_read)


def check(cells):
    """The per-term files must reproduce the totals MATRIX.json publishes."""
    worst = 0.0
    for e in cells:
        c = e["cell"]
        gap = abs(e["mosaics"]["total"] - c["mosaics"])
        worst = max(worst, gap)
        if gap > 1e-6:
            print(f"  warn {c['system']} {c['force_field']}: MOSAICS file total "
                  f"{e['mosaics']['total']:.9f} vs MATRIX {c['mosaics']:.9f}", file=sys.stderr)
        for eng, (terms, floor, path) in e["refs"].items():
            gap = abs(terms["total"] - c[eng])
            worst = max(worst, gap)
            if gap > 1e-6:
                print(f"  warn {c['system']} {c['force_field']} {eng}: file total "
                      f"{terms['total']:.9f} vs MATRIX {c[eng]:.9f}", file=sys.stderr)
    return worst


# --------------------------------------------------------------- drawing
def draw(cells):
    S.apply()
    fig = plt.figure(figsize=(3.4, 3.3))
    gs = gridspec.GridSpec(1, 1, left=0.20, right=0.985, top=0.90, bottom=0.27)
    axb = fig.add_subplot(gs[0])
    # The totals are drawn by figures/energy-vs-openmm.py. They are still
    # tallied below, on a throwaway axes, because main() reports on them.
    axa = fig.add_subplot(gs[0], frame_on=False)
    axa.set_visible(False)

    # ---- (a) totals, x positions grouped by class
    x, pos, ticks, groups = 0.0, {}, [], {}
    prev = None
    for e in cells:
        c = e["cell"]
        g = GROUP[c["klass"]]
        if prev is not None and g != prev:
            x += 1.3
        prev = g
        key = (c["system"], c["force_field"])
        pos[key] = x
        ticks.append((x, f"{SHORT[c['system']]} {c['force_field']}"))
        groups.setdefault(g, []).append(x)
        x += 1.0
    right = max(pos.values()) + 0.8
    left = -0.8

    bottom, top = 1.0e-13, 4.0
    offset = {"sander": -0.27, "openmm": 0.0, "gromacs": 0.27}
    floor_of = {"sander": FLOOR_SANDER, "gromacs": FLOOR_GROMACS}
    totals = []
    for e in cells:
        c = e["cell"]
        key = (c["system"], c["force_field"])
        for eng in ("sander", "openmm", "gromacs"):
            if c.get(eng) is None:
                continue
            d = abs(c["mosaics"] - c[eng])
            floor = floor_of.get(eng) or e["refs"][eng][1]
            totals.append({"key": key, "engine": eng, "d": d, "floor": floor})
            low = d < floor
            col = S.ENGINE_COLOUR[eng]
            axa.plot(pos[key] + offset[eng], max(d, bottom * 1.6),
                     marker=S.ENGINE_MARKER[eng], ms=5, ls="none", zorder=3,
                     mfc="none" if low else col, mec=col, mew=0.6, clip_on=False)

    for ax in (axb,):
        ax.set_yscale("log")
        ax.set_ylim(bottom, top)
        ax.set_yticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0])
        ax.axhline(KT300, color=S.LIGHT, lw=0.7, ls=(0, (4, 2.5)), zorder=1)
        ax.axhline(BAND, color=S.LIGHT, lw=0.7, ls=(0, (1, 2)), zorder=1)
        S.ygrid(ax)
        ax.tick_params(axis="x", length=0, pad=2)
    axb.text(len(TERMS) - 0.45, KT300 / 1.4, f"kT at 300 K, {KT300:.3f}",
             fontsize=6.5, color=S.TEXT, ha="right", va="top")
    axb.text(len(TERMS) - 0.45, BAND / 1.4, "0.001", fontsize=6.5,
             color=S.TEXT, ha="right", va="top")
    axb.set_ylabel("|$E_{\\mathrm{MOSAICS}}$ − $E_{\\mathrm{reference}}$|  (kcal/mol)")

    # ---- (b) per-term
    jitter = {"sander": -0.24, "openmm": 0.0, "gromacs": 0.24}
    terms = []
    for e in cells:
        c = e["cell"]
        key = (c["system"], c["force_field"])
        for eng, (ref, floor, path) in e["refs"].items():
            for t in TERMS:
                if t in ref and t in e["mosaics"]:
                    d = abs(e["mosaics"][t] - ref[t])
                    terms.append({"key": key, "engine": eng, "term": t, "d": d,
                                  "floor": floor_of.get(eng, floor)})
    for i, r in enumerate(terms):
        col = S.ENGINE_COLOUR[r["engine"]]
        low = r["d"] < r["floor"]
        # deterministic within-engine jitter so stacked markers stay separable
        wobble = (((i * 0.6180339887) % 1.0) - 0.5) * 0.16
        axb.plot(TERMS.index(r["term"]) + jitter[r["engine"]] + wobble,
                 max(r["d"], bottom * 1.6), marker=S.ENGINE_MARKER[r["engine"]],
                 ms=3.2, ls="none", zorder=3, alpha=0.9, clip_on=False,
                 mfc="none" if low else col, mec=col, mew=0.5)
    axb.set_xlim(-0.6, len(TERMS) - 0.4)
    axb.set_xticks(range(len(TERMS)))
    axb.set_xticklabels([TERM_LABEL[t] for t in TERMS], rotation=60, ha="right",
                        rotation_mode="anchor", fontsize=7)


    handles = [Line2D([], [], marker=S.ENGINE_MARKER[e], ls="none", ms=5,
                      mfc=S.ENGINE_COLOUR[e], mec=S.ENGINE_COLOUR[e], mew=0.6,
                      label=S.ENGINE_NAME[e]) for e in ("sander", "openmm", "gromacs")]
    handles.append(Line2D([], [], marker="o", ls="none", ms=5, mfc="none",
                          mec=S.TEXT, mew=0.6,
                          label="hollow: below the print resolution of the reference file"))
    fig.legend(handles=handles, loc="upper left", ncol=2, frameon=False,
               bbox_to_anchor=(0.02, 1.035), handletextpad=0.4,
               columnspacing=1.2, fontsize=6.5)
    return fig, totals, terms


def main():
    cells, files_read = load()
    worst = check(cells)
    fig, totals, terms = draw(cells)
    pdf, png = S.save(fig, OUT)
    plt.close(fig)

    ds = sorted(r["d"] for r in totals)
    within = {}
    for r in totals:
        within[r["key"]] = max(within.get(r["key"], 0.0), r["d"])
    n_band = sum(1 for v in within.values() if v <= BAND)
    lo_a = sum(1 for r in totals if r["d"] < r["floor"])
    lo_b = sum(1 for r in terms if r["d"] < r["floor"])
    split = {}
    for r in terms:
        split.setdefault(r["term"], 0)
        split[r["term"]] += 1
    unresolved = []
    for e in cells:
        for eng, (ref, _, path) in e["refs"].items():
            for t in ("lj", "elec"):
                if t not in ref:
                    unresolved.append((e["cell"]["system"], e["cell"]["force_field"], eng, t, path))

    print(f"wrote {pdf}\nwrote {png}")
    print(f"files read: {len(files_read)}")
    for f in files_read:
        print(f"   {f}")
    print(f"cells {len(cells)}   panel (a) markers {len(totals)}   panel (b) markers {len(terms)}")
    print(f"panel (b) markers per term: " + ", ".join(f"{t} {split[t]}" for t in TERMS))
    print(f"largest gap between a file-derived total and MATRIX.json: {worst:.3e}")
    print(f"|dE| over {len(ds)} comparisons: min {ds[0]:.3e}   median {statistics.median(ds):.3e}   max {ds[-1]:.3e}")
    print(f"cells with every available reference within {BAND}: {n_band} of {len(cells)}")
    for key, v in sorted(within.items(), key=lambda kv: -kv[1]):
        if v > BAND:
            print(f"   outside: {key[0]} {key[1]}  largest gap {v:.3e}")
    print(f"hollow markers: (a) {lo_a} of {len(totals)}   (b) {lo_b} of {len(terms)}")
    print(f"per-term LJ/electrostatics not resolvable: {len(unresolved)} engine-terms")
    for u in unresolved:
        print("   ", *u)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
