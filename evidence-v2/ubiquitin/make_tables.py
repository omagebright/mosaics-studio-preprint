#!/usr/bin/env python3
"""Assemble the ubiquitin energy tables from the raw engine outputs.

Reads only committed output files.  Nothing is typed in by hand.
"""
import re, sys
from pathlib import Path

B = Path("/Users/bright/Library/CloudStorage/OneDrive-Personal/Documents/MOSAICS_AllAtom_Paper/evidence-v2/ubiquitin")
KJ = 4.184
TARGET   = 332.0637133          # CODATA 2018, evidence/engines.txt
CONST = {"sander": 332.052200, "openmm": 332.063713,
         "gromacs": 332.063713, "mosaics": 332.063441}

def sander(ff):
    t = (B/ff/"sander"/"ubiquitin_sander_sp.out").read_text().split("FINAL RESULTS")[1]
    g = lambda k: float(re.search(re.escape(k)+r"\s*=\s*(-?[\d.]+)", t).group(1))
    out = {"bond": g("BOND"), "angle": g("ANGLE"), "dihedral": g("DIHED"),
           "lj": g("VDWAALS")+g("1-4 VDW"), "coulomb": g("EEL")+g("1-4 EEL")}
    if "CMAP" in t: out["cmap"] = g("CMAP")
    return out

def openmm(ff):
    t = (B/ff/"openmm"/"openmm_sp.out").read_text()
    g = lambda k: float(re.search(re.escape(k)+r"\s+(-?[\d.]+)", t).group(1))
    out = {"bond": g("HarmonicBondForce"), "angle": g("HarmonicAngleForce"),
           "dihedral": g("PeriodicTorsionForce"),
           "lj": g("Lennard-Jones (charges zeroed)"),
           "coulomb": g("Coulomb (epsilons zeroed)")}
    if "CMAPTorsionForce" in t: out["cmap"] = g("CMAPTorsionForce")
    return out

def gromacs(ff):
    lines = (B/ff/"gromacs"/"sp.xvg").read_text().splitlines()
    leg = [re.search(r'legend "(.*)"', l).group(1) for l in lines if l.startswith("@ s")]
    vals = [float(x) for x in [l for l in lines if l and l[0] not in "#@"][0].split()][1:]
    d = dict(zip(leg, vals))
    out = {"bond": d["Bond"]/KJ, "angle": d["Angle"]/KJ,
           "dihedral": (d["Proper Dih."]+d["Per. Imp. Dih."])/KJ,
           "lj": (d["LJ-14"]+d["LJ (SR)"])/KJ,
           "coulomb": (d["Coulomb-14"]+d["Coulomb (SR)"])/KJ}
    if "CMAP Dih." in d: out["cmap"] = d["CMAP Dih."]/KJ
    return out

def mosaics(path):
    t = Path(path).read_text()
    if "Report energies at step 1" not in t: return None
    b = t.split("Report energies at step 1")[1]
    g = lambda k: float(re.search(re.escape(k)+r"\s+(-?[\d.e+-]+)", b).group(1))
    out = {"bond": g("-Bond energy"), "angle": g("-Bend energy"),
           "dihedral": g("-Torsion energy"),
           "lj": g("-Total Lennard-Jones energy"),
           "coulomb": g("-Total Coulomb energy")}
    if "-Cmap energy" in b: out["cmap"] = g("-Cmap energy")
    return out

TERMS = ["bond", "angle", "dihedral", "cmap", "lj", "coulomb"]

for ff, label in (("ff14sb", "ff14SB"), ("ff19sb", "ff19SB")):
    cur = mosaics(B/ff/"mosaics-current"/"mosaics_sp.out")
    sta = mosaics(B/ff/"mosaics-3.9.1"/"mosaics_sp.out")
    sa, om, gm = sander(ff), openmm(ff), gromacs(ff)
    print("%s on ubiquitin (1UBQ)" % label)
    print("-" * len("%s on ubiquitin (1UBQ)" % label))
    print("  term          MOSAICS 3.9.1        MOSAICS current           sander"
          "           OpenMM          GROMACS")
    for k in TERMS:
        if k not in om: continue
        s = "  %-12s" % k
        s += "%20s" % ("no output" if sta is None else "%.8f" % sta[k])
        s += " %18.8f" % cur[k]
        s += " %16.8f" % sa[k]
        s += " %16.8f" % om[k]
        s += " %16.8f" % (gm[k] if k in gm else float("nan"))
        print(s.replace("nan", "  -  "))
    print("             (electrostatic values above are as printed)")
    corr = lambda v, e: v * TARGET / CONST[e]
    row = "  coulomb*    "
    row += "%20s" % ("no output" if sta is None else "%.8f" % corr(sta["coulomb"], "mosaics"))
    row += " %18.8f" % corr(cur["coulomb"], "mosaics")
    row += " %16.8f" % corr(sa["coulomb"], "sander")
    row += " %16.8f" % corr(om["coulomb"], "openmm")
    row += " %16.8f" % corr(gm["coulomb"], "gromacs")
    print(row + "   corrected to one constant")
    def total(d, e, has_gmx_cmap=True):
        t = d["bond"] + d["angle"] + d["dihedral"] + d["lj"] + corr(d["coulomb"], e)
        if "cmap" in d: t += d["cmap"]
        return t
    tr = "  total       "
    tr += "%20s" % ("no output" if sta is None else "%.8f" % total(sta, "mosaics"))
    tr += " %18.8f" % total(cur, "mosaics")
    tr += " %16.8f" % total(sa, "sander")
    tr += " %16.8f" % total(om, "openmm")
    if ff == "ff19sb":
        tr += "%16s" % "-"
    else:
        tr += " %16.8f" % total(gm, "gromacs")
    print(tr)
    print()
    print("  MOSAICS minus OpenMM, after correction:")
    for k in TERMS:
        if k not in om: continue
        v = corr(cur[k], "mosaics") - corr(om[k], "openmm") if k == "coulomb" else cur[k] - om[k]
        print("    current  %-9s %+.3e" % (k, v))
    print("    3.9.1    every term   no output: the run stops before any energy")
    print()
    print("  GROMACS minus OpenMM:")
    for k in TERMS:
        if k not in om or k not in gm: continue
        v = corr(gm[k], "gromacs") - corr(om[k], "openmm") if k == "coulomb" else gm[k] - om[k]
        print("    %-11s %+.3e" % (k, v))
    print()
    print("  sander minus OpenMM:")
    for k in TERMS:
        if k not in om: continue
        v = corr(sa[k], "sander") - corr(om[k], "openmm") if k == "coulomb" else sa[k] - om[k]
        print("    %-11s %+.3e" % (k, v))
    print()
