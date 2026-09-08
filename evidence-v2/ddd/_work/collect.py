#!/usr/bin/env python3
"""Read every engine's committed output file for ddd and print the tables.

Nothing is recomputed from another program's numbers; each value below is
parsed out of the file named beside it.
"""
import re, json
from pathlib import Path

BASE = Path("/Users/bright/Library/CloudStorage/OneDrive-Personal/Documents/MOSAICS_AllAtom_Paper/evidence-v2/ddd")
KJ = 4.184
TARGET = 332.0637133
K = {"sander": 332.052200, "openmm": 332.063713, "gromacs": 332.063713, "mosaics": 332.063441}

def sander_terms(p):
    blk = Path(p).read_text(errors="replace").split("Maximum number of minimization cycles reached")[0]
    pat = re.compile(r"(BOND|ANGLE|DIHED|VDWAALS|EEL|HBOND|1-4 VDW|1-4 EEL|RESTRAINT)\s*=\s*([-+0-9.Ee]+)")
    v = {k.strip(): float(x) for k, x in pat.findall(blk)}
    return {"bond": v["BOND"], "angle": v["ANGLE"], "torsion": v["DIHED"],
            "proper": None, "improper": None,
            "lj": v["VDWAALS"] + v["1-4 VDW"], "coul": v["EEL"] + v["1-4 EEL"]}

def openmm_terms(p):
    t = Path(p).read_text()
    g = lambda k: float(re.search(re.escape(k) + r"\s+([-0-9.]+)", t).group(1))
    return {"bond": g("HarmonicBondForce"), "angle": g("HarmonicAngleForce"),
            "torsion": g("PeriodicTorsionForce"), "proper": None, "improper": None,
            "lj": g("Lennard-Jones (charges zeroed)"), "coul": g("Coulomb (epsilons zeroed)")}

def gromacs_terms(p):
    lines = Path(p).read_text().splitlines()
    leg = {}
    for l in lines:
        m = re.match(r'@ s(\d+) legend "(.*)"', l)
        if m: leg[int(m.group(1))] = m.group(2)
    row = [float(x) for x in lines[-1].split()]
    v = {leg[i]: row[i + 1] / KJ for i in range(len(leg))}
    return {"bond": v["Bond"], "angle": v["Angle"],
            "torsion": v["Proper Dih."] + v["Per. Imp. Dih."],
            "proper": v["Proper Dih."], "improper": v["Per. Imp. Dih."],
            "lj": v["LJ (SR)"] + v["LJ-14"], "coul": v["Coulomb (SR)"] + v["Coulomb-14"]}

def mosaics_terms(p):
    t = Path(p).read_text()
    blk = t.split(">>Report energies at step 1<<")[1]
    g = lambda k: float(re.search(re.escape(k) + r"\s+([-0-9.e+]+)", blk).group(1))
    return {"bond": g("-Bond energy"), "angle": g("Bend-Bend energy"),
            "ureybradley": g("Bond-Bend energy"),
            "torsion": g("-Torsion energy"), "proper": g("Tors-Tors energy"),
            "improper": g("Tors-Improper energy"),
            "lj14": g("Onfo-Lennard-Jones energy"),
            "lj": g("-Total Lennard-Jones energy"), "coul": g("-Total Coulomb energy")}

def corr(v, eng):
    return None if v is None else v * (TARGET / K[eng])

FF = [("ol21_ol3", "OL21/OL3"), ("ol15_ol3", "OL15/OL3"), ("bsc1", "parmbsc1"), ("bs0", "parmbsc0")]
out = []
for key, name in FF:
    d = BASE / key
    row = {"ff": name, "key": key}
    row["sander"] = sander_terms(d / "sander" / f"{key}_sander_sp.out")
    row["openmm"] = openmm_terms(d / "openmm" / "openmm_sp.out")
    row["gromacs"] = gromacs_terms(d / "gromacs" / "sp.xvg")
    for build, tag in (("mosaics-current", "current"), ("mosaics-3.9.1", "stable")):
        f = d / build / "mosaics_sp.out"
        try:
            row[tag] = mosaics_terms(f)
        except Exception:
            row[tag] = None
    out.append(row)

TERMS = [("bond", "bond"), ("angle", "angle"), ("torsion", "torsion (proper+improper)"),
         ("proper", "  proper only"), ("improper", "  improper only"),
         ("lj", "lennard-jones (total)"), ("coul", "electrostatics (as printed)")]

for r in out:
    print("=" * 92)
    print(f"{r['ff']} on ddd (1BNA)")
    print("=" * 92)
    hdr = f"{'term':<28}{'MOSAICS 3.9.1':>17}{'MOSAICS current':>19}{'sander':>15}{'OpenMM':>19}{'GROMACS':>19}"
    print(hdr)
    for k, lab in TERMS:
        def f(src, eng=None):
            if src is None: return "n/a"
            v = src.get(k)
            return "-" if v is None else f"{v:.8f}"
        print(f"{lab:<28}{f(r['stable']):>17}{f(r['current']):>19}{f(r['sander']):>15}{f(r['openmm']):>19}{f(r['gromacs']):>19}")
    # corrected electrostatics
    def c(src, eng):
        if src is None or src.get("coul") is None: return "n/a"
        return f"{corr(src['coul'], eng):.8f}"
    print(f"{'electrostatics (corrected)':<28}{c(r['stable'],'mosaics'):>17}{c(r['current'],'mosaics'):>19}{c(r['sander'],'sander'):>15}{c(r['openmm'],'openmm'):>19}{c(r['gromacs'],'gromacs'):>19}")
    print()
    for tag, label in (("current", "CURRENT"), ("stable", "3.9.1")):
        s = r[tag]
        if s is None:
            print(f"  MOSAICS {label}: did not run (see RESULTS.md)"); continue
        print(f"  MOSAICS {label} minus OpenMM, electrostatics after correction:")
        o = r["openmm"]; g = r["gromacs"]
        for k, lab in (("bond","bond"),("angle","angle"),("torsion","torsion"),("lj","lj")):
            print(f"    {lab:<10} {s[k]-o[k]:+.4e}")
        print(f"    {'coulomb':<10} {corr(s['coul'],'mosaics')-corr(o['coul'],'openmm'):+.4e}")
        print(f"  MOSAICS {label} minus GROMACS, split torsion:")
        print(f"    {'proper':<10} {s['proper']-g['proper']:+.4e}")
        print(f"    {'improper':<10} {s['improper']-g['improper']:+.4e}")
        print()
