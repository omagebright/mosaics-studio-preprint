#!/usr/bin/env python3
"""How much of each MOSAICS deck the BPTI system actually touches.

Rows are matched the way MOSAICS matches them: bond/bend/vdw/1-4 rows by
atom-type tuple in either direction, torsion rows with exact rows taking
precedence over wildcard (X) rows.  The count reported is the number of
DISTINCT DECK ROWS this one structure exercises, over the number of rows the
deck contains.
"""
import re, sys, json
from pathlib import Path
from collections import Counter
import parmed as pmd

deckdir, prmtop, out = sys.argv[1], sys.argv[2], sys.argv[3]
tag = Path(deckdir).name.replace("db_", "")
D = Path(deckdir)
stem = "mosaics_openmm-%s" % tag

def blocks(path, opener):
    txt = Path(path).read_text()
    return re.findall(r"~%s\[(.*?)\]" % opener, txt, re.S)

def atoms_of(b, n):
    return tuple(re.search(r"\\atom%d\{([^}]*)\}" % i, b).group(1) for i in range(1, n + 1))

p = pmd.load_file(prmtop)
T = lambda a: "PROT-" + a.type

res = {}

# --- bonds
rows = [atoms_of(b, 2) for b in blocks(D / (stem + ".bond"), "bond_parm")]
used = {tuple(sorted((T(b.atom1), T(b.atom2)))) for b in p.bonds}
hit = {r for r in rows if tuple(sorted(r)) in used}
res["bond"] = (len(hit), len(rows))

# --- bends
rows = [atoms_of(b, 3) for b in blocks(D / (stem + ".bend"), "bend_parm")]
used = set()
for a in p.angles:
    k = (T(a.atom1), T(a.atom2), T(a.atom3))
    used.add(min(k, k[::-1]))
hit = {r for r in rows if min(r, r[::-1]) in used}
res["bend"] = (len(hit), len(rows))

# --- vdw and 1-4
for name, opener, pairs in (
        ("vdw", "inter_parm", {tuple(sorted((T(a), T(b))))
                               for i, a in enumerate(p.atoms) for b in p.atoms[i:]}),
        ("onfo", "onefour_parm", {tuple(sorted((T(d.atom1), T(d.atom4))))
                                  for d in p.dihedrals if not d.ignore_end and not d.improper})):
    rows = [atoms_of(b, 2) for b in blocks(D / (stem + "." + name), opener)]
    hit = {r for r in rows if tuple(sorted(r)) in pairs}
    res[name] = (len(hit), len(rows))

# --- torsions and impropers, with MOSAICS precedence
tb = blocks(D / (stem + ".tors_and_impr"), "torsion_parm")
prop, impr = [], []
for i, b in enumerate(tb):
    (impr if "\\label{improper}" in b else prop).append((i, atoms_of(b, 4)))

def match(rowset, quad, exact_only):
    out = []
    for i, r in rowset:
        for q in (quad, quad[::-1]):
            ok = all(rr == "X" or rr == qq for rr, qq in zip(r, q))
            if exact_only and "X" in r:
                ok = False
            if ok:
                out.append(i); break
    return out

for label, rowset, dihs in (
        ("torsion_proper", prop, [d for d in p.dihedrals if not d.improper]),
        ("torsion_improper", impr, [d for d in p.dihedrals if d.improper])):
    chosen = set()
    quads = {(T(d.atom1), T(d.atom2), T(d.atom3), T(d.atom4)) for d in dihs}
    for q in quads:
        m = match(rowset, q, exact_only=True) or match(rowset, q, exact_only=False)
        chosen.update(m)
    res[label] = (len(chosen), len(rowset))

# --- cmap
if p.cmaps:
    used = {c.type.idx for c in p.cmaps}
    total = len(blocks(D / (stem + ".cmap"), "cmap_parm")) if (D / (stem + ".cmap")).exists() else 0
    if total == 0:
        txt = (D / (stem + ".cmap")).read_text()
        total = len(re.findall(r"~cmap", txt))
    res["cmap"] = (len(used), total)

# --- composition
comp = Counter(r.name for r in p.residues)
meta = {
    "prmtop": prmtop,
    "atoms": len(p.atoms),
    "residues": len(p.residues),
    "net_charge": round(sum(a.charge for a in p.atoms), 6),
    "distinct_amber_atom_types": sorted({a.type for a in p.atoms}),
    "n_distinct_amber_atom_types": len({a.type for a in p.atoms}),
    "residue_counts": dict(sorted(comp.items())),
    "n_bonds": len(p.bonds), "n_angles": len(p.angles),
    "n_proper_dihedral_terms": sum(1 for d in p.dihedrals if not d.improper),
    "n_improper_dihedral_terms": sum(1 for d in p.dihedrals if d.improper),
    "n_one_four_pairs": sum(1 for d in p.dihedrals if not d.ignore_end and not d.improper),
    "n_cmap_terms": len(p.cmaps),
    "deck_rows_touched": {k: {"touched": v[0], "in_deck": v[1],
                              "percent": round(100.0 * v[0] / v[1], 1) if v[1] else None}
                          for k, v in res.items()},
}
Path(out).write_text(json.dumps(meta, indent=2) + "\n")
print(json.dumps(meta["deck_rows_touched"], indent=2))
print("atom types:", meta["n_distinct_amber_atom_types"])
