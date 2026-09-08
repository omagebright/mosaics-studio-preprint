#!/usr/bin/env python3
"""How much of a MOSAICS parameter deck a system actually touches.

The referee's objection is that the committed systems exercise an unstated
fraction of each deck.  This turns that into a stated fraction.

Method, and its limits, stated so the number can be argued with:

  * The system's bonds, angles, torsions and impropers are taken from the
    AMBER topology tleap built -- the same molecule MOSAICS reads, since the
    reference topology and the MOSAICS structure are built from one file.
  * AMBER atom types map onto MOSAICS deck types by the prefix the deck uses
    (`CX` -> `PROT-CX`).  That mapping is checked, not assumed: any type that
    fails to map is reported.
  * BOND and BEND deck rows are keyed on an exact type tuple, so the count of
    rows touched is exact.
  * TORSION rows may carry `X` wildcards, and where several rows match one
    dihedral the engine picks one.  Two numbers are therefore given: rows
    whose pattern is matched by at least one dihedral in the system (an upper
    bound on rows used) and wildcard-free rows matched exactly (a lower
    bound).
  * VDW / ONFO rows are keyed on a type pair, and every pair of types present
    interacts, so the count is the number of deck rows over the used types.
"""

import itertools
import re
import sys
from pathlib import Path

import parmed as pmd


def deck_rows(path, tag, natoms):
    """Read one deck file; return a list of atom-type tuples, one per row."""
    text = Path(path).read_text()
    rows = []
    for m in re.finditer(r"~%s\[(.*?)\]" % tag, text, re.S):
        body = m.group(1)
        atoms = []
        ok = True
        for i in range(1, natoms + 1):
            a = re.search(r"\\atom%d\{([^}]*)\}" % i, body)
            if a is None:
                ok = False
                break
            atoms.append(a.group(1))
        if ok:
            label = re.search(r"\\label\{([^}]*)\}", body)
            rows.append((tuple(atoms), label.group(1) if label else ""))
    return rows


def canon2(t):
    return tuple(sorted(t))


def canon3(t):
    return t if t[0] <= t[2] else (t[2], t[1], t[0])


def matches(pattern, actual):
    """Wildcard match of a type tuple, forward or reversed."""
    for cand in (actual, tuple(reversed(actual))):
        if all(p == "X" or p == c for p, c in zip(pattern, cand)):
            return True
    return False


def main():
    prmtop, deck_dir, prefix, tag = sys.argv[1:5]
    struct = pmd.load_file(prmtop)
    P = "PROT-"

    used_types = sorted({P + a.type for a in struct.atoms})

    bonds = {canon2((P + b.atom1.type, P + b.atom2.type))
             for b in struct.bonds}
    bends = {canon3((P + a.atom1.type, P + a.atom2.type, P + a.atom3.type))
             for a in struct.angles}

    propers, impropers = set(), set()
    for d in struct.dihedrals:
        q = (P + d.atom1.type, P + d.atom2.type,
             P + d.atom3.type, P + d.atom4.type)
        (impropers if d.improper else propers).add(q)

    d = Path(deck_dir)
    bond_rows = deck_rows(d / ("mosaics_%s.bond" % tag), "bond_parm", 2)
    bend_rows = deck_rows(d / ("mosaics_%s.bend" % tag), "bend_parm", 3)
    tors_rows = deck_rows(d / ("mosaics_%s.tors_and_impr" % tag),
                          "torsion_parm", 4)
    vdw_rows = deck_rows(d / ("mosaics_%s.vdw" % tag), "inter_parm", 2)
    onfo_rows = deck_rows(d / ("mosaics_%s.onfo" % tag), "onefour_parm", 2)

    bond_hit = {canon2(r) for r, _ in bond_rows} & bonds
    bend_hit = {canon3(r) for r, _ in bend_rows} & bends

    tors_deck = [(r, lab) for r, lab in tors_rows if lab != "improper"]
    impr_deck = [(r, lab) for r, lab in tors_rows if lab == "improper"]

    def hits(rows, actuals):
        upper, lower = 0, 0
        for r, _ in rows:
            hit = any(matches(r, a) for a in actuals)
            if hit:
                upper += 1
                if "X" not in r:
                    lower += 1
        return upper, lower

    tu, tl = hits(tors_deck, propers)
    iu, il = hits(impr_deck, impropers)

    pairs = {canon2(p) for p in itertools.combinations_with_replacement(
        used_types, 2)}
    vdw_hit = {canon2(r) for r, _ in vdw_rows} & pairs
    onfo_hit = {canon2(r) for r, _ in onfo_rows} & pairs

    print("topology            %s" % prmtop)
    print("atoms               %d" % len(struct.atoms))
    print("residues            %d" % len(struct.residues))
    print()
    print("term        deck rows   rows this system touches")
    print("-" * 62)
    print("atom types  %9d   %d" % (len({r[0] for r, _ in vdw_rows} |
                                        {r[1] for r, _ in vdw_rows}),
                                    len(used_types)))
    print("bond        %9d   %d" % (len(bond_rows), len(bond_hit)))
    print("bend        %9d   %d" % (len(bend_rows), len(bend_hit)))
    print("torsion     %9d   %d matched (%d of them wildcard-free)"
          % (len(tors_deck), tu, tl))
    print("improper    %9d   %d matched (%d of them wildcard-free)"
          % (len(impr_deck), iu, il))
    print("vdw pair    %9d   %d" % (len(vdw_rows), len(vdw_hit)))
    print("onfo pair   %9d   %d" % (len(onfo_rows), len(onfo_hit)))
    print()
    print("distinct type tuples present in the system")
    print("  bond %d   bend %d   proper %d   improper %d"
          % (len(bonds), len(bends), len(propers), len(impropers)))
    print("  bonds %d   angles %d   dihedral terms %d"
          % (len(struct.bonds), len(struct.angles), len(struct.dihedrals)))
    print()
    print("atom types used: %s" % ", ".join(t[len(P):] for t in used_types))


if __name__ == "__main__":
    main()
