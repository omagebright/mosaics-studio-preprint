#!/usr/bin/env python3
"""How much of a MOSAICS deck one structure actually touches.

Reads the RTF to get each residue's atom -> MOSAICS-type map and its
intra-residue bonds, walks the PDB, adds the O3'(i)-P(i+1) backbone links,
and enumerates the type tuples the structure's bonds, bends, 1-4 pairs and
non-bonded pairs present.  Then counts how many rows of each deck file at
least one of those tuples matches.

Torsion rows are matched on the four type slots with X as a wildcard.
\\mol_type and \\res_context are NOT used in the match -- the engine's own
selection rule for them is not reimplemented here -- so the torsion figure
is an upper bound on rows touched, and is reported as such.
"""
import re, sys, itertools
from collections import defaultdict

rtf_path, pdb_path, deck_stem = sys.argv[1], sys.argv[2], sys.argv[3]

# ---- RTF ----
res_atoms = {}      # resi -> {atomname: type}
res_bonds = {}      # resi -> [(a,b)]
cur = None
for raw in open(rtf_path):
    line = raw.split('!')[0].rstrip()
    if not line.strip():
        continue
    tok = line.split()
    if tok[0] == 'RESI':
        cur = tok[1]
        res_atoms[cur] = {}
        res_bonds[cur] = []
    elif cur and tok[0] == 'ATOM':
        res_atoms[cur][tok[1]] = tok[2]
    elif cur and tok[0] in ('BOND', 'DOUBLE'):
        rest = tok[1:]
        for i in range(0, len(rest) - 1, 2):
            res_bonds[cur].append((rest[i], rest[i + 1]))

# ---- PDB ----
residues = []       # list of (resname, {atomname: global index})
order = []
idx = 0
for line in open(pdb_path):
    if not line.startswith('ATOM'):
        continue
    name = line[12:16].strip(); rn = line[17:20].strip(); ri = int(line[22:26])
    if not residues or residues[-1][0] != ri:
        residues.append((ri, rn, {}))
    residues[-1][2][name] = idx
    order.append((rn, name))
    idx += 1

types = {}
for ri, rn, atoms in residues:
    if rn not in res_atoms:
        sys.exit(f"residue {rn} not in RTF")
    for an, gi in atoms.items():
        if an not in res_atoms[rn]:
            sys.exit(f"atom {rn}-{an} not in RTF")
        types[gi] = res_atoms[rn][an]

bonds = set()
for ri, rn, atoms in residues:
    for a, b in res_bonds[rn]:
        if a in atoms and b in atoms:
            bonds.add(frozenset((atoms[a], atoms[b])))
links = 0
for (ri, rn, a), (rj, rm, b) in zip(residues, residues[1:]):
    if "O3'" in a and 'P' in b:
        bonds.add(frozenset((a["O3'"], b['P']))); links += 1

adj = defaultdict(set)
for bd in bonds:
    x, y = tuple(bd)
    adj[x].add(y); adj[y].add(x)

n = idx
bond_pairs = {tuple(sorted((types[x], types[y]))) for x, y in (tuple(b) for b in bonds)}

angles = set(); angle_tris = set()
for c in range(n):
    for a, b in itertools.combinations(sorted(adj[c]), 2):
        angles.add((a, c, b))
        t = (types[a], types[c], types[b])
        angle_tris.add(min(t, t[::-1]))

torsions = set(); tors_quads = set()
for b, c in (tuple(x) for x in bonds):
    for a in adj[b] - {c}:
        for d in adj[c] - {b}:
            if a == d:
                continue
            torsions.add((a, b, c, d))
            q = (types[a], types[b], types[c], types[d])
            tors_quads.add(min(q, q[::-1]))

# 1-4 pairs: ends of a torsion that are not 1-2 or 1-3
onethree = {frozenset((a, b)) for a, c, b in angles}
onefour = set()
for a, b, c, d in torsions:
    p = frozenset((a, d))
    if p not in bonds and p not in onethree:
        onefour.add(p)
onefour_pairs = {tuple(sorted((types[x], types[y]))) for x, y in (tuple(p) for p in onefour)}

nb_pairs = set()
alltypes = sorted(set(types.values()))
for t1, t2 in itertools.combinations_with_replacement(alltypes, 2):
    nb_pairs.add(tuple(sorted((t1, t2))))

def deck_rows(path, kind):
    txt = open(path).read()
    rows = []
    for chunk in txt.split('~')[1:]:
        if not chunk.startswith(kind):
            continue
        ats = re.findall(r'\\atom(\d)\{([^}]*)\}', chunk)
        slots = [v for _, v in sorted(ats, key=lambda z: z[0])]
        rows.append((slots, chunk))
    return rows

def report(label, path, kind, touched, symmetric=True):
    rows = deck_rows(path, kind)
    hit = 0
    tset = set(touched)
    for slots, chunk in rows:
        key = tuple(slots)
        if key in tset or (symmetric and key[::-1] in tset):
            hit += 1
    print(f"{label:<26} deck rows {len(rows):>6}   touched {hit:>6}   "
          f"({100.0*hit/max(len(rows),1):5.1f}%)   distinct tuples in system {len(tset)}")
    return len(rows), hit, len(tset)

print(f"atoms {n}   residues {len(residues)}   distinct MOSAICS atom types {len(alltypes)}")
print(f"bonds {len(bonds)} (of which {links} inter-residue O3'-P links)   "
      f"bends {len(angles)}   torsions(1-4 paths) {len(torsions)}   1-4 pairs {len(onefour)}")
print()
report('bond', deck_stem + '.bond', 'bond_parm', bond_pairs)
report('bend', deck_stem + '.bend', 'bend_parm', angle_tris)
report('onefour (1-4 LJ)', deck_stem + '.onfo', 'onefour_parm', onefour_pairs)
report('inter (LJ)', deck_stem + '.vdw', 'inter_parm', nb_pairs)

# torsions, with X wildcards, ignoring mol_type / res_context
rows = deck_rows(deck_stem + '.tors_and_impr', 'torsion_parm')
def matches(slots, q):
    return all(s == 'X' or s == t for s, t in zip(slots, q))
hit = 0; hit_impr = 0; tot_impr = 0
for slots, chunk in rows:
    impr = '\\label{improper}' in chunk
    tot_impr += impr
    ok = any(matches(slots, q) or matches(slots, q[::-1]) for q in tors_quads)
    hit += ok
    hit_impr += ok and impr
print(f"{'torsion+improper':<26} deck rows {len(rows):>6}   touched {hit:>6}   "
      f"({100.0*hit/len(rows):5.1f}%)   distinct quadruplets in system {len(tors_quads)}")
print(f"{'  of which improper rows':<26} deck rows {tot_impr:>6}   touched {hit_impr:>6}")
molt = defaultdict(int); molt_hit = defaultdict(int)
for slots, chunk in rows:
    m = re.search(r'\\mol_type\{(\w+)\}', chunk)
    tag = m.group(1) if m else 'none'
    molt[tag] += 1
    if any(matches(slots, q) or matches(slots, q[::-1]) for q in tors_quads):
        molt_hit[tag] += 1
for tag in sorted(molt):
    print(f"    mol_type {tag:<6} rows {molt[tag]:>6}   touched {molt_hit[tag]:>6}")
print()
print("torsion figure ignores \\mol_type and \\res_context, so it is an upper bound")
