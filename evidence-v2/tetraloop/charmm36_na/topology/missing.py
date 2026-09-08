#!/usr/bin/env python3
"""List the type tuples this structure needs that the deck has no row for."""
import re, sys, itertools
from collections import defaultdict
rtf_path, pdb_path, deck_stem = sys.argv[1], sys.argv[2], sys.argv[3]
res_atoms, res_bonds, cur = {}, {}, None
for raw in open(rtf_path):
    line = raw.split('!')[0].rstrip()
    if not line.strip(): continue
    tok = line.split()
    if tok[0] == 'RESI':
        cur = tok[1]; res_atoms[cur] = {}; res_bonds[cur] = []
    elif cur and tok[0] == 'ATOM': res_atoms[cur][tok[1]] = tok[2]
    elif cur and tok[0] in ('BOND', 'DOUBLE'):
        r = tok[1:]
        for i in range(0, len(r) - 1, 2): res_bonds[cur].append((r[i], r[i+1]))
residues, idx = [], 0
for line in open(pdb_path):
    if not line.startswith('ATOM'): continue
    name, rn, ri = line[12:16].strip(), line[17:20].strip(), int(line[22:26])
    if not residues or residues[-1][0] != ri: residues.append((ri, rn, {}))
    residues[-1][2][name] = idx; idx += 1
types, owner = {}, {}
for ri, rn, atoms in residues:
    for an, gi in atoms.items():
        types[gi] = res_atoms[rn][an]; owner[gi] = (ri, rn, an)
bonds = set()
for ri, rn, atoms in residues:
    for a, b in res_bonds[rn]:
        if a in atoms and b in atoms: bonds.add(frozenset((atoms[a], atoms[b])))
junction = {}
for (ri, rn, a), (rj, rm, b) in zip(residues, residues[1:]):
    if "O3'" in a and 'P' in b:
        p = frozenset((a["O3'"], b['P'])); bonds.add(p)
        junction[p] = f"{rn}{ri}->{rm}{rj}"
adj = defaultdict(set)
for bd in bonds:
    x, y = tuple(bd); adj[x].add(y); adj[y].add(x)
def rows(path, kind):
    out = []
    for chunk in open(path).read().split('~')[1:]:
        if not chunk.startswith(kind): continue
        ats = re.findall(r'\\atom(\d)\{([^}]*)\}', chunk)
        out.append(tuple(v for _, v in sorted(ats, key=lambda z: z[0])))
    return out
bond_rows = set()
for r in rows(deck_stem + '.bond', 'bond_parm'):
    bond_rows.add(r); bond_rows.add(r[::-1])
missing_bonds = defaultdict(list)
for bd in bonds:
    x, y = tuple(bd)
    if (types[x], types[y]) not in bond_rows:
        lbl = junction.get(bd, 'intra-residue')
        missing_bonds[(types[x], types[y], lbl)].append((owner[x], owner[y]))
print(f"BOND rows missing: {len(missing_bonds)} distinct, "
      f"{sum(len(v) for v in missing_bonds.values())} instances of {len(bonds)} bonds")
for (t1, t2, lbl), inst in sorted(missing_bonds.items()):
    print(f"  {t1:<12} {t2:<12} {lbl:<16} x{len(inst)}   e.g. {inst[0][0]} - {inst[0][1]}")
bend_rows = set()
for r in rows(deck_stem + '.bend', 'bend_parm'):
    bend_rows.add(r); bend_rows.add(r[::-1])
angles, missing_ang = [], defaultdict(int)
for c in range(idx):
    for a, b in itertools.combinations(sorted(adj[c]), 2):
        angles.append((a, c, b))
        if (types[a], types[c], types[b]) not in bend_rows:
            missing_ang[(types[a], types[c], types[b])] += 1
print(f"BEND rows missing: {len(missing_ang)} distinct, "
      f"{sum(missing_ang.values())} instances of {len(angles)} bends")
for k in sorted(missing_ang)[:40]:
    print(f"  {k[0]:<12} {k[1]:<12} {k[2]:<12} x{missing_ang[k]}")
