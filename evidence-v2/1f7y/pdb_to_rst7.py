#!/usr/bin/env python3
"""Write an AMBER restart whose coordinates are exactly the MOSAICS PDB's.

A PDB carries three decimals; the tleap restart carries seven. Two engines
handed different files evaluate different structures (up to 5e-4 A here,
which is worth 0.6 kcal/mol of Lennard-Jones on 3321 atoms), so every
reference engine in this tree reads a restart written FROM the PDB. Atoms are
matched by (residue number, atom name) after undoing the MOSAICS renames,
never by position.

    python3 pdb_to_rst7.py <prmtop> <mosaics.pdb> <out.rst7>
"""
import sys
import numpy as np
import parmed as pmd

BACK = {"O1P": "OP1", "O2P": "OP2", "H2''": "H2'", "H2'": "HO2'"}
RNA = {"RA", "RC", "RG", "RU", "RA5", "RC5", "RG5", "RU5", "RA3", "RC3", "RG3", "RU3"}

prmtop, pdb, out = sys.argv[1:4]
parm = pmd.load_file(prmtop)
pos = {}
for line in open(pdb):
    if not line.startswith("ATOM"):
        continue
    res = line[17:20].strip()
    name = line[12:16].strip()
    if res in RNA and name in BACK:
        name = BACK[name]
    key = (int(line[22:26]), name)
    if key in pos:
        sys.exit(f"duplicate atom {key}")
    pos[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
if len(pos) != len(parm.atoms):
    sys.exit(f"{len(pos)} PDB atoms versus {len(parm.atoms)} prmtop atoms")
xyz = np.zeros((len(parm.atoms), 3))
for atom in parm.atoms:
    xyz[atom.idx] = pos[(atom.residue.idx + 1, atom.name)]
parm.coordinates = xyz
parm.save(out, format="rst7", overwrite=True)
check = pmd.load_file(prmtop, out)
worst = float(np.abs(check.coordinates - xyz).max())
if worst != 0.0:
    sys.exit(f"restart does not reproduce the PDB: {worst} A")
print(f"{out}: {len(parm.atoms)} atoms, round trip 0.0 A")
