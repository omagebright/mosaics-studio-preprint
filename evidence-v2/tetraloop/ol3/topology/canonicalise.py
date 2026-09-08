#!/usr/bin/env python3
"""Force the reference coordinates to be exactly the PDB the engines read.

tleap writes the rst7 to 1e-7 A and the PDB to 1e-3 A.  MOSAICS reads the
PDB.  Feeding the rst7 to sander/OpenMM/GROMACS and the PDB to MOSAICS puts
a 5e-4 A geometry difference into the comparison, which is worth ~1e-3
kcal/mol on a stiff bond -- larger than everything being measured.

The committed rule is that the structure MOSAICS reads is the one the
reference topology is built from.  This rewrites the rst7 from the PDB so
all four engines see one set of coordinates.
"""
import sys
import numpy as np
import parmed as pmd

prmtop, rst7_in, pdb, rst7_out = sys.argv[1:5]
p = pmd.load_file(prmtop, xyz=rst7_in)
xyz, names = [], []
for line in open(pdb):
    if line.startswith(("ATOM", "HETATM")):
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        names.append(line[12:16].strip())
xyz = np.array(xyz)
assert len(xyz) == len(p.atoms), (len(xyz), len(p.atoms))
assert names == [a.name for a in p.atoms], "atom order differs between prmtop and pdb"
print("max |rst7 - pdb| before =", np.abs(np.asarray(p.coordinates) - xyz).max())
p.coordinates = xyz
p.save(rst7_out, overwrite=True)
q = pmd.load_file(prmtop, xyz=rst7_out)
print("max |rst7 - pdb| after  =", np.abs(np.asarray(q.coordinates) - xyz).max())
