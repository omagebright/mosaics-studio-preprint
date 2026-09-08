#!/usr/bin/env python3
"""Put chain IDs back on a tleap-saved PDB.

tleap's savepdb writes a blank chain column.  1BNA is a two-chain duplex and
MOSAICS needs the two strands distinguished, so residues 1-12 are labelled A
and 13-24 are labelled B -- exactly the chain assignment of the cleaned
deposit the topology was built from.  Nothing else in the file is touched.
"""
import sys
src, dst = sys.argv[1], sys.argv[2]
out = []
for line in open(src):
    if line.startswith(("ATOM", "HETATM")):
        resseq = int(line[22:26])
        ch = "A" if resseq <= 12 else "B"
        line = line[:21] + ch + line[22:]
    out.append(line)
open(dst, "w").writelines(out)
