#!/usr/bin/env python3
"""Build the BPTI (5PTI) AMBER reference topology for a MOSAICS comparison.

Same contract as tools/build_protein_reference.py: the structure MOSAICS reads
is the one the reference topology is built from, and the restart file the
reference engines read is written FROM that PDB so both sides evaluate
identical coordinates to the last bit.

Difference from the committed builder: BPTI is a whole protein read from a PDB,
not a sequence, and it carries three disulfides that must be made explicitly in
leap.  Termini are charged (NH3+/COO-) per BUILD-DECISIONS.md.
"""
import json, subprocess, sys
from pathlib import Path
import numpy as np
import parmed as pmd

DISULFIDES = [(5, 55), (14, 38), (30, 51)]


def read_pdb_atoms(path):
    out = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out.append((line[21], line[22:27].strip(),
                        line[17:20].strip(), line[12:16].strip(), line))
    return out


def pdb_coordinate_map(pdb):
    order, seen, out = 0, None, {}
    for chain, resseq, resname, atom, line in read_pdb_atoms(pdb):
        if resseq != seen:
            seen = resseq
            order += 1
        out[(order, resname, atom)] = (float(line[30:38]), float(line[38:46]),
                                       float(line[46:54]))
    return out


def main(name, leaprc, outdir, struct):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    src = Path(struct)
    (out / src.name).write_text(src.read_text())

    leap_in = out / ("build_%s.in" % name)
    bonds = "\n".join("bond mol.%d.SG mol.%d.SG" % p for p in DISULFIDES)
    leap_in.write_text(
        "source %s\n"
        "mol = loadpdb %s\n"
        "%s\n"
        "saveamberparm mol %s.prmtop %s.rst7\n"
        "savepdb mol %s_leap.pdb\n"
        "quit\n" % (leaprc, src.name, bonds, name, name, name))
    proc = subprocess.run(["tleap", "-f", leap_in.name], cwd=out,
                          capture_output=True, text=True)
    (out / ("build_%s.log" % name)).write_text(
        proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0:
        print(proc.stdout[-4000:]); raise SystemExit("tleap failed")

    leap_pdb = out / ("%s_leap.pdb" % name)
    # MOSAICS reads the leap PDB with a chain id stamped on; no atom is renamed
    # (BPTI has no ACE/NME caps, and every standard residue already matches).
    mos = out / ("%s_mosaics.pdb" % name)
    lines = []
    for line in leap_pdb.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            line = line[:21] + "A" + line[22:]
        lines.append(line)
    mos.write_text("".join(lines))

    struct_t = pmd.load_file(str(out / ("%s.prmtop" % name)))
    lookup = pdb_coordinate_map(leap_pdb)
    if len(lookup) != len(struct_t.atoms):
        raise SystemExit("atom count mismatch %d vs %d"
                         % (len(struct_t.atoms), len(lookup)))
    coords = np.zeros((len(struct_t.atoms), 3), dtype=np.float64)
    for i, a in enumerate(struct_t.atoms):
        coords[i] = lookup[(a.residue.idx + 1, a.residue.name, a.name)]
    struct_t.coordinates = coords
    matched = out / ("%s_matched.rst7" % name)
    struct_t.save(str(matched), format="rst7", overwrite=True)

    check = pmd.load_file(str(out / ("%s.prmtop" % name)), xyz=str(matched))
    mlook = pdb_coordinate_map(mos)
    worst = 0.0
    for i, a in enumerate(check.atoms):
        want = mlook[(a.residue.idx + 1, a.residue.name, a.name)]
        got = check.coordinates[i]
        worst = max(worst, max(abs(got[k] - want[k]) for k in range(3)))
    if worst > 1e-9:
        raise SystemExit("coordinates disagree by %.3e A" % worst)

    # Confirm the three disulfides really are in the topology.
    ss = []
    for b in check.bonds:
        if b.atom1.name == "SG" and b.atom2.name == "SG":
            r1, r2 = b.atom1.residue, b.atom2.residue
            d = np.linalg.norm(check.coordinates[b.atom1.idx]
                               - check.coordinates[b.atom2.idx])
            ss.append({"res1": "%s%d" % (r1.name, r1.idx + 1),
                       "res2": "%s%d" % (r2.name, r2.idx + 1),
                       "distance_angstrom": round(float(d), 4)})
    if len(ss) != 3:
        raise SystemExit("expected 3 SG-SG bonds, topology has %d" % len(ss))

    resnames = [r.name for r in check.residues]
    meta = {"name": name, "leaprc": leaprc, "structure": src.name,
            "atoms": len(check.atoms), "residues": len(check.residues),
            "net_charge": round(sum(a.charge for a in check.atoms), 6),
            "first_residue": resnames[0], "last_residue": resnames[-1],
            "disulfides": ss,
            "max_coordinate_disagreement_angstrom": worst,
            "files": {"topology": "%s.prmtop" % name,
                      "reference_coordinates": "%s_matched.rst7" % name,
                      "mosaics_structure": "%s_mosaics.pdb" % name}}
    (out / "build_metadata.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:])
