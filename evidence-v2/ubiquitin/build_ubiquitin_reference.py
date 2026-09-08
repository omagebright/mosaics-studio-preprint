#!/usr/bin/env python3
"""Build a matched ubiquitin (1UBQ) reference topology for a MOSAICS comparison.

Same approach as tools/build_protein_reference.py in the CURRENT worktree, with
the one change 1UBQ forces: the sequence is not typed out, it is read from a
cleaned PDB.  Everything else is carried over unchanged --

  * tleap builds the topology and writes its own PDB;
  * the tleap PDB is CANONICAL and the restart file the reference engines read
    is written FROM it, so MOSAICS and the reference engines evaluate the same
    coordinates to the last bit;
  * the overlay is asserted, not assumed.

Chemistry follows evidence-v2/_structures/BUILD-DECISIONS.md exactly:
  HIS 68 -> HIE, charged NH3+/COO- termini (no caps), no waters, no ions,
  hydrogens rebuilt by tleap.

There is no cap rename here: 1UBQ has no ACE/NME.  Every residue and atom name
tleap writes is already the name the ff14SB / ff19SB RTF uses; that is checked
against the RTF rather than asserted.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run(argv, cwd, log_path):
    proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True)
    log_path.write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError("%s failed, see %s" % (argv[0], log_path))
    return proc.stdout


def read_pdb_atoms(path):
    out = []
    for line in path.read_text().splitlines():
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
        out[(order, resname, atom)] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def write_mosaics_pdb(source, dest, chain_id):
    """Stamp a chain id.  No atom or residue name is changed."""
    lines = []
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            line = line[:21] + chain_id + line[22:]
        lines.append(line)
    dest.write_text("".join(lines))


def rtf_templates(rtf_path):
    """RESI name -> set of ATOM names, read from a MOSAICS RTF."""
    out, current = {}, None
    for line in Path(rtf_path).read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "RESI":
            current = parts[1]
            out[current] = set()
        elif parts[0] == "ATOM" and current is not None:
            out[current].add(parts[1])
    return out


def check_against_rtf(mosaics_pdb, rtf_path):
    """Every residue must resolve to an RTF template holding all its atoms.

    MOSAICS picks the NXXX template from backbone H2/H3 and the CXXX template
    from OXT (params/FORCEFIELD_USAGE.md section 4), so the same rule is
    applied here to decide which template a residue must match.
    """
    templates = rtf_templates(rtf_path)
    residues, order, seen = [], 0, None
    for chain, resseq, resname, atom, line in read_pdb_atoms(mosaics_pdb):
        if resseq != seen:
            seen = resseq
            order += 1
            residues.append([resname, resseq, set()])
        residues[-1][2].add(atom)

    problems, resolved = [], []
    for index, (resname, resseq, atoms) in enumerate(residues):
        name = resname
        if {"H2", "H3"} <= atoms:
            name = "N" + resname
        elif "OXT" in atoms:
            name = "C" + resname
        if name not in templates:
            problems.append("residue %s %s -> no RTF template %s"
                            % (resname, resseq, name))
            continue
        missing = atoms - templates[name]
        if missing:
            problems.append("residue %s %s (template %s): atoms not in "
                            "template: %s" % (resname, resseq, name,
                                              sorted(missing)))
        resolved.append(name)
    return resolved, problems


def overlay_coordinates(prmtop, leap_pdb, mosaics_pdb, rst7_out):
    import numpy as np
    import parmed as pmd

    struct = pmd.load_file(str(prmtop))
    lookup = pdb_coordinate_map(leap_pdb)
    if len(lookup) != len(struct.atoms):
        raise RuntimeError("atom count mismatch: topology %d, pdb %d"
                           % (len(struct.atoms), len(lookup)))

    coords = np.zeros((len(struct.atoms), 3), dtype=np.float64)
    unmatched = []
    for index, atom in enumerate(struct.atoms):
        key = (atom.residue.idx + 1, atom.residue.name, atom.name)
        if key not in lookup:
            unmatched.append(key)
            continue
        coords[index] = lookup[key]
    if unmatched:
        raise RuntimeError("%d topology atoms have no match in %s, first few: %s"
                           % (len(unmatched), leap_pdb.name, unmatched[:5]))

    renamed = pdb_coordinate_map(mosaics_pdb)
    if len(renamed) != len(lookup):
        raise RuntimeError("the MOSAICS structure has a different atom count")
    if sorted(lookup.values()) != sorted(renamed.values()):
        raise RuntimeError("writing the MOSAICS structure changed a coordinate")

    struct.coordinates = coords
    struct.save(str(rst7_out), format="rst7", overwrite=True)
    return len(struct.atoms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, type=Path)
    ap.add_argument("--name", required=True)
    ap.add_argument("--leaprc", default="leaprc.protein.ff14SB")
    ap.add_argument("--rtf", required=True, type=Path)
    ap.add_argument("--chain-id", default="A")
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    name = args.name

    # HIS 68 -> HIE, per BUILD-DECISIONS.md.  Written out so the file tleap
    # reads is on disk beside its input and can be inspected.
    src = args.pdb.read_text().splitlines(keepends=True)
    his_to_hie = 0
    fixed = []
    for line in src:
        if line.startswith("ATOM") and line[17:20] == "HIS":
            line = line[:17] + "HIE" + line[20:]
            his_to_hie += 1
        fixed.append(line)
    leap_input_pdb = out / ("%s_leap_input.pdb" % name)
    leap_input_pdb.write_text("".join(fixed))

    leap_in = out / ("build_%s.in" % name)
    leap_in.write_text(
        "source %s\n"
        "mol = loadpdb %s\n"
        "check mol\n"
        "charge mol\n"
        "saveamberparm mol %s.prmtop %s.rst7\n"
        "savepdb mol %s_leap.pdb\n"
        "quit\n" % (args.leaprc, leap_input_pdb.name, name, name, name))
    leap_log = out / ("build_%s.log" % name)
    run(["tleap", "-f", leap_in.name], out, leap_log)

    leap_pdb = out / ("%s_leap.pdb" % name)
    mosaics_pdb = out / ("%s_mosaics.pdb" % name)
    write_mosaics_pdb(leap_pdb, mosaics_pdb, args.chain_id)

    resolved, problems = check_against_rtf(mosaics_pdb, args.rtf)

    matched_rst7 = out / ("%s_matched.rst7" % name)
    natom = overlay_coordinates(out / ("%s.prmtop" % name), leap_pdb,
                                mosaics_pdb, matched_rst7)

    import parmed as pmd
    check = pmd.load_file(str(out / ("%s.prmtop" % name)), xyz=str(matched_rst7))
    lookup = pdb_coordinate_map(mosaics_pdb)
    worst = 0.0
    for index, atom in enumerate(check.atoms):
        want = lookup[(atom.residue.idx + 1, atom.residue.name, atom.name)]
        got = check.coordinates[index]
        worst = max(worst, max(abs(got[k] - want[k]) for k in range(3)))
    if worst > 1e-9:
        raise RuntimeError("coordinates disagree by %.3e Angstrom after the "
                           "overlay" % worst)

    residues, seen = [], set()
    for _, resseq, resname, _, _ in read_pdb_atoms(mosaics_pdb):
        if resseq not in seen:
            seen.add(resseq)
            residues.append(resname)

    counts = {}
    for r in residues:
        counts[r] = counts.get(r, 0) + 1

    meta = {
        "name": name,
        "source_pdb": str(args.pdb),
        "leaprc": args.leaprc,
        "rtf": str(args.rtf),
        "his_atoms_renamed_to_hie": his_to_hie,
        "atoms": natom,
        "residue_count": len(residues),
        "residue_sequence": residues,
        "residue_counts": counts,
        "rtf_templates_resolved": resolved,
        "rtf_template_counts": {k: resolved.count(k) for k in sorted(set(resolved))},
        "rtf_problems": problems,
        "total_charge_e": round(sum(a.charge for a in check.atoms), 6),
        "max_coordinate_disagreement_angstrom": worst,
        "files": {
            "topology": "%s.prmtop" % name,
            "reference_coordinates": "%s_matched.rst7" % name,
            "mosaics_structure": "%s_mosaics.pdb" % name,
        },
    }
    (out / "build_metadata.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
