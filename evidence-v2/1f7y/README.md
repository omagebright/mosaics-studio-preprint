# 1F7Y: S15 protein + 16S rRNA fragment, one deck carrying both families

The first MOSAICS system run against a parameter deck that holds protein
(`PROT-`) and RNA (`RNA-`) atom types together. Decks and topology files:
`engine/params/db_ff19sb_ol3`, `db_ff14sb_ol3`,
`top_openmm-ff19sb_ol3_protein_rna.rtf`, `top_openmm-ff14sb_ol3_protein_rna.rtf`
(engine branch `bright/protein-nucleic-deck-1f7y`, commit 8832542).

System: 3321 atoms, chain B = 57 nt RNA (G5 ... C3), chain A = 86 aa protein
(PRO ... ILE, charged termini), net charge -49, vacuum, no cutoff. The
structure is `work/preprint-I/systems/builds/1f7y/` (shared by both force
fields; the leap PDBs are identical).

## Layout

```
ff19sb_ol3/            OL3 sourced BEFORE ff19SB in tleap (the correct ff19SB topology)
  topology/            prmtop, rst7 written from the MOSAICS PDB, build input + log, chain PDBs
  openmm/              openmm_sp.out   (scripts/openmm_single_point.py, Reference, NoCutoff)
  sander/              sander_sp.in/out (imin=1, maxcyc=1, cut=9999)
  mosaics-current/     sp_*.input, renamed PDB, mosaics_sp.out
  as-built/            the same four dirs on the ORIGINAL topology (ff19SB sourced first)
ff14sb_ol3/            same four dirs, the as-built topology (source order does not matter here)
pdb_to_rst7.py         writes the rst7 from the PDB, matched by (residue, atom name)
```

All three engines read the same coordinates: the rst7 is written FROM the
MOSAICS PDB (three decimals) and round-trips at 0.0 A. The tleap restart
differs from the PDB by up to 5e-4 A, which is worth 0.66 kcal/mol of
Lennard-Jones on this system, so the builds/ reference numbers (sander
13673.064686, OpenMM 13673.064702) are on a different geometry and are not
the comparison targets here.

## Result, kcal/mol, MOSAICS minus OpenMM

| term | ff19SB/OL3 | ff14SB/OL3 |
|---|---:|---:|
| bond | +1.9e-11 | +1.9e-11 |
| bend | -1.55e-03 | -1.55e-03 |
| torsion + improper | +7.8e-05 | +3.8e-04 |
| cmap | -1.6e-11 | 0 |
| nonbonded (LJ + Coulomb) | -3.67e-03 | -3.67e-03 |
| total, own constants | 13685.347922 vs 13685.353068 (-5.1e-03) | 14238.325033 vs 14238.329881 (-4.8e-03) |
| total, CODATA | 13685.351529 vs 13685.353068 (-1.5e-03) | 14238.328640 vs 14238.329881 (-1.2e-03) |
| sander, CODATA | 13685.353082 | 14238.329782 |

The bend residual is theta0 written as 109.5 in the OL3 deck where the
prmtop carries 109.50004693 (1.91113635 rad); the pure-RNA control has it
too. The non-bonded residual is the Coulomb constant (MOSAICS
332.06344122017, OpenMM 332.063713); LJ agrees with sander to 1e-4 on both
the plain and the 1-4 sum.

## The sulfur trap, ff19SB only

`ff19sb_ol3/as-built/` is the topology from
`source leaprc.protein.ff19SB` then `source leaprc.RNA.OL3`. OL3 loads
`parm10.dat`, which overwrites the S/SH Lennard-Jones that `frcmod.ff19SB`
installs (GAFF2: R* 1.9825, eps 0.2824) with parm10's 2.0000/0.2500. Four
Met sulfurs; every other parameter identical (`prmtop_diff.py` in the
session). OpenMM on that topology: 13686.013579, 0.66 above the deck, all
of it Lennard-Jones (plain +0.603, 1-4 +0.058). The deck carries the ff19SB
values, as the OpenMM XML does, so the reference was rebuilt OL3-first.
