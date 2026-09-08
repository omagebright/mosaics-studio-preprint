# F-23: the Tax nonamer number was computed on a structure not on `main`

Found 5 September 2026 while reconciling the deck with the four merges of 4 September.

## What was wrong

The tax9 rows in `work/preprint-I/energies/MATRIX.json` (208.442972 ff14SB,
149.508639 ff19SB) and the sentence in PR #65's body ("moves the Tax 9-mer from
240.983628 to 208.442972 under ff14SB") were computed on
`work/preprint-I/systems/normalized/tax_peptide_capped_rebuilt.pdb`. That file is a
rebuild from raw 1HHK chain C with hydrogens re-added, not the file #65 merged:

| comparison against `origin/main` `ea00b02` | atoms moved | X-H mean (A) | bond term, OpenMM |
|---|---|---|---|
| `main` before #65 (`30efc88`) | 12, all cap atoms | 1.076 | 5.27 |
| `tax_peptide_capped_rebuilt.pdb` | 77 | 1.009 | 201.00 |

So the deck quoted a value no committed structure reproduces. The 240.98 in the PR
text is the pre-#65 file and is right; the 208.44 is not the post-#65 file.

## Corrected numbers, on `examples/mhc_peptide_ff14sb/tax_peptide_capped.pdb` at `ea00b02`

`work/preprint-I/energies/rebuild_tax9_from_main.py`, output in `paper/evidence-v2/tax9/`.

| force field | MOSAICS | sander (CODATA) | OpenMM | MOSAICS-OpenMM | sander-OpenMM |
|---|---|---|---|---|---|
| ff14SB | -4.074209 | -4.074064 | -4.074077 | 1.3e-04 | 1.4e-05 |
| ff19SB | -63.456307 | -62.994264 | -62.994267 | **4.6e-01** | 3.5e-06 |

## The 0.46 on ff19SB is one improper at one nitrogen

Every term agrees to 1e-5 or better except Tors-Improper: MOSAICS 4.2438, AMBER 4.7057,
difference 0.4619. Recomputing every improper from the prmtop in both atom orders
(`compare_protein_impropers.py` plus the order check below) locates all of it at LEU2 N:

```
prmtop order   C(ACE1) CA(LEU2) N(LEU2) H(LEU2)   phi = +115.45   E = 1.7938
RTF order      -C       H        N       CA        phi = -128.92   E = 1.3319
```

An AMBER improper is `k(1+cos(2 phi - 180))` over the four atoms in a fixed order with
the third atom central. tleap writes `C CA N H`; the MOSAICS RTF record is
`IMPROPER -C H N CA`. Those are the same planarity restraint only when the centre is
near planar. At the other 62 backbone nitrogens in BPTI and tax9 the two orders differ
by under 0.05 each and by 0.03 in total. At LEU2 they differ by 0.46 because LEU2's amide
hydrogen is 60 degrees out of plane:

```
LEU2   C-N-H 99.1   H-N-CA 110.6   O=C-N-H dihedral -121.2     (planar amide: ~120, ~119, 0)
LEU3   C-N-H 118.5  H-N-CA 115.2   O=C-N-H dihedral -174.0
PHE4   C-N-H 118.0  H-N-CA 117.0   O=C-N-H dihedral -175.8
```

**That hydrogen is a second geometry defect in the committed structure, and #65 did
not touch it.** #65 moved exactly the twelve cap atoms. LEU2:H was placed against the
old 90-degree ACE cap and stayed where it was when the cap moved under it. The
1.79 kcal/mol in that one improper (typical site: 0.03) is the same kind of finding as
the cap defect: invisible to a cross-engine comparison, visible only by measuring the
structure.

Two things follow, kept separate:

1. **Structure:** LEU2:H needs re-placing in the peptide plane, as a follow-up to #65.
   One atom. `tools/rebuild_peptide_caps.py --check` should also measure the first
   amide H, or the same defect comes back with the next capped peptide.
2. **Engine/deck:** MOSAICS's per-residue `IMPROPER` atom order differs from tleap's at
   the backbone N. On a well-built structure the effect is ~1e-3 per residue; it is not
   the cause of any published number, but it is why an out-of-plane hydrogen shows as
   an "engine disagreement" instead of the geometry warning it is. Worth one line in
   the RTF README, not a code change.

## The tax9 atom names

`normalize_atom_names.py` renamed 7 atoms (ACE HH31-33 to H1-3, NME CH3 to C and
HH31-33 to H1-3) so tleap's ff14SB/ff19SB templates match. Geometry untouched.

## BPTI ff19SB 0.042, closed the same evening

`5PTI/ff19SB` MOSAICS-OpenMM is -4.2080e-02, all in the dihedral term. It is the same
mechanism as LEU2 above, at one residue: **ASP3 N**. The RTF record `-C H N CA` and
tleap's `C CA N H` read phi = +158.4 and -156.8 at that nitrogen, E = 0.2983 against
0.3403, difference -0.042041. Phase rounding (deck exactly pi, prmtop 3.141594 rad)
accounts for the rest: +3.1e-06 improper, -4.2e-05 proper. Residual -8.5e-08.

An earlier pass in this file said the OXT improper was missing and left 0.0119
unexplained. That was a bookkeeping error: residue 58 is `CALA`, not `ALA`, and
`CALA` carries `IMPROPER CA O C OXT` (RTF line 966) in tleap's order. Nothing is missing.

Closure scripts and their output: `paper/evidence-v2/bpti/ff19sb/closure/` (REPORT.md
has the exact commands). Bound: emulation of the RTF+deck reproduces MOSAICS's printed
improper to 6.5e-14 and proper to 1.7e-10; the engine's own improper list was not probed.

## The one standing item this leaves

The MOSAICS RTF writes the backbone-N improper as `-C H N CA` for every standard residue
and as `-C CA N H` for the C-terminal variants. tleap writes `C CA N H` everywhere. At a
planar nitrogen the two orders agree to under 0.05 kcal/mol per residue; at an out-of-plane
one they do not, and the gap reads as an engine disagreement instead of the geometry
warning it is. Whether the 78 standard records should be re-ordered to match AMBER is a
question for Peter: it is a deck change, it moves every ff14SB/ff19SB protein number by
about 1e-3, and it changes no physics. Not done.
