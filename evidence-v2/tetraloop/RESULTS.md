# tetraloop -- 2KOC, the UUCG hairpin

Run 2026-08-25.  Every number below was printed by a binary that was run here
and is in a file under `evidence-v2/tetraloop/`.  Nothing is estimated,
nothing is carried across from another system, and where an engine or a build
could not produce a number the cell says so and the raw output that says so is
committed beside it.

    STABLE   /Users/bright/Documents/MOSAICS/work/staging/bin/mosaics-3.9.1
             md5 6c268d2866adbcf189f03fa8afaa1b51
    CURRENT  /Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/examples/mosaics.x
             md5 1ba6ddfdd4e59e29a54bb15319c18f01

Both binaries print `>> version.3.9.1 <<`, so an output file cannot be
attributed to a build by reading it (gate report, Finding A).  Which build
wrote which file here is fixed by the directory it is in and by the md5 above.

One line, before the detail: **on this system 3.9.1 produced no energy at all
under either AMBER deck, and CURRENT produced no energy at all under
CHARMM36.  There is no force field on which both builds can be compared
number against number.**

---

## 1. The structure, and the chemistry as built

    source        evidence-v2/_structures/2KOC_A_model1.pdb
                  RCSB 2KOC, solution NMR, MODEL 1 of 20, chain A
    cleaned       294 heavy atoms, 14 residues, no HETATM, no altlocs
    built         447 atoms after the topology builder added 153 hydrogens
    sequence      G G C A C U U C G G U G C C
    net charge    -13

    residue    count   atoms each
    ---------------------------------------------------------------
    RG5          1     32     5-prime residue, free 5'-OH (HO5')
    RG           4     34     internal guanine  (2, 9, 10, 12)
    RC           4     31     internal cytosine (3, 5, 8, 13)
    RA           1     33     internal adenine  (4)
    RU           3     30     internal uracil   (6, 7, 11)
    RC3          1     32     3-prime residue, free 3'-OH (HO3')

    by base      A 1    C 5    G 5    U 3      -- all four RNA bases
    termini      5'-OH / 3'-OH.  The deposited bare 5'-phosphate on G1
                 (OP3, P, OP1, OP2) is removed, as BUILD-DECISIONS.md
                 requires and for the reason recorded there.

This is why the system is in the paper.  The committed evidence has RNA only
as one strand of the 1efs heteroduplex, whose RNA is `UUCUCUUCCUCUC` -- U and
C only, in a duplex.  2KOC carries A and G as well, and carries them in a
non-canonical loop rather than a helix.

**Hydrogens, and the one place this departs from BUILD-DECISIONS.md.**  Rule 3
of that file says hydrogens are tleap's for the AMBER decks and OpenMM's
template for the CHARMM path.  Here one hydrogen set -- tleap's, from
`leaprc.RNA.OL3` -- is used for all three force fields, and the same 447-atom
file is what MOSAICS reads and what the CHARMM36 OpenMM reference reads.  Two
reasons, both checkable:

  - the committed CHARMM36 controls are themselves AMBER leap geometry.
    `evidence/charmm36.txt` C-11 says so in as many words ("the oligo is
    AMBER leap geometry evaluated with CHARMM parameters"), and
    `tools/charmm_reference_energy.py` is written to take exactly such a file
    and rename `HO5'`/`HO3'` to `H5T`/`H3T` itself;
  - three force fields on one geometry is what makes a term-by-term
    force-field comparison mean anything.  Two geometries would put a
    coordinate difference into every CHARMM36-minus-OL3 difference.

The stronger committed rule -- the structure MOSAICS reads is the one the
reference topology is built from -- is obeyed exactly, and it needed work:
see section 2.

    every engine here reads     tetraloop_mosaics.pdb, 447 atoms
    MOSAICS naming              RG5/RG/RC/RA/RU/RC3, O1P/O2P, H2''/H2'
                                produced by tools/pdb_to_mosaics.py
                                --profile terminal, checked against the RTF

---

## 2. One correction to the protocol that had to be made first

tleap writes an rst7 to seven decimals and a PDB to three.  MOSAICS reads the
PDB.  Handing the rst7 to sander/OpenMM/GROMACS and the PDB to MOSAICS puts a
geometry difference of up to **4.963e-04 A** into the comparison.  On a first
pass that alone moved the bond term by 1.4e-03 kcal/mol and the torsion term
by 4.8e-03 kcal/mol -- three to four orders larger than the cross-engine
spread being measured, and it would have been read as a MOSAICS defect.

So the rst7 is rewritten from the PDB before any reference engine is run:
`topology/canonicalise.py`, committed, prints `max |rst7 - pdb| after = 0.0`.
All four engines then see one set of coordinates.  Both AMBER topology
directories carry the original tleap rst7 and the canonical one, so this is
auditable rather than asserted.

---

## 3. What the two builds did, per force field

    force field        3.9.1                          CURRENT
    ---------------------------------------------------------------------------
    OL3 (chiOL3)       REFUSES the deck, exit 1       runs, agrees with all
                       no topology, no energy         three references
    OL21/OL3           REFUSES the deck, exit 1       runs, agrees with all
                       no topology, no energy         three references
    CHARMM36 nucleic   runs to completion, exit 0     REFUSES at connectivity,
                       builds 14 separate molecules   exit 1, no energy
                       instead of 1 chain of 14;      -- the deck has no row
                       total out by +19640 kcal/mol   for 10 of the 13
                       CMAP: no line printed          backbone junctions

---

## 4. OL3 (chiOL3)

Reference built with `source leaprc.RNA.OL3` alone -- the force field a pure
RNA chain would actually be given.  `topology/build_ol3.in` is committed.
MOSAICS deck: `db_ol15_ol3_terminal` with
`top_openmm-ol15_ol3_terminal_hybrid_chidef.rtf`.  MOSAICS ships no OL3-only
deck; OL3 is the RNA half of every OL-family deck, and section 6 measures
whether that half really is shared rather than taking the claim on trust.

kcal/mol.  Electrostatics printed as each program prints it, then corrected to
332.0637133 kcal A/(mol e^2); each engine's own constant is in
`evidence/engines.txt`.  Nothing else is corrected.

    term          MOSAICS 3.9.1        MOSAICS CURRENT        sander            OpenMM           GROMACS
    bond          (no energy)             27.05879556      27.05880000      27.05879556      27.05879565
    angle         (no energy)            481.17851043     481.17860000     481.17861059     481.17861066
    dihedral      (no energy)            336.34612065     336.34610000     336.34612239     336.34611974
    lj            (no energy)            119.18151111     119.18150000     119.18151095     119.18147036
    coulomb       (no energy)           -616.31351163    -616.29270000    -616.31401225    -616.31400956
                 (electrostatic values above are as printed)
    coulomb*      (no energy)           -616.31401702    -616.31406882    -616.31401281    -616.31401012
    total         (no energy)            347.45092073     347.45093118     347.45102668     347.45098630

    MOSAICS CURRENT minus OpenMM, after correction:
      bond       -4.405e-13
      angle      -1.002e-04
      dihedral   -1.737e-06
      lj         +1.553e-07
      coulomb    -4.210e-06

    OpenMM minus GROMACS, the two references written independently of each
    other, on the same numbers:
      bond       -9.258e-08
      angle      -7.054e-08
      dihedral   +2.650e-06
      lj         +4.059e-05
      coulomb    -2.691e-06

    files  ol3/mosaics-current/mosaics_sp.out
           ol3/mosaics-3.9.1/mosaics_sp.out          (the refusal)
           ol3/sander/ol3_sander_sp.out
           ol3/openmm/openmm_sp.out
           ol3/gromacs/sp.xvg
           ol3/energies.tsv

The angle residual of -1.0e-04 is the same term and the same sign as the
committed 1efs tables, which carry -1.764e-04 on every force field.  It is not
new here.

### What 3.9.1 printed

    Reading in database from file .../db_ol15_ol3_terminal/mosaics_openmm-ol15_ol3.tors_and_impr
    ...
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Key word mol_type unknown in database line:
    ~torsion_parm
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1, while reading the third of the five deck files.  The bond deck (145
entries) and the bend deck (371 entries) were read without complaint; nothing
after the torsion deck is attempted -- no RTF, no topology, no energy, no
output file.

---

## 5. OL21/OL3

Reference built with `source leaprc.DNA.OL21` + `source leaprc.RNA.OL3`,
`topology/build_ol21_ol3.in`.  MOSAICS deck: `db_ol21_ol3_terminal` with
`top_openmm-ol21_ol3_terminal_hybrid_chidef.rtf`.

    term          MOSAICS 3.9.1        MOSAICS CURRENT        sander            OpenMM           GROMACS
    bond          (no energy)             27.05879556      27.05880000      27.05879556      27.05879565
    angle         (no energy)            481.17851043     481.17860000     481.17861059     481.17861066
    dihedral      (no energy)            336.34612065     336.34610000     336.34612239     336.34611974
    lj            (no energy)            119.18151111     119.18150000     119.18151095     119.18147036
    coulomb       (no energy)           -616.31351163    -616.29270000    -616.31401225    -616.31400956
                 (electrostatic values above are as printed)
    coulomb*      (no energy)           -616.31401702    -616.31406882    -616.31401281    -616.31401012
    total         (no energy)            347.45092073     347.45093118     347.45102668     347.45098630

    MOSAICS CURRENT minus OpenMM, after correction: identical to OL3 above,
    to every digit printed.

    files  ol21_ol3/mosaics-current/mosaics_sp.out
           ol21_ol3/mosaics-3.9.1/mosaics_sp.out     (the refusal)
           ol21_ol3/sander/ol21_ol3_sander_sp.out
           ol21_ol3/openmm/openmm_sp.out
           ol21_ol3/gromacs/sp.xvg
           ol21_ol3/energies.tsv

3.9.1 refuses this deck with the same error at the same file.

---

## 6. OL3 and OL21/OL3 are the same force field on this system, measured

`params/FORCEFIELD_USAGE.md` says "the OL-family decks differ primarily in the
DNA torsion set; the RNA half (OL3) is shared".  A pure RNA system is where
that can be tested, and it was, in two places:

  - **the reference.**  The two prmtops built by the two tleap inputs were
    compared parameter by parameter with ParmEd: identical atom types,
    identical charges (max difference 0.0), identical bond, angle and dihedral
    sets including phase, periodicity and force constant, identical
    Lennard-Jones A and B coefficient tables.  Both reference topologies are
    the same topology.  Consequently sander, OpenMM and GROMACS return the
    same numbers under both.  The two committed `sp.xvg` files differ only in
    their timestamp, their working-directory header and the GROMACS quote
    line; the single data line is character for character identical;
  - **the deck.**  MOSAICS CURRENT run on `db_ol15_ol3_terminal` and on
    `db_ol21_ol3_terminal` returns 347.4514261182977 both times, agreeing to
    all 16 printed digits on every term.

So the two force-field rows above are one measurement reported twice.  That is
the result, not a shortcut: on a pure RNA chain the OL21 half is never
reached, and the shared-RNA-half claim is confirmed rather than assumed.  The
two deck files are not identical -- `db_ol21_ol3_terminal` has 5 more bond
rows, 17 more bend rows and 14 more torsion rows -- and none of the extra rows
is touched by an RNA-only structure.

---

## 7. CHARMM36 nucleic

### 7a. The reference exists

`tools/charmm_reference_energy.py` on `charmm36_2024.xml`, 5TER/3TER, vacuum,
NoCutoff, Reference platform.  It built the system cleanly: 467 intra-residue
bonds from the CHARMM templates, 13 backbone O3'-P links, 0 deoxyribose
patches, no unbonded atoms, 447 particles.

    term                     kcal/mol      OpenMM force
    ---------------------------------------------------------------
    bond                 82.6094045415     HarmonicBondForce
    Urey-Bradley (1-3)  105.1030085698     HarmonicBondForce (split out)
    angle               551.8866392877     HarmonicAngleForce
    torsion             440.1552839637     PeriodicTorsionForce
    improper (harmonic)   0.1132546012     CustomTorsionForce
    1-4 Lennard-Jones   692.0249784181     CustomBondForce
    Lennard-Jones      -120.8488547469     CustomNonbondedForce
    electrostatics      355.5807512504     NonbondedForce
    correction map        0.0000000000     CMAPTorsionForce
    ---------------------------------------------------------------
    TOTAL              2106.6244658855

CMAP is 0.0000000000, which measures C-04 on a fourteen-residue mixed-sequence
RNA rather than asserting it: no nucleic residue in CHARMM36 uses a correction
map.

    file  charmm36_na/openmm/charmm36_openmm_reference.out

sander and GROMACS have no CHARMM36 path here.  The reason is written out in
`charmm36_na/sander/NOT-RUN.txt` and `charmm36_na/gromacs/NOT-RUN.txt`, with
the installed file list that supports it: `$AMBERHOME/dat/chamber/` ships
CHARMM36 protein and carbohydrate but **no `par_all36_na.prm` and no
`top_all36_na.rtf`**, so chamber cannot build a nucleic prmtop, and building
one out of the OpenMM System would make the "independent" engine a reader of
OpenMM's output.  The committed evidence takes the same position:
`evidence/coverage.txt` records `--` for sander and GROMACS on both CHARMM36
rows.

### 7b. CURRENT cannot run it, and the deck is why

    Setting up bonds
    ...
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Bond is not found in database.
    indices 32 33
    labels GUA-O3' GUA-P
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1.  The RTF loaded (1450 lines), all 24 residue templates resolved, the
residue topology and the velocities were set up; it stopped at the first
inter-residue bond, RG5 1 -> RG 2.

This is a deck coverage gap, not an engine defect, and it is exactly
enumerable.  The CHARMM36 deck types every atom as `<CHARMM residue>-<atom>`,
so a backbone link carries the identity of **both** residues it joins.  The
deck was generated from the ACGU/ACGT 4-mers in `work/charmm36-na/oligos`,
which contain only the junctions of that cycle.  The tetraloop needs thirteen
junctions and the deck has rows for three of them:

    junction   in the deck?     where in the tetraloop
    ------------------------------------------------------------
    A -> C     yes              A4  -> C5
    C -> G     yes              C8  -> G9
    G -> U     yes              G10 -> U11
    U -> A     yes              (not present in this sequence)
    G -> G     NO               G1 -> G2,  G9 -> G10
    G -> C     NO               G2 -> C3,  G12 -> C13
    C -> A     NO               C3 -> A4
    C -> U     NO               C5 -> U6
    U -> U     NO               U6 -> U7
    U -> C     NO               U7 -> C8
    U -> G     NO               U11 -> G12
    C -> C     NO               C13 -> C14

**10 of the 13 backbone bonds have no row in the deck**, and with them 40
bends across 32 distinct type triples.  The full list, with the residue pairs
that need each one, is in `charmm36_na/topology/deck_missing_rows.txt`.

Nothing about this is specific to PR #48 or to this build.  A homopolymer step
-- G-G -- is unrepresentable in the committed CHARMM36 nucleic deck, and so is
any dinucleotide step outside the four-cycle A->C->G->U->A.  Closing it means
regenerating the deck from oligos that carry every ordered base pair, which is
16 RNA steps and 16 DNA steps, not a change to the converter or the engine.
It is the same shape of problem as C-09/C-10, where a 4-mer could not define
an internal ADE.

### 7c. 3.9.1 runs it, and the reason it runs is a defect

3.9.1 completes with exit 0 on the same files.  It never asks for the
`GUA-O3'`-`GUA-P` row because it never makes that bond:

    Total number of molecules: 14
    Number of biopolymer molecules: 14

Fourteen molecules, one per residue, where CURRENT builds one chain of
fourteen.  This is the topology defect the gate report records as 3c, and on
this system it has the perverse effect of making the older build the only one
that produces a number.

kcal/mol.  MOSAICS electrostatics corrected to 332.0637133; OpenMM's constant
is already that value.

    term                  MOSAICS 3.9.1           OpenMM      3.9.1 - OpenMM
    bond                 140.3454002346    82.6094045415            +57.7360
    Urey-Bradley          57.0154605292   105.1030085698            -48.0875
    angle                463.9373994625   551.8866392877            -87.9492
    torsion              298.4024641562   440.1552839637           -141.7528
    improper               0.1050825835     0.1132546012             -0.0082
    1-4 Lennard-Jones    351.2086739952   692.0249784181           -340.8163
    Lennard-Jones      19986.9902810903  -120.8488547469         +20107.8391
    electrostatics       448.4079533106   355.5807512504
                     (as printed)
    electrostatics*      448.4083210159   355.5807512504            +92.8276
    correction map     no line printed      0.0000000000
    TOTAL              21746.4130830674  2106.6244658855         +19639.7886

Every term is wrong, and the Lennard-Jones term says why: with the residues
detached, every 1-2, 1-3 and 1-4 pair that crosses a backbone bond is counted
as a full non-bonded contact, and the phosphate-to-sugar contacts across a
severed link are at bonded distance.  19987 kcal/mol against -121.

Two things 3.9.1 does get right on this deck, and they are worth stating so
the "before" column is not a caricature: it parses the whole CHARMM36 deck
without complaint, and it evaluates both terms `evidence/charmm36.txt` said
needed no engine change -- Urey-Bradley comes back non-zero in the bend
record's own bond slot (`Bond-Bend energy 57.0154605292`) and the harmonic
impropers come back (`Tors-Improper energy 0.1050825835`, against OpenMM's
0.1132546012 on a molecule 3.9.1 has taken apart).  C-02 and C-03 hold against
the old binary here as they did on the 4-mer controls.

CMAP: 3.9.1's report has **no `-Cmap energy` line at all** -- not a zero, no
line.  The term is absent rather than evaluated.  On a nucleic system that
costs nothing, because the true value is 0.0000000000; on a protein it would
be the 1.92 kcal/mol the gate report measured on the ff19SB peptide.

    files  charmm36_na/mosaics-3.9.1/mosaics_sp.out
           charmm36_na/mosaics-current/mosaics_sp.out   (the refusal)
           charmm36_na/openmm/charmm36_openmm_reference.out
           charmm36_na/energies.tsv
           charmm36_na/topology/deck_missing_rows.txt

---

## 8. The honest control: 3.9.1 on its own era's files

3.9.1's refusal in sections 4 and 5 could be blamed on being handed decks from
a later branch, so it was handed the decks that ship beside it.  The staging
tree it lives in carries a full panel at
`work/staging/forcefields/panel/`, synchronised 2026-07-22 from
`pminary/MOSAICS` commit 5ad1e934, including `db_ol15_ol3_terminal`,
`db_ol21_ol3_terminal` and their RTFs, whose residues are named in 3.9.1's own
vocabulary (`ADE/GUA/CYT/URA`, `GR5`, `CR3`) rather than the current
`RG5/RG/RC/RA/RU/RC3`.

The structure was renamed into that vocabulary and both were run.  Both fail
identically:

    Reading in database from file .../staging/forcefields/panel/db_ol21_ol3_terminal/mosaics_openmm-ol21_ol3.tors_and_impr
    ...
    Key word mol_type unknown in database line:
    ~torsion_parm

exit 1.  The reason is that **all five deck files of both terminal profiles
are byte-identical between the staging tree and the current worktree**
(`cmp` clean on bond, bend, tors_and_impr, onfo, vdw for both OL15/OL3 and
OL21/OL3).  Only the RTFs differ, in 52 lines, of which 24 are `RESI` names.
So 3.9.1 cannot load the terminal deck that ships in its own staging tree.

There is no other route.  Surveying every RTF in that tree for a residue
template that actually declares a ribose 2'-hydroxyl and a free 5'-hydroxyl:

    RTF                                              ATOM O2'   ATOM HO5'
    ------------------------------------------------------------------------
    amber94/top_all94_prot_na_chidef.rtf                  1          0
    ff14sb_openmm/top_openmm-ff14sb_protein.rtf           0          0
    kb_3pt/top_3pt_prot_na.rtf                            0          0
    panel/top_openmm-ol15_ol3_terminal_hybrid...rtf      12          8
    panel/top_openmm-ol21_ol3_terminal_hybrid...rtf      12          8
    panel/top_all99-bsc1_prot_hybrid_chidef.rtf          12          8
    panel/top_amber-ol24_ol3_terminal_hybrid...rtf       12          8
    panel/top_all99-bs0_prot_hybrid_chidef.rtf            4          0

3.9.1's own `amber94` deck -- the one nucleic deck in the tree with no
`mol_type` in it -- declares `ADE` with `H2'` and `H2''` on C2' and no O2' and
no HO5': it is DNA chemistry with no terminal variants.  Every RTF in that
tree that can express a 5'-OH ribose belongs to the panel, and every panel
torsion deck carries the vocabulary 3.9.1 cannot parse.

    files  ol3/mosaics-3.9.1/native_naming_control/ctl_sp.out
           ol21_ol3/mosaics-3.9.1/native_naming_control/ctl_sp.out

---

## 9. Composition and deck coverage

The referee's objection is that the committed systems exercise an unstated
fraction of each deck.  Here is the fraction, measured.

### Topology

    atoms 447    residues 14    net charge -13
    bonds 480    of which 13 are inter-residue O3'-P backbone links
    bends 857
    torsion paths (bonded 1-4 routes) 1280
    1-4 pairs 1138

### The AMBER reference topology, as ParmEd reads it

    distinct AMBER atom types used            27
      C C4 C5 CA CB CI CP CQ CS CT H H1 H2 H4 H5 HA HO
      N* N2 NA NB NC O O2 OH OS P
    bond parameter rows used                  29 of 46 in the prmtop
    angle parameter rows used                 58 of 94
    dihedral terms                          1783   (1688 proper, 95 improper)
    distinct dihedral parameter rows          62 of 62   (58 proper, 4 improper)
    Lennard-Jones pair rows                  105   (14 LJ types, upper triangle)

### The MOSAICS decks, rows touched

Computed with `ol3/topology/coverage.py`: each atom is given its MOSAICS type
from the RTF, the structure's bonds/bends/1-4/non-bonded type tuples are
enumerated, and each deck row is asked whether any tuple matches it.

    OL15/OL3 terminal deck (db_ol15_ol3_terminal)
      term               deck rows   touched    of deck   distinct tuples
      bond                     145        46      31.7%        46
      bend                     371        94      25.3%        94
      onefour (1-4 LJ)        1596        97       6.1%        97
      inter (LJ)              1596       378      23.7%       378
      torsion + improper      1404       544      38.7%       164 quadruplets
        mol_type{rna}          548       485
        mol_type{dna}          553         0
        no mol_type            303        59

    OL21/OL3 terminal deck (db_ol21_ol3_terminal)
      bond                     150        46      30.7%        46
      bend                     388        94      24.2%        94
      onefour (1-4 LJ)        1596        97       6.1%        97
      inter (LJ)              1596       378      23.7%       378
      torsion + improper      1418       544      38.4%       164 quadruplets
        mol_type{rna}          548       485
        mol_type{dna}          553         0
        no mol_type            317        59

    CHARMM36 nucleic deck (db_charmm36_na)
      bond                     237       147      62.0%       155  (8 unmatched)
      bend                     562       266      47.3%       298 (32 unmatched)
      onefour (1-4 LJ)       15753       403       2.6%       403
      inter (LJ)             15753      9453      60.0%      9453
      torsion + improper      1036       406      39.2%       462 quadruplets

    Two limits on the torsion figure, stated rather than hidden.  The match is
    on the four type slots with X as a wildcard; \mol_type and \res_context are
    NOT used, because the engine's own selection rule for them is not
    reimplemented here.  The torsion percentages are therefore upper bounds.
    And the enumeration walks bonded 1-4 routes, so the improper rows -- whose
    quadruplets are not bonded paths -- are not counted as touched; impropers
    are read from the RTF's own IMPROPER lines (63 in the CHARMM36 RTF, 95
    improper terms in the AMBER prmtop).

    files  ol3/topology/deck_coverage.txt
           ol21_ol3/topology/deck_coverage.txt
           charmm36_na/topology/deck_coverage.txt
           charmm36_na/topology/deck_missing_rows.txt

Two things fall straight out.  A pure RNA structure touches **none** of the
553 `mol_type{dna}` torsion rows -- which is the point of the tag and is the
first direct confirmation that the tag partitions the deck as intended.  And
the CHARMM36 numbers are the only ones where the system needs more tuples than
the deck can supply: 155 bond tuples against 147 rows matched, 298 bend
triples against 266.  That difference is section 7b.

---

## 10. What failed, and why

    1  3.9.1 cannot load the OL3 deck.
       "Key word mol_type unknown in database line: ~torsion_parm", exit 1,
       while reading db_ol15_ol3_terminal/*.tors_and_impr.  1101 of the 1404
       torsion rows carry \mol_type and \res_context; neither string occurs
       anywhere in the 3.9.1 binary.  No topology, no energy, no output.

    2  3.9.1 cannot load the OL21/OL3 deck.  Same error, same file position,
       1101 of 1418 rows.

    3  3.9.1 cannot load the terminal decks that ship in its own staging
       tree either -- they are byte-identical to the current ones.  Its only
       mol_type-free nucleic deck, amber94, has no ribose O2' and no 5'-OH
       terminus in any residue template, so it cannot represent this molecule
       at all.

    4  CURRENT cannot run CHARMM36 on this system.  "Bond is not found in
       database. labels GUA-O3' GUA-P", exit 1.  The committed CHARMM36
       nucleic deck was generated from ACGU/ACGT 4-mers and carries only the
       four junctions of that cycle; 10 of this hairpin's 13 backbone
       junctions, and 40 bends across 32 type triples, have no row.  A
       homopolymer step is unrepresentable.  Deck, not engine: the same six
       files, unchanged, run the committed 4-mer controls.

    5  3.9.1 under CHARMM36 returns a wrong energy rather than refusing.  It
       builds 14 molecules instead of 1 chain, and the total is
       21746.413 against the OpenMM reference 2106.624, +19640 kcal/mol.
       Its Lennard-Jones term alone is +20108 out.

    6  3.9.1 prints no -Cmap energy line under CHARMM36.  The term is absent,
       not zero.  Costless here, because the true CMAP energy of an RNA
       strand is 0.0000000000 (measured, section 7a), but it is the same
       absence the gate report priced at 1.92 kcal/mol on the ff19SB peptide.

    7  sander and gmx_d have no CHARMM36 nucleic path in this environment.
       No par_all36_na.prm / top_all36_na.rtf in $AMBERHOME/dat/chamber, no
       .psf; the only source is OpenMM's own XML, so a "reference" from them
       would be reading OpenMM's output.  Recorded in NOT-RUN.txt in both
       directories.  Matches evidence/coverage.txt, which carries -- for both
       on CHARMM36.

    Not a failure, but it must be said: 3.9.1 and CURRENT are compared
    number against number on nothing.  Under both AMBER decks 3.9.1 has no
    number; under CHARMM36 CURRENT has none.  The before/after on this system
    is a capability comparison, not a residual comparison.

---

## 11. Files

    tetraloop/
      RESULTS.md                                 this file
      <ff>/energies.tsv                          machine-readable table
      <ff>/topology/                             tleap input, log, prmtop, both
                                                 rst7 files, tleap PDB, the
                                                 MOSAICS PDB and its mapping
                                                 report, canonicalise.py,
                                                 deck_coverage.txt
      <ff>/sander/                               sander_sp.in and mdout
      <ff>/openmm/                               openmm_sp.out (LJ/Coulomb split),
                                                 the committed tool's tsvs,
                                                 the scripts that produced them
      <ff>/gromacs/                              sp.mdp, sp.log, sp.xvg, grompp.log
      <ff>/mosaics-3.9.1/                        STABLE input, stdout, stderr,
                                                 exit code, and for the OL decks
                                                 native_naming_control/
      <ff>/mosaics-current/                      CURRENT input, stdout, stderr,
                                                 exit code, energy files
      charmm36_na/{sander,gromacs}/NOT-RUN.txt   why those two cells are empty

    GROMACS settings follow engines.txt: box (2*(d+1)+2)*10 A and
    rlist = rvdw = rcoulomb = d+1 nm, d = 3.3973354962 nm the largest
    interatomic distance, so cutoff 4.397335 nm and box 107.947 A.
    Energies are read from sp.xvg, never from the gmx_d screen output.
