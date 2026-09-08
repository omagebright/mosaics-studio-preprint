# ddd -- the Drew-Dickerson dodecamer, PDB 1BNA

Run 2026-08-25. Every number below was printed by one of the five programs and
is in a file under `evidence-v2/ddd/`. Nothing is estimated, and no number is
carried across from another system.

    STABLE   /Users/bright/Documents/MOSAICS/work/staging/bin/mosaics-3.9.1
             md5 6c268d2866adbcf189f03fa8afaa1b51   built 21 Jul 2026
    CURRENT  /Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/examples/mosaics.x
             md5 1ba6ddfdd4e59e29a54bb15319c18f01   built 17 Aug 2026

Both print `>> version.3.9.1 <<`, so an output file cannot be attributed to a
build by reading it (gate Finding A). Attribution here is by directory.

References: sander 24.0, OpenMM 8.5.2 Reference platform, GROMACS 2025.4
double precision (`gmx_d`). Vacuum, dielectric 1, no cutoff, no temperature,
no move, single point. Electrostatics corrected to 332.0637133 kcal A/(mol
e^2) using the per-engine constants in `evidence/engines.txt`.

---

## 1. Composition

    source            evidence-v2/_structures/1BNA_AB.pdb  (cleaned 1BNA, both chains)
    chains            A residues 1-12, B residues 13-24, CGCGAATTCGCG each
    heavy atoms       486 as deposited-and-cleaned
    hydrogens         272, all added by tleap; none from the depositor
    atoms built       758
    residues          24
    net charge        -22

    residue     count   note
    ---------------------------------------------------------------
    DC5           2     5-prime-OH terminus, chains A and B
    DG3           2     3-prime-OH terminus, chains A and B
    DC            6
    DG            6
    DA            4
    DT            4

    by base       DA 4   DC 8   DG 8   DT 4      all four DNA bases
    termini       2 x 5-prime-OH, 2 x 3-prime-OH  both true chain ends

This closes the two gaps `evidence/systems.txt` records for 1efs, which has
only DG and DA on its DNA strand and no true terminus anywhere: 1efs is padded
so that every compared residue is interior.

Bonded terms in the built topology (OL21/OL3; identical for OL15/OL3 and
parmbsc0, parmbsc1 differs only in torsion count):

    bonds 816   angles 1476   proper torsion terms 3024   improper terms 260

### Parameter rows this system touches

Counted from the AMBER topology, one row per distinct
(atom-type tuple, force constant, period, phase). Set against the same count
for the committed 1efs OL21/OL3 topology, rebuilt for this comparison with
`tools/build_1efs_amber_reference.py` into `_work/1efs_ref/`.

    term              ddd   1efs   shared   ddd only   1efs only   union
    ---------------------------------------------------------------------
    bond               50     52       38         12          14      64
    angle             103    110       76         27          34     137
    proper torsion    224    253      155         69          98     322
    improper           38     28       20         18           8      46
    atom types         29     30       27          2           3      32

    atom types in ddd and not in 1efs:   C1, CM
    atom types in 1efs and not in ddd:   C4, CI, CS

So ddd is not a smaller copy of 1efs. It exercises 12 bond rows, 27 angle
rows, 69 proper-torsion rows and 18 improper rows the committed system never
reaches, and adding it takes the union coverage of the OL21/OL3 deck from 253
to 322 proper-torsion rows.

Deck sizes for scale (rows in the MOSAICS deck, DNA and RNA together, so a
pure-DNA system can only ever touch a part):

    db_ol21_ol3_terminal    bond 150   bend 388   tors_and_impr 1418   onfo/vdw 1596
    db_ol15_ol3_terminal    bond 145   bend 371   tors_and_impr 1404   onfo/vdw 1596
    db_bsc1                 bond 144   bend 350   tors_and_impr  660   onfo/vdw 2211
    db_bs0                  bond 175   bend 567   tors_and_impr 1060   onfo/vdw 2145

CHARMM36 deck generated from this structure (`charmm36/topology/deck-from-ddd/`),
as reported by the converter reading the built System:

    bonds 136   angles 280   Urey-Bradley 66   harmonic impropers 11
    torsion types 374 (64 of them written as energy 0)   1-4 pairs 1904
    LJ classes 32 over 117 types   NBFix pairs 0
    round-trip worst disagreement 8.9e-16 Angstrom on sigma

---

## 2. How the topologies were built

The committed rule holds throughout: **the structure MOSAICS reads is the one
the reference topology is built from.** For every force field the chain is

    evidence-v2/_structures/1BNA_AB.pdb
      -> TER inserted between chain A and chain B (coordinates untouched)
      -> tleap  ->  prmtop + rst7 + PDB          <- sander, OpenMM, GROMACS
      -> tools/pdb_to_mosaics.py --profile terminal   <- MOSAICS, both builds

tleap inputs are committed as `<ff>/topology/build_<ff>.in`, with the log
beside them.

Two protocol points that change numbers and are therefore stated:

**The rst7 is rewritten from the PDB.** tleap places hydrogens at full
precision but `savepdb` writes three decimals, so the reference engines and
MOSAICS would otherwise read coordinates differing by up to 5e-4 Angstrom on
the 272 hydrogens. That is worth 0.044 kcal/mol on the bond term alone --
larger than everything being measured. The rst7 was therefore reloaded from
the saved PDB before any reference engine was run
(`_work/`, and `<ff>/topology/ddd_<ff>_leap.rst7` keeps the untouched tleap
output). After this step every engine reads the same coordinates and the bond
term agrees to 5e-13.

**The MOSAICS proposal is `tors` with all sigmas zero, not `cart`.** With
`prop_type{tors}` and zero sigmas the step-1 report equals the initial energy
to 1e-13, verified in every run here. With `prop_type{cart}` it does not: the
step-1 report describes a configuration that has already moved. See Finding D.

---

## 3. Energies, force field by force field

kcal/mol. `sander` prints four decimals, so a sander number is known to about
2e-04. GROMACS values are read from `sp.xvg` and converted at 4.184 kJ/kcal.
GROMACS is the only reference that separates proper from improper torsion.

### OL21/OL3

    term                     MOSAICS 3.9.1    MOSAICS current       sander           OpenMM          GROMACS
    bond                        REFUSED         94.35242574      94.35240000      94.35242574      94.35242567
    angle                       REFUSED        829.11679555     829.11710000     829.11713378     829.11713384
    torsion (proper+improper)   REFUSED        651.66162827     651.66010000     651.66008445     651.66063456
      proper only               REFUSED        651.02934460               -                -       651.02936353
      improper only             REFUSED          0.63228367               -                -         0.63127103
    lennard-jones (total)       REFUSED       -103.23432771    -103.23430000    -103.23432793    -103.23435540
    electrostatics as printed   REFUSED        669.59161828     669.56900000     669.59211450     669.59217878
    electrostatics corrected    REFUSED        669.59216736     669.59221607     669.59211510     669.59217938

    MOSAICS current minus OpenMM, electrostatics after correction:
      bond       +5.542e-13
      angle      -3.382e-04
      torsion    +1.544e-03
      lj         +2.119e-07
      coulomb    +5.226e-05
    MOSAICS current minus GROMACS, torsion split:
      proper     -1.893e-05
      improper   +1.013e-03

    files  ol21_ol3/mosaics-current/mosaics_sp.out
           ol21_ol3/mosaics-3.9.1/mosaics_sp.out
           ol21_ol3/sander/ol21_ol3_sander_sp.out
           ol21_ol3/openmm/openmm_sp.out
           ol21_ol3/gromacs/sp.xvg

### OL15/OL3

    term                     MOSAICS 3.9.1    MOSAICS current       sander           OpenMM          GROMACS
    bond                        REFUSED         94.35242574      94.35240000      94.35242574      94.35242567
    angle                       REFUSED        829.11679555     829.11710000     829.11713378     829.11713384
    torsion (proper+improper)   REFUSED        648.70607901     648.70450000     648.70453519     648.70508533
      proper only               REFUSED        648.07379534               -                -       648.07381429
      improper only             REFUSED          0.63228367               -                -         0.63127103
    lennard-jones (total)       REFUSED       -103.23432771    -103.23430000    -103.23432793    -103.23435540
    electrostatics as printed   REFUSED        669.59161828     669.56900000     669.59211450     669.59217878
    electrostatics corrected    REFUSED        669.59216736     669.59221607     669.59211510     669.59217938

    MOSAICS current minus OpenMM, electrostatics after correction:
      bond       +5.542e-13
      angle      -3.382e-04
      torsion    +1.544e-03
      lj         +2.119e-07
      coulomb    +5.226e-05
    MOSAICS current minus GROMACS, torsion split:
      proper     -1.895e-05
      improper   +1.013e-03

    files  ol15_ol3/{mosaics-current,mosaics-3.9.1,sander,openmm,gromacs}/

### parmbsc1

    term                     MOSAICS 3.9.1    MOSAICS current       sander           OpenMM          GROMACS
    bond                        REFUSED         94.35242366      94.35240000      94.35242369      94.35242376
    angle                       REFUSED        907.21237243     907.21270000     907.21274174     907.21274187
    torsion (proper+improper)   REFUSED        548.39218871     548.00690000     548.00695422     548.00750167
      proper only               REFUSED        547.76004711               -                -       547.37623064
      improper only             REFUSED          0.63214160               -                -         0.63127103
    lennard-jones (total)       REFUSED        -99.66702859     -99.66700000     -99.66703064     -99.66705832
    electrostatics as printed   REFUSED        669.46836750     669.44560000     669.46886362     669.46892806
    electrostatics corrected    REFUSED        669.46891648     669.46881180     669.46886422     669.46892866

    MOSAICS current minus OpenMM, electrostatics after correction:
      bond       -2.334e-08
      angle      -3.693e-04
      torsion    +3.852e-01        <-- see Finding A
      lj         +2.050e-06
      coulomb    +5.226e-05

    files  bsc1/{mosaics-current,mosaics-3.9.1,sander,openmm,gromacs}/

### parmbsc0

    NEITHER BUILD RUNS. The three references do:

    term                        sander           OpenMM          GROMACS
    bond                     94.35240000      94.35242574      94.35242567
    angle                   829.11710000     829.11713378     829.11713384
    torsion                 468.52210000     468.52209578     468.52209106
      proper only                     -                -       467.88980736
      improper only                   -                -         0.63228370
    lennard-jones          -103.23430000    -103.23432793    -103.23435540
    electrostatics corr.    669.59221607     669.59211510     669.59217938

    files  bs0/{sander,openmm,gromacs}/, bs0/mosaics-current/mosaics_sp.out,
           bs0/mosaics-3.9.1/mosaics_sp.out

### CHARMM36 nucleic

Compared against OpenMM only, as the committed route requires: sander and
GROMACS have no CHARMM36 topology here and `evidence/coverage.txt` records the
same for rna_acgu and dna_acgt. Reference is
`tools/charmm_reference_energy.py` on `charmm36_2024.xml`, 5TER/3TER. The deck
is generated from this structure's own built System with
`tools/openmm_system_to_mosaics.py`, which is the committed C-09/C-11 route.

    term                   MOSAICS 3.9.1    MOSAICS current            OpenMM
    bond                    123.23959455      123.43677555      123.43677555
    Urey-Bradley            209.15116252      283.72075760      283.72075760
    angle                   650.95888039      665.00808332      665.00808332
    torsion (proper)        617.75241199      825.68899872      825.68899872
    improper (harmonic)       0.75127317        0.75127317        0.75127317
    1-4 Lennard-Jones       487.55507048      518.76542418      518.76542418
    Lennard-Jones        312306.17292295     -371.89918926     -371.89918926
    electrostatics printed  1138.42896112     1463.54898620     1463.55018537
    electrostatics corr.    1138.43029400     1463.55018635     1463.55018669
    correction map                    -                 -         0.00000000
    TOTAL (corrected)    315534.01161027     3509.02230962     3509.02230997

    MOSAICS current minus OpenMM:
      bond       +2.649e-11
      Urey-Bradley +1.410e-11
      angle      -3.729e-11
      proper     +2.149e-11
      improper   -6.198e-11
      1-4 LJ     -3.558e-11
      lj         +6.594e-12
      coulomb    -3.482e-07   (after correction)
      total      -3.483e-07   (after correction)

    MOSAICS 3.9.1 minus OpenMM: every term. LJ is out by +3.13e+05.

    files  charmm36/mosaics-current/deck-from-ddd/mosaics_sp.out
           charmm36/mosaics-3.9.1/deck-from-ddd/mosaics_sp.out
           charmm36/openmm/openmm_charmm36_sp.out
           charmm36/topology/deck-from-ddd/                     the deck
           charmm36/topology/deck-from-ddd-build.log            its round trip

CHARMM36 is the one force field on this system where MOSAICS reproduces the
reference to the digits the reference prints: eleven significant figures on
every bonded and Lennard-Jones term, and 3.5e-07 on electrostatics after the
constant correction. It is also the only one of the five whose deck was
generated from this structure rather than from another.

---

## 4. Before and after, in one table

    force field   3.9.1                                     current
    ---------------------------------------------------------------------------
    OL21/OL3      REFUSES the deck: "Key word mol_type      runs; agrees with
                  unknown". 1101 mol_type rows in the       OpenMM to 5e-13
                  terminal deck. No energy of any kind.     (bond) and 3.4e-04
                                                            (angle); torsion
                                                            +1.5e-03, carried
                                                            by the improper
    OL15/OL3      REFUSES the deck: same error, same        as OL21/OL3
                  1101 rows
    parmbsc1      REFUSES the deck: same error, 207         runs; bond, angle,
                  mol_type rows                             LJ and Coulomb all
                                                            agree; torsion is
                                                            out by +0.385
                                                            (Finding A)
    parmbsc0      REFUSES the deck: same error, 191         REFUSES the RTF:
                  mol_type rows                             "Total 28 atoms in
                                                            1. residue(DC)
                                                            doesn't fit
                                                            topology"
    CHARMM36      runs, and builds 24 detached molecules    committed deck:
                  instead of 2 chains of 12. LJ +312306     REFUSES, "Bend is
                  against -371.9. Total 315534.0 against    not found in
                  3509.0. Urey-Bradley and the harmonic     database ...
                  impropers are evaluated -- the improper   GUA-O3' DEOX-P
                  is exact to 11 figures -- so the term     CYT-O5'".
                  algebra is fine and the topology is not.  deck generated from
                                                            this structure:
                                                            runs and reproduces
                                                            OpenMM to 1e-11

The honest control, 3.9.1 given its own era's files:

    3.9.1 ships its own panel under staging/forcefields/panel. Its OL21/OL3 and
    OL15/OL3 terminal decks carry the same 1101 mol_type rows, so 3.9.1 refuses
    its own decks too -- checked, not assumed
    (ol21_ol3/mosaics-3.9.1/native_naming/). Its parmbsc1 deck has no mol_type
    rows and its RTF carries the 3.9.1-era terminal templates CD5 and GD3, so
    parmbsc1 is the one force field this system can be put to 3.9.1 fairly.
    Renamed to 3.9.1's vocabulary (DC5 -> CD5, DG -> GUD, and so on; atom names
    and coordinates untouched) it:

      - builds the molecule correctly: 2 chains, charge -22;
      - computes initial energies IDENTICAL to CURRENT on the same files, all
        16 digits: bond 94.35242366497357, bend 907.2123724296151,
        torsion 548.0386364645709, total 2119.40489087859;
      - then DIES with signal 10 (exit 138) inside the torsional move, before
        it reports any energy. Reproducible.

    Under a cartesian proposal it survives, but then reports a configuration
    that has already moved (Finding D), so that run is evidence of nothing but
    the crash being specific to the torsional path.

    So on this system 3.9.1's energy algebra for parmbsc1 is exactly CURRENT's.
    What it cannot do is load any deck that says which sugar a torsion belongs
    to, connect residues that are named in the current vocabulary, or complete
    a torsional step on a 758-atom duplex.

    files  bsc1/mosaics-3.9.1/native_naming/mosaics_sp_3.9.1.out   (crash)
           bsc1/mosaics-3.9.1/native_naming/mosaics_sp_CURRENT.out (same files)
           bsc1/mosaics-3.9.1/native_naming/cart/                  (cart run)
           ol21_ol3/mosaics-3.9.1/native_naming/mosaics_sp_3.9.1.out

---

## 5. What failed, and why

    F-ddd-1  parmbsc0, both builds. The parmbsc0 RTF
             (top_all99-bs0_prot_hybrid_chidef.rtf) has 8 residue templates,
             all internal: RA RU RC RG DA DT DC DG. It has no DC5 and no DG3,
             so a chain with true 5-prime and 3-prime ends cannot be expressed.
             params/README.md states this ("The bs0 RTF currently has no
             terminal-variant templates"). tools/pdb_to_mosaics.py --rtf
             refuses first; run anyway, CURRENT prints
             "Total 28 atoms in 1. residue(DC) doesn't fit topology" and 3.9.1
             never gets past the deck. The three references all produce
             energies, so the row is a MOSAICS-side gap and not a chemistry
             one.
             Logs: bs0/topology/pdb_to_mosaics_bs0_rtfcheck.log,
                   bs0/mosaics-{current,3.9.1}/mosaics_sp.out

    F-ddd-2  parmbsc1 topology, stock recipe. The committed parmbsc1 build
             recipe (oldff/leaprc.ff99 + parmBSC1.lib) has no PDB-v3 DNA
             terminal residue map: leaprc.ff99's addPdbResMap knows GUA/ADE/
             CYT/THY and the single letters G/A/C/T, but not DG/DA/DC/DT. Given
             a PDB-v3 file it builds INTERNAL residues at both chain ends and
             tleap silently adds a 5-prime phosphate -- 760 atoms, net charge
             -24 instead of -22. This is a silent wrong answer, not an error.
             Fixed here by adding the four lines that oldff/leaprc.ff99bsc0
             already ships; no atom type is touched. The failing build is kept.
             Files: bsc1/topology/stock_resmap_failure/

    F-ddd-3  CHARMM36 with the committed deck, CURRENT. db_charmm36_na was
             generated from 4-mers that carry every base as 5-prime, internal
             and 3-prime, but only in the sequences A-C-G-T. 1BNA is
             CGCGAATTCGCG and contains base steps those oligos never had, so
             the deck has no row for the backbone bend at a G-to-C step:
             "Bend is not found in database. labels GUA-O3' DEOX-P CYT-O5'".
             The engine is right to stop. Generating the deck from 1BNA's own
             System removes the gap and the energy then matches OpenMM to
             1e-11. This is the referee's objection made concrete: a deck
             validated on one set of oligos is not a deck validated for a
             sequence.

    F-ddd-4  CHARMM36, 3.9.1. Runs, but builds 24 separate molecules from a
             24-residue duplex, so every inter-residue bond, bend and 1-4 is
             missing and every pair that should be excluded is a full
             non-bonded contact. The same defect the gate recorded on
             dna_acgt and rna_acgu, at duplex size.

    F-ddd-5  OL21/OL3, OL15/OL3, parmbsc1, parmbsc0, 3.9.1. All four decks
             refused at the same line, with the same message. The vocabulary
             that says which sugar a torsion row belongs to did not exist in
             3.9.1, and its own shipped OL decks use it too.

---

## 6. Findings for the authors

### Finding A -- parmbsc1 has no terminal profile, and on a true-terminal chain that is worth 0.385 kcal/mol

`params/FORCEFIELD_USAGE.md` says it plainly: "a true-terminal chain run under
the standard profile (or vice versa) will give the wrong torsion energy". OL21,
OL15 and OL24 each ship a `_terminal` deck. parmbsc1 and parmbsc0 do not. So
1BNA -- a true-terminal chain -- can only be run under parmbsc1's standard
deck, and its torsion is out by +0.385 kcal/mol while every other term agrees
to 1e-06.

Measured, not inferred. Running CURRENT on this identical structure with the
OL21/OL3 **standard** deck instead of the terminal one gives torsion
652.02926728 against the terminal profile's 651.66162827 and OpenMM's
651.66008445 -- an error of +0.369, the same size and sign. Bond, angle,
Lennard-Jones and Coulomb are bit-identical between the two profiles.

    file  _work/diag_ol21_standard_profile/out.txt

The paper should either say parmbsc1 and parmbsc0 are not validated for
true-terminal systems, or a terminal profile should be built for them.

### Finding B -- the improper torsion, 1.0e-03 kcal/mol, invisible on 1efs

MOSAICS returns 0.63228367 for the OL21/OL3 improper energy. GROMACS returns
0.63127103 for the same force field on the same coordinates, and 0.63228370
for parmbsc0. MOSAICS agrees with the parmbsc0 value to eight figures.

The cause is atom ordering, and it was located rather than guessed. tleap under
DNA.OL21 writes the adenine N9 improper as C8/C4/N9/C1'; under ff99bsc0 it
writes C4/C8/N9/C1'. Swapping the first two atoms is not a symmetry of the term
when the group is not exactly planar. Six instances differ, on DA 5, 16 and 17,
and they account for the whole 1.0e-03.

This could not be seen on 1efs, whose entire improper energy is 5.3e-04
kcal/mol. On 1BNA the improper term is 0.63 kcal/mol and the disagreement is
0.16 per cent of it.

### Finding C -- gate Finding B is a reading error, and the committed CHARMM36 DNA row is correct

The gate report says the committed `mosaics_total_kcal_dna 1954.691252416373`
does not reproduce, and that CURRENT gives 849.7216661162275 on the same
structure and deck.

Both numbers come out of the same run. `1954.691252416373` is the **initial**
energy; `849.72` is the **step-1 report**, and the committed CHARMM36 example
input uses `\prop_type{cart}`, under which the step-1 report describes a
configuration that has already moved. Rerun with `\prop_type{tors}` and all
sigmas zero, CURRENT on `dna_acgt.pdb` reproduces the committed C-11 table on
every single term:

    term                 committed C-11      CURRENT, this run
    bond                 18.9218490859       18.92184908593143
    Urey-Bradley        167.6513797824      167.6513797823018
    angle               412.7880943139      412.7880943138928
    torsion             122.9466154569      122.9466154568938
    improper              0.0010255532        0.00102555323518276
    1-4 Lennard-Jones  1412.2714009823     1412.27140098232
    Lennard-Jones       176.6882040552      176.6882040550464
    electrostatics     -356.5773168135     -356.5773168134513
    TOTAL              1954.6912524164     1954.691252416168

Nothing in the CHARMM36 evidence needs correcting. What needs correcting is the
example input, or a note next to it: on a strict single point the reported
energy must be read from the initial block, or the proposal must be `tors`.

    file  charmm36/_control_dna_acgt/run/mosaics_sp.out

### Finding D -- `prop_type{cart}` with zero sigmas is not a null move

Consequence of Finding C, stated separately because it affects any future run.
With `\prop_type{tors}` and `\prop_tors_sig{0}`, MOSAICS' step-1 report equals
its initial energy to 1e-13 -- verified on every run in this directory. With
`\prop_type{cart}` and the same zero sigmas it does not: on 1BNA under
CHARMM36 the initial total is 3509.021109 and the step-1 report is 3504.591605.
Every committed CHARMM36 example input uses `cart`.

### Finding E -- 3.9.1 crashes on a torsional step on this system

Given its own parmbsc1 deck, its own RTF and its own residue naming, 3.9.1
builds 1BNA correctly and computes initial energies identical to CURRENT's to
all 16 digits, then terminates with signal 10 (exit 138) inside the torsional
move, before printing a step-1 report. Reproducible. Not a deck problem: the
same files under `\prop_type{cart}` complete.

### Finding F -- `tools/openmm_system_to_mosaics.py` writes into `params/`

Lines 849-851 copy the generated RTF over
`params/top_charmm36_nucleic.rtf` on every run, whatever `--outdir` says. Two
deck generations during this work therefore overwrote the committed 24-residue
CHARMM36 RTF in the CURRENT worktree with a 6-residue and then a 4-residue one.
**It has been restored**, verified byte-identical to
`params/db_charmm36_na/top_charmm36_nucleic.rtf` and to
`stack4-rewrite/params/top_charmm36_nucleic.rtf`
(md5 `e2eed6fb7dc2ee3c88ea01a40cca6ccb`). No other file in either worktree was
modified. The side effect should be removed or made opt-in before anyone else
runs that tool.

### Finding G -- MOSAICS truncates a database path at 128 characters

Both builds. A `\mol_parm_file{...}` or `\*_database_file{...}` whose value is
longer than 128 characters is cut, and the remainder is reported as an unknown
keyword:

    Key word -ddd/top_charmm36_nucleic.rtf unknown in input file line:
    ~sim_mol_def

The path that triggered it was 152 characters. It is silent about the real
cause. All runs here that use long paths were repeated with relative ones.

---

## 7. Coverage row, in the form of `evidence/coverage.txt`

      force field    system    MOSAICS 3.9.1   MOSAICS current  sander  OpenMM  GROMACS   result
      -------------------------------------------------------------------------------------------
      OL21/OL3       ddd            no              yes           yes     yes     yes      pass
      OL15/OL3       ddd            no              yes           yes     yes     yes      pass
      parmbsc1       ddd            no              yes           yes     yes     yes      part   F-ddd-A
      parmbsc0       ddd            no              no            yes     yes     yes      fail   F-ddd-1
      CHARMM36 na    ddd            runs, wrong     yes           --      yes     --       pass   F-ddd-3

    "part" for parmbsc1: bond, angle, Lennard-Jones and electrostatics agree
    with all three references to 1e-06 or better; the torsion does not, and the
    reason is that no terminal profile exists for this force field.

---

## 8. Files

    _structures/1BNA_AB.pdb            the cleaned deposit (not written here)
    ddd/_work/1BNA_AB_ter.pdb          the same file with a TER between chains
    ddd/_work/openmm_sp_amber_decomposed.py   OpenMM AMBER single point, committed format
    ddd/_work/collect.py               reads the four AMBER decks' output files and prints
                                       their tables; output in _work/tables.txt
    ddd/_work/add_chain_ids.py         puts chain IDs back on a tleap PDB
    ddd/_work/1efs_ref/                1efs OL21/OL3 topology, for the coverage comparison
    ddd/_work/diag_ol21_standard_profile/     Finding A

    ddd/<ff>/topology/build_<ff>.in            tleap input, committed
    ddd/<ff>/topology/build_<ff>.log           tleap log
    ddd/<ff>/topology/ddd_<ff>.prmtop/.rst7    reference topology (rst7 from the PDB)
    ddd/<ff>/topology/ddd_<ff>_leap.rst7       untouched tleap coordinates
    ddd/<ff>/topology/ddd_<ff>_mosaics.pdb     what MOSAICS reads
    ddd/<ff>/topology/amber_energy_summary.tsv sander terms, parsed
    ddd/<ff>/sander/<ff>_sander_sp.out         sander output
    ddd/<ff>/openmm/openmm_sp.out              OpenMM output, per force + LJ/Coulomb split
    ddd/<ff>/openmm/openmm_vs_amber_comparison.tsv   tools/openmm_amber_single_point.py
    ddd/<ff>/gromacs/{sp.mdp,sp.log,sp.xvg}    GROMACS settings and output
    ddd/<ff>/mosaics-3.9.1/mosaics_sp.out      STABLE
    ddd/<ff>/mosaics-current/mosaics_sp.out    CURRENT

    ddd/charmm36/topology/deck-from-ddd/       deck generated from this structure
    ddd/charmm36/openmm/openmm_charmm36_sp.out reference
    ddd/charmm36/_control_dna_acgt/            Finding C
    ddd/bsc1/topology/stock_resmap_failure/    F-ddd-2
    ddd/bsc1/mosaics-3.9.1/native_naming/      the honest control, and Finding E
