# Ubiquitin (PDB 1UBQ) -- ff14SB and ff19SB, MOSAICS 3.9.1 against MOSAICS current

Run 2026-08-25. Every number below is read from a file in this directory. No
number is estimated, carried across from another system, or reconstructed from
a similar run. Where a build produced no number, the cell says so.

    STABLE   /Users/bright/Documents/MOSAICS/work/staging/bin/mosaics-3.9.1
             md5 6c268d2866adbcf189f03fa8afaa1b51        built 21 Jul 2026
    CURRENT  /Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/examples/mosaics.x
             md5 1ba6ddfdd4e59e29a54bb15319c18f01        built 17 Aug 2026

Both binaries print `>> version.3.9.1 <<`, so an output file here cannot be
attributed to a build by reading its banner (Finding A of the gate report).
The directory name is the attribution.

Protocol as committed: vacuum, dielectric 1, no cutoff, no temperature, no
move, single point. Electrostatics compared only after correcting every engine
to 332.0637133 kcal A/(mol e^2), using the per-engine constants measured in
`evidence/engines.txt`. Nothing else is corrected.

---

## Headline

**On this system MOSAICS 3.9.1 returns no energy at all, for either force
field.** It is not that it is inaccurate here; it stops before the first term.
1UBQ is a whole protein with charged NH3+/COO- termini, and 3.9.1 has no way
to select a terminal residue template. The current build runs both force fields
and agrees with sander, OpenMM and GROMACS to about 1e-4 kcal/mol on a
1231-atom system.

**ff19SB CMAP coverage: the deck holds 16 correction maps and this system
selects 13 of them, across 74 backbone terms.** The committed peptide selects
4 maps across 4 terms. This is the direct answer to the coverage objection.

---

## 1. Composition

    source          evidence-v2/_structures/1UBQ_A.pdb  (cleaned from RCSB)
    chain           A, one chain, one molecule
    residues        76
    heavy atoms     601 as deposited; 1231 atoms after tleap rebuilt hydrogen
    termini         CHARGED, not capped -- NMET at residue 1 (+1),
                    CGLY at residue 76 (-1), per BUILD-DECISIONS.md rule 5
    histidine       HIS 68 -> HIE (10 atom records renamed before tleap)
    net charge      0, as built (tleap reports 0.000000; MOSAICS reports
                    5.55112e-16)
    disordered tail LEU 73 / ARG 74 / GLY 75 / GLY 76 at occupancy
                    0.45 / 0.45 / 0.25 / 0.25; all heavy atoms present,
                    coordinates taken as deposited

Residue types present -- 18 of the 20, missing only CYS and TRP:

    ALA 2   ARG 4   ASN 2   ASP 5   GLN 6   GLU 6   GLY 6   HIE 1   ILE 7
    LEU 9   LYS 7   MET 1   PHE 2   PRO 3   SER 3   THR 7   TYR 1   VAL 4

RTF templates resolved by the current build (19 distinct templates; the two
terminal templates are selected by the engine, not named in the PDB):

    ALA 2  ARG 4  ASN 2  ASP 5  CGLY 1  GLN 6  GLU 6  GLY 5  HIE 1  ILE 7
    LEU 9  LYS 7  NMET 1  PHE 2  PRO 3  SER 3  THR 7  TYR 1  VAL 4

Connectivity, identical under both force fields:

    bonds 1237    angles 2257    dihedral terms 5499

The MOSAICS structure file is byte-identical for the two force fields
(md5 `d430aac3b6c3e9bf0d81873f5085d6bb`), so any difference between the ff14SB
and ff19SB rows is the force field and nothing else. The two MOSAICS input
files are byte-identical between the STABLE and CURRENT directories; this was
checked with `diff` before each pair of runs.

### 1a. How much of each deck this system touches

The referee's objection is that the committed systems exercise an unstated
fraction of each deck. Stated, then, and next to the committed peptide so the
change is visible. Method and its limits are in `deck_coverage.py`; raw output
in `<ff>/topology/deck_coverage.txt` and
`<ff>/topology/deck_coverage_peptide_reference.txt`.

ff14SB deck (`params/db_ff14sb`):

    term          deck rows   ubiquitin   committed peptide
    atom types           34          28                  12
    bond                 90          56                  18
    bend                258         130                  38
    torsion             465     258 (231 wildcard-free)   59 (51)
    improper             61      20 (19 wildcard-free)     4 (4)
    vdw pair            595         406                  78
    onfo pair           595         406                  78

ff19SB deck (`params/db_ff19sb`):

    term          deck rows   ubiquitin   committed peptide
    atom types           36          29                  12
    bond                102          60                  18
    bend                329         141                  38
    torsion             616     262 (232 wildcard-free)   51 (43)
    improper             74      22 (21 wildcard-free)     4 (4)
    vdw pair            666         435                  78
    onfo pair           666         435                  78
    cmap maps            16          13                   4

The torsion row is given twice because a deck torsion row may carry `X`
wildcards and several rows can match one dihedral; the first figure counts deck
rows whose pattern is matched by at least one dihedral in the system (an upper
bound on rows used) and the figure in brackets counts wildcard-free rows
matched exactly (a lower bound). Bond, bend and pair rows are keyed on an exact
type tuple, so those counts are exact.

### 1b. ff19SB CMAP -- the number the coverage objection asks for

The current build prints, in `ff19sb/mosaics-current/mosaics_sp.out` line 147:

    cmap: 16 map(s) loaded (N=24), 74 backbone term(s)

The committed peptide prints the same line with `4 backbone term(s)`.

"16 map(s) loaded" is the deck's contents in both cases, so the selected count
has to be derived from the deck's `~cmap_binding` table, which keys a map id on
a residue name. For 1UBQ's 74 interior residues:

    map  0   x5    GLY
    map  1   x2    ALA
    map  2   x3    SER
    map  3   x7    THR
    map  4  x11    ILE, VAL
    map  5   x2    ASN
    map  6   x5    ASP
    map  7   x6    GLU
    map  8   x4    ARG
    map  9   x7    LYS
    map 10  x13    HIE, LEU, PHE, TYR
    map 11   x3    PRO
    map 12   x6    GLN

**13 of the 16 maps are selected.** The three not selected are map 13 (CYS,
free cysteine), map 14 (GLH) and map 15 (ASH) -- chemistries 1UBQ does not
contain. The committed peptide selects 4: maps 0, 1, 2 and 4.

This is confirmed independently by the reference topology: tleap's ff19SB
prmtop for this system carries exactly 13 `CMAP_PARAMETER_` sections
(`CMAP_PARAMETER_01` through `_13`), the same count arrived at from the deck.

---

## 2. Energies

Every value in kcal/mol. `angle` is MOSAICS's `-Bend energy`; the
Urey-Bradley slot (`Bond-Bend energy`) is 0 for both AMBER force fields, as it
must be. `dihedral` is proper plus improper. `lj` and `coulomb` are MOSAICS's
`-Total Lennard-Jones energy` and `-Total Coulomb energy`, that is direct plus
1-4, matching the committed convention (sander `VDWAALS + 1-4 VDW` and
`EEL + 1-4 EEL`).

### ff14SB on ubiquitin

    term          MOSAICS 3.9.1        MOSAICS current           sander           OpenMM          GROMACS
    bond              no output           125.61189200     125.61190000     125.61189200     125.61189197
    angle             no output           275.36422591     275.36490000     275.36487329     275.36487285
    dihedral          no output           917.58633065     917.58630000     917.58649896     917.58649594
    lj                no output           534.14577725     534.14580000     534.14577715     534.14575741
    coulomb           no output         -2399.66666449   -2399.58550000   -2399.66863505   -2399.66864316
               (electrostatic values above are as printed)
    coulomb*          no output         -2399.66863227   -2399.66870122   -2399.66863722   -2399.66864533   corrected to one constant
    total             no output          -546.96040647    -546.95980122    -546.95959583    -546.95962717

  MOSAICS current minus OpenMM, after correction:

    bond       +6.963e-13
    angle      -6.474e-04
    dihedral   -1.683e-04
    lj         +1.046e-07
    coulomb    +4.945e-06

  MOSAICS 3.9.1 minus OpenMM: not defined. 3.9.1 produced no energy.

  files  ff14sb/mosaics-3.9.1/mosaics_sp.out        (error, exit 1)
         ff14sb/mosaics-current/mosaics_sp.out
         ff14sb/sander/ubiquitin_sander_sp.out
         ff14sb/openmm/openmm_sp.out
         ff14sb/gromacs/sp.xvg

### ff19SB on ubiquitin

    term          MOSAICS 3.9.1        MOSAICS current           sander           OpenMM          GROMACS
    bond              no output           125.61189200     125.61190000     125.61189200     125.61189197
    angle             no output           275.36422591     275.36490000     275.36487329     275.36487285
    dihedral          no output           372.75171815     372.75040000     372.75052718     372.75052438
    cmap              no output            69.95284751      69.95280000      69.95284751                -
    lj                no output           533.93769334     533.93770000     533.93769324     533.93767352
    coulomb           no output         -2399.66666449   -2399.58550000   -2399.66863505   -2399.66864316
               (electrostatic values above are as printed)
    coulomb*          no output         -2399.66863227   -2399.66870122   -2399.66863722   -2399.66864533   corrected to one constant
    total             no output         -1022.05025538   -1022.05100122   -1022.05080401                -

  The two blank GROMACS cells are the F-06 conversion limit, at full scale.
  tleap's ff19SB prmtop for this system holds 13 distinct CMAP maps; the ParmEd
  route writes a single `[ cmaptypes ]` entry, keyed on the five backbone atom
  types `C N XC C N`, and applies that one map to all 74 backbone terms. This
  was checked in the written file, not assumed. GROMACS then reports
  61.54332385 kcal/mol against the 69.95284751 the other three engines agree
  on -- an 8.41 kcal/mol gap that is a property of the conversion route as
  used, not a disagreement about physics. The cells are left blank rather than
  filled with that number, exactly as the committed peptide table does. Every
  other GROMACS value on this system agrees with the rest to about 1e-5.

  MOSAICS current minus OpenMM, after correction:

    bond       +6.963e-13
    angle      -6.474e-04
    dihedral   +1.191e-03
    cmap       -2.700e-13
    lj         +9.904e-08
    coulomb    +4.945e-06

  MOSAICS 3.9.1 minus OpenMM: not defined. 3.9.1 produced no energy.

  files  ff19sb/mosaics-3.9.1/mosaics_sp.out         (error, exit 1)
         ff19sb/mosaics-3.9.1-nocmap/mosaics_sp.out  (error, exit 1)
         ff19sb/mosaics-current/mosaics_sp.out
         ff19sb/sander/ubiquitin_sander_sp.out
         ff19sb/openmm/openmm_sp.out
         ff19sb/gromacs/sp.xvg

### The two residuals worth naming

Neither is a failure; both are stated because a reader should not have to
discover them.

**angle, -6.47e-04 kcal/mol.** That is -2.9e-07 per angle over 2257 angles.
The committed peptide's angle residual is -9.57e-06 over 97 angles, or
-9.9e-08 per angle, and 1efs's is -1.76e-04. The bend deck writes its
equilibrium angle to at most two decimal places; the largest disagreement
between a deck `eq_bend` and the corresponding tleap equilibrium angle over the
130 angle types this system uses is 5.4e-05 degrees, at `O2-CO-O2`
(deck 126.0, tleap 126.00005400053564). The residual scales with the number of
angles, which is what a per-row rounding does.

**ff19SB dihedral, +1.19e-03 kcal/mol.** Splitting it against GROMACS, which
reports proper and improper separately, puts essentially all of it in the
improper term:

    ff19SB   proper    MOSAICS 364.386335264   GROMACS 364.386429493   -9.4e-05
             improper  MOSAICS   8.365382884   GROMACS   8.364094885   +1.3e-03

    ff14SB   proper    MOSAICS 909.231470963   GROMACS 909.231650096   -1.8e-04
             improper  MOSAICS   8.354859686   GROMACS   8.354845841   +1.4e-05

So the ff19SB harmonic-improper set carries a residual about a hundred times
the ff14SB one, on a term worth 8.4 kcal/mol out of a 1022 kcal/mol total. It
is recorded here as measured; no cause is asserted, because none was measured.

### Cross-engine spread, for scale

The three reference engines agree with each other independently of MOSAICS:

    GROMACS minus OpenMM        sander minus OpenMM
    bond      -2.602e-08        bond      +8.005e-06
    angle     -4.412e-07        angle     +2.671e-05
    dihedral  -3.021e-06        dihedral  -1.990e-04
    lj        -1.974e-05        lj        +2.285e-05
    coulomb   -8.113e-06        coulomb   -6.400e-05

sander prints four decimals, so its column is known to about 2e-04; that is the
floor on any sander difference, not a disagreement.

---

## 3. What MOSAICS 3.9.1 could and could not do

### 3a. ff14SB -- refused at the topology, before any energy

3.9.1 reads the input, reads all five ff14SB deck files (90 bond, 258 bend,
526 torsion, 595 onfo, 595 vdw entries -- the same counts the current build
reports), reads the 1231-line PDB and reads the RTF. It then stops:

    ===================================================================================
    Setting up residue topology
    -----------------------------------------------------------------------------------
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Total 19 atoms in 1. residue(MET) doesn't fit topology
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1, stderr empty, no energy of any kind. Residue 1 is a charged N-terminal
methionine: 19 atoms, the 17 of MET plus the two extra ammonium hydrogens. The
RTF holds the right template, `NMET`, with those 19 atoms. 3.9.1 does not
select it.

The current build does. `params/FORCEFIELD_USAGE.md` section 4 states the rule
it uses -- a first residue carrying backbone `H2` and `H3` takes the `NXXX`
template, a last residue carrying `OXT` takes `CXXX` -- and the current run
resolves `NMET` at residue 1 and `CGLY` at residue 76 without either name
appearing in the PDB.

This is a capability, not a tuning difference, and it can be seen in the
binaries. `strings` over each:

    template name    in 3.9.1   in CURRENT
    NGLY / CGLY         no          yes
    NPRO / CPRO         no          yes
    NCYX / CCYX         no          yes
    NTER / CTER         yes         yes      (the nucleic-acid patches)

The protein terminal vocabulary is absent from 3.9.1. `cmap_binding` and
`cmap_database_file` are likewise absent from 3.9.1 and present in CURRENT,
consistent with the gate report.

### 3b. ff19SB -- refused one step earlier, at the input file

    ===================================================================================
    Reading simulation input file sp_ubiquitin_ff19sb_vac.input
    -----------------------------------------------------------------------------------
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Wrong setting with:
    energy_term
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1, while still reading the simulation input. `energy_term` is a keyword
3.9.1 knows, but `cmap` is not a value it accepts, so it never reaches the
`.cmap` deck. Same failure as gate section 3b on the committed peptide.

Stripping the two cmap directives, as the gate report did, does not get 3.9.1
any further on this system -- it simply moves the stop to the same place ff14SB
stops:

    Total 19 atoms in 1. residue(MET) doesn't fit topology

`ff19sb/mosaics-3.9.1-nocmap/`, exit 1. So on ubiquitin the CMAP gap is not
even reachable: 3.9.1 fails on the terminus first. The 69.95 kcal/mol of CMAP
correction that sander, OpenMM and the current MOSAICS all agree on is
unrepresentable in 3.9.1, but that is inferred from the peptide result, not
measured here, because here 3.9.1 never gets to a CMAP term.

### 3c. The control: 3.9.1 cannot be handed the terminal names either

3.9.1's failure in 3a is a missing selection rule, so the obvious question is
whether naming the templates explicitly in the PDB rescues it, the way
renaming 1efs to 3.9.1's own residue vocabulary rescued the nucleic case in
gate section 3d. It does not, and the reason is structural.

The current build's own output PDB writes the resolved templates:

    ATOM      1    N NMET A   1      27.340  24.430   2.614
    ATOM   1231  OXT CGLY A  76      40.816  39.621  36.302

`NMET` and `CGLY` are four characters. The PDB residue-name field the MOSAICS
reader uses is three. Feeding that file back to either build truncates the name
in the reader, and both fail:

    3.9.1     Total 19 atoms in 1. residue(NME) doesn't fit topology
    CURRENT   NNME 0 was not found in .rtf

(The CURRENT message is its own N-prefix rule applied to the truncated `NME`.)
There is therefore no PDB that names an `NXXX` or `CXXX` template in a way
MOSAICS reads: the selection has to be internal, which is what the current
build added and 3.9.1 does not have. `<ff>/mosaics-3.9.1-control-explicit-termini/`
holds both attempts and both builds' output.

Incidentally: the current build writes a protein PDB with terminal residues
that its own reader cannot read back. That round-trip is broken in both builds.
It is not in the way of anything here -- the structure MOSAICS reads is the one
the reference topology is built from, and that file uses standard three-letter
names -- but it is a real defect and it is recorded so it is not rediscovered.

### 3d. Per force field, as asked

    force field   3.9.1 result on ubiquitin
    ---------------------------------------------------------------------------
    ff14SB        REFUSES at topology setup: "Total 19 atoms in 1. residue(MET)
                  doesn't fit topology". Reads the input, all five decks, the
                  PDB and the RTF first. No energy of any kind. Cannot select
                  the NMET template; the protein terminal template names are
                  not in the binary. Explicit naming is impossible (3c).
    ff19SB        REFUSES at the input file: "Wrong setting with: energy_term".
                  With the cmap directives stripped it gets as far as ff14SB
                  and stops at the same N-terminus. No energy of any kind.

    force field   current build result
    ---------------------------------------------------------------------------
    ff14SB        Runs. One protein molecule, chain A, 76 residues, NMET/CGLY
                  resolved, charge 0. Agrees with OpenMM to 6e-4 kcal/mol worst
                  term on a 1231-atom system.
    ff19SB        Runs. Same topology. CMAP evaluated: 16 maps loaded, 13
                  selected, 74 backbone terms, 69.95284751 kcal/mol against
                  OpenMM's 69.95284751 (difference -2.7e-13) and sander's
                  69.9528.

---

## 4. What failed, plainly

    1  MOSAICS 3.9.1, ff14SB, ubiquitin       no energy. Topology setup error at
                                              residue 1 (charged N-terminus).
    2  MOSAICS 3.9.1, ff19SB, ubiquitin       no energy. Input-file error on
                                              \energy_term{cmap}.
    3  MOSAICS 3.9.1, ff19SB, cmap stripped   no energy. Same terminus error
                                              as 1.
    4  3.9.1 with explicit NMET/CGLY names    impossible: the reader's residue
                                              field is three characters wide and
                                              the template names are four.
                                              Both builds fail on such a file.
    5  GROMACS ff19SB CMAP                    numerically wrong (61.543 against
                                              69.953) because the ParmEd route
                                              collapses 13 distinct maps onto
                                              one [ cmaptypes ] entry. Cell left
                                              blank in the table. Conversion
                                              route, not physics. F-06 at scale.
    6  MOSAICS current, ff19SB improper       +1.3e-03 kcal/mol against GROMACS,
                                              a hundred times the ff14SB
                                              residual. Recorded, cause not
                                              measured.

Nothing else failed. Both topologies built cleanly (tleap: 0 errors; the seven
warnings are five close contacts between tleap-placed hydrogens on deposited
crystal coordinates, and the two routine N-/C-terminal PDB-naming notes), the
coordinate overlay agreed to 0.0 Angstrom, and every reference engine ran to
completion.

---

## 5. Files

    build_ubiquitin_reference.py            the builder, adapted from
                                            tools/build_protein_reference.py
    deck_coverage.py                        the coverage counter of section 1a
    make_tables.py                          assembles section 2 from raw output
    tables.txt                              its output

    <ff>/topology/build_ubiquitin_<ff>.in   the tleap input, committed
    <ff>/topology/build_ubiquitin_<ff>.log  tleap's output
    <ff>/topology/ubiquitin_<ff>_leap_input.pdb   what tleap read (HIS->HIE)
    <ff>/topology/ubiquitin_<ff>.prmtop     the reference topology
    <ff>/topology/ubiquitin_<ff>_matched.rst7     coordinates rebuilt from the
                                            PDB MOSAICS reads
    <ff>/topology/ubiquitin_<ff>_mosaics.pdb      what MOSAICS reads
    <ff>/topology/build_metadata.json       counts, checks, provenance
    <ff>/topology/deck_coverage*.txt        section 1a

    <ff>/sander/sander_sp.in                settings, committed
    <ff>/sander/ubiquitin_sander_sp.out     raw output
    <ff>/openmm/openmm_sp.out               raw output, tools/openmm_single_point.py
    <ff>/openmm/openmm_vs_amber_comparison.tsv    tools/openmm_amber_single_point.py
    <ff>/gromacs/sp.mdp                     settings, committed
    <ff>/gromacs/sp.xvg                     the energies read (never the screen)
    <ff>/gromacs/sp.log                     raw output

    <ff>/mosaics-3.9.1/                     STABLE, input + output + exit code
    <ff>/mosaics-current/                   CURRENT, input + output + exit code
    <ff>/mosaics-3.9.1-control-explicit-termini/  section 3c, both builds
    ff19sb/mosaics-3.9.1-nocmap/            section 3b

The committed convention is that the GROMACS topology is regenerated rather
than kept, by the recipe in `evidence/engines.txt`. Here `sys.gro` (72K),
`sys.top` (673K ff14SB, 684K ff19SB) and `sp.tpr` (387K) are left in place
anyway, because the `[ cmaptypes ]` collapse of section 2 is a claim about the
contents of `sys.top` and a reader should be able to open the file and see it
rather than regenerate it and hope for the same result.
