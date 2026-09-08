# BPTI (PDB 5PTI) — ff14SB and ff19SB, MOSAICS 3.9.1 against the current build

Vacuum, dielectric 1, no cutoff, no temperature, no move, single point.
Every number below is read from a file committed under this directory. Where a
build produced no number, that is stated and the error it printed is quoted.

Builds under test:

    STABLE   /Users/bright/Documents/MOSAICS/work/staging/bin/mosaics-3.9.1
             built 21 July 2026, before the force-field work
    CURRENT  /Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/examples/mosaics.x
             branch bright/fix-staging-nme-rna, PR #48, 19 August 2026

Both binaries print `version.3.9.1` in their banner. They are told apart by
path, not by what they say about themselves.

---

## 1. What this system is

    5PTI chain A, 58 residues, 892 atoms with hydrogen, net charge +6
    three disulfides, all six cysteines used, no free CYS
    termini charged: NARG at residue 1, CALA at residue 58 — not capped
    no histidine, so no protonation choice

Residue counts (identical for both force fields):

    ALA 6  ARG 6  ASN 3  ASP 2  CYX 6  GLN 1  GLU 2  GLY 6  ILE 2  LEU 2
    LYS 4  MET 1  PHE 4  PRO 4  SER 1  THR 3  TYR 4  VAL 1        = 18 kinds

Structure and chemistry follow `_structures/BUILD-DECISIONS.md` exactly:
hydrogens stripped from the neutron deposit and rebuilt by tleap, altloc A
kept, all HETATM dropped, OXT dropped and rebuilt, charged termini.

### The disulfides are built, not assumed

`5PTI_A_cyx.pdb` is `_structures/5PTI_A.pdb` with all six CYS renamed CYX
(36 atom records touched, 6 residues x 6 heavy atoms). The three bonds are
made explicitly in leap — `topology/build_bpti_<ff>.in`:

    bond mol.5.SG  mol.55.SG
    bond mol.14.SG mol.38.SG
    bond mol.30.SG mol.51.SG

and confirmed present in the resulting topology before anything was run.
`topology/build_metadata.json` records the three SG–SG bonds found in the
prmtop with their lengths, which reproduce the deposited SSBOND distances:

    CYX5  – CYX55   2.0430 A     (deposit 2.04)
    CYX14 – CYX38   2.0302 A     (deposit 2.03)
    CYX30 – CYX51   2.0218 A     (deposit 2.02)

The build asserts exactly three SG–SG bonds and fails if it finds any other
number. tleap reported `Errors = 0`; its three warnings are the +6 net charge
and the two PDB-format terminal renames NARG→ARG, CALA→ALA.

The structure MOSAICS reads (`bpti_<ff>_mosaics.pdb`) is the tleap output PDB
with a chain id stamped on it and no atom renamed — BPTI carries no ACE/NME
caps, so the cap-renaming step of the committed protein builder does not
apply. The restart file the reference engines read is written *from* that PDB;
the maximum coordinate disagreement after the overlay is 0.0 A exactly, for
both force fields.

---

## 2. What the two MOSAICS builds did

| | ff14SB | ff19SB |
|---|---|---|
| **MOSAICS 3.9.1 (STABLE)** | **FAILED** — reads the deck, dies on the structure | **FAILED** — dies on the input file |
| **MOSAICS current** | ran, all terms | ran, all terms |

### 3.9.1 on ff14SB — reads the deck, cannot build the residues

3.9.1 parses the ff14SB profile with no complaint. It reports exactly the same
deck sizes as the current build (90 bond, 258 bend, 526 torsion/improper, 595
1-4, 595 vdW rows), reads the 2564-line RTF, and reads all 892 PDB lines. It
then stops (`mosaics-3.9.1/mosaics_sp.out`, exit status 1):

    ===================================================================================
    Setting up residue topology
    -----------------------------------------------------------------------------------
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Total 26 atoms in 1. residue(ARG) doesn't fit topology
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

Residue 1 is a charged N-terminal arginine: 26 atoms, N + H1 + H2 + H3. The
RTF has the template it needs, `NARG`, but 3.9.1 has no rule that turns a
residue named `ARG` carrying H2 and H3 into `NARG`, so it looks up `ARG` (24
atoms) and the atom count does not match. The current build makes that
selection and prints the result in its own system table:

     1. residue is 1  NARG  natom 26, charge  2
    58. residue is 58  CALA natom 11, charge -1
    molecular charge 6

The deck is not the obstacle — 3.9.1 read it. Whole-protein terminal chemistry
is. This is the first residue in the file, so 3.9.1 never reaches a CYX and
never reaches a disulfide on this system.

### 3.9.1 on ff19SB — rejects the input file

3.9.1 does not get as far as the deck. It stops on the first line of the
`~sim_mol_def` block that ff19SB requires (`mosaics-3.9.1/mosaics_sp.out`):

    Reading simulation input file sp_bpti_ff19sb_vac.input
    -----------------------------------------------------------------------------------
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Wrong setting with:
    energy_term
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

`\energy_term` exists in 3.9.1 — its binary carries `\energy_term{cryo_em}`
messages — but `cmap` is not one of the values it accepts. Two labelled
diagnostic runs under `mosaics-3.9.1/` peel the input back one line at a time
to put the rest of the refusal on the record:

    diagnostic_no_cmap_keyword/   \energy_term{cmap} removed
                                  -> "Key word cmap_database_file unknown in
                                      input file line: ~sim_mol_def"
    diagnostic_no_cmap_at_all/    both cmap lines removed
                                  -> "Total 26 atoms in 1. residue(ARG)
                                      doesn't fit topology"

So 3.9.1 has no CMAP at all: not the energy term, not the deck keyword, not
the file format. With both cmap lines gone it reaches the same terminal-residue
wall as ff14SB. No ff19SB energy is reported from the second diagnostic; an
ff19SB number computed with the correction map silently absent is not an ff19SB
number.

### The disulfides themselves

3.9.1 stops at residue 1, so its handling of CYX is never exercised on this
structure. What can be said from the binaries themselves: `mosaics.x` (current)
contains the strings `NCYX`, `CCYX`, `Cannot form disulfide between residues
%s and %s`, `A CYX sulfur is within disulfide distance of more than one CYX
residue` and `Allocate more disulfide exclusions`. `mosaics-3.9.1` contains
none of them — zero matches for `disulf` or `CYX`. The disulfide path
(`Pos::setup_residue_bond` pairing CYX sulfurs within 2.5 A and
`Residue_bond::setup_disulfide`) is code the stable build does not have.

The current build's own residue table shows all six CYX at the expected
positions (5, 14, 30, 38, 51, 55) and a total system charge of +6, matching
tleap.

---

## 3. Energies

kcal/mol. Electrostatics shown as printed and again corrected to one common
constant, 332.0637133 kcal A/(mol e^2), using the per-engine constants in
`evidence/engines.txt` (sander 332.052200, MOSAICS 332.063441, OpenMM and
GROMACS 332.063713, left uncorrected). Nothing else is corrected.

The GROMACS `dihedral` cell is `Proper Dih.` + `Per. Imp. Dih.` summed, because
sander, OpenMM and MOSAICS all fold AMBER impropers into one dihedral total.
The split is printed underneath each table.

### ff14SB on BPTI

    term            MOSAICS-3.9.1    MOSAICS-current             sander             OpenMM            GROMACS
    bond                  FAILED       112.06814628       112.06810000       112.06814628       112.06814627
    angle                 FAILED       389.22411959       389.22460000       389.22463044       389.22463026
    dihedral              FAILED       628.40110371       628.40110000       628.40119305       628.40119097
    lj                    FAILED       -46.94246505       -46.94250000       -46.94246538       -46.94248757
    coulomb               FAILED     -1415.06190298     -1415.01410000     -1415.06306700     -1415.06307003
                (electrostatic values above are as printed)
    coulomb*              FAILED     -1415.06306336     -1415.06316301     -1415.06306700     -1415.06307003
    total                 FAILED      -332.31099845      -332.31186301      -332.31156262      -332.31159011

    MOSAICS-current minus OpenMM, after correction:
      bond       -7.9581e-13
      angle      -5.1085e-04
      dihedral   -8.9334e-05
      lj         +3.3351e-07
      coulomb    +3.6317e-06

    MOSAICS-3.9.1 minus OpenMM: there is no number. See section 2.

    sander minus OpenMM, after correction      GROMACS minus OpenMM
      bond       -4.6278e-05                     bond       -6.8125e-09
      angle      -3.0436e-05                     angle      -1.7755e-07
      dihedral   -9.3046e-05                     dihedral   -2.0803e-06
      lj         -3.4617e-05                     lj         -2.2189e-05
      coulomb    -9.6018e-05                     coulomb    -3.0321e-06

    dihedral split      MOSAICS-current      GROMACS
      proper              622.797188667    622.797287763
      improper              5.603915045      5.603903203

    files  ff14sb/mosaics-3.9.1/mosaics_sp.out     (the failure)
           ff14sb/mosaics-current/mosaics_sp.out
           ff14sb/sander/bpti_sander_sp.out
           ff14sb/openmm/openmm_sp.out
           ff14sb/gromacs/sp.xvg

sander's residuals are all at its own printing precision (four decimals, about
2e-4 on numbers this size), so sander is not in disagreement with anything
here. Every term of the current build matches OpenMM to better than 1e-3, the
worst being the angle term at -5.1e-4 on a 389 kcal/mol total — the same
small angle-term residual the committed evidence records on 1efs (-1.8e-4 on
221 kcal/mol) and on the peptide (-9.6e-6 on 0.99 kcal/mol), roughly in
proportion to the number of angles.

### ff19SB on BPTI

    term            MOSAICS-3.9.1    MOSAICS-current             sander             OpenMM            GROMACS
    bond                  FAILED       112.06814628       112.06810000       112.06814628       112.06814627
    angle                 FAILED       389.22411959       389.22460000       389.22463044       389.22463026
    dihedral              FAILED       230.48181304       230.52390000       230.52389269       230.52389054
    cmap                  FAILED        43.18683042        43.18680000        43.18683042      -59.28850741
    lj                    FAILED       -48.65177437       -48.65180000       -48.65177467       -48.65179708
    coulomb               FAILED     -1415.06190298     -1415.01410000     -1415.06306700     -1415.06307003
                (electrostatic values above are as printed)
    coulomb*              FAILED     -1415.06306336     -1415.06316301     -1415.06306700     -1415.06307003
    total                 FAILED      -688.75276802      -688.71156301      -688.71134184      -791.18670746

    MOSAICS-current minus OpenMM, after correction:
      bond       -7.9581e-13
      angle      -5.1085e-04
      dihedral   -4.2080e-02      <-- see below
      cmap       -5.7554e-13
      lj         +3.0155e-07
      coulomb    +3.6317e-06

    MOSAICS-3.9.1 minus OpenMM: there is no number. See section 2.

    sander minus OpenMM, after correction      GROMACS minus OpenMM
      bond       -4.6278e-05                     bond       -6.8125e-09
      angle      -3.0436e-05                     angle      -1.7755e-07
      dihedral   +7.3095e-06                     dihedral   -2.1551e-06
      cmap       -3.0423e-05                     cmap       -1.0248e+02   <-- see below
      lj         -2.5326e-05                     lj         -2.2410e-05
      coulomb    -9.6018e-05                     coulomb    -3.0321e-06

    dihedral split      MOSAICS-current      GROMACS
      proper              224.973266426    224.973306166
      improper              5.508546613      5.550584369

    files  ff19sb/mosaics-3.9.1/mosaics_sp.out     (the failure)
           ff19sb/mosaics-current/mosaics_sp.out
           ff19sb/sander/bpti_sander_sp.out
           ff19sb/openmm/openmm_sp.out
           ff19sb/gromacs/sp.xvg

The correction map matches OpenMM to 5.8e-13 and sander to its printing
precision. CMAP on a real protein — 56 maps, 13 of the deck's 16 distinct maps
— is reproduced exactly by the current build.

---

## 4. Two disagreements that are not about the MOSAICS build

### The ff19SB improper term, current build, -0.0420 kcal/mol

The ff19SB dihedral residual above is not spread across the torsions. Splitting
it (the table's last block) puts essentially all of it in the impropers:

    proper    MOSAICS-current 224.973266426  vs GROMACS 224.973306166   -4.0e-05
    improper  MOSAICS-current   5.508546613  vs GROMACS   5.550584369   -4.20e-02

On the same structure with ff14SB the improper term agrees to +1.2e-05. So the
defect is specific to the ff19SB deck, and it is confined to one term worth
0.76% of that term and 0.006% of the system's total energy. sander, OpenMM and
GROMACS all agree with each other on the ff19SB impropers, so the target is not
in doubt.

`ff19sb/improper_decomposition.txt` computes every improper of the reference
topology directly from the prmtop with E = k(1 + cos(n*phi - phase)),
independently of all four engines, and groups them by atom-type quadruple. What
that shows is that ff19SB splits the alpha-carbon type: `XC` on the 56 residues
that carry a CMAP, `CX` retained on the two charged termini (ARG 1 and ALA 58,
the only two residues in this structure with no CMAP). That split produces
backbone improper quadruples the capped committed peptide can never produce —
one `C-XC-N-H`, one `CX-N-C-O` and one `C-CX-N-H` at the CX/XC boundary,
against 56 `XC-N-C-O` and 51 `C-H-N-XC` in the interior. The ff19SB improper
rows in the deck are keyed `XC-N-C-O`, `C-H-N-XC`, `C-CX-N-H`, `CX-N-C-O`,
`C-CT-N-XC`, `N-X-CT-XC`. Which of the boundary quadruples MOSAICS resolves
differently was not isolated further; what is measured is the term-level
residual and its confinement to the impropers.

This is a finding of this system. The committed protein evidence is an
ACE/NME-capped peptide, which has no charged termini and therefore no CX/XC
boundary at all, and its ff19SB dihedral residual is +1.3e-06.

### The GROMACS CMAP number is wrong, and this time it is not blank

The committed evidence (`evidence/results.txt`, F-06) leaves the GROMACS CMAP
cell for the peptide blank, because the ParmEd conversion route keys
`[ cmaptypes ]` on the five backbone atom types spanning phi and psi, all four
peptide residues share those types, the four maps collapse onto one, and the
conversion keeps only the first.

On BPTI the same route does not go blank. It produces a number:

    GROMACS CMAP Dih.   -59.28850741 kcal/mol
    sander CMAP          43.18680000
    OpenMM CMAP          43.18683042
    MOSAICS CMAP         43.18683042

    GROMACS minus OpenMM  -102.48 kcal/mol

`sys.top` contains `[ cmaptypes ]` with 16 entries and a `[ cmap ]` section, so
nothing was dropped for GROMACS to notice; the maps are simply keyed on
backbone atom types that cannot tell 13 different residue maps apart, and 56
CMAP terms were evaluated against whichever map the key resolved to. All three
other engines agree to 3e-05. **The GROMACS CMAP cell for BPTI is a wrong
number produced by the conversion route, not a disagreement about physics, and
it should not be quoted as a GROMACS ff19SB CMAP energy.** It is recorded here
because a silently wrong number is more dangerous than a blank one, and the
peptide could not show it: with one residue type there was nothing to collapse
onto, so the failure showed up as a missing cell instead of a wrong cell.

Every other GROMACS value on this structure agrees with OpenMM to 2e-05 or
better, in both force fields.

---

## 5. Composition and coverage

The referee's objection is that the committed systems exercise an unstated
fraction of each deck. What BPTI touches, counted from the reference topology
and matched against the committed decks row by row (`tools/coverage.py`,
`<ff>/topology/deck_coverage.json`):

    quantity                       ff14SB    ff19SB
    atoms                             892       892
    residues                           58        58    (18 kinds)
    net charge                         +6        +6
    distinct AMBER atom types          21        22    (ff19SB adds XC)
    bonds                             906       906
    angles                           1626      1626
    proper torsion terms             3791      3791
    improper torsion terms            199       199
    1-4 pairs                        2347      2347
    CMAP terms                          0        56

    DECK ROWS TOUCHED             ff14SB              ff19SB
    bond                       48 /  90  53.3%     53 / 102  52.0%
    bend                      117 / 258  45.3%    131 / 329  39.8%
    torsion, proper           221 / 465  47.5%    233 / 616  37.8%
    torsion, improper          15 /  61  24.6%     17 /  74  23.0%
    Lennard-Jones pairs       231 / 595  38.8%    253 / 666  38.0%
    1-4 pairs                 103 / 595  17.3%    108 / 666  16.2%
    CMAP maps                       n/a           13 /  16  81.2%
    CMAP residue bindings           n/a           18 /  28  64.3%

The row counts in the denominator are the deck sizes MOSAICS itself prints as
it reads each file, so they are the engine's own numbers. The bond and bend
figures in the numerator are confirmed independently by sander, which prints
the count of distinct parameter rows in the topology it read: `NUMBND = 48,
NUMANG = 117` for ff14SB and `NUMBND = 53, NUMANG = 131` for ff19SB — the same
48/117 and 53/131. The CMAP figure is confirmed twice over: 13 distinct maps
from the prmtop, and 13 distinct `map_id` values from the deck's
`~cmap_binding` rows for the 18 residue kinds present.

Torsion rows are matched the way MOSAICS matches them, exact rows taking
precedence over wildcard (`X`) rows, and counted in either direction.

One protein of 58 residues therefore exercises between a sixth and a half of
each AMBER protein deck, and four fifths of the ff19SB correction maps. What it
does not exercise: HID/HIE/HIP, TRP, ASH/GLH/LYN/CYM, HYP, NHE, the ACE/NME
caps, and every terminal template other than NARG and CALA.

---

## 6. What failed, and why

| row | outcome |
|---|---|
| ff14SB, MOSAICS 3.9.1 | **failed.** `Total 26 atoms in 1. residue(ARG) doesn't fit topology`. 3.9.1 reads the ff14SB deck and the RTF without error, but has no rule mapping a charged N-terminal residue onto the `NARG` template. No energy. |
| ff19SB, MOSAICS 3.9.1 | **failed.** `Wrong setting with: energy_term`. `cmap` is not an accepted value in 3.9.1. With that line removed it then rejects `cmap_database_file` as an unknown keyword; with both removed it hits the same terminal-residue error. No energy. |
| disulfides, MOSAICS 3.9.1 | **not reached.** 3.9.1 stops at residue 1. Its binary contains no CYX or disulfide code at all, so it has no path to form an SG–SG bond, but on this structure that is inferred from the binary, not measured from a run. |
| ff14SB, MOSAICS current | ran. All five terms within 5.2e-4 of OpenMM. |
| ff19SB, MOSAICS current | ran. Bond, angle, CMAP, LJ and electrostatics within 5.2e-4; **dihedral -0.0420, all of it in the impropers**, at the ff19SB CX/XC type boundary created by the charged termini. |
| ff19SB, GROMACS | ran, and produced a **wrong CMAP energy**, -59.2885 against +43.1868 from the other three. ParmEd's `[ cmaptypes ]` keying cannot tell 13 residue maps apart. Every other GROMACS term on this structure is right to 2e-05. |
| ff14SB / ff19SB, sander | ran. All residuals at its four-decimal printing precision. |
| ff14SB / ff19SB, OpenMM | ran. Used as the reference column. |

---

## 7. Files

    5PTI_A_cyx.pdb                 the cleaned deposit with the six CYX renamed
    energies.json                  every number in section 3, machine readable
    tools/build_bpti.py            the builder (tleap + explicit disulfides + coordinate overlay)
    tools/coverage.py              the deck-row coverage count
    tools/tables.py                the table generator
    <ff>/energy_table.txt          section 3 as written by the generator
    <ff>/improper_decomposition.txt  every improper of the reference topology, by type
    <ff>/topology/                 tleap input, log, prmtop, rst7, matched rst7,
                                   the PDB MOSAICS reads, build_metadata.json,
                                   deck_coverage.json
    <ff>/sander/                   sander_sp.in and bpti_sander_sp.out
    <ff>/openmm/                   openmm_sp.out
    <ff>/gromacs/                  sp.mdp, sp.log, sp.xvg, sp.edr, sys.gro, sys.top,
                                   grompp and mdrun output, box_note.txt
    <ff>/mosaics-3.9.1/            the input, and the output in which it fails
    <ff>/mosaics-current/          the input, and the run that succeeded
    ff19sb/mosaics-3.9.1/diagnostic_no_cmap_keyword/    labelled, not a result row
    ff19sb/mosaics-3.9.1/diagnostic_no_cmap_at_all/     labelled, not a result row

Environment: sander 24.0, OpenMM 8.5.2 Reference platform, GROMACS 2025.4
double precision (`gmx_d`), ParmEd 4.3.1, tleap from AmberTools 24, all from
`/opt/homebrew/Caskroom/miniforge/base/envs/mosaics-ff-local`. Decks and tools
from the CURRENT worktree, `/Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna`.
