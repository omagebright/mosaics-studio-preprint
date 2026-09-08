# Gate report -- evidence-v2

Run 2026-08-25. Everything below was produced by running the two binaries.
Every number quoted here is in a file under `evidence-v2/_gate/`; nothing is
carried across from another system and nothing is estimated.

    STABLE   /Users/bright/Documents/MOSAICS/work/staging/bin/mosaics-3.9.1
             md5 6c268d2866adbcf189f03fa8afaa1b51   built 21 Jul 2026
    CURRENT  /Users/bright/Documents/MOSAICS/work/pr-fix-staging-nme-rna/examples/mosaics.x
             md5 1ba6ddfdd4e59e29a54bb15319c18f01   built 17 Aug 2026

Verdict: **PASS**, with two findings the authors need before anything
downstream is written up. Neither blocks the work.

---

## 1. Both binaries start

Run with no argument, each prints its banner, reports `Cannot open input file
(null)` and stops. Both start.

    STABLE   banner reads  >>  version.3.9.1  <<
    CURRENT  banner reads  >>  version.3.9.1  <<

**THE TWO BUILDS PRINT THE SAME VERSION STRING.** The banner is not evidence
of which binary produced an output file, and a reader of these outputs cannot
tell them apart by reading them. They are different binaries -- different md5,
different behaviour throughout this report -- but the version was never bumped
past 3.9.1. Anything downstream must record the path and md5, not the banner.
Worth a one-line fix before the force-field stack merges.

Raw: `_gate/run/*/mosaics_sp.out`

---

## 2. CURRENT reproduces the committed numbers -- exactly, and on all six

The gate asked for 1efs OL21/OL3. That reproduces, and so does every other
committed MOSAICS single point. `diff` of the printed energy block against the
committed file is empty in all six cases; on OL21/OL3 the *whole output file*
is identical apart from the timestamp line.

    system    force field   committed file                                    result
    1efs      OL21/OL3      1efs/ol21_ol3/mosaics/mosaics_sp.out              identical
    1efs      OL15/OL3      1efs/ol15_ol3/mosaics/mosaics_sp.out              identical
    1efs      parmbsc1      1efs/bsc1/mosaics/mosaics_sp.out                  identical
    1efs      parmbsc0      1efs/bs0/mosaics/mosaics_sp.out                   identical
    peptide   ff14SB        protein/ff14sb/mosaics/mosaics_sp.out             identical
    peptide   ff19SB        protein/ff19sb/mosaics/mosaics_sp.out             identical

The OL21/OL3 terms the gate named, as CURRENT prints them:

    term       asked for          CURRENT prints
    bond       70.99877741        70.99877740744154
    angle      220.56489609       220.5648960933117
    dihedral   596.17955421       596.1795542133353
    lj        -219.29643681      -219.2964368100409
    coulomb   1857.73806648      1857.738066483539

(The committed table's `total` of 2526.18637965 is the corrected-constant sum,
not a number MOSAICS prints; MOSAICS prints 2526.184857387645 using its own
constant, as the committed file also does. Same number, same convention.)

**PR #48 changed no energy.** Checked directly rather than inferred: a
recursive diff of `source/` between `stack4-rewrite` (PR #47, the branch the
paper's numbers came from) and `pr-fix-staging-nme-rna` (PR #48) returns
exactly one differing file,

    source/transform/stage/stage_defs/stage_defs_AA_protein.cpp   19 lines

and the change is confined to resolving ACE/NME methyl-hydrogen and
methyl-carbon atom indices (`H1/H2/H3` vs `HH31/HH32/HH33`; NME writing its
methyl carbon as `C`) inside the *staging* transform. Nothing in an energy
term, nothing on the nucleic path. Consistent with the six identical outputs.

Raw: `_gate/<forcefield>/current/mosaics_sp.out`

---

## 3. What 3.9.1 can and cannot do -- the "before" column

This is the substance. Four force fields 3.9.1 refuses outright; the rest it
runs and gets wrong. **In no case does 3.9.1 reproduce a committed number.**

### 3a. Refused at the deck: OL21/OL3, parmbsc1, parmbsc0

    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Key word mol_type unknown in database line:
    ~torsion_parm
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1, while reading `mosaics_*.tors_and_impr`. Nothing further is attempted:
no topology, no energy, no output.

`\mol_type{dna}` / `\mol_type{rna}` is how the deck says a torsion row belongs
to one sugar and not the other -- the DNA/RNA discrimination 1efs exists to
test. It appears 231 times in the OL21/OL3 torsion deck (122 dna, 109 rna),
207 times in parmbsc1, 191 in parmbsc0, and alongside it `\res_context{IIII}`.

Neither string occurs anywhere in the 3.9.1 binary:

    keyword              in 3.9.1   in CURRENT
    mol_type                no          yes
    res_context             no          yes
    cmap_binding            no          yes
    cmap_database_file      no          yes

So this is not a version 3.9.1 can be talked into. The vocabulary a deck needs
to distinguish DNA from RNA torsions did not exist in it.

Raw: `_gate/{ol21_ol3,bsc1,bs0}/stable/mosaics_sp.out`

### 3b. Refused at the input file: ff19SB

    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@
    Wrong setting with:
    energy_term
    @@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@

exit 1, while still reading the simulation input. `energy_term` exists as a
keyword in 3.9.1, but `cmap` is not one of the values it accepts, so it never
reaches the `.cmap` deck.

**Deleting the two cmap directives makes it run, and the term is then absent
rather than zero.** With `\energy_term{cmap}` and `\cmap_database_file{...}`
stripped, 3.9.1 loads the ff19SB deck and completes. Its report has no
`-Cmap energy` line at all -- not a zero, no line. CURRENT prints
`-Cmap energy -1.921240056940442`, matching the committed value and sander's
-1.92120000 and OpenMM's -1.92124000. So the CMAP correction on this peptide
is 1.92 kcal/mol of physics that 3.9.1 cannot represent.

Its torsion term, notably, is right: 6.885498028432928 against CURRENT's
6.885498028432925. It reads the ff19SB torsion deck perfectly well. What it
cannot do is the correction map, and what it cannot do besides is 3c.

Raw: `_gate/ff19sb/stable/mosaics_sp.out`, `_gate/ff19sb/stable_nocmap/mosaics_sp.out`

### 3c. Runs, and returns wrong energies: OL15/OL3, ff14SB, CHARMM36

Where no new keyword is in the way, 3.9.1 completes -- and builds the wrong
molecule. **It does not connect residues into a chain.**

    system / force field     3.9.1 builds                CURRENT builds
    1efs OL15/OL3            26 separate molecules       2 chains of 13
    peptide ff14SB           3 molecules, 1 biopolymer   1 chain of 6
                             (ACE and NME detached as
                              "other molecules")
    rna_acgu CHARMM36        4 separate molecules        1 chain of 4
    dna_acgt CHARMM36        4 separate molecules        1 chain of 4

Every inter-residue bond, bend and 1-4 is therefore missing, and every pair
that should have been excluded is counted as a full non-bonded contact. The
Lennard-Jones term is the symptom:

    system / force field    3.9.1 total LJ      CURRENT total LJ    committed
    1efs OL15/OL3           +134453.4459599     -219.4578376(*)     -219.29643683
    peptide ff14SB           +55142.5978396      134.5056784         134.50567845
    rna_acgu CHARMM36         +5586.6164456     2795.6454518        2795.6454518
    dna_acgt CHARMM36         +5311.1190227      446.1783543        (see finding B)

The ff14SB row is the cleanest statement of the defect in the whole set,
because there is nothing else it could be. The deck files, the RTF, the PDB
and the input are byte-identical for both builds -- `staging/forcefields/
ff14sb_openmm/*` and `params/db_ff14sb/*` diff clean on all five deck files
and on `top_openmm-ff14sb_protein.rtf`. Same files, same command, and:

    term        3.9.1                      CURRENT / committed
    bond          0.1474592554762238         0.1474918101971024
    angle         0.9872532462321912         0.9873252507945278
    dihedral     25.12266535615493          31.98193724561154
    lj        55142.59783957601            134.5056784450612
    coulomb     -58.45457876604456          -71.35989833880687
    TOTAL     55110.40063866777             96.2625344128578

55,014 kcal/mol, because ACE and NME are not attached to the peptide. This is
the cap defect: 3.9.1 detaches exactly the two residues that have one backbone
branch instead of two.

Raw: `_gate/{ol15_ol3,ff14sb,charmm36_rna}/stable/mosaics_sp.out`,
`_gate/charmm36_dna_oligo/stable/mosaics_sp.out`

### 3d. The honest control: 3.9.1 given its own naming still gets it wrong

3.9.1's failure to build a chain in 3c is partly a naming change, and the
paper must not claim more than it can. The panel RTFs were renamed in this
work: 3.9.1's shipped RTF calls RNA residues ADE/URA/CYT/GUA and DNA
ADD/THD/CYD/GUD, and the current RTF calls them RA/RU/RC/RG and DA/DT/DC/DG.
Apart from those 24 `RESI` names and 4 comment lines, the two OL15/OL3 RTFs
are identical -- same charges, same atom types, same bonds, same CHIGROUP.
The decks are identical too, save the 1-4 scale being written to more digits
(0.8333333333 against 0.8333333333333334).

So 1efs was renamed to 3.9.1's convention and rerun on 3.9.1's own files. It
then builds the chain -- 2 molecules of 13 -- and still returns wrong energies:

    term       3.9.1, its own naming    CURRENT / committed OL15/OL3   error
    bond          71.04187731389734        70.99877741121985        +0.0431
    angle        220.9686005245711        220.5648960758804         +0.4037
    dihedral     592.8214352955542        592.8489526607921        -0.0275
    lj          -219.4152430493324       -219.2964368273058        -0.1188
    coulomb     1858.327568345644        1857.738066483539         +0.5895
    TOTAL       2523.744764277621        2522.854781651514         +0.8900

For scale: sander, OpenMM and GROMACS agree with the committed MOSAICS number
to about 1e-04 on every term. 3.9.1 is out by 0.89 kcal/mol on the total and
0.59 on electrostatics -- four orders of magnitude outside the cross-engine
spread.

Running CURRENT on those same 3.9.1-era files separates the two causes. It
gets bond right (70.99877741122013, to 12 digits) but angle wrong
(220.957399169732), and a total of 2523.602055162067. So:

  - the **bond** term is 3.9.1's own defect: same files, both builds, and only
    CURRENT gets 70.9988;
  - the **angle** term is the deck: neither build recovers 220.5649 from the
    old residue names, so the RTF rename carries real chemistry and is not
    cosmetic.

Neither route lets 3.9.1 reach the committed answer. That is the claim the
evidence supports, and it is the claim to make.

Raw: `_gate/_control_391_native_naming/ctl_sp.out` (3.9.1),
`_gate/_control_391_native_naming/ctl_sp_CURRENT.out` (CURRENT, same files)

### 3e. What 3.9.1 does get right

Worth saying, so the "before" column is not a caricature. On the CHARMM36
nucleic deck 3.9.1 parses everything and evaluates both terms that
`evidence/charmm36.txt` said needed no engine change: Urey-Bradley comes back
non-zero in the bend's own bond slot (`Bond-Bend energy 24.60066542625691` on
rna_acgu) and the harmonic impropers come back
(`Tors-Improper 0.004089189089420226`). C-02 and C-03 are confirmed against
the old binary. The ff19SB torsion term is exact to 15 digits. The parser and
the term algebra are largely fine; the topology and the new deck vocabulary
are not.

### Per force field, as asked

    force field       3.9.1 result
    ---------------------------------------------------------------------------
    ff14SB            runs; ACE and NME detached; total 55110.401 against
                      96.2625. Byte-identical deck and RTF, so this is the
                      engine alone.
    ff19SB            REFUSES the input: "Wrong setting with: energy_term".
                      With cmap stripped it runs and prints NO Cmap line --
                      the term is absent, not zero. 1.92 kcal/mol unrepresent-
                      able. Caps also detached, total 55092.163 against
                      69.2449.
    parmbsc0          REFUSES the deck: "Key word mol_type unknown". 191
                      mol_type rows. No energy of any kind.
    parmbsc1          REFUSES the deck: same error. 207 mol_type rows.
    OL15/OL3          runs on the new deck; 26 detached residues, total
                      137879.984 against 2522.856. On its own 3.9.1-era files
                      it builds the chain and is still out by 0.89 kcal/mol
                      (3d).
    OL21/OL3          REFUSES the deck: same error. 231 mol_type rows. This
                      is the force field the whole 1efs comparison is built
                      on, and 3.9.1 cannot load it at all.
    CHARMM36 nucleic  runs; parses the whole deck; evaluates Urey-Bradley and
                      harmonic impropers correctly; 4 detached residues.
                      rna_acgu total 5852.768 against 3011.907.
                      dna_acgt total 5708.053 against 849.722.
                      CMAP: no Cmap line, as in ff19SB.
    CHARMM36 protein  not attempted, and cannot be: no protein CHARMM36 deck
                      exists in either tree. C-04 stands -- ~cmap_binding
                      cannot express a selector conditioned on the following
                      residue's class.

---

## Findings for the authors

### Finding A -- the two builds report the same version string
Both print `>> version.3.9.1 <<`. An output file cannot be attributed to a
build by reading it. Low effort, worth fixing before the stack merges.

### Finding B -- the committed CHARMM36 DNA row does not reproduce, and PR #48 is not why

`evidence/charmm36_dna_acgt.tsv`, `evidence/charmm36.txt` C-11 and
`params/db_charmm36_na/manifest.tsv` all carry

    mosaics_total_kcal_dna   1954.691252416373
    reference_structure_dna  work/charmm36-na/oligos/dna_acgt.pdb

CURRENT, run on that exact structure with that exact deck, gives
**849.7216661162275**, and disagrees on every term:

    term                committed          CURRENT (this run)
    bond                 18.92184908593162    53.21465006384997
    Urey-Bradley        167.6513797823571    157.6823163369957
    angle               412.7880943139457    423.4043055570307
    torsion             122.9466154568943    122.8276067764581
    improper              0.001025553235209    0.003590713275591
    1-4 Lennard-Jones  1412.271400982297     276.5865011018748
    Lennard-Jones       176.6882040551721    169.5918531581951
    Coulomb            -356.5773168134589   -353.5891575914526
    TOTAL              1954.691252416373     849.7216661162275

What has been ruled out:

  - **The engine.** PR #47 and PR #48 differ in one source file, protein
    staging only (section 2). It cannot move a CHARMM36 DNA single point.
  - **The deck.** All six files of `params/db_charmm36_na/` plus
    `top_charmm36_nucleic.rtf` are byte-identical between the two worktrees,
    manifest included.
  - **The structure.** `work/charmm36-na/oligos/dna_acgt.pdb` and
    `examples/pure_terminal_controls/dna_acgt.pdb` have identical coordinates
    (same md5 over columns 31-54); the 48 lines that differ are two-letter
    atom names justified differently in the column. CURRENT returns
    849.7216661162275 for both.
  - **The harness.** The RNA half of the same TSV reproduces *exactly*
    (3011.907304997552, all 16 digits) through the identical path, and the
    committed RNA run `stack4-rewrite/work/charmm36-na/rna4-compare/` confirms
    the protocol. Whatever is wrong is specific to the DNA row.
  - **Topology.** CURRENT builds one chain, DA5 30 / DC 30 / DG 33 / DT3 33 =
    126 atoms, charge -3. Sane.

There is no committed run output for the DNA case anywhere on disk -- there is
a `rna4-compare/` directory with `out.pot_energy` and `openmm.out`, and no
`dna4-compare/`. The number went into the TSV, the manifest and charmm36.txt
without its output file being kept, so it cannot be audited, only re-run.
The likeliest reading is that it was produced against an earlier state of the
converter or an earlier `dna_acgt.pdb` that was later overwritten (all the
oligo files carry one mtime, 13 Aug 2026 05:41:20). C-11's own remark that
"the DNA 1-4 Lennard-Jones is 1412 kcal/mol because the oligo is AMBER leap
geometry" points the same way: 1412 kcal/mol of 1-4 LJ is a clashing
structure, and today's file gives 276.6.

**This does not affect the gate's pass.** It touches one row of the CHARMM36
supplementary table, not the OL21/OL3, OL15/OL3, bsc1, bsc0, ff14SB or ff19SB
comparisons, all six of which reproduce to the last digit. But the row as
committed cannot be defended, and either the number or the provenance has to
be corrected before it is published.

---

## Structures fetched

Downloaded 2026-08-25 from `https://files.rcsb.org/download/XXXX.pdb` into
`evidence-v2/_structures/raw/`; cleaned copies beside them. Full composition,
protonation, disulfide and terminus decisions are in
`evidence-v2/_structures/BUILD-DECISIONS.md`, which is the file downstream
agents follow.

    1UBQ   ubiquitin, chain A, 76 res, 601 heavy atoms.  X-ray, one model,
           no gaps, no non-standard residues, no missing heavy atoms, no
           altlocs, 58 waters dropped. One histidine (HIS 68) -> HIE.
           LEU73/ARG74/GLY75/GLY76 are deposited at occupancy 0.45/0.45/
           0.25/0.25 -- all atoms present, coordinates taken as deposited.

    5PTI   BPTI, chain A, 58 res, 453 heavy atoms.  Joint X-ray/neutron, one
           model, no gaps, no non-standard residues, no missing heavy atoms.
           Dropped 189 DOD, 5 PO4, 1 UNX, and all 103 D + 344 H.  Altlocs at
           GLU 7 and MET 52 only; conformer A kept (9 atoms dropped, 462 ->
           453).  Three disulfides, exactly the deposited SSBONDs:
           CYS 5-CYS 55 (2.04 A), CYS 14-CYS 38 (2.03 A), CYS 30-CYS 51
           (2.02 A).  All six cysteines are in a bridge; none free.

    1BNA   Drew-Dickerson dodecamer, chains A (1-12) and B (13-24),
           CGCGAATTCGCG each, 24 res, 486 heavy atoms.  X-ray, one model, no
           gaps, no non-standard residues, 80 waters dropped.  BOTH CHAINS
           ARE KEPT -- it is a self-complementary duplex and one strand is
           not the dodecamer.  DC 1 and DC 13 carry no P/OP1/OP2 in the
           deposit: they are already free 5'-hydroxyls, so nothing is added
           or removed.

    2KOC   UUCG tetraloop hairpin, chain A, MODEL 1 of 20 (solution NMR),
           14 res, 294 heavy atoms after cleaning.  GGCACUUCGGUGCC -- A 1,
           C 5, G 5, U 3, all four RNA bases.  No gaps, no non-standard
           residues, no missing heavy atoms, no altlocs, no HETATM at all.

           Chosen over 1F7Y because 1F7Y is a ribosomal protein-RNA complex:
           its 57-nt RNA is chain B, and the file also carries an 89-residue
           protein chain with 32 MSE atoms, 9 MG, 1 K, 2 NA, 27 waters, 67
           alternate-location pairs and 10 residues missing from the model.
           2KOC is the tetraloop and nothing else.

           Its G1 is deposited with a bare 5'-terminal phosphate (OP3, P,
           OP1, OP2).  Those four atoms are removed, giving a 5'-OH.  This is
           forced, not preferred: evidence/charmm36.txt C-08 measured that a
           bare 5'-phosphate is the one terminus charmm36_2024.xml cannot
           express, which is why 1efs has no CHARMM36 reference.  Keeping it
           would have bought a second system CHARMM36 cannot run.  The
           committed nucleic controls are 5TER/3TER for the same reason, so
           2KOC now matches them and matches 1BNA.

---

## Files written

    _gate/GATE-REPORT.md                                  this file
    _gate/params                                          symlink to the CURRENT worktree params
    _gate/run/{stable,current}/                           item-1 startup and the 1efs OL21/OL3 rerun
    _gate/<ff>/{stable,current}/mosaics_sp.out            full matrix, both builds
    _gate/ff19sb/stable_nocmap/                           3.9.1 on ff19SB with the cmap directives removed
    _gate/charmm36_dna_oligo/{stable,current}/            CHARMM36 DNA on the manifest's own reference structure
    _gate/_control_391_native_naming/                     3.9.1 and CURRENT on 3.9.1-era RTF and deck
    _structures/raw/*.pdb                                 unmodified RCSB downloads
    _structures/{1UBQ_A,5PTI_A,1BNA_AB,2KOC_A_model1}.pdb cleaned
    _structures/BUILD-DECISIONS.md                        chemistry every downstream agent follows
    _structures/{clean_structures.py,survey_structures.py} the two scripts, so both are reproducible
