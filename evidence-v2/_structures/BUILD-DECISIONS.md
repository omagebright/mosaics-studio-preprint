# Build decisions for the evidence-v2 systems

Every downstream agent follows these exactly. Consistency across the four
systems matters more than any individual choice; where a choice is arguable
the reason is given so it can be argued with rather than guessed at.

Cleaned files live beside this note. The script that produced them is
`clean_structures.py`; the script that produced the composition tables is
`survey_structures.py`. Raw downloads are in `raw/`.

    curl https://files.rcsb.org/download/XXXX.pdb        fetched 2026-08-25

## Rules applied to all four

1. Only `ATOM` records survive. Every `HETATM` is dropped: waters (HOH, DOD),
   the 5 PO4 and 1 UNX of 5PTI, all ions. No ligands are kept.
2. Alternate locations: keep altloc `A` or blank, drop the rest, and blank the
   altloc column in the output.
3. All hydrogen and deuterium is stripped, whatever the source structure
   carried, and rebuilt by the topology builder. Hydrogen placement is
   therefore tleap's (`build_protein_reference.py` /
   `build_1efs_amber_reference.py`) for the AMBER decks and OpenMM's template
   for the CHARMM path -- never the depositor's. This is the only way the four
   systems get hydrogens from one source.
4. `OXT` is dropped and let the builder add its own C-terminal oxygens.
5. Protein termini are CHARGED, not capped: NH3+ / COO-. The committed peptide
   is ACE/NME-capped because it is a fragment; a whole protein is not a
   fragment and capping 1UBQ or 5PTI would misrepresent it. Both proteins get
   the same treatment.
6. Nucleic termini are 5'-OH and 3'-OH -- see the 2KOC entry for why this is
   forced rather than preferred.

## 1UBQ -- ubiquitin

    chain A, 76 residues, 601 heavy atoms          1UBQ_A.pdb
    X-ray 1.8 A, single model, no gaps, no non-standard residues,
    no missing heavy atoms, no altlocs, 58 waters dropped

    ALA 2  ARG 4  ASN 2  ASP 5  GLN 6  GLU 6  GLY 6  HIS 1  ILE 7  LEU 9
    LYS 7  MET 1  PHE 2  PRO 3  SER 3  THR 7  TYR 1  VAL 4

  HIS PROTONATION. One histidine, HIS 68. Build it as HIE (neutral, proton on
  NE2). Reasons: it is the AMBER default for a solvent-exposed histidine, it
  is the tautomer tleap picks when asked for HIS, and at pH 7 the neutral form
  dominates. HID and HIP are both wrong for this site and, more to the point,
  would make 1UBQ the only system whose protonation was chosen by hand.

  DISORDERED TAIL. LEU 73, ARG 74, GLY 75 and GLY 76 are deposited at
  occupancy 0.45 / 0.45 / 0.25 / 0.25. All heavy atoms are present, so nothing
  is rebuilt; the coordinates are taken as deposited. Note it in the paper --
  a single-point energy on a 0.25-occupancy tail is a real number about a
  poorly determined piece of the model.

  Net charge as built: 0 (7 LYS + 4 ARG + 1 N-terminus against 5 ASP + 6 GLU
  + 1 C-terminus).

## 5PTI -- BPTI, three disulfides

    chain A, 58 residues, 453 heavy atoms          5PTI_A.pdb
    joint X-ray / neutron, single model, no gaps, no non-standard residues,
    no missing heavy atoms
    dropped: 189 DOD, 5 PO4, 1 UNX, 103 D and 344 H

    ALA 6  ARG 6  ASN 3  ASP 2  CYS 6  GLN 1  GLU 2  GLY 6  ILE 2  LEU 2
    LYS 4  MET 1  PHE 4  PRO 4  SER 1  THR 3  TYR 4  VAL 1

  DISULFIDES. Three, exactly as the deposited SSBOND records give them, and
  all six cysteines are in one:

      CYS  5 - CYS 55     S-S 2.04 A
      CYS 14 - CYS 38     S-S 2.03 A
      CYS 30 - CYS 51     S-S 2.02 A

  All six are built as CYX and bonded in leap:

      bond mol.5.SG  mol.55.SG
      bond mol.14.SG mol.38.SG
      bond mol.30.SG mol.51.SG

  No free CYS remains, so there is no cysteine protonation question here.

  ALTLOCS. GLU 7 and MET 52 are the only two, both A/B. Conformer A is kept,
  which drops 9 atoms; that is the whole of the 462-to-453 heavy-atom change.

  DEUTERIUM. This is a neutron structure, so 103 of its hydrogen positions are
  deposited as D. They are stripped with the H under rule 3. Do not try to
  keep them: a D is not an H to any of the four engines and mixing them would
  put a mass difference into a comparison that has nothing to do with mass.

  No histidine, so no protonation choice.
  Net charge as built: +6 (4 LYS + 6 ARG + N-terminus against 2 ASP + 2 GLU +
  C-terminus). BPTI really is strongly cationic; do not neutralise it.

## 1BNA -- Drew-Dickerson dodecamer

    chains A and B, 24 residues, 486 heavy atoms   1BNA_AB.pdb
    X-ray 1.9 A, single model, no gaps, no non-standard residues,
    80 waters dropped

    chain A  residues  1-12   CGCGAATTCGCG    DA 2  DC 4  DG 4  DT 2
    chain B  residues 13-24   CGCGAATTCGCG    DA 2  DC 4  DG 4  DT 2

  BOTH CHAINS ARE KEPT. This is the one structure where "chain A unless the
  structure needs otherwise" does not apply. 1BNA is a self-complementary
  duplex and half of it is not the Drew-Dickerson dodecamer.

  TERMINI. DC 1 and DC 13 carry no P/OP1/OP2 in the deposit -- they are free
  5'-hydroxyls already, which is what the deposit means by a 5' end. Nothing
  is removed and nothing is added; the terminal residues build as DC5 and the
  3' ends as DG3. This matches the 5'-OH convention rule 6 imposes on 2KOC,
  so both nucleic systems have the same kind of chain end.

  Net charge as built: -22 (22 internal phosphates, 11 per strand).

## 2KOC -- UUCG tetraloop hairpin, the RNA choice

    chain A, model 1, 14 residues, 294 heavy atoms  2KOC_A_model1.pdb
    solution NMR, 20 models -- MODEL 1 IS TAKEN, as required
    no gaps, no non-standard residues, no missing heavy atoms, no altlocs,
    no HETATM at all

    GGCACUUCGGUGCC       A 1   C 5   G 5   U 3

  WHY 2KOC AND NOT 1F7Y. The brief allowed either. 2KOC is a clean single
  chain of 14 nucleotides containing all four RNA bases, with nothing else in
  the file -- no protein, no ions, no waters, no alternate locations, no
  modified residues. 1F7Y is a ribosomal protein-RNA complex: its RNA is chain
  B of 57 nucleotides, and the file also carries an 89-residue protein chain
  with 32 MSE (selenomethionine) atoms, 9 MG, 1 K, 2 NA, 27 waters, 67 pairs
  of alternate locations and 10 residues absent from the model. Every one of
  those is a decision that would have to be made and defended, and none of
  them is about the tetraloop. 2KOC is the tetraloop and nothing else.

  THE 5'-PHOSPHATE IS REMOVED, AND THIS IS NOT COSMETIC. As deposited, G1
  carries OP3, P, OP1 and OP2 -- a bare 5'-terminal phosphate. Those four
  atoms are dropped, leaving O5' to take a hydroxyl proton, so the built chain
  is 5'-OH / 3'-OH.

  The reason is recorded in evidence/charmm36.txt, C-08: a bare 5'-phosphate
  is the one terminus charmm36_2024.xml cannot express. Every patch in that
  file which can end a chain (5TER, 5MET) carries a RemoveExternalBond and so
  removes the phosphate; every patch that keeps the phosphate (5PHO, 5POM)
  carries none and so cannot end a chain. This is exactly why 1efs has no
  CHARMM36 reference. Keeping 2KOC's deposited 5'-phosphate would buy a second
  system that CHARMM36 cannot be run on, and the point of adding these four
  structures is to widen the force-field coverage, not to repeat the gap.

  The committed nucleic controls (`examples/pure_terminal_controls/`) are
  5TER/3TER for the same reason, so 2KOC now matches them and matches 1BNA.

  Net charge as built: -13 (13 internal phosphates; the 5' terminal phosphate
  is the one that was removed).

## Summary of the chemistry that has to be identical everywhere

    hydrogens        stripped from every deposit, rebuilt by the builder
    altlocs          A (or blank) only
    protein termini  charged NH3+ / COO-, never capped
    nucleic termini  5'-OH and 3'-OH, never 5'-phosphate
    HIS              HIE  (1UBQ HIS 68 -- the only histidine in the set)
    CYS              all six of 5PTI are CYX in three disulfides; no free CYS
    waters/ions      none kept
    vacuum, dielectric 1, no cutoff, no temperature, no move, single point
