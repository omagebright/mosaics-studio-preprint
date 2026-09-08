"""Build a CHARMM36 nucleic System in OpenMM, and say which terminus works.

WHY THE BONDS HAVE TO BE BUILT HERE

OpenMM's PDB reader creates intra-residue bonds from its own table of standard
residue names. It does not know RA5, RC, RG or RU3, so it reads those residues
as a set of atoms with no bonds at all -- every atom in the first residue of
examples/pure_terminal_controls/rna_acgu.pdb comes back with zero bonds.

Template matching then cannot succeed for any residue, and the message it
prints is about the nearest template by element count:

    No template found for residue 0 (ADE). The set of atoms is similar to
    OMA, but is missing 1 H atom, 1 C atom, 2 O atoms, and 1 P atom.

That message reads like a statement about the chemistry and is not one. It is
what an unbonded residue looks like. Anything concluded about which terminus
CHARMM can express, before the bonds exist, is concluded from noise.

So this script builds the connectivity itself, from the force field's own
residue templates, and only then asks OpenMM to build a System.

WHAT IT COMPARES

Two candidate chain ends, on the same strand:

    5TER   removes P, O1P and O2P and adds H5T on O5'   -> 5-prime-OH
    5PHO   keeps the phosphate and adds O5T and H5T     -> 5-prime-phosphate

5TER and 5MET each carry a RemoveExternalBond element and 5PHO does not, so on
a reading of the file 5PHO leaves P's external bond unsatisfied and cannot end
a chain. This script tests that rather than asserting it.

Usage:
    python3 charmm_reference_energy.py <structure.pdb> <5TER|5PHO>
"""

from __future__ import annotations

import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from openmm import app, unit


DATA = Path(app.__file__).parent / "data"
XML = "charmm36_2024.xml"

# PDB residue names -> CHARMM template names
CHARMM_NAME = {
    "A": "ADE", "A5": "ADE", "A3": "ADE",
    "DA": "ADE", "DA5": "ADE", "DA3": "ADE",
    "RA": "ADE", "RA5": "ADE", "RA3": "ADE", "ADE": "ADE",
    "G": "GUA", "G5": "GUA", "G3": "GUA",
    "DG": "GUA", "DG5": "GUA", "DG3": "GUA",
    "RG": "GUA", "RG5": "GUA", "RG3": "GUA", "GUA": "GUA",
    "C": "CYT", "C5": "CYT", "C3": "CYT",
    "DC": "CYT", "DC5": "CYT", "DC3": "CYT",
    "RC": "CYT", "RC5": "CYT", "RC3": "CYT", "CYT": "CYT",
    "U": "URA", "U5": "URA", "U3": "URA",
    "RU": "URA", "RU5": "URA", "RU3": "URA", "URA": "URA",
    "T": "THY", "T5": "THY", "T3": "THY",
    "DT": "THY", "DT5": "THY", "DT3": "THY", "THY": "THY",
}


def charmm_forcefield():
    """ForceField for charmm36_2024.xml, with thymine 3-prime DNA unambiguous.

    That residue matches both THY+DEOX+3TER and URA+DEOX+3TER+5-methyl.
    Prefer the template whose name starts with the CHARMM residue we assigned.
    """
    import types
    from openmm.app.forcefield import _createResidueSignature
    from openmm.app import forcefield as _ffmod
    ff = app.ForceField(XML)
    _orig = app.ForceField._getResidueTemplateMatches

    def _prefer_native(self, res, bondedToAtom, templateSignatures=None,
                       ignoreExternalBonds=False, ignoreExtraParticles=False):
        if templateSignatures is None:
            templateSignatures = self._templateSignatures
        signature = _createResidueSignature([atom.element for atom in res.atoms()])
        if signature in templateSignatures:
            all_matches = []
            for templ in templateSignatures[signature]:
                match = _ffmod.compiled.matchResidueToTemplate(
                    res, templ, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                if match is not None:
                    all_matches.append((templ, match))
            native = [m for m in all_matches
                      if m[0].name.split("-")[0] == res.name]
            if len(native) == 1:
                return [native[0][0], native[0][1]]
        return _orig(self, res, bondedToAtom, templateSignatures,
                     ignoreExternalBonds, ignoreExtraParticles)

    ff._getResidueTemplateMatches = types.MethodType(_prefer_native, ff)
    return ff


def load_definitions():
    root = ET.parse(DATA / XML).getroot()
    templates = {}
    for res in root.find("Residues").findall("Residue"):
        templates[res.get("name")] = {
            "atoms": [a.get("name") for a in res.findall("Atom")],
            "bonds": [(b.get("atomName1"), b.get("atomName2"))
                      for b in res.findall("Bond")],
            "external": [e.get("atomName") for e in res.findall("ExternalBond")],
        }
    patches = {}
    for p in root.find("Patches").findall("Patch"):
        patches[p.get("name")] = {
            "add_atoms": [a.get("name") for a in p.findall("AddAtom")],
            "remove_atoms": [a.get("name") for a in p.findall("RemoveAtom")],
            "add_bonds": [(b.get("atomName1"), b.get("atomName2"))
                          for b in p.findall("AddBond")],
            "remove_bonds": [(b.get("atomName1"), b.get("atomName2"))
                             for b in p.findall("RemoveBond")],
            "removes_external": len(p.findall("RemoveExternalBond")),
        }
    return templates, patches


def bond_set(template, applied=()):
    """The template's internal bonds, with each patch's additions and removals.

    Patches compose: a DNA chain end carries DEOX as well as a terminal patch,
    and CHARMM ships DEO5TER for exactly that pair. They are applied in the
    order given.
    """
    bonds = {frozenset(b) for b in template["bonds"]}
    for patch in applied:
        if patch is None:
            continue
        for b in patch["remove_bonds"]:
            bonds.discard(frozenset(b))
        for a in patch["remove_atoms"]:
            bonds = {b for b in bonds if a not in b}
        for b in patch["add_bonds"]:
            bonds.add(frozenset(b))
    return bonds


def build(pdb_path: Path, five_prime: str):
    templates, patches = load_definitions()
    pdb = app.PDBFile(str(pdb_path))
    top = pdb.topology

    residues = list(top.residues())
    for r in residues:
        if r.name not in CHARMM_NAME:
            raise SystemExit(f"no CHARMM name for residue {r.name!r}")
        r.name = CHARMM_NAME[r.name]

    # The committed controls spell the terminal hydroxyl hydrogens HO5' and
    # HO3'; CHARMM's 5TER and 3TER patches call them H5T and H3T. Renaming
    # them here means the control itself can be used, with no second copy of
    # a structure to keep in step with the first.
    renamed = 0
    # OpenMM's PDB reader canonicalises DNA to PDB names (OP1, C7) even
    # when the file already used CHARMM's O1P and C5M. Put them back so
    # template matching sees the same names as charmm36_2024.xml.
    for name_from, name_to in (
            ("HO5'", "H5T"), ("HO3'", "H3T"),
            ("OP1", "O1P"), ("OP2", "O2P"),
            ("C7", "C5M"), ("H71", "H51"), ("H72", "H52"), ("H73", "H53"),
    ):
        for r in residues:
            for a in r.atoms():
                if a.name == name_from:
                    a.name = name_to
                    renamed += 1
    if renamed:
        print(f"  renamed {renamed} atom(s) to CHARMM spelling")

    five_patch = f"{five_prime}_0"
    three_patch = "3TER_1"
    for needed in (five_patch, three_patch):
        if needed not in patches:
            raise SystemExit(f"{needed} is not a patch in {XML}")

    # OpenMM's PDB reader already bonds residues in its standard table.
    # DC and DG are in that table; RA5 and RU3 are not. The RNA path
    # therefore starts from an unbonded topology and the DNA path does
    # not, unless those bonds are dropped. CHARMM templates are the
    # source of truth either way -- keep one copy of each pair.
    top._bonds.clear()

    first_of_chain = {list(c.residues())[0] for c in top.chains()}
    last_of_chain = {list(c.residues())[-1] for c in top.chains()}
    added = 0
    deoxy = 0
    for res in residues:
        applied = []
        # Deoxyribose is a patch in CHARMM, not a separate template. A residue
        # with no O2' is DNA and carries DEOX on top of anything else.
        if not any(a.name == "O2'" for a in res.atoms()):
            applied.append(patches["DEOX_0"])
            deoxy += 1
        if res in first_of_chain:
            applied.append(patches[five_patch])
        elif res in last_of_chain:
            applied.append(patches[three_patch])
        wanted = bond_set(templates[res.name], applied)
        by_name = {a.name: a for a in res.atoms()}
        for pair in wanted:
            n1, n2 = tuple(pair)
            if n1 in by_name and n2 in by_name:
                top.addBond(by_name[n1], by_name[n2])
                added += 1

    # The backbone link: O3' of one residue to P of the next, WITHIN A CHAIN.
    #
    # Walking the flat residue list would bond the last residue of one strand
    # to the first of the next, which on a duplex silently makes two chains
    # into one and changes every energy that depends on connectivity.
    links = 0
    for chain in top.chains():
        chain_residues = list(chain.residues())
        for a, b in zip(chain_residues, chain_residues[1:]):
            left = {x.name: x for x in a.atoms()}
            right = {x.name: x for x in b.atoms()}
            if "O3'" in left and "P" in right:
                top.addBond(left["O3'"], right["P"])
                links += 1

    lonely = [a.name for a in residues[0].atoms()
              if not any(a in (x, y) for x, y in top.bonds())]
    print(f"  bonds built {added} intra-residue, {links} backbone links, "
          f"{deoxy} deoxyribose")
    print(f"  unbonded atoms in residue 1: {lonely if lonely else 'none'}")
    print(f"  {five_patch} removes an external bond: "
          f"{'yes' if patches[five_patch]['removes_external'] else 'NO'}")

    ff = charmm_forcefield()
    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff,
                             constraints=None, rigidWater=False)

    # Split the Urey-Bradley term out of the bond term.
    #
    # OpenMM has no Urey-Bradley force. CHARMM's 1-3 harmonic goes into
    # HarmonicBondForce alongside the real bonds, so the reported "bond"
    # energy is two different terms added together. MOSAICS carries the 1-3
    # part in the bend record's own bond slot, so a term-by-term comparison
    # needs them apart.
    #
    # A bond of this force whose two atoms are not bonded in the topology is
    # a 1-3 pair. On the RNA control that is 42 of 178.
    from openmm import HarmonicBondForce
    topology_bonds = {frozenset((b[0].index, b[1].index)) for b in top.bonds()}
    for index in range(system.getNumForces()):
        force = system.getForce(index)
        if not isinstance(force, HarmonicBondForce):
            continue
        urey = HarmonicBondForce()
        keep = []
        for j in range(force.getNumBonds()):
            a1, a2, length, k = force.getBondParameters(j)
            if frozenset((a1, a2)) in topology_bonds:
                keep.append((a1, a2, length, k))
            else:
                urey.addBond(a1, a2, length, k)
        if urey.getNumBonds() == 0:
            break
        while force.getNumBonds() > len(keep):
            force.setBondParameters(force.getNumBonds() - 1, 0, 0, 0.0, 0.0)
            break
        # rebuild the bond force with only the real bonds
        rebuilt = HarmonicBondForce()
        for a1, a2, length, k in keep:
            rebuilt.addBond(a1, a2, length, k)
        system.removeForce(index)
        bond_index = system.addForce(rebuilt)
        urey_index = system.addForce(urey)
        split = {bond_index: "bond", urey_index: "Urey-Bradley (1-3)"}
        break
    else:
        split = {}

    return system, top, pdb, split


def main() -> int:
    pdb_path = Path(sys.argv[1])
    five_prime = sys.argv[2] if len(sys.argv) > 2 else "5TER"
    print(f"{pdb_path.name}, 5-prime patch {five_prime}")
    try:
        system, top, pdb, split = build(pdb_path, five_prime)
    except SystemExit:
        raise
    except Exception as exc:
        print(f"  BUILD FAILED: {str(exc)[:300]}")
        return 1

    print(f"  BUILD OK   {system.getNumParticles()} particles, "
          f"{system.getNumForces()} forces")

    from openmm import LangevinIntegrator, Platform, Context

    # One force group per force, so the total can be reported term by term.
    # A single total is not comparable against MOSAICS, which reports bond,
    # bend, torsion, onefour and inter separately.
    for index in range(system.getNumForces()):
        system.getForce(index).setForceGroup(index)

    integrator = LangevinIntegrator(300, 1, 0.001)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions(pdb.positions)

    # What each OpenMM force IS under this force field. The class name is not
    # the term: CustomBondForce here is the 1-4 Lennard-Jones, not
    # Urey-Bradley, and its energy expression says so --
    # 4*epsilon*((sigma/r)^12-(sigma/r)^6). None of its 320 pairs is a bond in
    # the topology. Urey-Bradley has no force of its own in OpenMM; CHARMM's
    # 1-3 harmonic goes into HarmonicBondForce with the real bonds, which is
    # why build() separates them.
    MEANING = {
        "HarmonicBondForce": "bond",
        "HarmonicAngleForce": "angle",
        "PeriodicTorsionForce": "torsion",
        "CustomTorsionForce": "improper (harmonic)",
        "CustomBondForce": "1-4 Lennard-Jones",
        "NonbondedForce": "electrostatics",
        "CustomNonbondedForce": "Lennard-Jones",
        "CMAPTorsionForce": "correction map",
        "CMMotionRemover": "-",
    }
    kcal = unit.kilocalorie_per_mole
    total = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  {'term':<22} {'kcal/mol':>20}   {'openmm force':<22}")
    for index in range(system.getNumForces()):
        cls = system.getForce(index).__class__.__name__
        term = split.get(index) or MEANING.get(cls, cls)
        e = context.getState(getEnergy=True,
                             groups={index}).getPotentialEnergy()
        print(f"  {term:<22} {e.value_in_unit(kcal):>20.10f}   {cls:<22}")
    print(f"  {'TOTAL':<22} {total.value_in_unit(kcal):>20.10f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
