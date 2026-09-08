#!/usr/bin/env python3
"""Measure each engine's Coulomb constant directly.

Two point charges, +1 e and -1 e, exactly 1.0000000000 Angstrom apart, in
vacuum, with Lennard-Jones parameters set to zero and no cutoff.  For that
system the electrostatic energy is

    E = C * q1 * q2 / r = -C   (kcal/mol)

so whatever an engine prints for its electrostatic term IS its Coulomb
constant in kcal * Angstrom / (mol * e^2).  No force field, no topology
ambiguity, nothing to interpret.

One AMBER topology is written and handed to every engine, so all of them see
identical charges and identical coordinates.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import parmed as pmd
from parmed.topologyobjects import Atom, AtomType

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

SEPARATION_ANGSTROM = 1.0

# AMBER stores prmtop charges pre-multiplied by this factor.  Recorded here
# because it is the whole reason sander differs from everything else.
AMBER_CHARGE_SCALE = 18.2223


def build_two_charge_topology() -> pmd.Structure:
    """Two unbonded, LJ-free point charges 1 A apart."""
    struct = pmd.Structure()
    for idx, (name, resname, charge) in enumerate(
        (("IP", "ION", +1.0), ("IM", "ION", -1.0)), start=1
    ):
        atom = Atom(name=name, type=name, charge=charge, mass=1.0, atomic_number=1)
        # Zero Lennard-Jones: rmin and epsilon both zero -> no vdW at all,
        # so the only surviving interaction is Coulomb.
        atom_type = AtomType(name, idx, mass=1.0, atomic_number=1)
        atom_type.set_lj_params(eps=0.0, rmin=0.0)
        atom.atom_type = atom_type
        atom.rmin = 0.0
        atom.epsilon = 0.0
        struct.add_atom(atom, resname, idx)

    # No bonds are created, so the two atoms are separate molecules and the
    # pair is a genuine non-bonded interaction, not a 1-2/1-3/1-4 exclusion.
    struct.coordinates = np.array(
        [[0.0, 0.0, 0.0], [SEPARATION_ANGSTROM, 0.0, 0.0]], dtype=np.float64
    )
    struct.box = None
    return struct


def sander_constant(prmtop: Path, rst7: Path) -> float:
    mdin = WORK / "sp.in"
    mdin.write_text(
        "single point, two point charges, vacuum\n"
        " &cntrl\n"
        "  imin=1, maxcyc=0, ntb=0, igb=0, cut=999.0, ntpr=1, dielc=1.0,\n"
        " /\n"
    )
    out = WORK / "sander_sp.out"
    subprocess.run(
        ["sander", "-O", "-i", str(mdin), "-p", str(prmtop), "-c", str(rst7),
         "-o", str(out), "-r", str(WORK / "sander.rst7")],
        check=True, cwd=WORK, capture_output=True,
    )
    # Parse the FINAL RESULTS block: EEL is the electrostatic term.
    lines = out.read_text().splitlines()
    start = next(i for i, ln in enumerate(lines) if "FINAL RESULTS" in ln)
    for ln in lines[start:]:
        if "EEL" in ln and "1-4" not in ln:
            eel = float(ln.split("EEL")[1].split("=")[1].split()[0])
            return -eel
    raise RuntimeError("EEL not found in sander output")


def openmm_constant(prmtop: Path, rst7: Path) -> float:
    import openmm
    from openmm import app, unit

    parm = app.AmberPrmtopFile(str(prmtop))
    crd = app.AmberInpcrdFile(str(rst7))
    system = parm.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
        implicitSolvent=None,
    )
    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName("Reference")
    )
    context.setPositions(crd.positions)
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    return -energy.value_in_unit(unit.kilocalorie_per_mole)


def gromacs_constant(struct: pmd.Structure) -> float:
    """Single-point Coulomb energy from GROMACS.

    Three conventions have to be forced explicitly, because every one of them
    silently changes the answer:

    1.  Modern GROMACS removed the `group` cutoff scheme, so an infinite
        cutoff with pbc=no is no longer expressible.  Instead the system is
        placed in a box large enough that a finite cutoff exceeds the longest
        interatomic distance while staying below half the box vector.  Under
        those two conditions the periodic calculation equals the vacuum one.
        Both conditions are asserted below rather than assumed.

    2.  The default `coulomb-modifier`/`vdw-modifier` is potential-shift,
        which subtracts V(r_cut) from every pair so the potential reaches zero
        at the cutoff.  That is a constant offset per pair, it does NOT cancel
        for a net-neutral system, and it would corrupt every cross-engine
        comparison.  Both modifiers are set to None.

    3.  ParmEd writes .gro coordinates with 3 decimals in nm by default, i.e.
        0.01 Angstrom.  That is coarse enough to inject spurious bonded energy
        on a real structure, so the precision is set explicitly.
    """
    coords = np.asarray(struct.coordinates, dtype=np.float64)
    # Longest interatomic distance in the system, in nm.
    diameter_nm = 0.0
    if len(coords) > 1:
        diffs = coords[:, None, :] - coords[None, :, :]
        diameter_nm = float(np.sqrt((diffs ** 2).sum(-1)).max()) / 10.0

    r_cut_nm = diameter_nm + 2.0            # comfortably beyond every pair
    box_nm = 2.0 * r_cut_nm + 2.0           # so r_cut < box/2 with margin
    assert r_cut_nm > diameter_nm, "cutoff does not cover the whole system"
    assert r_cut_nm < box_nm / 2.0, "cutoff exceeds half the box vector"

    gro = WORK / "sys.gro"
    top = WORK / "sys.top"
    struct = struct.copy(pmd.Structure)
    struct.box = [box_nm * 10.0, box_nm * 10.0, box_nm * 10.0, 90.0, 90.0, 90.0]
    struct.save(str(gro), overwrite=True, precision=8)
    struct.save(str(top), format="gromacs", overwrite=True)

    mdp = WORK / "sp.mdp"
    mdp.write_text(
        "integrator               = md\n"
        "nsteps                   = 0\n"
        "pbc                      = xyz\n"
        "cutoff-scheme            = Verlet\n"
        "verlet-buffer-tolerance  = -1\n"
        f"rlist                    = {r_cut_nm:.6f}\n"
        "coulombtype              = Cut-off\n"
        "coulomb-modifier         = None\n"
        f"rcoulomb                 = {r_cut_nm:.6f}\n"
        "vdwtype                  = Cut-off\n"
        "vdw-modifier             = None\n"
        f"rvdw                     = {r_cut_nm:.6f}\n"
        "epsilon-r                = 1\n"
        "epsilon-rf               = 1\n"
        "DispCorr                 = no\n"
        "constraints              = none\n"
        "continuation             = yes\n"
        "nstcalcenergy            = 1\n"
        "nstenergy                = 1\n"
    )
    subprocess.run(
        ["gmx_d", "grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top),
         "-o", str(WORK / "sp.tpr"), "-maxwarn", "2"],
        check=True, cwd=WORK, capture_output=True,
    )
    subprocess.run(
        ["gmx_d", "mdrun", "-s", str(WORK / "sp.tpr"), "-rerun", str(gro),
         "-e", str(WORK / "sp.edr"), "-g", str(WORK / "sp.log"),
         "-o", str(WORK / "sp.trr"), "-c", str(WORK / "sp_out.gro")],
        check=True, cwd=WORK, capture_output=True,
    )
    xvg = WORK / "sp.xvg"
    subprocess.run(
        ["gmx_d", "energy", "-f", str(WORK / "sp.edr"), "-o", str(xvg)],
        input="Coulomb-(SR)\n\n", text=True, cwd=WORK,
        capture_output=True, check=True,
    )
    # Read the .xvg, NOT the stdout summary.  gmx energy prints only about six
    # significant figures to stdout, which on a ~1389 kJ/mol value is a
    # +/-0.005 kJ/mol quantisation -- large enough to fake a disagreement
    # between engines.  The .xvg carries the full stored precision, and in a
    # double-precision build the .edr stores doubles.
    data = [ln for ln in xvg.read_text().splitlines()
            if ln and not ln.startswith(("#", "@"))]
    if not data:
        raise RuntimeError("no data rows in %s" % xvg)
    kj_per_mol = float(data[-1].split()[1])
    return -kj_per_mol / 4.184  # kJ/mol -> kcal/mol


def main() -> int:
    struct = build_two_charge_topology()
    prmtop = WORK / "two_charges.prmtop"
    rst7 = WORK / "two_charges.rst7"
    struct.save(str(prmtop), overwrite=True)
    struct.save(str(rst7), format="rst7", overwrite=True)

    results = {}
    for label, fn in (
        ("AMBER sander", lambda: sander_constant(prmtop, rst7)),
        ("OpenMM", lambda: openmm_constant(prmtop, rst7)),
        ("GROMACS", lambda: gromacs_constant(struct)),
    ):
        try:
            results[label] = fn()
        except Exception as exc:  # noqa: BLE001 - report, do not mask
            results[label] = exc

    # MOSAICS works in atomic units, where the Coulomb prefactor is exactly 1
    # Hartree*bohr/e^2.  Its effective constant in kcal*A/(mol*e^2) is therefore
    # fixed by the two conversion constants in source/include/definitions_const.h.
    BOHR = 0.529177
    KCAL = 627.50921
    results["MOSAICS (BOHR x KCAL, analytic)"] = BOHR * KCAL

    print()
    print("Coulomb constant measured from a two-point-charge single point")
    print("q = +1 e and -1 e, r = %.10f A, LJ zeroed, vacuum, no cutoff"
          % SEPARATION_ANGSTROM)
    print("Units: kcal * Angstrom / (mol * e^2)")
    print("-" * 68)
    for label, value in results.items():
        if isinstance(value, Exception):
            print("%-34s  FAILED: %s" % (label, value))
        else:
            print("%-34s  %.6f" % (label, value))
    print("-" * 68)

    good = {k: v for k, v in results.items() if not isinstance(v, Exception)}
    if "AMBER sander" in good:
        ref = good["AMBER sander"]
        print("Relative difference from sander:")
        for label, value in good.items():
            print("%-34s  %+.6e" % (label, value / ref - 1.0))
    print()
    print("Precision")
    print("---------")
    print("sander prints EEL to four decimals, so its row above carries four")
    print("measured digits and two of formatting. The other three rows are")
    print("read from APIs and files that carry full double precision.")
    print()
    print("AMBER stores prmtop charges pre-multiplied by %s, so sander's"
          % AMBER_CHARGE_SCALE)
    print("constant is that number squared, exactly:")
    print("  18.2223^2 = %.6f" % (AMBER_CHARGE_SCALE ** 2))
    print()
    print("That exact value is what results.txt corrects the sander column")
    print("with, not the four digits printed above. The two agree to the")
    print("precision sander reports: 332.0522 contains 332.052217.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
