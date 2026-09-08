#!/usr/bin/env python3
"""Single-point AMBER prmtop energy in OpenMM, in the committed openmm_sp.out format.

Reproduces the layout of evidence/1efs/*/openmm/openmm_sp.out: per-force
energies, then the NonbondedForce split into Lennard-Jones and Coulomb by
re-evaluating the same System twice, once with every charge zeroed and once
with every Lennard-Jones epsilon zeroed.
"""
import argparse, copy
from pathlib import Path
import openmm
from openmm import Context, Platform, VerletIntegrator, unit
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, NoCutoff

KCAL = unit.kilocalories_per_mole

def group_energies(system, positions, platform):
    s = copy.deepcopy(system)
    names = []
    for i, f in enumerate(s.getForces()):
        f.setForceGroup(i)
        names.append(f.__class__.__name__)
    ctx = Context(s, VerletIntegrator(0.001*unit.picoseconds), platform)
    ctx.setPositions(positions)
    out = {}
    for i, n in enumerate(names):
        e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(KCAL)
        out[n] = out.get(n, 0.0) + e
    del ctx
    return out

def nb_variant(system, positions, platform, zero):
    s = copy.deepcopy(system)
    for f in s.getForces():
        if isinstance(f, openmm.NonbondedForce):
            for i in range(f.getNumParticles()):
                q, sig, eps = f.getParticleParameters(i)
                if zero == "charge":
                    f.setParticleParameters(i, 0.0, sig, eps)
                else:
                    f.setParticleParameters(i, q, sig, 0.0)
            for i in range(f.getNumExceptions()):
                p1, p2, qq, sig, eps = f.getExceptionParameters(i)
                if zero == "charge":
                    f.setExceptionParameters(i, p1, p2, 0.0, sig, eps)
                else:
                    f.setExceptionParameters(i, p1, p2, qq, sig, 0.0)
            f.setForceGroup(31)
        else:
            f.setForceGroup(0)
    ctx = Context(s, VerletIntegrator(0.001*unit.picoseconds), platform)
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True, groups={31}).getPotentialEnergy().value_in_unit(KCAL)
    del ctx
    return e

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--rst7", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    prmtop = AmberPrmtopFile(a.prmtop)
    inpcrd = AmberInpcrdFile(a.rst7)
    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=None,
                                 rigidWater=False, removeCMMotion=False)
    platform = Platform.getPlatformByName("Reference")
    pos = inpcrd.positions

    per_force = group_energies(system, pos, platform)
    lj = nb_variant(system, pos, platform, "charge")
    coul = nb_variant(system, pos, platform, "epsilon")
    nb_as_built = per_force.get("NonbondedForce", 0.0)

    L = []
    L.append("OpenMM single-point energy")
    L.append("=" * 74)
    L.append(f"{'OpenMM version':<22} {openmm.__version__}")
    L.append(f"{'Platform':<22} Reference (double precision)")
    L.append(f"{'nonbondedMethod':<22} NoCutoff")
    L.append(f"{'constraints':<22} None")
    L.append(f"{'rigidWater':<22} False")
    L.append(f"{'removeCMMotion':<22} False")
    L.append(f"{'implicitSolvent':<22} None")
    L.append(f"{'prmtop':<22} {Path(a.prmtop).name}")
    L.append(f"{'rst7':<22} {Path(a.rst7).name}")
    L.append(f"{'atoms':<22} {system.getNumParticles()}")
    L.append("")
    L.append("Energies in kcal/mol")
    L.append("-" * 74)
    for k in sorted(per_force):
        L.append(f"{k:<30} {per_force[k]:.12f}")
    L.append("")
    L.append("NonbondedForce decomposition")
    L.append("-" * 74)
    L.append("OpenMM reports one combined non-bonded energy.  The system was")
    L.append("evaluated twice more, once with every charge set to zero and once")
    L.append("with every Lennard-Jones epsilon set to zero.")
    L.append("")
    L.append(f"{'Lennard-Jones (charges zeroed)':<30} {lj:.12f}")
    L.append(f"{'Coulomb (epsilons zeroed)':<30} {coul:.12f}")
    L.append(f"{'sum of the two':<30} {lj+coul:.12f}")
    L.append(f"{'NonbondedForce as built':<30} {nb_as_built:.12f}")
    L.append(f"{'difference':<30} {abs((lj+coul)-nb_as_built):.3e}")
    L.append("")
    L.append("The last two lines must agree; the difference is printed so that a")
    L.append("reader can confirm it rather than take it on trust.")
    Path(a.out).write_text("\n".join(L) + "\n")
    print("\n".join(L))

main()
