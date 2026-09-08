#!/usr/bin/env python3
"""Single point with OpenMM on an AMBER prmtop, LJ and Coulomb separated.

Same protocol and same printed layout as the committed
evidence/1efs/*/openmm/openmm_sp.out: NoCutoff, no constraints, Reference
platform, one force group per force, then two extra evaluations -- charges
zeroed, and epsilons zeroed -- to split NonbondedForce.
"""
import sys
from pathlib import Path
import openmm
from openmm import Context, Platform, VerletIntegrator, unit
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, NoCutoff

prmtop_path, rst7_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

prmtop = AmberPrmtopFile(prmtop_path)
inpcrd = AmberInpcrdFile(rst7_path)

def build():
    return prmtop.createSystem(constraints=None, rigidWater=False,
                               removeCMMotion=False, nonbondedMethod=NoCutoff)

def energies(system):
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    ctx = Context(system, VerletIntegrator(0.001*unit.picoseconds),
                  Platform.getPlatformByName("Reference"))
    ctx.setPositions(inpcrd.positions)
    out = {}
    for i, f in enumerate(system.getForces()):
        e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
        name = f.__class__.__name__
        out[name] = out.get(name, 0.0) + e.value_in_unit(unit.kilocalories_per_mole)
    return out

base = build()
vals = energies(base)

# charges zeroed -> Lennard-Jones only
s_lj = build()
for f in s_lj.getForces():
    if isinstance(f, openmm.NonbondedForce):
        for i in range(f.getNumParticles()):
            q, sig, eps = f.getParticleParameters(i)
            f.setParticleParameters(i, 0.0*unit.elementary_charge, sig, eps)
        for i in range(f.getNumExceptions()):
            p1, p2, qq, sig, eps = f.getExceptionParameters(i)
            f.setExceptionParameters(i, p1, p2, 0.0, sig, eps)
lj = energies(s_lj)["NonbondedForce"]

# epsilons zeroed -> Coulomb only
s_q = build()
for f in s_q.getForces():
    if isinstance(f, openmm.NonbondedForce):
        for i in range(f.getNumParticles()):
            q, sig, eps = f.getParticleParameters(i)
            f.setParticleParameters(i, q, sig, 0.0*unit.kilojoule_per_mole)
        for i in range(f.getNumExceptions()):
            p1, p2, qq, sig, eps = f.getExceptionParameters(i)
            f.setExceptionParameters(i, p1, p2, qq, sig, 0.0)
coul = energies(s_q)["NonbondedForce"]

L = []
L.append("OpenMM single-point energy")
L.append("=" * 74)
L.append(f"OpenMM version         {openmm.__version__}")
L.append("Platform               Reference (double precision)")
L.append("nonbondedMethod        NoCutoff")
L.append("constraints            None")
L.append("rigidWater             False")
L.append("removeCMMotion         False")
L.append("implicitSolvent        None")
L.append(f"prmtop                 {Path(prmtop_path).name}")
L.append(f"rst7                   {Path(rst7_path).name}")
L.append(f"atoms                  {base.getNumParticles()}")
L.append("")
L.append("Energies in kcal/mol")
L.append("-" * 74)
for k in sorted(vals):
    L.append(f"{k:<30} {vals[k]:.12f}")
L.append("")
L.append("NonbondedForce decomposition")
L.append("-" * 74)
L.append("OpenMM reports one combined non-bonded energy.  The system was")
L.append("evaluated twice more, once with every charge set to zero and once")
L.append("with every Lennard-Jones epsilon set to zero.")
L.append("")
L.append(f"Lennard-Jones (charges zeroed) {lj:.12f}")
L.append(f"Coulomb (epsilons zeroed)      {coul:.12f}")
L.append(f"sum of the two                 {lj + coul:.12f}")
L.append(f"NonbondedForce as built        {vals['NonbondedForce']:.12f}")
L.append(f"difference                     {abs(lj + coul - vals['NonbondedForce']):.3e}")
L.append("")
L.append("The last two lines must agree; the difference is printed so that a")
L.append("reader can confirm it rather than take it on trust.")
Path(out_path).write_text("\n".join(L) + "\n")
print("\n".join(L))
