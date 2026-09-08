"""OpenMM single point on the tax9 topology built from the committed structure."""
import openmm, openmm.app as app, openmm.unit as u
prm = app.AmberPrmtopFile("../topology/tax9_ff19sb.prmtop")
crd = app.AmberInpcrdFile("../topology/tax9_ff19sb.rst7")
sys_ = prm.createSystem(nonbondedMethod=app.NoCutoff, constraints=None,
                        removeCMMotion=False)
for i, f in enumerate(sys_.getForces()):
    f.setForceGroup(i)
ctx = openmm.Context(sys_, openmm.VerletIntegrator(0.001*u.femtoseconds),
                     openmm.Platform.getPlatformByName("Reference"))
ctx.setPositions(crd.positions)
print("OpenMM", openmm.version.version, "Reference, NoCutoff, no constraints")
print(f"{'term':38s} kcal/mol")
tot = 0.0
for i, f in enumerate(sys_.getForces()):
    e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
    e = e.value_in_unit(u.kilocalories_per_mole)
    tot += e
    print(f"{f.__class__.__name__:38s} {e:16.10f}")
print(f"{'TOTAL':38s} {tot:16.10f}")
