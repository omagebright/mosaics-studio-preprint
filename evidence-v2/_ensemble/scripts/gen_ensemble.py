#!/usr/bin/env python3
"""Generate a decorrelated conformational ensemble for a MOSAICS cross-engine test.

The comparison protocol used everywhere else in evidence-v2 is a single point at
one conformation.  A parameter that is wrong but inactive at that one geometry is
invisible to it.  This script produces the conformations that make the same
comparison an ensemble measurement (shirts2017engines).

Method
------
Short Langevin MD in OpenMM from the same AMBER topology the reference engines
read.  Two deliberate departures from the single-point protocol, both stated
here because they change what conformations come out:

  * the generating MD uses the OBC2 implicit solvent, not vacuum.  A vacuum
    trajectory of a -22e DNA duplex or of a charged protein is dominated by
    unscreened Coulomb and produces collapsed or exploded geometries rather
    than the range of conformations a force field is actually used over.  The
    ENERGY COMPARISON that consumes these frames is still strict vacuum,
    dielectric 1, no cutoff -- the solvent term exists only to steer the
    trajectory.
  * no bond constraints and a 1 fs step, so the bond and angle terms are
    sampled rather than frozen.  A constrained trajectory would leave the bond
    term at one value and defeat half the test.

Frames are written only every --spacing-ps picoseconds so that successive
frames are not the same conformation twice.  The spread is measured, not
asserted: pairwise heavy-atom RMSD over the saved frames is reported.

Output
------
  positions.npy      (nframes, natoms, 3) float64, Angstrom, as OpenMM produced
  ensemble.json      every setting used, plus the RMSD statistics
  rmsd_pairwise.txt  the full pairwise heavy-atom RMSD matrix summary
"""

from __future__ import annotations

import argparse
import json
import platform as pyplatform
import time
from pathlib import Path

import numpy as np
import openmm
from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, NoCutoff, OBC2, Simulation


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Heavy-atom RMSD after optimal superposition."""
    p = P - P.mean(axis=0)
    q = Q - Q.mean(axis=0)
    cov = p.T @ q
    v, s, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    dm = np.diag([1.0, 1.0, d])
    rot = v @ dm @ wt
    diff = p @ rot - q
    return float(np.sqrt((diff * diff).sum() / len(p)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True, type=Path)
    ap.add_argument("--rst7", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--spacing-ps", type=float, default=50.0)
    ap.add_argument("--equil-ps", type=float, default=200.0)
    ap.add_argument("--temperature-K", type=float, default=300.0)
    ap.add_argument("--friction-per-ps", type=float, default=5.0)
    ap.add_argument("--timestep-fs", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=20260825)
    ap.add_argument("--platform", default="OpenCL")
    ap.add_argument("--precision", default="mixed")
    ap.add_argument("--minimize-steps", type=int, default=0)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    prm = AmberPrmtopFile(str(args.prmtop))
    inp = AmberInpcrdFile(str(args.rst7))
    system = prm.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=True,
        implicitSolvent=OBC2,
    )
    integ = LangevinMiddleIntegrator(
        args.temperature_K * unit.kelvin,
        args.friction_per_ps / unit.picosecond,
        args.timestep_fs * unit.femtosecond,
    )
    integ.setRandomNumberSeed(args.seed)
    plat = Platform.getPlatformByName(args.platform)
    sim = Simulation(prm.topology, system, integ, plat)
    actual_props = {n: sim.context.getPlatform().getPropertyValue(sim.context, n)
                    for n in sim.context.getPlatform().getPropertyNames()}
    sim.context.setPositions(inp.positions)
    if args.minimize_steps:
        sim.minimizeEnergy(maxIterations=args.minimize_steps)
    sim.context.setVelocitiesToTemperature(args.temperature_K * unit.kelvin, args.seed)

    steps_per_ps = int(round(1000.0 / args.timestep_fs))
    equil_steps = int(round(args.equil_ps * steps_per_ps))
    spacing_steps = int(round(args.spacing_ps * steps_per_ps))

    t0 = time.time()
    sim.step(equil_steps)
    frames = []
    for _ in range(args.frames):
        sim.step(spacing_steps)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        frames.append(np.asarray(pos, dtype=np.float64))
    wall = time.time() - t0

    X = np.stack(frames)
    np.save(args.out_dir / "positions.npy", X)

    masses = np.array(
        [system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles())]
    )
    heavy = np.where(masses > 2.5)[0]
    ref = np.asarray(inp.positions.value_in_unit(unit.angstrom), dtype=np.float64)

    n = len(X)
    pair = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = kabsch_rmsd(X[i][heavy], X[j][heavy])
            pair[i, j] = pair[j, i] = r
    to_ref = np.array([kabsch_rmsd(X[i][heavy], ref[heavy]) for i in range(n)])
    iu = np.triu_indices(n, 1)
    pv = pair[iu]

    with (args.out_dir / "rmsd_pairwise.txt").open("w") as fh:
        fh.write("Heavy-atom RMSD after optimal superposition, Angstrom\n")
        fh.write("nframes %d, nheavy %d\n\n" % (n, len(heavy)))
        fh.write("pairwise over all %d frame pairs\n" % len(pv))
        fh.write("  min  %.3f\n  mean %.3f\n  max  %.3f\n" % (pv.min(), pv.mean(), pv.max()))
        fh.write("  adjacent-frame (i, i+1) min %.3f mean %.3f\n"
                 % (np.diag(pair, 1).min(), np.diag(pair, 1).mean()))
        fh.write("\nRMSD to the reference (single-point) conformation\n")
        fh.write("  min  %.3f\n  mean %.3f\n  max  %.3f\n"
                 % (to_ref.min(), to_ref.mean(), to_ref.max()))
        fh.write("\nper frame RMSD to reference\n")
        for i, r in enumerate(to_ref):
            fh.write("  %3d  %.3f\n" % (i, r))

    meta = {
        "prmtop": str(args.prmtop),
        "rst7": str(args.rst7),
        "openmm_version": openmm.version.version,
        "python": pyplatform.python_version(),
        "platform": args.platform,
        "platform_properties": actual_props,
        "generating_potential": "AMBER topology + OBC2 implicit solvent, NoCutoff",
        "integrator": "LangevinMiddle",
        "temperature_K": args.temperature_K,
        "friction_per_ps": args.friction_per_ps,
        "timestep_fs": args.timestep_fs,
        "constraints": "none",
        "seed": args.seed,
        "equil_ps": args.equil_ps,
        "spacing_ps": args.spacing_ps,
        "frames": args.frames,
        "total_md_ps": args.equil_ps + args.spacing_ps * args.frames,
        "wall_seconds": round(wall, 1),
        "natoms": int(X.shape[1]),
        "nheavy": int(len(heavy)),
        "rmsd_pairwise_min": round(float(pv.min()), 4),
        "rmsd_pairwise_mean": round(float(pv.mean()), 4),
        "rmsd_pairwise_max": round(float(pv.max()), 4),
        "rmsd_adjacent_min": round(float(np.diag(pair, 1).min()), 4),
        "rmsd_to_reference_min": round(float(to_ref.min()), 4),
        "rmsd_to_reference_mean": round(float(to_ref.mean()), 4),
        "rmsd_to_reference_max": round(float(to_ref.max()), 4),
    }
    (args.out_dir / "ensemble.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(json.dumps(meta, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
