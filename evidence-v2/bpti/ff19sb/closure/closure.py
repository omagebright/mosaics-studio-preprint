#!/usr/bin/env python3
"""Close the BPTI ff19SB dihedral-term gap between MOSAICS-current and OpenMM.

Every number is recomputed here from the prmtop, the rst7, the MOSAICS RTF and
the MOSAICS improper deck; nothing is taken from a previous analysis.

Sign convention: IUPAC/AMBER/OpenMM (phi = atan2(|b2| b1.(b2xb3), (b1xb2).(b2xb3))).
The sign only matters because the prmtop stores phase = 3.141594 rad, not pi.

Usage:
  closure.py <prmtop> <rst7> <rtf> <deck> <mosaics_sp.out>
"""
from __future__ import annotations
import math, pathlib, re, sys
import numpy as np
import parmed
import openmm, openmm.app as app, openmm.unit as u

BOLTZ = 315777.0      # definitions_const.h:16
KCAL = 627.50921      # definitions_const.h:17
KJ = 4.184


def dihedral(c):
    b1, b2, b3 = c[1] - c[0], c[2] - c[1], c[3] - c[2]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return math.atan2(np.linalg.norm(b2) * np.dot(b1, n2), np.dot(n1, n2))


def read_rtf(path):
    out, cur = {}, None
    for line in pathlib.Path(path).read_text().splitlines():
        s = line.strip()
        if s.startswith("RESI "):
            cur = s.split()[1]
            out[cur] = {"atoms": {}, "impr": []}
        elif cur and s.startswith("ATOM "):
            t = s.split()
            out[cur]["atoms"][t[1]] = t[2]
        elif cur and (s.startswith("IMPROPER") or s.startswith("IMPR ")):
            toks = s.split()[1:]
            for i in range(0, len(toks) - 3, 4):
                out[cur]["impr"].append(tuple(toks[i:i + 4]))
    return out


def read_deck(path):
    txt = pathlib.Path(path).read_text()
    return [dict(re.findall(r"\\(\w+)\{([^}]*)\}", b))
            for b in re.findall(r"~torsion_parm\[(.*?)\]", txt, re.S)]


def wild_eq(a, b):
    return a.upper() == "X" or b.upper() == "X" or a.upper() == b.upper()


def deck_lookup(ents, types):
    """tors.cpp:186-210 (last exact match wins, else last wildcard match)."""
    exact = wild = None
    for e in ents:
        if e.get("label", "").lower() != "improper":
            continue
        lbl = [e["atom1"], e["atom2"], e["atom3"], e["atom4"]]
        if all(types[i].upper() == lbl[i].upper() for i in range(4)) or \
           all(types[i].upper() == lbl[3 - i].upper() for i in range(4)):
            exact = e
        if all(wild_eq(types[i], lbl[i]) for i in range(4)) or \
           all(wild_eq(types[i], lbl[3 - i]) for i in range(4)):
            wild = e
    return exact if exact is not None else wild


def mosaics_energy(e, phi):
    """tors.cpp:756-779, converted to kcal/mol."""
    g = lambda k: float(e.get(k, 0))
    c, s = math.cos(phi), math.sin(phi)
    vc = sum(g(f"p{i}") * c ** i for i in range(7))
    sisum = sum(g(f"s{i}") * s ** (i - 1) for i in range(2, 7))
    vs = sisum * c + g("s1") * s
    vsp = sum(g(f"sp{i}") * s ** i for i in range(2, 7))
    return (vc + vs + vsp) / BOLTZ * KCAL


def main():
    prm, rst, rtf_p, deck_p, out_p = sys.argv[1:6]
    p = parmed.load_file(prm, xyz=rst)
    rtf = read_rtf(rtf_p)
    deck = read_deck(deck_p)
    res = list(p.residues)
    xyz = {a.idx: np.array([a.xx, a.xy, a.xz]) for a in p.atoms}
    txt = pathlib.Path(out_p).read_text()
    mos_imp = float(re.search(r"Tors-Improper energy\s+([-\d.eE+]+)", txt).group(1))
    mos_prop = float(re.search(r"Tors-Tors energy\s+([-\d.eE+]+)", txt).group(1))

    # ---- OpenMM reference, whole PeriodicTorsionForce ------------------------
    top = app.AmberPrmtopFile(prm)
    crd = app.AmberInpcrdFile(rst)
    sysm = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
    for i, f in enumerate(sysm.getForces()):
        f.setForceGroup(i)
    ctx = openmm.Context(sysm, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(crd.positions)
    ig = [i for i, f in enumerate(sysm.getForces()) if isinstance(f, openmm.PeriodicTorsionForce)][0]
    omm_dihed = ctx.getState(getEnergy=True, groups={ig}).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole) / KJ

    # ---- analytic recompute of every prmtop term, prmtop phase ---------------
    def E_prm(d, phi, phase_rad):
        return d.type.phi_k * (1 + math.cos(d.type.per * phi - phase_rad))

    amber, prop_prm, prop_pi = {}, 0.0, 0.0
    phases = set()
    for d in p.dihedrals:
        a = [d.atom1, d.atom2, d.atom3, d.atom4]
        phi = dihedral([xyz[x.idx] for x in a])
        ph = math.radians(d.type.phase)
        phases.add(round(ph, 9))
        ph_exact = math.pi if abs(ph - math.pi) < 1e-3 else 0.0
        if d.improper:
            amber[frozenset(x.idx for x in a)] = dict(
                names=tuple(x.name for x in a), types=tuple(x.type for x in a),
                res=f"{a[2].residue.name}{a[2].residue.idx + 1}", k=d.type.phi_k, per=d.type.per,
                phi=phi, e=E_prm(d, phi, ph), e_pi=E_prm(d, phi, ph_exact))
        else:
            prop_prm += E_prm(d, phi, ph)
            prop_pi += E_prm(d, phi, ph_exact)
    imp_prm = sum(v["e"] for v in amber.values())

    # ---- MOSAICS build: terminal names, RTF order, RTF types, deck lookup ------
    def rtf_name(i):
        n = res[i].name
        return ("N" + n) if i == 0 else ("C" + n) if i == len(res) - 1 else n

    built = {}
    for i, r in enumerate(res):
        rn = rtf_name(i)
        for rec in rtf[rn]["impr"]:
            quad, types = [], []
            for tok in rec:
                if tok.startswith("-"):
                    j, nm = i - 1, tok[1:]
                elif tok.startswith("+"):
                    j, nm = i + 1, tok[1:]
                else:
                    j, nm = i, tok
                if j < 0 or j >= len(res):
                    quad = None
                    break
                hit = next((a for a in res[j].atoms if a.name == nm), None)
                if hit is None:
                    quad = None
                    break
                quad.append(hit)
                types.append(rtf[rtf_name(j)]["atoms"][nm])
            if not quad:
                continue
            ent = deck_lookup(deck, types)
            assert ent is not None, (rn, rec, types)
            phi = dihedral([xyz[a.idx] for a in quad])
            built[frozenset(a.idx for a in quad)] = dict(
                res=f"{r.name}{i + 1}", rtf_res=rn, rec=rec, types=tuple(types),
                k=float(ent["p0"]) / BOLTZ * KCAL / 2, phi=phi, e=mosaics_energy(ent, phi),
                names=tuple(a.name for a in quad))
    imp_mos = sum(v["e"] for v in built.values())

    W = 78
    print("=" * W)
    print("BPTI ff19SB dihedral term: MOSAICS-current vs OpenMM, closure")
    print("=" * W)
    print(f"prmtop phases present (rad)          : {sorted(phases)}")
    print(f"OpenMM PeriodicTorsionForce          : {omm_dihed:.9f}")
    print(f"analytic recompute, prmtop phase     : {prop_prm + imp_prm:.9f}   (diff {prop_prm + imp_prm - omm_dihed:+.3e})")
    print(f"   proper  {prop_prm:.9f}   improper {imp_prm:.9f}")
    print(f"MOSAICS printed                      : proper {mos_prop:.9f}   improper {mos_imp:.9f}")
    print(f"MOSAICS emulation (RTF order, deck)  : improper {imp_mos:.9f}   (emulation - printed {imp_mos - mos_imp:+.3e})")
    print()
    common = set(amber) & set(built)
    print(f"impropers: prmtop {len(amber)}, MOSAICS-built {len(built)}, common {len(common)}, "
          f"prmtop-only {len(set(amber) - set(built))}, built-only {len(set(built) - set(amber))}")
    kmis = [k for k in common if abs(amber[k]['k'] - built[k]['k']) > 1e-9]
    print(f"deck k != prmtop k                   : {len(kmis)}")
    order_diff = [k for k in common if amber[k]['names'] != built[k]['names']
                  and amber[k]['names'] != built[k]['names'][::-1]]
    print(f"atom order differs (RTF vs prmtop)   : {len(order_diff)}")
    for k in order_diff:
        a, b = amber[k], built[k]
        print(f"   {a['res']:6} prmtop {'-'.join(a['names']):14} types {'-'.join(a['types']):14} phi {math.degrees(a['phi']):9.4f}  E {a['e']:.9f}")
        print(f"   {'':6} RTF    {'-'.join(b['names']):14} types {'-'.join(b['types']):14} phi {math.degrees(b['phi']):9.4f}  E {b['e']:.9f}   k both {a['k']}")

    # ---- contributions ---------------------------------------------------------
    c_order = sum(built[k]['e'] - amber[k]['e'] for k in order_diff)
    c_phase_imp = sum(built[k]['e'] - amber[k]['e'] for k in common if k not in order_diff)
    c_phase_imp_check = sum(amber[k]['e_pi'] - amber[k]['e'] for k in common if k not in order_diff)
    c_phase_prop = prop_pi - prop_prm
    d_imp = mos_imp - imp_prm
    d_prop = mos_prop - prop_prm
    print()
    print("contribution                                                    kcal/mol")
    print("-" * W)
    print(f"I1  ASP3 N improper: RTF order (-C,H,N,CA) vs prmtop (C,CA,N,H)  {c_order:+.9f}")
    print(f"I2  phase pi (deck) vs 3.141594 rad (prmtop), other {len(common) - len(order_diff)} impropers {c_phase_imp:+.9f}")
    print(f"    (same, from prmtop k with exact-pi phase, cross-check)        {c_phase_imp_check:+.9f}")
    print(f"P1  phase pi/0 (deck) vs prmtop phase, 3791 proper terms          {c_phase_prop:+.9f}")
    print("-" * W)
    print(f"improper: MOSAICS - OpenMM = {d_imp:+.9f};  I1+I2 = {c_order + c_phase_imp:+.9f};  residual {d_imp - c_order - c_phase_imp:+.3e}")
    print(f"proper  : MOSAICS - OpenMM = {d_prop:+.9f};  P1    = {c_phase_prop:+.9f};  residual {d_prop - c_phase_prop:+.3e}")
    tot = d_imp + d_prop
    print(f"dihedral: MOSAICS - OpenMM = {tot:+.9f};  sum   = {c_order + c_phase_imp + c_phase_prop:+.9f};  residual {tot - c_order - c_phase_imp - c_phase_prop:+.3e}")

    # ---- reconcile the earlier analysis's 0.011867 -------------------------------
    print()
    print("reconciliation of the earlier 'unexplained 0.011867' (residue 58 resolved as ALA instead of CALA):")
    last = len(res) - 1
    r = res[last]
    byname = {a.name: a for a in r.atoms}
    prevC = next(a for a in res[last - 1].atoms if a.name == "C")
    oxt_key = frozenset(byname[n].idx for n in ("CA", "O", "C", "OXT"))
    e_oxt = built[oxt_key]["e"]
    q_ala = [prevC, byname["H"], byname["N"], byname["CA"]]     # ALA record  -C H N CA
    q_cala = [prevC, byname["CA"], byname["N"], byname["H"]]    # CALA record -C CA N H
    ent_ala = deck_lookup(deck, ["PROT-C", "PROT-H", "PROT-N", "PROT-CX"])
    k_prm = amber[frozenset(a.idx for a in q_cala)]["k"]
    phi_ala, phi_cala = dihedral([xyz[a.idx] for a in q_ala]), dihedral([xyz[a.idx] for a in q_cala])
    e_ala = k_prm * (1 + math.cos(2 * phi_ala - math.pi))
    e_cala = mosaics_energy(deck_lookup(deck, ["PROT-C", "PROT-CX", "PROT-N", "PROT-H"]), phi_cala)
    print(f"   ALA58 CA-O-C-OXT: present in RTF CALA and built by MOSAICS       {e_oxt:+.9f}  (earlier: counted as missing)")
    print(f"   ALA58 N improper with ALA record order  (-C H N CA), k={k_prm}   {e_ala:+.9f}  (earlier value; deck lookup for C-H-N-CX -> {ent_ala})")
    print(f"   ALA58 N improper with CALA record order (-C CA N H)              {e_cala:+.9f}  (what MOSAICS builds; = prmtop order)")
    print(f"   net error of the earlier bookkeeping = -OXT + (ALA - CALA)       {-e_oxt + (e_ala - e_cala):+.9f}")
    print(f"   earlier 'residual' 5.520414 - {mos_imp:.9f}                        {5.520414 - mos_imp:+.9f}")


if __name__ == "__main__":
    main()
