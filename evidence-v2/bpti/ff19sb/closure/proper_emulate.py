#!/usr/bin/env python3
"""Emulate MOSAICS proper-torsion energy from the deck to attribute the 8e-8 proper residual.
For every unique proper quad in the prmtop (MOSAICS autogenerates the same quads from bonds),
look up the deck entry by RTF atom types (tors.cpp last-match rule) and evaluate tors.cpp:756-779.
Usage: proper_emulate.py <prmtop> <rst7> <rtf> <deck>
"""
import sys, math, re, pathlib
import numpy as np, parmed
sys.path.insert(0, "/tmp/bpti19")
from closure import dihedral, read_rtf, read_deck, wild_eq, mosaics_energy, BOLTZ, KCAL

prm, rst, rtf_p, deck_p = sys.argv[1:5]
p = parmed.load_file(prm, xyz=rst)
rtf = read_rtf(rtf_p); deck = read_deck(deck_p)
res = list(p.residues)
xyz = {a.idx: np.array([a.xx, a.xy, a.xz]) for a in p.atoms}
def rtf_name(i):
    n = res[i].name
    return ("N"+n) if i == 0 else ("C"+n) if i == len(res)-1 else n
atype = {}
for i, r in enumerate(res):
    for a in r.atoms:
        atype[a.idx] = rtf[rtf_name(i)]["atoms"][a.name]

def lookup_proper(types):
    exact = wild = None
    for e in deck:
        if e.get("label", "").lower() == "improper":
            continue
        lbl = [e["atom1"], e["atom2"], e["atom3"], e["atom4"]]
        if all(types[i].upper() == lbl[i].upper() for i in range(4)) or all(types[i].upper() == lbl[3-i].upper() for i in range(4)):
            exact = e
        if all(wild_eq(types[i], lbl[i]) for i in range(4)) or all(wild_eq(types[i], lbl[3-i]) for i in range(4)):
            wild = e
    return exact if exact is not None else wild

quads = {}
for d in p.dihedrals:
    if d.improper: continue
    key = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
    if key[::-1] in quads: key = key[::-1]
    quads.setdefault(key, []).append(d)
print("unique proper quads:", len(quads), " prmtop proper terms:", sum(len(v) for v in quads.values()))
e_mos = e_prm = e_pi = 0.0
miss = 0
worst = []
for key, ds in quads.items():
    phi = dihedral([xyz[i] for i in key])
    types = [atype[i] for i in key]
    ent = lookup_proper(types)
    if ent is None:
        miss += 1; continue
    em = mosaics_energy(ent, phi)
    ep = sum(d.type.phi_k*(1+math.cos(d.type.per*phi-math.radians(d.type.phase))) for d in ds)
    epi = sum(d.type.phi_k*(1+math.cos(d.type.per*phi-(math.pi if d.type.phase > 90 else 0.0))) for d in ds)
    e_mos += em; e_prm += ep; e_pi += epi
    worst.append((abs(em-epi), key, types, em, epi, ep))
worst.sort(reverse=True)
print(f"deck-emulated MOSAICS proper : {e_mos:.9f}   (unmatched quads {miss})")
print(f"prmtop, prmtop phase         : {e_prm:.9f}")
print(f"prmtop, exact pi/0 phase     : {e_pi:.9f}")
print(f"emulated - printed 224.973266426 : {e_mos-224.973266426:+.3e}")
print(f"emulated - exact-pi analytic     : {e_mos-e_pi:+.3e}   <- deck coefficient rounding")
print("worst 5 |emulated - exact-pi| quads:")
for w, key, types, em, epi, ep in worst[:5]:
    print(f"   {'-'.join(types):40} em={em:.9f} epi={epi:.9f} d={em-epi:+.2e}")
