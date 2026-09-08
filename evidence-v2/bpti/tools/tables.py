import json,re
from pathlib import Path
BASE=Path("/Users/bright/Library/CloudStorage/OneDrive-Personal/Documents/MOSAICS_AllAtom_Paper/evidence-v2/bpti")
KJ=4.184; TARGET=332.0637133
C={"sander":332.052200,"openmm":332.063713,"gromacs":332.063713,"mosaics":332.063441}
def fixc(v,e): return v*(TARGET/C[e]) if e in ("sander","mosaics") else v
def mosaics(p):
    d={}
    for line in Path(p).read_text().splitlines():
        s=line.strip()
        for key,tag in (("-Bond energy","bond"),("-Bend energy","angle"),
                        ("-Torsion energy","dihedral"),("-Total Lennard-Jones energy","lj"),
                        ("-Total Coulomb energy","coulomb"),("-Cmap energy","cmap"),
                        ("Tors-Tors energy","tors_proper"),("Tors-Improper energy","tors_improper"),
                        ("Onfo-Coulomb energy","onfo_coul"),("Onfo-Lennard-Jones energy","onfo_lj")):
            if s.startswith(key) and tag not in d: d[tag]=float(s.split()[-1])
    return d
def sander(p):
    blk=Path(p).read_text().split("NSTEP")[1]; v={}
    for m in re.finditer(r"(BOND|ANGLE|DIHED|VDWAALS|EEL|1-4 VDW|1-4 EEL|CMAP)\s*=\s*(-?[\d.]+)",blk):
        v.setdefault(m.group(1),float(m.group(2)))
    o={"bond":v["BOND"],"angle":v["ANGLE"],"dihedral":v["DIHED"],
       "lj":v["VDWAALS"]+v["1-4 VDW"],"coulomb":v["EEL"]+v["1-4 EEL"]}
    if "CMAP" in v: o["cmap"]=v["CMAP"]
    return o
def openmm(p):
    v={}
    for line in Path(p).read_text().splitlines():
        for key,tag in (("HarmonicBondForce","bond"),("HarmonicAngleForce","angle"),
                        ("PeriodicTorsionForce","dihedral"),("CMAPTorsionForce","cmap"),
                        ("Lennard-Jones (charges zeroed)","lj"),("Coulomb (epsilons zeroed)","coulomb")):
            if line.startswith(key): v[tag]=float(line.split()[-1])
    return v
def gmx(p):
    L=Path(p).read_text().splitlines()
    leg=[l.split('"')[1] for l in L if l.startswith("@ s")]
    v=dict(zip(leg,[float(x) for x in L[-1].split()][1:]))
    o={"bond":v["Bond"]/KJ,"angle":v["Angle"]/KJ,
       "dihedral":(v["Proper Dih."]+v["Per. Imp. Dih."])/KJ,
       "lj":(v["LJ-14"]+v["LJ (SR)"])/KJ,"coulomb":(v["Coulomb-14"]+v["Coulomb (SR)"])/KJ,
       "_proper":v["Proper Dih."]/KJ,"_improper":v["Per. Imp. Dih."]/KJ}
    if "CMAP Dih." in v: o["cmap"]=v["CMAP Dih."]/KJ
    return o
out={}
for ff in ("ff14sb","ff19sb"):
    d={"mosaics-current":mosaics(BASE/ff/"mosaics-current/mosaics_sp.out"),
       "sander":sander(BASE/ff/"sander/bpti_sander_sp.out"),
       "openmm":openmm(BASE/ff/"openmm/openmm_sp.out"),
       "gromacs":gmx(BASE/ff/"gromacs/sp.xvg")}
    terms=["bond","angle","dihedral"]+(["cmap"] if "cmap" in d["openmm"] else [])+["lj","coulomb"]
    L=[]
    L.append("%-12s %18s %18s %18s %18s"%("term","MOSAICS-current","sander","OpenMM","GROMACS"))
    for t in terms:
        L.append("%-12s %18.8f %18.8f %18.8f %18.8f"%(t,d["mosaics-current"][t],d["sander"][t],d["openmm"][t],d["gromacs"][t]))
    L.append("%12s(electrostatic values above are as printed)"%"")
    L.append("%-12s %18.8f %18.8f %18.8f %18.8f   corrected to one constant"%("coulomb*",
        fixc(d["mosaics-current"]["coulomb"],"mosaics"),fixc(d["sander"]["coulomb"],"sander"),
        d["openmm"]["coulomb"],d["gromacs"]["coulomb"]))
    tot=lambda e: sum(fixc(d[e][t],e) if t=="coulomb" else d[e][t] for t in terms)
    L.append("%-12s %18.8f %18.8f %18.8f %18.8f"%("total",tot("mosaics-current"),tot("sander"),tot("openmm"),tot("gromacs")))
    L.append("")
    L.append("MOSAICS-current minus OpenMM, after correction:")
    for t in terms:
        a=fixc(d["mosaics-current"][t],"mosaics") if t=="coulomb" else d["mosaics-current"][t]
        L.append("  %-10s %+.4e"%(t,a-d["openmm"][t]))
    L.append("MOSAICS-3.9.1 minus OpenMM: no run. 3.9.1 produced no energy on this system.")
    L.append("")
    L.append("sander minus OpenMM, after correction:")
    for t in terms:
        a=fixc(d["sander"][t],"sander") if t=="coulomb" else d["sander"][t]
        L.append("  %-10s %+.4e"%(t,a-d["openmm"][t]))
    L.append("GROMACS minus OpenMM:")
    for t in terms:
        L.append("  %-10s %+.4e"%(t,d["gromacs"][t]-d["openmm"][t]))
    L.append("")
    L.append("dihedral split (MOSAICS reports proper and improper separately;")
    L.append("GROMACS reports Proper Dih. and Per. Imp. Dih.; the table above sums them)")
    L.append("  %-22s MOSAICS-current %18.9f   GROMACS %18.9f"%("proper",d["mosaics-current"]["tors_proper"],d["gromacs"]["_proper"]))
    L.append("  %-22s MOSAICS-current %18.9f   GROMACS %18.9f"%("improper",d["mosaics-current"]["tors_improper"],d["gromacs"]["_improper"]))
    txt="\n".join(L)+"\n"
    (BASE/ff/"energy_table.txt").write_text(txt)
    out[ff]=d
    print("#"*30,ff); print(txt)
(BASE/"energies.json").write_text(json.dumps(out,indent=2)+"\n")
