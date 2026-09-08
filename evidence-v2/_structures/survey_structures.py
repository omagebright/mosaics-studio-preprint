import sys, collections
STD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
          'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
# heavy-atom templates (no OXT, no H)
AA_HEAVY = {
 'ALA':'N CA C O CB','ARG':'N CA C O CB CG CD NE CZ NH1 NH2','ASN':'N CA C O CB CG OD1 ND2',
 'ASP':'N CA C O CB CG OD1 OD2','CYS':'N CA C O CB SG','GLN':'N CA C O CB CG CD OE1 NE2',
 'GLU':'N CA C O CB CG CD OE1 OE2','GLY':'N CA C O','HIS':'N CA C O CB CG ND1 CD2 CE1 NE2',
 'ILE':'N CA C O CB CG1 CG2 CD1','LEU':'N CA C O CB CG CD1 CD2','LYS':'N CA C O CB CG CD CE NZ',
 'MET':'N CA C O CB CG SD CE','PHE':'N CA C O CB CG CD1 CD2 CE1 CE2 CZ','PRO':'N CA C O CB CG CD',
 'SER':'N CA C O CB OG','THR':'N CA C O CB OG1 CG2','TRP':'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
 'TYR':'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH','VAL':'N CA C O CB CG1 CG2'}
NUC_BASE = {'DA':'N9 C8 N7 C5 C6 N6 N1 C2 N3 C4','DG':'N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4',
 'DC':'N1 C2 O2 N3 C4 N4 C5 C6','DT':'N1 C2 O2 N3 C4 O4 C5 C7 C6',
 'A':'N9 C8 N7 C5 C6 N6 N1 C2 N3 C4','G':'N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4',
 'C':'N1 C2 O2 N3 C4 N4 C5 C6','U':'N1 C2 O2 N3 C4 O4 C5 C6'}
DNA_SUG="P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' C1'"
RNA_SUG="P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1'"

def survey(path, chain, model=None):
    res=collections.OrderedDict(); inmodel = model is None
    for L in open(path):
        if L.startswith('MODEL'):
            inmodel = (model is not None and int(L.split()[1])==model)
        elif L.startswith('ENDMDL') and model is not None:
            inmodel=False
        elif L.startswith('ATOM') and inmodel:
            if L[21]!=chain: continue
            alt=L[16]
            if alt not in (' ','A'): continue
            key=(L[22:27].strip(), L[17:20].strip())
            res.setdefault(key,[]).append(L[12:16].strip())
    return res

for pid,ch,mdl in [('1UBQ','A',None),('5PTI','A',None),('1BNA','A',None),('1BNA','B',None),('2KOC','A',1)]:
    p=f'/Users/bright/Library/CloudStorage/OneDrive-Personal/Documents/MOSAICS_AllAtom_Paper/evidence-v2/_structures/raw/{pid}.pdb'
    r=survey(p,ch,mdl)
    print(f'===== {pid} chain {ch}' + (f' model {mdl}' if mdl else '') + f'  residues={len(r)} =====')
    comp=collections.Counter(k[1] for k in r)
    print('  composition:', ' '.join(f'{k}:{v}' for k,v in sorted(comp.items())))
    print('  n_atoms(altloc A/blank):', sum(len(v) for v in r.values()))
    nums=[int(''.join(c for c in k[0] if c.isdigit() or c=='-')) for k in r]
    gaps=[(nums[i],nums[i+1]) for i in range(len(nums)-1) if nums[i+1]!=nums[i]+1]
    print('  numbering:', nums[0],'to',nums[-1], '| gaps:', gaps if gaps else 'none')
    nonstd=[k for k in r if k[1] not in STD_AA and k[1] not in NUC_BASE]
    print('  non-standard residues:', nonstd if nonstd else 'none')
    miss=[]
    for (num,name),atoms in r.items():
        if name in AA_HEAVY: tmpl=AA_HEAVY[name].split()
        elif name in ('DA','DG','DC','DT'): tmpl=(DNA_SUG+' '+NUC_BASE[name]).split()
        elif name in ('A','G','C','U'): tmpl=(RNA_SUG+' '+NUC_BASE[name]).split()
        else: continue
        have=set(a for a in atoms if not a.startswith('H') and not a[0].isdigit() and a!='OXT')
        # tolerate 5'-terminal phosphate absence
        m=[a for a in tmpl if a not in have]
        if m: miss.append((num,name,m))
    print('  residues missing heavy atoms:', miss if miss else 'none')
    hyd=sum(1 for v in r.values() for a in v if a.startswith('H') or (a[0].isdigit() and 'H' in a))
    print('  hydrogens present:', hyd)
