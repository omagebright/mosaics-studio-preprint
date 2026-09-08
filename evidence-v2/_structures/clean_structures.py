import os
BASE='/Users/bright/Library/CloudStorage/OneDrive-Personal/Documents/MOSAICS_AllAtom_Paper/evidence-v2/_structures'
RAW=os.path.join(BASE,'raw')

def clean(pid,out,chains,model=None,drop_atoms=None,strip_h=True,drop_oxt=True):
    drop_atoms = drop_atoms or {}
    lines=[]; inm = model is None; n=0
    for L in open(os.path.join(RAW,pid+'.pdb')):
        if L.startswith('MODEL'):
            inm = (model is not None and int(L.split()[1])==model); continue
        if L.startswith('ENDMDL'):
            if model is not None: inm=False
            continue
        if not L.startswith('ATOM'): continue          # drops every HETATM: waters, DOD, PO4, UNX, ions
        if not inm: continue
        if L[21] not in chains: continue
        if L[16] not in (' ','A'): continue            # altloc: keep A / blank
        name=L[12:16].strip(); elem=L[76:78].strip()
        if strip_h and (elem in ('H','D') or (name[0] in 'HD' and elem=='')): continue
        if drop_oxt and name=='OXT': continue
        key=(L[21], L[22:27].strip())
        if name in drop_atoms.get(key,()): continue
        n+=1
        lines.append(L[:16]+' '+L[17:21]+L[21]+L[22:]) # blank the altloc flag
    with open(os.path.join(BASE,out),'w') as f:
        f.writelines(lines); f.write('END\n')
    return n

print('1UBQ_A.pdb      ', clean('1UBQ','1UBQ_A.pdb',{'A'}), 'atoms')
print('5PTI_A.pdb      ', clean('5PTI','5PTI_A.pdb',{'A'}), 'atoms')
print('1BNA_AB.pdb     ', clean('1BNA','1BNA_AB.pdb',{'A','B'}), 'atoms')
# 2KOC: model 1, and remove the 5'-terminal phosphate on G1 so the 5' end is an OH
print('2KOC_A_model1.pdb', clean('2KOC','2KOC_A_model1.pdb',{'A'},model=1,
      drop_atoms={('A','1'):{'OP3','P','OP1','OP2'}}), 'atoms')
