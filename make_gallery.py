"""Compose the eleven ray-traced renders into one labelled gallery figure.
Each render is cropped to its content and centred on a fixed 4:3 canvas so
that every panel has the same geometry."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.family']='DejaVu Sans'
import matplotlib.pyplot as plt, matplotlib.image as mpimg
MAN=json.load(open('data/MANIFEST.json'))
order=['1efs','1BNA','dna_acgt','2KOC','rna_acgu','1UBQ','5PTI','peptide','tax9','1HHK','1F7Y']
title={'1efs':('1efs','DNA:RNA hybrid duplex'),'1BNA':('1BNA','B-DNA duplex'),'dna_acgt':('d(ACGT)','DNA 4-mer, free ends'),
       '2KOC':('2KOC','UUCG tetraloop hairpin'),'rna_acgu':('r(ACGU)','RNA 4-mer, free ends'),'1UBQ':('1UBQ','ubiquitin'),
       '5PTI':('5PTI','BPTI, three disulfides'),'peptide':('ACE-AGSV-NME','capped tetrapeptide'),
       'tax9':('Tax 11-19','capped nonamer'),'1HHK':('1HHK','HLA-A*02:01 / β2m / Tax'),'1F7Y':('1F7Y','S15 + 16S rRNA')}
def framed(path, aspect=4/3):
    img=mpimg.imread(path)
    if img.shape[2]==4:
        mask=img[:,:,3]>0.02
        rgb=img[:,:,:3]*img[:,:,3:4]+(1-img[:,:,3:4])
    else:
        rgb=img[:,:,:3]; mask=rgb.sum(2)<2.97
    ys,xs=np.where(mask)
    pad=20
    y0,y1=max(ys.min()-pad,0),min(ys.max()+pad,img.shape[0]); x0,x1=max(xs.min()-pad,0),min(xs.max()+pad,img.shape[1])
    crop=rgb[y0:y1,x0:x1]
    h,w=crop.shape[:2]
    if w/h>aspect: W=w; H=int(round(w/aspect))
    else: H=h; W=int(round(h*aspect))
    canvas=np.ones((H,W,3))
    oy=(H-h)//2; ox=(W-w)//2
    canvas[oy:oy+h,ox:ox+w]=crop
    return canvas
fig=plt.figure(figsize=(7.0,5.6))
gs=fig.add_gridspec(3,4,wspace=0.05,hspace=0.32,left=0.01,right=0.99,top=0.96,bottom=0.02)
cells=[(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(2,2)]
for (r,c),s in zip(cells,order):
    ax=fig.add_subplot(gs[r,c])
    ax.imshow(framed(f'figures/systems/{s}.png'))
    ax.set_axis_off()
    nm,ds=title[s]
    ax.text(0.0,1.10,nm,transform=ax.transAxes,fontsize=7.5,fontweight='bold',color='#0b0b0b',va='bottom')
    ax.text(0.0,1.02,ds,transform=ax.transAxes,fontsize=6.3,color='#52514e',va='bottom')
    m=MAN[s]
    ax.text(0.0,-0.03,f"{m['atoms']:,} atoms, {m['residues']} residues, {len(m['chains'])} chain{'s' if len(m['chains'])>1 else ''}",
            transform=ax.transAxes,fontsize=6,color='#52514e',va='top')
ax=fig.add_subplot(gs[2,3]); ax.set_axis_off()
ax.text(0.02,0.98,'Colour by polymer',fontsize=7.5,fontweight='bold',va='top',color='#0b0b0b',transform=ax.transAxes)
for i,(lab,col) in enumerate([('DNA','#eb6834'),('RNA','#1baf7a'),('protein','#2a78d6'),('peptide','#4a3aa7')]):
    ax.add_patch(plt.Rectangle((0.04,0.74-i*0.16),0.11,0.10,color=col,transform=ax.transAxes,clip_on=False))
    ax.text(0.21,0.79-i*0.16,lab,fontsize=6.8,va='center',color='#0b0b0b',transform=ax.transAxes)
ax.text(0.02,0.02,'Coordinates as MOSAICS reads them;\nhydrogens present, not drawn.',fontsize=5.8,color='#52514e',transform=ax.transAxes,va='bottom')
fig.savefig('figures/systems-gallery.pdf',bbox_inches='tight',pad_inches=0.02)
fig.savefig('figures/systems-gallery.png',bbox_inches='tight',pad_inches=0.02,dpi=250)
print('ok')
