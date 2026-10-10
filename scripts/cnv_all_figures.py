# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import csv,gzip,os
from pathlib import Path
os.environ['MPLCONFIGDIR']=str(Path('work/mplcache').resolve())
import numpy as np,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(2,4,figsize=(13,7));rng=np.random.default_rng(20261010)
for ax,n in zip(axs.flat,range(3,11)):
 p=f'HCC{n:02}';d=Path('work/cnv')/p;lab=list(csv.DictReader((d/'cell_labels.tsv').open(),delimiter='\t'));lab.sort(key=lambda x:(x['tissue'],x['label'],x['cnv_group']))
 # Preserve tissue/label/CNV-group ordering, deterministic sampling only for display.
 display=[]
 for tissue in ['Adjacent','Tumor']:
  members=[x for x in lab if x['tissue']==tissue]
  chosen=set(rng.choice(len(members),min(150,len(members)),replace=False));display.extend(x for i,x in enumerate(members) if i in chosen)
 cells=[x['cell'] for x in display]
 with gzip.open(d/'hepatocyte_reference_bins.tsv.gz','rt') as f:
  hdr=next(f).strip().split('\t');rows=list(csv.reader(f,delimiter='\t'))
 columns=hdr if len(hdr)==len(rows[0])-1 else hdr[1:];order=[columns.index(c) for c in cells];bins=[r[0] for r in rows];M=np.array([[float(v) for v in r[1:]] for r in rows])[:,order];bi=sorted(range(len(bins)),key=lambda i:tuple(map(int,bins[i].replace('chr','').split(':'))));bb=[bins[i] for i in bi]
 im=ax.imshow(np.log2(np.maximum(M[bi].T,1e-6)),aspect='auto',cmap='RdBu_r',vmin=-.4,vmax=.4,interpolation='nearest');starts=[i for i,b in enumerate(bb) if (i==0 or b.split(':')[0]!=bb[i-1].split(':')[0]) and int(b.split(':')[0][3:]) in [1,4,8,13,17,21]];ax.set_xticks(starts,[bb[i].split(':')[0].replace('chr','') for i in starts],fontsize=8);a=sum(x['tissue']=='Adjacent' for x in display);ax.axhline(a-.5,color='black',lw=.7);ax.set_yticks([max(0,(a-1)/2),a+(len(display)-a-1)/2],['Adjacent','Tumor'],fontsize=8);nt=sum(x['tissue']=='Tumor' for x in lab);nc=sum(x['tissue']=='Tumor' and x['label']=='RNA_CNV_supported_candidate' for x in lab);ax.set_title(f'{p}  |  candidates {nc}/{nt}',fontsize=10)
 for i in range(1,len(display)):
  if (display[i]['tissue'],display[i]['label'])!=(display[i-1]['tissue'],display[i-1]['label']):ax.axhline(i-.5,color='black',lw=.5,linestyle='--')
 ax.set_xlabel('Autosome',fontsize=8)
fig.subplots_adjust(left=.06,right=.9,top=.94,bottom=.07,hspace=.32,wspace=.3);cax=fig.add_axes([.93,.25,.012,.5]);fig.colorbar(im,cax=cax,label='log2 relative RNA-CNV signal');fig.savefig('outputs/GSE149614_inferCNV_sekiz_hasta.png',dpi=180);fig.savefig('outputs/GSE149614_inferCNV_sekiz_hasta.pdf');plt.close(fig)
