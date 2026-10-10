# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,csv,os
from pathlib import Path
os.environ['MPLCONFIGDIR']=str(Path('work/mplcache').resolve())
import numpy as np,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
p='HCC06';lab=list(csv.DictReader(open(f'work/cnv/{p}/cell_labels.tsv'),delimiter='\t'));lab.sort(key=lambda x:(x['tissue'],x['label'],x['cnv_group']));cells=[x['cell'] for x in lab]
fig,axs=plt.subplots(1,2,figsize=(13,5),sharey=True)
for ax,mode in zip(axs,['primary','hepatocyte_reference']):
 with gzip.open(f'work/cnv/{p}/{mode}_bins.tsv.gz','rt') as f:
  hdr=next(f).strip().split('\t');rows=list(csv.reader(f,delimiter='\t'))
 # First column is empty in header then cell identifiers.
 if len(hdr)==len(rows[0])-1:columns=hdr
 else:columns=hdr[1:]
 order=[columns.index(c) for c in cells];bins=[r[0] for r in rows];M=np.array([[float(v) for v in r[1:]] for r in rows])[:,order]
 b_order=sorted(range(len(bins)),key=lambda i:tuple(map(int,bins[i].replace('chr','').split(':'))));M=M[b_order];bb=[bins[i] for i in b_order]
 im=ax.imshow(np.log2(np.maximum(M.T,1e-6)),aspect='auto',cmap='RdBu_r',vmin=-.4,vmax=.4,interpolation='nearest');starts=[i for i,b in enumerate(bb) if i==0 or b.split(':')[0]!=bb[i-1].split(':')[0]];selected=[i for i in starts if int(bb[i].split(':')[0][3:])<=12 or int(bb[i].split(':')[0][3:])%2==1];ax.set_xticks(selected,[bb[i].split(':')[0].replace('chr','') for i in selected],fontsize=7);ax.set_xlabel('Autosome (10-Mb gene bins)');ax.set_title('Non-hepatocyte references' if mode=='primary' else 'Adjacent hepatocyte reference');ax.axhline(sum(x['tissue']=='Adjacent' for x in lab)-.5,color='black',lw=.8)
 for i in range(1,len(lab)):
  if (lab[i]['tissue'],lab[i]['label'])!=(lab[i-1]['tissue'],lab[i-1]['label']):ax.axhline(i-.5,color='black',lw=.8,linestyle='--')
axs[0].set_ylabel('HCC06 hepatocytes; adjacent first, tumor second');fig.subplots_adjust(left=.06,right=.87,bottom=.13,top=.9,wspace=.12);cax=fig.add_axes([.9,.25,.014,.5]);fig.colorbar(im,cax=cax,label='log2 relative RNA-CNV signal');fig.savefig('outputs/HCC06_inferCNV_referans_duyarliligi.png',dpi=170);plt.close(fig)
