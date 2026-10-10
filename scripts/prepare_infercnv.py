# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,csv,json,re
from pathlib import Path
import numpy as np
from scipy import sparse,io
root=Path('work/cnv');root.mkdir(exist_ok=True)
with gzip.open('work/data/GSE149614_HCC.metadata.updated.txt.gz','rt') as f:meta={r['Cell']:r for r in csv.DictReader(f,delimiter='\t')}
rng=np.random.default_rng(20261010);selected=[];groups=[];manifest={}
for n in range(3,11):
 p=f'HCC{n:02}';heps=[c for c,r in meta.items() if r['sample'] in [p+'N',p+'T'] and r['celltype']=='Hepatocyte'];refs=[]
 for ct in ['T/NK','B','Endothelial']:
  candidates=[c for c,r in meta.items() if r['sample']==p+'N' and r['celltype']==ct]
  chosen=rng.choice(candidates,min(100,len(candidates)),replace=False).tolist();refs+=chosen
 selected+=heps+refs;groups+=[p]*len(heps+refs);manifest[p]={'hepatocytes':len(heps),'reference_cells':len(refs),'references_by_type':{ct:sum(meta[c]['celltype']==ct for c in refs) for ct in ['T/NK','B','Endothelial']}}
selected=set(selected)
with gzip.open('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rt') as f:
 cells=f.readline().strip().split('\t');idx=np.array([i for i,c in enumerate(cells) if c in selected]);cells=np.array(cells)[idx];parts=[];cols=[];indptr=[0];genes=[]
 for k,line in enumerate(f):
  gene,rest=line.rstrip().split('\t',1);a=np.fromstring(rest,sep='\t',dtype=np.int32)[idx];nz=np.flatnonzero(a);parts.append(a[nz]);cols.append(nz.astype(np.int32));indptr.append(indptr[-1]+len(nz));genes.append(gene)
  if (k+1)%5000==0:print('raw genes',k+1,flush=True)
 M=sparse.csr_matrix((np.concatenate(parts),np.concatenate(cols),np.array(indptr)),shape=(len(genes),len(cells)))
np.savez_compressed(root/'cells_genes.npz',cells=cells,genes=np.array(genes));sparse.save_npz(root/'raw_counts.npz',M)
for p in manifest:
 mask=np.array([meta[c]['patient']==p for c in cells]);d=root/p;d.mkdir(exist_ok=True)
 with gzip.open(d/'counts.mtx.gz','wb') as f:io.mmwrite(f,M[:,mask])
 (d/'genes.txt').write_text('\n'.join(genes)+'\n');(d/'cells.txt').write_text('\n'.join(cells[mask])+'\n')
 annotations=[]
 for c in cells[mask]:
  r=meta[c];label=('hep_'+r['site'] if r['celltype']=='Hepatocyte' else 'ref_'+r['celltype'].replace('/','_'))
  annotations.append(c+'\t'+label)
 (d/'annotations.tsv').write_text('\n'.join(annotations)+'\n')
Path('work/cnv/input_manifest.json').write_text(json.dumps(manifest,indent=2));print(manifest,flush=True)
