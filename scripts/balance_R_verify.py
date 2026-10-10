# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import numpy as np,json
w=np.load('work/balance_weights.npz');r=json.loads(Path('work/balance_results.json').read_text())
for cohort,rows in r.items():
 if cohort=='GSE149614':z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz');y=z['signature']
 else:z=np.load('work/GSE189175_cluster_inputs.npz');c=np.load('work/GSE189175_scores.npz');y=c['repair']
 C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
 for row in rows:
  k=row['patient'] if cohort=='GSE149614' else 'sn'+row['patient'];ix=w[k+'_index'];t=((z['sample'][ix]==k+'T') if cohort=='GSE149614' else z['tissue'][ix]=='Tumor').astype(float)
  X=(C[ix]-C[ix].mean(0))/C[ix].std(0)
  np.savetxt('work/overlap_'+k+'.tsv',np.c_[t,y[ix],X],delimiter='\t',header='T\ty\tx1\tx2\tx3\tx4\tx5',comments='')
