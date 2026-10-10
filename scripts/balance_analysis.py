# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from balance_core import analyze
import numpy as np,json
from pathlib import Path
z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz')
C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
results={'GSE149614':[],'GSE189175':[]};weights={}
for n in range(3,11):
 patient=f'HCC{n:02}';ix=np.flatnonzero(np.isin(z['sample'],[patient+'N',patient+'T']));t=(z['sample'][ix]==patient+'T').astype(float)
 r,w,p=analyze(C[ix],z['signature'][ix],t);r['patient']=patient;results['GSE149614'].append(r);weights[patient+'_index']=ix;weights[patient+'_weight']=w;weights[patient+'_ps']=p
z=np.load('work/GSE189175_cluster_inputs.npz');c=np.load('work/GSE189175_scores.npz');C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
for patient in ['1','2','3']:
 ix=np.flatnonzero(c['coremask']&(z['patients']==patient));t=(z['tissue'][ix]=='Tumor').astype(float)
 r,w,p=analyze(C[ix],c['repair'][ix],t);r['patient']=patient;results['GSE189175'].append(r);weights['sn'+patient+'_index']=ix;weights['sn'+patient+'_weight']=w;weights['sn'+patient+'_ps']=p
Path('work/balance_results.json').write_text(json.dumps(results,indent=2));np.savez_compressed('work/balance_weights.npz',**weights)
for cohort,rr in results.items():
 print(cohort)
 for r in rr:print(r['patient'],r['overlap']['difference'],round(r['overlap']['ESS_N'],1),round(r['overlap']['ESS_T'],1),max(abs(v) for v in r['overlap']['SMD'].values()),'entropy',r.get('entropy',{}).get('difference'),r['logistic_gradient_max'],flush=True)
