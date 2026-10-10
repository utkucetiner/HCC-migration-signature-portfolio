# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,numpy as np,csv
# Use the same fixed balancing definitions; no repair-dependent CNV label selection.
from balance_core import analyze
z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz');assert np.array_equal(z['cells'],c['cells']);lookup={x:i for i,x in enumerate(z['cells'])};C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']];results=[]
for n in range(3,11):
 p=f'HCC{n:02}';f=Path('work/cnv')/p/'cell_labels.tsv'
 if not f.exists():continue
 rows=list(csv.DictReader(f.open(),delimiter='\t'));assert len({r['cell'] for r in rows})==len(rows);assert all(r['cell'] in lookup for r in rows)
 a=np.array([lookup[r['cell']] for r in rows if r['tissue']=='Adjacent' and r['label']!='RNA_CNV_supported_candidate'],dtype=int);b=np.array([lookup[r['cell']] for r in rows if r['tissue']=='Tumor' and r['label']=='RNA_CNV_supported_candidate'],dtype=int)
 result={'patient':p,'comparison':'RNA-CNV-supported tumor hepatocytes versus adjacent hepatocytes excluding RNA-CNV-supported candidates; adjacent cells not DNA-verified benign','n_adjacent':len(a),'n_candidate_T':len(b),'labels':{t:{lab:sum(r['tissue']==t and r['label']==lab for r in rows) for lab in ['RNA_CNV_supported_candidate','CNV_signal_indeterminate','CNV_low_or_undetected']} for t in ['Adjacent','Tumor']}}
 if min(len(a),len(b))>=20:
  ix=np.r_[a,b];t=np.r_[np.zeros(len(a)),np.ones(len(b))];r,w,ps=analyze(C[ix],z['signature'][ix],t);result['balancing']=r
  weights=np.r_[np.ones(len(a))/len(a),np.ones(len(b))/len(b)];cov=C[ix]-np.average(C[ix],axis=0,weights=weights);sd=np.sqrt(np.average(cov**2,axis=0,weights=weights));valid=sd>1e-10;X=np.c_[np.ones(len(ix)),t,cov[:,valid]/sd[valid]];beta,_,rank,sv=np.linalg.lstsq(X*np.sqrt(weights)[:,None],z['signature'][ix]*np.sqrt(weights),rcond=None);result['regression']={'rank':int(rank),'columns':int(X.shape[1]),'constant_covariates_omitted':np.array(['log_depth','S','G2M','identity','zone'])[~valid].tolist(),'condition':float(sv[0]/sv[-1])};
  if rank==X.shape[1]:result['regression']['difference']=float(beta[1])
 results.append(result)
Path('work/cnv_signature_results.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
