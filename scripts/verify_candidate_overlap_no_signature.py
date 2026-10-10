# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import csv,json,numpy as np
from pathlib import Path
z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz');lookup={x:i for i,x in enumerate(z['cells'])};C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']];results=json.loads(Path('work/cnv_no_signature_signature_results.json').read_text());expected=[]
for r in results:
 if 'balancing' not in r:continue
 rows=list(csv.DictReader((Path('work/cnv_no_signature')/r['patient']/'cell_labels.tsv').open(),delimiter='\t'));a=[lookup[x['cell']] for x in rows if x['tissue']=='Adjacent' and x['label']!='RNA_CNV_supported_candidate'];b=[lookup[x['cell']] for x in rows if x['tissue']=='Tumor' and x['label']=='RNA_CNV_supported_candidate'];ix=np.r_[a,b];t=np.r_[np.zeros(len(a)),np.ones(len(b))];x=C[ix];x=(x-x.mean(0))/np.where(x.std(0)>1e-12,x.std(0),1)
 np.savetxt('work/candidate_overlap_no_signature_'+r['patient']+'.tsv',np.c_[t,z['signature'][ix],x],delimiter='\t',header='T\ty\tx1\tx2\tx3\tx4\tx5',comments='');expected.append((r['patient'],r['balancing']['overlap']['difference'],r['regression']['difference']))
Path('work/candidate_overlap_no_signature_expected.tsv').write_text('patient\texpected\tregression_expected\n'+''.join(f'{p}\t{v}\t{beta}\n' for p,v,beta in expected))
