# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import numpy as np,json
from pathlib import Path
from scipy.optimize import minimize
from scipy.special import expit
r=json.loads(Path('work/balance_results.json').read_text());weights=np.load('work/balance_weights.npz')
for cohort,rows in r.items():
 if cohort=='GSE149614':
  z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz');outcome=z['signature']
 else:z=np.load('work/GSE189175_cluster_inputs.npz');c=np.load('work/GSE189175_scores.npz');outcome=c['repair']
 C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
 for row in rows:
  key=row['patient'] if cohort=='GSE149614' else 'sn'+row['patient'];ix=weights[key+'_index'];t=((z['sample'][ix]==key+'T') if cohort=='GSE149614' else z['tissue'][ix]=='Tumor').astype(float)
  X=(C[ix]-C[ix].mean(0))/C[ix].std(0);Z=np.c_[X,X**2];Z=(Z-Z.mean(0))/Z.std(0);D=np.c_[np.ones(len(t)),Z]
  def fg(b):
   e=D@b;return np.logaddexp(0,e).sum()-t@e,D.T@(expit(e)-t)
  f=minimize(fg,np.zeros(11),jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':3000});p=expit(D@f.x);w=np.where(t==1,1-p,p);a=t==0;b=t==1;wa=w[a]/w[a].sum();wb=w[b]/w[b].sum();err=abs(wb@Z[b]-wa@Z[a]).max();d={'difference':float(wb@outcome[ix][b]-wa@outcome[ix][a]),'ESS_N':float(1/(wa@wa)),'ESS_T':float(1/(wb@wb)),'max_standardized_moment_error':float(err),'gradient_max':float(abs(fg(f.x)[1]).max()),'optimizer_success':bool(f.success),'coefficient_norm':float(np.linalg.norm(f.x))};row['quadratic_overlap']=d
  print(cohort,key,d,flush=True)
Path('work/balance_results.json').write_text(json.dumps(r,indent=2))
md='\n## İkinci moment duyarlılığı\n\nBeş standardize kovaryatın kareleri eğilim skoruna eklendi. Böylece doğrusal ve karesel momentler dengelendi. Etkileşimler, ortak çok boyutlu destek ve ölçülmemiş bileşenler yine sınanmış değildir. Küçük ESS, daha fazla denge için bilgi kaybını gösterir.\n\n| Kohort | Hasta | Karesel örtüşme T−N | ESS N/T | En büyük moment farkı |\n|---|---|---:|---:|---:|\n'
for cohort,rows in r.items():
 for row in rows:
  q=row['quadratic_overlap'];md+=f"| {cohort} | {row['patient']} | {q['difference']:+.4f} | {q['ESS_N']:.1f}/{q['ESS_T']:.1f} | {q['max_standardized_moment_error']:.2e} |\n"
p=Path('outputs/kovaryat-dengeleme-v02.md');p.write_text(p.read_text()+md)
