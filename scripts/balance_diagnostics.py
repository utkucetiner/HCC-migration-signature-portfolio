# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import numpy as np,json
from pathlib import Path
r=json.loads(Path('work/balance_results.json').read_text());w=np.load('work/balance_weights.npz')
def calc(C,t,w):
 a=t==0;b=t==1;wa=w[a]/w[a].sum();wb=w[b]/w[b].sum();out=[]
 for j in range(C.shape[1]):
  ma=wa@C[a,j];mb=wb@C[b,j];va=wa@((C[a,j]-ma)**2);vb=wb@((C[b,j]-mb)**2)
  v=C[:,j];order=np.argsort(v);aa=np.where(t==0,w/w[a].sum(),0)[order];bb=np.where(t==1,w/w[b].sum(),0)[order]
  # Evaluate cumulative distribution at ends of ties, not arbitrary order inside tied scores.
  ends=np.r_[v[order][:-1]!=v[order][1:],True];ks=np.max(abs(np.cumsum(aa)[ends]-np.cumsum(bb)[ends]))
  out.append({'variance_ratio_T_N':float(vb/va),'weighted_KS':float(ks)})
 return out
names=['log_depth','S','G2M','identity','portal_minus_central']
for cohort,rows in r.items():
 if cohort=='GSE149614':
  z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz');C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
 else:
  z=np.load('work/GSE189175_cluster_inputs.npz');c=np.load('work/GSE189175_scores.npz');C=np.c_[np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']]
 for row in rows:
  key=row['patient'] if cohort=='GSE149614' else 'sn'+row['patient'];ix=w[key+'_index'];ww=w[key+'_weight'];t=(z['sample'][ix]==key+'T') if cohort=='GSE149614' else z['tissue'][ix]=='Tumor'
  row['distribution_diagnostics']=dict(zip(names,calc(C[ix],t,ww)))
Path('work/balance_results.json').write_text(json.dumps(r,indent=2))
p=Path('outputs/kovaryat-dengeleme-v02.md');md=p.read_text()+'\n## Ortalamaların ötesinde dağılım kontrolü\n\nAğırlıklı Kolmogorov–Smirnov uzaklığı betimleyici dağılım farkıdır; p değeri değildir. Eşit ortalamalar tam dağılım örtüşmesi sağlamaz.\n\n| Kohort | Hasta | En büyük ağırlıklı KS | En küçük / büyük varyans oranı T/N |\n|---|---|---:|---:|\n'
for cohort,rows in r.items():
 for row in rows:
  dd=list(row['distribution_diagnostics'].values());vr=[x['variance_ratio_T_N'] for x in dd];md+=f"| {cohort} | {row['patient']} | {max(x['weighted_KS'] for x in dd):.3f} | {min(vr):.2f}/{max(vr):.2f} |\n"
p.write_text(md)
