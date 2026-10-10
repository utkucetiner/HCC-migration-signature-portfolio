# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,itertools
import numpy as np
root=_WORKSPACE;o=root/'outputs';d=json.loads((o/'hasta_duzeyi_belirsizlik.json').read_text());rng=np.random.default_rng(20261010)
checks=[]
for r in d['results']:
 x=np.array(r['patient_differences']);n=len(x)
 mc=x[rng.integers(0,n,size=(200000,n))].mean(axis=1)
 approximate=np.quantile(mc,[.025,.975]); exact=np.array(r['exact_patient_bootstrap_percentile_95'])
 assert np.max(abs(approximate-exact))<.0015
 for a in r['leave_one_patient_out']:
  keep=[v for p,v in zip(r['patients'],x) if p!=a['excluded_patient']]
  assert np.isclose(sum(keep)/len(keep),a['mean'],atol=1e-14)
 if n==3:
  full=np.array([x[list(ids)].mean() for ids in itertools.product(range(3),repeat=3)])
  assert np.allclose(np.quantile(full,[.025,.975],method='inverted_cdf'),exact)
 checks.append({'cohort':r['cohort'],'method':r['method'],'monte_carlo_max_endpoint_difference':float(max(abs(approximate-exact))),'loo_sign_flip_patients':[a['excluded_patient'] for a in r['leave_one_patient_out'] if np.sign(a['mean'])!=np.sign(r['mean'])]})
(o/'hasta_duzeyi_hesaplama_kontrolu.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
import os
os.environ['MPLCONFIGDIR']=str(root/'work/mplcache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
selected=[r for r in d['results'] if r['method'] in ['zone_adjusted_regression','overlap']]
fig,axes=plt.subplots(1,2,figsize=(12,6.5),gridspec_kw={'width_ratios':[1.1,1]})
labels=[]
for i,r in enumerate(selected):
 label=r['cohort'].replace('GSE149614_scRNA_CNV_candidates','scRNA — CNV adayları').replace('GSE149614_scRNA','scRNA — tüm hepatositler').replace('GSE189175_snRNA','snRNA — 3 hasta')+'\n'+('Düzeltilmiş regresyon' if r['method']=='zone_adjusted_regression' else 'Overlap ağırlığı')
 labels.append(label);lo,hi=r['exact_patient_bootstrap_percentile_95']; axes[0].plot([lo,hi],[i,i],color='#245b80',lw=2);axes[0].scatter([r['mean']],[i],color='#245b80',zorder=4)
 vals=[a['mean'] for a in r['leave_one_patient_out']];axes[1].scatter(vals,[i]*len(vals),color='#ad5931',s=30,alpha=.8);axes[1].scatter([r['mean']],[i],marker='|',color='black',s=170)
for ax in axes:
 ax.axvline(0,color='gray',ls='--');ax.set_ylim(len(selected)-.5,-.5);ax.set_xlabel('İmza puanı farkı (tümör − komşu)');ax.spines[['top','right']].set_visible(False)
axes[0].set_yticks(range(len(selected)),labels,fontsize=9);axes[1].set_yticks(range(len(selected)),['']*len(selected));axes[0].set_title('Hasta bootstrap %95 aralıkları');axes[1].set_title('Bir hasta çıkarılınca ortalamalar')
fig.suptitle('Hasta düzeyindeki belirsizlik ve duyarlılık',fontsize=14)
fig.text(.5,.015,'Hastalar eşit ağırlıklı; kohortlar ayrı. CNV ve hasta içi model belirsizliği yeniden hesaplanmadı.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.045,1,.95]);fig.savefig(o/'hasta_duzeyi_belirsizlik.png',dpi=180);fig.savefig(o/'hasta_duzeyi_belirsizlik.pdf');plt.close(fig)
