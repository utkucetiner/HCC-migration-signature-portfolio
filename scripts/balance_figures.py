# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import os
from pathlib import Path
os.environ['MPLCONFIGDIR']=str(Path('work/mplcache').resolve())
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,json,numpy as np
r=json.loads(Path('work/balance_results.json').read_text());fig,axs=plt.subplots(1,2,figsize=(10,4.5),gridspec_kw={'width_ratios':[2,1]})
for ax,(cohort,rows) in zip(axs,r.items()):
 x=np.arange(len(rows));ax.axhline(0,color='gray',lw=.7)
 for offset,method,label,color in [(-.16,'raw','Raw','#999999'),(0,'overlap','Overlap weights','#2166ac'),(.16,'entropy','Entropy balance','#b2182b')]:
  vals=[a.get(method,{}).get('difference',np.nan) for a in rows];ax.scatter(x+offset,vals,label=label,color=color,s=30)
 ax.set_xticks(x,[a['patient'].replace('HCC','') for a in rows]);ax.set_xlabel('Patient');ax.set_title(cohort+(' (scRNA-seq)' if cohort=='GSE149614' else ' (snRNA-seq)'));ax.set_ylabel('Weighted signature difference (T − N)');ax.spines[['top','right']].set_visible(False)
axs[0].legend(fontsize=8);fig.tight_layout();fig.savefig('outputs/kovaryat_dengeleme.png',dpi=180);fig.savefig('outputs/kovaryat_dengeleme.pdf');plt.close(fig)
