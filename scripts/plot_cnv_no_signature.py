# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,os
import numpy as np
root=_WORKSPACE;o=root/'outputs';os.environ['MPLCONFIGDIR']=str(root/'work/mplcache')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
labels=json.loads((o/'CNV_47_gen_disarida_etiket_karsilastirmasi.json').read_text())['results'];old={r['patient']:r for r in json.loads((root/'work/cnv_signature_results.json').read_text())};new={r['patient']:r for r in json.loads((root/'work/cnv_no_signature_signature_results.json').read_text())};assert len(labels)==8
fig,axs=plt.subplots(1,2,figsize=(11,5));y=np.arange(8)
axs[0].barh(y,[100*r['candidate_jaccard'] for r in labels],color='#245b80');axs[0].set_yticks(y,[r['patient'] for r in labels]);axs[0].invert_yaxis();axs[0].set_xlim(0,105);axs[0].set_xlabel('Aday kümelerinin Jaccard uyumu (%)');axs[0].set_title('Eski–yeni RNA-CNV adayları')
for i,r in enumerate(labels):
 p=r['patient'];a=old[p]['balancing']['overlap']['difference'];b=new[p]['balancing']['overlap']['difference'];axs[1].plot([a,b],[i,i],color='gray');axs[1].scatter(a,i,color='#245b80',label='Özgün CNV' if i==0 else None);axs[1].scatter(b,i,color='#ad5931',label='47 gen çıkarılmış' if i==0 else None)
axs[1].set_yticks(y,[r['patient'] for r in labels]);axs[1].invert_yaxis();axs[1].axvline(0,color='gray',ls='--');axs[1].set_xlabel('Overlap ağırlıklı imza farkı (tümör − komşu)');axs[1].set_title('İmza farkının duyarlılığı');axs[1].legend(fontsize=9)
for ax in axs:ax.spines[['top','right']].set_visible(False)
fig.suptitle('47 imza geni çıkarılmış inferCNV: sekiz hasta',fontsize=14);fig.text(.5,.02,'RNA temelli sağlamlık kontrolü; bağımsız DNA doğrulaması değildir.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(o/'CNV_47_gen_disarida_karsilastirma.png',dpi=180);fig.savefig(o/'CNV_47_gen_disarida_karsilastirma.pdf');plt.close(fig)
