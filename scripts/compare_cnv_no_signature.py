# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import csv,json
root=_WORKSPACE;results=[]
for i in range(3,11):
 p=f'HCC{i:02d}';old=root/'work/cnv'/p/'cell_labels.tsv';new=root/'work/cnv_no_signature'/p/'cell_labels.tsv'
 if not new.exists():continue
 read=lambda f:{r['cell']:r for r in csv.DictReader(f.open(),delimiter='\t')}
 a=read(old);b=read(new);assert set(a)==set(b)
 transition={}
 for cell,r in a.items():
  key=r['label']+' -> '+b[cell]['label'];transition[key]=transition.get(key,0)+1
 old_c={c for c in a if a[c]['label']=='RNA_CNV_supported_candidate'};new_c={c for c in b if b[c]['label']=='RNA_CNV_supported_candidate'}
 results.append({'patient':p,'hepatocytes':len(a),'unchanged_labels':sum(a[c]['label']==b[c]['label'] for c in a),'changed_labels':sum(a[c]['label']!=b[c]['label'] for c in a),'old_candidates':len(old_c),'new_candidates':len(new_c),'candidate_jaccard':len(old_c&new_c)/len(old_c|new_c) if old_c|new_c else None,'transitions':transition})
(root/'outputs/CNV_47_gen_disarida_etiket_karsilastirmasi.json').write_text(json.dumps({'completed_patients':len(results),'expected_patients':8,'results':results},indent=2));print(json.dumps(results,indent=2))
