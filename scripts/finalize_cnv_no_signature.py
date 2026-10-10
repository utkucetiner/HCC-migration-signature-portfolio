# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,csv,gzip,hashlib,importlib.util
import numpy as np
root=_WORKSPACE;target=root/'work/cnv_no_signature';o=root/'outputs'
patients=[f'HCC{i:02d}' for i in range(3,11)]
manifest=json.loads((target/'input_manifest.json').read_text());audit=[];all_labels=[]
for p in patients:
 d=target/p
 for mode in ['primary','hepatocyte_reference']:
  assert (d/mode/'COMPLETED').exists() and (d/mode/'final_infercnv.rds').stat().st_size>0
 rows=list(csv.DictReader((d/'cell_labels.tsv').open(),delimiter='\t'))
 old=list(csv.DictReader((root/'work/cnv'/p/'cell_labels.tsv').open(),delimiter='\t'))
 assert len(rows)==len(old) and {r['cell'] for r in rows}=={r['cell'] for r in old}
 for mode in ['primary','hepatocyte_reference']:
  with gzip.open(d/(mode+'_bins.tsv.gz'),'rb') as f:
   while f.read(1024*1024):pass
 audit.append({'patient':p,'hepatocytes':len(rows),'label_sha256':hashlib.sha256((d/'cell_labels.tsv').read_bytes()).hexdigest()});all_labels.extend(rows)
assert len(all_labels)==13762
with (o/'GSE149614_CNV_47_gen_disarida_hucre_etiketleri.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(all_labels[0]),delimiter='\t');w.writeheader();w.writerows(all_labels)
spec=importlib.util.spec_from_file_location('patient_uncertainty',_SOURCE_ROOT/'scripts/patient_uncertainty.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
new=json.loads((root/'work/cnv_no_signature_signature_results.json').read_text());old={r['patient']:r for r in json.loads((root/'work/cnv_signature_results.json').read_text())};assert len(new)==8
summary=[]
for name,func in [('overlap',lambda r:r.get('balancing',{}).get('overlap',{}).get('difference')),('zone_adjusted_regression',lambda r:r.get('regression',{}).get('difference'))]:
 valid=[r for r in new if func(r) is not None];summary.append(module.summarize('GSE149614_scRNA_CNV_without_signature',name,[r['patient'] for r in valid],[func(r) for r in valid],[r['patient'] for r in new if func(r) is None]))
comparison=json.loads((o/'CNV_47_gen_disarida_etiket_karsilastirmasi.json').read_text());assert comparison['completed_patients']==8
result={'status':'complete','patients':8,'hepatocytes':len(all_labels),'changed_labels':sum(r['changed_labels'] for r in comparison['results']),'old_candidates':sum(r['old_candidates'] for r in comparison['results']),'new_candidates':sum(r['new_candidates'] for r in comparison['results']),'summaries':summary,'verification':audit,'limits':'RNA-only robustness check, not DNA confirmation. Fixed within-patient effect bootstrap; classification/model uncertainty not fully propagated.'}
(o/'CNV_47_gen_disarida_nihai_ozet.json').write_text(json.dumps(result,indent=2))
p=o/'CNV_47_gen_disarida_ara_rapor.md';text=p.read_text();text=text.replace('Bu sayı sekize ulaşmadan rapor ara sonuçtur.','Sekiz hastanın tamamı sonuçlanmıştır.');text+='\n\n## Tam kohort özeti\n\n'
text+=f"13.762 hepatositin {result['changed_labels']} tanesinin etiketi değişti. Aday sayısı {result['old_candidates']} → {result['new_candidates']}.\n\n"
for s in summary:
 text+=f"{s['method']}: eşit hasta ağırlıklı ortalama {s['mean']:+.5f}; hasta bootstrap %95 aralığı {s['exact_patient_bootstrap_percentile_95']}; tek hasta çıkarılınca aralık [{s['loo_min']:+.5f}, {s['loo_max']:+.5f}].\n\n"
text+='Yeni aralıklar da sabit hasta içi farklara koşulludur; CNV, referans ve model seçiminin bütün belirsizliğini kapsamaz. CopyKAT bu aşamada çalıştırılmadı.\n';(o/'CNV_47_gen_disarida_nihai_rapor.md').write_text(text)
print(json.dumps({k:v for k,v in result.items() if k not in ['verification']},indent=2))
