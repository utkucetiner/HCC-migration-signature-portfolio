# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,csv
root=_WORKSPACE;o=root/'outputs'
a=json.loads((o/'CNV_47_gen_disarida_etiket_karsilastirmasi.json').read_text());old={r['patient']:r for r in json.loads((root/'work/cnv_signature_results.json').read_text())};f=root/'work/cnv_no_signature_signature_results.json';new={r['patient']:r for r in json.loads(f.read_text())} if f.exists() else {}
lines=['# 47 imza geni çıkarılmış inferCNV: karşılaştırma raporu','',f"Tamamlanan hasta: {a['completed_patients']}/8. Bu sayı sekize ulaşmadan rapor ara sonuçtur.",'','47 gen sayım matrisinden ve gen sırası girdisinden çıkarıldı. Hücreler, referans seçimleri ve CNV ayarları korundu. İmza puanı özgün 47 genle hesaplanmaya devam eder. Böylece imza değişimi değil, CNV seçiminin duyarlılığı sınanır. Gen çıkarma CNV normalizasyon paydasını da değiştirir.','','| Hasta | Hepatosit | Aynı etiket | Değişen etiket | Eski aday | Yeni aday | Aday Jaccard | Eski overlap farkı | Yeni overlap farkı | Yeni ESS N / T |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
for r in a['results']:
 p=r['patient'];b=new.get(p,{}).get('balancing',{}).get('overlap',{});previous=old[p]['balancing']['overlap']['difference'];v=b.get('difference');delta='bekliyor' if v is None else f'{v:+.5f}';ess='bekliyor' if v is None else f"{b['ESS_N']:.1f} / {b['ESS_T']:.1f}"
 lines.append(f"| {p} | {r['hepatocytes']} | {r['unchanged_labels']} | {r['changed_labels']} | {r['old_candidates']} | {r['new_candidates']} | {r['candidate_jaccard']:.4f} | {previous:+.5f} | {delta} | {ess} |")
lines+=['','## Sınırlar','','Etiket uyumu bağımsız DNA doğrulaması değildir. Referans doku etkileri, genler arasındaki korelasyon ve RNA temelli sınıflamanın belirsizliği devam eder. Özgün imza puanının normalizasyonunda tüm transkriptom sayımları kullanılmıştır; CNV için gen dışlama ortak ölçüm kaynağının bütün bağımlılıklarını ortadan kaldırmaz.','','Etkin örneklem büyüklüğü 20 altına düşen gruplar özellikle dikkatle yorumlanmalıdır. İşaret korunması istatistiksel kesinlik değildir. Tam kohort özetini yalnızca sekiz hastanın tamamı sonuçlandığında hesaplamak gerekir.']
(o/'CNV_47_gen_disarida_ara_rapor.md').write_text('\n'.join(lines)+'\n');print(f"Report updated: {a['completed_patients']}/8")
