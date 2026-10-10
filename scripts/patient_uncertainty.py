"""Exact patient-resampling bootstrap of fixed within-patient contrasts.
No cell independence, no pooled cohorts, no CNV/selection refit implied.
"""
# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json, math, csv, hashlib
import numpy as np
from patient_core import compositions, summarize
ROOT=_WORKSPACE
W=ROOT/'work'; O=ROOT/'outputs'; O.mkdir(exist_ok=True)





expanded=json.loads((W/'expanded_control_results.json').read_text())['patients']
bal=json.loads((W/'balance_results.json').read_text())
cnv=json.loads((W/'cnv_signature_results.json').read_text())
sn=json.loads((W/'GSE189175_score_results.json').read_text())['patients']
results=[]
def add(cohort,method,rows,fn):
    available=[];excluded=[]
    for r in rows:
        v=fn(r)
        (excluded if v is None else available).append(r['patient'] if v is None else (r['patient'],v))
    results.append(summarize(cohort,method,[p for p,v in available],[v for p,v in available],excluded))
add('GSE149614_scRNA','raw',expanded,lambda r:r['raw'])
add('GSE149614_scRNA','zone_adjusted_regression',expanded,lambda r:r['zone']['difference'])
for key in ['overlap','quadratic_overlap','entropy']:
    add('GSE149614_scRNA',key,bal['GSE149614'],lambda r,k=key:r.get(k,{}).get('difference'))
add('GSE149614_scRNA_CNV_candidates','zone_adjusted_regression',cnv,lambda r:r['regression']['difference'])
for key in ['overlap','entropy']:
    add('GSE149614_scRNA_CNV_candidates',key,cnv,lambda r,k=key:r['balancing'].get(k,{}).get('difference'))
add('GSE189175_snRNA','raw',sn,lambda r:r['selections']['core']['raw'])
add('GSE189175_snRNA','zone_adjusted_regression',sn,lambda r:r['selections']['core']['models']['zone']['difference'])
for key in ['overlap','quadratic_overlap','entropy']:
    add('GSE189175_snRNA',key,bal['GSE189175'],lambda r,k=key:r.get(k,{}).get('difference'))
meta={'target':'Equal-patient mean of fixed within-patient tumour minus adjacent signature-score contrasts; method-specific populations. Not a causal effect.','bootstrap':'Exact multinomial patient-resampling percentile distribution, all count vectors, no Monte Carlo error.','limits':['Small patient samples, especially n=3, limit interval calibration and generalisability.','Patient resampling treats patient contrasts as fixed: does not refit CNV classification, clustering or within-patient weights, and does not quantify their uncertainty.','Student-t intervals assume approximately normally distributed independent patient contrasts; reported as sensitivity, not evidence of normality.','No multiplicity-adjusted significance claims or equivalence conclusions.','Entropy summaries omit infeasible patients and target a different subset.','LOO is influence sensitivity, not prediction/cross-validation.','CNV labels remain RNA-derived candidates; adjacent tissue is not DNA-verified benign.'],'input_sha256':{f:hashlib.sha256((W/f).read_bytes()).hexdigest() for f in ['expanded_control_results.json','balance_results.json','cnv_signature_results.json','GSE189175_score_results.json']}}
(O/'hasta_duzeyi_belirsizlik.json').write_text(json.dumps({'methods':meta,'results':results},indent=2))
with (O/'hasta_duzeyi_disarida_birakma.tsv').open('w') as f:
    writer=csv.writer(f,delimiter='\t');writer.writerow(['cohort','method','excluded_patient','remaining_n','full_mean','loo_mean','change_from_full'])
    for r in results:
        for a in r['leave_one_patient_out']:writer.writerow([r['cohort'],r['method'],a['excluded_patient'],r['n_patients']-1,r['mean'],a['mean'],a['change_from_full']])
lines=['# Hasta düzeyindeki belirsizlik ve tek hasta çıkarma analizi','','Her hastaya eşit ağırlık verilmiştir. Fark, tümör eksi komşu doku yönündedir; ölçü birimi 47 genin ortalama log-normalize ifade puanıdır. Kohortlar havuzlanmamıştır.','','## Yöntem','','Hasta içi kayıtlı farklar sabit tutularak hastalar yerine koymalı yeniden örneklenmiştir. Bootstrap dağılımı multinom sayım vektörlerinin tamamı üzerinden tam olarak hesaplanmıştır; rastgele örnekleme hatası yoktur. Sekiz hasta için 6.435 farklı sayım vektörü vardır. Aralıklar yüzde 2,5 ve 97,5 ağırlıklı dağılım eşikleridir. Ek olarak hasta farklarının yaklaşık normal dağıldığı varsayımına dayanan Student-t aralıkları verilmiştir.','','Bu aralıklar CNV etiketlerinin, hücre anotasyonunun veya hasta içi model tahminlerinin tüm belirsizliğini kapsamaz. Küçük hasta sayısında nominal %95 kapsam garanti değildir. Üç hastalı kohortta bootstrap dağılımı özellikle kesiklidir; güven aralığı biçimi güçlü doğrulama olarak yorumlanamaz.','','## Özet','','| Kohort / yöntem | n | Ortalama | Hasta bootstrap %95 aralığı | Student-t %95 aralığı | Tek hasta çıkarılınca ortalama aralığı | İşaret değişimi |','|---|---:|---:|---|---|---|---|']
for r in results:
    fmt=lambda a:f'[{a[0]:.5f}, {a[1]:.5f}]'
    lines.append(f"| {r['cohort']} / {r['method']} | {r['n_patients']} | {r['mean']:.5f} | {fmt(r['exact_patient_bootstrap_percentile_95'])} | {fmt(r['student_t_95'])} | {fmt([r['loo_min'],r['loo_max']])} | {'Evet' if r['loo_sign_flip'] else 'Hayır'} |")
lines+=['','## Hasta seçimi ve yorum sınırları','']
for r in results:
    if r['excluded_unavailable_patients']: lines.append(f"- {r['cohort']} / {r['method']}: hesaplanabilir dengeleme hedefi bulunmadığı için dışarıda kalan hastalar: {', '.join(r['excluded_unavailable_patients'])}. Bu özet tüm hastaların ortalamasıyla aynı hedef değildir.")
lines+=['','Tek hasta çıkarma, tek bir hastanın ortalama üzerindeki etkisini sınar. Kalan ortalamanın işaretinin korunması bütün hastaların aynı yönde olduğu veya biyolojik ilişkinin doğrulandığı anlamına gelmez. Aralığın sıfırı içermesi etkisizlik ya da eşdeğerlik kanıtı değildir. Düşük etkin örneklem büyüklükleri ve doku seçimiyle ilişkili yanlılıklar sürer.','','Ham farklar bağlam sağlamak için eklenmiştir; düzeltilmiş modellerin yerine kullanılmamalıdır. Çoklu yöntemler keşifsel duyarlılık değerlendirmesidir; anlamlı sonuç seçmek için tarama yapılmamıştır.','','## Hesaplama kontrolü','','Tüm bootstrap olasılıkları 1 toplamı, tam dağılım varyansı Var(hasta farkı)/n eşitliği ve her kohortta benzersiz hasta kimlikleri denetlenmiştir. Leave-one-out sonuçları ayrı TSV dosyasında kayıtlıdır.']
(O/'hasta_duzeyi_belirsizlik_raporu.md').write_text('\n'.join(lines)+'\n')
print(json.dumps([{k:r[k] for k in ['cohort','method','n_patients','mean','exact_patient_bootstrap_percentile_95','loo_min','loo_max','loo_sign_flip']} for r in results],indent=2))
