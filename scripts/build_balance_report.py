# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,numpy as np
r=json.loads(Path('work/balance_results.json').read_text());md=['# Hasta içinde kovaryat dengeleme: ayrı kohort sonuçları','','10 Ekim 2026','','Önceki beş ayrı %5–%95 aralığının kesişimi filtresi ana dengeleme yöntemi olarak kullanılmadı. Hasta içinde tümör dokusunda bulunma olasılığı, beş kovaryatlı cezasız lojistik modelle tahmin edildi: log1p toplam sayım, S, G2M, hepatosit kimliği ve portal−santral zonasyon. Onarım puanı eğilim skoruna alınmadı.','','Tümör hücrelerinin ağırlığı 1−e(X), komşu doku hücrelerininki e(X) olarak tanımlandı; her dokunun ağırlıkları kendi içinde toplamı bire gelecek şekilde normalize edildi. Hesaplanan fark, bu ağırlıklarla ortalama 47 genlik imza puanının T−N farkıdır. Bütün hücreler hesaplamaya girdi. Bu yöntem birebir PSM değildir; eğilim skoruna dayalı örtüşme ağırlıklandırmasıdır. Hedef, bütün dokunun ortalama etkisinden kovaryatları örtüşen hücrelerin betimleyici farkına değişir.','','Entropi duyarlılığında her dokunun beş kovaryat ortalaması, iki dokunun ağırlıksız ortalamalarının orta noktasına taşındı. Her grup için hedefin konveks zarf içinde olup olmadığı doğrusal programlama ile sınandı; mümkünse entropi ikilisi optimize edildi. Başarısız dengeleme için etki raporlanmadı.','','SMD paydası iki dokunun ağırlıksız varyanslarının ortalamasının kareköküdür. ESS=(Σw)²/Σw² etkin hücre sayısıdır; bağımsız hasta sayısı değildir. Küçük SMD yalnızca modele alınan birinci momentlerin dengesi anlamına gelir. Dağılımların, etkileşimlerin veya ölçülmemiş kovaryatların dengelendiğini göstermez.']
for cohort,rows in r.items():
 md+=['',f'## {cohort}: '+('scRNA-seq, sekiz hasta' if cohort=='GSE149614' else 'snRNA-seq, üç hasta'),'','| Hasta | İlk N/T | Ham T−N | Örtüşme T−N | ESS N/T | En büyük |SMD| | Entropi T−N |','|---|---:|---:|---:|---:|---:|---:|']
 # Escape pipe in column header for Markdown.
 md[-2]=md[-2].replace('|SMD|','SMD (mutlak)')
 for x in rows:
  a=x['overlap'];eb=x.get('entropy');ef=f"{eb['difference']:+.4f}" if eb else 'uygun hedef yok'
  md.append(f"| {x['patient']} | {x['n_N']}/{x['n_T']} | {x['raw']['difference']:+.4f} | {a['difference']:+.4f} | {a['ESS_N']:.1f}/{a['ESS_T']:.1f} | {max(abs(v) for v in a['SMD'].values()):.2e} | {ef} |")
 md+=['',f"Örtüşme ağırlıklı farkların eşit hasta ağırlıklı ortalaması {np.mean([x['overlap']['difference'] for x in rows]):+.4f}. Bu ortalama ayrı kohort içindedir; güven aralığı veya hipotez testi üretilmedi."]
md+=['','## Yorum sınırları','','HCC07 için orta nokta entropi hedefi uygulanabilir değildir. Örtüşme ağırlıkları yine hesaplanabilir. Bazı hastalarda hücrelerin büyük bölümü çok düşük ağırlık alır; sekiz hastanın hesaplamada kalması güçlü bilgi sağlandığı anlamına gelmez. Özellikle normal doku ESS değerleri HCC05 için yaklaşık 27, HCC07 için 30, HCC08 için 42 düzeyindedir.','','GSE189175 hasta 2 için ilk hücre sayıları 4053/4666 iken ESS yaklaşık 252/325 olur. Önceki çok katı filtrede 3/36 hücre kalması yönteme duyarlıdır, fakat ağırlıklandırma da ciddi ayrışmaya işaret eder. Yeni +0.0381 değeri teknik veya biyolojik kökeni belirlenmiş bir etki değildir.','','Kovaryatlar tümör biyolojisinin sonuçları da olabilir. Bu yüzden dengeleme nedensel tümör etkisini tanımlamaz, gerçek biyolojik programın bir kısmını çıkarabilir. Hücreler hasta içinde bağımlıdır; hücre sayısıyla p değeri üretilemez. Malignite etiketleri bu dengeleme aşamasında henüz kullanılmadı; burada doku karşılaştırması vardır.','','Kohortlar birleştirilmedi. scRNA-seq ve snRNA-seq katsayıları ve ağırlıklı farkları ayrı tutuldu.','','Yöntem kaynağı: [Li, Morgan ve Zaslavsky, Balancing Covariates via Propensity Score Weighting](https://doi.org/10.1080/01621459.2016.1260466).']
Path('outputs/kovaryat-dengeleme-v02.md').write_text('\n'.join(md)+'\n')
Path('outputs/Sharma_kohortu_eleme_eki.md').write_text('''# Ek veri: Sharma kohortunun elenmesi

Sharma/Cell 2020 HCC veri kaynağı (Mendeley 6wmzcskt6k/1), eşleşmiş hepatosit analizi kriterlerini karşılamadığı için elenmiştir. Ana analiz ve doğrulama sonuçlarına dahil edilmemiştir.

73589 hücrelik dosyada `patientno` gerçek hasta kimliği, `patient_id` örnek/batch sayacıdır. Komşu dokularda aday hepatosit kümeleri kullanıldığında hasta başına yalnızca 0–7 aday hücre bulundu; en az 20 hücre/doku gerektiren eşleşmiş karşılaştırma yapılamadı. Bu aday kümeler kaynak yazarların doğrulanmış barkod etiketleri değildir. Dolayısıyla eleme, karaciğerde hepatosit bulunmadığı anlamına gelmez; indirilen dosyanın bu analiz için uygun olmadığı anlamına gelir.

Onarım sonucu veya eşleşmiş model bu kohortta hesaplanmadı. Daha fazla anotasyon çalışması ana proje kapsamında sürdürülmeyecek. GSE175111 erişim kodunun bu indirilen dosyayla birebir ilişkisi bu aşamada doğrulanmadığından, veri kimliği DOI ve Mendeley arşiviyle tutuldu.

[Makale](https://doi.org/10.1016/j.cell.2020.08.040) · [İndirilen veri arşivi](https://data.mendeley.com/datasets/6wmzcskt6k/1)
''')
print('reports saved')
