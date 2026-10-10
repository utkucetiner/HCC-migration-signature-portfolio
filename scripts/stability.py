# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import numpy as np,json
from pathlib import Path
z=np.load('work/hepatocyte_scoring_inputs.npz')
s=z['sample'];y=z['signature'];depth=np.log1p(z['total']);cycle=z['cycling']
rng=np.random.default_rng(20261010);out=[]
for i in range(3,11):
 p=f'HCC{i:02}';a=np.flatnonzero(s==p+'N');b=np.flatnonzero(s==p+'T')
 ix=np.r_[a,b];t=np.r_[np.zeros(len(a)),np.ones(len(b))]
 # Equal sample weighting prevents the sample with more cells dominating the fit.
 weights=np.r_[np.full(len(a),1/len(a)),np.full(len(b),1/len(b))]
 def fit(columns):
  X=np.column_stack([np.ones(len(ix)),t]+columns)
  w=np.sqrt(weights)
  return float(np.linalg.lstsq(X*w[:,None],y[ix]*w,rcond=None)[0][1])
 raw=float(y[b].mean()-y[a].mean())
 dep=fit([depth[ix]-np.average(depth[ix],weights=weights)])
 both=fit([depth[ix]-np.average(depth[ix],weights=weights),cycle[ix]-np.average(cycle[ix],weights=weights)])
 draws=np.array([rng.choice(y[b],len(b),replace=True).mean()-rng.choice(y[a],len(a),replace=True).mean() for _ in range(2000)])
 lo,hi=np.quantile(draws,[.025,.975])
 out.append({'patient':p,'raw':raw,'depth_adjusted_descriptive':dep,'depth_cycle_adjusted_descriptive':both,'cell_bootstrap_lo':float(lo),'cell_bootstrap_hi':float(hi),'n_nontumor':len(a),'n_tumor':len(b)})
d=np.array([r['raw'] for r in out])
patient_draws=np.mean(rng.choice(d,(10000,len(d)),replace=True),axis=1)
pinterval=np.quantile(patient_draws,[.025,.975])
md=['# Derinlik ve yeniden örnekleme duyarlılığı','','10 Ekim 2026','','Bu çalışma sekiz hastadaki Hepatocyte etiketli hücrelerin betimleyici kontrolüdür. Maligniteyi, tümör mekanizmasını veya çoğalmadan bağımsız rejenerasyonu doğrulamaz.','','Her hasta içinde iki dokuya eşit toplam ağırlık verilerek imza puanı doku göstergesi ve log1p toplam sayımla doğrusal olarak modellendi. İkinci betimleyici model dört genlik döngü puanını da içeriyor. Katsayılar hücre düzeyindeki ilişkileri özetler; hasta düzeyinde nedensel veya doğrulanmış düzeltilmiş etki değildir. P değerleri hesaplanmadı.','','| Hasta | İlk fark | Derinlik kontrollü | Derinlik ve 4 gen döngü kontrollü | Hücre yeniden örnekleme alt | Üst |','|---|---:|---:|---:|---:|---:|']
for r in out:md.append(f"| {r['patient']} | {r['raw']:+.4f} | {r['depth_adjusted_descriptive']:+.4f} | {r['depth_cycle_adjusted_descriptive']:+.4f} | {r['cell_bootstrap_lo']:+.4f} | {r['cell_bootstrap_hi']:+.4f} |")
md+=['','## Belirsizliğin anlamı','','Her dokunun gözlenen hücreleri kendi içinde 2.000 kez yeniden örneklendi. Tablo aralıkları mevcut hücre örnekleminin kararlılığını özetler; bağımsız doku tekrarları, yeni hastalar, hazırlama etkileri veya biyolojik popülasyon belirsizliğini kapsamaz. Hücre bağımlılığı bu aralıkları dar gösterebilir.','','Sekiz hastanın ham farklarına eşit ağırlık verilerek 10.000 hasta yeniden örneklemesi yapıldı. Küçük hasta sayısı ve tek kohort nedeniyle yalnızca keşif amaçlıdır.',f"",'Hasta ortalama ham farkı: '+f'{d.mean():+.4f}; hasta bootstrap yüzdelik aralığı: [{pinterval[0]:+.4f}, {pinterval[1]:+.4f}].','','Doğrusal derinlik kontrolü toplam sayımla biyolojik süreçlerin ilişkisini de çıkarabilir. Dört gen bütün hücre döngüsünü temsil etmez. Model, kovaryat aralıkları örtüşmediğinde ekstrapolasyon yapabilir; bu nedenle sonuç nihai çıkarım değildir.','','## Sonraki gereklilik','','Kovaryat örtüşmesi incelenmeli; daha geniş hücre döngüsü ve alternatif imza puanlaması tanımlanmalı. Dokuya özgü malignite anotasyonu ve bağımsız kohort doğrulaması tamamlanmalı. Araştırmanın tek yönlü artış hipotezi sonuçlara uydurulmayacak.']
Path('outputs/derinlik-ve-kararlilik-raporu.md').write_text('\n'.join(md)+'\n')
Path('work/stability_results.json').write_text(json.dumps({'patients':out,'patient_mean':float(d.mean()),'patient_bootstrap_interval':pinterval.tolist(),'seed':20261010},indent=2))
print(json.dumps(out,indent=2));print('Patient mean and interval',d.mean(),pinterval)
