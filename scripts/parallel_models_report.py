# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import json
from pathlib import Path
sc=json.loads(Path('work/expanded_control_results.json').read_text());sn=json.loads(Path('work/GSE189175_score_results.json').read_text());bal=json.loads(Path('work/balance_results.json').read_text());md=['# scRNA-seq ve snRNA-seq: ayrı model sonuçları','','Bu tablolar iki ayrı kohortun hasta içi betimleyici sonuçlarıdır. Hücreler, katsayılar veya hastalar ortak bir etki tahmininde havuzlanmadı. Skor, ortak 47 genin log1p(10000×sayım/toplam sayım) ortalamasıdır; aynı sayısal tanım ölçüm platformlarını eşdeğer yapmaz.','','## GSE149614: scRNA-seq','','| Hasta | Derinlik+S+G2M | +Kimlik | +Zonasyon | Örtüşme ağırlıklı fark | Karesel örtüşme farkı |','|---|---:|---:|---:|---:|---:|']
print(sc['patients'][0].keys())
for row in sc['patients']:
 b=next(x for x in bal['GSE149614'] if x['patient']==row['patient']);md.append(f"| {row['patient']} | {row['cycle']['difference']:+.4f} | {row['identity']['difference']:+.4f} | {row['zone']['difference']:+.4f} | {b['overlap']['difference']:+.4f} | {b['quadratic_overlap']['difference']:+.4f} |")
md+=['','## GSE189175: snRNA-seq','','| Hasta | Derinlik+S+G2M | +Kimlik | +Zonasyon | Örtüşme ağırlıklı fark | Karesel örtüşme farkı |','|---|---:|---:|---:|---:|---:|']
for row in sn['patients']:
 a=row['selections']['core']['models'];b=next(x for x in bal['GSE189175'] if x['patient']==row['patient']);md.append(f"| {row['patient']} | {a['cycle']['difference']:+.4f} | {a['identity']['difference']:+.4f} | {a['zone']['difference']:+.4f} | {b['overlap']['difference']:+.4f} | {b['quadratic_overlap']['difference']:+.4f} |")
md+=['','Lineer model sütunları tümör doku göstergesinin katsayısıdır. Örtüşme sütunları yeniden ağırlıklandırılmış ortalama skor farklarıdır; lineer regresyon katsayıları değildir. Hedef hücre popülasyonları farklı olabilir. Hepsi T−N yönündedir.','','GSE189175’in soy etiketleri yayımlanmış barkod anotasyonları değildir; çoklu belirteçle yeniden yapılan keşif amaçlı anotasyona dayanır. CNV destekli aday malignite altkümesinin GSE149614 sonucu ayrı raporlanır; bu tablolar tüm hepatosit soyundaki doku karşılaştırmalarıdır.','','Hasta 2’nin snRNA-seq sonucundaki değişimin intron/ekzon oranıyla veya SoupX düzeltmesiyle ilişkisi ölçülemedi; gerekli ham/katmanlı dosyalar yayımlanmış bu veri paketinde yok. Bir teknik açıklama ileri sürülmedi.','','[Kovaryat dengeleme ayrıntıları](kovaryat-dengeleme-v02.md) · [Teknik dosya sınırları](sc_sn_teknik_sinirlar.md)']
Path('outputs/iki_kohort_ayri_modeller.md').write_text('\n'.join(md)+'\n')
