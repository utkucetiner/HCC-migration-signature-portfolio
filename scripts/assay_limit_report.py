# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import gzip
series=Path('work/data/GSE189175_series.txt').read_text();features=[]
with gzip.open('work/data/GSE189175_counts_features.tsv.gz','rt') as f:
 for i in range(3):features.append(next(f).rstrip())
md='''# Ayrı kohortlar ve teknik bileşenlerin ölçülebilirliği

GSE149614 tek-hücre RNA-seq, GSE189175 tek-çekirdek RNA-seq olarak ayrı analiz edildi. Ortak bir hücre matrisi veya ortak bir model katsayısı oluşturulmadı. Katsayıların karşılaştırılması keşif amaçlıdır; ölçüm platformu ve hasta etiyolojisi birbirinden ayrıştırılamaz.

GSE189175 GEO seri kaydı ham verilerin insan katılımcı gizliliği nedeniyle sağlanmadığını açıkça bildiriyor: “raw data is not provided due to human subject privacy concerns”. Yayınlanan dört ek dosya gen sayım matrisi, gen özellikleri, barkodlar ve metadatadır. Gen matrisi intron/ekzon ayrımı içermiyor. Bu yüzden intronik/ekzonik oran ve hasta 2 katsayısının bu oranla ilişkisi hesaplanamadı. Toplam UMI veya saptanan gen sayısı bu oranların yerine geçirilemez.

SoupX için mevcut tutulan hücre/çekirdek matrisi yeterli değildir. Standart ortam RNA profilini tahmin etmek için boş damlacıkları da içeren filtresiz sayım matrisi gerekir; eldeki GSE189175 matrisi DIEM sonrası tutulmuş çekirdekleri içerir. Ham diziler, barkodlu hizalama dosyası ve filtresiz damlacık matrisi bu yerel dosyalar arasında yok. SoupX çalıştırılmış gibi bir çıktı oluşturulmadı.

İntron/ekzon oranını sınamak için barkodla ilişkilendirilebilen BAM veya eşdeğer intronik/ekzonik sayım katmanları gerekir. Ambient RNA düzeltmesi için örnek bazında filtresiz ve filtreli droplet matrisleri gerekir. Bunlar olmadan hasta 2’deki yön ve büyüklük değişimi ambient RNA’ya veya çekirdek yakalama oranına bağlanamaz. Mevcut derinlik ve çoklu soy belirteci kontrolleri yalnızca dolaylı duyarlılık incelemeleridir.

[Veri kaydı](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189175) · [Kaynak makale](https://doi.org/10.1186/s13073-022-01055-5) · [SoupX resmî uygulaması](https://github.com/constantAmateur/SoupX)
'''
assert 'raw data is not provided' in series
Path('outputs/sc_sn_teknik_sinirlar.md').write_text(md)
