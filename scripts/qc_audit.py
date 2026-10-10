# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,csv,json
from pathlib import Path
import numpy as np
with gzip.open('work/data/GSE149614_HCC.metadata.updated.txt.gz','rt',newline='') as f:
 meta={r['Cell']:r for r in csv.DictReader(f,delimiter='\t')}
prev=json.loads(Path('work/matrix_audit.json').read_text())
mapping={g:g for g in prev['present_signature_genes']}
mapping.update({'CYTOR':'LINC00152','PRAG1':'SGK223'})
assert len(mapping)==47 and len(set(mapping.values()))==47
targets=set(mapping.values())
with gzip.open('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rt') as f:
 cells=f.readline().strip().split('\t')
 n=len(cells);total=np.zeros(n,dtype=np.int64);det=np.zeros(n,dtype=np.int32);mt=np.zeros(n,dtype=np.int64);sig={};mtgenes=[]
 for i,line in enumerate(f,1):
  gene,rest=line.rstrip('\r\n').split('\t',1)
  v=np.fromstring(rest,sep='\t',dtype=np.int64)
  assert len(v)==n and np.all(v>=0)
  total+=v;det+=v>0
  if gene.startswith('MT-'):mt+=v;mtgenes.append(gene)
  if gene in targets:sig[gene]=v>0
  if i%5000==0:print('Checked',i,'gene rows',flush=True)
assert set(sig)==targets
assert np.all(total>0)
pct=100*mt/total
samples=sorted({r['sample'] for r in meta.values()})
result=[]
for sample in samples:
 mask=np.array([meta[c]['sample']==sample and meta[c]['celltype']=='Hepatocyte' for c in cells])
 count=int(mask.sum());v=total[mask];d=det[mask];p=pct[mask]
 cov={g:int(sig[m][mask].sum()) for g,m in mapping.items()}
 result.append({'sample':sample,'hepatocytes':count,'median_counts':float(np.median(v)),'median_genes':float(np.median(d)),'median_mt_pct':float(np.median(p)),'p95_mt_pct':float(np.percentile(p,95)),'cells_lt200_genes':int((d<200).sum()),'cells_gt20pct_mt':int((p>20).sum()),'signature_genes_detected_any':sum(z>0 for z in cov.values()),'signature_genes_detected_at_least10pct':sum(z/count>=.1 for z in cov.values()),'signature_detection':cov})
out={'verified_signature_mapping':mapping,'n_cells':n,'n_genes':i,'mt_genes':mtgenes,'samples':result,'scope':'Descriptive QC; no filtering or tumour association performed.'}
Path('work/qc_audit.json').write_text(json.dumps(out,indent=2))
md=['# HCC hepatosit kalite ve gen saptanma raporu','','10 Ekim 2026','','71.915 hücre ve 25.712 gen satırındaki tüm sayımlar negatif olmayan tam sayı olarak okundu. Her hücrede toplam sayım sıfırdan büyük. Analiz yalnızca yayımlanmış Hepatocyte etiketli hücreleri örnek bazında özetler; maligniteyi doğrulamaz.','','47 genlik geçici liste kullanıldı: 45 doğrudan eşleşme, CYTOR için LINC00152 ve PRAG1 için SGK223. İki aday eşleşme eklenmedi.','','| Örnek | Hepatosit | Medyan sayım | Medyan gen | Medyan MT % | İmza saptanan gen /47 | En az %10 hücrede /47 |','|---|---:|---:|---:|---:|---:|---:|']
for r in result:
 md.append(f"| {r['sample']} | {r['hepatocytes']} | {r['median_counts']:.0f} | {r['median_genes']:.0f} | {r['median_mt_pct']:.2f} | {r['signature_genes_detected_any']} | {r['signature_genes_detected_at_least10pct']} |")
md+=['','## Yorum sınırları','','MT yüzdesi matristeki MT- ile başlayan genlerin toplam sayım içindeki oranıdır. Dosya önceden işlenmiş olduğu için bu ölçümler ham droplet kalite kontrolünün yerini tutmaz. Doublet veya ambient RNA değerlendirmesi tamamlanmadı.','','İmza saptanması en az bir hücrede sıfırdan büyük sayım bulunmasıdır; skor veya etki değildir. %10 sütunu yalnızca saptanma kararlılığına dair betimleyici bir özettir; gen eleme kuralı değildir. Hücreler veya hastalar bu rapora dayanarak elenmedi.','','Az hepatosit içeren örnekler için hasta düzeyindeki belirsizlik ve yeniden örneklemeyle kararlılık incelenmeli. Gen sayıları veya MT yüzdesi için sonuçlara göre eşik seçilmeyecek. Önce kalite kuralları, malignite kapsamı ve puanlama yöntemi belirlenecek.','','[Veri kaynağı](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)','','Henüz tümör etkisi, onarım skoru veya p değeri hesaplanmadı.']
Path('outputs/hepatosit-kalite-raporu.md').write_text('\n'.join(md)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k!='signature_detection'} for r in result],indent=2),flush=True)
