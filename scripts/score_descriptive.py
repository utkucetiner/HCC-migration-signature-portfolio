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
from types import SimpleNamespace
def ranks(x):
 order=np.argsort(x,kind='stable');v=x[order]
 starts=np.r_[0,np.flatnonzero(v[1:]!=v[:-1])+1];ends=np.r_[starts[1:],len(x)]
 ranked=np.empty(len(x),dtype=float)
 for a,b in zip(starts,ends):ranked[order[a:b]]=(a+b-1)/2
 return ranked
def spearmanr(x,y):
 a,b=ranks(x),ranks(y)
 return SimpleNamespace(statistic=np.corrcoef(a,b)[0,1])
qc=json.loads(Path('work/qc_audit.json').read_text())
mapping=qc['verified_signature_mapping'];cycle=['MKI67','CCNA2','CCNB2','STMN1']
targets=set(mapping.values())|set(cycle)
with gzip.open('work/data/GSE149614_HCC.metadata.updated.txt.gz','rt',newline='') as f:meta={r['Cell']:r for r in csv.DictReader(f,delimiter='\t')}
with gzip.open('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rt') as f:
 cells=f.readline().strip().split('\t')
 indices=np.array([i for i,c in enumerate(cells) if meta[c]['celltype']=='Hepatocyte'])
 selected=[cells[i] for i in indices]
 total=np.zeros(len(indices),dtype=np.int64);det=np.zeros(len(indices),dtype=np.int32);data={}
 for i,line in enumerate(f,1):
  g,s=line.rstrip('\r\n').split('\t',1)
  v=np.fromstring(s,sep='\t',dtype=np.int64)[indices]
  total+=v;det+=v>0
  if g in targets:data[g]=v
  if i%5000==0:print('Processed',i,flush=True)
assert targets<=set(data)
normalized={g:np.log1p(v/total*10000) for g,v in data.items()}
score=np.mean([normalized[g] for g in mapping.values()],axis=0)
cycling=np.mean([normalized[g] for g in cycle],axis=0)
samples=np.array([meta[c]['sample'] for c in selected])
np.savez_compressed('work/hepatocyte_scoring_inputs.npz',cells=np.array(selected),sample=samples,total=total,detected=det,signature=score,cycling=cycling)
summary={}
for s in sorted(set(samples)):
 mask=samples==s
 summary[s]={'cells':int(mask.sum()),'mean_signature':float(score[mask].mean()),'median_signature':float(np.median(score[mask])),'mean_cycling':float(cycling[mask].mean()),'rho_signature_counts':float(spearmanr(score[mask],total[mask]).statistic),'rho_signature_detected':float(spearmanr(score[mask],det[mask]).statistic),'rho_signature_cycling':float(spearmanr(score[mask],cycling[mask]).statistic)}
pairs=[]
for i in range(3,11):
 p=f'HCC{i:02}';t=summary[p+'T'];n=summary[p+'N']
 pairs.append({'patient':p,'tumor_mean':t['mean_signature'],'nontumor_mean':n['mean_signature'],'difference':t['mean_signature']-n['mean_signature'],'cycle_difference':t['mean_cycling']-n['mean_cycling']})
out={'score':'Mean of log1p(count/total*10000) for 47 matched genes per hepatocyte; then mean within sample','cycle':'Same normalization, four published cycling markers','samples':summary,'pairs':pairs,'inferential_test_performed':False}
Path('work/descriptive_scores.json').write_text(json.dumps(out,indent=2))
md=['# İlk betimleyici hepatosit imza puanları','','10 Ekim 2026','','Bu çıktı keşifseldir. Hepatocyte etiketli hücrelerin doku kökenini karşılaştırır; maligniteyi doğrulamaz. Puan, yayımlanmış 50 genlik listenin 47 eşleşmiş geniyle uyarlanmıştır. Eksik üç gen sıfır sayılmadı.','','Hücre toplam sayımı 10.000 ölçeğine normalize edildi, log1p dönüşümü yapıldı ve 47 gen ortalaması alındı. Ardından hücre puanları örnek içinde ortalandı. Bu basit puan gen ifadesi ve saptanma özelliklerinden etkilenebilir; özgül rejenerasyon ölçüsü olarak doğrulanmış değildir.','','| Hasta | Tümör ortalaması | Tümör dışı ortalaması | Eşleşmiş fark | Hücre döngüsü farkı |','|---|---:|---:|---:|---:|']
for r in pairs:md.append(f"| {r['patient']} | {r['tumor_mean']:.4f} | {r['nontumor_mean']:.4f} | {r['difference']:+.4f} | {r['cycle_difference']:+.4f} |")
md+=['','## Hücre içi karıştırıcı ilişkiler','','Aşağıdaki Spearman korelasyonları her örneğin hücreleri içinde betimleyici olarak hesaplandı. Hücre sayısını hasta örneklemi gibi kullanan p değerleri hesaplanmadı. Döngü puanı MKI67, CCNA2, CCNB2, STMN1 ortalamasıdır; kapsamlı bir hücre döngüsü düzeltmesi değildir.','','| Örnek | İmza ve toplam sayım rho | İmza ve saptanan gen rho | İmza ve döngü rho |','|---|---:|---:|---:|']
for s,r in summary.items():md.append(f"| {s} | {r['rho_signature_counts']:.3f} | {r['rho_signature_detected']:.3f} | {r['rho_signature_cycling']:.3f} |")
md+=['','## Yorum sınırı','','Bir puan farkı, çoğalmadan bağımsız onarım, malign hücre mekanizması veya nedensellik göstermez. Sayım derinliği, hepatosit durumu ve düşük hücre sayıları etkili olabilir. Hücre döngüsü listesi genişletilmeli; gen ifadesine göre eşleştirilmiş rastgele kümeler, alternatif puanlama ve hasta düzeyinde kararlılık analizleri yapılmalı. İstatistiksel bağımsız birim hastadır.','','Hiçbir hasta veya hücre bu puanlara göre elenmedi. Hasta düzeyinde p değeri, doğrulama kohortu sonucu veya düzeltilmiş tümör etkisi henüz yok.','','[Veri kaynağı](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)','','[İmza ve dört döngü belirtecinin kaynağı](https://www.nature.com/articles/s41586-024-07376-2)']
Path('outputs/ilk-betimleyici-puanlar.md').write_text('\n'.join(md)+'\n')
print(json.dumps(pairs,indent=2))
