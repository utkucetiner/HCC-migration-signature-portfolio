# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,json
from pathlib import Path
import numpy as np
z=np.load('work/hepatocyte_scoring_inputs.npz')
cells=z['cells'];samples=z['sample'];total=z['total']
sample_names=np.array(sorted(set(samples)))
masks=[np.flatnonzero(samples==s) for s in sample_names]
paired=[f'HCC{i:02}{t}' for i in range(3,11) for t in ('N','T')]
pidx=np.array([list(sample_names).index(s) for s in paired])
pairmask=np.isin(samples,paired)
names=[];means=[];detections=[]
with gzip.open('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rt') as f:
 header=f.readline().strip().split('\t');lookup={c:i for i,c in enumerate(header)};ix=np.array([lookup[c] for c in cells])
 for i,line in enumerate(f,1):
  g,s=line.rstrip('\r\n').split('\t',1);v=np.fromstring(s,sep='\t',dtype=np.int64)[ix]
  log=np.log1p(v/total*10000)
  names.append(g);means.append([log[m].mean() for m in masks]);detections.append(float((v[pairmask]>0).mean()))
  if i%5000==0:print('Summarized',i,flush=True)
means=np.asarray(means);detections=np.asarray(detections);names=np.asarray(names)
np.savez_compressed('work/gene_sample_summaries.npz',genes=names,samples=sample_names,means=means,detection=detections)
q=json.loads(Path('work/qc_audit.json').read_text());targets=list(q['verified_signature_mapping'].values())
idx=np.array([list(names).index(g) for g in targets])
avg=means[:,pidx].mean(axis=1)
# Two-dimensional quantile strata based only on expression and detection, without tissue labels.
def bins(v):
 edges=np.unique(np.quantile(v,np.linspace(0,1,6)))
 return np.searchsorted(edges[1:-1],v,side='right')
eb,db=bins(avg),bins(detections)
excluded=set(targets)|{'MKI67','CCNA2','CCNB2','STMN1'}
pool=np.array([i for i,g in enumerate(names) if g not in excluded and not g.startswith('MT-') and avg[i]>0])
rng=np.random.default_rng(20261010);random_sets=[];fallback=0
for b in range(1000):
 chosen=[]
 for t in idx:
  candidates=pool[(eb[pool]==eb[t])&(db[pool]==db[t])]
  candidates=np.setdiff1d(candidates,chosen)
  if len(candidates)==0:
   fallback+=1;candidates=np.setdiff1d(pool,chosen)
   distance=abs(eb[candidates]-eb[t])+abs(db[candidates]-db[t])
   candidates=candidates[distance==distance.min()]
  chosen.append(int(rng.choice(candidates)))
 random_sets.append(chosen)
sets=np.asarray(random_sets)
signature=means[idx][:,pidx].mean(axis=0)
control=means[sets][:,:,pidx].mean(axis=1)
d=signature[1::2]-signature[::2];rd=control[:,1::2]-control[:,::2]
med=np.median(rd,axis=0);lo,hi=np.quantile(rd,[.025,.975],axis=0)
percentile=np.mean(rd<=d,axis=0)*100
md=['# İfade düzeyi eşleştirilmiş rastgele gen kontrolü','','10 Ekim 2026','','47 genlik uyarlanmış imza, 1.000 rastgele 47 genlik kümeyle karşılaştırıldı. Rastgele kümeler eşleşmiş örneklerdeki ortalama normalize ifade ve hücrelerde saptanma oranının beşerli dilimlerini eşleştiriyor. Seçim tümör farklarına göre yapılmadı. İmza genleri, dört döngü belirteci ve MT- genleri kontrol havuzundan çıkarıldı; küme içinde genler tekrarlanmıyor. Sabit rastgelelik tohumu 20261010.','','Bu basit dilim eşleştirmesi bazal ifade özelliklerini yaklaşık eşleştirir; korelasyon yapısını, gen uzunluğunu veya bütün biyolojik özellikleri eşleştirmez. Kümeler birbirinden bağımsız biyolojik deney değildir. Yüzdelikler klinik anlamlılık veya düzeltilmiş p değeri değildir.','','| Hasta | İmza farkı | Kontrol medyan farkı | Kontrol %2.5 | Kontrol %97.5 | İmza kontrol yüzdeliği |','|---|---:|---:|---:|---:|---:|']
for j in range(8):
 md.append(f"| HCC{j+3:02} | {d[j]:+.4f} | {med[j]:+.4f} | {lo[j]:+.4f} | {hi[j]:+.4f} | {percentile[j]:.1f} |")
md+=['','Kontrol aralığı rastgele gen kümelerinin dağılımıdır; hasta etkisinin güven aralığı değildir. Pozitif ve negatif uçlar ayrı ayrı keşif amacıyla görülebilir; bunlardan sonuç sonrası eşik veya gen seçimi yapılmayacak.','','Hasta başına fark, önce hücre düzeyinde normalize log ifadenin örnek ortalaması, sonra T eksi N olarak hesaplandı. Hücrelerin örnekler arasında farklı derinliğe sahip olması hâlâ etkili olabilir. Malignite sınıflaması, kapsamlı döngü kontrolü ve bağımsız kohort doğrulaması tamamlanmadı.']
Path('outputs/rastgele-gen-kontrolu.md').write_text('\n'.join(md)+'\n')
Path('work/random_control_result.json').write_text(json.dumps({'seed':20261010,'n_sets':1000,'fallback_draws':fallback,'signature_difference':d.tolist(),'control_median':med.tolist(),'lo':lo.tolist(),'hi':hi.tolist(),'percentile':percentile.tolist()},indent=2))
print('Fallback draws:',fallback)
print('Percentiles:',percentile.tolist())
print('Signature differences:',d.tolist())
print('Random median:',med.tolist())
