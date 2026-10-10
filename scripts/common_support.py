# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import numpy as np,json
from pathlib import Path
z=np.load('work/hepatocyte_scoring_inputs.npz');c=np.load('work/expanded_control_scores.npz')
assert np.array_equal(z['cells'],c['cells'])
C=np.column_stack([np.log1p(z['total']),c['S'],c['G2M'],c['identity'],c['portal']-c['central']])
names=['log_depth','S','G2M','identity','portal_minus_central']
original=json.loads(Path('work/expanded_control_results.json').read_text())
rows=[]
for n in range(3,11):
 p=f'HCC{n:02}';a=np.flatnonzero(z['sample']==p+'N');b=np.flatnonzero(z['sample']==p+'T')
 bounds=[];keep=np.ones(len(z['cells']),bool)
 for j in range(5):
  an=np.quantile(C[a,j],[.05,.95]);bt=np.quantile(C[b,j],[.05,.95]);lo=max(an[0],bt[0]);hi=min(an[1],bt[1])
  bounds.append([float(lo),float(hi)])
  keep&=(C[:,j]>=lo)&(C[:,j]<=hi)
 aa=a[keep[a]];bb=b[keep[b]]
 row={'patient':p,'original_N':len(a),'original_T':len(b),'retained_N':len(aa),'retained_T':len(bb),'bounds':dict(zip(names,bounds)),'full_adjusted':next(x['zone']['difference'] for x in original['patients'] if x['patient']==p)}
 if min(len(aa),len(bb))>=20:
  ix=np.r_[aa,bb];wgt=np.r_[np.full(len(aa),1/len(aa)),np.full(len(bb),1/len(bb))];w=np.sqrt(wgt)
  centered=C[ix]-np.average(C[ix],axis=0,weights=wgt);sd=np.sqrt(np.average(centered**2,axis=0,weights=wgt))
  if np.all(sd>0):
   X=np.column_stack([np.ones(len(ix)),np.r_[np.zeros(len(aa)),np.ones(len(bb))],centered/sd])
   beta,res,rank,sv=np.linalg.lstsq(X*w[:,None],z['signature'][ix]*w,rcond=None)
   row.update(rank=int(rank),condition=float(sv[0]/sv[-1]),restricted_raw=float(z['signature'][bb].mean()-z['signature'][aa].mean()))
   if rank==7:row['restricted_adjusted']=float(beta[1])
 rows.append(row)
Path('work/common_support_results.json').write_text(json.dumps({'rule':'Intersection of within-tissue 5th–95th percentile ranges on all five covariates; minimum 20 cells per tissue for descriptive fit; fixed before restricted scores. This is marginal overlap, not joint support or matching.','patients':rows},indent=2))
md=['# Kovaryat örtüşmesi duyarlılığı','','## Amaç ve kural','','Genişletilmiş modelin örtüşmeyen hücrelere duyarlılığını incelemek için her hastanın iki dokusundaki beş kovaryatın ayrı ayrı %5–%95 aralıklarının kesişimi alındı. Bir hücre ancak beş aralıkta da kalıyorsa tutuldu. Bu kural sınırlı örneklemde keşif amaçlıdır; ortak çok boyutlu destek, eşleştirme veya kovaryat dengesi sağlamaz. Eşikler kısıtlı farklar hesaplanmadan önce sabitlendi.','','En az 20 hücre/doku kaldığında aynı eşit doku ağırlıklı doğrusal model tekrar kuruldu; bu eşik istatistiksel yeterlilik güvencesi değildir.','','| Hasta | Kalan N / ilk N | Kalan T / ilk T | Tam model | Kısıtlı ham | Kısıtlı model |','|---|---:|---:|---:|---:|---:|']
for r in rows:
 def f(k):return f"{r[k]:+.4f}" if k in r else 'hesaplanmadı'
 md.append(f"| {r['patient']} | {r['retained_N']}/{r['original_N']} | {r['retained_T']}/{r['original_T']} | {f('full_adjusted')} | {f('restricted_raw')} | {f('restricted_adjusted')} |")
md+=['','## Yorum sınırı','','Hücreleri bu özelliklere göre seçmek hedef popülasyonu değiştirir ve biyolojik sinyali de çıkarabilir. Bu nedenle kısıtlı sonuç tam modelin yerine geçmez. Az hücre kalan hastalar kanıt yokluğu veya fark yokluğu olarak yorumlanamaz. Hasta sayısı hâlâ sekizdir; hücre düzeyinde anlamlılık testi yapılmadı. Pozitif değer tümör dokusunda daha yüksek onarım imzası puanı anlamına gelir.','','[Önceki genişletilmiş kontrol raporu](genisletilmis-kontroller-raporu.md)']
Path('outputs/kovaryat-ortusmesi-duyarliligi.md').write_text('\n'.join(md)+'\n')
print(json.dumps(rows,indent=2))
