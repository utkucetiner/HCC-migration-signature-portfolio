# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import numpy as np,json,gzip
from pathlib import Path
z=np.load('work/hepatocyte_scoring_inputs.npz')
sets={'S':Path('work/data/s_genes.txt').read_text().splitlines(),'G2M':Path('work/data/g2m_genes.txt').read_text().splitlines(),'identity':'ALB TTR TF HP HNF4A CYP2A6'.split(),'central':'CYP2E1 CYP3A4 ADH1B CXCL2 ADH1A RBP4 ALAS1 ADH4 ADH1C APOH AOX1'.split(),'portal':'CYP2A6 CRP SDS CYP2A7 SLPI FGA A2M HAL TF LBP ACTG1 CP KRT18 IGFBP2'.split()}
available=set(Path('work/matrix_gene_ids.txt').read_text().splitlines())
coverage={k:{'requested':len(v),'present':[g for g in v if g in available],'missing':[g for g in v if g not in available]} for k,v in sets.items()}
mapping=json.loads(Path('work/qc_audit.json').read_text())['verified_signature_mapping']
plan={'sets':sets,'coverage':coverage,'models':['raw','depth + S + G2M','depth + S + G2M + identity','depth + S + G2M + identity + portal-minus-central','depth + S + G2M + identity + portal + central'],'signature_overlap':{k:sorted(set(v)&set(mapping.values())) for k,v in sets.items()},'source_cycle':'https://github.com/satijalab/seurat/blob/master/data/cc.genes.rda','source_atlas':'https://www.nature.com/articles/s41586-024-07376-2','scoring':'mean log1p(count/total*10000), not Seurat module score','note':'Exploratory plan fixed before expanded outcome computation; original gene symbols, no unverified alias replacement.'}
Path('work/expanded_control_plan.json').write_text(json.dumps(plan,indent=2))
print(json.dumps(coverage),flush=True)
targets=set().union(*(set(v['present']) for v in coverage.values()))
data={}
with gzip.open('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz','rt') as f:
 header=f.readline().strip().split('\t'); lookup={c:i for i,c in enumerate(header)};idx=np.array([lookup[c] for c in z['cells']])
 for line in f:
  gene,values=line.rstrip('\r\n').split('\t',1)
  if gene in targets:data[gene]=np.fromstring(values,sep='\t',dtype=np.int64)[idx]
assert set(data)==targets
scores={k:np.mean([np.log1p(data[g]/z['total']*10000) for g in v['present']],axis=0) for k,v in coverage.items()}
np.savez_compressed('work/expanded_control_scores.npz',**scores,cells=z['cells'])
out=[]
for n in range(3,11):
 p=f'HCC{n:02}';a=np.flatnonzero(z['sample']==p+'N');b=np.flatnonzero(z['sample']==p+'T');ix=np.r_[a,b];t=np.r_[np.zeros(len(a)),np.ones(len(b))]
 weights=np.r_[np.full(len(a),1/len(a)),np.full(len(b),1/len(b))];w=np.sqrt(weights)
 controls=[np.log1p(z['total'][ix]),scores['S'][ix],scores['G2M'][ix],scores['identity'][ix],(scores['portal']-scores['central'])[ix]]
 def fit(cols):
  standardized=[]
  for c in cols:
   centered=c-np.average(c,weights=weights);sd=np.sqrt(np.average(centered**2,weights=weights));assert sd>0;standardized.append(centered/sd)
  X=np.column_stack([np.ones(len(ix)),t]+standardized);Xw=X*w[:,None]
  beta,res,rank,singular=np.linalg.lstsq(Xw,z['signature'][ix]*w,rcond=None)
  return {'difference':float(beta[1]),'rank':int(rank),'columns':X.shape[1],'condition':float(singular[0]/singular[-1])}
 row={'patient':p,'n_N':len(a),'n_T':len(b),'raw':float(z['signature'][b].mean()-z['signature'][a].mean()),'cycle':fit(controls[:3]),'identity':fit(controls[:4]),'zone':fit(controls),'separate_zones':fit(controls[:4]+[scores['portal'][ix],scores['central'][ix]])}
 row['control_mean_differences']={k:float(v[b].mean()-v[a].mean()) for k,v in scores.items()};out.append(row)
result={'plan':plan,'patients':out}
Path('work/expanded_control_results.json').write_text(json.dumps(result,indent=2))
md=['# Genişletilmiş kontrol analizi','','10 Ekim 2026','','GSE149614, sekiz eşleşmiş hasta; yalnızca Hepatocyte etiketli hücreler. Pozitif değer tümör dokusundaki 47 genlik onarım imzasının daha yüksek olduğunu gösterir.','','## Kaynak ve kapsam','','Hücre döngüsü: Seurat resmî cc.genes listesi (43 S, 54 G2/M). Hepatosit kimliği ve zonasyon: Matchett ve arkadaşlarının insan karaciğer atlası, ek Tablo 2. Listeler hesaplamadan önce kaydedildi; sonuçlara göre gen seçilmedi.']
for k,v in coverage.items():md.append(f"- {k}: {len(v['present'])}/{v['requested']} gen mevcut; eksik: {', '.join(v['missing']) or 'yok'}.")
md+=['','Tüm gen listeleri:']
for k,v in sets.items():md.append(f"- {k}: "+', '.join(v)+'.')
md+=['','## Hasta bazında sonuç','','| Hasta | Ham | Derinlik + S + G2/M | + Kimlik | + Zonasyon farkı | İki zon ayrı |','|---|---:|---:|---:|---:|---:|']
for r in out:md.append(f"| {r['patient']} | {r['raw']:+.4f} | {r['cycle']['difference']:+.4f} | {r['identity']['difference']:+.4f} | {r['zone']['difference']:+.4f} | {r['separate_zones']['difference']:+.4f} |")
md+=['','## Yöntem ve yorum sınırları','','Her gen için log1p(sayım / hücre toplam sayımı × 10.000) hesaplandı; her liste puanı mevcut genlerin eşit ağırlıklı ortalamasıdır. Bu işlem Seurat CellCycleScoring veya AddModuleScore değildir. S ve G2/M ayrı sürekli kovaryatlardır; hücre evresi atanmadı. Her hastanın iki dokusuna eşit toplam ağırlık veren doğrusal model kuruldu. Kovaryatlar ağırlıklı olarak standartlaştırıldı; tablo doku katsayısını gösterir.','','Kimlik listesindeki TF ve CYP2A6 portal listesinde de vardır. Ana zonasyon kontrolü portal eksi merkezi puandır; iki zonun ayrı eklendiği model ek duyarlılık kontrolüdür. Tüm modellerde tam matris sırası ve koşul sayısı kontrol edildi. Bu, biyolojik örtüşmeyi veya ekstrapolasyonu ortadan kaldırmaz.','','Onarım imzası ile kontrol listelerinin doğrudan gen örtüşmesi: '+json.dumps(plan['signature_overlap'])+'.','','| Hasta | Ana model sırası / sütun | Koşul sayısı | İki zon model koşul sayısı |','|---|---:|---:|---:|']
for r in out:md.append(f"| {r['patient']} | {r['zone']['rank']}/{r['zone']['columns']} | {r['zone']['condition']:.2f} | {r['separate_zones']['condition']:.2f} |")
for key in ['raw','cycle','identity','zone','separate_zones']:
 vals=np.array([r[key] if key=='raw' else r[key]['difference'] for r in out]);md.append(f"\n{key}: {sum(vals>0)} pozitif, {sum(vals<0)} negatif; eşit hasta ağırlıklı ortalama {vals.mean():+.4f}.")
md+=['','Katsayılar keşif amaçlı betimleyici ilişkilerdir; bağımsız biyolojik tekrar, nedensel etki veya p değeri değildir. Kimlik ve zonasyon gerçek hastalık biyolojisinin parçası olabilir; bunları çıkarmak aşırı düzeltme yaratabilir. Önceki derinlik örtüşmesi sorunları sürer. Daha çok kovaryat eklemek hücreleri bağımsız hasta yapmaz. Hepatocyte etiketi maligniteyi doğrulamaz.','','[Seurat yöntem kaynağı](https://satijalab.org/seurat/articles/cell_cycle_vignette.html) · [Resmî gen listesi](https://github.com/satijalab/seurat/blob/master/data/cc.genes.rda) · [İnsan atlası](https://www.nature.com/articles/s41586-024-07376-2)']
Path('outputs/genisletilmis-kontroller-raporu.md').write_text('\n'.join(md)+'\n')
print(json.dumps(out,indent=2))
