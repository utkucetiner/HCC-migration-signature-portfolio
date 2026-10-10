# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import numpy as np,json,gzip,csv
from scipy import io,sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
import igraph as ig,leidenalg
from scipy.stats import rankdata
seed=20261010
with gzip.open('work/data/GSE189175_meta_data.tsv.gz','rt') as f:meta=list(csv.DictReader(f,delimiter='\t'))
with gzip.open('work/data/GSE189175_counts_features.tsv.gz','rt') as f:feat=[x.rstrip().split('\t') for x in f]
genes=np.array([x[1] for x in feat]);cells=np.array([r[''] for r in meta]);samples=np.array([r['SampleID'] for r in meta]);patients=np.array([r['Patient'] for r in meta]);tissue=np.array([r['Tissue'] for r in meta])
alias={'LINC00152':'CYTOR','SGK223':'PRAG1'}
signature=[alias.get(g,g) for g in json.loads(Path('work/qc_audit.json').read_text())['verified_signature_mapping'].values()]
controls=json.loads(Path('work/expanded_control_plan.json').read_text())['sets']
excluded=set(signature).union(*map(set,controls.values()))
plan={'seed':seed,'qc':'Existing DIEM-filtered nuclei plus detected genes >=200 and mt fraction <=0.20; assay-specific sensitivity needed','cluster_features':'Top 2000 normalized log expression dispersion mean rank across 6 samples; genes expressed >=20 nuclei; mean normalized log expression >0.02; unique gene symbols; no MT-/RPL/RPS or signature/control genes','embedding':'Centered PCA on clipped standardized 2000 features; 30 PCs; 20 nearest neighbors; Leiden resolution 1.0, no batch correction','annotations':'Canonical marker profiles assessed without examining repair outcomes; inferred lineage, not author barcode annotation or malignancy','signature':signature,'excluded':sorted(excluded)}
Path('work/GSE189175_analysis_plan.json').write_text(json.dumps(plan,indent=2))
print('Reading sparse matrix',flush=True)
with threadpool_limits(limits=4):
 with gzip.open('work/data/GSE189175_counts_matrix.mtx.gz','rb') as f:A=io.mmread(f)
assert A.shape==(58100,39995) and A.nnz==54475478
assert np.all(A.data>0) and np.all(A.data==np.rint(A.data))
X=A.T.tocsr().astype(np.float32);del A
total=np.asarray(X.sum(axis=1)).ravel();det=np.diff(X.indptr);mt=np.asarray(X[:,np.char.startswith(genes,'MT-')].sum(axis=1)).ravel()/total
assert np.all(total>0)
keep=(det>=200)&(mt<=.2)
qc={s:{'before':int(sum(samples==s)),'after':int(sum((samples==s)&keep)),'median_counts':float(np.median(total[samples==s])),'median_genes':float(np.median(det[samples==s]))} for s in np.unique(samples)}
X=X[keep];cells=cells[keep];samples=samples[keep];patients=patients[keep];tissue=tissue[keep];total=total[keep];det=det[keep];mt=mt[keep]
N=X.multiply((10000/total)[:,None]).tocsr();N.data=np.log1p(N.data)
symbols,count=np.unique(genes,return_counts=True);unique=set(symbols[count==1]);detection=np.asarray((N>0).sum(axis=0)).ravel()
means=np.asarray(N.mean(axis=0)).ravel()
eligible=np.array([g in unique and g not in excluded and not g.startswith(('MT-','RPL','RPS')) for g in genes])&(detection>=20)&(means>.02)
ix=np.flatnonzero(eligible);ranks=[]
for s in np.unique(samples):
 a=N[samples==s][:,ix];mu=np.asarray(a.mean(axis=0)).ravel();var=np.asarray(a.power(2).mean(axis=0)).ravel()-mu**2
 dispersion=var/(mu+.01);ranks.append(rankdata(dispersion)/len(ix))
best=ix[np.argsort(np.mean(ranks,axis=0))[-2000:]]
D=N[:,best].toarray();mu=D.mean(axis=0);sd=D.std(axis=0);D=(D-mu)/np.maximum(sd,1e-6);np.clip(D,-10,10,out=D)
print('PCA',len(cells),flush=True)
with threadpool_limits(limits=4):
 model=PCA(n_components=30,svd_solver='randomized',random_state=seed)
 pcs=model.fit_transform(D)
 del D
 knn=NearestNeighbors(n_neighbors=21,n_jobs=4).fit(pcs);dist,neighbors=knn.kneighbors(pcs)
edges=np.column_stack([np.repeat(np.arange(len(cells)),20),neighbors[:,1:].ravel()])
g=ig.Graph(n=len(cells),edges=edges,directed=False);g.simplify()
part=leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,resolution_parameter=1.0,seed=seed)
labels=np.array(part.membership)
sparse.save_npz('work/GSE189175_normalized.npz',N)
np.savez_compressed('work/GSE189175_cluster_inputs.npz',genes=genes,cells=cells,samples=samples,patients=patients,tissue=tissue,total=total,detected=det,mt=mt,labels=labels,pcs=pcs,features=genes[best])
markers={'Hepatocyte':['ALB','TTR','HNF4A','ASGR1','FABP1','ALDOB','APOA1','APOA2'],'Immune':['PTPRC','CD3D','NKG7','LST1','CD79A'],'Endothelial':['PECAM1','VWF','KDR','EMCN'],'Stromal':['COL1A1','COL1A2','DCN','RGS5'],'Cholangiocyte':['KRT19','KRT7','EPCAM','SOX9']}
bygene={g:j for j,g in enumerate(genes) if g in unique};summaries=[]
for k in sorted(np.unique(labels)):
 mask=labels==k;n=int(mask.sum());row={'cluster':int(k),'n':n,'sample_counts':{s:int(sum(mask&(samples==s))) for s in np.unique(samples)},'markers':{}}
 for group,gs in markers.items():
  for gn in gs:
   if gn in bygene:
    v=N[mask,bygene[gn]]
    row['markers'][gn]={'mean':float(v.mean()),'detection':float(v.nnz/n)}
 summaries.append(row)
result={'qc':qc,'clusters':summaries,'n_clusters':len(summaries),'marker_sets':markers,'explained_variance':model.explained_variance_ratio_.tolist()}
Path('work/GSE189175_cluster_audit.json').write_text(json.dumps(result,indent=2))
print('Completed',len(summaries),'clusters',json.dumps(qc),flush=True)
for r in summaries:print(r['cluster'],r['n'],{k:round(np.mean([r['markers'][g]['mean'] for g in v if g in r['markers']]),3) for k,v in markers.items()},r['sample_counts'])
