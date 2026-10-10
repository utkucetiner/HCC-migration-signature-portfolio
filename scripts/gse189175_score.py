# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import numpy as np,json
from scipy import sparse
z=np.load('work/GSE189175_cluster_inputs.npz')
N=sparse.load_npz('work/GSE189175_normalized.npz')
genes=z['genes'];unique,counts=np.unique(genes,return_counts=True);dupes=set(unique[counts>1])
lookup={g:i for i,g in enumerate(genes)}
plan=json.loads(Path('work/GSE189175_analysis_plan.json').read_text())
source=json.loads(Path('work/expanded_control_plan.json').read_text())
core=[0,1,3,4,5,9,10,11,12];mixed=[15,16];broad=core+mixed
annotation={'core_clusters':core,'mixed_sensitivity_clusters':mixed,'excluded':{'2':'endothelial-like','6':'neural-like / ambiguous hepatic ambient signal','7':'myeloid-like','8':'lymphoid-like','13':'stromal-like','14':'cholangiocyte-like','17':'TTR+/HTR2C+, ambiguous lineage'},'basis':'Multi-marker cluster profiles, before repair outcome contrasts. Inferred labels, not supplied author annotations. ALB alone not used. Mixed clusters not assumed doublets or malignant.'}
Path('work/GSE189175_annotation_decision.json').write_text(json.dumps(annotation,indent=2))
scoring_sets={'repair':plan['signature'],**source['sets']}
scores={}
for k,gs in scoring_sets.items():
 present=[g for g in gs if g in lookup]
 assert not set(present)&dupes
 scores[k]=np.asarray(N[:,[lookup[g] for g in present]].mean(axis=1)).ravel()
markers=['HNF4A','ASGR1','FABP1','ALDOB']
positive=np.column_stack([np.asarray(N[:,lookup[g]].toarray()).ravel()>0 for g in markers]).sum(axis=1)
negative=np.column_stack([np.asarray(N[:,lookup[g]].toarray()).ravel()>0 for g in ['PTPRC','PECAM1','VWF','KRT19','COL1A1']]).sum(axis=1)
coremask=np.isin(z['labels'],core)
masks={'core':coremask,'mixed_included':np.isin(z['labels'],broad),'strict_markers':coremask&(positive>=2)&(negative==0)}
# Selection decisions and sensitivity rules stored before fit outcomes.
annotation['strict_rule']='Within core clusters, >=2 of HNF4A/ASGR1/FABP1/ALDOB detected and none of PTPRC/PECAM1/VWF/KRT19/COL1A1 detected. Alters target population and dropout/depth selection.'
Path('work/GSE189175_annotation_decision.json').write_text(json.dumps(annotation,indent=2))
C=np.column_stack([np.log1p(z['total']),scores['S'],scores['G2M'],scores['identity'],scores['portal']-scores['central']])
out=[]
def fit(ix,tn):
 a=np.flatnonzero(~tn);b=np.flatnonzero(tn)
 weights=np.where(tn,1/len(b),1/len(a));w=np.sqrt(weights)
 ret={}
 for model,nc in [('cycle',3),('identity',4),('zone',5)]:
  c=C[ix,:nc];c=c-np.average(c,axis=0,weights=weights);sd=np.sqrt(np.average(c**2,axis=0,weights=weights))
  if np.any(sd<1e-10):ret[model]={'error':'constant_covariate'};continue
  X=np.column_stack([np.ones(len(ix)),tn.astype(float),c/sd]);coef,res,rank,sv=np.linalg.lstsq(X*w[:,None],scores['repair'][ix]*w,rcond=None)
  ret[model]={'difference':float(coef[1]),'rank':int(rank),'columns':int(X.shape[1]),'condition':float(sv[0]/sv[-1])}
 return ret
for patient in ['1','2','3']:
 row={'patient':patient,'selections':{}}
 for label,mask in masks.items():
  sel=mask&(z['patients']==patient);a=np.flatnonzero(sel&(z['tissue']=='AdjNonTumor'));b=np.flatnonzero(sel&(z['tissue']=='Tumor'))
  rr={'n_N':len(a),'n_T':len(b)}
  if min(len(a),len(b))>=20:
   rr['raw']=float(scores['repair'][b].mean()-scores['repair'][a].mean());ix=np.r_[a,b];tn=np.r_[np.zeros(len(a),bool),np.ones(len(b),bool)]
   rr['models']=fit(ix,tn)
   # Marginal common-range sensitivity, same five-control criterion as discovery.
   retain=np.ones(len(ix),bool)
   for j in range(5):
    qa=np.quantile(C[a,j],[.05,.95]);qb=np.quantile(C[b,j],[.05,.95]);lo=max(qa[0],qb[0]);hi=min(qa[1],qb[1])
    retain&=(C[ix,j]>=lo)&(C[ix,j]<=hi)
   rr['common_range_n_N']=int(sum(retain&~tn));rr['common_range_n_T']=int(sum(retain&tn))
   if min(rr['common_range_n_N'],rr['common_range_n_T'])>=20:rr['common_range_models']=fit(ix[retain],tn[retain])
  row['selections'][label]=rr
 out.append(row)
np.savez_compressed('work/GSE189175_scores.npz',**scores,coremask=coremask,strict_mask=masks['strict_markers'])
results={'annotation':annotation,'patients':out,'cell_count':int(sum(coremask)),'limitations':['n=3 independent patients','reannotated lineage, not author barcode labels','no DNA-confirmed malignancy','no batch correction in exploratory clustering','ambient RNA/doublet risk unresolved','scRNA versus snRNA assay difference','descriptive within-patient regressions, not causal or inferential effects']}
Path('work/GSE189175_score_results.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
