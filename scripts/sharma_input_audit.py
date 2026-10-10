# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import h5py,numpy as np,json
from pathlib import Path
meta=json.loads(Path('work/sharma_metadata_audit.json').read_text())
plan=json.loads(Path('work/expanded_control_plan.json').read_text())
mapping=json.loads(Path('work/qc_audit.json').read_text())['verified_signature_mapping']
with h5py.File('work/data/Sharma_HCC.h5ad') as f:
 obs=f['obs'][:];genes=[x.decode() for x in f['raw.var']['index']];lookup={g:i for i,g in enumerate(genes)};assert len(genes)==len(set(genes))
 coverage={'repair47':{'present':[g for g in mapping.values() if g in lookup],'missing':[g for g in mapping.values() if g not in lookup]}}
 for k,v in plan['sets'].items():coverage[k]={'present':[g for g in v if g in lookup],'missing':[g for g in v if g not in lookup]}
 markers='ALB TTR TF HP HNF4A CYP2A6 APOA1 APOA2 ASGR1 FABP1 PTPRC CD3D NKG7 LST1 PECAM1 VWF COL1A1 KRT19 EPCAM'.split()
 target=set(markers)|set().union(*(set(v['present']) for v in coverage.values()))
 selected=sorted(target&set(lookup));out=np.zeros((len(obs),len(selected)));gmap=np.full(len(genes),-1,int)
 for j,g in enumerate(selected):gmap[lookup[g]]=j
 ptr=f['raw.X/indptr'][:];total=np.zeros(len(obs));det=np.diff(ptr);maxerr=0.
 for start in range(0,len(obs),2000):
  end=min(start+2000,len(obs));a,b=int(ptr[start]),int(ptr[end])
  d=f['raw.X/data'][a:b].astype(float);idx=f['raw.X/indices'][a:b]
  recovered=np.expm1(d);assert np.isfinite(recovered).all() and np.all(recovered>=0)
  maxerr=max(maxerr,float(np.max(abs(recovered-np.rint(recovered)))))
  cellidx=np.repeat(np.arange(start,end),np.diff(ptr[start:end+1]))
  np.add.at(total,cellidx,recovered)
  cols=gmap[idx];mask=cols>=0;out[cellidx[mask],cols[mask]]=recovered[mask]
 assert np.all(total>0)
 relative=np.abs(total-obs['n_counts'])/np.maximum(1,obs['n_counts'])
 normalized=np.log1p(out/total[:,None]*10000)
 np.savez_compressed('work/sharma_selected_inputs.npz',genes=np.array(selected),counts=out,total=total,normalized=normalized,obs=obs)
 marker_summary=[]
 names=f['uns/rank_genes_groups/names'][:]
 for k in range(29):
  ix=obs['louvain']==k
  marker_summary.append({'cluster':str(k),'n':int(ix.sum()),'top8':[x.decode() for x in names[str(k)][:8]],'marker_means':{g:float(normalized[ix,selected.index(g)].mean()) for g in markers if g in selected}})
 meta.update(coverage=coverage,raw_recovery={'max_integer_error':maxerr,'max_relative_total_error':float(relative.max()),'median_relative_total_error':float(np.median(relative))},cluster_markers=marker_summary)
 Path('work/sharma_metadata_audit.json').write_text(json.dumps(meta,indent=2))
 print(json.dumps({'coverage':coverage,'raw_recovery':meta['raw_recovery']},indent=2))
