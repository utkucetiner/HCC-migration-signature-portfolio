# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
import gzip,csv,json,hashlib,time
from pathlib import Path
matrix=Path('data/external/GSE149614_HCC.scRNAseq.S71915.count.txt.gz')
meta=Path('work/data/GSE149614_HCC.metadata.updated.txt.gz')
with gzip.open(meta,'rt',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
ids={r['Cell'] for r in rows}
genes='AKAP12 ITGB1-DT FMNL2 ITGA2 TNFAIP8 PLOD2 ITGAV BICC1 TPM1 KLHL29 CREB5 DTNA MIR4435-2HG GMDS ANXA2 NDRG1 GNA13 DGKH AL359313.1 RAB3B MKL1 PRKCA AL138828.1 PRAG1 SLC38A1 CYTOR HKDC1 PTPRM ASAP1 FLNB DUSP5 SAMD4A FAT1 UBASH3B MAP4K4 TXNRD1 LPGAT1 LUCAT1 MYH9 BCL2L1 SERPINE1 PALLD CBLB PKN2-AS1 HIF1A ASPH ARHGAP12 FAM13A RELN ITGB5'.split()
assert len(genes)==len(set(genes))==50
seen=set();dupes=[];bad=[];found={};n=0
with gzip.open(matrix,'rt') as f:
 h=f.readline().rstrip('\r\n').split('\t');print('Header first fields:',h[:4],flush=True)
 # Matrix may have an empty gene-label field or header containing only cell IDs.
 cells=h[1:] if h[0] not in ids else h
 assert all(c in ids for c in cells),'Unexpected cell header format'
 for line in f:
  v=line.rstrip('\r\n').split('\t');g=v[0];n+=1
  if len(v)!=len(cells)+1:bad.append((n,g,len(v)))
  if g in seen:dupes.append(g)
  seen.add(g)
  if g in genes:
   counts=[int(z) for z in v[1:]]
   assert min(counts)>=0
   found[g]={'cells_with_nonzero_counts':sum(z>0 for z in counts),'total_counts':sum(counts)}
  if n%5000==0:print('Rows checked:',n,flush=True)
out={'matrix_gene_rows':n,'matrix_cells':len(cells),'unique_matrix_cells':len(set(cells)),'metadata_cells':len(rows),'matrix_only_cells':sorted(set(cells)-ids),'metadata_only_cells':sorted(ids-set(cells)),'malformed_rows':bad,'duplicate_gene_ids':dupes,'signature_size':len(genes),'present_signature_genes':sorted(found),'missing_signature_genes':[g for g in genes if g not in found],'signature_detection':found,'matrix_sha256':hashlib.sha256(matrix.read_bytes()).hexdigest(),'scope':'Structural and signature coverage audit only; no tumour association tested.'}
Path('work/matrix_audit.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='signature_detection'},indent=2),flush=True)
