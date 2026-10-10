# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import json,hashlib,gzip,shutil
import numpy as np
from scipy import io,sparse
root=_WORKSPACE;source=root/'work/cnv';target=root/'work/cnv_no_signature';target.mkdir(exist_ok=True)
mapping=json.loads((root/'work/qc_audit.json').read_text())['verified_signature_mapping'];excluded=sorted(set(mapping.values()));assert len(excluded)==47
(target/'excluded_signature_genes.txt').write_text('\n'.join(excluded)+'\n')
order=(source/'gene_order.tsv').read_text().splitlines();(target/'gene_order.tsv').write_text('\n'.join(l for l in order if l.split('\t')[0] not in excluded)+'\n')
manifest={}
for patient in ['HCC07','HCC06','HCC05','HCC10','HCC09','HCC08','HCC04','HCC03']:
 src=source/patient;dest=target/patient;dest.mkdir(exist_ok=True)
 genes=(src/'genes.txt').read_text().splitlines();keep=np.array([g not in excluded for g in genes]);present=sorted(set(genes)&set(excluded));assert len(present)==47
 with gzip.open(src/'counts.mtx.gz','rb') as f: counts=io.mmread(f).tocsc()
 filtered=counts[keep,:]; assert filtered.shape[1]==counts.shape[1] and filtered.min()>=0
 assert np.asarray(filtered.sum(axis=0)).min()>=100, 'Removing genes would alter cell inclusion'
 with gzip.open(dest/'counts.mtx.gz','wb') as f: io.mmwrite(f,filtered,field='integer')
 (dest/'genes.txt').write_text('\n'.join(g for g,k in zip(genes,keep) if k)+'\n')
 for name in ['cells.txt','annotations.tsv']:shutil.copy2(src/name,dest/name)
 manifest[patient]={'before_genes':len(genes),'after_genes':filtered.shape[0],'cells':filtered.shape[1],'excluded_present':present,'minimum_remaining_counts':float(np.asarray(filtered.sum(axis=0)).min()),'source_matrix_sha256':hashlib.sha256((src/'counts.mtx.gz').read_bytes()).hexdigest(),'filtered_matrix_sha256':hashlib.sha256((dest/'counts.mtx.gz').read_bytes()).hexdigest()}
 print(patient, filtered.shape,flush=True)
(target/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
print('Input checks passed: same cells, removed 47 genes, originals preserved.')
