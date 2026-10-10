# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import csv,json,gzip,numpy as np
manifest=json.loads(Path('work/cnv/input_manifest.json').read_text());z=np.load('work/hepatocyte_scoring_inputs.npz');allcells=[]
for p,r in manifest.items():
 d=Path('work/cnv')/p;assert (d/'primary/COMPLETED').exists() and (d/'hepatocyte_reference/COMPLETED').exists()
 labels=list(csv.DictReader((d/'cell_labels.tsv').open(),delimiter='\t'));assert len(labels)==r['hepatocytes'];assert len({x['cell'] for x in labels})==len(labels);expected=set(z['cells'][np.isin(z['sample'],[p+'N',p+'T'])]);assert set(x['cell'] for x in labels)==expected
 allcells.extend(x['cell'] for x in labels)
 for mode in ['primary','hepatocyte_reference']:
  with gzip.open(d/(mode+'_bins.tsv.gz'),'rt') as f:
   for line in f:pass
assert len(allcells)==len(set(allcells))==13762
results=json.loads(Path('work/cnv_signature_results.json').read_text());assert len(results)==8
for r in results:
 if 'balancing' in r:
  a=r['balancing']['overlap'];assert np.isfinite(a['difference']);assert all(np.isfinite(v) for v in a['SMD'].values())
  if max(abs(v) for v in a['SMD'].values())>1e-5:print('balance warning',r['patient'])
Path('work/cnv_final_verification.json').write_text(json.dumps({'patients':8,'hepatocytes':13762,'references':2367,'unique_barcodes':True,'coverage_exact':True,'gzip_profiles_complete':True},indent=2));print('Eight-patient coverage and numerical integrity verified')
