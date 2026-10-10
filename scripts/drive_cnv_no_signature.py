# Runtime paths: source and working data are separate.
import os as _os
from pathlib import Path as _RuntimePath
_SOURCE_ROOT = _RuntimePath(__file__).resolve().parents[1]
_WORKSPACE = _RuntimePath(_os.environ.get("HCC_WORKSPACE", str(_SOURCE_ROOT))).resolve()
_os.chdir(_WORKSPACE)
(_WORKSPACE / "outputs").mkdir(exist_ok=True)
from pathlib import Path
import subprocess,json,time,sys,shutil
ROOT=_WORKSPACE; target=ROOT/'work/cnv_no_signature'
patients=['HCC07','HCC06','HCC05','HCC10','HCC09','HCC08','HCC04','HCC03']
status={'analysis':'inferCNV excluding 47 signature genes','started':time.strftime('%Y-%m-%d %H:%M:%S'),'patients':{}}
for p in patients:
    status['patients'][p]={}
    for mode,script,marker in [('primary','run_infercnv_no_signature.R','primary/COMPLETED'),('hepatocyte_reference','run_infercnv_hepref_no_signature.R','hepatocyte_reference/COMPLETED'),('classification','classify_cnv_no_signature.R','cell_labels.tsv')]:
        if (target/p/marker).exists():
            status['patients'][p][mode]={'status':'already_complete'};continue
        status['patients'][p][mode]={'status':'running','start':time.time()}; (target/'run_status.json').write_text(json.dumps(status,indent=2))
        log=target/p/(mode+'.log')
        print(p,mode,'START',flush=True)
        with log.open('w') as f:
            result=subprocess.run([shutil.which(_os.environ.get('HCC_RSCRIPT', 'Rscript')) or 'Rscript',str(_SOURCE_ROOT/'scripts'/script),p],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
        status['patients'][p][mode].update({'exit_code':result.returncode,'elapsed_seconds':time.time()-status['patients'][p][mode]['start'],'status':'complete' if result.returncode==0 else 'failed'})
        (target/'run_status.json').write_text(json.dumps(status,indent=2))
        if result.returncode!=0 or not (target/p/marker).exists():
            print(log.read_text()[-4000:],flush=True);raise SystemExit(1)
        print(p,mode,'COMPLETE',flush=True)
status['finished']=time.strftime('%Y-%m-%d %H:%M:%S');(target/'run_status.json').write_text(json.dumps(status,indent=2))
print('ALL EIGHT PATIENTS COMPLETE',flush=True)

with (target/'signature_sensitivity.log').open('w') as f:
    result=subprocess.run([sys.executable,str(_SOURCE_ROOT/'scripts/cnv_signature_sensitivity_no_signature.py')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
if result.returncode: raise SystemExit(result.returncode)

subprocess.run([sys.executable,str(_SOURCE_ROOT/"scripts/compare_cnv_no_signature.py")],cwd=ROOT,check=True)
