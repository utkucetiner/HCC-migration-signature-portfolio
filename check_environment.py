"""Check interpreter and dependencies without opening study data."""
import importlib.metadata, json, shutil, sys
from pathlib import Path
root=Path(__file__).resolve().parent
expected={line.split('==')[0]:line.split('==')[1] for line in (root/'environment/requirements.txt').read_text().splitlines() if '==' in line}
packages={};ok=sys.version_info[:2]==(3,11)
for name,version in expected.items():
    try: actual=importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError: actual=None
    packages[name]={'expected':version,'installed':actual,'matches':actual==version}
    ok=ok and actual==version
print(json.dumps({'python':sys.version.split()[0],'tested_python_minor':'3.11','packages':packages,'Rscript_available':shutil.which('Rscript') is not None,'python_environment_matches':ok,'scope':'Does not validate biological data, R package dependencies or a full CNV rerun.'},indent=2))
sys.exit(0 if ok else 1)
