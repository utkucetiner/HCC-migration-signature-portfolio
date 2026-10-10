"""Run one documented stage; never silently run the full study."""
import argparse, os, subprocess, sys, shutil
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("script", help="Filename under scripts/; see docs/CODE_INDEX.md")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--patient", action="append", default=[], help="R patient argument; repeat as needed")
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    script = source / "scripts" / args.script
    if Path(args.script).name != args.script or not script.is_file() or script.suffix not in (".py", ".R"):
        parser.error("Choose an existing .py or .R filename under scripts/.")
    workspace = args.workspace.resolve()
    if not workspace.is_dir(): parser.error("Workspace does not exist.")
    env = dict(os.environ, HCC_WORKSPACE=str(workspace), HCC_SOURCE_ROOT=str(source))
    if script.suffix == ".py": executable = sys.executable
    else:
        executable = shutil.which(env.get("HCC_RSCRIPT", "Rscript"))
        if not executable: parser.error("Rscript is not available; install R or set HCC_RSCRIPT.")
    return subprocess.run([executable, str(script), *args.patient], cwd=source, env=env).returncode

if __name__ == "__main__": sys.exit(main())
