"""Run the five-repetition simulation pilot in a temporary copy."""
import argparse,shutil,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--rscript',default=shutil.which('Rscript'));args=p.parse_args()
if not args.rscript:p.error('Rscript not found; provide --rscript')
root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='hpv_simulation_pilot_') as temp:
 work=Path(temp);shutil.copytree(root/'simulation',work/'simulation')
 subprocess.run([args.rscript,'simulation/simulation_master.R','--stage','1','--no-summarize','--no-tables'],cwd=work,check=True)
print('PASS: simulation pilot ran in a temporary copy; frozen simulation files were not changed')
