"""Run the portable analysis from the frozen, documented cohort boundary."""
import argparse, os, shutil, subprocess, sys
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--rscript',default=shutil.which('Rscript'),help='Path to Rscript if it is not on PATH')
p.add_argument('--workers',type=int,default=1)
p.add_argument('--full-bootstrap',action='store_true',help='Run all 1000 replicates per marker; default is a 10-replicate smoke test')
p.add_argument('--simulation-stage4',action='store_true',help='Rerun all 24 x 1000 simulation scenarios; lengthy')
args=p.parse_args()
if not args.rscript: p.error('Rscript not found; provide --rscript')
if args.workers<1:p.error('--workers must be positive')
root=Path(__file__).resolve().parent;os.chdir(root)
def run(cmd):
 print('\n>>>',subprocess.list2cmdline([str(x) for x in cmd]),flush=True)
 subprocess.run([str(x) for x in cmd],check=True)
run([sys.executable,'restore_inputs.py'])
run([sys.executable,'prepare_inputs.py'])
run([sys.executable,'verify_cohort_rules.py'])
for script in ['refit_alcohol_four_level_20261102.R','hpv_extra_covariates_table1_20261102.R','hpv_sync_clr_20261102.R','audit_clr_early.R']:
 run([args.rscript,script])
run([args.rscript,'hpv_extra_covariates_performance_20261102.R','10',str(args.workers)])
run([args.rscript,'recompute_dca_baseline.R'])
run([args.rscript,'verify_simulation.R'])
run([sys.executable,'verify_results.py'])
if args.full_bootstrap:
 run([args.rscript,'hpv_extra_covariates_performance_20261102.R','1000',str(args.workers)])
 run([sys.executable,'verify_full_bootstrap.py'])
else:
 shutil.copy2('reference/HPV_performance.csv','hpv_sex_covariates_audit_20261102/full/HPV_performance.csv')
 shutil.copy2('reference/HPV_decision_curves.csv','hpv_sex_covariates_audit_20261102/full/HPV_decision_curves.csv')
run([args.rscript,'hpv_extra_covariates_figures_20261102.R'])
if args.simulation_stage4:
 run([args.rscript,'simulation/simulation_master.R','--stage','4','--workers',str(args.workers)])
print('\nPASS: requested steps completed. Pilot bootstrap intervals are not final estimates.')
