"""Compare portable analysis output with frozen numerical reference values."""
from pathlib import Path
import numpy as np
import pandas as pd

p=Path('.'); ref=p/'reference'; failures=[]
def compare(label,a,b,cols,tol=1e-9):
 assert len(a)==len(b),(label,len(a),len(b))
 for col in cols:
  x=pd.to_numeric(a[col]).to_numpy();y=pd.to_numeric(b[col]).to_numpy()
  diff=np.nanmax(np.abs(x-y));print(label,col,'max abs diff',diff)
  if not np.allclose(x,y,rtol=tol,atol=tol,equal_nan=True):failures.append((label,col,diff))

models=pd.read_csv(p/'alcohol_four_level_20261102/models/model_status.csv');frozen=pd.read_csv(ref/'constraint.csv')
assert len(models)==16 and models.status.eq('fitted').all()
m=models.merge(frozen,on=['outcome','ratio'],validate='1:1',suffixes=('','_ref'))
for left,right in [('n','n_ref'),('events','events_ref'),('D','constraint_deviation'),('se','se_ref'),('CI_low','CI_low_ref'),('CI_high','CI_high_ref'),('p','working_LRT_p'),('Holm_7','Holm_p_7_unique_constraints')]:
 x=m[left].to_numpy();y=m[right].to_numpy();diff=np.max(np.abs(x-y));print('constraint',left,'max abs diff',diff)
 if not np.allclose(x,y,rtol=1e-9,atol=1e-9):failures.append(('constraint',left,diff))
for name,path,columns in [
 ('CLR_complete_case','hpv_sync_20261102/CLR_complete_case_results.csv',['N','events','D','se','CI_low','CI_high','working_LRT_p']),
 ('CLR_MI','hpv_sync_20261102/CLR_MI_results.csv',['N','events','D','se','CI_low','CI_high','Wald_p']),
 ('CLR_2003_2010','clr_unit_audit_20261103/CLR_2003_2010_complete_case.csv',['N','events','D','se','CI_low','CI_high','working_LRT_p'])]:
 compare(name,pd.read_csv(p/path),pd.read_csv(ref/(name+'.csv')),columns)
for name,outfile in [('table1_descriptive','Table1_descriptive_preliminary.csv'),('table1_tests','Table1_design_tests_preliminary.csv')]:
 a=pd.read_csv(p/'hpv_sex_covariates_audit_20261102'/outfile);b=pd.read_csv(ref/(name+'.csv'))
 compare(name,a,b,['n','weighted_percent'] if name.endswith('descriptive') else ['p'])
pilot=pd.read_csv(p/'hpv_sex_covariates_audit_20261102/pilot_10/HPV_performance.csv');f=pd.read_csv(ref/'HPV_performance.csv')
assert list(pilot.marker)==list(f.marker) and (pilot.B_valid==10).all()
compare('performance pilot point',pilot,f,['N','events','auc_composite','auc_components','delta_auc','nri','idi'])
compare('DCA covariate baseline',pd.read_csv(p/'hpv_sex_covariates_audit_20261102/HPV_dca_covariates_only.csv'),pd.read_csv(ref/'HPV_dca_covariates_only.csv'),['threshold','covariates_only'])
compare('DCA 14 pilot point curves',pd.read_csv(p/'hpv_sex_covariates_audit_20261102/pilot_10/HPV_decision_curves.csv'),pd.read_csv(ref/'HPV_decision_curves.csv'),['threshold','composite','components','treat_all','treat_none'])
rep=pd.read_csv(ref/'HPV_bootstrap_replicates.csv');assert len(rep)==14000
for row in f.itertuples():
 z=rep[rep.marker==row.marker];assert len(z)==1000 and z[['delta_auc','nri','idi']].notna().all().all()
 for metric in ['delta_auc','nri','idi']:
  lo,hi=z[metric].quantile([.025,.975]).to_numpy()
  assert np.allclose([lo,hi],[getattr(row,metric+'_low'),getattr(row,metric+'_high')],rtol=1e-9,atol=1e-9)
assert not failures,failures
print('PASS: 16 models, Table 1, CLR/MI/early sensitivity, 14 pilot point estimates and saved 1000-replicate percentile intervals')
