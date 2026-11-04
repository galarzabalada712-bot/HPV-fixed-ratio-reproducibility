"""Rebuild sex-history and CRP inputs from the included public NHANES XPT files.

The 11,884-person frozen candidate cohort and complete sampling design are inputs;
this script does not claim to recover the original historical selection process.
Run from the package root.
"""
from pathlib import Path
import numpy as np
import pandas as pd

root=Path('.')
d=pd.read_csv(root/'alcohol_four_level_20261102/analysis_input.csv')
membership=pd.read_csv(root/'alcohol_four_level_20261102/fixed_model_membership.csv')
assert len(d)==11884 and d.SEQN.is_unique
parts=[]
for i,letter in enumerate('CDEFGHI'):
 year=2003+2*i; raw=pd.read_sas(root/f'source_cache/SXQ_{letter}.xpt',format='xport')
 age='SXD030' if letter=='C' else 'SXD031'
 partner='SXQ100' if letter=='C' else ('SXQ101' if letter in 'DE' else 'SXD101')
 x=d.loc[d.cycle_start.eq(year),['SEQN']].merge(raw[['SEQN',age,partner]],on='SEQN',how='left',validate='1:1')
 a=x[age];v=x[partner].mask(x[partner].eq(5.397605346934028e-79),0)
 av=a.between(9,69)&a.eq(np.floor(a))&~a.isin([77,99,777,999])
 vv=v.between(0,2000)&v.eq(np.floor(v))&~v.isin([77777,99999])
 x['sexagenew']=np.select([av&a.lt(18),av&a.ge(18)],[1.,2.],default=np.nan)
 x['sextotal']=np.select([vv&v.eq(0),vv&v.eq(1),vv&v.between(2,5),vv&v.between(6,10),vv&v.gt(10)],[0.,1.,2.,3.,4.],default=np.nan)
 parts.append(x[['SEQN','sexagenew','sextotal']])
sex=pd.concat(parts,ignore_index=True); assert len(sex)==11884 and sex.SEQN.is_unique
expected=pd.read_csv(root/'reference/sex_covariates.csv')
assert sex.SEQN.equals(expected.SEQN)
for v in ['sexagenew','sextotal']: assert np.allclose(sex[v],expected[v],equal_nan=True)
sex.to_csv(root/'hpv_sex_covariates_audit_20261102/model_covariates.csv',index=False)
ids=membership.loc[(membership.outcome=='HPV')&(membership.ratio=='PLR'),'SEQN']
q=d[d.SEQN.isin(ids)][['SEQN','cycle_start','HPVnew','activity']].copy();assert len(q)==10230
crp=[]
for year,name,var in [(2003,'L11_C.xpt','LBXCRP'),(2005,'CRP_D.xpt','LBXCRP'),(2007,'CRP_E.xpt','LBXCRP'),(2009,'CRP_F.xpt','LBXCRP'),(2015,'HSCRP_I.xpt','LBXHSCRP')]:
 x=pd.read_sas(root/'nhanes_raw'/name,format='xport')[['SEQN',var]].rename(columns={var:'CRP_raw'})
 x['CRP_mg_L']=x.CRP_raw*(10 if year<2015 else 1);x['source']=name;crp.append(x)
q=q.merge(pd.concat(crp,ignore_index=True),on='SEQN',how='left',validate='1:1')
ref=pd.read_csv(root/'reference/CRP_preflight.csv'); check=q.merge(ref[['SEQN','CRP_mg_L']],on='SEQN',validate='1:1',suffixes=('','_expected'))
assert np.allclose(check.CRP_mg_L,check.CRP_mg_L_expected,equal_nan=True)
q.to_csv(root/'hpv_sync_20261102/fixed_HPV_sample_preflight_by_person.csv',index=False)
print('PASS: 11884 sex-history records and 10230 CRP input records match frozen reference')
