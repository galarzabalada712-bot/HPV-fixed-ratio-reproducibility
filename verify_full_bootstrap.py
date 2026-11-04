"""Verify a newly computed full 1000-replicate run against frozen results."""
from pathlib import Path
import numpy as np
import pandas as pd
root=Path('.')
full=root/'hpv_sex_covariates_audit_20261102/full';ref=root/'reference'
for name,key in [('HPV_performance.csv',['marker']),('HPV_bootstrap_replicates.csv',['marker','replicate']),('HPV_decision_curves.csv',['marker','threshold'])]:
 a=pd.read_csv(full/name);b=pd.read_csv(ref/name)
 assert list(a.columns)==list(b.columns) and len(a)==len(b),name
 assert a[key].equals(b[key]),name
 for col in a.select_dtypes(include='number').columns:
  assert np.allclose(a[col],b[col],atol=1e-9,rtol=1e-9,equal_nan=True),(name,col)
 print('PASS full',name,len(a),'rows')
