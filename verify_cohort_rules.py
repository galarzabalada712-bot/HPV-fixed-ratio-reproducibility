"""Check recoverable covariate derivations against the frozen cohort input."""
import numpy as np
import pandas as pd
from cohort_rules import derive, derive_adjustments, classify

d=pd.read_csv('alcohol_four_level_20261102/analysis_input.csv')
early=derive(d); later=derive_adjustments(d); alcohol=classify(d)
comparisons=[('diabetes',early.diabetes_analysis,d.diabetesnew),
 ('leisure_activity',early.activity_leisure_analysis,d.activity),
 ('education',later.education,d.education),('PIR',later.PIR,d.PIR),
 ('smoking',later.smokingnew,d.smokingnew),('hypertension',later.hypertentionnew,d.hypertentionnew),
 ('followup_years',later.FOLLOWUP_YEARS,d.FOLLOWUP_YEARS),('alcohol_four_levels',alcohol.Alcoholuse,d.Alcoholuse)]
for name,a,b in comparisons:
 same=np.isclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),equal_nan=True,rtol=1e-12,atol=1e-12)
 print(name,'matched',int(same.sum()),'of',len(same))
 assert same.all(),name
marital=np.select([d.DMDMARTL.eq(5),d.DMDMARTL.isin([1,6]),d.DMDMARTL.isin([2,3,4])],[1,2,3],default=4)
assert np.array_equal(marital,d.marital)
assert np.allclose(d.WTfinal,d.WTMEC2YR/7,atol=1e-9,rtol=1e-12)
print('PASS: current marital category and seven-cycle weight match frozen data')
