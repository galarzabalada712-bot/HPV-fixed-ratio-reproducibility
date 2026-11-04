"""Portable historical derivation functions applied to the frozen study input.
Final alcohol and marital rules were changed later; see verify_cohort_rules.py.
"""
import numpy as np
import pandas as pd

COMPONENTS = ['LBXRDW', 'LBXWBCSI', 'LBDNENO', 'LBDLYMNO', 'LBDMONO', 'LBXPLTSI', 'LBXHGB', 'LBXSAL']

def normalize_xport_zeros(d):
    # pandas XPORT zero-decoding artifact, independently checked with haven.
    # Exact sentinel only: never blanket-round small genuine measurements.
    return d.replace(5.397605346934028e-79, 0.0)

def derive(d):
    d = d.copy()
    def col(name):
        return pd.to_numeric(d.get(name, pd.Series(np.nan, index=d.index)), errors='raise')
    diagnosed = col('DIQ010').eq(1)
    medication = col('DIQ050').eq(1) | col('pills').eq(1)
    gh = col('LBXGH').where(col('LBXGH').gt(0) & col('LBXGH').le(100))
    lab_positive = gh.ge(6.5)
    med_unknown = col('DIQ050').notna() & ~col('DIQ050').isin([1, 2])
    med_unknown |= col('pills').notna() & ~col('pills').isin([1, 2])
    negative = col('DIQ010').eq(2) & gh.lt(6.5) & ~med_unknown
    d['diabetes_corrected'] = np.select([diagnosed | medication | lab_positive, negative], [1, 0], default=2)
    d['diabetes_reason'] = np.select([diagnosed, medication, lab_positive, negative, col('DIQ010').eq(3), gh.isna()], ['reported_diagnosis', 'reported_medication', 'hba1c_ge_6_5', 'explicit_no_and_hba1c_below_6_5', 'borderline_without_positive', 'hba1c_unavailable_without_positive'], default='insufficient_or_unknown_questionnaire')
    d['diabetes_borderline_reported'] = col('DIQ010').eq(3)
    d['diabetes_medication_positive'] = medication
    d['diabetes_hba1c_positive'] = lab_positive

    lifetime = col('ALQ101').eq(1) | col('ALQ110').eq(1)
    low_lifetime = col('ALQ110').eq(2)
    conflict = col('ALQ101').eq(1) & low_lifetime
    q, u = col('ALQ120Q'), col('ALQ120U')
    valid_q = (u.eq(1) & q.between(1, 7)) | (u.eq(2) & q.between(1, 31)) | (u.eq(3) & q.between(1, 366))
    valid_q &= q.eq(np.floor(q))
    d['alcohol_status'] = np.select([conflict, low_lifetime, lifetime & q.eq(0), lifetime & valid_q], [3, 0, 1, 2], default=3)
    d['alcohol_reason'] = np.select([conflict, low_lifetime, lifetime & q.eq(0), lifetime & valid_q, col('_ALQ_present').ne(1)], ['conflicting_lifetime_answers', 'lifetime_under_12_not_never', 'lifetime_ge_12_past_year_zero', 'lifetime_ge_12_past_year_positive', 'not_in_public_ALQ_file'], default='insufficient_or_invalid_response')
    d['alcohol_ALQ110_structural_skip'] = col('ALQ101').eq(1) & col('ALQ110').isna()
    # Keep intensity as a separate descriptor, not a heavy-drinking diagnosis.
    drinks = col('ALQ130')
    cap = pd.Series(np.where(col('cycle_start').eq(2015), 15, 95), index=d.index)
    drinks = drinks.where(drinks.ge(1) & drinks.le(cap) & ~drinks.isin([77, 99, 777, 999]))
    d['drinks_per_drinking_day'] = drinks.where(d.alcohol_status.eq(2))

    early = col('cycle_start').isin([2003, 2005])
    vigorous = col('PAD200').where(early, col('PAQ650'))
    moderate = col('PAD320').where(early, col('PAQ665'))
    activity_yes = vigorous.eq(1) | moderate.eq(1)
    activity_no = vigorous.eq(2) & moderate.eq(2)
    d['activity_leisure_any'] = np.select([activity_yes, activity_no], [1, 0], default=2)
    d['activity_reason'] = np.select([activity_yes, activity_no, vigorous.eq(3) | moderate.eq(3)], ['any_leisure_yes', 'both_leisure_no', 'unable_response'], default='insufficient_or_unknown_response')
    d['activity_recall_period'] = np.where(early, 'past_30_days', 'typical_week')
    domains = []
    for name, question, days, minutes in [('vigorous_work', 'PAQ605', 'PAQ610', 'PAD615'), ('moderate_work', 'PAQ620', 'PAQ625', 'PAD630'), ('vigorous_recreation', 'PAQ650', 'PAQ655', 'PAD660'), ('moderate_recreation', 'PAQ665', 'PAQ670', 'PAD675')]:
        response, freq, duration = col(question), col(days), col(minutes)
        valid = response.eq(1) & freq.between(1, 7) & freq.eq(np.floor(freq)) & duration.ge(10) & duration.lt(1440)
        total = (freq * duration).where(valid)
        total = total.mask(response.eq(2), 0).mask(early)
        output = name + '_minutes_week'
        d[output] = total
        domains.append(output)
    v = d[domains[0]] + d[domains[2]]
    m = d[domains[1]] + d[domains[3]]
    d['work_recreation_equivalent_minutes'] = m + 2 * v
    total = d.work_recreation_equivalent_minutes
    weekly = pd.Series(np.select([total.eq(0), total.gt(0) & total.lt(150), total.ge(150)], [0, 1, 2], default=3), index=d.index, dtype=float)
    d['activity_weekly_category'] = weekly.mask(early)
    # Explicit NA versions prevent unknown codes from becoming real categories.
    d['diabetes_analysis'] = d.diabetes_corrected.where(d.diabetes_corrected.ne(2))
    d['alcohol_status_analysis'] = d.alcohol_status.where(d.alcohol_status.ne(3))
    d['activity_leisure_analysis'] = d.activity_leisure_any.where(d.activity_leisure_any.ne(2))
    return d

def eligibility(d):
    d = d.copy()
    d['female_18_59'] = d.RIAGENDR.eq(2) & d.RIDAGEYR.between(18, 59)
    d['valid_hpv'] = d.LBDRPCR.isin([1, 2])
    d['valid_components'] = (np.isfinite(d[COMPONENTS]) & d[COMPONENTS].gt(0)).all(axis=1) & d.LBDMONO.ge(1e-10)
    d['available_components'] = (np.isfinite(d[COMPONENTS]) & d[COMPONENTS].ge(0)).all(axis=1)
    d['valid_bmi'] = np.isfinite(d.BMXBMI) & d.BMXBMI.gt(0)
    d['valid_design'] = np.isfinite(d.WTMEC2YR) & d.WTMEC2YR.gt(0) & d[['SDMVSTRA', 'SDMVPSU']].notna().all(axis=1)
    d['candidate_eligible'] = d[['female_18_59', 'valid_hpv', 'valid_components', 'valid_bmi', 'valid_design']].all(axis=1)
    d['candidate_available_components'] = d[['female_18_59', 'valid_hpv', 'available_components', 'valid_bmi', 'valid_design']].all(axis=1)
    d['exclusion_reason'] = np.select([~d.female_18_59, ~d.valid_hpv, ~d.valid_components, ~d.valid_bmi, ~d.valid_design], ['outside_female_18_59', 'hpv_not_valid_1_or_2', 'component_missing_or_nonpositive', 'bmi_missing_or_nonpositive', 'invalid_weight_or_design'], default='candidate_eligible_not_final_model_sample')
    return d

def derive_adjustments(d):
    d = d.copy()
    def c(name):
        return pd.to_numeric(d.get(name, pd.Series(np.nan, index=d.index)))
    youth = c('RIDAGEYR').between(18, 19)
    ed_y = c('DMDEDUC3')
    ed_a = c('DMDEDUC2')
    d['education'] = np.select([
        youth & (ed_y.between(0, 12) | ed_y.isin([55, 66])),
        youth & ed_y.isin([13, 14]), youth & ed_y.eq(15),
        ~youth & ed_a.isin([1, 2]), ~youth & ed_a.eq(3),
        ~youth & ed_a.isin([4, 5])], [1, 2, 3, 1, 2, 3], default=np.nan)
    pir = c('INDFMPIR')
    d['PIR'] = np.select([pir.between(0, 1.3), pir.gt(1.3) & pir.le(3.5),
                          pir.gt(3.5) & pir.le(5)], [1, 2, 3], default=np.nan)
    d['smokingnew'] = np.select([c('SMQ020').eq(2),
        c('SMQ020').eq(1) & c('SMQ040').eq(3),
        c('SMQ020').eq(1) & c('SMQ040').isin([1, 2])], [0, 1, 2], default=np.nan)
    d['marital'] = c('DMDMARTL').map({1:1, 6:1, 5:2, 2:3, 3:3, 4:3})
    qpos = c('BPQ020').eq(1) | c('BPQ040A').eq(1)
    qneg = c('BPQ020').eq(2) | c('BPQ040A').eq(2)
    d['hypertension_questionnaire'] = np.select([qpos, qneg], [1, 0], default=2)
    for stem, label in [('BPXSY', 'Systolic'), ('BPXDI', 'Diastolic')]:
        readings = pd.concat([c(stem + str(i)) for i in range(1, 5)], axis=1)
        d[label] = readings.mean(axis=1)
    bpos = d.Systolic.ge(140) | d.Diastolic.ge(90)
    bneg = d.Systolic.lt(140) & d.Diastolic.lt(90)
    d['hypertension_measured'] = np.select([bpos, bneg], [1, 0], default=2)
    d['hypertension_code'] = np.select([
        d.hypertension_questionnaire.eq(1) | d.hypertension_measured.eq(1),
        d.hypertension_questionnaire.eq(0) | d.hypertension_measured.eq(0)], [1, 0], default=2)
    d['hypertentionnew'] = d.hypertension_code.where(d.hypertension_code.ne(2))
    d['MORTSTAT'] = c('MORTSTAT')
    d['mortality_eligible'] = c('ELIGSTAT').eq(1) & c('MORTSTAT').isin([0, 1]) & c('PERMTH_INT').gt(0)
    d['FOLLOWUP_YEARS'] = c('PERMTH_INT') / 12
    return d

def classify(d):
    q=d.ALQ130
    special=(d.cycle_start.eq(2003)&q.isin([77,99]))|(~d.cycle_start.eq(2003)&q.isin([777,999]))
    cap=np.where(d.cycle_start.eq(2015),15,95)
    valid=q.ge(1)&q.le(cap)&q.eq(np.floor(q))&~special
    category=np.select([special,d.ALQ110.eq(2),d.ALQ110.eq(1)&valid&q.eq(1),
                        d.ALQ110.eq(1)&valid&q.gt(1)],[3,0,1,2],default=3).astype(int)
    reason=np.select([special,d.ALQ110.eq(2),d.ALQ110.eq(1)&valid&q.eq(1),
                      d.ALQ110.eq(1)&valid&q.gt(1),d.ALQ101.eq(1)&d.ALQ110.isna(),
                      d.ALQ110.eq(1)&d.ALQ120Q.eq(0),d.ALQ110.eq(1)],
                     ['ALQ130_special_cleaned','ALQ110_below_12','ALQ130_one_drink',
                      'ALQ130_above_one','ALQ110_structural_skip_ALQ101_yes',
                      'ALQ110_yes_past_year_zero','ALQ110_yes_intensity_unavailable'],
                     default='ALQ110_missing_refused_or_other')
    return pd.DataFrame({'Alcoholuse':category,'alcohol_four_reason':reason,
                         'special_130':special,'valid_130':valid},index=d.index)

