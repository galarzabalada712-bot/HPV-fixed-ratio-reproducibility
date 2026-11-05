# HPV fixed-ratio study: reproducibility package

## Scope and provenance

This directory contains the corrected post-hoc analysis of NHANES 2003–2016, starting from the **frozen reconstructed candidate cohort of 11,884 participants** and the full survey design frame. It is independent of the original project directory and its historical analysis directory. The initial historical eligibility/filtering sequence was not recovered: this package therefore does **not** claim that it reconstructs the original historical cohort from raw NHANES files. The frozen input is included, alongside its field-level dictionary, original-source URLs and a description of its derivation. This package is publicly available in the GitHub repository. Users must comply with applicable NHANES data-use conditions described in `DATA_USE.md`. No repository DOI or repository-wide licence is specified.

The 12 included public NHANES XPT files are the seven sexual-behavior files and five CRP files actually read by `prepare_inputs.py` and `audit_clr_early.R`. The `reference/official_source_list_103.csv` file lists the 103 study sources with official URLs and recorded SHA-256 checksums; the other 91 source files are not required to run the frozen-input pipeline and are not bundled. Whether every supporting historical derivation dependency is listed has not been independently established. Do not imply that all 103 raw files are present in this package.

## Software

Tested on Windows 11 with R 4.6.0, survey 4.5, survival 3.8-6, mice 3.19.0, haven 2.5.5 and Python 3.12.14, pandas 3.0.1, NumPy 2.3.5. Plotting additionally uses ggplot2, tidyr and svglite. Simulation reruns require MASS, pROC, dplyr, tidyr, purrr and ggplot2; future/furrr are optional. R package versions observed in the CLR run are in `reference/CLR_sessionInfo.txt`; simulation versions are in `simulation/output/summaries/session_info.txt`. Package versions may affect numerical reproducibility. Some packages emitted R version-build warnings, and R emitted locale startup warnings in testing; these were not model-fitting warnings.

## Quick reproducibility run

Make a writable copy of the directory. Before running, check the saved hashes with `python verify_manifest.py` (the analysis overwrites generated output files and will then change those hashes). From any shell directory, invoke:

```powershell
python path/to/package/run_all.py --rscript "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" --workers 2
```

Omit `--rscript` if Rscript is on PATH. Python must provide pandas and NumPy; `read_sas` must support XPT. The driver locates its own root and runs all steps in that directory. If running the component scripts individually, run `python restore_inputs.py` from the package root first. It overwrites generated output files within the copy; frozen references remain under `reference/`. It performs:

1. Decompress the two frozen input CSV files, verify their SHA-256 hashes, then derive and check sexual-history categories and CRP input from the included XPT files (`prepare_inputs.py`). The original cohort and survey frame are inputs, not recomputed. `cohort_rules.py` and `verify_cohort_rules.py` independently reapply and check current diabetes, activity, education, PIR, smoking, hypertension, follow-up, alcohol, marital and weight derivations against all 11,884 frozen records.
2. Refit eight named HPV constraints and eight exploratory mortality constraints with paired samples, survey design, marital strata for mortality, and seven-constraint Holm correction (`refit_alcohol_four_level_20261102.R`).
3. Recompute Table 1 statistics, CLR complete cases, 20-imputation measured-cycle analysis and post-hoc 2003–2010 sensitivity (`hpv_extra_covariates_table1_20261102.R`, `hpv_sync_clr_20261102.R`, `audit_clr_early.R`).
4. Run 14 HPV performance comparisons with 10 paired bootstrap replicates as a smoke test, refit the covariate-only decision curve and compare point estimates/curves to frozen results. **The 10-replicate intervals are not the manuscript intervals.**
5. Independently reaggregate the saved 24 simulation RDS files and check their published summaries (`verify_simulation.R`), check the model estimates, output samples, saved 1,000-replicate intervals and table statistics (`verify_results.py`), then regenerate the HPV plots from frozen performance results and decision-curve sources.

The full quick pipeline was successfully run from an isolated copy under a different parent directory, with 16 constraints, all CLR runs, Table 1, 14 performance point estimates, DCA, 24 simulation raw files and the manuscript bootstrap percentile intervals checked against frozen output. A separate simulation stage-1 smoke test (one scenario, five repetitions) also passed. No new 1,000-replicate full HPV bootstrap or complete stage-4 simulation was run during packaging.

## Full optional reruns

Add `--full-bootstrap` to rerun 1,000 paired survey-bootstrap replicates for all 14 HPV comparisons; the driver runs `verify_full_bootstrap.py` against the frozen replicate and summary files. This can take substantial computation time. For a 10-replicate pilot, the pseudorandom sequence differs from a full 1,000-replicate run because the replicate weights are generated in different sized batches; therefore only point estimates and DCA, not the pilot intervals or replicate-by-replicate values, are compared to the full run.

To run the five-repetition simulation pilot safely, use `python run_simulation_pilot.py --rscript "C:/Program Files/R/R-4.6.0/bin/Rscript.exe"`; it operates in a temporary copy. Running `simulation_master.R --stage 1` directly inside the package overwrites a same-named stage-4 RDS even with `--no-summarize --no-tables`.

Add `--simulation-stage4` to rerun the 24 scenarios with 1,000 repetitions each. This is a long operation and writes into `simulation/output/`, including raw simulation data and summaries; use a separate copy to preserve the original packaged simulations. The original stage-4 run was previously recorded at approximately 11 hours 43 minutes with four workers, but runtime depends on hardware and the current package was not timed for stage 4. `verify_simulation.R` checks existing raw RDS and summaries; it does not itself rerun the simulation.

## Data and result locations

- `alcohol_four_level_20261102/analysis_input.csv.gz`: compressed frozen reconstructed cohort, restored to `.csv` by the driver; 133 columns including public source fields and derived audit fields.
- `alcohol_missing_reanalysis/survey_design_frame.csv.gz`: compressed full 68,303-record frame with weights, strata and PSU, used before selecting analysis domains. The directory name is a historical label; the current analysis uses the four-level alcohol coding in `analysis_input.csv`.
- `alcohol_four_level_20261102/fixed_model_membership.csv`: frozen SEQN membership for each outcome and ratio; identifiers are NHANES public sequence IDs.
- `source_cache/SXQ_*.xpt` and `nhanes_raw/`: the included public NHANES files used for the current sex and CRP preparations.
- `data_dictionary.csv`: variables, data types, missing counts and descriptions for the frozen cohort. Core analysis variables and derived audit fields are described; original NHANES source fields refer to the corresponding cycle-specific CDC codebooks for complete response labels and skip patterns.
- `reference/`: frozen numerical results, 14,000 bootstrap records, official-file list, original CRP and sex-category comparison inputs, baseline decision curve, Figure 1 historical PDF and HPV figure PDFs.
- `simulation/output/raw/`: original 24 RDS records; `simulation/output/summaries/` and `simulation/output/tables/` contain associated summaries. The obsolete historical manifest was deliberately not included.
- `hpv_sex_covariates_audit_20261102/full/`: initially seeded with frozen performance and decision curve CSVs so figure generation does not require a full bootstrap; running `--full-bootstrap` replaces these with fresh outputs.
- `hpv_sex_covariates_audit_20261102/pilot_10/`, `hpv_sync_20261102/`, `clr_unit_audit_20261103/` and `alcohol_four_level_20261102/models/`: generated analysis outputs; retain the `reference/` copy for comparisons. To keep the distribution compact, 34 pre-generated model RDS files are omitted and are recreated by the driver. Their paths are listed in `reference/regenerable_model_objects_omitted.csv`. The 24 original simulation RDS and frozen bootstrap replicates remain included.

## Analysis definitions and limitations

The primary endpoint is overall HPV infection. Constraint models use 10,230 participants (4,496 HPV-positive), with 8 named ratios but 7 unique constraints because LMR and MLR are reciprocal. The secondary exploratory mortality models use 10,179 participants and 262 events; they use a Cox working model with baseline hazards stratified by the four-level marital variable. The marital not-recorded stratum has 211 participants and no deaths. Violation of proportional hazards and incomplete covariate data limit interpretation. Do not interpret a non-significant Holm-adjusted test as proof of the constraint.

HPV performance uses 8,781 participants and 3,935 positives, additionally adjusting for leisure activity, age at first sex and lifetime male partners. Comparisons and decision curves are apparent fitting-sample performance, not externally validated or optimism-corrected. The sexual-behavior wording differs by survey period; missing values are not treated as valid responses.

Alcohol is four-level: ALQ110=2 denotes fewer than 12 lifetime drinks (not necessarily abstinence); ALQ110=1 and a valid ALQ130=1 is Mild; ALQ110=1 and valid ALQ130>1 is Heavy; remaining cases are Not recorded. ALQ130 is drinks per drinking day; the final category includes structural skips and is not simply a nonresponse group.

For CLR, CRP from 2003–2010 is multiplied by 10 to convert mg/dL to mg/L; 2015–2016 hsCRP is already mg/L. Assay methods and detection limits differ. The two unmeasured 2011–2014 cycles are excluded; only three missing 2015–2016 concentrations are imputed. The early-period analysis is post hoc and cannot demonstrate assay equivalence.

Simulation and the HPV performance outputs do not establish clinical utility. Mortality performance, NHANES III and the old restricted cubic spline results were not adapted and are excluded. The historical Figure 1 PDF is supplied as a conceptual artifact; the complete original drawing source was not recovered, so the figure cannot be regenerated from source here. Recoverable covariate rules are included and checked against the frozen records, but their raw-source joins and the original cohort-selection code are not yet packaged as an independently executable end-to-end chain, and the complete response labels and skip patterns for raw NHANES fields require the cycle-specific official codebooks. These limitations apply to the publicly available package; source-to-cohort reproduction is not complete.
