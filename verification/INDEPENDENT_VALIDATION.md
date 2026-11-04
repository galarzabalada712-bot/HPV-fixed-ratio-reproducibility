# Independent verification record

The package was copied to a separate temporary directory under a different parent, and `run_all.py` was launched by absolute path from outside that copy. It did not require the original project or historical release paths. It completed the following checks with the recorded local R/Python versions:

- Seven-cycle sex categories (11,884) and five-cycle CRP inputs for 10,230 fixed HPV participants matched frozen reference values.
- Diabetes, leisure activity, education, PIR, smoking, hypertension, follow-up time and four-level alcohol each matched all 11,884 cohort records; current marital coding and weight rule also matched.
- Sixteen constraint models, their matched IDs and 7-constraint Holm adjustment matched the frozen numerical results. Table 1, complete-case CLR, 20-imputation CLR, and 2003–2010 CLR sensitivity estimates matched.
- All 14 HPV models completed a 10-replicate paired-bootstrap smoke test; point estimates and decision curves matched frozen full-run values. The saved 14,000 bootstrap results reproduce frozen 1,000-replicate percentile intervals. A pilot's bootstrap values are not expected to equal the first 10 replicates of the full run because batching changes the pseudorandom sequence.
- Existing 24 simulation RDS files each contain 1,000 successful runs and independently reproduce the summary statistics. A separate stage-1 five-repetition simulation pilot passed in a temporary copy.

During testing a direct stage-1 invocation overwrote one same-named simulation RDS in the disposable test copy; the original RDS were restored to that copy before the final successful verification. The included `run_simulation_pilot.py` prevents this by working in a temporary copy. This packaging cycle did not rerun the complete 1,000-replicate HPV performance bootstrap or the complete stage-4 simulation. The original historical source-to-cohort chain and Figure 1 drawing source remain unavailable.
