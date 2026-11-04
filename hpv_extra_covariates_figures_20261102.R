#!/usr/bin/env Rscript
root<-'.'
s<-readLines(file.path(root,'hpv_sync_dca_figure_20261102.R'),warn=FALSE)
old<-"root<-'./hpv_sync_20261102'"
i<-which(s==old);stopifnot(length(i)==1)
s[i]<-"root<-'./hpv_sex_covariates_audit_20261102'"
eval(parse(text=s),envir=.GlobalEnv)
