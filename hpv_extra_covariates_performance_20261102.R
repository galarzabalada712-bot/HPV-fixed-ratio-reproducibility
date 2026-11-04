#!/usr/bin/env Rscript
# Reuse the verified performance engine with the user-selected extra covariates.
root<-'.'
s<-readLines(file.path(root,'hpv_sync_performance_20261102.R'),warn=FALSE)
replace_exact<-function(old,new){i<-which(s==old);stopifnot(length(i)==1);s[i]<<-new}
replace_exact("root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sync_20261102')",
              "root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sex_covariates_audit_20261102')")
replace_exact("d<-read.csv(file.path(base,'analysis_input.csv'));frame<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))",
              "d<-read.csv(file.path(base,'analysis_input.csv'));sex<-read.csv(file.path(out,'model_covariates.csv'));stopifnot(!anyDuplicated(sex$SEQN));d<-merge(d,sex,by='SEQN',all.x=TRUE,sort=FALSE);frame<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))")
replace_exact("cv<-c('age','RIDRETH1','education','PIR','smokingnew','Alcoholuse','hypertentionnew','diabetesnew','BMXBMI','marital')",
              "cv<-c('age','RIDRETH1','education','PIR','smokingnew','Alcoholuse','hypertentionnew','diabetesnew','BMXBMI','marital','activity','sexagenew','sextotal')")
replace_exact("keep<-dat$in_domain;stopifnot(sum(keep)==10230,setequal(dat$SEQN[keep],ids))",
              "keep<-dat$in_domain & complete.cases(dat[,c('activity','sexagenew','sextotal')]);stopifnot(sum(keep)==8781,sum(dat$HPVnew[keep])==3935,all(dat$SEQN[keep] %in% ids))")
replace_exact("y<-z$HPVnew;w<-weights(di,type='sampling');stopifnot(sum(y)==4496,all(levels(z$Alcoholuse)==as.character(0:3)))",
              "y<-z$HPVnew;w<-weights(di,type='sampling');stopifnot(length(y)==8781,sum(y)==3935,all(levels(z$Alcoholuse)==as.character(0:3)))")
eval(parse(text=s),envir=.GlobalEnv)
