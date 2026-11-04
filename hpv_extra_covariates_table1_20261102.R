#!/usr/bin/env Rscript
root<-'.';s<-readLines(file.path(root,'hpv_sync_table1_20261102.R'),warn=FALSE)
replace_exact<-function(old,new){i<-which(s==old);stopifnot(length(i)==1);s[i]<<-new}
replace_exact("root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sync_20261102')",
              "root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sex_covariates_audit_20261102')")
replace_exact("d<-read.csv(file.path(base,'analysis_input.csv'));fr<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))",
              "d<-read.csv(file.path(base,'analysis_input.csv'));sex<-read.csv(file.path(out,'model_covariates.csv'));d<-merge(d,sex,by='SEQN',all.x=TRUE,sort=FALSE);fr<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))")
replace_exact(" activity=c('No reported leisure activity','Any reported leisure activity'),hypertentionnew=c('No','Yes'),diabetesnew=c('No','Yes'))",
              " activity=c('No reported leisure activity','Any reported leisure activity'),hypertentionnew=c('No','Yes'),diabetesnew=c('No','Yes'),sexagenew=c('Before 18 years','At least 18 years'),sextotal=c('0','1','2–5','6–10','More than 10'))")
replace_exact("levels0<-c('smokingnew','Alcoholuse','activity','hypertentionnew','diabetesnew')",
              "levels0<-c('smokingnew','Alcoholuse','activity','hypertentionnew','diabetesnew','sextotal')")
eval(parse(text=s),envir=.GlobalEnv)
