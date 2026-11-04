#!/usr/bin/env Rscript
root<-getwd()
base<-file.path(root,"alcohol_four_level_20261102")
fixed<-read.csv(file.path(base,"fixed_model_membership.csv"))
d<-read.csv(file.path(base,"analysis_input.csv"))
stopifnot(nrow(d)==11884,all(d$Alcoholuse %in% 0:3),!anyNA(d$Alcoholuse))
# Reuse the verified model runner, redirecting both input and output in memory.
script<-readLines(file.path(root,"marital_stratified_all.R"),warn=FALSE)
old<-'"marital_stratified_all"';stopifnot(sum(grepl(old,script,fixed=TRUE))==1)
script<-gsub(old,'"alcohol_four_level_20261102/models"',script,fixed=TRUE)
idx<-which(script=='eval(parse(text=lines),envir=.GlobalEnv)');stopifnot(length(idx)==1)
script[idx]<-paste(
 'subline(\'base <- file.path(root,"alcohol_missing_reanalysis")\',\'base <- file.path(root,"alcohol_four_level_20261102")\')',
 'subline(\'frame <- read.csv(file.path(base,"survey_design_frame.csv"))\',\'frame <- read.csv(file.path(root,"alcohol_missing_reanalysis/survey_design_frame.csv"))\')',
 'eval(parse(text=lines),envir=.GlobalEnv)',sep='\n')
eval(parse(text=script),envir=.GlobalEnv)
base<-file.path(root,"alcohol_four_level_20261102")
actual<-read.csv(file.path(base,"models/model_membership.csv"))
res<-read.csv(file.path(base,"models/model_status.csv"))
stopifnot(nrow(res)==16,nrow(actual)==nrow(fixed))
for(i in seq_len(nrow(res))) {
 row<-res[i,]; a<-actual$SEQN[actual$outcome==row$outcome & actual$ratio==row$ratio]
 b<-fixed$SEQN[fixed$outcome==row$outcome & fixed$ratio==row$ratio]
 stopifnot(setequal(a,b))
 if(row$status %in% c("fitted","fitted_with_warning")) {
  z<-readRDS(file.path(base,"models",paste0(row$outcome,"_",row$ratio,".rds")))
  tt<-readRDS(file.path(base,"models",paste0(row$outcome,"_",row$ratio,"_tests.rds")))
  stopifnot(identical(as.numeric(a),as.numeric(z$SEQN)),
            identical(z$SEQN,z$free$survey.design$variables$SEQN),
            identical(z$SEQN,z$constrained$survey.design$variables$SEQN),
            tt$lrt$df==1,tt$lrt$test.terms=="s")
  if(row$outcome=="Mortality")stopifnot(length(attr(terms(z$free),"specials")$strata)==1)
 }
}
cat("Completed 16 fits; fixed memberships and saved paired model objects checked.\n")
