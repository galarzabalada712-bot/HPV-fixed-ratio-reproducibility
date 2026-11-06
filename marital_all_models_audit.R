#!/usr/bin/env Rscript
# Post-hoc screenshot coding. Eight HPV fits and eight Cox estimability checks.
suppressPackageStartupMessages({library(survey);library(survival)})
options(survey.lonely.psu="adjust")
root <- getwd()
base <- file.path(root,"alcohol_missing_reanalysis")
out <- file.path(root,"marital_all_models_audit")
dir.create(out,showWarnings=FALSE)
d <- read.csv(file.path(base,"analysis_input.csv"))
# Post-hoc amendment: examination is the default; interview is a sensitivity analysis.
origin<-Sys.getenv('MORTALITY_ORIGIN',unset='exam');stopifnot(origin %in% c('exam','interview'))
months<-if(origin=='exam') d$PERMTH_EXM else d$PERMTH_INT
stopifnot(all(is.na(months) | months>0))
d$FOLLOWUP_YEARS <- months/12
frame <- read.csv(file.path(base,"survey_design_frame.csv"))
stopifnot(nrow(d)==11884,!anyDuplicated(d$SEQN),!anyDuplicated(frame$SEQN))
raw <- d$DMDMARTL
d$marital <- ifelse(raw %in% 5,1,ifelse(raw %in% c(1,6),2,ifelse(raw %in% c(2,3,4),3,4)))
stopifnot(all(d$marital[raw %in% 5]==1),all(d$marital[raw %in% c(1,6)]==2),
          all(d$marital[raw %in% c(2,3,4)]==3),all(d$marital[is.na(raw)|raw %in% c(77,99)]==4))
d$marital_source <- ifelse(is.na(raw) & d$RIDAGEYR<20 & d$cycle_start>=2007,
                           "structural_not_public",ifelse(raw %in% c(77,99),"refused_or_unknown",
                           ifelse(is.na(raw),"other_missing","observed")))
write.csv(d[,c("SEQN","RIDAGEYR","cycle_start","DMDMARTL","marital","marital_source")],
          file.path(out,"marital_classification.csv"),row.names=FALSE,na="")
extra<-setdiff(names(d),names(frame)); dat<-cbind(frame,d[match(frame$SEQN,d$SEQN),extra,drop=FALSE])
rownames(dat)<-NULL; dat$in_cohort<-dat$SEQN %in% d$SEQN
dat$mortality_eligible<-as.character(dat$mortality_eligible)%in%c("True","TRUE")
covs<-c("age","RIDRETH1","education","PIR","smokingnew","Alcoholuse","hypertentionnew","diabetesnew","BMXBMI","marital")
for(v in setdiff(c(covs,"activity"),"BMXBMI")) dat[[v]]<-factor(dat[[v]])
ratios<-data.frame(name=c("RAR","NPAR","MAR","NAR","NLR","PLR","LMR","MLR"),
 A=c("RDW","NEU_PCT","MONOCYTE","NEUTROPHIL","NEUTROPHIL","PLATELET","LYMPHOCYTE","MONOCYTE"),
 B=c("ALBUMIN","ALBUMIN","ALBUMIN","ALBUMIN","LYMPHOCYTE","LYMPHOCYTE","MONOCYTE","LYMPHOCYTE"))
results<-list();groups<-list();members<-list()
for(outcome in c("HPV","Mortality")) for(j in seq_len(nrow(ratios))) {
  rr<-ratios[j,]
  a<-dat[[rr$A]];b<-dat[[rr$B]]
  dat$r<-log(ifelse(a>0,a,NA))-log(ifelse(b>0,b,NA));dat$s<-log(ifelse(a>0,a,NA))+log(ifelse(b>0,b,NA))
  factors<-if(outcome=="HPV") covs else c(covs,"activity")
  event<-if(outcome=="HPV") "HPVnew" else "MORTSTAT"
  keep<-dat$in_cohort & complete.cases(dat[,c("r","s",event,factors)])
  if(outcome=="Mortality") keep<-keep & dat$mortality_eligible
  stopifnot(!anyNA(keep))
  design<-svydesign(ids=~SDMVPSU,strata=~design_stratum,weights=~WTfinal,nest=TRUE,data=dat)[keep,]
  design$variables<-droplevels(design$variables);z<-design$variables
  members[[length(members)+1]]<-data.frame(outcome=outcome,ratio=rr$name,SEQN=z$SEQN)
  for(i in 1:4){q<-z[z$marital==i,,drop=FALSE];groups[[length(groups)+1]]<-data.frame(outcome=outcome,ratio=rr$name,marital=i,n=nrow(q),events=sum(q[[event]]==1),structural=sum(q$marital_source=="structural_not_public"),refused_unknown=sum(q$marital_source=="refused_or_unknown"))}
  response<-if(outcome=="HPV") "HPVnew" else "Surv(FOLLOWUP_YEARS,MORTSTAT)"
  f<-as.formula(paste(response,"~",paste(c("r","s",factors),collapse=" + ")))
  f0<-as.formula(paste(response,"~",paste(c("r",factors),collapse=" + ")))
  warnings<-character(); result<-data.frame(outcome=outcome,ratio=rr$name,n=nrow(z),events=sum(z[[event]]==1),design_df=degf(design),D=NA_real_,se=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p=NA_real_,status="failed",warning="",error="")
  tryCatch(withCallingHandlers({
    fit<-if(outcome=="HPV") svyglm(f,design=design,family=quasibinomial()) else svycoxph(f,design=design)
    fit0<-if(outcome=="HPV") svyglm(f0,design=design,family=quasibinomial()) else svycoxph(f0,design=design)
    stopifnot(nrow(model.frame(fit))==nrow(z),nrow(model.frame(fit0))==nrow(z))
    # Zero-event level implies separation; save object for audit, never report inferential P.
    unreliable<-outcome=="Mortality" && any(table(z$marital,z$MORTSTAT)[,"1"]==0)
    saveRDS(list(free=fit,constrained=fit0,warning=warnings,formula=f,formula0=f0,SEQN=z$SEQN),file.path(out,paste0(outcome,"_",rr$name,".rds")))
    if(unreliable) {result$status<-"unreliable_zero_event_category"} else {
      stopifnot(all(is.finite(coef(fit))),all(is.finite(vcov(fit))),all(is.finite(coef(fit0))))
      test<-regTermTest(fit,~s,method="LRT"); w<-regTermTest(fit,~s,method="Wald")
      result$D<-2*coef(fit)[["s"]];result$se<-2*sqrt(vcov(fit)["s","s"])
      crit<-qt(.975,w$ddf);result$CI_low<-result$D-crit*result$se;result$CI_high<-result$D+crit*result$se
      result$p<-as.numeric(test$p);result$status<-"fitted"
      saveRDS(list(lrt=test,wald=w),file.path(out,paste0(outcome,"_",rr$name,"_tests.rds")))
    }
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e){result$error<<-conditionMessage(e)})
  result$warning<-paste(unique(warnings),collapse=" | ")
  if(result$status=="fitted" && nzchar(result$warning))result$status<-"fitted_with_warning"
  results[[length(results)+1]]<-result
  cat(outcome,rr$name,result$status,"N",nrow(z),"events",result$events,"warning",result$warning,"error",result$error,"\n")
}
r<-do.call(rbind,results);r$Holm_7<-NA_real_;i<-which(r$outcome=="HPV" & r$ratio!="MLR")
if(all(r$status[i]=="fitted")){r$Holm_7[i]<-p.adjust(r$p[i],"holm",n=7);r$Holm_7[r$ratio=="MLR" & r$outcome=="HPV"]<-r$Holm_7[r$ratio=="LMR" & r$outcome=="HPV"]}
write.csv(r,file.path(out,"model_status.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,groups),file.path(out,"marital_group_support.csv"),row.names=FALSE)
write.csv(do.call(rbind,members),file.path(out,"model_membership.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
stopifnot(nrow(r)==16)
