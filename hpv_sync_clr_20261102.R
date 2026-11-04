#!/usr/bin/env Rscript
# CLR sensitivity in the five survey cycles with observed CRP. Results are
# exploratory; the two wholly unmeasured cycles cannot supply within-cycle donors.
suppressPackageStartupMessages({library(survey);library(mice)})
options(survey.lonely.psu='adjust')
root <- '.'
base <- file.path(root,'alcohol_four_level_20261102')
out <- file.path(root,'hpv_sync_20261102')
d <- read.csv(file.path(base,'analysis_input.csv'))
fixed <- read.csv(file.path(base,'models/model_membership.csv'))
crp <- read.csv(file.path(out,'fixed_HPV_sample_preflight_by_person.csv'))
fr <- read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))
ids <- fixed$SEQN[fixed$outcome=='HPV' & fixed$ratio=='PLR']
stopifnot(length(ids)==10230,!anyDuplicated(ids),!anyDuplicated(crp$SEQN),!anyDuplicated(fr$SEQN))
ii <- match(fr$SEQN,d$SEQN)
dat <- cbind(fr,d[ii,setdiff(names(d),names(fr)),drop=FALSE])
dat$CRP_mg_L <- crp$CRP_mg_L[match(dat$SEQN,crp$SEQN)]
dat$in_fixed <- dat$SEQN %in% ids
dat$measured_cycle <- dat$cycle_start %in% c(2003,2005,2007,2009,2015)
cv <- c('age','RIDRETH1','education','PIR','smokingnew','Alcoholuse','hypertentionnew','diabetesnew','BMXBMI','marital')
for(v in setdiff(cv,'BMXBMI'))dat[[v]] <- factor(dat[[v]])
dat$r <- log(dat$CRP_mg_L)-log(dat$LYMPHOCYTE)
dat$s <- log(dat$CRP_mg_L)+log(dat$LYMPHOCYTE)
design <- svydesign(ids=~SDMVPSU,strata=~design_stratum,weights=~WTfinal,nest=TRUE,data=dat)
keep <- dat$in_fixed & dat$measured_cycle & is.finite(dat$r) & is.finite(dat$s)
stopifnot(sum(keep)==7089)
sub <- design[keep,];sub$variables <- droplevels(sub$variables)
free <- as.formula(paste('HPVnew ~ r + s +',paste(cv,collapse=' + ')))
fixedform <- as.formula(paste('HPVnew ~ r +',paste(cv,collapse=' + ')))
messages <- character()
fits <- withCallingHandlers({list(free=svyglm(free,design=sub,family=quasibinomial()),
  constrained=svyglm(fixedform,design=sub,family=quasibinomial()))},
  warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
stopifnot(length(messages)==0,identical(fits$free$survey.design$variables$SEQN,fits$constrained$survey.design$variables$SEQN))
tt <- regTermTest(fits$free,~s,method='LRT');ww <- regTermTest(fits$free,~s,method='Wald')
est <- 2*coef(fits$free)['s'];se <- 2*sqrt(vcov(fits$free)['s','s']); crit<-qt(.975,ww$ddf)
res <- data.frame(analysis='complete_case_measured_cycles',N=sum(keep),events=sum(sub$variables$HPVnew==1),
  design_df=degf(sub),D=unname(est),se=unname(se),CI_low=unname(est-crit*se),CI_high=unname(est+crit*se),
  working_LRT_p=as.numeric(tt$p),warning=paste(messages,collapse=' | '))
write.csv(res,file.path(out,'CLR_complete_case_results.csv'),row.names=FALSE)
write.csv(data.frame(SEQN=sub$variables$SEQN),file.path(out,'CLR_complete_case_membership.csv'),row.names=FALSE)
saveRDS(list(models=fits,lrt=tt,wald=ww),file.path(out,'CLR_complete_case_models.rds'))
cat('CLR complete case:');print(res)
# Impute only the three missing values in the cycle where hsCRP was actually
# measured. There are no within-cycle donors in 2011–2014, which remain out of
# this measured-cycle sensitivity analysis.
mi_domain <- dat$in_fixed & dat$measured_cycle
stopifnot(sum(mi_domain)==7092,sum(mi_domain & is.na(dat$CRP_mg_L))==3,
          all(dat$cycle_start[mi_domain & is.na(dat$CRP_mg_L)]==2015))
z <- dat[mi_domain & dat$cycle_start==2015,]
pred <- data.frame(logCRP=log(z$CRP_mg_L),HPVnew=z$HPVnew,logLYM=log(z$LYMPHOCYTE),
                   logWeight=log(z$WTfinal),stratum=factor(z$design_stratum),PSU=factor(z$SDMVPSU),
                   z[,cv,drop=FALSE],check.names=FALSE)
pred <- droplevels(pred)
stopifnot(sum(is.na(pred$logCRP))==3,all(complete.cases(pred[,-1])))
method <- mice::make.method(pred);method[] <- '';method['logCRP'] <- 'pmm'
predmat <- mice::make.predictorMatrix(pred);predmat[,] <- 0
predmat['logCRP',setdiff(names(pred),'logCRP')] <- 1
set.seed(20261102)
miwarn <- character()
imputed <- withCallingHandlers(mice(pred,m=20,maxit=5,method=method,predictorMatrix=predmat,printFlag=FALSE),
  warning=function(w){miwarn<<-c(miwarn,conditionMessage(w));invokeRestart('muffleWarning')})
stopifnot(all(vapply(seq_len(20),function(i)all(is.finite(complete(imputed,i)$logCRP)),logical(1))))
qs<-us<-numeric(20);mi_events<-integer(20);impute_records<-list()
for(i in seq_len(20)) {
  zz<-dat
  zz$CRP_mg_L[match(z$SEQN,zz$SEQN)]<-exp(complete(imputed,i)$logCRP)
  zz$r<-log(zz$CRP_mg_L)-log(zz$LYMPHOCYTE);zz$s<-log(zz$CRP_mg_L)+log(zz$LYMPHOCYTE)
  desi<-svydesign(ids=~SDMVPSU,strata=~design_stratum,weights=~WTfinal,nest=TRUE,data=zz)[mi_domain,]
  desi$variables<-droplevels(desi$variables)
  f<-withCallingHandlers(svyglm(free,design=desi,family=quasibinomial()),
     warning=function(w){miwarn<<-c(miwarn,conditionMessage(w));invokeRestart('muffleWarning')})
  stopifnot(length(coef(f))>0,all(is.finite(coef(f))),all(is.finite(vcov(f))))
  qs[i]<-2*coef(f)['s'];us[i]<-4*vcov(f)['s','s'];mi_events[i]<-sum(desi$variables$HPVnew==1)
  impute_records[[i]]<-data.frame(imputation=i,SEQN=z$SEQN[is.na(z$CRP_mg_L)],
                                  imputed_CRP_mg_L=exp(complete(imputed,i)$logCRP[is.na(z$CRP_mg_L)]))
}
pooled<-mice::pool.scalar(qs,us,n=degf(desi)+1,k=length(coef(f)))
pm<-2*pt(-abs(pooled$qbar/sqrt(pooled$t)),df=pooled$df);ci<-qt(.975,df=pooled$df)*sqrt(pooled$t)
resmi<-data.frame(analysis='MI_20_within_2015_measured_cycles_only',N=sum(mi_domain),events=mi_events[1],
  missing_imputed=3,D=pooled$qbar,se=sqrt(pooled$t),df=pooled$df,CI_low=pooled$qbar-ci,CI_high=pooled$qbar+ci,
  Wald_p=pm,warning=paste(unique(miwarn),collapse=' | '))
write.csv(resmi,file.path(out,'CLR_MI_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,impute_records),file.path(out,'CLR_MI_imputed_values.csv'),row.names=FALSE)
write.csv(data.frame(imputation=seq_along(qs),D=qs,variance=us),file.path(out,'CLR_MI_estimates.csv'),row.names=FALSE)
saveRDS(list(imputation=imputed,pooled=pooled,warning=miwarn,predictors=predmat,seed=20261102),file.path(out,'CLR_MI_diagnostics.rds'))
cat('CLR measured-cycle MI:');print(resmi)
capture.output(sessionInfo(),file=file.path(out,'CLR_sessionInfo.txt'))
