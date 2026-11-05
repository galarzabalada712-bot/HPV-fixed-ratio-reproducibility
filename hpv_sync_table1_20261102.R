#!/usr/bin/env Rscript
# Descriptive preliminary Table 1 for both the candidate cohort and model domain.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sync_20261102')
d<-read.csv(file.path(base,'analysis_input.csv'));fr<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))
m<-read.csv(file.path(base,'models/model_membership.csv'));ids<-m$SEQN[m$outcome=='HPV' & m$ratio=='PLR']
dat<-cbind(fr,d[match(fr$SEQN,d$SEQN),setdiff(names(d),names(fr)),drop=FALSE]);dat$in_cohort<-dat$SEQN %in% d$SEQN;dat$in_model<-dat$SEQN %in% ids
des<-svydesign(ids=~SDMVPSU,strata=~design_stratum,weights=~WTfinal,nest=TRUE,data=dat)
groups<-list(age=c('18–29','30–39','40–49','50–59'),RIDRETH1=c('Mexican American','Other Hispanic','Non-Hispanic White','Non-Hispanic Black','Other'),
 education=c('Below high school','High school','Above high school'),marital=c('Never married','Married or living with partner','Previously married','Not recorded'),
 PIR=c('≤1.3','>1.3 to 3.5','>3.5'),smokingnew=c('Never','Former','Current'),Alcoholuse=c('Never','Mild','Heavy','Not recorded'),
 activity=c('No reported leisure activity','Any reported leisure activity'),hypertentionnew=c('No','Yes'),diabetesnew=c('No','Yes'))
levels0<-c('smokingnew','Alcoholuse','activity','hypertentionnew','diabetesnew')
outputs<-list()
tests<-list()
for(domain in c('candidate_cohort','HPV_model')) {
 ss<-if(domain=='candidate_cohort')dat$in_cohort else dat$in_model
 s<-des[ss,];z<-s$variables;stopifnot(sum(ss)==if(domain=='candidate_cohort')11884 else 10230)
 s$variables$HPVgroup<-factor(s$variables$HPVnew)
 for(v in names(groups)) {
  values<-as.character(if(v %in% levels0)seq_along(groups[[v]])-1 else seq_along(groups[[v]]))
  z$group<-factor(ifelse(is.na(z[[v]]),'NA',as.character(z[[v]])),levels=c(values,'NA'))
  s$variables$group<-droplevels(z$group)
  # Table 1: education and hypertension tests use known values only; descriptive counts retain missing values.
  test_design<-if(domain=='candidate_cohort' && v %in% c('education','hypertentionnew')) s[!is.na(z[[v]]),] else s
  test_design$variables$group<-droplevels(test_design$variables$group)
  warns<-character();pv<-tryCatch(withCallingHandlers(svychisq(~group+HPVgroup,test_design,statistic='F')$p.value,
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){warns<<-c(warns,conditionMessage(e));NA_real_})
  tests[[length(tests)+1]]<-data.frame(domain=domain,variable=v,method='Rao-Scott design-adjusted F',p=as.numeric(pv),warning=paste(warns,collapse=' | '))
  for(j in seq_along(levels(z$group)))for(h in 0:1) {
   sel<-z$HPVnew==h; lev<-levels(z$group)[j]; n<-sum(sel & z$group==lev)
   pp<-sum(weights(s,type='sampling')[sel & z$group==lev])/sum(weights(s,type='sampling')[sel])
   outputs[[length(outputs)+1]]<-data.frame(domain=domain,variable=v,code=lev,label=if(lev=='NA')'Missing' else groups[[v]][match(lev,values)],HPV=h,n=n,weighted_percent=100*pp)
  }
 }
 for(h in 0:1){sh<-s[s$variables$HPVnew==h,];est<-svymean(~BMXBMI,sh,na.rm=TRUE);outputs[[length(outputs)+1]]<-data.frame(domain=domain,variable='BMXBMI',code='mean',label='Mean BMI (standard error)',HPV=h,n=sum(!is.na(sh$variables$BMXBMI)),weighted_percent=as.numeric(est))
  outputs[[length(outputs)+1]]<-data.frame(domain=domain,variable='BMXBMI',code='SE',label='BMI standard error',HPV=h,n=sum(!is.na(sh$variables$BMXBMI)),weighted_percent=as.numeric(SE(est)))}
 fit<-svyglm(BMXBMI~factor(HPVnew),s);tests[[length(tests)+1]]<-data.frame(domain=domain,variable='BMXBMI',method='survey weighted linear regression',p=summary(fit)$coefficients[2,4],warning='')
}
tab<-do.call(rbind,outputs);write.csv(tab,file.path(out,'Table1_descriptive_preliminary.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,tests),file.path(out,'Table1_design_tests_preliminary.csv'),row.names=FALSE,na='')
cat('PASS: candidate cohort 11884; HPV fixed model 10230; design weighted percentages; labels and missing recorded\n')
