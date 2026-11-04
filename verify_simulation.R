#!/usr/bin/env Rscript
root<-'.';out<-file.path(root,'verification')
dir.create(out,showWarnings=FALSE)
sim<-file.path(root,'simulation/output');paths<-sort(list.files(file.path(sim,'raw'),pattern='\\.rds$',full.names=TRUE))
inf<-read.csv(file.path(sim,'summaries/summary_inferential.csv'));pred<-read.csv(file.path(sim,'summaries/summary_predictive.csv'))
stopifnot(length(paths)==24,nrow(inf)==24,nrow(pred)==24,!anyDuplicated(inf$scenario_label))
results<-lapply(paths,function(path){x<-readRDS(path);r<-x$results;s<-x$scenario;label<-as.character(s$scenario_label[1]);a<-inf[inf$scenario_label==label,];b<-pred[pred$scenario_label==label,]
 stopifnot(nrow(r)==1000,nrow(a)==1,nrow(b)==1,all(r$fit_success==1),all(r$ratio_converged==1),all(r$component_converged==1),all(r$covariate_converged==1))
 data.frame(scenario_label=label,n_repetitions=nrow(r),failed=sum(r$fit_success!=1),warning_count=sum(nzchar(r$warning_message)),
    rejection=mean(r$reject_lrt),rejection_diff=mean(r$reject_lrt)-a$rejection_rate,
    bias=mean(r$D_hat-r$D_true),bias_diff=mean(r$D_hat-r$D_true)-a$bias_D,
    coverage=mean(r$cover_D),coverage_diff=mean(r$cover_D)-a$coverage,
    delta_auc=mean(r$delta_auc),delta_auc_diff=mean(r$delta_auc)-b$delta_auc,
    test_seed=x$test_seed,n_test=s$n_test[1],rho=s$rho[1],D=s$D_true[1])})
z<-do.call(rbind,results);stopifnot(max(abs(z$rejection_diff),abs(z$bias_diff),abs(z$coverage_diff),abs(z$delta_auc_diff))<1e-10,length(unique(z$test_seed))==24)
write.csv(z,file.path(out,'simulation_raw_reaggregation_check.csv'),row.names=FALSE)
cat('PASS: 24 x 1000 successful fits, original RDS vs summaries; 24 scenario-specific test seeds. Total warning-bearing repetitions:',sum(z$warning_count),'\n')

