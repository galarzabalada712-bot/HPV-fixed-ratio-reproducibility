#!/usr/bin/env Rscript
# Apply the same prespecified stratified Cox specification to all ratios.
root<-getwd()
lines<-readLines(file.path(root,"marital_all_models_audit.R"),warn=FALSE)
subline<-function(old,new) {
 idx<-which(lines==old);stopifnot(length(idx)==1)
 lines[idx]<<-new
}
subline('out <- file.path(root,"marital_all_models_audit")',
        'out <- file.path(root,"marital_stratified_all")')
subline('  f0<-as.formula(paste(response,"~",paste(c("r",factors),collapse=" + ")))',
        paste('  f0<-as.formula(paste(response,"~",paste(c("r",factors),collapse=" + ")))',
              '  if(outcome=="Mortality") {f<-update(f,.~. -marital+strata(marital)); f0<-update(f0,.~. -marital+strata(marital))}',sep='\n'))
subline('    unreliable<-outcome=="Mortality" && any(table(z$marital,z$MORTSTAT)[,"1"]==0)',
        '    unreliable<-outcome=="Mortality" && length(warnings)>0')
subline('    # Zero-event level implies separation; save object for audit, never report inferential P.',
        '    # Stratified baseline avoids a coefficient for the zero-event level; retain warnings.')
subline('    if(unreliable) {result$status<-"unreliable_zero_event_category"} else {',
        '    if(unreliable) {result$status<-"unreliable_fitting_warning"} else {')
subline('write.csv(r,file.path(out,"model_status.csv"),row.names=FALSE,na="")',
        paste('ii<-which(r$outcome=="Mortality" & r$ratio!="MLR")',
              'if(all(r$status[ii]=="fitted")){r$Holm_7[ii]<-p.adjust(r$p[ii],"holm",n=7);r$Holm_7[r$outcome=="Mortality" & r$ratio=="MLR"]<-r$Holm_7[r$outcome=="Mortality" & r$ratio=="LMR"]}',
              'write.csv(r,file.path(out,"model_status.csv"),row.names=FALSE,na="")',sep='\n'))
eval(parse(text=lines),envir=.GlobalEnv)
