#!/usr/bin/env Rscript
# HPV performance with a fixed primary sample and paired survey bootstrap.
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust')
root<-'.';base<-file.path(root,'alcohol_four_level_20261102');out<-file.path(root,'hpv_sync_20261102')
args<-commandArgs(trailingOnly=TRUE);B<-if(length(args))as.integer(args[1]) else 1000L
workers<-if(length(args)>1)as.integer(args[2]) else 1L
stopifnot(B>=10,workers>=1)
tag<-if(B==1000L)'full' else paste0('pilot_',B)
dir.create(file.path(out,tag),showWarnings=FALSE)
d<-read.csv(file.path(base,'analysis_input.csv'));frame<-read.csv(file.path(root,'alcohol_missing_reanalysis/survey_design_frame.csv'))
membership<-read.csv(file.path(base,'models/model_membership.csv'))
ids<-membership$SEQN[membership$outcome=='HPV' & membership$ratio=='PLR']
ix<-match(frame$SEQN,d$SEQN);dat<-cbind(frame,d[ix,setdiff(names(d),names(frame)),drop=FALSE]);dat$in_domain<-dat$SEQN %in% ids
cv<-c('age','RIDRETH1','education','PIR','smokingnew','Alcoholuse','hypertentionnew','diabetesnew','BMXBMI','marital')
for(v in setdiff(cv,'BMXBMI'))dat[[v]]<-factor(dat[[v]])
dat$HEMOGLOBIN<-dat$LBXHGB
with(dat,{stopifnot(all(is.finite(RDW[in_domain])),all(is.finite(HEMOGLOBIN[in_domain])))})
dat$RAR<-dat$RDW/dat$ALBUMIN;dat$NPAR<-dat$NEU_PCT/dat$ALBUMIN
dat$MAR<-dat$MONOCYTE/dat$ALBUMIN;dat$NAR<-dat$NEUTROPHIL/dat$ALBUMIN
dat$PNI<-10*dat$ALBUMIN+5*dat$LYMPHOCYTE;dat$NLR<-dat$NEUTROPHIL/dat$LYMPHOCYTE
dat$PLR<-dat$PLATELET/dat$LYMPHOCYTE;dat$LMR<-dat$LYMPHOCYTE/dat$MONOCYTE;dat$MLR<-dat$MONOCYTE/dat$LYMPHOCYTE
dat$SII<-dat$PLATELET*dat$NEUTROPHIL/dat$LYMPHOCYTE;dat$SIRI<-dat$NEUTROPHIL*dat$MONOCYTE/dat$LYMPHOCYTE
dat$PPN<-dat$PLATELET*dat$NEUTROPHIL;dat$PIV<-dat$PLATELET*dat$NEUTROPHIL*dat$MONOCYTE/dat$LYMPHOCYTE
dat$HALP<-dat$HEMOGLOBIN*dat$ALBUMIN*dat$LYMPHOCYTE/dat$PLATELET
markers<-list(RAR=c('RDW','ALBUMIN'),NPAR=c('NEU_PCT','ALBUMIN'),MAR=c('MONOCYTE','ALBUMIN'),
 NAR=c('NEUTROPHIL','ALBUMIN'),PNI=c('ALBUMIN','LYMPHOCYTE'),NLR=c('NEUTROPHIL','LYMPHOCYTE'),
 PLR=c('PLATELET','LYMPHOCYTE'),LMR=c('LYMPHOCYTE','MONOCYTE'),MLR=c('MONOCYTE','LYMPHOCYTE'),
 SII=c('PLATELET','NEUTROPHIL','LYMPHOCYTE'),SIRI=c('NEUTROPHIL','MONOCYTE','LYMPHOCYTE'),
 PPN=c('PLATELET','NEUTROPHIL'),PIV=c('PLATELET','NEUTROPHIL','MONOCYTE','LYMPHOCYTE'),
 HALP=c('HEMOGLOBIN','ALBUMIN','LYMPHOCYTE','PLATELET'))
keep<-dat$in_domain;stopifnot(sum(keep)==10230,setequal(dat$SEQN[keep],ids))
for(v in unique(c(names(markers),unlist(markers))))stopifnot(all(is.finite(dat[[v]][keep])),all(dat[[v]][keep]>0))
design<-svydesign(ids=~SDMVPSU,strata=~design_stratum,weights=~WTfinal,nest=TRUE,data=dat)
RNGkind("L'Ecuyer-CMRG");set.seed(20261102)
# Generate full-design replicates in bounded batches, then retain domain rows.
# This preserves full sampling strata/PSUs without allocating a 68k x B matrix.
repweights<-matrix(NA_real_,sum(keep),B)
for(start in seq.int(1L,B,by=25L)) {
  end<-min(B,start+24L)
  repdesign<-as.svrepdesign(design,type='bootstrap',replicates=end-start+1L)
  repweights[,start:end]<-weights(repdesign[keep,],type='analysis')
  rm(repdesign);gc(FALSE)
}
di<-design[keep,];z<-droplevels(di$variables)
y<-z$HPVnew;w<-weights(di,type='sampling');stopifnot(sum(y)==4496,all(levels(z$Alcoholuse)==as.character(0:3)))
auc<-function(y,p,w){o<-order(p);y<-y[o];p<-p[o];w<-w[o];rp<-rle(p);en<-cumsum(rp$lengths);nw<-cumsum(w*(y==0))[en];before<-c(0,head(nw,-1));ties<-nw-before;sum(w*(y==1)*(rep(before,rp$lengths)+.5*rep(ties,rp$lengths)))/(sum(w[y==1])*sum(w[y==0]))}
metrics<-function(y,a,b,w){ev<-y==1;up<-b>a;dn<-b<a;wm<-function(x,k)sum(w[k]*x[k])/sum(w[k]);
 c(auc_a=auc(y,a,w),auc_b=auc(y,b,w),nri=wm(up,ev)-wm(dn,ev)+wm(dn,!ev)-wm(up,!ev),idi=(wm(b,ev)-wm(a,ev))-(wm(b,!ev)-wm(a,!ev)))}
calc<-function(nm){
 terms<-markers[[nm]]; composite<-if(nm=='PNI')nm else paste0('log(',nm,')');components<-if(nm=='PNI')terms else paste0('log(',terms,')')
 f0<-as.formula(paste('HPVnew ~',paste(c(composite,cv),collapse=' + ')));f1<-as.formula(paste('HPVnew ~',paste(c(components,cv),collapse=' + ')))
 X0<-model.matrix(f0,z);X1<-model.matrix(f1,z);Y<-z$HPVnew
 runfit<-function(X,ww){warning<-character();f<-tryCatch(withCallingHandlers(glm.fit(X,Y,weights=ww/mean(ww),family=quasibinomial(),control=glm.control(maxit=100)),warning=function(e){warning<<-c(warning,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)NULL)
 if(is.null(f)||!f$converged||any(!is.finite(f$coefficients))||length(warning))return(NULL);f$fitted.values}
 p0<-runfit(X0,w);p1<-runfit(X1,w);stopifnot(!is.null(p0),!is.null(p1))
 ref<-svyglm(f0,design=di,family=quasibinomial());stopifnot(max(abs(as.numeric(predict(ref,type='response'))-p0))<1e-7)
 pt<-metrics(Y,p0,p1,w);th<-seq(.05,.80,.01);prev<-weighted.mean(Y,w)
 nb<-function(p,t)sum(w*(Y==1)*(p>=t))/sum(w)-sum(w*(Y==0)*(p>=t))/sum(w)*t/(1-t)
 dca<-data.frame(marker=nm,threshold=th,composite=vapply(th,function(t)nb(p0,t),numeric(1)),components=vapply(th,function(t)nb(p1,t),numeric(1)),treat_all=prev-(1-prev)*th/(1-th),treat_none=0)
 boots<-matrix(NA_real_,B,3,dimnames=list(NULL,c('delta_auc','nri','idi')))
 for(i in seq_len(B)){wb<-repweights[,i];a<-runfit(X0,wb);b<-runfit(X1,wb);if(!is.null(a)&&!is.null(b)){v<-metrics(Y,a,b,wb);boots[i,]<-c(v['auc_b']-v['auc_a'],v['nri'],v['idi'])}}
 valid<-complete.cases(boots);if(sum(valid)<.9*B)stop(nm,': fewer than 90% valid replicates')
 vals<-boots[valid,,drop=FALSE];delta<-pt['auc_b']-pt['auc_a'];qs<-apply(vals,2,quantile,probs=c(.025,.975),names=FALSE)
 result<-data.frame(marker=nm,N=length(Y),events=sum(Y),auc_composite=pt['auc_a'],auc_components=pt['auc_b'],delta_auc=delta,
 delta_auc_low=qs[1,1],delta_auc_high=qs[2,1],nri=pt['nri'],nri_low=qs[1,2],nri_high=qs[2,2],idi=pt['idi'],idi_low=qs[1,3],idi_high=qs[2,3],B=B,B_valid=sum(valid),stringsAsFactors=FALSE)
 cat('Completed',nm,'valid replicates',sum(valid),'of',B,'\n');flush.console()
 list(result=result,dca=dca,bootstrap=data.frame(marker=nm,replicate=seq_len(B),boots),formulas=c(as.character(f0)[3],as.character(f1)[3]))
}
if(workers>1){cl<-parallel::makeCluster(workers);on.exit(parallel::stopCluster(cl),add=TRUE);parallel::clusterExport(cl,c('markers','cv','z','y','w','repweights','B','auc','metrics','calc','di'),envir=environment());parallel::clusterEvalQ(cl,{library(survey);NULL});result<-parallel::parLapply(cl,names(markers),calc)}else result<-lapply(names(markers),calc)
tab<-do.call(rbind,lapply(result,`[[`,'result'));write.csv(tab,file.path(out,tag,'HPV_performance.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(result,`[[`,'dca')),file.path(out,tag,'HPV_decision_curves.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(result,`[[`,'bootstrap')),file.path(out,tag,'HPV_bootstrap_replicates.csv'),row.names=FALSE)
write.csv(data.frame(SEQN=z$SEQN),file.path(out,tag,'HPV_model_membership.csv'),row.names=FALSE)
saveRDS(list(result=result,seed=20261102,B=B,workers=workers,df=degf(di),covariates=cv,session=sessionInfo()),file.path(out,tag,'HPV_analysis_audit.rds'))
print(tab);cat('PASS: paired weights and fixed sample for 14 models, outputs:',tag,'\n')
