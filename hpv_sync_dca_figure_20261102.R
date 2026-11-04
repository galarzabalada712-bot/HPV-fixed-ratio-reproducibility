#!/usr/bin/env Rscript
# Figure contract: compare apparent in-sample net benefit of paired component
# and composite models, with three reference strategies in a 14-panel grid.
# Source: full fixed-sample models; common 0.05–0.80 threshold and survey weights.
suppressPackageStartupMessages({library(ggplot2);library(tidyr);library(svglite)})
root<-'./hpv_sync_20261102'
x<-read.csv(file.path(root,'full/HPV_decision_curves.csv'))
b<-read.csv(file.path(root,'HPV_dca_covariates_only.csv'))
stopifnot(length(unique(x$marker))==14,nrow(x)==14*76,nrow(b)==76)
x<-merge(x,b,by='threshold',all.x=TRUE,sort=FALSE);stopifnot(!anyNA(x$covariates_only))
names(x)[names(x)=='treat_all']<-'Treat all';names(x)[names(x)=='treat_none']<-'Treat none'
names(x)[names(x)=='composite']<-'Composite + covariates';names(x)[names(x)=='components']<-'Components + covariates'
names(x)[names(x)=='covariates_only']<-'Covariates only'
long<-pivot_longer(x,-c(marker,threshold),names_to='model',values_to='net_benefit')
long$marker<-factor(long$marker,levels=c('RAR','NPAR','MAR','NAR','PNI','NLR','PLR','LMR','MLR','SII','SIRI','PPN','PIV','HALP'))
palette<-c('Composite + covariates'='#20639B','Components + covariates'='#CE7E00','Covariates only'='#555555','Treat all'='#31846A','Treat none'='#A54C4C')
p<-ggplot(long,aes(threshold,net_benefit,colour=model,linetype=model))+
 geom_line(linewidth=.45)+facet_wrap(~marker,ncol=4)+coord_cartesian(ylim=c(-.05,.45))+
 scale_colour_manual(values=palette)+scale_linetype_manual(values=c('solid','dashed','dotted','longdash','dotdash'))+
 labs(x='Threshold probability',y='Survey-weighted net benefit',colour=NULL,linetype=NULL)+
 theme_bw(base_size=9)+theme(legend.position='bottom',legend.box='horizontal',panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey95'))+
 guides(colour=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE))
ggsave(file.path(root,'full/HPV_DCA_candidate.svg'),p,width=11.5,height=10,device=svglite::svglite)
ggsave(file.path(root,'full/HPV_DCA_candidate.pdf'),p,width=11.5,height=10,device=cairo_pdf)
ggsave(file.path(root,'full/HPV_DCA_candidate.png'),p,width=11.5,height=10,dpi=180)
write.csv(long,file.path(root,'full/HPV_DCA_plot_source.csv'),row.names=FALSE)
performance<-read.csv(file.path(root,'full/HPV_performance.csv'))
stopifnot(nrow(performance)==14,all(performance$B==1000),all(performance$B_valid>=900))
performance$marker<-factor(performance$marker,levels=rev(c('RAR','NPAR','MAR','NAR','PNI','NLR','PLR','LMR','MLR','SII','SIRI','PPN','PIV','HALP')))
p2<-ggplot(performance,aes(x=delta_auc,y=marker))+
 geom_vline(xintercept=0,colour='grey55',linetype='dashed',linewidth=.4)+
 geom_errorbar(aes(xmin=delta_auc_low,xmax=delta_auc_high),width=.16,orientation='y',colour='#20639B')+
 geom_point(size=1.8,colour='#20639B')+theme_bw(base_size=10)+
 labs(x='Apparent difference in AUC (components minus composite)',y=NULL)+theme(panel.grid.minor=element_blank())
ggsave(file.path(root,'full/HPV_delta_AUC_candidate.svg'),p2,width=6.5,height=6.5,device=svglite::svglite)
ggsave(file.path(root,'full/HPV_delta_AUC_candidate.pdf'),p2,width=6.5,height=6.5,device=cairo_pdf)
ggsave(file.path(root,'full/HPV_delta_AUC_candidate.png'),p2,width=6.5,height=6.5,dpi=180)
cat('Candidate figure produced; visual review and manuscript-specific checks remain.\n')
