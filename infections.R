gc();rm(list=ls())
library(ggplot2)
library(data.table)
library(dplyr)
library(lme4)
library(influence.ME)
require(multcomp)

###########
#load data#
###########

#meta analysis for combining intervention effects
get_meta_estimate<-function(summaries){
  estimates<-summaries[,'Estimate']
  std_errors<-summaries[,'Std. Error']
  W<-1/(std_errors^2)  
  meta.est<-sum(W*estimates)/sum(W)
  meta.se<-sqrt(1/sum(W))
  meta.cilow<-meta.est+1.96*meta.se
  meta.cihigh<-meta.est-1.96*meta.se  
  out<-data.table(Estimate=round(meta.est,3),CI95low=round(meta.cilow,3),CI95high=round(meta.cihigh,3))
  return(out)
}


#create new data.table with counts for SABs combined - impacts estimation of intercept
#only need to included relevant variables for model fitting
#also include Start_date just in case need to map back later
HAI.bySite_1<-HAI.bySite[,list(N=sum(N),OBD.adj=unique(OBD.adj)),by=list(Start_Date,Site,Infection_type,Week_no,IntStart_4wks,IntStart_8wks)]
#scale week no to be between 0 and 1
HAI.bySite_1[,'Week_no_scaled':=Week_no/max(Week_no)]
#Null model
model.0.output<-list()
model.0.formula<-formula(N~Week_no_scaled+(1|Site))
model.0.output[[1]]<-glmer(model.0.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='SAB',])
model.0.output[[2]]<-glmer(model.0.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='CDI',])
model.0.output[[3]]<-glmer(model.0.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='VRE',])
names(model.0.output)<-c('SAB','CDI','VRE')

#MODEL 1 - binary intervention switch 4 weeks into intervention phase
model.1.output<-list()
model.1.formula<-formula(N~Week_no_scaled+IntStart_4wks+(1|Site))
model.1.output[[1]]<-glmer(model.1.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='SAB',])
model.1.output[[2]]<-glmer(model.1.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='CDI',])
model.1.output[[3]]<-glmer(model.1.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='VRE',])
names(model.1.output)<-c('SAB','CDI','VRE')

#MODEL 2 - binary intervention switch 8 weeks into intervention phase (ie one month delay)
model.2.formula<-formula(N~Week_no_scaled+IntStart_8wks+(1|Site))
model.2.output<-list()
model.2.output[[1]]<-glmer(model.2.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='SAB',])
model.2.output[[2]]<-glmer(model.2.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='CDI',])
model.2.output[[3]]<-glmer(model.2.formula,family = 'poisson',offset=log(OBD.adj/10000),data=HAI.bySite_1[Infection_type=='VRE',])
names(model.2.output)<-c('SAB','CDI','VRE')

#extract AICs
model.1.AIC<-unlist(lapply(model.1.output,AIC))
model.2.AIC<-unlist(lapply(model.2.output,AIC))
model.AIC<-rbind(model.1.AIC,model.2.AIC)

#meta analysis - Model 1
model.1.estimates<-do.call('rbind',lapply(model.1.output,function(x) summary(x)$coefficients['IntStart_4wks',c('Estimate','Std. Error')]))
std_errors<-model.1.estimates[,'Std. Error']
estimates<-(model.1.estimates[,'Estimate'])

#forest plot
forest.plot.data_1<-setDT(data.frame(model.1.estimates),keep.rownames = T)
setnames(forest.plot.data_1,'rn','HAI')

forest.plot.data_1$HAI<-factor(forest.plot.data_1$HAI,levels=c('CDI','SAB','VRE'))
forest.plot.data_1$HAI<-factor(forest.plot.data_1$HAI,levels=rev(levels(forest.plot.data_1$HAI)))


g<-ggplot(forest.plot.data_1,aes(y=exp(Estimate),ymin=exp(Estimate-1.96*Std..Error),ymax=exp(Estimate+1.96*Std..Error),x=HAI))+geom_point(size=2)+geom_errorbar(width=0.5)+coord_flip()+geom_hline(aes(yintercept=1))
g2<-g+xlab('Healthcare associated infection')+ylab('Intervention effect (RR)')+annotate("text", y = c(0.8,1.2), x = c(0.5,0.5), label = c("Intervention reduces infections", "Intervention increases infections"))+theme_classic()+theme(text=element_text(size=14))
g2

#Meta analysis estimate- model 1
W<-1/(std_errors^2)
meta.est<-sum(W*estimates)/sum(W)
meta.se<-sqrt(1/sum(W))
meta.cilow<-meta.est+1.96*meta.se
meta.cihigh<-meta.est-1.96*meta.se
meta.pvalue<-2*pnorm(-abs(meta.est/meta.se))

#meta analysis - model 2
model.2.estimates<-do.call('rbind',lapply(model.2.output,function(x) summary(x)$coefficients['IntStart_8wks',c('Estimate','Std. Error')]))
std_errors<-model.2.estimates[,'Std. Error']
estimates<-(model.2.estimates[,'Estimate'])
W<-1/(std_errors^2)
meta.est<-sum(W*estimates)/sum(W)
meta.se<-sqrt(1/sum(W))
meta.cilow<-meta.est+1.96*meta.se
meta.cihigh<-meta.est-1.96*meta.se
meta.pvalue<-2*pnorm(-abs(meta.est/meta.se))

#pre/post infection rates
CDI.PrePostRate<-exp(cumsum(fixef(model.1.output$CDI)[c(1,3)]))
SAB.PrePostRate<-exp(cumsum(fixef(model.1.output$SAB)[c(1,3)]))
VRE.PrePostRate<-exp(cumsum(fixef(model.1.output$VRE)[c(1,3)]))


#bootstrapping for prediction intervals - model 1 - 2 year horizon from baseline
new.data<-data.table(Week_no=0:103,OBD.adj=10000)
new.data[,c('Week_no_scaled','Weeks_IntStart'):=list(Week_no/116,Week_no-52)]
new.data[,'IntStart_4wks':=ifelse(Weeks_IntStart>=4,1,0)]

new.data_null<-data.table(Week_no=0:103,OBD.adj=10000,Week_no_scaled=(0:103)/116,IntStart_4wks=0)

boot.PI <- function(.) {
  exp(predict(.,newdata=new.data,re.form=~0,type='link'))
}
boot.PI_null<-function(.){
  exp(predict(.,newdata=new.data_null,re.form=~0,type='link'))
}
boot.PI_AbsRisk<-function(.) {
  summary(.)$coefficients[,'Estimate']
}



#bootstrap model 1
SAB.boot<-bootMer(model.1.output$SAB, boot.PI, nsim=1000, use.u=TRUE, type="parametric")
CDI.boot<-bootMer(model.1.output$CDI, boot.PI, nsim=1000, use.u=TRUE, type="parametric")
VRE.boot<-bootMer(model.1.output$VRE, boot.PI, nsim=1000, use.u=TRUE, type="parametric")

#bootstrap null
SAB.boot_null<-bootMer(model.1.output$SAB, boot.PI_null, nsim=1000, use.u=TRUE, type="parametric")
CDI.boot_null<-bootMer(model.1.output$CDI, boot.PI_null, nsim=1000, use.u=TRUE, type="parametric")
VRE.boot_null<-bootMer(model.1.output$VRE, boot.PI_null, nsim=1000, use.u=TRUE, type="parametric")

#bootstrap absolute risk
SAB.boot_AbsRisk<-bootMer(model.1.output$SAB, boot.PI_AbsRisk, nsim=1000, use.u=TRUE, type="parametric")
CDI.boot_AbsRisk<-bootMer(model.1.output$CDI, boot.PI_AbsRisk, nsim=1000, use.u=TRUE, type="parametric")
VRE.boot_AbsRisk<-bootMer(model.1.output$VRE, boot.PI_AbsRisk, nsim=1000, use.u=TRUE, type="parametric")


HAI.boot<-rbindlist(list(melt(CDI.boot$t,varnames=c('Iteration','Time')), melt(SAB.boot$t,varnames=c('Iteration','Time')), melt(VRE.boot$t,varnames=c('Iteration','Time'))),idcol='HAI')
HAI.boot$HAI<-factor(HAI.boot$HAI,levels=1:3,labels=c('CDI','SAB','VRE'))

HAI.boot_null<-rbindlist(list(melt(CDI.boot_null$t,varnames=c('Iteration','Time')), melt(SAB.boot_null$t,varnames=c('Iteration','Time')), melt(VRE.boot_null$t,varnames=c('Iteration','Time'))),idcol='HAI')
HAI.boot_null$HAI<-factor(HAI.boot_null$HAI,levels=1:3,labels=c('CDI','SAB','VRE'))

HAI.boot_summary<-HAI.boot[,list(boot.mean=mean(value),boot.lower=quantile(value,0.025),boot.upper=quantile(value,0.975)),by=list(HAI,Time)]
HAI.boot_summary[,'Weeks_IntStart':=(Time-1)-52]

HAI.boot_null_summary<-HAI.boot_null[,list(boot.mean_null=mean(value),boot.lower_null=quantile(value,0.025),boot.upper_null=quantile(value,0.975)),by=list(HAI,Time)]
HAI.boot_null_summary[,'Weeks_IntStart':=(Time-1)-52]

setkeyv(HAI.boot_summary,c('HAI','Time'))
setkeyv(HAI.boot_null_summary,c('HAI','Time'))

HAI.boot_summary_all<-HAI.boot_summary[HAI.boot_null_summary,]

#build up trend plot with null model projections first
g.boot_null<-ggplot(HAI.boot_summary_all,aes(x=Weeks_IntStart,y=boot.mean_null,ymin=boot.lower_null,ymax=boot.upper_null,group=HAI))+geom_ribbon(alpha=.2)+geom_line(size=1,linetype='dashed')
g.boot_null<-g.boot_null+scale_y_continuous(breaks=seq(-3,4,0.25))+scale_x_continuous(breaks=seq(-52,52,8))  +xlab('Weeks since intervention start')+ylab('Rate per 10,000 occupied bed days')+theme_classic()+theme(text=element_text(size=14))+geom_vline(aes(xintercept=0),size=1,linetype='dotted')
#add model 1 trend lines
g.boot<-g.boot_null+geom_line(aes(x=i.Weeks_IntStart,y=boot.mean,group=HAI,colour=HAI),size=2)+geom_ribbon(aes(x=i.Weeks_IntStart,ymin=boot.lower,ymax=boot.upper,fill=HAI),alpha=.5)

g.boot

