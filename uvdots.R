gc();rm(list=ls())
require(data.table)
require(lme4)
require(ggplot2)
require(multcomp)
require(XLConnect)
require(influence.ME)

#load data#
source('FP_functions.R')


################################################
#FIGURE 1

#calc.percentage of dot locations clean - from start of control period up to a max. of 50 weeks
percent.clean<-audit.dt[,list(y=100*mean(Clean,na.rm=T),TimefromControlStart=min(Time.from.Start.Control),N.Ward=length(unique(Ward)),N.Room=length(unique(Room))), by = list(Site,Period,Audit_no,Dot_recode)]
percent.clean[,'Audit_label':=paste(Period,Audit_no,sep='_')]

#plot by Site as time trend
g<-ggplot(percent.clean,aes(x=TimefromControlStart,y=y,colour=Dot_recode))+geom_line(size=1.2)+geom_vline(aes(xintercept=8),linetype='dashed')+xlim(c(0,60))
g1<-g+facet_wrap(~Site)+theme(text=element_text(size=14),legend.title=element_blank())+xlab('Weeks since start of Control period')+ylab('% Dot locations clean')


################################################
#TABLE 1 - summary stats with respect to changes in cleaning performance: by period, and site
#number of dots
#number of dots declared clean
#percent clean
#define max number of audit weeks to truncate at
audit.end<-50
percent.clean.tab<-audit.dt[Time.from.Start.Control<=audit.end & !is.na(Clean),list(N=.N,N.Clean=sum(Clean),prop.Clean=round(100*mean(Clean),2)),by=list(Site,Period)]


#write table to text file for paper
write.table(percent.clean.tab,file='U:/Research/Projects/ihbi/aushsi/hospital_cleaning/Publications/Effectiveness paper/Tables and Figures/Table1.txt',col.names=T,row.names=F)


#descriptive statistics of FTP % clean by hospital
audit.clean_summary<-dcast(audit.clean[,mean(Clean_N/(Clean_N+NotClean_N)),by=list(Period,Site,Room)],Site+Room~Period,value.var='V1')
audit.clean_summary[,'CleanImproved':=ifelse(Intervention>Control,1,0)]
audit.clean_summary[,'IntDiff':=Intervention-Control]
#change labels for plotting
audit.clean_summary[,'Hospital':=factor(Site,labels=1:11)]
ggplot(audit.clean_summary,aes(x=Hospital,y=100*IntDiff,fill=Room))+
  geom_bar(stat='identity',position='dodge')+coord_flip()+
  xlab('Hospital')+ylab('% Frequent touch points cleaned (Intervention-Control)')+
  scale_y_continuous(breaks=seq(-30,70,10))+theme_bw()+
  theme(text=element_text(size=14))

#################################################
#GLMM analysis

#take audit.dt and determine the number of clean/not clean sites by Site, Period and Audit_no
#inclue up to ten intervention audits (described in Period and Audit_no)
audit.clean<-audit.dt[Time.from.Start.Control<=audit.end & !is.na(Clean),list(Clean_N=sum(Clean==1,na.rm=T),NotClean_N=sum(Clean==0,na.rm=T)),by=list(Site,Period,Audit_no,Dot_recode,Time.from.Start.Control,Time.from.Start.Int)]
setnames(audit.clean,'Dot_recode','Room')
#set up scaling of time covariates - control and intervention
scale.control<-compute_scale_x(audit.clean$Time.from.Start.Control)
scale.int<-compute_scale_x(audit.clean$Time.from.Start.Int)  


#scale.control<-max(audit.clean$Time.from.Start.Control)
#scale.int<-max(audit.clean$Time.from.Start.Int)
#add scaled versions to audit.clean (in place of using scale() in each model formula)
audit.clean[,c('Time.from.Start.Control_scaled','Time.from.Start.Int_scaled'):=list(Time.from.Start.Control/scale.control,Time.from.Start.Int/scale.int)]

#create generic list to store AIC for each model
model.AIC<-list()


#null model 
model.A0<-glmer(cbind(Clean_N, NotClean_N) ~ Room + (1 | Site),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000)))
model.AIC$A0<-summary(model.A0)$AIC['AIC']

#models from protocol paper
#1. a simple binary intervention (yes/no)
model.A1 <- glmer(cbind(Clean_N, NotClean_N) ~ Period + Room + Period:Room + (1 | Site),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000)))
model.AIC$A1<-summary(model.A1)$AIC['AIC']


#2. a linear intervention using the time since intervention
model.A2 <- glmer(cbind(Clean_N, NotClean_N) ~ Time.from.Start.Int_scaled + Room + Time.from.Start.Int_scaled:Room + (1 | Site),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000)))
model.AIC$A2<-summary(model.A2)$AIC['AIC']

#additional models

#interrupted time series model: trend from start of intervention period (all interactions included)  
model.A3 <- glmer(cbind(Clean_N, NotClean_N) ~ Period + Room + Time.from.Start.Int_scaled + Room:(Period+Time.from.Start.Int_scaled) + (1 | Site),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000)))
model.AIC$A3<-summary(model.A3)$AIC['AIC']

#TABLE 2: extract AIC for each model
A.AIC<-data.table(Model=paste0('A',0:3),AIC=unlist(model.AIC,use.names=F))
#########################################################
#TABLE 3: parameter estimate summary for best fitting model by AIC

best.model.AIC<-A.AIC[which.min(AIC),]$Model
best.model.summary<-summary(get(paste('model',best.model.AIC,sep='.')))$coefficients

# #########################################################
#TABLE 4: Predictions and log ORs
#predict under best fitting model (Model A3)
# want predictions of % clean by dot location for (i) contorl period, (ii), start of intervention period, (iii) end of intervention period
Int.last<-10

IntTimeStart_mean<-audit.clean[Period=='Intervention' & Audit_no==1,mean(Time.from.Start.Int_scaled)]
IntTimeEnd_mean<-audit.clean[Period=='Intervention' & Audit_no==Int.last,mean(Time.from.Start.Int_scaled)]

#in IntTime, include 0 to better graph step change
IntTime<-c(10e-3,seq(IntTimeStart_mean,IntTimeEnd_mean,length.out=Int.last))
ControlTimes<-IntTimeStart_mean-((IntTimeEnd_mean-IntTimeStart_mean)/Int.last)*(2:1)

new.data<-data.table(Period=rep(rep(c('Control','Intervention'),c(2,length(IntTime))),2),Room=rep(c('Bathroom','Bedroom'),each=2+length(IntTime)),Time.from.Start.Int_scaled=rep(c(0,0,IntTime),2),Time=rep(c(ControlTimes,IntTime),2))
new.data[,'pred':=predict(model.A3,newdata=new.data,re.form=NA,type='response')]

audit.clean[,'pred.modelA3':=predict(model.A3,re.form=NA,type='response')]

#use ghlt to find log OR of step change by Dot location

c1<-rbind('Bathroom: Control vs Intervention (Step change)'=c(0,1,0,0,0,0),
          'Bedroom: Control vs Intervention (Step change)'=c(0,1,0,0,1,0),
          'Bedroom vs. Bath: Control vs intervention (Step change)'=c(0,0,0,0,1,0),
          'Bathroom: Intervention Time'=c(0,0,0,1,0,0),
          'Bedroom: Intervention Time'=c(0,0,0,1,0,1),
          'Bedroom vs Bath: Intervention Time'=c(0,0,0,0,0,1)
)

OR.summary<-summary(glht(model.A3,c1))

tmp<-rbindlist(list(summary(OR.summary)$test[c('coefficients','sigma','pvalues')]))
OR.tab<-tmp[,list(hypotheis=rownames(c1),coefficients,coefficients-1.96*sigma,coefficients+1.96*sigma,pvalues)]
rm(tmp)

#add relative risks compared with control
#model based predictions for control period

control.pred<-predict(model.A3,newdata=data.frame(Period='Control',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=0),re.form=NA,type='response')
int0.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=0),re.form=NA,type='response')
int4.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=0.4),re.form=NA,type='response')

int10.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=1),re.form=NA,type='response')
int20.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=2),re.form=NA,type='response')
int30.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=3),re.form=NA,type='response')
int40.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=4),re.form=NA,type='response')



#plot over raw data
percent.clean.dot<-audit.dt[Time.from.Start.Control<=audit.end,list(y=100*mean(Clean,na.rm=T),TimefromControlStart=min(Time.from.Start.Control),N.Ward=length(unique(Ward)),N.Room=length(unique(Room))), by = list(Site,Period,Audit_no,Dot_recode)]
setnames(percent.clean.dot,'Dot_recode','Room')
percent.clean.dot[,'Audit_label':=paste(Period,Audit_no,sep='_')]
#change dot_recode to room, to be consistent with main paper

g<-ggplot(percent.clean.dot,aes(x=TimefromControlStart,y=y,colour=Room))+geom_point(size=2,alpha=0.5)
g2<-g+geom_vline(aes(xintercept=8),size=1,linetype='dotted')+geom_line(data=audit.clean,aes(x=Time.from.Start.Control,y=100*pred.modelA3,linetype=Room),size=1)+theme(text=element_text(size=14),legend.title=element_blank())+xlab('Weeks since start of control period')+ylab('% Frequent touch points cleaned')+theme_classic()

####################################################
control.pred<-predict(model.A3,newdata=data.frame(Period='Control',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=0),re.form=NA,type='response')
int0.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=0),re.form=NA,type='response')
int10.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=1),re.form=NA,type='response')
int20.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=2),re.form=NA,type='response')
int30.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=3),re.form=NA,type='response')
int40.pred<-predict(model.A3,newdata=data.frame(Period='Intervention',Room=c('Bathroom','Bedroom'),Time.from.Start.Int_scaled=4),re.form=NA,type='response')

model.pred_tab<-data.table(do.call('rbind',list(control.pred,int10.pred,int20.pred,int30.pred,int40.pred)))
model.pred_tab[,'Period':=c('Control','Int_10wks','Int_20wks','Int_30wks','Int_40wks')]
setnames(model.pred_tab,1:2,c('Bathroom','Bedroom'))

#add extra columns to model.pred_tab for plot annotations
model.pred_tab[,'X':=c(6,18,28,38,48)]


g<-ggplot(percent.clean.dot,aes(x=TimefromControlStart,y=y,colour=Room))+geom_point(size=2,alpha=0.5)
g2<-g+annotate("segment",x=8,xend=8,y=0,yend=100,linetype='dotted')+geom_line(data=new.data,aes(x=Time*scale.int+8,y=100*pred,linetype=Room),size=1)+xlab('Weeks since start of control period')+ylab('% Frequent touch points cleaned')+theme_classic()+theme(text=element_text(size=14))+annotate("text",x=c(model.pred_tab$X),y=rep(110,5),label=100*round(model.pred_tab$Bathroom,2),colour="#F8766D",fontface=2)+annotate("text",x=model.pred_tab$X,y=rep(105,5),label=100*round(model.pred_tab$Bedroom,2),colour="#00BFC4",fontface=2)+scale_y_continuous(breaks=seq(0,100,25))
g2


#plot over raw data
percent.clean.dot<-audit.dt[!Audit_label %in% exclude.audit,list(y=100*mean(Clean,na.rm=T),TimefromControlStart=min(Time.from.Start.Control),N.Ward=length(unique(Ward)),N.Room=length(unique(Room))), by = list(Site,Period,Audit_no,Dot_recode)]
percent.clean.dot[,'Audit_label':=paste(Period,Audit_no,sep='_')]

setnames(percent.clean.dot,'Dot_recode','Room')
setnames(new.data,'Dot_recode','Room')

g<-ggplot(percent.clean.dot,aes(x=TimefromControlStart,y=y,colour=Room))+geom_point(size=2,alpha=0.5)
g2<-g+annotate("segment",x=8,xend=8,y=0,yend=100,linetype='dotted')+geom_line(data=new.data,aes(x=Time*scale.int+8,y=100*pred,linetype=Room),size=1)+xlab('Weeks since start of control period')+ylab('% Frequent touch points cleaned')+theme_classic()+theme(text=element_text(size=14))+annotate("text",x=c(model.pred_tab$X),y=rep(110,5),label=100*round(model.pred_tab$Bathroom,2),colour="#F8766D",fontface=2)+annotate("text",x=model.pred_tab$X,y=rep(105,5),label=100*round(model.pred_tab$Bedroom,2),colour="#00BFC4",fontface=2)+scale_y_continuous(breaks=seq(0,100,25))


##############################################
IntTime_scaled<-seq(0,40,10)/scale.int
new.data.Int<-data.table(Time.from.Start.Int_scaled=rep(IntTime_scaled,2),Period='Intervention',Room=rep(c('Bathroom','Bedroom'),each=length(IntTime_scaled)))
new.data.Control<-data.table(Time.from.Start.Int_scaled=c(0,0),Period='Control',Room=c('Bathroom','Bedroom'))
new.data<-rbindlist(list(new.data.Control,new.data.Int))
new.data[,'pred':=predict(model.A3,newdata=new.data,re.form=NA,type='response')]

mySumm_bath <- function(.) {
  #predict(., newdata=sleepstudy, re.form=NULL)
  predict(.,newdata=new.data[Room=='Bathroom',],re.form=~0,type='response')
}
mySumm_bed <- function(.) {
  #predict(., newdata=sleepstudy, re.form=NULL)
  predict(.,newdata=new.data[Room=='Bedroom',],re.form=~0,type='response')
}
Bath.boot<-bootMer(model.A3, mySumm_bath, nsim=1000, use.u=TRUE, type="parametric")
Bed.boot<-bootMer(model.A3, mySumm_bed, nsim=1000, use.u=TRUE, type="parametric")

Bath.boot_summary<-data.table(melt(Bath.boot$t,varnames=c('Iteration','Time')))
Bed.boot_summary<-data.table(melt(Bed.boot$t,varnames=c('Iteration','Time')))

Bath.boot_summary<-Bath.boot_summary[,list(mean.boot=mean(value),lower95.boot=quantile(value,0.025),upper95.boot=quantile(value,0.975)),by=Time]
Bed.boot_summary<-Bed.boot_summary[,list(mean.boot=mean(value),lower95.boot=quantile(value,0.025),upper95.boot=quantile(value,0.975)),by=Time]


FTP.boot_summary<-rbindlist(list(Bath.boot_summary,Bed.boot_summary),idcol=T)
setnames(FTP.boot_summary,'.id','Room')
FTP.boot_summary[,'Room':=factor(Room,levels=1:2,labels=c('Bathroom','Bedroom'))]
FTP.boot_summary[,'Time_1':=rep(c('Control',paste0('Int',0:4)),2)]

#plot
g<-ggplot(FTP.boot_summary[Time_1 %in% c('Control','Int0','Int2','Int4')],aes(x=Time_1,y=100*mean.boot,ymin=100*lower95.boot,ymax=100*upper95.boot,group=Room))
g+geom_bar(aes(fill=Room),position=position_dodge(), stat="identity")+  
  geom_errorbar(colour="black", width=0.1, position=position_dodge(0.9))+
  xlab('')+ylab('% Frequent touch points cleaned')+
  theme_bw() + theme(text=element_text(size=14),axis.text.x = element_text(angle=45,hjust=1))+
  scale_x_discrete(labels= c('Control','First Intervention audit','20 weeks implementation','40 weeks implementation'))+
  scale_y_continuous(breaks=seq(0,90,10))


################################################
#leave one out analysis for model A3
model.A3 <- glmer(cbind(Clean_N, NotClean_N) ~ Period + Room + Time.from.Start.Int_scaled + Room:(Period+Time.from.Start.Int_scaled) + (1 | Site),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000)))

sens_model.A3<-influence(model.A3,"Site")
cooksd_model.A3<-cooks.distance(sens_model.A3)
dfbetas_model.A3<-setDT(data.frame(dfbetas(sens_model.A3)),keep.rownames = T)
colnames(dfbetas_model.A3)<-c('Site','Intercept','Intervention','Bedroom','IntTime','IntbyRoom','IntTimebyRoom')
dfbetas_model.A3.long<-melt(dfbetas_model.A3,id.vars='Site',value.name='DFBETA')

pchange_model.A3<-setDT(data.frame(pchange(sens_model.A3)),keep.rownames = T)
colnames(pchange_model.A3)<-c('Site','Intercept','Intervention: Step change','Bedroom','Intervention: Time','IntbyRoom','IntTimebyRoom')
pchange_model.A3.long<-melt(pchange_model.A3,id.vars='Site',value.name='percent.change')

g<-ggplot(pchange_model.A3.long[!variable %in% c('IntbyRoom','IntTimebyRoom'),],aes(x=Site,y=percent.change,group=variable))+geom_histogram(stat='identity')+facet_wrap(~variable)+coord_flip()

g+ylab('% change in parameter estimate')+xlab('Hospital excluded')+theme_bw()+theme(text=element_text(size=14),strip.background = element_rect(fill='white'))

#################################################
#Extension to fractional polynomials (FP1)
#################################################
###########################################################
#Fit FP1 models - including interactino with dot locations
#define set of FP powers and corresponding labels for data.table
FP.powers<-c(-2,-1,-.5,0,.5,1,2,3)
FP.cols<-paste('FP',gsub('-','minus',FP.powers),sep='_')

#add FP transforms as variables to audit.clean - names in 'FP.cols'
#based on time from intervention start
#for control period, variable = 0
audit.clean[,c(FP.cols):=0]
audit.clean[Period=='Intervention',c(FP.cols):=lapply(FP.powers,function(p) FP_transform_x_1(Time.from.Start.Int_scaled[Period=='Intervention'],p,shift=1))]

# #add additional columns for FPs with repeated powers
# audit.clean[,c(paste(FP.cols,'rp',sep='_')):=0]
# audit.clean[Period=='Intervention',c(paste(FP.cols,'rp',sep='_')):=lapply(FP.powers,function(p) FP_transform_x_2(Time.from.Start.Control_scaled[Period=='Intervention'],p))]

#define model formula for FP1 model
FP1.formula<-function(p){formula(paste0('cbind(Clean_N, NotClean_N) ~ Period + Room + ',p,'+ Room:(Period + ', p,') + (1|Site)'))}

#fit all FP1 models and store as a list
model.FP1<-lapply(FP.cols,function(p) glmer(FP1.formula(p),data=audit.clean, family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl = list(maxfun=50000))))
names(model.FP1)<-paste(FP.cols,'1',sep='_')
#extract AICs for each FP model
FP1.AIC<-data.table(Model=names(model.FP1),AIC=extract_AIC_FP(model.FP1))
FP1.AIC[,'Model':=factor(Model,levels=paste(FP.cols,'1',sep='_'),labels=FP.powers)]
############################################################################
#ORs

c1<-rbind('Bathroom: Control vs Intervention (Step change)'=c(0,1,0,0,0,0),
          'Bedroom: Control vs Intervention (Step change)'=c(0,1,0,0,1,0),
          'Bedroom vs. Bath: Control vs intervention (Step change)'=c(0,0,0,0,1,0),
          'Bathroom: Intervention'=c(0,0,0,40/scale.int,0,0),
          'Bedroom: Intervention'=c(0,0,0,40/scale.int,0,40/scale.int),
          'Bedroom vs Bath: Intervention'=c(0,0,0,0,0,40/scale.int)
)

OR.summary.FP1<-summary(glht(model.FP1[[3]],c1))
tmp<-rbindlist(list(summary(OR.summary.FP1)$test[c('coefficients','sigma','pvalues')]))
OR.tab.FP1<-tmp[,list(hypotheis=rownames(c1),coefficients,coefficients-1.96*sigma,coefficients+1.96*sigma,pvalues)]
rm(tmp)
##############################################################################
#Predict under best fitting  FP1 - add to current new.data object so scaled is the same
new.data[,'FP_minus0.5':=0]
new.data[Period=='Intervention','FP_minus0.5':=FP_transform_x_1(Time.from.Start.Int_scaled[Period=='Intervention'],-0.5,shift=1)]
new.data[,'pred.FPminus05':=predict(model.FP1[[3]],newdata=new.data,re.form=NA,type='response')]
#Plot predictions
plot.pred<-melt(new.data,id.vars=c('Period','Room','Time.from.Start.Int_scaled','FP_minus0.5','Time'),variable.name = 'Model',value.name = 'Prediction')
plot.pred$Model<-factor(plot.pred$Model,levels=c('pred','pred.FPminus05'),labels=c('Model A3','FP1'))
plot.pred[,'Time.from.Start.Control_scaled':=((scale.int*Time)+8)/scale.control]

g<-ggplot(percent.clean.dot,aes(x=TimefromControlStart,y=y,colour=Room))+geom_point(size=2,alpha=0.5)

g4<-g+geom_vline(aes(xintercept=8),size=1,linetype='dashed')+geom_line(data=plot.pred,aes(x=(Time.from.Start.Control_scaled)*scale.control,y=100*Prediction,colour=Room,linetype=Model),size=1)+theme_bw()+xlab('Weeks since start of control period')+ylab('% Frequent touch points cleaned')+theme(text=element_text(size=14),legend.title=element_blank())
g4

#bootstrapepd predictions under fp1(-0.5) model
new.data_boot<-data.table(Period=rep(c('Control','Intervention'),c(2,10)),Room=rep(c('Bathroom','Bedroom'),6),Time.from.Start.Int_scaled=rep(c(0,0,10,20,30,40)/scale.int,each=2))
new.data_boot[,'FP_minus0.5':=0]
new.data_boot[Period=='Intervention','FP_minus0.5':=FP_transform_x_1(Time.from.Start.Int_scaled[Period=='Intervention'],-0.5,shift=1)]



mySumm_bath <- function(.) {
  #predict(., newdata=sleepstudy, re.form=NULL)
  predict(.,newdata=new.data_boot[Room=='Bathroom',],re.form=~0,type='response')
}
mySumm_bed <- function(.) {
  #predict(., newdata=sleepstudy, re.form=NULL)
  predict(.,newdata=new.data_boot[Room=='Bedroom',],re.form=~0,type='response')
}




Bath.boot_fp1<-bootMer(model.FP1[[3]], mySumm_bath, nsim=1000, use.u=TRUE, type="parametric")
Bed.boot_fp1<-bootMer(model.FP1[[3]], mySumm_bed, nsim=1000, use.u=TRUE, type="parametric")

Bath.boot_fp1_summary<-data.table(melt(Bath.boot_fp1$t,varnames=c('Iteration','Time')))
Bed.boot_fp1_summary<-data.table(melt(Bed.boot_fp1$t,varnames=c('Iteration','Time')))

Bath.boot_fp1_summary<-Bath.boot_fp1_summary[,list(mean.boot=mean(value),lower95.boot=quantile(value,0.025),upper95.boot=quantile(value,0.975)),by=Time]
Bed.boot_fp1_summary<-Bed.boot_fp1_summary[,list(mean.boot=mean(value),lower95.boot=quantile(value,0.025),upper95.boot=quantile(value,0.975)),by=Time]
################################################
