######################################################
## Code to generate constituent models for super-   ##
##  learner, given a fully constructed data matrix, ##
##  and run all principal analyses used for paper.  ##
######################################################
##
## James Liley, 2020
##


######################################################
## Directories, packages and scripts                ##
######################################################

# Working directory and location (windows or linux machine)
if (file.exists("C:/Users/jliley_1718-0370/Documents")) {
	location="Windows"
	setwd("//farr-fs1/Study Data/1718-0370/Research/Code/")
} else {
	location="Linux"
	setwd("/Study_Data/Research/Code")
}

# Scripts
lx=list.files("James/SPARRAv4/",pattern="*.R",full.names=TRUE)
for (i in 1:length(lx)) source(lx[i])

# Packages
if (location=="Windows") {
  library(pROC)
  library(PRROC)
  library(xgboost)
  library(glmnet)
  library(Metrics)
  library(data.table)
  library(cvAUC)
  library(latex2exp)
} else {
  library(h2o)
  library(tidyverse)
  library(fst)
}

# Directories
dir_cleanData = "../Data/Data_clean/"
model_dir="James/Analysis/Data/full_model/"
out_dir="James/Analysis/Output/full_model/"
plot_dir="James/Analysis/Output/full_model/"

# Specific file locations
fullmatrix=paste0(model_dir,"full_data_matrix.RDS") # Full training matrix location
design_exclusions_file=paste0(model_dir,"design_exclusions.RData")

# In text outputs, everything after this marker should be readable by R.
sink_marker="\n\n\n***************\n\n"
options(width=1e4)






######################################################
## File paths for data matrices and models          ##
######################################################


loc1=paste0(model_dir,"models_cv/model_fold1") # Model for prediction on fold 1 stored in various files {loc1}.RDS, {loc1}.pred.RDS, {loc1}.model..
loc2=paste0(model_dir,"models_cv/model_fold2")
loc3=paste0(model_dir,"models_cv/model_fold3")
output1=paste0(loc1,".output") # Predictions for fold 1 stored in {output1}.RDS
output2=paste0(loc2,".output")
output3=paste0(loc3,".output")


# Models fitted to design matrices with no topic features
loc1_notopic=paste0(model_dir,"models_cv/model_fold1_notopic") 
output1_notopic=paste0(loc1_notopic,".output")



######################################################
## Ensure determinacy                               ##
######################################################

seed=220
set.seed(seed)



######################################################
## Collect predictions, type and time of events     ##
######################################################

all_pred_file=paste0(out_dir,"all_pred.RDS")
if (!file.exists(all_pred_file)) {
  
  # Read in training matrix
  nhs=readRDS(fullmatrix)
  load(design_exclusions_file)
  nhs=nhs[v34,] # redefine nhs for memory considerations 
  rm(list=c("v34","v3x"))
  gc()
  
  # Load patients, episodes, list_of_data_tables
  load_cleanData(partition = "all",
                 subset_frac=1, subset_seed = 1234,
                 load_to_global=TRUE)
  
  # Gather predictions 
  pred1=readRDS(paste0(output1,".RDS"))
  pred2=readRDS(paste0(output2,".RDS"))
  pred3=readRDS(paste0(output3,".RDS"))
  
  # Predictions from fold 1 with no topic features
  pred1_notopic=readRDS(paste0(output1_notopic,".RDS"))
  
  
  ## Collect predictions from various constituents and super-learner
  cpred=grep("pred.",names(pred1),value=T)
  all_pred=data.frame(matrix(0,dim(nhs)[1],length(cpred)+1))
  colnames(all_pred)=c(cpred,"super")
  
  for (i in 1:3) {
    w=which(nhs$cv==i)
    for (j in 1:length(cpred)) {
      all_pred[w,j]=get(paste0("pred",i))[[j]]
    }
    all_pred$super[w]=get(paste0("pred",i))$final.full
  }
  
  ## Collect predictions from super-learner for predictor with no topic features
  all_pred$super_notopic=rep(NA,dim(all_pred)[1])
  w1=which (nhs$cv==1)
  all_pred$super_notopic[w1]=pred1_notopic$final.full
  
  ## Affix various other useful data
  all_pred=cbind(id=nhs$id,cv=patients$cv[match(nhs$id,patients$id)],
                 time=as.numeric(nhs$time),
                 cv=nhs$cv,
                 v3=as.numeric(nhs$v3score)/100,
                 age=nhs$age,sexM=(patients$gender[match(nhs$id,patients$id)]=="Male"),
                 simd=nhs$SIMD_DECILE_2016_SCT,ae2=nhs$num_ae2_attendances,
                 all_pred,target=as.numeric(nhs$target),reason=nhs$reason)
  
  wpos=which(nhs$target==TRUE)
  
  cutoffs=sort(as_datetime(unique(nhs$time)))
  SMR01M=list_of_data_tables$SMR01M

  # Initialise columns
  admission_type=rep(NA,dim(all_pred)[1])
  main_condition=rep(NA,dim(all_pred)[1])
  cis_marker=rep(NA,dim(all_pred)[1])
  time_to_first_event=rep(NA,dim(all_pred)[1])
  n_event=rep(NA,dim(all_pred)[1])
  

  for (it in 1:length(cutoffs)) {
    
    t0=cutoffs[it] # as_datetime(nhs$time[1]) # time cutoff for SPARRA score
    tsub=which(nhs$time==t0)
    idsub=nhs$id[tsub]
    
    tdiff=round(as.integer(SMR01M$time-t0)/(24*3600)) # time since time cutoff for each episode
    
    w1=which(SMR01M$id %in% intersect(idsub,nhs$id[wpos])) # individuals who had an AE
    w2=which(tdiff> -1 & tdiff<367)
    subep=SMR01M[intersect(w1,w2),] # all episodes in the right time frame, involving the right samples
    
    subep=subep[which(subep$emergency_admin==1),] # now only valid events, as per v4
    
    # Number of admissions. Sometimes events are recorded more than once,
    #  and conditions etc are only recorded for one of them, so we have to 
    #  be careful. We can sort this out using 'order' since it puts NAs last
    subep=subep[order(subep$main_condition),]
    idt=paste0(subep$id,"_",as.integer(subep$time))
    subep=subep[match(unique(idt),idt),]
    tn=table(subep$id)
    
    # only consider first episode for each sample
    subep=subep[order(subep$time),]
    subep=subep[match(unique(subep$id),subep$id),]
    
    # Admission type. We can match IDs because subep only contains each ID once.
    admission_type0=rep(NA,length(tsub))
    admission_type0[match(subep$id,idsub)]=subep$admission_type
    admission_type[tsub]=admission_type0
    
    # Main condition
    main_condition0=rep(NA,length(tsub))
    main_condition0[match(subep$id,idsub)]=subep$main_condition
    main_condition[tsub]=main_condition0
    
    # CIS marker
    cis_marker0=rep(NA,length(tsub))
    cis_marker0[match(subep$id,idsub)]=subep$cis_marker
    cis_marker[tsub]=cis_marker0
    
    # Time to first event
    time_to_first_event0=rep(NA,length(tsub))
    time_to_first_event[match(subep$id,all_pred$id)]=as.integer(subep$time-t0)/(24*3600)
    time_to_first_event[tsub]=time_to_first_event0
    
    # Number of events
    n_event0=rep(NA,length(tsub))
    n_event0[match(as.numeric(names(tn)),idsub)]=tn
    n_event[tsub]=n_event0
    
  }
  
  # Affix to all_pred
  all_pred=cbind(all_pred,admission_type,main_condition,cis_marker,time_to_first_event,n_event)
  all_pred=all_pred[which(!is.na(all_pred$cv)),]
  
  saveRDS(all_pred,file=all_pred_file)

} else all_pred=readRDS(all_pred_file)

if (exists("nhs")) rm(list=c("nhs"))
if (exists("patients")) rm(list=c("patients","episodes","list_of_data_tables"))




######################################################################################
######################################################################################
## Topic inclusion. The following scripts compare performance of models fitted to   ##
##  design matrices which either contain or do not contain topic-model derived      ##
##  features.                                                                       ##
######################################################################################
######################################################################################

# Relevant data; subset of rows of all_pred corresponding to fold 1
sub1=which(all_pred$cv==1)

y=all_pred$target[sub1]
yt=all_pred$super[sub1]
ynt=all_pred$super_notopic[sub1]


## Test of discrimination (difference in AUCs) and practical difference
library(pROC)
roc_topic=roc(y,yt)
roc_notopic=roc(y,ynt)

sink(paste0(plot_dir,"topics/topic_comparison.txt"))
r0=roc.test(roc_topic,roc_notopic)
r0$roc1=NULL; r0$roc2=NULL
print(r0)
cat("\n\n")
cat("P-value: ", r0$p.value)
cat("\n\n\n")

# Practical difference: number of true positives amongst top-N
topn=c(1000,2000,5000,1e4,5e4,1e5,5e5,1e6)
out=c()
ord_topic=order(-yt)
ord_notopic=order(-ynt)
for (i in 1:length(topn)) {
  n_top=sum(y[ord_topic[1:topn[i]]])
  n_notop=sum(y[ord_notopic[1:topn[i]]])
  out=cbind(out,c(n_top,n_notop,n_top-n_notop))
}
out=rbind(topn,out)
rownames(out)=c("TopN","N_adm_topics","N_adm_notopics","Difference")
print(out)

cat("\n\n")
cat("")

cat(sink_marker)
cat("r0="); dput(r0); cat("\n\n")
cat("out="); dput(out); cat("\n\n")

sink()




## ROC curves

xrocT=getroc(y,yt) # Topics included
xrocNT=getroc(y,ynt) # Topics not included

pdf(paste0(plot_dir,"topics/topic_roc.pdf"),width=3,height=3.5)
aucT=signif(xrocT$auc,digits=4); seT=signif(xrocT$se,digits=2);
aucNT=signif(xrocNT$auc,digits=4); seNT=signif(xrocNT$se,digits=2);
labs=c(paste0("Topics: ",aucT," (",seT,")"),
       paste0("No topics: ",aucNT," (",seNT,")"))
xcol=c("red","black")
roc_2panel(list(xrocT,xrocNT),labels=labs,col=xcol,title="AUROC (SE)",text.col=xcol,title.col="black")
dev.off()

sink(paste0(plot_dir,"topics/topic_roc.txt"))
cat("Topics","\n")
cat("Sensitivity\n",xrocT$sens,"\n\nSpecificity\n",xrocT$spec,"\n\nAUROCs\n",xrocT$auc,"\n\nSE\n",xrocT$se,"\n\n\n")
cat("\n\n")
cat("NoTopics","\n")
cat("Sensitivity\n",xrocNT$sens,"\n\nSpecificity\n",xrocNT$spec,"\n\nAUROCs\n",xrocNT$auc,"\n\nSE\n",xrocNT$se,"\n\n\n")
cat("\n\n")
cat(sink_marker)
cat("xrocT="); dput(xrocT); cat("\n\n")
cat("xrocNT="); dput(xrocNT); cat("\n\n")
sink()












### PR curves

xprcT=getprc(y,yt)
xprcNT=getprc(y,ynt) 

pdf(paste0(plot_dir,"topics/topic_prc.pdf"),width=3,height=3.5)
aucT=signif(xprcT$auc,digits=4); seT=signif(xprcT$se,digits=2);
aucNT=signif(xprcNT$auc,digits=4); seNT=signif(xprcNT$se,digits=2);
labs=c(paste0("Topics: ",aucT," (",seT,")"),
       paste0("No topics: ",aucNT," (",seNT,")"))
xcol=c("red","black")
prc_2panel(list(xprcT,xprcNT),labels=labs,col=xcol,title="AUROC (SE)",text.col=xcol,title.col="black")
dev.off()

sink(paste0(plot_dir,"topics/topic_prc.txt"))
cat("Topics","\n")
cat("Sensitivity\n",xprcT$sens,"\n\nPPV\n",xprcT$ppv,"\n\nAUROCs\n",xprcT$auc,"\n\nSE\n",xprcT$se,"\n\n\n")
cat("\n\n")
cat("NoTopics","\n")
cat("Sensitivity\n",xprcNT$sens,"\n\nPPV\n",xprcNT$ppv,"\n\nAUROCs\n",xprcNT$auc,"\n\nSE\n",xprcNT$se,"\n\n\n")
cat("\n\n")
cat(sink_marker)
cat("xprcT="); dput(xprcT); cat("\n\n")
cat("xprcNT="); dput(xprcNT); cat("\n\n")
sink()




### Calibration curves
xcalT=plotcal(y,yt,plot=FALSE,kernel=TRUE)
xcalNT=plotcal(y,ynt,plot=FALSE,kernel=TRUE)

pdf(paste0(plot_dir,"topics/topic_cal.pdf"),width=3,height=3.5)
labs=c("Topics","No topics")
xcol=c("red","black")
cci=rgb(0.5,0.5,0.5,alpha=0.5) # colour for confidence envelope
cal_2panel(list(xcalT,xcalNT),labels=labs,col=xcol,text.col=xcol,ci_col=c(NA,cci))
dev.off()

sink(paste0(plot_dir,"topics/topic_cal.txt"))
cat("Topics\n\n")
cat("Obs\n",xcalT$x,"\n\nPred\n",xcalT$y,"\n\nLower\n",xcalT$lower,"\n\nUpper\n",xcalT$upper)
cat("\n\n\nNoTopics\n\n")
cat("Obs\n",xcalNT$x,"\n\nPred\n",xcalNT$y,"\n\nLower\n",xcalNT$lower,"\n\nUpper\n",xcalNT$upper)
cat(sink_marker)
cat("xcalT="); dput(xcalT); cat("\n\n")
cat("xcalNT="); dput(xcalNT); cat("\n\n")
sink()




# Calibration at disagreement points

pt=0.05; # Level of disagreement

w_t=which(yt-ynt> pt); # Higher score for topics-included
w_nt=which(yt-ynt< -pt); # Higher score for topics-not-included

# True and predicted values for higher score with topics-included
y_t=y[w_t]; 
ypn_t=ynt[w_t]; yp_t=yt[w_t]

# True and predicted values for higher score with topics-not-included
y_nt=y[w_nt]; 
ypn_nt=ynt[w_nt]; yp_nt=yt[w_nt]

# Calibration for higher scores with topics-included
cal_t=plotcal(y_t,yp_t,plot=F); # Calibration for topics-included
caln_t=plotcal(y_t,ypn_t,plot=F); # Calibration for topics-not-included

# Calibration for higher scores with topics-not-included
cal_nt=plotcal(y_nt,yp_nt,plot=F); # Calibration for topics-included
caln_nt=plotcal(y_nt,ypn_nt,plot=F); # Calibration for topics-not-included

# Plot calibration curves
pdf(paste0(plot_dir,"topics/topic_cal_dif.pdf"),width=3,height=3.5)

labs=c("T, T higher","NT, T higher","T, NT higher","NT, NT higher")
xcol=c("blue","black","blue","black")
xlty=c(1,1,2,2)
cci1=rgb(0,0,1,alpha=0.3) # colour for confidence envelope
cci2=rgb(0.5,0.5,0.5,alpha=0.3) # colour for confidence envelope
cal_2panel(list(cal_t,caln_t,cal_nt,caln_nt),labels=labs,col=xcol,lty=xlty,text.col=xcol,ci_col=c(cci1,cci2,cci1,cci2))
dev.off()

sink(paste0(plot_dir,"topics/topic_cal_dif.txt"))
cat("topicsIncluded_topicsHigher\n\n")
cat("Obs\n",cal_t$x,"\n\nPred\n",cal_t$y,"\n\nLower\n",cal_t$lower,"\n\nUpper\n",cal_t$upper)
cat("\n\n\ntopicsNotIncluded_topicsHighers\n\n")
cat("Obs\n",caln_t$x,"\n\nPred\n",caln_t$y,"\n\nLower\n",caln_t$lower,"\n\nUpper\n",caln_t$upper)
cat("\n\n\ntopicsIncluded_NoTopicsHigher\n\n")
cat("Obs\n",cal_t$x,"\n\nPred\n",cal_t$y,"\n\nLower\n",cal_t$lower,"\n\nUpper\n",cal_t$upper)
cat("\n\n\ntopicsNotIncluded_NoTopicsHighers\n\n")
cat("Obs\n",caln_nt$x,"\n\nPred\n",caln_nt$y,"\n\nLower\n",caln_nt$lower,"\n\nUpper\n",caln_nt$upper)
cat("\n\n\n")
cat(sink_marker)
cat("cal_t="); dput(cal_t); cat("\n\n")
cat("caln_t="); dput(caln_t); cat("\n\n")
cat("cal_nt="); dput(cal_nt); cat("\n\n")
cat("caln_nt="); dput(caln_nt); cat("\n\n")
sink()


sink(paste0(plot_dir,"topics/topic_cal_dif_test.txt"))

# Test: are calibration differences closer to 0 for topic-included models?
# Null hypothesis: calibration difference is equally poor in each bin for topics included/not included.
caldif_t=abs(c(cal_t$y-cal_t$x,cal_nt$y-cal_nt$x))
caldif_nt=abs(c(caln_t$y-caln_t$x,caln_nt$y-caln_nt$x))
print(wilcox.test(caldif_t,caldif_nt))

sink()


######################################################################################
######################################################################################
## Admission causes. Proportion of 'admissions' which were actually deaths in each  ##
##  age group                                                                       ##
######################################################################################
######################################################################################

# Age cutoffs: aggregate groups
age_cuts=c(0,5,20 + 5*(0:14))

# 5-year age brackets 20-90, 10 SIMD deciles
n_tot=matrix(0,length(age_cuts),10)
n_d=n_tot
n_b=n_tot


for (a in 1:(dim(n_tot)[1]-1)) {
  for (s in 1:dim(n_tot)[2]) {
    sub=which((all_pred$age >= age_cuts[a]) & (all_pred$age < age_cuts[a+1]) & (all_pred$simd==s))
    n_tot[a,s]=sum(all_pred$target[sub])
    n_d[a,s]=length(which((all_pred$target[sub]==1) & (all_pred$reason[sub]=="D")))
    n_b[a,s]=length(which((all_pred$target[sub]==1) & (all_pred$reason[sub]=="B")))
  }
  print(a)
}
# >=90
amax=dim(n_tot)[1]
for (s in 1:dim(n_tot)[2]) {
  sub=which((all_pred$age >= age_cuts[amax]) & (all_pred$simd==s))
  n_tot[amax,s]=sum(all_pred$target[sub])
  n_d[amax,s]=length(which((all_pred$target[sub]==1) & (all_pred$reason[sub]=="D")))
  n_b[amax,s]=length(which((all_pred$target[sub]==1) & (all_pred$reason[sub]=="B")))
}


# IMPORTANT: privacy
n_tot[which(n_tot<5)]=0
n_d[which(n_d<5)]=0
n_b[which(n_b<5)]=0


# Column and row names
rownames(n_tot)=c(paste0(age_cuts[1:(dim(n_tot)[1]-1)],"_",age_cuts[2:dim(n_tot)[1]]),">90"); 
colnames(n_tot)=paste0("SIMD",1:10)
rownames(n_b)=c(paste0(age_cuts[1:(dim(n_tot)[1]-1)],"_",age_cuts[2:dim(n_tot)[1]]),">90"); 
colnames(n_d)=paste0("SIMD",1:10)
rownames(n_d)=c(paste0(age_cuts[1:(dim(n_tot)[1]-1)],"_",age_cuts[2:dim(n_tot)[1]]),">90"); 
colnames(n_b)=paste0("SIMD",1:10)


nt=rowSums(n_tot)[1:amax]
nd=rowSums(n_d)[1:amax]
nb=rowSums(n_b)[1:amax]
x=(c(age_cuts,100)+c(0,age_cuts))[2:(amax+1)]/2

# General plot
pdf(paste0(plot_dir,"Analytics/admission_type_by_age_cumulative.pdf"),width=5,height=5)
plot(0,type="n",xlim=range(x),ylim=c(0,1),xlab="Age",ylab="Cumulative proportion",xaxs="i",yaxs="i")
polygon(c(min(x),min(x),max(x),max(x)),c(0,1,1,0),col="gray",border=NA)
polygon(c(x,max(x),min(x)),c((nd+nb)/nt,0,0),col="red",border=NA)
polygon(c(x,max(x),min(x)),c((nd/nt),0,0),col="black",border=NA)
legend("topleft",c("Adm.","Death","Both"),pch=16,col=c("gray","black","red"))
dev.off()
sink(paste0(plot_dir,"Analytics/admission_type_by_age_cumulative.txt"))
cat("Number of total admissions by age and SIMD\n\n")
print(n_tot)
cat("\n\n\nDeaths prior to admission in target year by age and SIMD\n\n")
print(n_d)
cat("\n\n\nAdmission and death in target year by age and SIMD\n\n")
print(n_b)
cat("\n\n\n")
cat(sink_marker)
cat("n_tot="); dput(n_tot); cat("\n\n")
cat("n_d="); dput(n_d); cat("\n\n")
cat("n_b="); dput(n_b); cat("\n\n")
sink()





# Breakdown by low/high SIMD
nt3=cbind(rowSums(n_tot[,1:3]),rowSums(n_tot[,4:7]),rowSums(n_tot[,8:10]))[1:amax,]
nd3=cbind(rowSums(n_d[,1:3]),rowSums(n_d[,4:7]),rowSums(n_d[,8:10]))[1:amax,]
nb3=cbind(rowSums(n_b[,1:3]),rowSums(n_b[,4:7]),rowSums(n_b[,8:10]))[1:amax,]

pdf(paste0(plot_dir,"Analytics/admission_type_by_age_simd.pdf"),width=5,height=5)
plot(0,type="n",xlim=range(x),ylim=c(0,0.3),xlab="Age",ylab="Proportion",xaxs="i",yaxs="i")
for (i in 1:3) {
  lines(x,(nd3/nt3)[,i],col="black",lty=i)
  lines(x,(nb3/nt3)[,i],col="red",lty=i)
  #  lines(x,((nb3+nd3)/nt3)[,i],col="red",lty=i)
}
legend("topleft",c("SIMD 1-3","SIMD 4-7","SIMD 8-10",
                   "Death","Both"),
       ncol=1,col=c(rep("black",3),"black","red"),lty=c(1:3,NA,NA),pch=c(rep(NA,3),"|","|"))
dev.off()



######################################################################################
######################################################################################
## Performance. The following scripts evaluate the performance of SPARRAv4 and      ##
##  compare it against SPARRAv3.                                                    ##
######################################################################################
######################################################################################


######################################################
## Generate matrices for each constituent predictor ##
######################################################

prednames=grep("pred.",colnames(all_pred),value=TRUE)
for (i in 1:length(prednames)) {
  predmat = data.frame(id = 1:nrow(all_pred),
    time = 1,
    crossval_group = all_pred$cv,
    risk_score = all_pred[,prednames[i]],
    target = (all_pred$target==max(as.numeric(all_pred$target))))
  if (max(predmat$risk_score)>2) predmat$risk_score=predmat$risk_score/100
  
  assign(paste0(prednames[i],".pred"),predmat)
  rm(predmat)
  gc()
}

super.pred = data.frame(id = 1:nrow(all_pred),
  time = 1,
  crossval_group = all_pred$cv,
  risk_score = all_pred$super,
  target = (all_pred$target==max(as.numeric(all_pred$target))))

longprednames=case_when(
  prednames=="pred.ann.h2o" ~"ANN (deep learning)",
  prednames=="pred.glm.h2o"~"GLM, elastic net",
  prednames=="pred.nb.h2o"~"Naive Bayes",
  prednames=="pred.rf1.h2o"~"RF, max. depth 20", 
  prednames=="pred.rf2.h2o"~"RF, max. depth 40",
  prednames=="pred.v3"~"SPARRA version 3",
  prednames=="pred.xgb1.xgb"~"XG-boost, max. depth 4",
  prednames=="pred.xgb2.xgb"~"XG-boost, max. depth 8",
  prednames=="pred.xgb3.xgb"~"XG-boost, max. depth 3")



######################################################
## Draw ROC and CAL plots for each predictor        ##
######################################################

# Note - predictors include RECONSTRUCTED v3 but not ACTUAL v3. Comparison with v3 is done later, and can't be compared directly with these, since some samples don't have a valid v3 score.

# Separate
for (i in 1:length(prednames)) {
  pred=get(paste0(prednames[i],".pred"))
 
  # ROC 
  pdf(paste0(plot_dir,"Performance/constituent_predictors/roc_",prednames[i],".pdf"),width=3,height=3.5)
  xroc=getroc(pred$target,pred$risk_score,cv=pred$crossval_group)
  lab=paste0(1:3,"; ",signif(xroc$auc[1:3],digits=3)," (",signif(xroc$se[1:3],digits=2),")")
  roc_2panel(xroc,labels=lab,col=c("black","red","blue"),title="CV fold; AUROC (SE)")
  dev.off()
  
  sink(paste0(plot_dir,"Performance/constituent_predictors/roc_",prednames[i],".txt"))
  cat("Sensitivity\n",xroc$sens,"\n\nSpecificity\n",xroc$spec,"\n\nAUROCs\n",xroc$auc,"\n\nSE\n",xroc$se)
  cat(sink_marker) 
  cat("xroc=")
  dput(xroc)
  sink()
  
  # PRC
  pdf(paste0(plot_dir,"Performance/constituent_predictors/prc_",prednames[i],".pdf"),width=3,height=3.5)
  xprc=getprc(pred$target,pred$risk_score,cv=pred$crossval_group)
  lab=paste0(1:3,"; ",signif(xprc$auc[1:3],digits=3)," (",signif(xprc$se[1:3],digits=2),")")
  prc_2panel(xprc,labels=lab,col=c("black","red","blue"),title="CV fold; AUPRC (SE)")
  dev.off()

  sink(paste0(plot_dir,"Performance/constituent_predictors/prc_",prednames[i],".txt"))
  cat("Sensitivity\n",xprc$sens,"\n\nPPV\n",xprc$ppv,"\n\nAUPRCs\n",xprc$auc,"\n\nSE\n",xprc$se)
  cat(sink_marker) 
  cat("xprc=")
  dput(xprc)
  sink()
  

  
  # Calibration
  pdf(paste0(plot_dir,"Performance/constituent_predictors/cal_",prednames[i],".pdf"),width=3,height=3.5)
  for (f in 1:3) 
    assign(paste0("xcal",f),
      plotcal(pred$target[which(pred$crossval==f)],pred$risk_score[which(pred$crossval==f)],
        n=20,plot=FALSE,kernel=TRUE));
  cal_2panel(list(xcal1,xcal2,xcal3),labels=1:3,col=c("black","red","blue"),title="CV fold")
  dev.off()
  
  sink(paste0(plot_dir,"Performance/constituent_predictors/cal_",prednames[i],".txt"))
  cat("Obs\n",xcal1$x,"\n\nPred\n",xcal1$y,"\n\nLower\n",xcal1$lower,"\n\nUpper\n",xcal1$upper)
  cat("\n\nObs\n",xcal2$x,"\n\nPred\n",xcal2$y,"\n\nLower\n",xcal2$lower,"\n\nUpper\n",xcal2$upper)
  cat("\n\nObs\n",xcal3$x,"\n\nPred\n",xcal3$y,"\n\nLower\n",xcal3$lower,"\n\nUpper\n",xcal3$upper)
  cat(sink_marker) 
  cat("xcal1="); dput(xcal1); cat("\n\n")
  cat("xcal2="); dput(xcal2); cat("\n\n")
  cat("xcal3="); dput(xcal3); cat("\n\n")
  sink()
  
  
  # Other metrics  
  sink(paste0(plot_dir,"Performance/constituent_predictors/metrics_",prednames[i],".txt"))
  performance_metrics=get_performance_metrics(pred)
  cat("Performance metrics:\n\n")
  print(performance_metrics)
  cat(sink_marker)
  cat("performance_metrics=")
  dput(performance_metrics)
  sink()
}


pred=super.pred

# ROC 
pdf(paste0(plot_dir,"Performance/v4_final/roc_super.pdf"),width=3,height=3.5)
xroc=getroc(pred$target,pred$risk_score,cv=pred$crossval_group)
lab=paste0(1:3,"; ",signif(xroc$auc[1:3],digits=3)," (",signif(xroc$se[1:3],digits=2),")")
roc_2panel(xroc,labels=lab,col=c("black","red","blue"),title="CV fold; AUROC (SE)")
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/roc_super.txt"))
cat("Sensitivity\n",xroc$sens,"\n\nSpecificity\n",xroc$spec,"\n\nAUROCs\n",xroc$auc,"\n\nSE\n",xroc$se)
cat(sink_marker)
cat("xroc=")
dput(xroc)
sink()

# PRC
pdf(paste0(plot_dir,"Performance/v4_final/prc_super.pdf"),width=3,height=3.5)
xprc=getprc(pred$target,pred$risk_score,cv=pred$crossval_group)
lab=paste0(1:3,"; ",signif(xprc$auc[1:3],digits=3)," (",signif(xprc$se[1:3],digits=2),")")
prc_2panel(xprc,labels=lab,col=c("black","red","blue"),title="CV fold; AUPRC (SE)")
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/prc_super.txt"))
cat("Sensitivity\n",xprc$sens,"\n\nPPV\n",xprc$ppv,"\n\nAUPRCs\n",xprc$auc,"\n\nSE\n",xprc$se)
cat(sink_marker)
cat("xprc=")
dput(xprc)
sink()



# Calibration
pdf(paste0(plot_dir,"Performance/v4_final/cal_super.pdf"),width=3,height=3.5)
for (f in 1:3) 
  assign(paste0("xcal",f),
    plotcal(pred$target[which(pred$crossval==f)],pred$risk_score[which(pred$crossval==f)],
      n=20,plot=FALSE,kernel=TRUE));
cal_2panel(list(xcal1,xcal2,xcal3),labels=1:3,col=c("black","red","blue"),title="CV fold")
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/cal_super.txt"))
cat("Obs\n",xcal1$x,"\n\nPred\n",xcal1$y,"\n\nLower\n",xcal1$lower,"\n\nUpper\n",xcal1$upper)
cat("\n\nObs\n",xcal2$x,"\n\nPred\n",xcal2$y,"\n\nLower\n",xcal2$lower,"\n\nUpper\n",xcal2$upper)
cat("\n\nObs\n",xcal3$x,"\n\nPred\n",xcal3$y,"\n\nLower\n",xcal3$lower,"\n\nUpper\n",xcal3$upper)
cat(sink_marker) 
cat("xcal1="); dput(xcal1); cat("\n\n")
cat("xcal2="); dput(xcal2); cat("\n\n")
cat("xcal3="); dput(xcal3); cat("\n\n")
sink()


# Other metrics  
sink(paste0(plot_dir,"Performance/v4_final/metrics_super.txt"))
get_performance_metrics(pred)
performance_metrics=get_performance_metrics(pred)
cat("Performance metrics:\n\n")
print(performance_metrics)
cat(sink_marker)
cat("performance_metrics=")
dput(performance_metrics)
sink()




######################################################
## Draw performance plots combined                  ##
######################################################

pdf(paste0(plot_dir,"Performance/constituent_predictors/roc_all.pdf"),width=6,height=7)
xroc0=getroc(super.pred$target,super.pred$risk_score)
for (i in 1:length(prednames)) {
  pred=get(paste0(prednames[i],".pred"))
  assign(paste0("xroc",i),getroc(pred$target,pred$risk_score))
}
labs=c("Final score",longprednames)
xcol=rep(c("black","blue","red"),10)[1:length(labs)]; xty=rep(1:10,each=3)[1:length(labs)];
roc_2panel(list(xroc0,xroc1,xroc2,xroc3,xroc4,xroc5,xroc6,xroc7,xroc8,xroc9),labels=labs,
  col=xcol,lty=xty,xy_col="gray",xy_lty=2)
dev.off()

sink(paste0(plot_dir,"Performance/constituent_predictors/roc_all.txt"))
for (i in 1:length(prednames)) {
  cat(prednames[i],"\n\n")
  xrocn=get(paste0("xroc",i))
  cat("Sensitivity\n",xrocn$sens,"\n\nSpecificity\n",xrocn$spec,"\n\nAUROCs\n",xrocn$auc,"\n\nSE\n",xrocn$se,"\n\n\n")
}
cat("Super\n\n")
cat("Sensitivity\n",xroc$sens,"\n\nSpecificity\n",xroc$spec,"\n\nAUROCs\n",xroc$auc,"\n\nSE\n",xroc$se)
cat("\n\n")
cat(sink_marker)
for(i in 1:length(prednames)) {
  xrocn=get(paste0("xroc",i))
  cat(paste0("xroc",i,"="))
  dput(xrocn)
  cat("\n\n")
}
cat("xroc0=")
dput(xroc0)
sink()


pdf(paste0(plot_dir,"Performance/constituent_predictors/prc_all.pdf"),width=6,height=7)
xprc0=getprc(super.pred$target,super.pred$risk_score)
for (i in 1:length(prednames)) {
  pred=get(paste0(prednames[i],".pred"))
  assign(paste0("xprc",i),getprc(pred$target,pred$risk_score))
}
labs=c("Final score",longprednames)
xcol=rep(c("black","blue","red"),10)[1:length(labs)]; xty=rep(1:10,each=3)[1:length(labs)];
prc_2panel(list(xprc0,xprc1,xprc2,xprc3,xprc4,xprc5,xprc6,xprc7,xprc8,xprc9),labels=labs,
  col=xcol,lty=xty)
dev.off()

sink(paste0(plot_dir,"Performance/constituent_predictors/prc_all.txt"))
for (i in 1:length(prednames)) {
  cat(prednames[i],"\n\n")
  xprcn=get(paste0("xprc",i))
  cat("Sensitivity\n",xprcn$sens,"\n\nPPV\n",xprcn$ppv,"\n\nAUPRCs\n",xprcn$auc,"\n\nSE\n",xprcn$se)
}
cat("Super\n\n")
cat("Sensitivity\n",xprc$sens,"\n\nPPV\n",xprc$ppv,"\n\nAUPRCs\n",xprc$auc,"\n\nSE\n",xprc$se)
cat(sink_marker)
for(i in 1:length(prednames)) {
  xprcn=get(paste0("xprc",i))
  cat(paste0("xprc",i,"="))
  dput(xprcn)
  cat("\n\n")
}
cat("xprc0=")
dput(xprc0)
sink()



pdf(paste0(plot_dir,"Performance/constituent_predictors/cal_all.pdf"),width=6,height=7)
for (i in 1:length(prednames)) {
  pred=get(paste0(prednames[i],".pred"))
  assign(paste0("xcal",i),plotcal(pred$target,pred$risk_score,plot=FALSE,kernel=TRUE))
}
xcal0=plotcal(super.pred$target,super.pred$risk_score,kernel=TRUE,plot=FALSE)
labs=c("Final score",longprednames)
xcol=rep(c("black","blue","red"),10)[1:length(labs)]; xty=rep(1:10,each=3)[1:length(labs)];
cal_2panel(list(xcal0,xcal1,xcal2,xcal3,xcal4,xcal5,xcal6,xcal7,xcal8,xcal9),col=xcol,
  lty=xty,xy_col="gray",xy_lty=2,labels=labs)
dev.off()

sink(paste0(plot_dir,"Performance/constituent_predictors/cal_all.txt"))
for (i in 1:length(prednames)) {
  cat(prednames[i],"\n\n")
  xcaln=get(paste0("xcal",i))
  cat("Obs\n",xcaln$x,"\n\nPred\n",xcaln$y,"\n\nLower\n",xcaln$lower,"\n\nUpper\n",xcaln$upper)
}
cat("Super\n\n")
cat("Obs\n",xcal0$x,"\n\nPred\n",xcal0$y,"\n\nLower\n",xcal0$lower,"\n\nUpper\n",xcal0$upper)
cat(sink_marker)
for(i in 1:length(prednames)) {
  xcaln=get(paste0("xcal",i))
  cat(paste0("xcal",i,"="))
  dput(xcaln)
  cat("\n\n")
}
cat("xcal0=")
dput(xcal0)
sink()








######################################################
## Maximum of v3 and v4                             ##
######################################################

# We will compare with v3 and v4
target=(all_pred$target==max(all_pred$target))
v4=all_pred$super
v3=all_pred$v3
cv=all_pred$cv

max34=pmax(v3,v4)

# ROC 
xroc3=getroc(target,v3,cv=cv)
xroc3b=getroc(target,v3)
xroc4=getroc(target,v4,cv=cv)
xroc4b=getroc(target,v4)
xrocm=getroc(target,max34,cv=cv)
xrocmb=getroc(target,max34)

pdf(paste0(plot_dir,"Performance/v4_final/roc_max.pdf"),width=3,height=3.5)
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2); 
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2); 
aucm=signif(xrocm$auc[4],digits=3); sem=signif(xrocm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
  paste0("Max: ",aucm," (",sem,")"))
roc_2panel(list(xroc3b,xroc4b,xrocmb),labels=labs,
  col=c("blue","black","red"),xy_lty=2,xy_col="gray",
  text.col=c("blue","red","black"),title="AUROC (SE)",
  title.col="black")
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/roc_max.txt"))
cat("V3\n")
cat("Sensitivity\n",xroc3b$sens,"\n\nSpecificity\n",xroc3b$spec,"\n\nAUROCs\n",xroc3b$auc,"\n\nSE\n",xroc3b$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xroc3$sens,"\n\nSpecificity\n",xroc3$spec,"\n\nAUROCs\n",xroc3$auc,"\n\nSE\n",xroc3$se)
cat("\n\nV4\n")
cat("Sensitivity\n",xroc4b$sens,"\n\nSpecificity\n",xroc4b$spec,"\n\nAUROCs\n",xroc4b$auc,"\n\nSE\n",xroc4b$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xroc4$sens,"\n\nSpecificity\n",xroc4$spec,"\n\nAUROCs\n",xroc4$auc,"\n\nSE\n",xroc4$se)
cat("\n\nMax\n")
cat("Sensitivity\n",xrocmb$sens,"\n\nSpecificity\n",xrocmb$spec,"\n\nAUROCs\n",xrocmb$auc,"\n\nSE\n",xrocmb$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xrocm$sens,"\n\nSpecificity\n",xrocm$spec,"\n\nAUROCs\n",xrocm$auc,"\n\nSE\n",xrocm$se)
cat(sink_marker)
cat("xroc3=")
dput(xroc3)
cat("\n\nxroc3b=")
dput(xroc3b)
cat("\n\nxroc4=")
dput(xroc4)
cat("\n\nxroc4b=")
dput(xroc4b)
cat("\n\nxrocm=")
dput(xrocm)
cat("\n\nxrocmb=")
dput(xrocmb)
sink()




# PRC
xprc3=getprc(target,v3,cv=cv)
xprc3b=getprc(target,v3)
xprc4=getprc(target,v4,cv=cv)
xprc4b=getprc(target,v4)
xprcm=getprc(target,max34,cv=cv)
xprcmb=getprc(target,max34)


pdf(paste0(plot_dir,"Performance/v4_final/prc_max.pdf"),width=3,height=3.5)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2); 
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2); 
aucm=signif(xprcm$auc[4],digits=3); sem=signif(xprcm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
  paste0("Max: ",aucm," (",sem,")"))
prc_2panel(list(xprc3b,xprc4b,xprcmb),labels=labs,col=c("blue","red","black"),
  title="AUPRC (SE)",title.col="black",text.col=c("blue","red","black"))
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/prc_max.txt"))
cat("v3\n")
cat("Sensitivity\n",xprc3b$sens,"\n\nPPV\n",xprc3b$ppv,"\n\nAUPRCs\n",xprc3b$auc,"\n\nSE\n",xprc3b$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xprc3$sens,"\n\nPPV\n",xprc3$ppv,"\n\nAUPRCs\n",xprc3$auc,"\n\nSE\n",xprc3$se)
cat("\n\nV4\n")
cat("Sensitivity\n",xprc4b$sens,"\n\nPPV\n",xprc4b$ppv,"\n\nAUPRCs\n",xprc4b$auc,"\n\nSE\n",xprc4b$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xprc4$sens,"\n\nPPV\n",xprc4$ppv,"\n\nAUPRCs\n",xprc4$auc,"\n\nSE\n",xprc4$se)
cat("\n\nMax\n")
cat("Sensitivity\n",xprcmb$sens,"\n\nPPV\n",xprcmb$ppv,"\n\nAUPRCs\n",xprcmb$auc,"\n\nSE\n",xprcmb$se)
cat("\n\nXV\n")
cat("Sensitivity\n",xprcm$sens,"\n\nPPV\n",xprcm$ppv,"\n\nAUPRCs\n",xprcm$auc,"\n\nSE\n",xprcm$se)
cat(sink_marker)
cat("xprc3="); dput(xprc3); cat("\n\n")
cat("xprc3b="); dput(xprc3b); cat("\n\n")
cat("xprc4="); dput(xprc4); cat("\n\n")
cat("xprc4b="); dput(xprc4b); cat("\n\n")
cat("xprcm="); dput(xprcm); cat("\n\n")
cat("xprcmb="); dput(xprcmb); cat("\n\n")
sink()






# Calibration
xcal3=plotcal(target,v3,n=20,plot=FALSE,kernel=TRUE);
xcal4=plotcal(target,v4,n=20,plot=FALSE,kernel=TRUE);
xcalm=plotcal(target,max34,n=20,plot=FALSE,kernel=TRUE);

pdf(paste0(plot_dir,"Performance/v4_final/cal_max.pdf"),width=3,height=3.5)
labs=c("v3", "v4","Max")
cic=c(rgb(0,0,1,alpha=0.5),rgb(1,0,1,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(xcal3,xcal4,xcalm),labels=labs,col=c("blue","red","black"),
  text.col=c("blue","red","black"),ci_col=cic,xy_col="gray",xy_lty=2)
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/cal_max.txt"))
cat("v3\n")
cat("Obs\n",xcal3$x,"\n\nPred\n",xcal3$y,"\n\nLower\n",xcal3$lower,"\n\nUpper\n",xcal3$upper)
cat("\n\nV4\n")
cat("Obs\n",xcal4$x,"\n\nPred\n",xcal4$y,"\n\nLower\n",xcal4$lower,"\n\nUpper\n",xcal4$upper)
cat("\n\nMax\n")
cat("Obs\n",xcalm$x,"\n\nPred\n",xcalm$y,"\n\nLower\n",xcalm$lower,"\n\nUpper\n",xcalm$upper)
cat(sink_marker)
cat("xcal3="); dput(xcal3); cat("\n\n")
cat("xcal4="); dput(xcal4); cat("\n\n")
cat("xcalm="); dput(xcalm); cat("\n\n")
sink()


max.pred= data.frame(id = 1:nrow(all_pred),
    time = 1,
    crossval_group = all_pred$cv,
    risk_score = pmax(all_pred$super,all_pred$v3),
    target = (all_pred$target==max(as.numeric(all_pred$target))))


# Other metrics  
sink(paste0(plot_dir,"Performance/v4_final/metrics_max.txt"))
performance_metrics=get_performance_metrics(max.pred)
cat("Performance metrics:\n\n")
print(performance_metrics)
cat(sink_marker)
cat("performance_metrics=")
dput(performance_metrics)
sink()



######################################################
## Density plot v3 vs v4                            ##
######################################################

# Kernel density estimates normalised to unit diagonal
xcol=colorRampPalette(c("white","darkgray","red","yellow"))(100)
xx=1; nn=100; sp=7; scmax=100
qplot= ((1:1000)/1000)^(1/4) # Add q-q plots at these quantiles. >500 samples in each bin.
library(matrixStats)


pi=all_pred$v3; pj=all_pred$super;  

w=which(pi<scmax/100 & pj<scmax/100)

pi=pi[w]; pj=pj[w]
x=round(nn*pi); y=round(nn*pj)
xout=xtabs(~xi+xj,cbind(xi=c(x,1:scmax),xj=c(y,1:scmax)))
xout[which(xout<4)]=0; xout[which(xout>3 & xout<7)]=6 # DO NOT DELETE: important for privacy
xout0=xout

# Marginals
tab1=table(round(pj*100)/100); d1=list(x=c(0,as.numeric(names(tab1))),y=c(0,tab1/sum(tab1))); 
tab2=table(round(pj*100)/100); d2=list(x=c(0,as.numeric(names(tab2))),y=c(0,tab2/sum(tab2)));

# Medians
q0=quantile(pi,qplot)
q1=quantile(pj,qplot)
q0a=1:scmax; q1a=q0a; for (ii in 1:scmax) q1a[ii]=median(y[which(x==ii)],na.rm=T); q0a=q0a/scmax; q1a=q1a/scmax
q0b=1:scmax; q1b=q0b; for (ii in 1:scmax) q1b[ii]=median(x[which(y==ii)],na.rm=T); q0b=q0b/scmax; q1b=q1b/scmax



pdf(paste0(plot_dir,"Performance/v4_final/v3_v4_density.pdf"),width=7,height=7)
xout=xout0
# Normalise columns, then rows
for (ii in 1:dim(xout)[2]) xout[ii,]=xout[ii,]/sum(xout[ii,],na.rm=T);
for (ii in 1:dim(xout)[1]) xout[,ii]=xout[,ii]/sum(xout[,ii],na.rm=T);  
xmax=min(colMaxs(xout[,10:(scmax-15)],na.rm=T));
xmin=min(xout[which(xout>0)],na.rm=T); ncol0=5

image((1:scmax)/100,(1:scmax)/100,xout,xlab=paste0("SPARRA v3"),
  ylab=paste0("SPARRA v4"),
  col=xcol,breaks=c(0,rep(xmin/2,ncol0),seq(xmin,4*xmax/3,length=length(xcol)-ncol0-1),1+max(xout,na.rm=T)),
  xaxs="i",yaxs="i",xlim=c(-1/sp,scmax/100),ylim=c(-1/sp,scmax/100),bty="n")

polygon(d1$x,d1$y/(sp*max(d1$y))-(1/sp),col="gray"); 
polygon(d2$y/(sp*max(d2$y))-(1/sp),d2$x,col="gray")

abline(0,1,col="black",lty=2,lwd=2)
#abline(0,1,col="white",lty=1,lwd=2)

#lines(q0,q1,col="black",lwd=4)
lines(q0,q1,col="blue",lwd=2)

#lines(q0a,q1a,col="black",lwd=4)
lines(q0a,q1a,col="darkgreen",lwd=2)
lines(q1b,q0b,col="darkgreen",lwd=2,lty=2)

legend("bottomright",c("Density (low)","", "Density (high)", "Marginal", "X-Y","Q-Q","Median (col)","Median (row)"),
  pch=c(16,16,16,NA,NA,NA,NA,NA),lty=c(NA,NA,NA,1,2,1,1,2),lwd=c(NA,NA,NA,1,2,2,2,2),
  col=c("darkgray","red","yellow","black","black","blue","darkgreen","darkgreen"))
dev.off()

sink(paste0(plot_dir,"Performance/v4_final/v3_v4_density.txt"))
options(width=2000)
print(xout0)
cat("\n\n")
print(d1$x)
cat("\n\n")
print(d1$y)
cat("\n\n\n")
print(d2$x)
cat("\n\n")
print(d2$y)
cat("\n\n")
print(q0)
cat("\n\n")
print(q1)     
cat("\n\n")
print(q0a)
cat("\n\n")
print(q1a)
cat(sink_marker)
cat("xout="); dput(xout); cat("\n\n")
cat("d1="); dput(d1); cat("\n\n")
cat("d2="); dput(d2); cat("\n\n")
cat("q0="); dput(q0); cat("\n\n")
cat("q1="); dput(q1); cat("\n\n")
cat("q0a="); dput(q0a); cat("\n\n")
cat("q1a="); dput(q1a); cat("\n\n")
cat("q0b="); dput(q0b); cat("\n\n")
cat("q1b="); dput(q1b); cat("\n\n")
sink()




######################################################
## SPARRA v3 vs v4 in various subcohorts            ##
######################################################

super=all_pred$super
v3=all_pred$v3
target=(all_pred$target==max(all_pred$target))
cv=all_pred$cv

times=sort(unique(all_pred$time))
  
### Definition of some subcohorts
subcohorts=list(
  which(all_pred$age>80), ### High risk: age >80
  which(all_pred$age>20 & all_pred$age<70 & all_pred$ae2<1), ### Low risk: age 20-70, no previous A&E admissions
  setdiff(1:dim(all_pred)[1],c(which(all_pred$age>80),which(all_pred$age>20 & all_pred$age<70 & all_pred$ae2<1))), # complement of previous two
  which(all_pred$time==times[1]), # final year
  which(all_pred$time==times[2]), # middle year
  which(all_pred$time==times[3]) # final year
)


### Names to index subcohort analysis files
names=c("age_gt80","age_20_70_noAE","complement","final_year","middle_year","first_year")

### Verbose names to appear on plots
vnames=c("Age > 80","Age 20-70, no prev. A&E","Mid-risk","May 2017 cutoff","Jan 2017 cutoff","May 2016 cutoff")



for (cht in 1:length(subcohorts)) {
  
sub=subcohorts[[cht]]
suffix=names[cht]


xroc3=getroc(target[sub],v3[sub],cv=cv[sub])
xroc3b=getroc(target[sub],v3[sub])
xroc4=getroc(target[sub],super[sub],cv=cv[sub]) 
xroc4b=getroc(target[sub],super[sub]) # Ignore cv structure when plotting for visual simplicity

pdf(paste0(plot_dir,"Performance/v3_vs_v4/roc_v3_vs_v4_",suffix,".pdf"),width=3,height=3.5)
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2); 
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2); 
labs=c(paste0("v3: ",auc3," (",se3,")"),  paste0("v4: ",auc4," (",se4,")"))
roc_2panel(list(xroc3b,xroc4b),labels=labs,col=c("blue","red"),xy_col="gray",xy_lty=2,
  title="AUROC (SE)",title.col="black",text.col=c("blue","red"))
dev.off()

sink(paste0(plot_dir,"Performance/v3_vs_v4/roc_v3_vs_v4_",suffix,".txt"))
cat("v3_byfold","\n")
cat("Sensitivity\n",xroc3$sens,"\n\nSpecificity\n",xroc3$spec,"\n\nAUROCs\n",xroc3$auc,"\n\nSE\n",xroc3$se,"\n\n\n")
cat("\n\n")
cat("v3","\n")
cat("Sensitivity\n",xroc3b$sens,"\n\nSpecificity\n",xroc3b$spec,"\n\nAUROCs\n",xroc3b$auc,"\n\nSE\n",xroc3b$se,"\n\n\n")
cat("\n\n")
cat("v4_byfold","\n")
cat("Sensitivity\n",xroc4$sens,"\n\nSpecificity\n",xroc4$spec,"\n\nAUROCs\n",xroc4$auc,"\n\nSE\n",xroc4$se,"\n\n\n")
cat("\n\n")
cat("v4_combined","\n")
cat("Sensitivity\n",xroc4b$sens,"\n\nSpecificity\n",xroc4b$spec,"\n\nAUROCs\n",xroc4b$auc,"\n\nSE\n",xroc4b$se,"\n\n\n")
cat(sink_marker)
cat("xroc3=")
dput(xroc3)
cat("\n\n")
cat("xroc3b=")
dput(xroc3b)
cat("\n\n")
cat("xroc4=")
dput(xroc4)
cat("\n\n")
cat("xroc4b=")
dput(xroc4b)
cat("\n\n")
sink()


### PR curves

xprc3=getprc(target[sub],v3[sub],cv=cv[sub])
xprc3b=getprc(target[sub],v3[sub])
xprc4=getprc(target[sub],super[sub],cv=cv[sub]) 
xprc4b=getprc(target[sub],super[sub]) # Ignore cv structure when plotting for visual simplicity

pdf(paste0(plot_dir,"Performance/v3_vs_v4/prc_v3_vs_v4_",suffix,".pdf"),width=3,height=3.5)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2); 
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2); 
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"))
prc_2panel(list(xprc3b,xprc4b),labels=labs,col=c("blue","red"),
  title="AUPRC (SE)",title.col="black",text.col=c("blue","red"))
dev.off()

sink(paste0(plot_dir,"Performance/v3_vs_v4/prc_v3_vs_v4_",suffix,".txt"))
cat("v3_byfold","\n")
cat("Sensitivity\n",xprc3$sens,"\n\nPPV\n",xprc3$ppv,"\n\nAUROCs\n",xprc3$auc,"\n\nSE\n",xprc3$se,"\n\n\n")
cat("\n\n")
cat("v3","\n")
cat("Sensitivity\n",xprc3b$sens,"\n\nPPV\n",xprc3b$ppv,"\n\nAUROCs\n",xprc3b$auc,"\n\nSE\n",xprc3b$se,"\n\n\n")
cat("\n\n")
cat("v4_byfold","\n")
cat("Sensitivity\n",xprc4$sens,"\n\nPPV\n",xprc4$ppv,"\n\nAUROCs\n",xprc4$auc,"\n\nSE\n",xprc4$se,"\n\n\n")
cat("\n\n")
cat("v4_combined","\n")
cat("Sensitivity\n",xprc4b$sens,"\n\nPPV\n",xprc4b$ppv,"\n\nAUROCs\n",xprc4b$auc,"\n\nSE\n",xprc4b$se,"\n\n\n")
cat(sink_marker)
cat("xprc3=")
dput(xprc3)
cat("\n\n")
cat("xprc3b=")
dput(xprc3b)
cat("\n\n")
cat("xprc4=")
dput(xprc4)
cat("\n\n")
cat("xprc4b=")
dput(xprc4b)
cat("\n\n")
sink()




### Calibration curves
xcal3=plotcal(target[sub],v3[sub],plot=FALSE,kernel=TRUE)
xcal4=plotcal(target[sub],super[sub],plot=FALSE,kernel=TRUE)


pdf(paste0(plot_dir,"Performance/v3_vs_v4/cal_v3_vs_v4_",suffix,".pdf"),width=3,height=3.5)
cal_2panel(list(xcal3,xcal4),label=c("v3","v4"),col=c("blue","red"),xy_col="gray",xy_lty=2,
  ci_col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,1,alpha=0.5)))
dev.off()

sink(paste0(plot_dir,"Performance/v3_vs_v4/cal_v3_vs_v4_",suffix,".txt"))
cat("v3\n\n")
cat("Obs\n",xcal3$x,"\n\nPred\n",xcal3$y,"\n\nLower\n",xcal3$lower,"\n\nUpper\n",xcal3$upper)
cat("\n\n\nv4_combined\n\n")
cat("Obs\n",xcal3$x,"\n\nPred\n",xcal3$y,"\n\nLower\n",xcal3$lower,"\n\nUpper\n",xcal3$upper)
cat(sink_marker)
cat("xcal3=")
dput(xcal3)
cat("\n\n")
cat("xcal4=")
dput(xcal4)
cat("\n\n")
sink()


}




######################################################
## SPARRA v3 vs v4 - other plots                    ##
######################################################

super=all_pred$super
v3=all_pred$v3
target=(all_pred$target==max(all_pred$target))

topn=1000*(1:50)
tp_super=0*topn
tp_v3=0*topn
cut_super=0*topn
cut_v3=0*topn
meancut_super=0*topn
meancut_v3=0*topn


for (i in 1:length(topn)) {
  cut_super[i]=-sort(-super)[topn[i]]
  cut_v3[i]=-sort(-v3)[topn[i]]
  
  w_v3=which(v3>cut_v3[i])
  w_super=which(super>cut_super[i])
  
  tp_v3[i]=sum(target[w_v3])
  tp_super[i]=sum(target[w_super])
  
  meancut_v3[i]=mean(v3[w_v3])
  meancut_super[i]=mean(super[w_super])
}


# Save file
pdf(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_num_in_topN.pdf"),width=7,height=4)

#par(mfrow=c(2,1))

par(mar=c(0.5,4,4,2))
plot(0,0,type="n",ylim=c(0,max(c(tp_v3,tp_super))),xlim=range(topn),#main="Admissions amongst highest-risk patients by v3 and v4",
  xaxt="n",xlab="",ylab="Number admitted")
points(topn,tp_v3,pch=16,col="blue")
points(topn,tp_super,pch=16,col="red")
lines(topn,topn*meancut_v3,col="blue")
lines(topn,topn*meancut_super,col="red")
legend("bottomright",c("v3","v4"),pch=16,col=c("blue","red"))

dev.off()


sink(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_num_in_topN.txt"))
options(width=3000)
print("The following table shows: highest risk [n] patients by super-learner or SPARRAv3 followed by the number of EDAs in those cohorts of patients and risk score cutoffs and the mean risk scores amongst those cohorts of patients")
print(rbind(topn,tp_super,tp_v3,cut_super,cut_v3, meancut_super,meancut_v3))
cat("\n\n")
cat(sink_marker)
cat("topn="); dput(topn); cat("\n\n")
cat("tp_super="); dput(tp_super); cat("\n\n")
cat("tp_v3="); dput(tp_v3); cat("\n\n")
cat("cut_super="); dput(cut_super); cat("\n\n")
cat("cut_v3="); dput(cut_v3); cat("\n\n")
cat("meancut_super="); dput(meancut_super); cat("\n\n")
cat("meancut_v3="); dput(meancut_v3); cat("\n\n")
sink()






pdf(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_num_in_topN_difference.pdf"),width=7,height=4)

par(mar=c(4,4,0.5,2))
plot(0,0,type="n",ylim=c(0,max(tp_super-tp_v3)),xlim=range(topn),
  xlab="Number of patients",ylab="Difference")
points(topn,tp_super-tp_v3,pch=16,col="black")
lines(loess.smooth(topn,abs(topn*meancut_v3-tp_v3)), col="blue")
lines(topn,abs(topn*meancut_super-tp_super), col="red")
if(FALSE) segments(topn,tp_super-tp_v3,topn,0*topn,lty=2,col="black")

legend("topleft",c("v3","v4"),lty=1,col=c("blue","red"))

dev.off()


sink(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_num_in_topN_difference.txt"))
options(width=3000)
print("The following table shows: highest risk [n] patients by super-learner or SPARRAv3 followed by the number of EDAs in those cohorts of patients and risk score cutoffs and the mean risk scores amongst those cohorts of patients")
print(rbind(topn,tp_super,tp_v3,cut_super,cut_v3, meancut_super,meancut_v3))
cat("\n\n")
cat(sink_marker)
cat("topn="); dput(topn); cat("\n\n")
cat("tp_super="); dput(tp_super); cat("\n\n")
cat("tp_v3="); dput(tp_v3); cat("\n\n")
cat("cut_super="); dput(cut_super); cat("\n\n")
cat("cut_v3="); dput(cut_v3); cat("\n\n")
cat("meancut_super="); dput(meancut_super); cat("\n\n")
cat("meancut_v3="); dput(meancut_v3); cat("\n\n")
sink()





# Number of patients needed to treat to target a number of avoidable admissions. This comparison requires
#  the same patient set for both v3 and v4, so we use the representative subsample derived earlier.

nadm=500*(1:10) # in order to prevent this many preventable admissions
tadm=round(nadm/0.20) # screen this total number of admissions

s_v3=0*nadm #  screen this many patients using v3
s_super=0*nadm # screen this many patients using the super learner

super=all_pred$super
v3=all_pred$v3
target=(all_pred$target==max(all_pred$target))

rv3=-sort(-v3[which(target)])
rsuper=-sort(-super[which(target)])

for (i in 1:length(nadm)) {
  cv3=rv3[tadm[i]]
  tadm2=length(which(v3>cv3 & target==1))
  csuper=rsuper[1+tadm2]
  s_v3[i]=length(which(v3>cv3))
  s_super[i]=length(which(super>csuper))
}


pdf(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_ntreat.pdf"),width=7,height=4)

plot(nadm,s_v3-s_super,col="red", #main="Number of patients needed to treat to target avoidable admissions",
  xlab="To target this many avoidable admissions",ylab="Treat this many fewer patients using v4",pch=16)

dev.off()

sink(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_ntreat.txt"))
print("The following table shows rows a),b),c) which indicate: to avoid a) admissions, treat the b) most-at-risk patients by v3, or c) most-at-risk patients by v4")
print(rbind(nadm,s_v3,s_super))
cat(sink_marker)
cat("nadm="); dput(nadm); cat("\n\n")
cat("s_v3="); dput(s_v3); cat("\n\n")
cat("s_super="); dput(s_super); cat("\n\n")
sink()



# Calibration of v3 and v4 amongst individuals for which score differs substantially

v4=all_pred$super
v3=all_pred$v3
target=(all_pred$target==max(all_pred$target))

del=0.1
h4=which(v4-v3>del)
h3=which(v3-v4>del)

cv3h3=plotcal(target[h3],v3[h3],kernel=T,plot=F)
cv4h3=plotcal(target[h3],v4[h3],kernel=T,plot=F)
cv3h4=plotcal(target[h4],v3[h4],kernel=T,plot=F)
cv4h4=plotcal(target[h4],v4[h4],kernel=T,plot=F)


pdf(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_calibration_differential.pdf"),width=3,height=3.5)
cic=c(rgb(1,0,0,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(cv3h3,cv3h4,cv4h3,cv4h4),labels=c("v3, v3>v4","v3, v4>v3","v4, v3>v4", "v4, v4>v3"),
  title=paste0("|v3-v4| > ",del),xy_col="gray",xy_lty=2,col=c("red","red","black","black"),
  lty=c(1,2,1,2),ci_col=cic)
dev.off()

sink(paste0(plot_dir,"Performance/v3_vs_v4/v3_vs_v4_calibration_differential.txt"))
print("Calibration amongst individuals with different v3 and v4 scores, |v3-v4|>0.1. In order: v3, v3>v4; v3, v4>v3; v4, v3>v4; v4, v4>v3\n")
cx=cv3h3; cat(cx$x,"\n",cx$y,"\n",cx$lower,"\n",cx$upper,"\n\n")
cx=cv3h4; cat(cx$x,"\n",cx$y,"\n",cx$lower,"\n",cx$upper,"\n\n")
cx=cv4h3; cat(cx$x,"\n",cx$y,"\n",cx$lower,"\n",cx$upper,"\n\n")
cx=cv4h4; cat(cx$x,"\n",cx$y,"\n",cx$lower,"\n",cx$upper,"\n\n")
cat(sink_marker)
cat("cv3h3="); dput(cv3h3); cat("\n\n")
cat("cv3h4="); dput(cv3h4); cat("\n\n")
cat("cv4h3="); dput(cv4h3); cat("\n\n")
cat("cv4h4="); dput(cv4h4); cat("\n\n")
sink()












######################################################################################
######################################################################################
## Analytics. The following scripts assess the behaviour of SPARRAv4 in terms of    ##
##  patient sub-cohorts.                                                            ##
######################################################################################
######################################################################################

######################################################
## Distribution of times-to-first-event             ##
######################################################

xpred=all_pred$super
xvec=all_pred$time_to_first_event

xmax=floor(10*max(xpred,na.rm=T))-1; mh=0
tmax=max(xvec,na.rm=T)
kx=function(x) density(c(-x,x,2*tmax - x),na.rm=T,from=1,to=tmax)
for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=kx(xvec[sub])
  mh=max(c(mh,dx$y))
  assign(paste0("d",i),dx)
}
d0=kx(xvec)


# We mostly predict imminent events
pdf(paste0(plot_dir,"Analytics/time_to_first_event.pdf"),width=5,height=5)

# To recreate: first set 
# x=c(1,1); tmax=1
xmax=8
mh=0; for (i in 1:xmax) mh=max(mh,get(paste0("d",i))$y)

plot(0,type="n",xlim=c(0,365),xlab="Days after time cutoff",ylim=c(0,mh), 
  ylab="Density", #main="Density of time-to-first-admission",
  yaxt="n")

gcol=gray((1:xmax)/(1.5*xmax))
for (i in 1:xmax) lines(get(paste0("d",i)),col=gcol[i],lty=1,lwd=2)
lines(d0,col="red",lwd=2)

legend("topright",lty=1,lwd=2,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"))

dev.off()

sink(paste0(plot_dir,"Analytics/time_to_first_event.txt"))
cat(d0$x,"\n",d0$y,"\n\n")
for (i in 1:xmax) {
  di=get(paste0("d",i))
  cat(di$x,"\n",di$y,"\n\n")
}
cat(sink_marker)
cat("d0="); dput(d0); cat("\n\n")
for (i in 1:xmax) {
  di=get(paste0("d",i))
  cat(paste0("d",i,"=")); dput(di); cat("\n\n")
  }
sink()




######################################################
## Distribution of number of events                 ##
######################################################

xpred=all_pred$super
xvar=all_pred$n_event

xmax=8; mh=0
nmax=min(which(table(xvar)<20)) -1 

xfreq=table(xvar); 
xfreq[nmax]=sum(xfreq[nmax:length(xfreq)])
xfreq=xfreq[1:nmax]; xfreq=100*xfreq/sum(xfreq)

for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=table(xvar[sub])
  if (length(dx)>nmax) {
    dx[nmax]=sum(dx[nmax:length(dx)])
    dx=dx[1:nmax]
  }
  dx=100*dx/sum(dx)
  mh=max(c(mh,dx))
  xfreq=rbind(xfreq,dx)
}


# We predict patients with multiple admissions better
pdf(paste0(plot_dir,"Analytics/n_events.pdf"),width=5,height=3)

xmax=8; mh=0
gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("red",gcol)

# Truncate for ease of viewing
xfreq1=xfreq
xfreq1[,6]=rowSums(xfreq[,6:dim(xfreq)[2]])
xfreq1=xfreq1[,1:6]

m0=dim(xfreq1)[1]; n0=dim(xfreq1)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq1)),bty="n",xaxt="n",
  xlab="N admissions",ylab="% patients")
for (i in 1:dim(xfreq1)[2]) 
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq1[,i],col=xcol,lty=1,lwd=1)
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=c(1:(n0-1),paste0(">",n0-1)))

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
  lwd=2)

dev.off()
sink(paste0(plot_dir,"Analytics/n_events.txt"))
cat("Frequencies:\n\n")
print(xfreq)
cat(sink_marker)
cat("xfreq="); dput(xfreq)
sink()




######################################################
## Distribution of times-to-only-event              ##
######################################################


# As for plot of time-to-first-event, but this time only including patients with one admission.
w=which(all_pred$n_event==1)
xpred=all_pred$super[w]
xvec=all_pred$time_to_first_event[w]

xmax=floor(10*max(xpred,na.rm=T))-1; mh=0
tmax=max(xvec,na.rm=T)
kx=function(x) density(c(-x,x,2*tmax - x),na.rm=T,from=1,to=tmax)
for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=kx(xvec[sub])
  mh=max(c(mh,dx$y))
  assign(paste0("d",i),dx)
}
d0=kx(xvec)


# We mostly predict imminent events
pdf(paste0(plot_dir,"Analytics/time_to_only_event.pdf"),width=5,height=5)

# To recreate: first set 
# x=c(1,1); tmax=1
xmax=8
mh=0; for (i in 1:xmax) mh=max(mh,get(paste0("d",i))$y)


plot(0,type="n",xlim=c(0,365),xlab="Days after time cutoff",ylim=c(0,mh), 
  ylab="Density",#main="Density of time-to-only-admission",
  yaxt="n")

gcol=gray((1:xmax)/(1.5*xmax))
for (i in 1:xmax) lines(get(paste0("d",i)),col=gcol[i],lty=1,lwd=2)
lines(d0,col="red",lwd=2)

legend("topright",lty=1,lwd=2,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"))

dev.off()
sink(paste0(plot_dir,"Analytics/time_to_only_event.txt"))
cat(paste0("All: \n"))
cat(d0$x,"\n",d0$y,"\n\n")
for (i in 1:xmax) {
  cat(paste0("Score bracket ",i,":")); cat("\n")
  di=get(paste0("d",i))
  cat(di$x,"\n",di$y,"\n\n")
}
cat(sink_marker)
cat(paste0("d0=")); dput(d0); cat("\n\n\n")
for (i in 1:xmax) {
  di=get(paste0("d",i))
  cat(paste0("d",i,"=")); dput(di); cat("\n\n")
}
sink()



######################################################
## Accuracy of predictions by diagnostic class      ##
######################################################

xpred=all_pred$super
xvar=as.factor(icd2sick(all_pred$main_condition)) # this function is found in auxiliary.R
xl0=levels(xvar)

xfreq=table(xvar);
nfreq=xfreq;
xfreq=100*xfreq/sum(xfreq)

for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=table(xvar[sub])
  dx=100*dx/sum(dx)
  mh=max(c(mh,dx))
  xfreq=rbind(xfreq,dx)
  nfreq=rbind(nfreq,table(xvar[sub]))
}

xfreq[which(nfreq<6)]=0


ox=order(-xfreq[1,])
xfreq=xfreq[,ox]; 

# Privacy
for (i in 1:dim(xfreq)[2]) {
  w=which(nfreq[,i]<5); xfreq[w,i]=5/sum(nfreq[,i]) # THIS IS IMPORTANT TO ENSURE NO PLOT ELEMENTS CONCERN <5 INDIVIDUALS
}

# We predict respiratory, endocrine/metabolic and mental health admissions disproportionately well
pdf(paste0(plot_dir,"Analytics/disease_class.pdf"),width=8,height=5)
par(mar=c(8.1,4.1,4.1,2.1))
xmax=8; mh=0
gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("red",gcol)

# Truncate xfreq for plotting simplicity
xfreq1=xfreq
xfreq1[,14]=rowSums(xfreq[,14:dim(xfreq)[2]])
xfreq1=xfreq1[,1:14]

xlabs=colnames(xfreq1)

m0=dim(xfreq1)[1]; n0=dim(xfreq1)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq1)),bty="n",xaxt="n",
  xlab="",ylab="% patients")
for (i in 1:dim(xfreq1)[2]) {
  xcol2=xcol; 
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq1[,i],col=xcol,lty=1,lwd=1)
}
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=xlabs,las=2,cex.axis=0.5)

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
  lwd=2)

dev.off()
sink(paste0(plot_dir,"Analytics/disease_class.txt"))
cat("Matrix of class frequencies: \n\n")
print(xfreq)
cat(sink_marker)
cat("xfreq="); dput(xfreq); cat("\n\n")
sink()


############################################################
## Distribution of predictions across admittees by type 1 ##
############################################################

xpred=log10(100*all_pred$super)
xvar=as.factor(icd2sick(all_pred$main_condition)) # this function is found in auxiliary.R
xl0=levels(xvar)

lscore0=xpred[which(all_pred$target==0)]
d0=density(lscore0)
xmed0=median(lscore0)

lscorex=xpred[which(all_pred$target==1)]
dx=density(lscorex)
xmedx=median(lscorex)


xmed=c(); dmax=c()
for (i in 1:length(xl0)) {
  w=which(xvar==xl0[i])
  lscore=xpred[w]
  dz=density(lscore)
  dmax=c(dmax,max(dx$y))
  xmed=c(xmed,median(lscore))
  assign(paste0("d",i),dz)
}

pdf(paste0(plot_dir,"Analytics/disease_class_violin.pdf"),width=6,height=5)
## If running from loaded data, first run
# lscore=c(0,0); lscore0=c(0,0)
par(mar=c(8.1,4.1,4.1,2.1))
dmax=0
for (i in 1:length(xmed)) {
  dz=get(paste0("d",i))
  dmax=c(dmax,max(dz$y))
}
sc=0.3/max(dmax)

plot(0,type="n",xlim=c(0,2+length(xl0)),ylim=c(0,2),bty="n",xaxt="n",yaxt="n",
  xlab="",ylab="Score (%), log scale")

polygon(c(sc*d0$y,rev(-sc*d0$y))+1/3,c(d0$x,rev(d0$x)),col="red",border=NA)
points(1/3,xmed0,pch=16)

polygon(c(sc*dx$y,rev(-sc*dx$y))+5/3,c(dx$x,rev(dx$x)),col="red",border=NA)
points(5/3,xmedx,pch=16)

ox=order(-xmed)
for (i0 in 1:length(xl0)) {
  i=ox[i0]
  dz=get(paste0("d",i))
  polygon(c(sc*dz$y,rev(-sc*dz$y))+i0+2,c(dz$x,rev(dz$x)),col="gray",border=NA)
  points(i0+2,xmed[i],pch=16)
}
axis(1,at=c(1/3,5/3,3:(2+length(xl0))),labels=c("Not admitted","All admitted",xl0[ox]),las=2,cex.axis=0.5)
yp=pretty(c(0,2),n=5)
axis(2,at=yp,labels=round((10^yp)),las=2,cex.axis=0.5)
dev.off()
sink(paste0(plot_dir,"Analytics/disease_class_violin.txt"))
cat("No_admission\n")
print(d0$x)
cat("\n")
print(d0$y)
cat("\n\n")
cat("All_admissions\n")
print(dx$x)
cat("\n")
print(dx$y)
cat("\n\n")
for (i in 1:length(xl0)) {
  cat(xl0[i],"\n")
  dx=get(paste0("d",i))
  print(dx$x)
  cat("\n")
  print(dx$y)
  cat("\n\n")
}
cat("Medians\n")
print(xmed)
cat(sink_marker)
cat("d0="); dput(d0); cat("\n\n")
cat("dx="); dput(dx); cat("\n\n")
cat("xl0="); dput(xl0); cat("\n\n")
cat("xmed0="); dput(xmed0); cat("\n\n")
cat("xmedx="); dput(xmedx); cat("\n\n")
for (i in 1:length(xl0)) {
  dx=get(paste0("d",i))
  cat(paste0("d",i,"=")); dput(dx); cat("\n\n")
}
cat("xmed="); dput(xmed); cat("\n\n")
sink()



######################################################
## Accuracy of predictions by diag. cl., PPV scale  ##
######################################################

xpred=all_pred$super
xvar=as.factor(icd2sick(all_pred$main_condition)) # this function is found in auxiliary.R
xl0=levels(xvar)

xmax=floor(10*max(xpred,na.rm=T))-1; mh=0

xfreq=table(xvar); 
nfreq=xfreq
xfreq=100*xfreq/length(xvar)

for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=table(xvar[sub])
  dx=100*dx/length(sub)
  mh=max(c(mh,dx))
  xfreq=rbind(xfreq,dx)
  nfreq=rbind(nfreq,table(xvar[sub]))
}

xfreq[which(nfreq<6)]=0

ox=order(-xfreq[1,])
xfreq=xfreq[,ox]; xlabs=xl0[ox]

gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("red",gcol)

# We predict respiratory, endocrine/metabolic and mental health admissions disproportionately well
pdf(paste0(plot_dir,"Analytics/disease_class_ppv.pdf"),width=8,height=5)
par(mar=c(8.1,4.1,4.1,2.1))

m0=dim(xfreq)[1]; n0=dim(xfreq)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq)),bty="n",xaxt="n",
  xlab="",ylab="PPV for condition (%)")
for (i in 1:dim(xfreq)[2]) {
  w=which(nfreq[,i]<5); xfreq[w,i]=5/sum(nfreq[,i]) # THIS IS IMPORTANT TO ENSURE NO PLOT ELEMENTS CONCERN <5 INDIVIDUALS
  xcol2=xcol; xcol2[w]="white"
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq[,i],col=xcol,lty=1,lwd=1)
}
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=xlabs,las=2,cex.axis=0.5)

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
  lwd=2)

dev.off()
sink(paste0(plot_dir,"Analytics/disease_class_ppv.txt"))
cat("Matrix of frequencies: \n\n")
print(xfreq)
cat(sink_marker)
cat("xfreq="); dput(xfreq); cat("\n\n")
sink()



######################################################
## Accuracy of predictions by admission type        ##
######################################################

xpred=all_pred$super
xvar=as.factor(adcode2lit(all_pred$admission_type)) # this function is found in auxiliary.R
xl0=levels(xvar)

xmax=floor(10*max(xpred,na.rm=T))-1; mh=0

xfreq=table(xvar); nfreq=xfreq
xfreq=100*xfreq/sum(xfreq)

for (i in 1:xmax) {
  sub=which(xpred> i/10)
  dx=table(xvar[sub])
  dx=100*dx/sum(dx)
  mh=max(c(mh,dx))
  xfreq=rbind(xfreq,dx)
  nfreq=rbind(nfreq,table(xvar[sub]))
}

xfreq[which(nfreq<6)]=0

ox=order(-xfreq[1,])
xfreq=xfreq[,ox]; xlabs=xl0[ox]

gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("red",gcol)

# We predict emergency admission better than urgent, and self inflicted injury better than external injury
pdf(paste0(plot_dir,"Analytics/admission_type.pdf"),width=6,height=5)
par(mar=c(6.8,4.1,4.1,2.1))

m0=dim(xfreq)[1]; n0=dim(xfreq)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq)),bty="n",xaxt="n",
  xlab="",ylab="% patients")
for (i in 1:dim(xfreq)[2]) {
  w=which(nfreq[,i]<5); xfreq[w,i]=5/sum(nfreq[,i]) # THIS IS IMPORTANT TO ENSURE NO PLOT ELEMENTS CONCERN <5 INDIVIDUALS
  xcol2=xcol; xcol2[w]="white"
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq[,i],col=xcol,lty=1,lwd=1)
}
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=xlabs,las=2,cex.axis=0.5)

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"red"),
  c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
  lwd=2)

dev.off()
sink(paste0(plot_dir,"Analytics/admission_type.txt"))
cat("Matrix of admission types: \n\n")
print(xfreq)
cat(sink_marker)
cat("xfreq="); dput(xfreq); cat("\n\n")
sink()


######################################################
## Predictive power in various patient cohorts      ##
######################################################

# Age
age_split=c(0,20,30,40,50,60,70,80,120)


perf=c(); xrate=c()
for (i in 1:(length(age_split)-1)) {
  sub=which(all_pred$age>age_split[i] & all_pred$age <= age_split[i+1] & is.finite(all_pred$super+all_pred$v3))
  psp=getroc(all_pred$target[sub],all_pred$super[sub])
  p3=getroc(all_pred$target[sub],all_pred$v3[sub])
  perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
  xrate=c(xrate,sum(all_pred$target[sub])/length(sub))
}

# 
pdf(paste0(plot_dir,"Analytics/performance_by_age.pdf"),width=6,height=5)
labs=c(); 
age_split=c(0,20,30,40,50,60,70,80,120)
for (i in 1:(length(age_split)-1)) labs=c(labs,paste0("(",age_split[i],",",age_split[i+1],"]"))

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=12; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
  xlab="",ylab="AUROC") #,main="ROC and EA frequency by age band")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}  
lines((m0+2)*(-0.5 + 1:n0),perfmin + xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=2)
axis(4,at=perfmin+seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomright",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()
sink(paste0(plot_dir,"Analytics/performance_by_age.txt"))
print(perf)
cat("\n\n")
print(xrate)
cat("\n\n")
cat(sink_marker)
cat("perf="); dput(perf); cat("\n\n")
cat("xrate="); dput(xrate); cat("\n\n")
sink()




# SIMD
simd=sort(unique(all_pred$simd))

perf=c(); xrate=c()
for (i in 1:(length(simd))) {
  sub=which(all_pred$simd == simd[i]  & is.finite(all_pred$super+all_pred$v3))
  psp=getroc(all_pred$target[sub],all_pred$super[sub])
  p3=getroc(all_pred$target[sub],all_pred$v3[sub])
  perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
  xrate=c(xrate,sum(all_pred$target[sub])/length(sub))
  labs=c(labs,simd[i])
}

# 
pdf(paste0(plot_dir,"Analytics/performance_by_simd.pdf"),width=6,height=5)
labs=1:10

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=25; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
  xlab="SIMD",ylab="AUROC") #,main="ROC and EA frequency by SIMD")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}
lines((m0+2)*(-0.5 + 1:n0),perfmin + xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=1)
axis(4,at=perfmin + seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomright",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()
sink(paste0(plot_dir,"Analytics/performance_by_simd.txt"))
print(perf)
cat("\n\n")
print(xrate)
cat(sink_marker)
cat("perf="); dput(perf); cat("\n\n")
cat("xrate="); dput(xrate); cat("\n\n")
sink()






# v3 cohorts
v3c=rep("",dim(all_pred)[1])
v3c[which(all_pred$age >= 75  & is.finite(all_pred$super+all_pred$v3))] = "FEC"
v3c[which(all_pred$age >= 16 & all_pred$age <= 74  & is.finite(all_pred$super+all_pred$v3))] = "LTC"
v3c[which(all_pred$age >= 16 & all_pred$age <= 55 & all_pred$ae2 >= 1  & is.finite(all_pred$super+all_pred$v3))] = "YED"
v3c[which(all_pred$age >= 16 & all_pred$age <= 55  & is.finite(all_pred$super+all_pred$v3))] = "LTC_YED"
v3c[which(all_pred$age < 16)] = "U16"

perf=c(); labs=c(); xrate=c()
sub=which(v3c=="FEC")
psp=getroc(all_pred$target[sub],all_pred$super[sub])
p3=getroc(all_pred$target[sub],all_pred$v3[sub])
perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
xrate=c(xrate,sum(all_pred$target[sub])/length(sub))

sub=which(v3c %in% c("LTC","LTC_YED"))
psp=getroc(all_pred$target[sub],all_pred$super[sub])
p3=getroc(all_pred$target[sub],all_pred$v3[sub])
perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
xrate=c(xrate,sum(all_pred$target[sub])/length(sub))

sub=which(v3c %in% c("YED","LTC_YED"))
psp=getroc(all_pred$target[sub],all_pred$super[sub])
p3=getroc(all_pred$target[sub],all_pred$v3[sub])
perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
xrate=c(xrate,sum(all_pred$target[sub])/length(sub))

sub=which(v3c %in% c("U16"))
psp=getroc(all_pred$target[sub],all_pred$super[sub])
p3=getroc(all_pred$target[sub],all_pred$v3[sub])
perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
xrate=c(xrate,sum(all_pred$target[sub])/length(sub))



## All of YED are in LTC_YED so this is unnecessary
#sub=which(v3c == "LTC_YED")
#psp=getroc(all_pred$target[sub],all_pred$super[sub])
#p3=getroc(all_pred$target[sub],all_pred$v3[sub])
#perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
#xrate=c(xrate,sum(all_pred$target[sub])/length(sub))

# Row and column names
colnames(perf)=c("FEC","LTC","YED","U16")
rownames(perf)=c("v4_auc","v3_auc","v4_se","v3_se")
names(xrate)=colnames(perf)


# 
pdf(paste0(plot_dir,"Analytics/performance_by_v3_cohort.pdf"),width=6,height=5)
labs=c("FEC","LTC","YED","U16") #,"LTC,YED")

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=20; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
  xlab="Cohort",ylab="AUROC") #,main="ROC and EA frequency in v3 cohorts")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
    (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}
lines((m0+2)*(-0.5 + 1:n0),perfmin+xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=1)
axis(4,at=perfmin + seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomleft",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()
sink(paste0(plot_dir,"Analytics/performance_by_v3_cohort.txt"))
print(perf)
cat("\n\n")
print(xrate)
cat(sink_marker)
cat("perf="); dput(perf); cat("\n\n")
cat("xrate="); dput(xrate); cat("\n\n")
sink()







######################################################################################
######################################################################################
## Shapley values. The following scripts analyse individual contributions to risk   ##
##  in SPARRAv4 and overall variable importance using Shapley values.               ##
######################################################################################
######################################################################################

######################################################################################
## Read in data. These are generated in estimate_shapley_values.R                   ##
######################################################################################

# Fold 1
shapleyx=readRDS(file=paste0("James/Analysis/Data/full_model/shapley_values.RDS"))
xpred=shapleyx$id_matrix
shapley=shapleyx$values

# Further details not in training matrix
if (!exists("all_pred")) all_pred=readRDS(all_pred_file)
dpred=all_pred[shapleyx$ids_in_full,]

# Normalise to non-topic values
vx=which(!grepl("topic",colnames(shapley)))
shapley=shapley[,vx]

if ("MASS" %in% (.packages())) detach("package:MASS",unload=TRUE)




######################################################################################
## Variable importance                                                              ##
######################################################################################

vx=colMeans(abs(shapley %>% select(-"v3score")),na.rm=T)
names(vx)=colnames(shapley %>% select(-"v3score"))
vx=vx[order(-vx)] # descending order

pdf(paste0(plot_dir,"Shapley_values/variable_importance.pdf"),width=6,height=5)
par(mar=c(5.1,12.1,4.1,2.1))

# For table:
# cbind(longvarnames(names(vx)[1:30]), signif(as.numeric(100*vx[1:30]),digits=2))

nvx=10 # plot this many top values
plot(0,xlim=c(0,max(vx)),ylim=c(0,nvx+1),type="n",bty="n",
  xlab="Mean absolute Shap. Val.",ylab="",yaxt="n",main="Variable importance")
axis(2,at=1:nvx,labels=longvarnames(names(vx)[nvx:1]),las=2,cex.axis=0.5)
segments(0,nvx:1,vx[1:nvx],nvx:1,lwd=2)
dev.off()
sink(paste0(plot_dir,"Shapley_values/variable_importance.txt"))
cat("Variables by importance\n\n")
print(vx)
cat(sink_marker)
cat("vx="); dput(vx); cat("\n\n")
sink()




######################################################################################
## Age vs Shapley value for age                                                     ##
######################################################################################

library(MASS)

rcex=1.2 # range expansion 
xx=xpred$age
yy=shapley$age
sexM=dpred$sexM
xxm=xx[which(sexM==TRUE)]; yym=yy[which(sexM==TRUE)]
xxf=xx[which(sexM==FALSE)]; yyf=yy[which(sexM==FALSE)]

res=50; kern=0.3
w=which(is.finite(xx+yy))
xx=xx[w]; yy=yy[w]
wm=which(is.finite(xxm+yym)); wf=which(is.finite(xxf + yyf))
xxm=xxm[wm]; yym=yym[wm]; xxf=xxf[wf]; yyf=yyf[wf]

# X co-ordinates to evaluate at
xgrid=seq(min(xx,na.rm=T),max(xx,na.rm=T),length.out=res)

# elegant but impractical
# tmat=outer(xgrid,xx,function(x,y) dnorm(x-y,sd=kern))
xytrue=mkernm(xgrid,xx,yy,kern=kern)
xtrue=xytrue[,1]; ytrue=xytrue[,2]; yvar=xytrue[,3]

# Split by sex. No difference
xytruem=mkernm(xgrid,xxm,yym,kern=kern)
xtruem=xytruem[,1]; ytruem=xytruem[,2]; yvarm=xytruem[,3]
xytruef=mkernm(xgrid,xxf,yyf,kern=kern)
xtruef=xytruef[,1]; ytruef=xytruef[,2]; yvarf=xytruef[,3]

# trim for privacy
w=which(xtrue<85)
xtrue=xtrue[w]; ytrue=ytrue[w]; yvar=yvar[w]

# Background image: can't plot actual points for privacy reasons, density might be OK
xrange=c(0,85) # range(xx); 
yrange=mean(ytrue) + rcex*(range(ytrue)-mean(ytrue)); 
kres=100; kx=(max(xrange)-min(xrange))/20; ky=(max(yrange)-min(yrange))/20 # resolution, x kernel width, y kernel width
scx=kde2d(xx,yy,c(kx,ky),n=kres,lims=c(xrange,yrange))
scx$z=scx$z/outer(rowSums(scx$z),rep(1,kres))
icol=gray(1-(1:100)/200)


pdf(paste0(plot_dir,"Shapley_values/age_effect.pdf"),width=6,heigh=5)
xrange=c(0,85)
icol=gray(1-(1:100)/200)

plot(0,type="n",xlim=xrange,ylim=yrange,bty="n",
  xlab="Age",ylab="Shapley value: Age") #,main="Importance of age as predictor")
image(scx,add=T,col=icol) 
#plot(xx,yy,col="gray",cex=0.5) # leave this out for the moment for privacy reasons

lines(xtrue,ytrue)
lines(xtrue,ytrue + sqrt(yvar),lty=2)
lines(xtrue,ytrue - sqrt(yvar),lty=2)

# split by sex: uninteresting
#lines(xtruem,ytruem,col="red")
#lines(xtruef,ytruef,col="blue")

legend("topleft",c("Mean","SD","Density"),lty=c(1,2,NA),
  pch=c(NA,NA,16),pt.cex=c(NA,NA,2),col=c("black","black","gray"),bty="n")

dev.off()
sink(paste0(plot_dir,"Shapley_values/age_effect.txt"))
cat(xtrue)
cat("\n")
cat(ytrue)
cat("\n")
cat(yvar)
cat("\n\n")
cat("Male_only\n")
cat(xtruem)
cat("\n")
cat(ytruem)
cat("\n")
cat("Female_only\n")
cat(xtruef)
cat("\n")
cat(ytruef)
cat("\n\n")
cat("Background density\n")
print(scx)
cat(sink_marker)
cat("xtrue="); dput(xtrue); cat("\n\n")
cat("ytrue="); dput(ytrue); cat("\n\n")
cat("yvar="); dput(yvar); cat("\n\n")
cat("xtruem="); dput(xtruem); cat("\n\n")
cat("ytruem="); dput(ytruem); cat("\n\n")
cat("xtruef="); dput(xtruef); cat("\n\n")
cat("ytruef="); dput(ytruef); cat("\n\n")
cat("scx="); dput(scx); cat("\n\n")
sink()




######################################################################################
## SIMD vs Shapley value for SIMD                                                   ##
######################################################################################

rcex=3 # range expansion 
dbw=0.0003 # bandwidth for KDE: must be wide enough for anonymisation
xx=xpred$SIMD_DECILE_2016_SCT
yy=shapley$SIMD_DECILE_2016_SCT
w=which(is.finite(xx+yy))
xx=xx[w]; yy=yy[w]
res=100; kern=0.5

dxx=matrix(0,10,512); dxy=dxx
mx=rep(0,10); sdx=mx
alpha=0.005
for (i in 1:10) {
  w=which(xx==i)
  yx=yy[w]
  mx[i]=mean(yx,na.rm=T)
  sdx[i]=sd(yx,na.rm=T)
  dx=density(yx,bw=dbw,na.rm=T,from=max(min(yx),mx[i]-3*sdx[i]),to=min(max(yx),mx[i]+3*sdx[i]))
  # Privacy - flatten small bits of density
  w=which.max(dx$y); dy=dx$y
  dy[1:w]=cummax(dy[1:w])
  dy[(w+1):length(dy)]=rev(cummax(rev(dy[(w+1):length(dy)])))
  dxx[i,]=dx$x; dxy[i,]=dy
}
# scale densities to max height 2/3
dxy=dxy/(1.5*max(dxy))

pdf(paste0(plot_dir,"Shapley_values/simd_effect.pdf"),width=6,heigh=5)

yrange=mean(mx) + rcex*(range(mx)-mean(mx))
plot(0,type="n",xlim=c(1,11),ylim=yrange,bty="n",
  xlab="SIMD",ylab="Shapley value: SIMD") #,main="Importance of SIMD as predictor")

for (i in 1:10) polygon(i+dxy[i,],dxx[i,],col="gray",border=NA)
points(1:10,mx,pch=16)

legend("topright",c("Mean","Density"),pch=c(16,16),col=c("black","gray"),bty="n")

dev.off()
sink(paste0(plot_dir,"Shapley_values/simd_effect.txt"))
print(mx)
cat("\n")
print(dxy)
cat(sink_marker)
cat("mx="); dput(mx); cat("\n\n")
cat("dxy="); dput(dxy); cat("\n\n")
sink()


######################################################################################
## Age SIMD equivalents                                                             ##
######################################################################################

dsimd=mx[1]-mx[10]
ytd=rev(cummin(rev(ytrue))) 

aget=seq(18,80) # if you are this old...
age_d=aget

for (i in 1:length(aget)) {
  s0=suppressWarnings(approx(xtrue,ytd,aget[i])$y) # average Shapley value for age for people of this age
  age_d[i]=suppressWarnings(approx(ytd,xtrue,s0+dsimd)$y) # equivalent age for Shapley value + SIMD contribution
}

pdf(paste0(plot_dir,"Shapley_values/equivalent_years_older_shapley.pdf"),width=5,height=5)

plot(0,type="n",xlim=range(aget),ylim=range(age_d-aget),
  xlab="For individuals of this age",ylab="SIMD1 vs 10 equiv. to this many years older")
lines(aget,age_d-aget)
dev.off()
sink(paste0(plot_dir,"Shapley_values/equivalent_years_older_shapley.txt"))
cat(aget)
cat("\n\n")
cat(age_d)
cat(sink_marker)
cat("aget="); dput(aget); cat("\n\n")
cat("age_d="); dput(age_d); cat("\n\n")
sink()



######################################################################################
## N adm vs Shapley value for N adm                                                 ##
######################################################################################

rcex=3 # range expansion 
dbw=0.001 # bandwidth for KDE: must be wide enough for anonymisation
xx1=xpred$num_emergency_admissions
xx2=xpred$num_elective_admissions
yy1=shapley$num_emergency_admissions
yy2=shapley$num_elective_admissions
xx1=case_when(
  xx1 == 0 ~ 1,
  xx1 %in% 1:3 ~ 2,
  xx1 %in% 4:6 ~ 3,
  xx1 %in% 7:9 ~ 4,
  xx1 > 9 ~ 5
)
xx2=case_when(
  xx2 == 0 ~ 1,
  xx2 %in% 1:3 ~ 2,
  xx2 %in% 4:6 ~ 3,
  xx2 %in% 7:9 ~ 4,
  xx2 > 9 ~ 5
)
w=which(is.finite(xx1+xx2+yy1+yy2))
xx1=xx1[w]; xx2=xx2[w]; yy1=yy1[w]; yy2=yy2[w]
res=100; kern=0.5



dxx1=matrix(0,max(xx1),512); dxy1=dxx1
dxx2=matrix(0,max(xx2),512); dxy2=dxx2

mx1=rep(0,max(xx1)); sdx1=mx1; mx2=mx1; sdx2=sdx1
for (i in  1:max(xx1)) {
  w=which(xx1==i)
  yx=yy1[w]
  mx1[i]=mean(yx,na.rm=T)
  sdx1[i]=sd(yx,na.rm=T)
  dx=density(yx,na.rm=T,bw=(max(yx)-min(yx))/5)
  dxx1[i,]=dx$x; dxy1[i,]=dx$y
  
  w=which(xx2==i)
  yx=yy2[w]
  mx2[i]=mean(yx,na.rm=T)
  sdx2[i]=sd(yx,na.rm=T)
  dx=density(yx,na.rm=T,bw=(max(yx)-min(yx))/5)
  dxx2[i,]=dx$x; dxy2[i,]=dx$y
}
# scale densities to max height 2/3
dxy1=dxy1/(2.1*max(dxy1))
dxy2=dxy2/(2.1*max(dxy2))

pdf(paste0(plot_dir,"Shapley_values/n_adm_effect.pdf"),width=6,heigh=5)

yrange=mean(c(mx1,mx2)) + rcex*(range(c(mx1,mx2)-mean(c(mx1,mx2))))
plot(0,type="n",xlim=c(min(xx1)-1,1+max(xx1)),ylim=yrange,bty="n",xaxt="n",
  xlab="Num. admissions",ylab="Shapley value: num. admissions") #,main="Importance of num. admissions as predictor")
axis(1,at=1:max(xx1),label=c("0","1:3","4:6","7:9","10+"))
abline(h=0,lty=2)

for (i in 1:max(xx1)) polygon(i+0.01+dxy1[i,],dxx1[i,],col="lightgray",border=NA)
for (i in 1:max(xx1)) polygon(i-0.01-dxy2[i,],dxx2[i,],col="darkgray",border=NA)

points(1:max(xx1),mx1,pch=16)
points(1:max(xx1),mx2,pch=16,col="red")


legend("topleft",c("Em. mean","Em. density","Elec. mean","Elec. density"),
  pch=c(16,16),col=c("black","lightgray","red","darkgray"),bty="n")

dev.off()
sink(paste0(plot_dir,"Shapley_values/n_adm_effect.txt"))
print(mx)
cat("\n")
print(dxx1)
cat("\n")
print(dxy1)
cat("\n")
print(dxx2)
cat("\n")
print(dxy2)
cat(sink_marker)
cat("mx="); dput(mx); cat("\n\n")
cat("dxx1="); dput(dxx1); cat("\n\n")
cat("dxy1="); dput(dxy1); cat("\n\n")
cat("dxx2="); dput(dxx1); cat("\n\n")
cat("dxy2="); dput(dxy2); cat("\n\n")
sink()


######################################################################################
## SIMD x Age vs Shapley value for age                                              ##
######################################################################################

library(MASS)

rcex=1.8 # range expansion 
xx=xpred$age
yy=shapley$age
zz=xpred$SIMD_DECILE_2016_SCT
res=50; kern=0.5
w=which(is.finite(xx+yy+zz))
xx=xx[w]; yy=yy[w]; zz=zz[w]

xgrid=seq(min(xx,na.rm=T),max(xx,na.rm=T),length.out=res)

ytrue=matrix(0,res,10); xtrue=matrix(0,res,10); yvar=ytrue
for (j in 1:10) {
 sub=which(zz==j)  
 for (i in 1:res) {
  wt=dnorm(xgrid[i]-xx[sub],sd=kern)
  s1=sum(wt); s2=sum(wt^2)
  ytrue[i,j]=sum(wt*yy[sub]) / s1
  xtrue[i,j]=sum(wt*xx[sub]) / s1
  yvar[i,j]= (s1/(s1^2 - s2))*sum(wt*((yy[sub]-ytrue[i])^2))
  if (sum(wt/s1)^2 < 1/5) { # THIS IS IMPORTANT: no density estimate is equivalent to a mean of <5 elements
    xtrue[i,j]=NA
    ytrue[i,j]=NA
  }
 }
}

# trim for privacy
library(matrixStats)
w=which(rowMaxs(xtrue)<85)
xtrue=xtrue[w,]; ytrue=ytrue[w,]; yvar=yvar[w,]


pdf(paste0(plot_dir,"Shapley_values/age_by_simd_effect.pdf"),width=6,heigh=5)

plot(0,type="n",xlim=c(0,90),ylim=range(ytrue),bty="n",
  xlab="Age",ylab="Shapley value: Age",main="Importance of age as predictor")
# image(scx,add=T,col=icol) # leave this out for the moment for privacy reasons
# plot(xx,yy,col="gray",cex=0.5)

for (i in 1:10) {
lines(xtrue[,i],ytrue[,i],col=rep(c("darkred","red","orange","yellow"),each=3)[i],lty=rep(1:3,4)[i])
}

legend(20,max(ytrue),title="SIMD",legend=1:10,col=(rep(c("darkred","red","orange","yellow"),each=3))[1:10],lty=rep(1:3,4)[1:10],bty="n")

dev.off()
sink(paste0(plot_dir,"Shapley_values/age_by_simd_effect.txt"))
print(xtrue)
cat("\n")
print(ytrue)
cat("\n")
print(yvar)
cat(sink_marker)
cat("xtrue="); dput(xtrue); cat("\n\n")
cat("ytrue="); dput(ytrue); cat("\n\n")
cat("yvar="); dput(yvar); cat("\n\n")
sink()



######################################################################################
## Direct effect of age and SIMD on admission rate                                  ##
######################################################################################

xx=seq(2,80); nn=xx; yy=xx; 
ysimd=matrix(0,length(xx),10); nsimd=ysimd
for (i in 1:length(xx)) {
  w=which(all_pred$age==xx[i])
  yy[i]=mean(all_pred$target[w])
  nn[i]=length(w)
  for (s in 1:10) { 
    ws=which(all_pred$age==xx[i] & all_pred$simd==s)
    ysimd[i,s]=mean(all_pred$target[ws])
    nsimd[i]=length(ws)
  }
}

pdf(paste0(plot_dir,"Shapley_values/age_simd_direct_effect.pdf"),width=5,height=5)

plot(0,type="n",xlim=range(xx),xlab="Age",ylab="EA freq.",col="red",ylim=range(ysimd))
for (i in 1:10) lines(xx,ysimd[,i],col=gray(i/20))
lines(xx,yy,col="red")
#lines(xx,yyf,col="red")
legend("topleft",c("All","1","...","10"),lty=c(1,1,NA,1),col=c("red",gray(1/20),NA,gray(1/2)))

dev.off()
sink(paste0(plot_dir,"Shapley_values/age_simd_direct_effect.txt"))
cat("Ages:\n")
cat(xx)
cat("\n\nMean\n")
cat(yy)
cat("\n\nBy_simd\n")
cat(ysimd)
cat(sink_marker)
cat("xx="); dput(xx); cat("\n\n")
cat("yy="); dput(yy); cat("\n\n")
cat("ysimd="); dput(ysimd); cat("\n\n")
sink()



######################################################################################
## Equivalent ages for different SIMD deciles                                       ##
######################################################################################

yeq=0*ysimd
ydsimd=ysimd
yd=rev(cummin(rev(yy)))
for (s in 1:10) ydsimd[,s]=rev(cummin(rev(ysimd[,s])))
for (i in 1:length(xx)) {
  for (s in 1:10) {
    yeq[i,s]=approx(yd,xx,ydsimd[i,s])$y # xx[min(which(ydsimd[,s]>yd[i]))]
  }
}
yeq[which(xx<18),]=NA #
yeq[which(yeq<18)]=NA # risk essentially as low as it ever goes


xcol=colorRampPalette(c("blue","gray","red"))(10)
pdf(paste0(plot_dir,"Shapley_values/age_simd_equivalent.pdf"),width=5,height=5)
plot(0,type="n",xlab="Chron. age",ylab="Effective age",xlim=c(15,max(xx)),ylim=c(20,max(xx)))
for (i in 1:10) lines(xx,yeq[,i],col=xcol[i])
abline(0,1,col="black",lty=2)
legend("bottomright",c("1","...","10"),title="SIMD",lty=1,col=c(xcol[1],xcol[5],xcol[10]))
dev.off()
sink(paste0(plot_dir,"Shapley_values/age_simd_equivalent.txt"))
print(yeq)
cat(sink_marker)
cat("yeq="); dput(yeq); cat("\n\n")
sink()


######################################################################################
## All Shapley value plots                                                          ##
######################################################################################

for (ii in 1:dim(shapley)[2]) {# ,which(colnames(shapley) %in% c("target","v3score")))) {

xx=xpred[[which(colnames(xpred)==colnames(shapley)[ii])]]; yy=shapley[[ii]]

if (length(unique(xx))>2) {
w=which(!is.na(xx+yy)); xx=xx[w]; yy=yy[w]
xr=sort(xx)[c(6,length(xx)-6)] # range
w=which(xx>=xr[1] & xx<= xr[2]) # privacy
xx=xx[w]; yy=yy[w]

nq=10
qx=quantile(unique(xx),(0:nq)/nq)
if (length(unique(qx))<3) nq=0
if (nq>3) while(nq>3 & length(which(table(cut(xx,breaks=unique(qx),include.lowest=TRUE))>5))<3) {
  nq=nq-1
  qx=quantile(unique(xx),(0:nq)/nq)
}


if (nq>3) {
qx[nq]=0.01+qx[nq]
mx=rep(0,nq)
sdx=mx
x0=mx
lab=rep("",nq)
for (i in 1:nq) {
  sub=which(xx>=qx[i] & xx<qx[i+1])
  if (length(sub)>5) {
  mx[i]=mean(yy[sub])
  x0[i]=mean(xx[sub])
  sdx[i]=sd(yy[sub])
  xx1=suppressWarnings(paste0(signif(min(xx[sub]),digits=2)))
  xx2=suppressWarnings(paste0(signif(max(xx[sub]),digits=2)))
  if (xx1==xx2) lab[i]=xx1 else lab[i]=paste(xx1,"-",xx2)
  } else {
    mx[i]=NA; x0[i]=NA; sdx[i]=NA; lab[i]=NA
  }
}
} else {
  mx=NA; x0=NA; sdx=NA; lab=NA
}

w=which(is.finite(x0))
x0=x0[w]; mx=mx[w]; sdx=sdx[w]; lab=lab[w]

} else {

x0=range(xx); 
mx=c(mean(yy[which(xx==min(xx))]),mean(yy[which(xx==max(xx))]))
sdx=c(sd(yy[which(xx==min(xx))]),sd(yy[which(xx==max(xx))]))
lab=c("No","Yes")
}

assign(paste0("x0_",ii),x0)
assign(paste0("mx_",ii),mx)
assign(paste0("sdx_",ii),sdx)
assign(paste0("lab_",ii),lab)

}

xmax=c(); xmin=c()
for (ii in 1:dim(shapley)[2]) {
  x0=get(paste0("x0_",ii))
  mx=get(paste0("mx_",ii))
  sdx=get(paste0("sdx_",ii))
  if (length(mx)>0) {
    xmax=c(xmax,max(mx+sdx,na.rm=T))
    xmin=c(xmin,min(mx-sdx,na.rm=T))
  } else {
    xmax=c(xmax,0)
    xmin=c(xmin,0)
  }
}
names(xmax)=colnames(shapley)
names(xmin)=colnames(shapley)
xmax=xmax[which(!names(xmax) %in% c("v3score","target"))]
xmin=xmin[which(!names(xmin) %in% c("v3score","target"))]

ecut=0.07
large_effect=names(xmax)[which(xmax-xmin>ecut)]; 
small_effect=names(xmin)[which(xmax-xmin <= ecut)]

large_range=range(c(xmin[large_effect],xmax[large_effect]))
small_range=range(c(xmin[small_effect],xmax[small_effect]))



for (ii in setdiff(1:dim(shapley)[2],which(colnames(shapley) %in% c("v3score","target")))) {

x0=get(paste0("x0_",ii))
mx=get(paste0("mx_",ii))
sdx=get(paste0("sdx_",ii))
lab=get(paste0("lab_",ii))


if (!all(is.na(x0))) {

if (colnames(shapley)[ii] %in% large_effect) {
  yvar=large_range 
  suffix="large_range"
} else {
  yvar=small_range  
  suffix="small_range"
}

pdf(paste0(plot_dir,"Shapley_values/All/shapley_",colnames(shapley)[ii],"_",suffix,".pdf"),width=7,height=6)
par(mar=c(7,4.1,4.1,2.1))
plot(0,type="n",xlim=range(x0),ylim=yvar,main=longvarnames(colnames(shapley)[ii]),xaxt="n",
  xlab="",ylab="Shapley value") # note fixed ylim
lines(x0,mx,lty=2,col="red")
segments(x0,mx-sdx,x0,mx+sdx,col="red")
points(x0,mx,pch=16,col="red")
axis(1,at=x0,labels=lab,las=2)
abline(h=0,lty=2)
dev.off()

sink(paste0(plot_dir,"Shapley_values/All/shapley_",colnames(shapley)[ii],"_",suffix,".txt"))
cat("X\n")
cat(x0)
cat("\n\n")
cat("Mean\n")
cat(mx)
cat("\n\nSD\n")
cat(sdx)
cat(sink_marker)
cat("x0="); dput(x0); cat("\n\n")
cat("mx="); dput(mx); cat("\n\n")
cat("sdx="); dput(sdx); cat("\n\n")
sink()
}

}


######################################################################################
## Topic analysis                                                                   ##
######################################################################################

# Need to look at fold 1 only
T1=readRDS(paste0(model_dir,"topic_model_fit_fold_23.rds"))
T2=readRDS(paste0(model_dir,"topic_model_fit_fold_13.rds"))
T3=readRDS(paste0(model_dir,"topic_model_fit_fold_12.rds"))

# Lookup tables for ICD10 and BNF
icd_lookup=read.csv("James/SPARRAv4/Lookup/icd10.csv",stringsAs=F)
bnf_lookup=read.csv("James/SPARRAv4/Lookup/BNF_lookup.csv",stringsAs=F)

# Report whenever probability of topic including word is > threshold
pthresh=0.01

# Topics with highest mean absolute Shapley values in fold 1
shapleyx=readRDS(file=paste0("James/Analysis/Data/full_model/shapley_values.RDS"))
xpred=shapleyx$id_matrix
shapley=shapleyx$values
dpred=all_pred[shapleyx$ids_in_full,]
vxt=colMeans(abs(shapley),na.rm=T)
names(vxt)=colnames(shapley)
vxt=vxt[order(-vxt)] # descending order, including topics this time

# Names and indices of topics in descending order of importance
txo=names(vxt)[grep("topic_",names(vxt))]
ixo=as.numeric(gsub("_","",substring(txo,7,8)))





# Privacy check - IMPORTANT

nfile="James/Analysis/Data/full_model/topic_check.RData"

if (!file.exists(nfile)) {

T1M=readRDS("James/Analysis/Data/full_model/doc_term_sparsematrix.rds")
T1M1=(T1M>0)

for (f in 1:3) {
  xf=c("23","13","12")[f]
  T1=readRDS(paste0("James/Analysis/Data/full_model/topic_model_fit_fold_",xf,".rds"))
  
  B1=T1@beta
  
  n1=rep(0,30)
  for (i in 1:30) {
    w=which(exp(B1[i,])>0.01)
    if (length(w)>1) {
      T1A=rowSums(T1M1[,w])
    } else T1A=T1M1[,w]
    T1B=rowSums(T1M1[,-w])
    n1[i]=length(which(T1A==length(w) & T1B==0))
    print(i)
  }
  assign(paste0("N",f),n1)
}

save(N1,N2,N3,file=nfile)

} else load(nfile)



# Look at betas; beta[i,j]=log prob that word j is in topic i
for (f in 1:3) {

# Get topic model
TM=get(paste0("T",f))

# Get numbers for fold f
NN=get(paste0("N",f))

# Order for topics
if (f==1) ix=ixo else ix=1:30

# Beta matrix
beta=TM@beta


topic_summary=list()
for (i in 1:30) { # ten most important topics
  if (NN[ix[i]]==0 | NN[ix[i]]>5) {
    w=which(exp(beta)[ix[i],]>pthresh)
    beta_w=exp(beta)[ix[i],w]
    xterms=TM@terms[w[order(-beta_w)]]
    terms=c()
    for (j in 1:length(xterms)) {
     if (xterms[j] %in% icd_lookup$Code.1) {
       term_j=paste0("(ICD10) ",
         icd_lookup$Full.Description[which(icd_lookup$Code.1==xterms[j])])
     } else if (xterms[j] %in% bnf_lookup$BNF_Section2) {
       term_j=paste0("(BNF) ",
         bnf_lookup$Description[which(bnf_lookup$BNF_Section2==xterms[j])])
     } else term_j=xterms[j]
     terms=c(terms,term_j)
    }
  topic_summary[[i]]=list(terms=terms,prob=sort(beta_w,dec=T),N=NN[ix[i]])
  }  else topic_summary[[i]]=list(terms=NULL,prob=NULL,N=NN[ix[i]])
}

assign(paste0("topic_summary",f),topic_summary)

}

sink(paste0(plot_dir,"Shapley_values/topic_breakdown.txt"))
cat("topic_summary1="); dput(topic_summary1); cat("\n\n")
cat("topic_summary2="); dput(topic_summary2); cat("\n\n")
cat("topic_summary3="); dput(topic_summary3); cat("\n\n")
sink()




######################################################################################
######################################################################################
## Miscellaneous                                                                    ##
######################################################################################
######################################################################################

if ("MASS" %in% (.packages())) detach("package:MASS",unload=TRUE)
load_cleanData(partition = "all",
  subset_frac=1, subset_seed = 1234,
  load_to_global=TRUE)

t0=as.Date("2may2016", "%d%b%Y")

inames=c("AE2","PIS","SMR00","SMR01","SMR04","SMR01E_SystemWatch")
longnames=c("A_E","Prescr","Outpt","Inp_day","MH_inp_day","Ger_systwatch")
for (i in 1:length(inames)) {
  if (i %in% 1:5) tabid=list_of_data_tables[[inames[i]]]$id else tabid=c(list_of_data_tables$SMR01E$id,list_of_data_tables$SystemWatch$id)
  pref=match(intersect(tabid,patients$id),patients$id)
  w=which(patients$date_of_death[pref]<t0)
  pref=pref[-w]
  
  pdf(paste0(plot_dir,"Description/age_",inames[i],"_",longnames[i],".pdf"),width=3,height=3)
  par(mar=c(2,2,0.1,0.1))
  xage=round(as.numeric(t0-patients$date_of_birth[pref])/365)
  xsex=(patients$gender[pref]=="Male")
  dm=density(xage[which(xsex)],na.rm=T,from=0,to=90,bw=3)
  df=density(xage[which(!xsex)],na.rm=T,from=0,to=90,bw=3)
  plot(0,type="n",xlim=c(0,90),ylim=c(0,max(c(dm$y,df$y,0.015))),
    xaxt="n",yaxt="n",xlab="",main="",ylab="",bty="n")
  c1=rgb(0,0,1,0.2); c2=rgb(1,0,0,0.2)
  polygon(c(0,dm$x,90,0),c(0,dm$y,0,0),col=c1,border=NA)
  polygon(c(0,df$x,90,0),c(0,df$y,0,0),col=c2,border=NA)
  axis(1,at=c(0,45,90),labels=F)
  axis(2,at=c(0,0.015),labels=F)
  dev.off()
  sink(paste0(plot_dir,"Description/age_",inames[i],"_",longnames[i],".txt"))
  cat("Male\n")
  cat(dm$x,"\n")
  cat(dm$y,"\n\n")
  cat("Female\n")
  cat(df$x,"\n")
  cat(df$y)
  cat(sink_marker)
  cat("dm="); dput(dm); cat("\n\n")
  cat("df="); dput(df); cat("\n\n")
  sink()
  
  pdf(paste0(plot_dir,"Description/simd_",inames[i],"_",longnames[i],".pdf"),width=3,height=3)
  par(mar=c(2,2,0.1,0.1))
  xsimd0=patients$SIMD_DECILE_2016_SCT[pref]
  xsimd0=xsimd0[which(!is.na(xsimd0))]
  xsimd=table(xsimd0)
  xsimd=xsimd/sum(xsimd,na.rm=T)
  plot(0,type="n",xlim=c(1,10),ylim=c(0,max(c(xsimd,0.15))),
    xaxt="n",yaxt="n",xlab="",main="",ylab="",bty="n")
  segments(1:10,0,1:10,xsimd,lwd=2)
  axis(1,at=c(1:10),labels=F)
  axis(2,at=c(0,0.1),labels=F)
  dev.off()
  sink(paste0(plot_dir,"Description/simd_",inames[i],"_",longnames[i],".txt"))
  cat(xsimd)
  cat(sink_marker)
  cat("xsimd="); dput(xsimd); cat("\n\n")
  sink()
  
}


xid=all_pred$id
xage=all_pred$age
xsex=all_pred$sex
xsimd0=all_pred$simd

pdf(paste0(plot_dir,"Description/age_all.pdf"),width=3.5,height=3.5)
par(mar=c(4,4,0.1,0.1))
dm=density(xage[which(xsex)],na.rm=T,from=0,to=90,bw=3)
df=density(xage[which(!xsex)],na.rm=T,from=0,to=90,bw=3)
plot(0,type="n",xlim=c(0,90),ylim=c(0,max(c(dm$y,df$y))),
  xaxt="n",yaxt="n",xlab="Age",main="",ylab="Density",bty="n")
c1=rgb(0,0,1,0.2); c2=rgb(1,0,0,0.2)
polygon(c(0,dm$x,90,0),c(0,dm$y,0,0),col=c1,border=NA)
polygon(c(0,df$x,90,0),c(0,df$y,0,0),col=c2,border=NA)
axis(1,at=c(0,45,90))
axis(2,at=c(0,0.015))
legend("topright",c("M","F"),pch=16,cex=1.5,col=c(c1,c2),bty="n")
dev.off()
sink(paste0(plot_dir,"Description/age_all.txt"))
cat("Male\n")
cat(dm$x,"\n")
cat(dm$y,"\n\n")
cat("Female\n")
cat(df$x,"\n")
cat(df$y)
cat(sink_marker)
cat("dm="); dput(dm); cat("\n\n")
cat("df="); dput(df); cat("\n\n")
sink()

pdf(paste0(plot_dir,"Description/simd_all.pdf"),width=3.5,height=3.5)
par(mar=c(4,4,0.1,0.1))
xsimd0=xsimd0[which(!is.na(xsimd0))]
xsimd=table(xsimd0)
xsimd=xsimd/sum(xsimd)
plot(0,type="n",xlim=c(1,10),ylim=c(0,max(xsimd)),
  xaxt="n",yaxt="n",xlab="SIMD decile",main="",ylab="Frequency",bty="n")
segments(1:10,0,1:10,xsimd,lwd=2)
axis(1,at=1:10)
axis(2,at=c(0,0.1))
dev.off()
sink(paste0(plot_dir,"Description/simd_all.txt"))
print(xsimd)
cat(sink_marker)
cat("xsimd="); dput(xsimd); cat("\n\n")
sink()






# Load training matrix for metadata
fullmatrix=paste0("James/Analysis/Data/full_model/full_data_matrix.RDS")
nhs=readRDS(fullmatrix)
exclusions_file=paste0("James/Analysis/Data/full_model/design_exclusions.RData")
load(exclusions_file)

# Generate table of demographics
xtab=c()
tx=sort(unique(nhs$time)); 
amax=9 # count individuals in 10y brackets up to age 10*amax, then >amax
varnames=setdiff(colnames(nhs),c(colnames(nhs)[grep("topic",colnames(nhs))],"target","id","age","SIMD_DECILE_2016_SCT","time","cv","v3score"))
for (i in 1:4) {
  if (i==1) sub=1:dim(nhs)[1] # all individual/time pairs under consideration
  if (i==2) sub=v34 # all individual/time pairs included
  if (i==3) sub=intersect(v34, which(nhs$target==TRUE)) # admissions
  if (i==4) sub=intersect(v34, which(nhs$target==FALSE)) # non-admissions
  
  stat=length(sub) # total number
  stat=c(stat,length(which(patients$gender[match(nhs$id[sub],patients$id)]=="Male"))) # number male
  stat=c(stat,length(which(patients$gender[match(nhs$id[sub],patients$id)]=="Female"))) # number female
  for (ii in 1:amax) stat=c(stat,length(which(nhs$age[sub]>= 10*(ii-1) & nhs$age[sub]<10*ii))) # number in 10 year age brackets up to amax
  stat=c(stat,length(which(nhs$age[sub]>=10*amax))) # number older than amax*10
  for (ii in 1:10) stat=c(stat,length(which(nhs$SIMD_DECILE_2016_SCT[sub]==ii))) # number with an A&E admission in lookback period
  stat=c(stat,length(which(nhs$ltc_total_count[sub]>0))) # number with a long-term condition
  stat=c(stat,length(which(nhs$target[sub]==TRUE))) # number with target positive
  for (ii in 1:length(varnames)) {
    X=nhs[[which(colnames(nhs)==varnames[ii])]][sub]
    if (grepl("since",varnames[ii])) w=which(X<= max(X,na.rm=T)-1) else w=which(X>0)
    stat=c(stat,length(w),mean(X[w],na.rm=T),median(X[w],na.rm=T),sd(X[w],na.rm=T))
  }
  xtab=rbind(xtab,stat)
  print(i)
}
rownames(xtab)=c("All","All_inc","adm","not_adm")
colnames(xtab)=c("Total","M","F",paste0("age_",10*(0:(amax-1)),"_",-1+10*(1:amax)),paste0("age_gt_",10*amax),
  paste0("SIMD",1:10),"LTC","Target",
  paste0(c("npos_","meanpos_","medpos_","sdpos_"),rep(varnames,each=4)))

sink(paste0(plot_dir,"Description/demographic_table.txt"))
cat(xtab)
cat(sink_marker)
cat("xtab="); dput(xtab); cat("\n\n")
sink()



# Exclusions
sink(paste0(plot_dir,"Description/exclusions.txt"))
cat("Total patient-time pairs under consideration\n")
cat(dim(nhs)[1])
cat("\n\n")
cat("Total patients recorded anywhere in data tables\n")
cat(length(unique(patients$id))); 
cat("\n\n")
cat("Total patients excl. all those without records prior to final time cutoff\n")
cat(length(unique(nhs$id)))
cat("\n\n")
cat("Total patients in study\n")
cat(length(unique(nhs$id[v34])))
cat("\n\n")
cat("Total records as above (patient,time pairs)\n")
cat(dim(nhs)[1])
cat("\n\n")
cat("Total records excl. individuals deceased prior to time cutoffs\n")
cat(length(which(!is.na(nhs$target))))
cat("\n\n")
cat("Total records excl. individuals with no v3 score\n")
cat(length(v34))
cat("\n\n")
cat("Total records for which indivudual was admitted and survived through predicion year\n")
cat(length(which(all_pred$reason=="E")))
cat("\n\n")
cat("Total records for which indivudual was admitted and later died in predicion year\n")
cat(length(which(all_pred$reason=="B")))
cat("\n\n")
cat("Total records for which indivudual died in prediction year prior to any admission\n")
cat(length(which(all_pred$reason=="D")))
cat("\n\n")
exclude_death=which(is.na(nhs$target))
exclude_v3=which(is.na(nhs$v3score))
exclude_simd=which(!(nhs$SIMD_DECILE_2016_SCT %in% 1:10))
exclude_list=which(nhs$id %in% read.csv(paste0("../../Linked Data/Unmatched_UPIs_without_UPI.csv.gz"))$UNIQUE_STUDY_ID)
cat("A: Number of records for which individual died before time cutoff\n")
cat(length(exclude_death))
cat("\n\n")
cat("B: Number of records for which individual had NA v3 score\n")
cat(length(exclude_v3))
cat("\n\n")
cat("C: Number of records for which individual had NA SIMD\n")
cat(length(exclude_simd))
cat("\n\n")
cat("D: Number of records for which individual was on exclusion list\n")
cat(length(exclude_list))
cat("\n\n")
cat("A and B\n")
cat(length(intersect(exclude_death,exclude_v3)))
cat("\n\n")
cat("A and C\n")
cat(length(intersect(exclude_death,exclude_simd)))
cat("\n\n")
cat("A and D\n")
cat(length(intersect(exclude_death,exclude_list)))
cat("\n\n")
cat("B and C\n")
cat(length(intersect(exclude_v3,exclude_simd)))
cat("\n\n")
cat("B and D\n")
cat(length(intersect(exclude_v3,exclude_list)))
cat("\n\n")
cat("C and D\n")
cat(length(intersect(exclude_simd,exclude_list)))
cat("\n\n")
cat("ABC\n")
cat(length(intersect(exclude_death,intersect(exclude_v3,exclude_simd))))
cat("\n\n")
cat("ABD\n")
cat(length(intersect(exclude_death,intersect(exclude_v3,exclude_list))))
cat("\n\n")
cat("ACD\n")
cat(length(intersect(exclude_death,intersect(exclude_simd,exclude_list))))
cat("\n\n")
cat("BCD\n")
cat(length(intersect(exclude_v3,intersect(exclude_simd,exclude_list))))
cat("\n\n")
cat("Number of unique individuals amongst admissions/deaths\n")
cat(length(unique(all_pred$id[which(all_pred$target==TRUE)])))
cat("\n\n")
cat("Number of unique individuals amongst non-admissions/non-deaths\n")
cat(length(unique(all_pred$id[which(all_pred$target==FALSE)])))
cat("\n\n")

sink()


nhs=nhs[v34,] # redefine nhs for memory considerations 

# Print stuff
sink(paste0(plot_dir,"Description/misc.txt"))
cat("Total patients\n")
cat(dim(patients)[1])
cat("\n\nTotal recorded episodes of any type\n")
cat(dim(episodes)[1])
cat("\n\nTotal patients used for training and testing\n")
cat(dim(nhs)[1])
cat("\n\nNumber of events in training and testing\n")
cat(sum(nhs$target))
cat("\n\nVariable names\n")
vx=setdiff(colnames(nhs),colnames(nhs)[grep("topic",colnames(nhs))])
for (i in 1:length(vx)) {
  cat(vx[i])
  cat("\n")
}
cat("\n\nSources\n")
for (i in 1:length(list_of_data_tables)) {
  cat("\n")
  cat(names(list_of_data_tables)[i])
  cat("\n")
  cat(dim(list_of_data_tables[[i]])[1])
  cat("\n")
  cat(length(unique(list_of_data_tables[[i]]$id)))
  cat("\n\n")
}

cat("\n\nSuper-learner coefficients\n")
m1=readRDS("James/Analysis/Data/full_model/models_cv/model_fold1.RDS")
m2=readRDS("James/Analysis/Data/full_model/models_cv/model_fold2.RDS")
m3=readRDS("James/Analysis/Data/full_model/models_cv/model_fold3.RDS")
g1=m1$slfit
g2=m2$slfit
g3=m3$slfit

cat("\n\nFold1\n")
cat("Lambda_min\n")
print(g1$glmnet.fit$beta[,which.min(abs(g1$lambda-g1$lambda.min))])

cat("\n\nFold2\n")
cat("Lambda_min\n")
print(g2$glmnet.fit$beta[,which.min(abs(g2$lambda-g2$lambda.min))])

cat("\n\nFold3\n")
cat("Lambda_min\n")
print(g3$glmnet.fit$beta[,which.min(abs(g3$lambda-g3$lambda.min))])



cat("\n\nCalibrating transforms\n")
x1=m1$ctransform$transform
x2=m2$ctransform$transform
x3=m3$ctransform$transform

cat("\n\nFold1\n")
cat("x_e="); dput(m1$x_e); cat("\n")
cat("y_e="); dput(m1$y_e); cat("\n\n\n")

cat("\n\nFold2\n")
cat("x_e="); dput(m2$x_e); cat("\n")
cat("y_e="); dput(m2$y_e); cat("\n\n\n")

cat("\n\nFold3\n")
cat("x_e="); dput(m3$x_e); cat("\n")
cat("y_e="); dput(m3$y_e); cat("\n\n\n")





### new
cat("\n\n\nMissingness\n")
g0=c(grep("days_since_last",colnames(nhs)),grep("yearssincediag",colnames(nhs)))
for (i in 1:dim(nhs)[2]) {
  if (i %in% g0) nmiss=0 else nmiss=length(which(is.na(nhs[,i]))) # different coding here; NA relevant
  cat(colnames(nhs)[i],", ",nmiss/dim(nhs)[1],"\n")
}
### end new

sink()




