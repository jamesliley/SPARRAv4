######################################################
## Code to demonstrate attrition of accuracy of     ##
##  SPARRA model with time                          ##
######################################################
##
## James Liley, 2020
##
## Run on windows, then linux, then windows again.



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
  library(ranger)
  library(glmnet)
  library(Metrics)
  library(tidyverse)
  library(glue)
  library(varhandle)
} else {
  library(h2o)
  library(fst)
}

# Directories
model_dir="James/Analysis/Data/time_attenuation/"
plot_dir="James/Analysis/Output/time_attenuation/"
data_dir="James/Analysis/Data/"

# File locations
topicfile=paste0(model_dir,"topic_model_time.rds")
trainmatfile=paste0(model_dir,"training_matrix_time_attenuation.rds")


# In text outputs, everything after this marker should be readable by R.
sink_marker="\n\n\n***************\n\n"
options(width=1e4)
options(max.print=1e7)

######################################################
## Flags                                            ##
######################################################

force_redo=FALSE


######################################################
## Time specifications                              ##
######################################################

# Training data uses this time cutoff
time_train=dmy_hm("1-5-2014 00:00")


######################################################
## Load data                                        ##
######################################################

if (!(file.exists(topicfile) & file.exists(trainmatfile)) | force_redo) {
load_cleanData(partition = "all",
	subset_frac=1, subset_seed = 1234,
	load_to_global=TRUE,
	load_sparrav3_scores = FALSE)
}


######################################################
## Fit topic model to prescription/diagnosis data   ##
##  collected before time_train                     ##
######################################################

if ((!file.exists(topicfile) | force_redo) & (location=="Windows")) {

library(topicmodels)
  
topic_train_time=topic_training_matrix(patients,episodes,list_of_data_tables,time_train,
	save_file=paste0(model_dir,"doc_term_long_table_training"))

# Restrict to individuals with valid v3 score and nonzero rows
sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
tx=as.numeric(sparraV3Scores[[2]])
sw=which(tx %in% as.numeric(time_train))
xv3=sparraV3Scores[sw,]; rm(list=c("sparraV3Scores","tx","sw"))
xv3=xv3[which(!is.na(xv3[[3]])),]
id=xv3[[1]]; ix=as.numeric(xv3[[2]])
idt=paste0(id,"_",ix) # slow

topic_train_time=topic_train_time[which(rownames(topic_train_time) %in% idt),]
topic_train_time=topic_train_time[which(rowSums(topic_train_time)>0),]

topic_model_time <- LDA(x = topic_train_time, k = 30)
saveRDS(topic_model_time,file=topicfile)

}

######################################################
## Transformation to training matrix                ##
######################################################


if ((!file.exists(trainmatfile) | force_redo) & (location=="Windows")) {

library(topicmodels)  
topic_model_time=readRDS(topicfile)
  
# Topics
topic_features <- transformer(patients, episodes, list_of_data_tables,
	list(
		all_codes_topic = list(topic_model_fit = topic_model_time)
	),
	as.integer(time_train),
	hard_max_lookback = 366)
topic_file=paste0(data_dir,"topic_temp.RDS")
saveRDS(topic_features,file=topic_file)
rm(topic_features)
gc()

# Long-term condition counts
ltc_features <- transformer(patients, episodes, list_of_data_tables,
	list(
		ltcs = list(output_type = "rawdata_NUMBEROFLTCs")
	),
	as.integer(time_train),
	hard_max_lookback = 366)
ltc_file=paste0(data_dir,"ltc_temp.RDS")
saveRDS(ltc_features,file=ltc_file)
rm(ltc_features)
gc()

# Time-since-event features
time_features <- transformer(patients, episodes, list_of_data_tables,
	list(
		last_episode_days_ago = list(source_table_names = c("AE2", "SMR00", "SMR01", "SMR01E", "SMR04")),
		last_emergency_admission_days_ago = list(),
		ltcs = list(output_type = "years_since_diag")
	),
	as.integer(time_train),
	hard_max_lookback = 366)
time_file=paste0(data_dir,"time_temp.RDS")
saveRDS(time_features,file=time_file)
rm(time_features)
gc()


# Features from SPARRAv3. We can only use these with a one-year lookback
v3_features <- transformer(patients, episodes, list_of_data_tables,
	list(
		sparrav3 = list()
	),
	as.integer(time_train),
	hard_max_lookback = 366)
v3_file=paste0(data_dir,"v3_temp.RDS")
saveRDS(v3_features,file=v3_file)
rm(v3_features)
gc()

# Some memory cleanup
#rm(list=c("list_of_data_tables","topic_model_time","patients","episodes"))
gc()

# Join
train_matrix = combine_training_matrices(v3_file,topic_file,time_file,ltc_file,
  patients,episodes,list_of_data_tables,
  keep_id_and_time=TRUE) 

## Remove extraneous variables relating to time cutoff
train_matrix = train_matrix %>% select(-c("topic_cutoff","ltc_cutoff"))

## Processing

# Add sex as predictor
train_matrix$sexM=as.integer((patients$gender=="Male")[match(train_matrix$id,patients$id)])

dsince=colnames(train_matrix)[grep("days_since",colnames(train_matrix))]
for (i in 1:length(dsince)) {
  x=train_matrix[[dsince[i]]]
  w=which(!is.finite(x)|(x==0))
  x[w]=2*366 # If the person has never had an episode, set value to 2* hardmax
  train_matrix[[dsince[i]]]=x
}
ysince=colnames(train_matrix)[grep("yearssince",colnames(train_matrix))]
for (i in 1:length(ysince)) {
  x=train_matrix[[ysince[i]]]
  w=which(!is.finite(x)|(x==0))
  x[w]=100 # If the person has never had an episode, set value to 100
  train_matrix[[ysince[i]]]=x
}





# Adjoin v3 scores
sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
tx=as.numeric(sparraV3Scores[[2]])
w3=which(tx %in% unique(as.numeric(train_matrix$time)))
id3=paste0(sparraV3Scores$id[w3],"_",tx[w3])
id4=paste0(train_matrix$id,"_",as.numeric(train_matrix$time))
train_matrix$v3score=sparraV3Scores$SPARRAv3_score[w3][match(id4,id3)]

rm(list=c("sparraV3Scores","tx","w3","id3","id4"))
rm(list=c("v3_features","time_features","ltc_features","topic_features"))
gc()

# Need to reload to ensure all individuals are in patient list
load_cleanData(partition = "all",
  subset_frac=1, subset_seed = 1234,
  load_to_global=TRUE,
  load_sparrav3_scores = FALSE)

train_matrix=train_matrix[which(train_matrix$id %in% patients$id),]
train_matrix=train_matrix[which(!is.na(train_matrix$v3score)),]
train_matrix=train_matrix[which(!is.na(train_matrix$target)),]

# Memory cleanup
rm(list=c("patients","episodes","list_of_data_tables"))
gc()


# Save
saveRDS(train_matrix,file=trainmatfile)

rm(list=c("train_matrix"))

}


######################################################
## Fit model                                        ##
######################################################

# sparrav4.fit.split saves predictions for Xpred to this file
time_model="James/Analysis/Data/time_attenuation/Model/time_model"

if (length(list.files(dirname(time_model),pattern="time_model\\."))<11) {
 train_matrix=readRDS(trainmatfile)
 train_id=train_matrix$id
 mod=SPARRAv4.fit.split(train_matrix %>% select(-c("id","time")),
	  seed=220,model_id=time_model,linux=(location=="Linux"))
}



######################################################
## Fit model with no topic features                 ##
######################################################


# sparrav4.fit.split saves predictions for Xpred to this file
time_model_notopics="James/Analysis/Data/time_attenuation/Model/time_model_notopics"


if (length(list.files(dirname(time_model_notopics),pattern="time_model_notopics\\."))<11) {
 train_matrix=readRDS(trainmatfile)
 train_id=train_matrix$id
 mod_notopics=SPARRAv4.fit.split(train_matrix %>% select(-c("id","time",matches("topic*"))),
    seed=220,model_id=time_model_notopics,linux=(location=="Linux"))
}




######################################################
## Test matrices                                    ##
######################################################

test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)

for (i in 1:length(test_times)) {

time_test=test_times[i]

testmatfile=paste0(data_dir,"time_attenuation/training_matrix_time_attenuation_",gsub(" ","_",time_test),".R")
if ((!file.exists(testmatfile) | force_redo) & (location=="Windows")) {
  
  if (!exists("patients")) {
  load_cleanData(partition = "all",
    subset_frac=1, subset_seed = 1234,
    load_to_global=TRUE,
    load_sparrav3_scores = FALSE)
  }
  if (!exists("topic_model_time")) {
    library(topicmodels)
    topic_model_time=readRDS(topicfile)
  }
  
  # Topics
  topic_features <- transformer(patients, episodes, list_of_data_tables,
    list(
      all_codes_topic = list(topic_model_fit = topic_model_time)
    ),
    as.integer(time_test),hard_max_lookback = 366)
  # Save: memory issues
  topic_file=paste0(data_dir,"topic_temp.RDS")
  saveRDS(topic_features,file=topic_file)
  rm(topic_features)
  gc()
  
  # Long-term condition counts
  ltc_features <- transformer(patients, episodes, list_of_data_tables,
    list(
      ltcs = list(output_type = "rawdata_NUMBEROFLTCs")
    ),
    as.integer(time_test),
    hard_max_lookback = 366)
  # Save: memory issues
  ltc_file=paste0(data_dir,"ltc_temp.RDS")
  saveRDS(ltc_features,file=ltc_file)
  rm(ltc_features)
  gc()
  
  
  # Time-since-event features
  time_features <- transformer(patients, episodes, list_of_data_tables,
    list(
      last_episode_days_ago = list(source_table_names = c("AE2", "SMR00", "SMR01", "SMR01E", "SMR04")),
      last_emergency_admission_days_ago = list(),
      ltcs = list(output_type = "years_since_diag")
    ),
    as.integer(time_test),
    hard_max_lookback = 366)
  # Save: memory issues
  time_file=paste0(data_dir,"time_temp.RDS")
  saveRDS(time_features,file=time_file)
  rm(time_features)
  gc()
  
  
  # Features from SPARRAv3. We can only use these with a one-year lookback
  v3_features <- transformer(patients, episodes, list_of_data_tables,
    list(sparrav3 = list()),
    as.integer(time_test),
    hard_max_lookback = 366)
  # Save: memory issues
  v3_file=paste0(data_dir,"v3_temp.RDS")
  saveRDS(v3_features,file=v3_file)
  rm(v3_features)
  gc()
  
  # Some memory cleanup
  rm(list=c("topic_model_time"))
  gc()
  
  # Join
  test_matrix = combine_training_matrices(v3_file,topic_file,time_file,ltc_file,
    patients,episodes,list_of_data_tables,
    keep_id_and_time=TRUE) 
  
  ## Remove extraneous variables relating to time cutoff
  test_matrix = test_matrix %>% select(-c("topic_cutoff","ltc_cutoff"))

  ## Processing
  
  # Add sex as predictor
  test_matrix$sexM=as.integer((patients$gender=="Male")[match(test_matrix$id,patients$id)])
  
  dsince=colnames(test_matrix)[grep("days_since",colnames(test_matrix))]
  for (ii in 1:length(dsince)) {
    x=test_matrix[[dsince[ii]]]
    w=which(!is.finite(x)|(x==0))
    x[w]=2*366 # If the person has never had an episode, set value to 2* hardmax
    test_matrix[[dsince[ii]]]=x
  }
  ysince=colnames(test_matrix)[grep("yearssince",colnames(test_matrix))]
  for (ii in 1:length(ysince)) {
    x=test_matrix[[ysince[ii]]]
    w=which(!is.finite(x)|(x==0))
    x[w]=100 # If the person has never had an episode, set value to 100
    test_matrix[[ysince[ii]]]=x
  }
  
    
  # Adjoin v3 scores
  sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
  tx=as.numeric(sparraV3Scores[[2]])
  w3=which(tx %in% unique(as.numeric(test_matrix$time)))
  id3=paste0(sparraV3Scores$id[w3],"_",tx[w3])
  id4=paste0(test_matrix$id,"_",as.numeric(test_matrix$time))
  test_matrix$v3score=sparraV3Scores$SPARRAv3_score[w3][match(id4,id3)]
  
  rm(list=c("sparraV3Scores","tx","w3","id3","id4"))
  gc()
  
  # Need to reload to ensure all individuals are in patient list
  #load_cleanData(partition = "all",
  #  subset_frac=1, subset_seed = 1234,
  #  load_to_global=TRUE,
  #  load_sparrav3_scores = FALSE)
  
  test_matrix=test_matrix[which(test_matrix$id %in% patients$id),]
  test_matrix=test_matrix[which(!is.na(test_matrix$v3score)),]
  test_matrix=test_matrix[which(!is.na(test_matrix$target)),]
  
  # Memory cleanup
  #rm(list=c("patients","episodes","list_of_data_tables"))
  gc()

  # Save
  saveRDS(test_matrix,file=testmatfile)
  rm(list=c("test_matrix"))
  gc()
  
  #assign(paste0("test_id",i),test_matrix$id)
  #assign(paste0("test_matrix",i),test_matrix %>% select(-c("id")))
  
  print(paste0("Constructed test matrix for time ",test_times[i]))

} 
}

for (i in 1:length(test_times)) {
  time_test=test_times[i]
  testmatfile=paste0(data_dir,"time_attenuation/training_matrix_time_attenuation_",gsub(" ","_",time_test),".R")
  testx=readRDS(testmatfile)
  assign(paste0("test_matrix",i),testx)
  assign(paste0("test_id",i),testx$id)
  rm(list=c("testx"))
}



######################################################
## Combine test matrices                            ##
######################################################

# 'Test' matrix is combined test matrices for various times
test_matrix_combined=test_matrix1
test_matrix_index=rep(1,dim(test_matrix1)[1])
rm(test_matrix1)
for (i in 2:length(test_times)) {
  tmx=get(paste0("test_matrix",i))
  test_matrix_combined=rbind(test_matrix_combined,tmx)
  test_matrix_index=c(test_matrix_index,rep(i,dim(tmx)[1]))
  rm(list=c("tmx",paste0("test_matrix",i)))
}
test_id=test_matrix_combined$id
test_time=test_matrix_combined$time
test_age=test_matrix_combined$age
gc()



######################################################
## Make predictions                                 ##
######################################################

pred_file="James/Analysis/Data/time_attenuation/Model/time_model.output.RDS"
pred_file_notopics="James/Analysis/Data/time_attenuation/Model/time_model_notopics.output.RDS"

if (!file.exists(paste0(plot_dir,"topics/topic_comparison.txt"))) {

output_name=paste0(time_model,".output")
pred=predict.sparra(test_matrix_combined %>% select(-c("id","time")),
  time_model,save_id=output_name,verbose=T,
  linux=(location=="Linux"))

output_name_notopics=paste0(time_model_notopics,".output")
  pred_notopics=predict.sparra(test_matrix_combined %>% select(-c("id","time",matches("topic*"))),
    time_model_notopics,save_id=output_name_notopics,verbose=T,
  linux=(location=="Linux"))

}


######################################################
## Read in predictions                              ##
######################################################

Ypred=readRDS(pred_file)
ypredtest=data.frame(Ypred)
colnames(ypredtest)[dim(ypredtest)[2]]="super"
ytest=test_matrix_combined$target
rm(list=c("Ypred","test_matrix_combined"))

Ypredn=readRDS(pred_file_notopics)
ypredtestn=data.frame(Ypredn)
colnames(ypredtestn)[dim(ypredtestn)[2]]="super"
rm(list=c("Ypredn"))



######################################################
## Topics vs no topics                              ##
######################################################

## Use this subset of test matrix
sub=which(test_matrix_index==1)


## Test and practical difference
library(pROC)
roc_topic=roc(ytest[sub],ypredtest$super[sub])
roc_notopic=roc(ytest[sub],ypredtestn$super[sub])

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
ord_topic=order(-ypredtest$super[sub])
ord_notopic=order(-ypredtestn$super[sub])
for (i in 1:length(topn)) {
  n_top=sum(ytest[sub][which(ord_topic <= topn[i])])
  n_notop=sum(ytest[sub][which(ord_notopic <= topn[i])])
  out=cbind(out,c(n_top,n_notop,n_top-n_notop))
}
out=rbind(topn,out)
rownames(out)=c("TopN","N_adm_topics","N_adm_notopics","Difference")
print(out)

cat(sink_marker)
cat("r0="); dput(r0); cat("\n\n")
cat("out="); dput(out); cat("\n\n")
sink()


## ROC curves

xrocT=getroc(ytest[sub],ypredtest$super[sub])
xrocNT=getroc(ytest[sub],ypredtestn$super[sub]) 

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

xprcT=getprc(ytest[sub],ypredtest$super[sub])
xprcNT=getprc(ytest[sub],ypredtestn$super[sub]) 

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
xcalT=plotcal(ytest[sub],ypredtest$super[sub],plot=FALSE,kernel=TRUE)
xcalNT=plotcal(ytest[sub],ypredtestn$super[sub],plot=FALSE,kernel=TRUE)

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





######################################################
## Plot of change in accuracy with time             ##
######################################################

# There are systematically different admission rates year by year, so we resample to an equal admission rate
mt=mean(ytest) # Grand mean of target probability
nt=1000000 # use 1 million samples
nt1=round(nt*mt) # sample this many individuals with AE
nt0=nt-nt1 # and this many without
seed=220 # make sure these are consistent



for (i in 1:dim(ypredtest)[2]) {
xname=colnames(ypredtest)[i]

## Get plot co-ordinates
for (time in 1:max(test_matrix_index)) {
  set.seed(seed)
  ind=which(test_matrix_index==time)
  s1=sample(ind[which(ytest[ind]==1)],nt1)
  s0=sample(ind[which(ytest[ind]==0)],nt0)
  indx=c(s0,s1)
  xy=getroc(ytest[indx],ypredtest[indx,i])
  px=getprc(ytest[indx],ypredtest[indx,i])
  cx=plotcal(ytest[indx],ypredtest[indx,i],n=20,kernel=T,plot=F)
  xroc=roc(ytest[indx],ypredtest[indx,i])
  assign(paste0("i",time),ind)
  assign(paste0("xy",time),xy)
  assign(paste0("px",time),px)
  assign(paste0("cx",time),cx)
  assign(paste0("xroc",time),xroc)
}
ntx=length(test_times)



# Number of time points
ntx=5

## ROC curves for predictor at various time points
pdf(paste0(plot_dir,"recalculate/ROC/roc.",xname,".pdf"),width=3,height=3.5)
roc_2panel(list(xy1,xy2,xy3,xy4,xy5),labels=as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

sink(paste0(plot_dir,"recalculate/ROC/roc.",xname,".txt"))
for (i in 1:ntx) {
 xy=get(paste0("xy",i))
 cat("Time",i,"\n","Sensitivity\n",xy$sens,"\n\nSpecificity\n",xy$spec,"\n\nAUROCs\n",xy$auc,"\n\nSE\n",xy$se,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  xy=get(paste0("xy",i))
  cat(paste0("xy",i,"=")); dput(xy); cat("\n\n")
}
sink()
 





 

## PR curves for predictor at various time points
pdf(paste0(plot_dir,"recalculate/PRC/prc.",xname,".pdf"),width=3,height=3.5)
prc_2panel(list(px1,px2,px3,px4,px5),as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

sink(paste0(plot_dir,"recalculate/PRC/prc.",xname,".txt"))
for (i in 1:ntx) {
  px=get(paste0("px",i))
  cat("Time",i,"\n","Sensitivity\n",px$sens,"\n\nPPV\n",px$ppv,"\n\nAUROCs\n",px$auc,"\n\nSE\n",px$se,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  px=get(paste0("px",i))
  cat(paste0("px",i,"=")); dput(px); cat("\n\n")
}
sink()


    
  
## Calibration curves for predictor at various time points
pdf(paste0(plot_dir,"recalculate/CAL/cal.",xname,".pdf"),width=3,height=3.5)
cci=rgb((ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),alpha=0.5) # colour for confidence envelope
cal_2panel(list(cx1,cx2,cx3,cx4,cx5),as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)),
  ci_col=c(NA,NA,NA,NA,cci))
dev.off()

sink(paste0(plot_dir,"recalculate/CAL/cal.",xname,".txt"))
for (i in 1:ntx) {
  cx=get(paste0("cx",i))
  cat("Time",i,"\n\n","Obs\n",cx$x,"\n\nPred\n",cx$y,"\n\nLower\n",cx$lower,"\n\nUpper\n",cx$upper,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  cx=get(paste0("cx",i))
  cat(paste0("cx",i,"=")); dput(cx); cat("\n\n")
}
sink()

}


## Compare AUROCs for super-learner
sink(paste0(plot_dir,"recalculate/comparisons.txt"))

for (i in 1:max(test_matrix_index)) {
 for (j in 1:max(test_matrix_index)) {
   if (i<j) {
     roc_t1=get(paste0("xroc",i))
     roc_t2=get(paste0("xroc",j))
     cat("Difference_timepoint_",i,"_",as.character(test_times[i]),"_vs_timepoint_",j,"_",as.character(test_times[j]),"\n\n")
     rt=roc.test(roc_t1,roc_t2)
     print(rt)
     rt$roc1=NULL; rt$roc2=NULL
     assign(paste0("rt_",i,"_",j),rt)
     cat("\n\n\n")
   }
 }
}
cat(sink_marker)
for (i in 1:max(test_matrix_index)) {
  for (j in 1:max(test_matrix_index)) {
    if (i<j) {
      rt=get(paste0("rt_",i,"_",j))
      cat(paste0("rt_",i,"_",j,"=")); dput(rt); cat("\n\n\n")
    }
  }
}
sink()




######################################################
## Performance of the same predictor over time      ##
######################################################

# We need a cohort of patients at each time point with similar age distributions

# Auxiliary function: subsample to match age distribution and proportion of admissions in each age.
#  Assume a reasonable number of admissions for each age.
match_sub=function(age,target,agesub,targetsub) {
  t1a=table(age[which(target==0)]); t1b=table(age[which(target==1)])
  t2a=table(agesub[which(targetsub==0)]); t2b=table(agesub[which(targetsub==1)])
  t1=as.vector(c(t1a,t1b)); t2=as.vector(c(t2a,t2b)) 
  t1x=t1/sum(t1); t2x=t2/sum(t2)
  nw=floor(min(t2/t1x)) # find the size of the largest possible matched sample; w is the 'weak point'
  sub=c()
  min_age=min(age); max_age=max(age)
  for (i in min(age):max(age)) sub=c(sub,sample(which(agesub==i & targetsub==0),floor(t1x[i+1-min_age]*nw)))
  for (i in min(age):max(age)) sub=c(sub,sample(which(agesub==i & targetsub==1),round(t1x[max_age-2*min_age+i+2]*nw)))
  return(sub)
}

# Age distribution across whole daaset (excluding very young and old). We will match to this.
age=test_age; target=ytest
w=which(age>10 & age<90); age=age[w]; target=target[w]


# Lookup
for (i in 1:length(test_times)) assign(paste0("id",i),which(test_matrix_index==i))
seed=220 # make sure these are consistent

ntx=length(test_times)
for (time in 1:length(test_times)) {
  set.seed(seed)
  eid=test_id1; for (ii in 1:time) eid=intersect(eid,get(paste0("test_id",ii)))  # Eligible IDS
  eind=get(paste0("id",time))[match(eid,get(paste0("test_id",time)))] # Indices of eligible IDs at time point 'time'
  
  e_age0=test_age[eind]
  w=which(e_age0>10 & e_age0<90)
  
  eind=eind[w] # trim to 10<age<90
  
  e_age=test_age[eind]; e_target=ytest[eind]
  indx=eind[match_sub(age,target,e_age,e_target)] # match by age and target
  
  ind1=id1[match(test_id[indx],test_id1)]
  
  xy=getroc(ytest[indx],ypredtest$super[ind1])
  px=getprc(ytest[indx],ypredtest$super[ind1])
  cx=plotcal(ytest[indx],ypredtest$super[ind1],n=20,kernel=T,plot=F)
  xroc=roc(ytest[indx],ypredtest$super[ind1])
  assign(paste0("xy",time),xy)
  assign(paste0("px",time),px)
  assign(paste0("cx",time),cx)
  assign(paste0("xroc",time),xroc)
  assign(paste0("indx",time),indx)
}

# Number of time points
ntx=5

## ROC curves for predictor at various time points
pdf(paste0(plot_dir,"cohort/drift/drift_roc.pdf"),width=3,height=3.5)
roc_2panel(list(xy1,xy2,xy3,xy4,xy5),labels=as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

sink(paste0(plot_dir,"cohort/drift/drift_roc.txt"))
for (i in 1:ntx) {
  xy=get(paste0("xy",i))
  cat("Time",i,"\n","Sensitivity\n",xy$sens,"\n\nSpecificity\n",xy$spec,"\n\nAUROCs\n",xy$auc,"\n\nSE\n",xy$se,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  xy=get(paste0("xy",i))
  cat(paste0("xy",i,"=")); dput(xy); cat("\n\n")
}
sink()


## PRC curves for predictor at various time points
pdf(paste0(plot_dir,"cohort/drift/drift_prc.pdf"),width=3,height=3.5)
prc_2panel(list(px1,px2,px3,px4,px5),labels=as.character(test_times),
  col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()
sink(paste0(plot_dir,"cohort/drift/drift_prc.txt"))
for (i in 1:ntx) {
  px=get(paste0("px",i))
  cat("Time",i,"\n","Sensitivity\n",px$sens,"\n\nPPV\n",px$ppv,"\n\nAUROCs\n",px$auc,"\n\nSE\n",px$se,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  px=get(paste0("px",i))
  cat(paste0("px",i,"=")); dput(px); cat("\n\n")
}
sink()


## Calibration curves for predictor at various time points
pdf(paste0(plot_dir,"cohort/drift/drift_cal.pdf"),width=5,height=5)
cci=rgb((ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),alpha=0.5)
cal_2panel(list(cx1,cx2,cx3,cx4,cx5),labels=as.character(test_times),
  col=gray(((1:ntx)-1)/(1.3*ntx)),ci_col=c(NA,NA,NA,NA,cci))
dev.off()
sink(paste0(plot_dir,"cohort/drift/drift_cal.txt"))
for (i in 1:ntx) {
  cx=get(paste0("cx",i))
  cat("Time",i,"\n\n","Obs\n",cx$x,"\n\nPred\n",cx$y,"\n\nLower\n",cx$lower,"\n\nUpper\n",cx$upper,"\n\n\n")
}
cat(sink_marker)
for (i in 1:ntx) {
  cx=get(paste0("cx",i))
  cat(paste0("cx",i,"=")); dput(cx); cat("\n\n")
}
sink()




# Recent vs earlier predictions
pdf(paste0(plot_dir,"cohort/auroc_comparison.pdf"),width=5,height=5)
ntx=5
cols=gray(((1:ntx)-1)/(1.3*ntx))
ltyx=c(1,2,3,4,5)

plot(0,type="n",xlim=c(1,ntx),ylim=c(0.755,0.791),xlab="Time point",ylab=expression(paste("AUROC")))
dx=matrix(0,ntx,ntx)
for (i in 1:ntx) {
  for (j in i:ntx) {
    # Prediction made for timepoint j using prediction from time point i
    indxj=get(paste0("indx",j))
    indxi=get(paste0("id",i))[match(test_id[indxj],get(paste0("test_id",i)))]
    dx[i,j]=getroc(ytest[indxj],ypredtest$super[indxi])$auc
  }
  if (i<ntx) lines(i:ntx,dx[i,i:ntx],col=cols[i],lty=ltyx[i]) else points(ntx,dx[ntx,ntx],col=cols[i],pch=16)
} 
legend("bottomleft",title="Predictor for:", legend=1:5, #as.character(test_times),
  col=cols,lty=c(1:4,NA),pch=c(rep(NA,ntx-1),16))
dev.off()
sink(paste0(plot_dir,"cohort/auroc_comparison.txt"))
print(dx)
cat(sink_marker)
cat("dx="); dput(dx); cat("\n\n")
sink()


######################################################
## Scores for fixed individuals over time           ##
######################################################

# Individuals with predictions at all time points
idx=test_id1
for (i in 2:length(test_times)) idx=intersect(idx,get(paste0("test_id",i)))
for (i in 1:length(test_times)) assign(paste0("id",i),which(test_matrix_index==i))
predmat=c()  
for (i in 1:length(test_times)) {
  idi=get(paste0("id",i)); txi=get(paste0("test_id",i))
  predmat=cbind(predmat,ypredtest$super[idi][match(idx,txi)]) # attach predictions
}
for (i in 1:length(test_times)) {
  idi=get(paste0("id",i)); txi=get(paste0("test_id",i))
  predmat=cbind(predmat,ytest[idi][match(idx,txi)]) # attach targets
}
colnames(predmat)=c(paste0("predicted",1:length(test_times)),paste0("target",1:length(test_times)))
predmat=as.data.frame(predmat)



######################################################
## Summary statistics                               ##
######################################################

sink(paste0(plot_dir,"cohort/summary.txt"))

# Basic
cat("Number of individuals: ",length(idx),"\n\n")
for (i in 1:length(test_times)) {
  cat("Number of admissions for time point",i,"(",as.character(test_times[i]),"): ",sum(predmat[,length(test_times)+i]),"\n\n")
}

# AUC for five-year target
cat("AUROC for predicting at least 1 admission in next 3y: ")
xtarget=pmax(predmat$target1,predmat$target2,predmat$target3,predmat$target4,predmat$target5)
roc3=(getroc(xtarget,predmat$predicted5))
print(roc3$auc)
cat("\n\n")
cat(sink_marker)
cat("roc3="); dput(roc3)
sink()




######################################################
## Cohort-based changes in risk                     ##
######################################################

nquant=100
ntx=length(test_times)
qmatc=matrix(0,nquant,ntx)
for (i in 1:ntx) qmatc[,i]=quantile(predmat[,i],((1:nquant)/nquant) - (1/(2*nquant)))



pdf(paste0(plot_dir,"cohort/change/cohort_change_over_time.pdf"),width=4,height=4)
par(mar=c(4,4,0.5,0.5))
nquant=100
ntx=length(test_times)
plot(0,type="n",xlim=c(1,ntx),ylim=c(0,0.05+max(qmatc)),
  xaxt="n",xlab="Time cutoff",
  ylab="P(EA)")
for (i in 0:10) abline(h=i/10,col="lightgray",lty=2)
axis(1,labels=c(expression("t"[1]),expression("t"[2]),expression("t"[3]),expression("t"[4]),expression("t"[5])),
  at=1:length(test_times))
for (i in 1:nquant) lines(1:ntx,qmatc[i,],col="lightblue",lty=1)
for (i in 1:10) lines(1:ntx,qmatc[round(nquant*i/10),],col="blue",lty=1)
dev.off()
sink(paste0(plot_dir,"cohort/change/cohort_change_over_time.txt"))
print(qmatc)
cat(sink_marker)
cat("qmatc="); dput(qmatc); cat("\n\n")
sink()

# Personal risk change with time by risk decile at first cutoff
qx=quantile(predmat[,1],(0:nquant)/nquant)
qmat=matrix(0,nquant+1,ntx)
for (i in 1:nquant) {
  w=which(predmat[,1]>qx[i] & predmat[,1]<=qx[i+1])
  for(j in 1:ntx) {
    qmat[i,j]=median(predmat[w,j])
  }
}


pdf(paste0(plot_dir,"cohort/change/individual_change_over_time.pdf"),width=4,height=4)
par(mar=c(4,4,0.5,0.5))
plot(0,type="n",xlim=c(1,ntx),ylim=c(0,0.05+max(qmat)),
  xaxt="n",xlab="Time cutoff",
  ylab="P(EA)")
for (i in 0:10) abline(h=i/10,col="lightgray",lty=2)
axis(1,labels=c(expression("t"[1]),expression("t"[2]),expression("t"[3]),expression("t"[4]),expression("t"[5])),
  at=1:length(test_times))
for (i in 1:nquant) lines(1:ntx,qmat[i,],col="lightblue",lty=1)
for (i in 1:10) lines(1:ntx,qmat[round(nquant*i/10),],col="blue",lty=1)
dev.off()
sink(paste0(plot_dir,"cohort/change/individual_change_over_time.txt"))
print(qmat)
cat(sink_marker)
cat("qmat="); dput(qmat); cat("\n\n")
sink()


######################################################
## How similar are predictions year-to-year?        ##
######################################################

# Kernel density estimates normalised to unit diagonal
xcol=colorRampPalette(c("white","darkgray","red","yellow"))(100)
xx=1; nn=100; sp=7; scmax=100
qplot= ((1:1000)/1000)^(1/4) # Add q-q plots at these quantiles. >500 samples in each bin.
library(matrixStats)


for (i in 1:length(test_times)) {
  for (j in 1:length(test_times)) {
   if (i<j) {
     pi=predmat[,i]; pj=predmat[,j]; 
     
     w=which(pi<scmax/100 & pj<scmax/100)
     
     pi=pi[w]; pj=pj[w]
     x=round(nn*pi); y=round(nn*pj)
     xout=xtabs(~xi+xj,cbind(xi=c(x,1:scmax),xj=c(y,1:scmax)))
     xout[which(xout<4)]=0; xout[which(xout>3 & xout<7)]=6 # DO NOT DELETE: important for privacy
     xout0=xout

     
     # Marginals
     d1=density(pi); d2=density(pj)
     
     pdf(paste0(plot_dir,"cohort/crosstime/time",i,"_vs_time",j,".pdf"),width=7,height=7)

     
     xout=xout0
     
     # Normalise columns, then rows
     for (ii in 1:dim(xout)[2]) xout[ii,]=xout[ii,]/sum(xout[ii,],na.rm=T);
     for (ii in 1:dim(xout)[1]) xout[,ii]=xout[,ii]/sum(xout[,ii],na.rm=T);  
     xmax=min(colMaxs(xout[,10:(scmax-15)],na.rm=T));
     xmin=min(xout[which(xout>0)],na.rm=T); ncol0=5
     
     image((1:scmax)/100,(1:scmax)/100,xout,xlab=paste0("Prediction for timepoint ",i,", ",as.character(test_times[i])),
       ylab=paste0("Prediction for timepoint ",j,", ",as.character(test_times[j])),
       col=xcol,breaks=c(0,rep(xmin/2,ncol0),seq(xmin,4*xmax/3,length=length(xcol)-ncol0-1),1+max(xout,na.rm=T)),
       xaxs="i",yaxs="i",xlim=c(-1/sp,scmax/100),ylim=c(-1/sp,scmax/100),bty="n")
     
     polygon(d1$x,d1$y/(sp*max(d1$y))-(1/sp),col="gray"); 
     polygon(d2$y/(sp*max(d2$y))-(1/sp),d2$x,col="gray")
     
     q0=quantile(pi,qplot)
     q1=quantile(pj,qplot)
     q0a=1:scmax; q1a=q0a; for (ii in 1:scmax) q1a[ii]=median(y[which(x==ii)],na.rm=T); q0a=q0a/scmax; q1a=q1a/scmax
     abline(0,1,col="black",lty=2,lwd=2)
     #abline(0,1,col="white",lty=1,lwd=2)
     
     #lines(q0,q1,col="black",lwd=4)
     lines(q0,q1,col="blue",lwd=2)
     
     #lines(q0a,q1a,col="black",lwd=4)
     lines(q0a,q1a,col="darkgreen",lwd=2)
     
     legend("bottomright",c("Density (low)","", "Density (high)", "Marginal", "X-Y","Q-Q","Median"),
       pch=c(16,16,16,NA,NA,NA,NA),lty=c(NA,NA,NA,1,2,1,1),lwd=c(NA,NA,NA,1,2,2,2),
       col=c("darkgray","red","yellow","black","black","blue","darkgreen"))
     dev.off()
     sink(paste0(plot_dir,"cohort/crosstime/time",i,"_vs_time",j,".txt"))
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
     cat("xout0="); dput(xout0); cat("\n\n")
     cat("d1="); dput(d1); cat("\n\n")
     cat("d2="); dput(d2); cat("\n\n")
     cat("q0="); dput(q0); cat("\n\n")
     cat("q1="); dput(q1); cat("\n\n")
     cat("q0a="); dput(q0a); cat("\n\n")
     cat("q1a="); dput(q1a); cat("\n\n")
     sink()
   }
}
}


