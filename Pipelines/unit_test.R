######################################################
## Code to run unit test. Generates mockup EHR      ##
##  tables, and fits a model to them.               ##
######################################################
##
## James Liley, 2020
##
## INSTRUCTIONS TO RUN ON NSH
## 1. Run up to end of 'Transformation to training and test matrices' on Windows machine
## 2. Run up to end of 'Check R and package versions' and then from start to end of  'Fit models and predict' on Linux machine
## 3. Run up to end of 'Packages and scripts', then 'Fit models and predict', then remainder of script on Windows machine.


######################################################
## Clear all variables and packages, set seed       ##
######################################################

rm(list=ls()) # remove all variables

# Random seed, used throughout
seed=220

######################################################
## Directories                                      ##
######################################################
##
## This section can be modified across settings

# Working directory and location (windows or linux machine)
if (file.exists("C:/Users/jliley_1718-0370/Documents")) {
	location="Windows"
	setwd("//farr-fs1/Study Data/1718-0370/Research/Code/")
	RNGkind(sample.kind="Rounding")
} else {
	location="Linux"
	setwd("/Study_Data/Research/Code/")
}

# Directories. These directories should all exist
dir_rawData = "James/Analysis/Testing/raw_data/" # directory in which mockup EHR tables will go
dir_cleanData = "James/Analysis/Testing/clean_data/" # directory where 'patients', 'episodes', 'list_of_data_tables' will go
dir_model="James/Analysis/Testing/model/" # directory where model files will go
dir_output="James/Analysis/Testing/output/" # directory where output and training/test matrices will go
dir_sparra="James/SPARRAv4/" # Directory containing scripts for SPARRAv4 functions

if (!all(file.exists(c(dir_rawData,dir_cleanData,dir_model,dir_output,dir_sparra)))) stop("Ensure appropriate directories exist")

# Hard max lookback
hardmax=366

######################################################
## Packages and scripts                             ##
######################################################

# Scripts
lx=list.files(dir_sparra,pattern="*.R",full.names=TRUE)
for (i in 1:length(lx)) source(lx[i])

# Packages
if (location=="Windows") {
	library(xgboost)
	library(glmnet)
	library(data.table)
	library(glue)
	library(ranger)
  library(topicmodels)
  library(varhandle)
} else {
	library(h2o)
	library(tidyverse)
	library(fst)
}


######################################################
## Check R and package versions                     ##
######################################################


# Check R versionon
if ((location=="Windows" & !(R.Version()$version.string=="R version 3.6.1 (2019-07-05)"))|
    (location=="Linux" & !(R.Version()$version.string=="R version 3.5.1 (2018-07-02)")))
  stop("R on Windows should be version 3.6.1") else message("R version correct")

# Check packages
pkgcheck=check_package_versions(linux=(location=="Linux"))
if (!pkgcheck[[1]]) stop(paste0("Packages have incorrect versions: packages ",paste0(pkgcheck$failed,collapse=", ")," should be versions ",paste0(pkgcheck$shouldbe,collapse=", ")," respectively")) else message("All installed package versions correct")


######################################################
## Mockup EHR tables                                ##
######################################################

if (location=="Windows") {

ehr_mockup(dir_rawData,seed=seed)

## Check reproducibility
set.seed(seed)
AE2 = haven::read_sav(file.path(dir_rawData, "Final_AE2_extract_incl_UniqueStudyID.zsav"))
rep_check=as.matrix(AE2)[cbind(sample(dim(AE2)[1],5),sample(dim(AE2)[2],5))]
rep_correct=c( "01Z","Z665T","28-02-2014","22:27","5:16"  )

if (!all(rep_check==rep_correct)) stop("EHR mockup reproducibility failed.") else message("EHR mockup correct")

}


######################################################
## Clean and load raw data                          ##
######################################################

if (location=="Windows") {

suppressWarnings(clean_rawData(force_redo=TRUE,dir_rawData=dir_rawData,dir_cleanData=dir_cleanData))

## Check reproducibility
set.seed(seed)
pts = read.fst(file.path(dir_cleanData, "patients.fst"))
rep_check=as.matrix(pts)[cbind(sample(dim(pts)[1],5),sample(dim(pts)[2],5,rep=TRUE))]
rep_correct=c(" 8","2","4213","0",NA)

w=which(is.na(rep_correct)); w2=setdiff(1:5,w)
if (!(all(is.na(rep_check[w])) & all(rep_check[w2]==rep_correct[w2]))) 
  stop("Data cleaning reproducibility failed.") else message("Data cleaning correct")

}

######################################################
## Load clean data                                  ##
######################################################

if (location=="Windows") {

load_cleanData(partition = "all",
	subset_frac=1, subset_seed = 1234,
	load_to_global=TRUE,
	load_sparrav3_scores = FALSE,
	dir_cleanData = dir_cleanData)


# Time cutoff
time_cutoff=dmy_hm("1-5-2015 00:00")

}

######################################################
## Topic model                                      ##
######################################################

if (location=="Windows") {

topic_matrix=topic_training_matrix(patients,episodes,list_of_data_tables,time_cutoff,
	save_file=paste0(dir_output,"topic_matrix"))
topic_matrix=topic_matrix[which(rowSums(topic_matrix)>0),]

topic_id=unlist(lapply(strsplit(rownames(topic_matrix),"_"),function(x) x[1]))
topic_cv=patients$cv[match(topic_id,patients$id)]

topic_train=topic_matrix[which(topic_cv==1),]
topic_test1=topic_matrix[which(topic_cv==2),]
topic_test2=topic_matrix[which(topic_cv==2),]

topic_model <- LDA(x = topic_train, k = 30,control=list(seed=seed))
saveRDS(topic_model,file=paste0(dir_output,"topic_model.RDS"))


# Check reproducibility
set.seed(seed)
tb=topic_model@beta
topic_check=tb[cbind(sample(1:dim(tb)[1],5,rep=T),sample(1:dim(tb)[2],5,rep=T))]
topic_correct=c(-100.0000, -100.0000, -141.3287, -180.5888, -705.0263)
if (!all(abs(topic_check-topic_correct)<1e-4)) stop("Reproducibility failed for topic model") else message("Topic model correct")

}

######################################################
## Transformation to training and test matrices     ##
######################################################


if (location=="Windows") {


# Topics
topic_features <- transformer(patients, episodes, list_of_data_tables,
  list(
    all_codes_topic = list(topic_model_fit = topic_model)
  ),
  as.integer(time_cutoff),
  hard_max_lookback = hardmax)


# Long-term condition counts
ltc_features <- transformer(patients, episodes, list_of_data_tables,
  list(
    ltcs = list(output_type = "rawdata_NUMBEROFLTCs")
  ),
  as.integer(time_cutoff),
  hard_max_lookback = hardmax)

# Time-since-event features
time_features <- transformer(patients, episodes, list_of_data_tables,
  list(
    last_episode_days_ago = list(source_table_names = c("AE2", "SMR00", "SMR01", "SMR01E", "SMR04")),
    last_emergency_admission_days_ago = list(),
    ltcs = list(output_type = "years_since_diag")
  ),
  as.integer(time_cutoff),
  hard_max_lookback = hardmax)


# Features from SPARRAv3. We can only use these with a one-year lookback
v3_features <- transformer(patients, episodes, list_of_data_tables,
  list(
    sparrav3 = list()
  ),
  as.integer(time_cutoff),
  hard_max_lookback = hardmax)


# Join
data_matrix = combine_training_matrices(v3_features,topic_features,time_features,ltc_features,
  patients,episodes,list_of_data_tables,
  keep_id_and_time=TRUE) 

## Remove extraneous variables relating to time cutoff
data_matrix = data_matrix %>% select(-c("topic_cutoff","ltc_cutoff"))

# Add sex as predictor
data_matrix$sexM=as.integer((patients$gender=="Male")[match(data_matrix$id,patients$id)])

# Fix 'Days since/Years since' style predictors'
dsince=colnames(data_matrix)[grep("days_since",colnames(data_matrix))]
for (i in 1:length(dsince)) {
  x=data_matrix[[dsince[i]]]
  w=which(!is.finite(x)|(x==0))
  x[w]=2*hardmax # If the person has never had an episode, set value to 2* hardmax
  data_matrix[[dsince[i]]]=x
}
ysince=colnames(data_matrix)[grep("yearssince",colnames(data_matrix))]
for (i in 1:length(ysince)) {
  x=data_matrix[[ysince[i]]]
  w=which(!is.finite(x)|(x==0))
  x[w]=100 # If the person has never had an episode, set value to 100
  data_matrix[[ysince[i]]]=x
}

# Adjoin (random) v3 scores
set.seed(seed)
v3score=sample(1:99,dim(data_matrix)[1],replace=T)
data_matrix$v3score=v3score

# Cleanup
data_matrix=data_matrix[which(!is.na(data_matrix$target)),]


# Split into training and test matrices
cv_data=patients$cv[match(data_matrix$id,patients$id)]
train_matrix=data_matrix[which(cv_data==1),]
test_matrix1=data_matrix[which(cv_data==2),]
test_matrix2=rbind(test_matrix1,data_matrix[which(cv_data==3),])

saveRDS(train_matrix,file=paste0(dir_output,"train_matrix.RDS"))
saveRDS(test_matrix1,file=paste0(dir_output,"test_matrix1.RDS"))
saveRDS(test_matrix2,file=paste0(dir_output,"test_matrix2.RDS"))

}

######################################################
## Fit models and predict                           ##
######################################################

## Note - on the NSH, this needs to be done carefully between the linux and windows machines.

# Read in data
train_matrix=readRDS(paste0(dir_output,"train_matrix.RDS"))
test_matrix1=readRDS(paste0(dir_output,"test_matrix1.RDS"))
test_matrix2=readRDS(paste0(dir_output,"test_matrix2.RDS"))

model_id=paste0(dir_model,"test_model")
output1=paste0(model_id,".output1")
output2=paste0(model_id,".output2")

lx=list.files(dir_model)
if ((length(lx)>0 & location=="Linux") |(length(lx)>15 & location=="Windows"))
  warning(paste0("Model directory ",dir_model," not empty. This can mean previously-fitted models are being used. Empty directory to ensure reproducibility"))
  

mod=SPARRAv4.fit.split(train_matrix %>% select(-c("id","time")),
  seed=seed,model_id=model_id,linux=(location=="Linux")) # Ignore warnings: there is no signal, so h2o has trouble.
pred1=predict.sparra(test_matrix1 %>% select(-c("id","time")),
  model_id,save_id=output1,verbose=T,linux=(location=="Linux"))
pred2=predict.sparra(test_matrix2 %>% select(-c("id","time")),
  model_id,save_id=output2,verbose=T,linux=(location=="Linux"))



######################################################
## Compare predictions                              ##
######################################################

pred1=readRDS(paste0(output1,".RDS"))
pred2=readRDS(paste0(output2,".RDS"))

p1=pred1$final.full
p2=pred2$final.full

# Check independence to training set
n1=length(p1)
n2=length(p2)

# The first n1 samples in test_matrix1 and test_matrix2 are the same. 
#  Identity of predictions for these samples indicates independence
#  of predictions from the the training set.
all(p1[1:n1]==p2[1:n1])

# Check reproducibility
set.seed(seed)
pred_check=p1[sample(n1,5)]
pred_correct=c(0.06233488, 0.06793719, 0.06407466, 0.04654139, 0.04739927)
if (all(abs(pred_check-pred_correct)<1e-4)) message("Output successfully reproduced") else stop("Output of prediction failed to reproduce")
