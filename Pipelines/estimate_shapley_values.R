######################################################
## Code to estimate Shapley values for all          ##
##  individuals
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
  library(ranger)
  library(glmnet)
  library(Metrics)
  library(varhandle)
} else {
  library(h2o)
  library(tidyverse)
  library(fst)
}

# Directories
dir_cleanData = "../Data/Data_clean/" # where raw data is 
output_dir="James/Analysis/Data/full_model/" # models saved here
fullmatrix=paste0(output_dir,"full_data_matrix.RDS") # full data matrix
exclusions_file=paste0(output_dir,"exclusions.RData")


####################################################################################
## Read data                                                                      ##
####################################################################################

nhs=readRDS(fullmatrix)
load(exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 


cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))


######################################################
## Cross-validation related variables               ##
######################################################

# Column names of training matrix for each fold 
sub12=grep("_fit12",colnames(nhs))
sub13=grep("_fit13",colnames(nhs))
sub23=grep("_fit23",colnames(nhs))

gensub=setdiff(1:ncol(nhs),c(sub12,sub13,sub23))



######################################################
## Subset of individuals to calculate SVs for       ##
######################################################

nsub=2e4 # Sample this many individuals. We will only use fold 1 for simplicity.
isub=(1:nsub)*floor(length(cv[[1]])/nsub)
xsub=cv[[1]][isub] # Indices of chosen individuals in cv[[1]]: effectively random, but deterministic between different R versions and versions of sample()




######################################################
## Fold 1                                           ##
######################################################

# Test matrix for fold 1
nhs_tst1_id=nhs[intersect(cv[[1]],xsub),sort(c(gensub,sub23))]  %>% select(c("id","time","cv")); 
nhs_tst1s=nhs[intersect(cv[[1]],xsub),sort(c(gensub,sub23))]  %>% select(-c("id","time","cv")); 

# Clean up
rm(list=c("nhs","v34"))
gc()


loc1=paste0(output_dir,"models_cv/model_fold1")
shapley1=get_shap.split(nhs_tst1s,loc1,proc_id="shapley_fold1",seed=220,ni=500,
  verbose=T,linux=(location=="Linux"),
  remove_when_done=FALSE)


######################################################
## Save  results                                    ##
######################################################

if (!(location=="Linux")) {
shapley1$id_matrix=cbind(nhs_tst1_id,nhs_tst1s)
shapley1$ids_in_full=intersect(cv[[1]],xsub)
saveRDS(shapley1,file=paste0(output_dir,"shapley_values.RDS"))
}


