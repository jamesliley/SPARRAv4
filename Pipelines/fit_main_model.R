######################################################
## Code to generate constituent models for super-   ##
##  learner, given a fully constructed data matrix, ##
##  and run all principal analyses used for paper.  ##
######################################################
##
## James Liley, 2020
##
##
##
## TEMP FILE DIRECTORY: may need to clear before starting
## /var/lib/docker/overlay2/ef44b903500484de40f907224c970661dd5a5e1358e4a2ad6eeae5662ca42608/diff/tmp/


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
  library(glue)
	library(cvAUC)
  library(ranger)
  library(varhandle)
  library(topicmodels)
} else {
	library(h2o)
	library(tidyverse)
	library(fst)
}

# Directories
dir_cleanData = "../Data/Data_clean/"
output_dir="James/Analysis/Data/full_model/"



######################################################
## Flags and general settings                       ##
######################################################

force_redo=FALSE # Set to TRUE to recalculate everything

time_cutoffs=c( # Time cutoffs for training matrix
  #  dmy_hm("1-5-2014 00:00"), # Earliest time in dataset is "2013-05-01 UTC" (except SPARRALTC!)
  #dmy_hm("1-12-2014 00:00"),
  #dmy_hm("1-5-2015 00:00"),
  #dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"), 
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")) # Latest date in dataset is 2018-05-01

hardmax=round(365.25*3) # Hard maximum lookback for training matrix
hardmax_ltc=round(365.25*10) # We have longer-term data for LTCs and can use a ten-year lookback



######################################################
## Clean raw data                                   ##
######################################################


if (!(location=="Linux")) clean_rawData(force_redo=force_redo)



######################################################
## Load raw data                                    ##
######################################################

fullmatrix=paste0(output_dir,"full_data_matrix.RDS")
if (!file.exists(fullmatrix) | force_redo) {
	load_cleanData(partition = "all",
		subset_frac=1, subset_seed = 1234,
		load_to_global=TRUE,
		load_sparrav3_scores = FALSE)
}


######################################################
## Fit appropriate topic models                     ##
######################################################

doc_term_loc=paste0(output_dir,"doc_term_sparsematrix")
topicfile=paste0(output_dir,"topic_model_fit")
if (location=="Windows" & (!file.exists(paste0(topicfile,"_fold_12.rds")) | force_redo)) {

  # It would take around 10 days to fit a topic model to the full topic training matrix,
  #  which is impractically long. We use only the final time cutoff to train. 
  time_cutoffs_topic_train=time_cutoffs[length(time_cutoffs)]
  time_min_topic_train=time_cutoffs[length(time_cutoffs)]- (hardmax*24*3600)
	topic_train=topic_training_matrix(patients,episodes,list_of_data_tables,
	  time_cutoff=time_cutoffs_topic_train,time_min=time_min_topic_train,
		save_file=doc_term_loc)
	
	print("Generated topic training matrix")
	
	# Restrict to individuals with valid v3 score and nonzero rows
	sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
	tx=as.numeric(sparraV3Scores[[2]])
	sw=which(tx %in% as.numeric(time_cutoffs_topic_train))
	xv3=sparraV3Scores[sw,]; rm(list=c("sparraV3Scores","tx","sw"))
	xv3=xv3[which(!is.na(xv3[[3]])),]
	id=xv3[[1]]; ix=as.numeric(xv3[[2]])
	idt=paste0(id,"_",ix) # slow
	
	topic_train=topic_train[which(rownames(topic_train) %in% idt),]
	topic_train=topic_train[which(rowSums(topic_train)>0),]
	
	print("Trimmed topic training matrix")
	
	# Get CV folds for samples
	topic_id=unlist(lapply(strsplit(rownames(topic_train),"_"),function(x) x[1]))
	cv_sparsematrix=patients$cv[match(topic_id,patients$id)]
	nfold=max(cv_sparsematrix)
	
	# About 1 day (24h) to run per fold
	for (i in 1:nfold) {
	  set.seed(i)
		submat=topic_train[which(cv_sparsematrix!=i),]
		topic_model_fit <- LDA(x = submat, k = 30)
		suffix=paste((1:nfold)[-i],collapse="")
		saveRDS(topic_model_fit,file=paste0(topicfile,"_fold_",suffix,".rds"))
		rm(submat,topic_model_fit)  
		print(suffix)
		print(paste0("Fitted topics for fold ",i))
	}

	rm(list=c("topic_train","cv_sparsematrix","id","idt","ix","xv3","topic_id"))
	gc()
} 




######################################################
## Transformation to various features               ##
######################################################

if (!file.exists(fullmatrix)) {
  ## Subdivisions of patient dataset in order to manage memory
  set.seed(221)
  nsplit=10
  patients_subgroups=rep(1:nsplit,nsplit + nrow(patients)/nsplit)[1:nrow(patients)][order(runif(nrow(patients)))]
  
	## Topic model - derived features
	topic_matrix_files=list.files(output_dir, pattern="topics_final_matrix_complete_training_*",full.names=T)
	topic_model_files=list.files(output_dir,pattern="topic_model_fit_fold_*",full.names=T)
	nfold=length(topic_model_files)
	print("Starting")
	if (location=="Windows" & ((length(topic_matrix_files)<nfold) | force_redo)) {
		
	  if (exists("patients")) rm(list=c("patients","episodes","list_of_data_tables"))
	  for (pid in 1:nsplit) {

	    # load all data
	    load_cleanData(partition = "all",
	      subset_frac=1, subset_seed = 1234,
	      load_to_global=TRUE,
	      load_sparrav3_scores = FALSE)
	    print(paste0("Loaded data, split ",pid))
	    
	    # Subsets of large data tables in order to manage memory   
	    id_sub=patients$id[which(patients_subgroups==pid)]
	    patients_sub=patients[which(patients$id %in% id_sub),]
	    episodes_sub=episodes[which(episodes$id %in% id_sub),]
	    list_of_data_tables_sub=list_of_data_tables
	    for (i in 1:length(list_of_data_tables)) 
	      list_of_data_tables_sub[[i]]=list_of_data_tables[[i]][which(list_of_data_tables[[i]]$id %in% id_sub),]
	  
	    # remove large objects; reload later
	    rm(list=c("patients","episodes","list_of_data_tables"))
	    gc()
	    
	    # Fit models
	    for (i in 1:nfold) {
	      topic_model_fit <- readRDS(file=topic_model_files[i])
	      topic_model_suffix=strsplit(rev(strsplit(topic_model_files[i],"_")[[1]])[1],"[.]")[[1]][1]
	      topic_features <- transformer(patients_sub, episodes_sub, list_of_data_tables_sub,
	        list(all_codes_topic = list(topic_model_fit = topic_model_fit)),
	        as.integer(time_cutoffs),
	        hard_max_lookback = hardmax)
	      print(paste0("Generated topic features for split ",pid,", fold ",i))
	      rm(topic_model_fit)
	      gc()
	      # Fix names: they would otherwise notbe fold-specific
	      txi=colnames(topic_features)
	      gt=grep("topic_",txi)
	      txi[gt]=paste0(txi[gt],"_fit",topic_model_suffix)
	      colnames(topic_features)=txi
	      # Save the resulting final_matrix_topics
	      write.fst(topic_features, paste0(output_dir,"topics_sub/topics_final_matrix_complete_training",topic_model_suffix,"_part",pid,".fst"))
	      # Clean up
	      rm(list=c("topic_features"))
	      gc()
	    }
	  }  

  	# Reassemble topic matrices
  	topic_matrix_files=list.files(paste0(output_dir,"/topics_sub/"), pattern="topics_final_matrix_complete_training_*",full.names=T)
  	suffixes=gsub("training","",unlist(lapply(strsplit(basename(topic_matrix_files),"_"),function(x) x[5])))
  	indices=gsub("part","",unlist(lapply(strsplit(basename(topic_matrix_files),"_"),function(x) x[6])))
  	indices=gsub(".fst","",indices)
	  sx=unique(suffixes)
	  for (i in 1:length(sx)) {
	    tx=c()
	    for (j in 1:length(unique(indices))) {
	      w=which(suffixes==sx[i] & indices==j)
	      mat=read.fst(topic_matrix_files[w])
	      tx=rbind(tx,mat)
	    }
	    write.fst(tx, paste0(output_dir,"topics_final_matrix_complete_training",sx[i],".fst"))
	  }
	}	
	

	## Reload data
	if (!exists("patients")) {
	load_cleanData(partition = "all",
	  subset_frac=1, subset_seed = 1234,
	  load_to_global=TRUE,
	  load_sparrav3_scores = FALSE)
	}
	
	## Features used in SPARRAv3
	v3file=paste0(output_dir,"v3_matrix.RDS")
	if (location=="Windows" & (!file.exists(v3file) | force_redo)) {
		v3_features <- transformer(patients, episodes, list_of_data_tables,
			list(sparrav3 = list()),
			as.integer(time_cutoffs),hard_max_lookback = hardmax)
		saveRDS(v3_features,file=v3file)
	} 
	
	## Long-term condition related features
	ltcfile=paste0(output_dir,"ltc_matrix.RDS")
	if (location=="Windows" & (!file.exists(ltcfile) | force_redo)) {
		ltc_features <- transformer(patients, episodes, list_of_data_tables,
			list(ltcs = list(output_type = "rawdata_NUMBEROFLTCs")),
			as.integer(time_cutoffs),hard_max_lookback = hardmax_ltc)
		saveRDS(ltc_features,file=ltcfile)
	} 
	
	
	## Time-to-event related features
	timefile=paste0(output_dir,"time_matrix.RDS")
	if (location=="Windows" & (!file.exists(timefile) | force_redo)) {
		time_features <- transformer(patients, episodes, list_of_data_tables,
			list(
				last_episode_days_ago = list(source_table_names = c("AE2", "SMR00", "SMR01", "SMR01E", "SMR04")),
				last_emergency_admission_days_ago = list(),
				ltcs = list(output_type = "years_since_diag")
			),
			as.integer(time_cutoffs),hard_max_lookback = hardmax)
		saveRDS(time_features,file=timefile)
	} 

		
	######################################################
	## Form final matrix                                ##
	######################################################
	
	## We want topic model predictions from all three folds, so we adjoin topic matrices
	final_topic_matrix_files=list.files(output_dir, pattern="topics_final_matrix_complete_training[0-9]{1,}.fst",full.names=T)
  topics_all=read.fst(final_topic_matrix_files[1])
  if ("topics_missing_code_total" %in% colnames(topics_all)) topics_all=topics_all %>% select(-c("topics_missing_code_total","topics_missing_code_unique"))
	for (i in 2:length(final_topic_matrix_files)) {
	  subtopic=read.fst(final_topic_matrix_files[i])
	  if ("topics_missing_code_total" %in% colnames(subtopic)) subtopic=subtopic %>% select(-c("topics_missing_code_total","topics_missing_code_unique"))
	  topics_all= topics_all %>% 
			left_join(subtopic,
				by=c("id","time_cutoff"))
	}
  if (exists("v3_features")) rm(list=c("v3_features","time_features","ltc_features"))
  rm(list=c("subtopic"))
  gc()
  
  
	## Final matrix
	nhs = combine_training_matrices(v3file,topics_all,timefile,ltcfile,patients,episodes,list_of_data_tables,keep_id_and_time=TRUE) 
	nhs=nhs[which(nhs$id %in% patients$id),]

	
	## Remove extraneous variables relating to time cutoff
	nhs = nhs %>% select(-c("topic_cutoff","ltc_cutoff"))
	
	
		
#########################################################
## Add sex and CV and do post-processing               ##
#########################################################
	
	# Add sex as predictor
	nhs$sexM=as.integer((patients$gender=="Male")[match(nhs$id,patients$id)])
	
	# CV folds are preset at the time of generating the 'patients' table. This allows splitting of the topic model fit by CV.
	nhs$cv= patients$cv[match(nhs$id,patients$id)]
	
  dsince=colnames(nhs)[grep("days_since",colnames(nhs))]
  for (i in 1:length(dsince)) {
    x=nhs[[dsince[i]]]
    w=which(!is.finite(x)|(x==0))
    x[w]=2*hardmax # If the person has never had an episode, set value to 2* hardmax
    nhs[[dsince[i]]]=x
  }
  ysince=colnames(nhs)[grep("yearssince",colnames(nhs))]
  for (i in 1:length(ysince)) {
    x=nhs[[ysince[i]]]
    w=which(!is.finite(x)|(x==0))
    x[w]=100 # If the person has never had an episode, set value to 100
    nhs[[ysince[i]]]=x
  }
  
	
	
######################################################
## Adjoin v3 scores                                 ##
######################################################
	
# Individual must have a valid v3 score
sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
tx=as.numeric(sparraV3Scores[[2]])
w3=which(tx %in% unique(as.numeric(nhs$time)))
id3=paste0(sparraV3Scores$id[w3],"_",tx[w3])
id4=paste0(nhs$id,"_",as.numeric(nhs$time))

nhs$v3score=sparraV3Scores$SPARRAv3_score[w3][match(id4,id3)]



######################################################
## Save or load                                     ##
######################################################

saveRDS(nhs,file=fullmatrix)
	
} else nhs=readRDS(fullmatrix)



######################################################
## Cross-validation related variables               ##
######################################################

# Column names of training matrix for each fold 
sub12=grep("_fit12",colnames(nhs))
sub13=grep("_fit13",colnames(nhs))
sub23=grep("_fit23",colnames(nhs))

gensub=setdiff(1:ncol(nhs),c(sub12,sub13,sub23))



######################################################
## Exclusions                                       ##
######################################################

exclusions_file=paste0(output_dir,"exclusions.RData")

if (!file.exists(exclusions_file)) {

# Individual must be alive at time cutoff
sub1=which(is.finite(nhs$target))

# Individual should have a valid v3 score
sub2=which(is.finite(nhs$v3score))

v34=intersect(sub1,sub2) # Inclusions for training set
v3x=setdiff(sub1,sub2) # samples with NA v3 scores for testing

save(v34,v3x,file=exclusions_file)

} else load(exclusions_file)



######################################################
## Predict on fold 1                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training and testing sets
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn1=nhs[c(cv[[2]],cv[[3]]),sort(c(gensub,sub23))] %>% select(-c("id","time","cv"))
nhs_tst1=nhs[cv[[1]],sort(c(gensub,sub23))]  %>% select(-c("id","time","cv")); 

# Clean up
rm(list=c("nhs","v34"))
gc()

loc1=paste0(output_dir,"models_cv/model_fold1")
output1=paste0(loc1,".output")
#if (!file.exists(paste0(output1,".RDS"))) {
	m1=SPARRAv4.fit.split(nhs_trn1,train=1:length(cv[[2]]),seed=220,model_id=loc1, 
		linux=(location=="Linux"))
	p1=predict.sparra(nhs_tst1,loc1,save_id=output1,verbose=T,
		linux=(location=="Linux"))
#}
rm(list=c("nhs_trn1","nhs_tst1"))
gc()



######################################################
## Predict on fold 2                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training set
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn2=nhs[c(cv[[3]],cv[[1]]),sort(c(gensub,sub13))]  %>% select(-c("id","time","cv"));
nhs_tst2=nhs[cv[[2]],sort(c(gensub,sub13))]  %>% select(-c("id","time","cv"));

# Clean up
rm(list=c("nhs","v34"))
gc()

loc2=paste0(output_dir,"models_cv/model_fold2")
output2=paste0(loc2,".output")
#if (!file.exists(paste0(output2,".RDS"))) {
	m2=SPARRAv4.fit.split(nhs_trn2,train=1:length(cv[[3]]),seed=220,model_id=loc2, 
		linux=(location=="Linux"))
	p2=predict.sparra(nhs_tst2,loc2,save_id=output2,verbose=T,
		linux=(location=="Linux"))
#}

rm(list=c("nhs_trn2","nhs_tst2"))
gc()



######################################################
## Predict on fold 3                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training and testing sets
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn3=nhs[c(cv[[1]],cv[[2]]),sort(c(gensub,sub12))]  %>% select(-c("id","time","cv"));
nhs_tst3=nhs[cv[[3]],sort(c(gensub,sub12))]  %>% select(-c("id","time","cv"));

# Clean up
rm(list=c("nhs","v34"))
gc()

loc3=paste0(output_dir,"models_cv/model_fold3")
output3=paste0(loc3,".output")
#if (!file.exists(paste0(output3,".RDS"))) {
	m3=SPARRAv4.fit.split(nhs_trn3,train=1:length(cv[[1]]),seed=220,model_id=loc3, 
		linux=(location=="Linux"))
	p3=predict.sparra(nhs_tst3,loc3,save_id=output3,verbose=T,
		linux=(location=="Linux"))
#}

rm(list=c("nhs_trn3","nhs_tst3"))
gc()
