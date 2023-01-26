######################################################
## Code to generate constituent models for super-   ##
##  learner, given a fully constructed data matrix, ##
##  and run all principal analyses used for paper.  ##
######################################################
##
## James Liley, 2020
##
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
dir_rawData = "../../Linked Data/"
dir_cleanData = "../Data/Data_clean/"
output_dir="James/Analysis/Data/full_model/"



######################################################
## Flags and general settings                       ##
######################################################

force_redo=FALSE # Set to TRUE to recalculate everything

time_cutoffs=c( # Time cutoffs for training matrix
  dmy_hm("1-5-2016 00:00"), 
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")) # Latest date in dataset is 2018-05-01

hardmax=round(365.25*3) # Hard maximum lookback for training matrix
hardmax_ltc=round(365.25*5000) # We have longer-term data for LTCs and can use no lookback.


######################################################
## Clean raw data                                   ##
######################################################


if (!(location=="Linux")) clean_rawData(force_redo=force_redo)
gc()


######################################################
## Load raw data                                    ##
######################################################

fullmatrix=paste0(output_dir,"full_data_matrix.RDS")
if (!file.exists(fullmatrix)) {
  load_cleanData(partition = "all",
               subset_frac=1, subset_seed = 1234,
               load_to_global=TRUE)
}


######################################################
## Sample exclusions for topic models               ##
######################################################

topic_exclusions_file=paste0(output_dir,"topic_exclusions.RData")

if (!file.exists(topic_exclusions_file)) {
  
  # Individuals with no SIMD
  exclude_simd=patients$id[which(is.na(patients$SIMD_DECILE_2016_SCT))]
  
  # Individuals who died before the time cutoff
  exclude_death=list()
  for (i in 1:length(time_cutoffs)) {
    exclude_death=patients$id[which(patients$date_of_death < time_cutoffs[i])]
  }
  
  # Individuals on exclusion list
  exclude_list=read.csv(paste0(dir_rawData,"Unmatched_UPIs_without_UPI.csv.gz"))$UNIQUE_STUDY_ID
  
  save(exclude_simd,exclude_death,exclude_list,file=topic_exclusions_file)
  
} else load(topic_exclusions_file)



######################################################
## Fit appropriate topic models                     ##
######################################################

# Save locations
doc_term_loc=paste0(output_dir,"doc_term_sparsematrix")
topicfile=paste0(output_dir,"topic_model_fit")

# Exclusions for topic model (final time cutoff)
topic_exclusions = c(exclude_simd,exclude_death[length(time_cutoffs)],exclude_list)

# Generate training matrices and fit models
if (location=="Windows" & (!file.exists(paste0(topicfile,"_fold_12.rds")) | force_redo)) {

  # It would take around 10 days to fit a topic model to the full topic training matrix,
  #  which is impractically long. We use only the final time cutoff to train. 
  time_cutoffs_topic_train=time_cutoffs[length(time_cutoffs)]
  time_min_topic_train=time_cutoffs[length(time_cutoffs)]- (hardmax*24*3600)
	topic_train=topic_training_matrix(patients,episodes,list_of_data_tables,
	  time_cutoff=time_cutoffs_topic_train,time_min=time_min_topic_train,
		save_file=doc_term_loc)
	
	if(TRUE) { # Check against saved
	  load("../Code/James/Analysis/Testing/check583.RData")
	  ixx=intersect(rownames(sub_topic),rownames(topic_train))
	  tt2=topic_train[ixx,]
	  ix2=c()
	  for (i in 1:length(ixx)) {
	     s1=sort(colnames(sub_topic)[which(sub_topic[ixx[i],]>0)])
	     s2=sort(colnames(tt2)[which(tt2[ixx[i],]>0)])
	     if (!all(s1==s2)) ix2=c(ix2,i)
	  }
	  if (length(ix2)==0) print("Check of topic training matrix successful")
	 }
	print("Generated topic training matrix")
	
	# Restrict to individuals with valid v3 score and nonzero rows
	sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
	tx=as.numeric(sparraV3Scores[[2]])
	sw=which(tx %in% as.numeric(time_cutoffs_topic_train))
	xv3=sparraV3Scores[sw,]; rm(list=c("sparraV3Scores","tx","sw"))
	xv3=xv3[which(!is.na(xv3[[3]])),]
	id=xv3[[1]]; ix=as.numeric(xv3[[2]])
	idt=paste0(id,"_",ix) # slow
	
	# Exclusions for non-valid v3 score and 0 rows
	topic_train=topic_train[which(rownames(topic_train) %in% idt),]
	topic_train=topic_train[which(rowSums(topic_train)>0),]
	
	# Other exclusions
	topic_id_preremoval = unlist(lapply(strsplit(rownames(topic_train),"_"),function(x) x[1]))
	topic_train=topic_train[which(!(topic_id_preremoval %in% topic_exclusions)),]
	
	print("Trimmed topic training matrix")
	print(dim(topic_train))
	
	# Get CV folds for samples
	topic_id=unlist(lapply(strsplit(rownames(topic_train),"_"),function(x) x[1]))
	cv_sparsematrix=patients$cv[match(topic_id,patients$id)]
	nfold=max(cv_sparsematrix)
	
	
	# About 1 day (24h) to run per fold
	for (i in 1:nfold) {
	  set.seed(i)
		submat=topic_train[which(cv_sparsematrix!=i),]
		topic_model_fit <- LDA(x = submat, k = 30, control=list(seed=i+225,verbose=1))
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
	      load_to_global=TRUE)
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
	    
	    # Topic features
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
	  load_to_global=TRUE)
	}
	
	
	### Approximate SPARRA v3 features
	for (t in 1:length(time_cutoffs)) {
	 v3temp=paste0(output_dir,"v3_temp",t,".RDS")
	 v3_features0=transformer_v3like(patients,episodes,list_of_data_tables,time_cutoff=time_cutoffs[t])
	 v3_features0$time_cutoff=as.numeric(time_cutoffs[t])
	 saveRDS(v3_features0,file=v3temp)
	 rm(list=c("v3_features0")); gc()
	 print(paste0("Completed v3 generation for time ", t, " of ",length(time_cutoffs)))
	}
	v3_features=c()
	for (t in 1:length(time_cutoffs)) {
	  v3temp=paste0(output_dir,"v3_temp",t,".RDS")
	  v3_features0=readRDS(v3temp)
	  v3_features=rbind(v3_features,v3_features0)
	  rm(list=c("v3_features0")); gc()
	}
	# Save v3 features
	v3file=paste0(output_dir,"v3_matrix.RDS")
	saveRDS(v3_features,file=v3file)
	rm(list=c("v3_features")); gc()
	
	### Days since/years since type features; rename some
	time_features=c();
	for (t in 1:length(time_cutoffs)) {
	  time_features0= transformer(patients, episodes, list_of_data_tables,
	                              list(
	                                last_episode_days_ago = list(
	                                  source_table_names = c("AE2", "SMR00", "SMR01M", "SMR04")),
	                                last_emergency_admission_days_ago = list(
	                                  source_table_names = c("SMR01M")
	                                ),
	                                last_elective_admission_days_ago = list(
	                                  source_table_names = c("SMR01M")
	                                )
	                              ),
	                              as.integer(time_cutoffs[t]),hard_max_lookback = hardmax) %>%
	    mutate(days_since_last_acute=days_since_last_SMR01M,
	           days_since_last_emergency_admission=days_since_last_SMR01M_emergency_only,
	           days_since_last_elective_admission=days_since_last_SMR01M_elective_only) %>%
	    select(-c(days_since_last_SMR01M,
	              days_since_last_SMR01M_emergency_only,
	              days_since_last_SMR01M_elective_only))
	  time_features=rbind(time_features,time_features0)
	}
	rm(list=c("time_features0")); gc()
	# Save time features
	timefile=paste0(output_dir,"time_matrix.RDS")
	saveRDS(time_features,file=timefile)
	
	
	
	### Long term condition related features. Fill 40 years if no diagnosis.
	ltc_features = c()
	for (t in 1:length(time_cutoffs)) {
	  ltc_features0 = transformer(patients, episodes, list_of_data_tables,
	                           list(ltcs = list(output_type = "total_count"),
	                                ltcs = list(output_type = "years_since_diag",
	                                            time_since_fill_value=40)),
	                           as.integer(time_cutoffs[t]),hard_max_lookback = hardmax_ltc)
	  ltc_features=rbind(ltc_features,ltc_features0)
	}
	rm(list=c("ltc_features0")); gc()
	# Save ltc features
	ltcfile=paste0(output_dir,"ltc_matrix.RDS")
	saveRDS(ltc_features,file=ltcfile)
	

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
  v3file=paste0(output_dir,"v3_matrix.RDS")
  timefile=paste0(output_dir,"time_matrix.RDS")
  ltcfile=paste0(output_dir,"ltc_matrix.RDS")
  nhs = combine_training_matrices(v3file,topics_all,timefile,ltcfile,patients,episodes,list_of_data_tables,keep_id_and_time=TRUE) 
	nhs=nhs[which(nhs$id %in% patients$id),]

	
	## Remove extraneous variables relating to time cutoff
	if ("topic_cutoff" %in% colnames(nhs)) nhs = nhs %>% select(-c("topic_cutoff"))
	if ("ltc_cutoff" %in% colnames(nhs)) nhs = nhs %>% select(-c("ltc_cutoff"))
	
	## Shorten some very long variable names
	nhs = nhs %>% rename(
	  # num_outpt_fu_gen=num_outpatient_appointment_followup_general,
	  # num_outpt_fu_psych=num_outpatient_appointment_followup_psych,  
	  # numLTCs_emerg_admin=numLTCs_resulting_in_emergency_admin,
	  # pis_antiinf_steroids=pis_antiinflammatory_corticosteroids,
	  # numLTCs_elec_admin=numLTCs_resulting_in_elective_admin,
	  # num_outpt_gen=num_outpatient_appointment_general,
	  # num_alc_sub_admin=num_alcohol_substance_admissions,
	  # numLTCs_other_admin=numLTCs_resulting_in_other_admin,  
	  # num_outpt_psych=num_outpatient_appointment_psych,
	  # pis_resp_steroids=pis_respiratory_corticosteroids, 
	  # pis_hyp_heart_fail=pis_hypertensive_heart_failure,
	  # num_alc_drug_ae=num_alcohol_drug_attendances, 
	  # pis_sexhorm_antag=pis_sex_hormone_antagonists,
	  # numLTCs_admin=numLTCs_resulting_in_admin,
	  ltc_FIRST_DIS_BLOOD_EPISODE_yearssincediag=ltc_FIRST_DISEASES_OF_THE_BLOOD_AND_BLOOD_FORMING_ORGANS_EPISODE_yearssincediag,
	  ltc_FIRST_OTHER_DIGESTIVE_EPISODE_yearssincediag=ltc_FIRST_OTHER_DISEASES_OF_DIGESTIVE_SYSTEM_EPISODE_yearssincediag,
	  ltc_FIRST_ENDOCRINE_MET_EPISODE_yearssincediag=ltc_FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE_yearssincediag
	)

#########################################################
## Add CV and do post-processing                       ##
#########################################################
	
	# CV folds are preset at the time of generating the 'patients' table. This allows splitting of the topic model fit by CV.
	nhs$cv= patients$cv[match(nhs$id,patients$id)]
	
	# Code 'days since...' type variables as high rather than low if the individual had no such episode
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
    w=which(!is.finite(x)) # Years can realistically be 0
    x[w]=rep(40,length(w)) # If the person has never had an episode, set value to 40
    nhs[[ysince[i]]]=x
  }

######################################################
## Adjoin v3 scores                                 ##
######################################################


	## Former - add 'reported' v3 scores. These use somewhat different EHRs to our records, so we use 'reconstructed' scores instead
  # Individual must have a valid v3 score
	sparraV3Scores=as_tibble(read.fst("../Results/Modeling/SPARRAv3_raw_scores/SPARRAv3_raw_scores_all.fst"))
	tx=as.numeric(sparraV3Scores[[2]])
	w3=which(tx %in% unique(as.numeric(nhs$time)))
	id3=paste0(sparraV3Scores$id[w3],"_",tx[w3])
	id4=paste0(nhs$id,"_",as.numeric(nhs$time))

	# Affix v3 scores
  nhs$v3score=sparraV3Scores$SPARRAv3_score[w3[match(id4,id3)]]



######################################################
## Save or load                                     ##
######################################################

saveRDS(nhs,file=fullmatrix)
	



######################################################
## Check against saved                              ##
######################################################


load("../Code/James/Analysis/Testing/check583.RData")

# Check final cutoff
sub3=nhs[which(as.numeric(nhs$time)==as.numeric(dmy_hm("1-5-2017 00:00"))),]

# Check v3-like features
xid=intersect(v3_features$id,sub3$id)
testv3=sub3[match(xid,sub3$id),colnames(v3_features)]
checkv3=v3_features[match(xid,v3_features$id),colnames(v3_features)]

itrack=c()
for (i in 1:dim(testv3)[2]) {
  if (!("character" %in% class(checkv3[[i]]))) {
   x1=as.numeric(checkv3[[i]]);
   x2=as.numeric(testv3[[i]])
   x1[which(is.na(x1))]=-1
   x2[which(is.na(x2))]=-1
  } else {
    x1=checkv3[[i]]
    x2=testv3[[i]]
    x1[which(is.na(x1))]="XXX"
    x2[which(is.na(x2))]="XXX"
  }
  if (!all(x1==x2)) {
    print(i)
    itrack=c(itrack,i)
  }
}
if (length(itrack)==0) print("Reproduced v3-like features from online matrix successfully")



# Check time-like features
xid=intersect(time_features$id,sub3$id)
time2=time_features %>% rename(time=time_cutoff)
testtime=sub3[match(xid,sub3$id),colnames(time2)]
checktime=time2[match(xid,time2$id),colnames(time2)]

itrack=c()
for (i in 1:dim(testtime)[2]) {
  if (!("character" %in% class(checktime[[i]]))) {
    x1=as.numeric(checktime[[i]]);
    x2=as.numeric(testtime[[i]])
    x1[which(is.na(x1))]=-1
    x2[which(is.na(x2))]=-1
    x1[which(x1==0)]=2192 # Post-processing recasts 0's as 2192s
  } else {
    x1=checktime[[i]]
    x2=testtime[[i]]
    x1[which(is.na(x1))]="XXX"
    x2[which(is.na(x2))]="XXX"
    x1[which(x1==0)]=2192 # Post-processing recasts 0's as 2192s
  }
  if (!all(x1==x2)) {
    print(i)
    itrack=c(itrack,i)
  }
}
if (length(itrack)==0) print("Reproduced time-since-event features from online matrix successfully")



# Check ltc-related features
xid=intersect(ltc_features$id,sub3$id)
ltc2=ltc_features %>% rename(time=time_cutoff)

## Shorten some very long variable names
ltc2 = ltc2 %>% rename(
  ltc_FIRST_DIS_BLOOD_EPISODE_yearssincediag=ltc_FIRST_DISEASES_OF_THE_BLOOD_AND_BLOOD_FORMING_ORGANS_EPISODE_yearssincediag,
  ltc_FIRST_OTHER_DIGESTIVE_EPISODE_yearssincediag=ltc_FIRST_OTHER_DISEASES_OF_DIGESTIVE_SYSTEM_EPISODE_yearssincediag,
  ltc_FIRST_ENDOCRINE_MET_EPISODE_yearssincediag=ltc_FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE_yearssincediag
)

testltc=sub3[match(xid,sub3$id),colnames(ltc2)]
checkltc=ltc2[match(xid,ltc2$id),colnames(ltc2)]


itrack=c()
for (i in 1:dim(testltc)[2]) {
  if (!("character" %in% class(checkltc[[i]]))) {
    x1=as.numeric(checkltc[[i]]);
    x2=as.numeric(testltc[[i]])
    x1[which(is.na(x1))]=-1
    x2[which(is.na(x2))]=-1
  } else {
    x1=checkltc[[i]]
    x2=testltc[[i]]
    x1[which(is.na(x1))]="XXX"
    x2[which(is.na(x2))]="XXX"
  }
  if (!all(x1==x2)) {
    print(i)
    itrack=c(itrack,i)
  }
}
if (length(itrack)==0) print("Reproduced long-term condition features from online matrix successfully")

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
## Sample exclusions for design matrix              ##
######################################################

design_exclusions_file=paste0(output_dir,"design_exclusions.RData")

if (!file.exists(design_exclusions_file)) {

  # Load data
  load_cleanData(partition = "all",
                 subset_frac=1, subset_seed = 1234,
                 load_to_global=TRUE)
  
  # Individuals with no SIMD
  exclude_simd=patients$id[which(is.na(patients$SIMD_DECILE_2016_SCT))]
  index_simd=which(!(nhs$id %in% exclude_simd))

  # Individuals on exclusion list
  exclude_list=read.csv(paste0(dir_rawData,"Unmatched_UPIs_without_UPI.csv.gz"))$UNIQUE_STUDY_ID
  index_list=which(!(nhs$id %in% exclude_list))
  
  # Individual must be alive at time cutoff
  index_death=which(is.finite(nhs$target))
  
  # Individual should have a valid v3 score
  index_v3=which(is.finite(nhs$v3score))
  
  # All inclusions
  v34=intersect(intersect(index_list,index_simd),intersect(index_death,index_v3)) # Inclusions for training set
  v3x=setdiff(intersect(intersect(index_list,index_simd),index_death),index_v3) # samples with NA v3 scores, but otherwise OK
  
  save(v34,v3x,file=design_exclusions_file)
  
} else load(design_exclusions_file)



######################################################
## Predict on fold 1                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(design_exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training and testing sets
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn1=nhs[c(cv[[2]],cv[[3]]),sort(c(gensub,sub23))] %>% select(-c("id","time","cv","reason"))
nhs_tst1=nhs[cv[[1]],sort(c(gensub,sub23))]  %>% select(-c("id","time","cv","reason")); 

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
if (location=="Linux") tryCatch(h2o.shutdown(prompt=FALSE),error=function(e) {})
gc()



## Fold 1 - no topics

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(design_exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training and testing sets
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn1_notopic=nhs[c(cv[[2]],cv[[3]]),sort(c(gensub,sub23))] %>% select(-c("id","time","reason",matches("topic*")))
nhs_tst1_notopic=nhs[cv[[1]],sort(c(gensub,sub23))]  %>% select(-c("id","time","reason",matches("topic*")))

# Clean up
rm(list=c("nhs","v34"))
gc()

loc1_notopic=paste0(output_dir,"models_cv/model_fold1_notopic")
output1_notopic=paste0(loc1_notopic,".output")
#if (!file.exists(paste0(output1,".RDS"))) {
	m1=SPARRAv4.fit.split(nhs_trn1_notopic,train=1:length(cv[[2]]),seed=220,model_id=loc1_notopic, 
		linux=(location=="Linux"))
	p1=predict.sparra(nhs_tst1_notopic,loc1_notopic,save_id=output1_notopic,verbose=T,
		linux=(location=="Linux"))
#}
rm(list=c("nhs_trn1_notopic","nhs_tst1_notopic"))
if (location=="Linux") tryCatch(h2o.shutdown(prompt=FALSE),error=function(e) {})
gc()





######################################################
## Predict on fold 2                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(design_exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training set
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn2=nhs[c(cv[[3]],cv[[1]]),sort(c(gensub,sub13))]  %>% select(-c("id","time","cv","reason"));
nhs_tst2=nhs[cv[[2]],sort(c(gensub,sub13))]  %>% select(-c("id","time","cv","reason"));

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
if (location=="Linux") tryCatch(h2o.shutdown(prompt=FALSE),error=function(e) {})
gc()



######################################################
## Predict on fold 3                                ##
######################################################

# Regererate nhs, for memory sakes
nhs=readRDS(fullmatrix)
load(design_exclusions_file)
nhs=nhs[v34,] # redefine nhs for memory considerations 

# Training and testing sets
cv=list(fold1=which(nhs$cv==1),fold2=which(nhs$cv==2),fold3=which(nhs$cv==3))
nhs_trn3=nhs[c(cv[[1]],cv[[2]]),sort(c(gensub,sub12))]  %>% select(-c("id","time","cv","reason"));
nhs_tst3=nhs[cv[[3]],sort(c(gensub,sub12))]  %>% select(-c("id","time","cv","reason"));

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
if (location=="Linux") tryCatch(h2o.shutdown(prompt=FALSE),error=function(e) {})
gc()
