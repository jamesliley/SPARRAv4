#' Returns posterior assignment probabilities of a patient belonging to each of the diagnosis and prescription clusters (=topics)
#' @param topic_model_fit a fitted topicmodels::LDA() object, trained on ICD10 codes representing medical history and prescription summaries
#' @param return_extra_codes binary, if TRUE, also return a column with the number of total and number of unique codes that are unrecognised by the trained topic model
#' 
#' 
transformer_all_codes_topic <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  topic_model_fit
  
){
  # Get a doc_term matrix to infer on
  topic_matrix=topic_training_matrix(patients,episodes,list_of_data_tables,
    time_cutoff=time_cutoff)
  
  # Ensure variable names match
  ctest=colnames(topic_matrix)
  ctrain=colnames(posterior(topic_model_fit)$terms)
  cplus=setdiff(ctrain,ctest)
  xmat=as(matrix(0,dim(topic_matrix)[1],length(cplus)),"dgCMatrix")
  colnames(xmat)=cplus
  topic_matrix=cbind(topic_matrix,xmat)
  topic_matrix=topic_matrix[,ctrain]
  rownames(topic_matrix)=unlist(lapply(strsplit(rownames(topic_matrix),"_"),function(x) x[1]))
  rm(list=c("xmat","ctest","ctrain","cplus"))
  
  print("Diagnosis/prescription matrix generated")
  gc()

  # Do the inference given the fit and the new doc_term matrix
  xpos=which(rowSums(topic_matrix)>0)
  lda_post <- posterior(
    topic_model_fit, 
    newdata = topic_matrix[xpos,])
  
  topicprob=matrix(0,dim(topic_matrix)[1],dim(lda_post$topics)[2])
  topicprob[xpos,]=lda_post$topics
  rownames(topicprob)=rownames(topic_matrix)
  colnames(topicprob)=colnames(lda_post$topics)
  
  rm(topic_matrix)
  rm(lda_post)
  
  print("Posterior calculated")
  gc()
    
  # Create the inference matrix 
  out_matrix <- as_tibble(
      topicprob,
      rownames = "id"
    ) %>% 
    rename_at(.vars = 2:ncol(.), ~ paste0("topic_", .)) %>%
    mutate(
      id = as.integer(id)
    )
  rm(topicprob)
  
  print("Inference matrix generated")
  

  out_missing_default = as.list(rep(0, length=(ncol(out_matrix)-1)))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  print("Post-processing completed")
  
  # Return in the required transformer output
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}




#' transformer_count
#' Returns the count of observations in the given column(s)
#'
#' @param column_names_flat (required) names of columns in format <TABLE_NAME.COLUMN_NAME>
#' @param unique_vals (optional, default = FALSE) 
#' 
#' Usage example
#' transformer_count(
#'       patients, episodes, list_of_data_tables, 
#'       time_cutoff,
#'       column_names_flat = c("AE2.DIAGNOSIS_1", "SMR00.DAYS_WAITING", "SMR01.DAYS_WAITING", "SMR01.time"),
#'       unique_vals = TRUE
#'  )
#' @return The counts of occurencies of the requested variables. If unique_vals=TRUE, only unique occurances are counted.
#' Note that when unique_vals==FALSE, we only really get the number of records from a table, 
#' ie SMR00.DAYS_WAITING and SMR00.time are going to return the exact same numbers (ie the number of rows in SMR00 per patient id)
#' 
#' @export
transformer_count <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  column_names_flat,
  
  # Further optional parameters (these have defaults)
  unique_vals = FALSE
){
  
  # Define the function to apply
  if (unique_vals){
    function_to_apply <- function(x)length(unique(x))
  } else {
    function_to_apply <- function(x)length(x)
  }
  
  # Create the out_matrix
  out_matrix = summarise_per_patient_from_list_of_tables(
    list_of_data_tables,
    column_names_flat,
    summary_func = function_to_apply,
    table_key = "id"
  )
  
  # Set missing defaults
  out_missing_default = as.list(rep(0, length=length(names(out_matrix)[-1])))
  names(out_missing_default) <- names(out_matrix)[-1]

  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}




#' Returns how long ago a patient last had an emergency admission
#'
#' @param source_table_names (required) names of the tables where one wishes to know the last interaction per patient 
#' (could be used to filter out emergency admissions to geriatric, SMR01E, or mental, SMR04, wards)
#' (By default we keep geriateric and do not keep mental emergency admissions)
#'     
#' Importantly this outputs 0 when there was no visit of a particular type, or at least 1 if there was a visit (even if within less than a day)
#' 
transformer_last_emergency_admission_days_ago <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  
  # # Further optional parameters (these have defaults)
  source_table_names = c("SMR01M")
  
){
  
  print("Check (last emergency admission days) 1/3")
  
  # Filter the appropriate tables
  data_inpatient <- list_of_data_tables$SMR01M %>% 
      select(id, time,  emergency_admin ) %>%
      mutate(source_table="SMR01M",
             source_row=1:dim(list_of_data_tables$SMR01M)[1]) %>%
      filter(emergency_admin==1)

  print("Check (last emergency admission days) 2/3")
  
  # Call the helper function "transformer_time_since ..."
  out <- transformer_last_episode_days_ago(
    patients = NULL,
    episodes = episodes %>%
      select(id, time, source_table) %>%
      semi_join(data_inpatient),
    list_of_data_tables = NULL,
    time_cutoff = time_cutoff,
    source_table_names = source_table_names
  )
  
  print("Check (last emergency admission days) 3/3")
  
  # Rename the columns adding "_emergency_only" to indicate filtering for emergency admissions only
  colnames(out$matrix) <- c("id", sapply(colnames(out$matrix)[-1], function(x)paste0(x, "_emergency_only")))
  names(out$missing_default) <- sapply(names(out$missing_default), function(x)paste0(x, "_emergency_only"))
  
  
  # Return
  out
           
}



#' Returns how long ago a patient last had an elective admission
#'
#' @param source_table_names (required) names of the tables where one wishes to know the last interaction per patient 
#' (could be used to filter out emergency admissions to geriatric, SMR01E, or mental, SMR04, wards)
#'     
#' Importantly this outputs 0 when there was no visit of a particular type, or at least 1 if there was a visit (even if within less than a day)
#' 
transformer_last_elective_admission_days_ago <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  
  # # Further optional parameters (these have defaults)
  source_table_names = c("SMR01M")
  
){
  
  print("Check (last elective admission days) 1/3")
  
  # Filter the appropriate tables
  data_inpatient <- list_of_data_tables$SMR01M %>% 
    select(id, time, elective_admin) %>%
    mutate(source_table="SMR01M",
           source_row=1:dim(list_of_data_tables$SMR01M)[1]) %>%
    filter(elective_admin==1)
  
  print("Check (last elective admission days) 2/3")
  
  # Call the helper function "transformer_time_since ..."
  out <- transformer_last_episode_days_ago(
    patients = NULL,
    episodes = episodes %>%
      select(id, time, source_table) %>%
      semi_join(data_inpatient),
    list_of_data_tables = NULL,
    time_cutoff = time_cutoff,
    source_table_names = source_table_names
  )
  
  print("Check (last elective admission days) 3/3")
  
  # Rename the columns adding "_emergency_only" to indicate filtering for emergency admissions only
  colnames(out$matrix) <- c("id", sapply(colnames(out$matrix)[-1], function(x)paste0(x, "_elective_only")))
  names(out$missing_default) <- sapply(names(out$missing_default), function(x)paste0(x, "_elective_only"))
  
  
  # Return
  out
  
}







#' Returns how long ago a patient last had an interaction with NHS of a particular type (sorted by tables)
#' Conditional filtering should be done prior to calling this helper funciton
#'
#' @param source_table_names (required) names of the tables where one wishes to know the last interaction per patient
#'     
#' Importantly this outputs 0 when there was no visit of a particular type, or at least 1 if there was a visit (even if within less than a day)
#' 
transformer_last_episode_days_ago <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  source_table_names
  
  # # Further optional parameters (these have defaults)
){
  
  if (("SMR01M" %in% source_table_names) & !is.null(list_of_data_tables$SMR01)) {
    # Correct SMR01M
    SMR01X=affix_smr01(
      list_of_data_tables$SMR01M,
      list_of_data_tables$SMR01,
      t-three_year_lookback,
      filtered=TRUE
    ); gc()
    if (!is.null(SMR01X)) {
      SMR01X$source_table="SMR01M"
      episodes=episodes %>%
        bind_rows(SMR01X %>%select(colnames(episodes)[colnames(episodes) %in% colnames(.)]))
    }
      
    SMR01XE=affix_smr01e(
      list_of_data_tables$SMR01M,
      list_of_data_tables$SMR01E,
      t-three_year_lookback,
      filtered=TRUE
    ); gc()
    if (!is.null(SMR01XE)) {
      SMR01XE$source_table="SMR01M"
      episodes=episodes %>%
        bind_rows(SMR01XE %>%select(colnames(episodes)[colnames(episodes) %in% colnames(.)]))
    }
    
  }
  
  print("Check (last episode days ago) 1/3")
  # Get the days since last visit for each table in "source_table_names"
  out_matrix <- episodes %>%
    ungroup() %>%
    mutate(time = as.integer(time)) %>% 
    select(id, time, source_table) %>%
    filter(source_table %in% source_table_names) %>% # Get rid of unnecessary information (conditional filtering should be done prior to calling this helper funciton)
    mutate(source_table = paste0("days_since_last_", source_table)) %>%
    mutate(source_table = factor(source_table, levels = paste0("days_since_last_", source_table_names))) %>% # ensure that returns as many columns as entries in source_table_names
    group_by(id, source_table)
  out_matrix=out_matrix %>%
    summarise(last_time = 24*3600*floor(max(time)/(24*3600))) %>% 
    ungroup() %>%
    mutate(last_time = time_cutoff - as.integer(last_time)) %>% # Get time since last visit (in seconds)
    mutate(last_time = last_time %/% as.integer(60*60*24)) %>% # convert to days 
    mutate(last_time = as.integer(ifelse(last_time==0, 1, last_time))) %>% # min 1 day if the entry exists, 0 means no entry
    spread(
       key = source_table,
       value = last_time,
       drop = FALSE
     )
           
  print("Check (last episode days ago) 2/3")
  
  # Set missing defaults
  out_missing_default = as.list(rep(0, length=length(names(out_matrix)[-1])))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  print("Check (last episode days ago) 3/3")
  
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
           
}

#' Returns the latest observation in a given column
#'
#' @param column_names_flat (required) names of columns in format <TABLE_NAME.COLUMN_NAME>
#' 
#' Usage example
#' transformer_latest(
#'       patients, episodes, list_of_data_tables, time_cutoff
#'       column_names_flat = c("AE2.DIAGNOSIS_1", "SMR00.DAYS_WAITING", "SMR01.DAYS_WAITING", "SMR01.time")
#'  )
#' 
transformer_latest <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  column_names_flat
  
  # Further optional parameters (these have defaults)
){

  # Define the function to apply
  function_to_apply <- function(x){
    ifelse(length(x)>0,
           x[length(x)],
           NA)
  }
  
  # Create the out_matrix
  out_matrix = summarise_per_patient_from_list_of_tables(
    list_of_data_tables,
    column_names_flat,
    summary_func = function_to_apply,
    table_key = "id"
  )
  
  out_missing_default = as.list(rep(NA, length=length(column_names_flat)))
  names(out_missing_default) <- column_names_flat

  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}




#' Returns various statistics about the Long Term Conditons stored in the "SPARRALTC" data table
#' @param output_type one of c("total_count", "years_since_diag", "binary", "rawdata_NUMBEROFLTCs")
#' Where "total_count" and "rawdata_NUMBEROFLTCs" refer to the total number of LTCs a person has been diagnoised with (and should be mostly equal)
#' "years_since_diag" returns one column per LTC with an integer representing years since diagnosis
#' "binary" returns one column per LTC, indicating whether or not the person has that LTC
#' @param time_since_fill_value for output_type="years_since_diag", uses this value for individuals not diagnosed with the relevant condition.
#'
#' @export
transformer_ltcs <- function(
  # Input data
  patients, episodes, list_of_data_tables,
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  output_type, # %in% c("total_count", "years_since_diag", "binary", "rawdata_NUMBEROFLTCs")
  
  time_since_fill_value=NA
  # Further optional parameters (these have defaults)
  
){

  print("Check (LTC) 1/5")
    
  if (!output_type %in% c("total_count", "years_since_diag", "binary", "rawdata_NUMBEROFLTCs")){
    stop("Non-implemented \"output_type\" input in transformer_ltc()")
  }
  
  if (output_type=="rawdata_NUMBEROFLTCs"){
    out_matrix <- list_of_data_tables$SPARRALTC %>%
      select(id, NUMBEROFLTCs) %>%
      distinct() %>%
      rename(ltc_rawdata_NUMBEROFLTCS:=NUMBEROFLTCs)
  }
  
  if (output_type=="total_count"){
    out_matrix <- list_of_data_tables$SPARRALTC %>%
      group_by(id) %>%
      summarise(ltc_total_count = n())
  }
  
  print("Check (LTC) 2/5")
  
  ltc_names <-c("FIRST_ARTHRITIS_EPISODE", "FIRST_ASTHMA_EPISODE", "FIRST_ATRIAL_FIBRILLATION_EPISODE", 
                "FIRST_COPD_EPISODE", "FIRST_CANCER_EPISODE", "FIRST_CEREBROVASCULAR_DISEASE_EPISODE", 
                "FIRST_CHRONIC_LIVER_DISEASE_EPISODE", "FIRST_DEMENTIA_EPISODE", 
                "FIRST_DIABETES_EPISODE", "FIRST_EPILEPSY_EPISODE", "FIRST_HEART_DISEASE_EPISODE", 
                "FIRST_HEART_FAILURE_EPISODE", "FIRST_MULTIPLE_SCLEROSIS_EPISODE", 
                "FIRST_PARKINSON_DISEASE_EPISODE", "FIRST_RENAL_FAILURE_EPISODE", 
                "FIRST_CONGENITAL_PROBLEMS_EPISODE", "FIRST_DISEASES_OF_THE_BLOOD_AND_BLOOD_FORMING_ORGANS_EPISODE", 
                "FIRST_OTHER_DISEASES_OF_DIGESTIVE_SYSTEM_EPISODE", "FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE"
  )
  
  if (output_type=="binary"){
    out_matrix <- list_of_data_tables$SPARRALTC %>%
      select(id, LTC_TYPE) %>%
      mutate(LTC_TYPE = factor(LTC_TYPE, levels = ltc_names)) %>%
      mutate(value = as.integer(1)) %>%
      spread(key = LTC_TYPE, value = value, fill=as.integer(0), drop=FALSE)
    
    colnames(out_matrix) = c(c("id"), sapply(colnames(out_matrix)[-1], function(x)stringr::str_c("ltc_", x)))
  }

  print("Check (LTC) 3/5")
  
  
  if (output_type=="years_since_diag"){
    out_matrix <- list_of_data_tables$SPARRALTC %>%
      select(id, LTC_TYPE, time) %>%
      mutate(LTC_TYPE = factor(LTC_TYPE, levels = ltc_names)) %>%
      mutate(value = as.integer(time)/(60*60*24*365.25)) %>% # Divide here to avoid overflow issues
      select(-time) %>%
      mutate(value = sapply(as.integer((time_cutoff/(60*60*24*365.25)) - value), function(x)max(x, as.integer(0)))) %>% # Subtract for time_cutoff (in years)
      spread(key = LTC_TYPE, value = value, fill=time_since_fill_value, drop=FALSE)
    
    colnames(out_matrix) = c(c("id"), sapply(colnames(out_matrix)[-1], function(x)stringr::str_c("ltc_", x, "_yearssincediag")))
  }
  
  
  print("Check (LTC) 4/5")
  
  
  out_missing_default = as.list(rep(0, length=(ncol(out_matrix)-1)))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  print("Check (LTC) 5/5")
  
  # Return in the required transformer output
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}







#' Returns the sum of observations in a given column
#'
#' @param column_names_flat (required) names of columns in format <TABLE_NAME.COLUMN_NAME>
#' 
#' Usage example
#' transformer_sum(
#'       patients, episodes, list_of_data_tables,  time_cutoff
#'       column_names_flat = c("SMR00.DAYS_WAITING", "SMR01.DAYS_WAITING", "PIS.NUMBER_OF_PAID_ITEMS", "PIS.PAID_GIC_INCL_BB")
#'  )
#' 
transformer_sum <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  column_names_flat
  
  # Further optional parameters (these have defaults)
  
){

  function_to_apply <- quo(sum(., na.rm=TRUE))

  # Create the out_matrix
  out_matrix = summarise_per_patient_from_list_of_tables(
    list_of_data_tables,
    column_names_flat,
    summary_func = function_to_apply,
    table_key = "id"
  )
  
  out_missing_default = as.list(rep(0, length=length(column_names_flat)))
  names(out_missing_default) <- column_names_flat

  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}





## New v3 transformer
## Simon Rogers, Oct 21


#' transformer_sparrav3approx_features
#' This transformer approximately replicates the features used by SPARRAv3. An exact replication of 
#' v3 features is given in the function v3_feature_creation().
#' 
#' The changes from the exact v3 specificiations are:
#'  1. A universal 3-year lookback on all tables except LTC, and an arbitrary lookback (30y) there
#'  2. Identical BNF section counts across cohorts
#'  3. Removed binarised versions of bnf_vitamins or bnf_bandages
#'  4. epilepsy_indicated and num_ltc identical across cohorts
#'  5. Removal of constant columns blood_indicated, endocrine_indicated, congenital_indicated and digestive_indicated
#'  6. Added pis_antiepileptics, pis_parkinsonism, pis_respiratory_corticosteroids, pis_minerals
#'  7. Added pis_bandages and pis_vitamins as multi-BNF categories
#'  8. Removed redundant pis_substance and corrected spelling of pis_anticoagulants
#'  9. Added 'indicated' covariates for other LTC-listed covariates, and removed BNF/LTC dual features epilepsy_indicated etc
#'  10. Used a single unified BNF grouping, adding multiple pis_ features.
#'  11. PIS cutoff is now the same as other cutoffs rather than one month earlier
#'  
#' This function was adapted from that written by Simon Rogers, Oct 21.
#' 
#' 
#' Usage example
#' transformer_sparrav3(patients, episodes, list_of_data_tables, time_cutoff)
#' @param patients patients table from load_cleanData()
#' @param episodes episodes table from load_cleanData()
#' @param list_of_data_tables list_of_data_tables object from load_cleanData
#' @param lookback lookback in days
#' @param time_cutoff time cutoff
#' @returns a list with components matrix and missing_default which can be parsed by the function transformer()
transformer_sparrav3approx_features = function(patients,episodes,list_of_data_tables,time_cutoff,
  lookback=round(3*365.25)) {
  # Cutoff data and mutate to unixtime
  episodes = episodes %>% mutate(unixtime = as.numeric(time))
  for (i in 1:length(list_of_data_tables)) {
    list_of_data_tables[[i]] = list_of_data_tables[[i]] %>% mutate(unixtime = as.numeric(time))
  }
  
  # Lookback in seconds
  general_lookback = lookback * 3600 * 24
  ltc_lookback = 100 * 366 * 3600 * 24 # 100 years, ie as long as possible
  lookback_table=data.frame(
    source_table=factor(c(
      "SMR01",
      "SMR00",
      "SMR01E",
      "SystemWatch",
      "SPARRALTC",
      "AE2",
      "SMR04",
      "PIS",
      "deaths"),levels=levels(episodes$source_table)),
    max_lookback = c(
      general_lookback, # These differ from the exact v3 covariates, which used a one-year lookback.
      general_lookback,
      general_lookback,
      general_lookback,
      ltc_lookback,
      general_lookback,
      general_lookback,
      general_lookback,
      ltc_lookback
    )
  ) %>% mutate(max_lookback = as.numeric(max_lookback))
  
  episodes = episodes %>% left_join(lookback_table, by =c("source_table"))
  
  # Filter episodes and data tables
  filtered_episodes = episodes %>%
    filter(
      unixtime < time_cutoff,
      unixtime > time_cutoff - max_lookback
    ) %>%
    select(-unixtime)
  
  filtered_list_of_data_tables = list_of_data_tables
  for (table_name in names(filtered_list_of_data_tables)) {
    max_l = lookback_table %>% filter(source_table == table_name)
    tc = time_cutoff
    print(paste0(table_name," ",max_l$max_lookback))
    filtered_list_of_data_tables[[table_name]] = filtered_list_of_data_tables[[table_name]] %>%
      filter(
        unixtime < tc,
        unixtime >= tc - max_l$max_lookback
      )
  }
  filtered_patients = patients %>% filter(id %in% filtered_episodes$id)
  
  # Output data structure
  out = NULL
  
  
  
  # SMR01 + Systemwatch merge
  # This merges data and adds flags for
  # emergency, alcohol, substance abuse, daycase, elective
  print(paste0("Merging SMR01"))
  SMR01_merged=merge_smr01(filtered_list_of_data_tables,time_cutoff=time_cutoff)
  print(paste0("Completed merge"))
  
  # Emergency / elective admissions and bed days
  print(paste0("Emergency / elective admissions and bed days"))
  tmp = admissions_and_bed_days(SMR01_merged)
  out = join_transformer_outputs(out,tmp)
  
  # Features that indicate a LTC based on either an LTC or a PIS entry
  print(paste0("Indicated features"))
  tmp = indicated_features(filtered_list_of_data_tables)
  out=join_transformer_outputs(out,tmp)
  
  # LTC features
  # - counts of LTCs are different per cohort
  # - this makes a general_ltc_count (for LTC, FE) and a yed_count (for YED): TODO (U16?)
  print(paste0("LTC_features"))
  tmp = ltc_features(filtered_list_of_data_tables$SPARRALTC)
  out = join_transformer_outputs(out, tmp)
  
  # General outpatient attendances
  print(paste0("General outpatient"))
  tmp = general_outpatient(filtered_list_of_data_tables$SMR00)
  out= join_transformer_outputs(out,tmp)
  
  # General outpatient attendances
  print(paste0("Psychiatric outpatient"))
  tmp = psych_outpatient(filtered_list_of_data_tables$SMR00)
  out= join_transformer_outputs(out,tmp)
  
  
  # Psych admissions
  print(paste0("Psych admissions"))
  tmp = psych_admissions(filtered_list_of_data_tables$SMR04)
  out = join_transformer_outputs(out,tmp)
  
  
  # AE2 attendances
  # Note - had to request additional attendance code column from Fraser
  print(paste0("AE2"))
  tmp = ae2_attendances(filtered_list_of_data_tables$AE2)
  out = join_transformer_outputs(out,tmp)
  
  # PIS section counts. These are the individual PIS features.
  # e.g. pis_diabetes
  print(paste0("PIS section counts"))
  for (pis_feature in names(bnf_groupings)) {
    print(paste0(pis_feature))
    tmp = bnf_count(
      filtered_list_of_data_tables$PIS,
      bnf_groupings[[pis_feature]],
      pis_feature
    )
    out=join_transformer_outputs(out,tmp)
  }
  
  
  # Number of BNF sections
  print(paste0("NUM BNF SECTIONS"))
  tmp = bnf_num_sections(filtered_list_of_data_tables$PIS)
  out = join_transformer_outputs(out,tmp)
  
  
  # Age, sex and SIMD (added by JL for NSH)
  print(paste0("Age, sex and SIMD"))
  tmp=transformer_patient_info(patients,episodes,list_of_data_tables,time_cutoff=time_cutoff)
  tmp$matrix$gender=as.integer(tmp$matrix$gender=="Male")
  tmp$matrix =tmp$matrix %>% rename(sexM=gender)
  names(tmp$missing_default)=gsub("gender","sexM",names(tmp$missing_default))
  out = join_transformer_outputs(out,tmp)
  
  
  #	out$matrix[is.na(out$matrix)] = 0
  return(out)
}



###################################################################################
###################################################################################
###                                                                             ###
###                              UTILITY METHODS                                ###
###                                                                             ###
###################################################################################
###################################################################################


##################################################################
##                      Feature joiner                          ##
##################################################################

join_transformer_outputs = function(out,new_inp) {
  if (is.null(out)) {
    out = new_inp
  } else {
    out$matrix = out$matrix %>%
      full_join(new_inp$matrix,by="id")
    out$missing_default = c(out$missing_default,new_inp$missing_default)
  }
  return(out)
}

##################################################################
##                     BNF section count                        ##
##################################################################

bnf_num_sections = function(PIS) {
  tmp = list()
  # Extract relevant columns, binarise number sold, keep just one of each section
  pis_sub = PIS %>%
    select(id, BNF_section,NUM_sold) %>%
    mutate(NUM_sold = ifelse(NUM_sold > 0 , 1 ,0)) %>%
    group_by(id, BNF_section) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  # Add the individual (all, we remove the excluded ones later)
  temp0 = pis_sub %>%
    filter(BNF_section %in% as.vector(unlist(bnf_groupings))) #%>%
  #	group_by(id) #%>%
  #summarise(all_sections = n_distinct(BNF_section)) # This is very slow on the NSH: rewritten below
  count_data=as_tibble(data.frame(
    id=unique(temp0$id),
    all_sections=as.integer(
      table(temp0$id[which(!duplicated(
        paste0(temp0$id,as.character(temp0$BNF_section))))]
      ))))
  
  
  # Do the total - i.e. all sections in the table
  temp0 = pis_sub #%>%
  #	group_by(id) #%>%
  #	summarise(num_bnf_total = n_distinct(BNF_section)) # This is very slow on the NSH: rewritten below
  temp=as_tibble(data.frame(
    id=unique(temp0$id),
    num_bnf_total=as.integer(
      table(temp0$id[which(!duplicated(
        paste0(temp0$id,as.character(temp0$BNF_section))))]
      ))))
  
  
  
  count_data = count_data %>%
    full_join(temp, by=c("id"))
  
  # Fill in NAs with 0
  count_data[is.na(count_data)] = 0
  
  # make the features
  tmp$matrix = count_data %>%
    transmute(
      id=id,
      num_bnf_total = num_bnf_total,
      num_bnf_sections = all_sections,
    )
  tmp$missing_default = list(
    num_bnf_total = 0,
    num_bnf_sections = 0
  )
  return(tmp)
}

##################################################################
##                   Binary count BNF features                  ##
##################################################################

bnf_binary_count = function(PIS, section_list,feature_name) {
  # Performs an addition on presence / absence of things in section_list
  # i.e., the number of unique things in section_list that are present
  # for this ID
  # see e.g. specs for BNF_VITAMINS and BNF_BANDAGES on p45 of spec
  tmp = list()
  tmp$matrix = PIS %>%
    select(id, BNF_section, NUM_sold) %>%
    filter(BNF_section %in% section_list) %>%
    group_by(id) %>%
    summarise(!!feature_name := n_distinct(BNF_section))
  tmp$missing_default = list()
  tmp$missing_default[feature_name]=0
  return(tmp)
}

##################################################################
##                  Standard count BNF features                 ##
##################################################################

bnf_count = function(PIS, section_list, feature_name) {
  tmp = list()
  tmp$matrix = PIS %>%
    select(id, BNF_section, NUM_sold) %>%
    filter(BNF_section %in% section_list) %>%
    group_by(id) %>%
    summarise(!!feature_name := sum(NUM_sold))
  tmp$missing_default = list()
  tmp$missing_default[feature_name] = 0
  return(tmp)
}

##################################################################
##                  AE2 attendance                              ##
##################################################################

## Note:: ATTEDANCE_CATEGORY might not be present in extract

ae2_attendances = function(AE2) {
  tmp = list()
  tmp$matrix = AE2 %>% 
    # filter(ATTENDANCE_CATEGORY ==1 | ATTENDANCE_CATEGORY == 3) %>% ### COMMENTED OUT BY JAMES ON NSH: We do not have this field
    group_by(id) %>%
    summarise(num_ae2_attendances = n())
  tmp$missing_default=list(
    num_ae2_attendances = 0
  )
  return(tmp)
}

##################################################################
##                  Psych admissions                            ##
##################################################################

psych_admissions = function(SMR04) {
  tmp = list()
  tmp$matrix = SMR04 %>%
    group_by(id) %>%
    summarise(
      num_psych_admissions = n_distinct(CIS_MARKER)
    )
  tmp$missing_default = list(
    num_psych_admissions = 0
  )
  return(tmp)
}

##################################################################
##                Admission and bed days features               ##
##################################################################

admissions_and_bed_days = function(SMR01_merged) {
  tmp = list()
  tmp$matrix = SMR01_merged %>%
    group_by(id) %>%
    summarise(
      num_emergency_admissions = sum(emergency_admin),
      emergency_bed_days = sum(emergency_admin * LENGTH_OF_STAY),
      num_elective_admissions = sum(elective_admin),
      elective_bed_days = sum(elective_admin * LENGTH_OF_STAY),
      alcohol_substance = sum(emergency_alcohol_substance),
      selfharm = sum(emergency_selfharm),
      num_dc_admissions = sum(dc_admin),
      num_alcohol_admissions = sum(emergency_alcohol),
      num_elective_dc_admissions = sum(elective_admin) + sum(dc_admin)
    )
  tmp$missing_default = list(
    num_emergency_admissions = 0,
    emergency_bed_days = 0,
    num_elective_admissions = 0,
    elective_bed_days = 0,
    alcohol_substance = 0,
    num_dc_admissions = 0,
    num_alcohol_admissions = 0,
    num_elective_dc_admissions = 0
  )
  return(tmp)
}

##################################################################
##         Features indicated by LTC or prescription            ##
##################################################################

indicated_features = function(list_of_data_tables) {
  # MS
  indicated_table = indicated(
    list_of_data_tables$SPARRALTC,
    list_of_data_tables$PIS,
    "FIRST_MULTIPLE_SCLEROSIS_EPISODE",
    NULL
  ) %>% transmute(id = id, MS_indicated = indicated)
  
  # Parkinsons
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_PARKINSONS_DISEASE_EPISODE",
        NULL
      ) %>% transmute(id = id, parkinsons_indicated = indicated),
      by = c("id")
    )
  
  # Epilepsy
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_EPILEPSY_EPISODE",
        NULL
      ) %>% transmute(id = id, epilepsy_indicated = indicated),
      by = c("id")
    )
  
  
  # Dementia
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_DEMENTIA_EPISODE",
        NULL
      ) %>% transmute(id = id, dementia_indicated = indicated),
      by = c("id")
    )
  
  # Asthma
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_ASTHMA_EPISODE",
        NULL
      ) %>% transmute(id = id, asthma_indicated = indicated),
      by = c("id")
    )
  
  # Diabetes
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_DIABETES_EPISODE",
        NULL
      ) %>% transmute(id = id, diabetes_indicated = indicated),
      by = c("id")
    )
  
  # Cancer
  indicated_table = indicated_table %>%
    full_join(
      indicated(
        list_of_data_tables$SPARRALTC,
        list_of_data_tables$PIS,
        "FIRST_CANCER_EPISODE",
        NULL
      ) %>% transmute(id = id, cancer_indicated = indicated),
      by = c("id")
    )
  
  
  indicated_table[is.na(indicated_table)] = 0
  
  tmp = list()
  tmp$matrix = indicated_table
  tmp$missing_default = list(
    MS_indcated = 0,
    parkinsons_indicated = 0,
    epilepsy_indicated = 0,
    epilepsy_indicated_yed = 0,
    asthma_indicated = 0,
    diabetes_indicated = 0,
    cancer_indicated = 0,
    blood_indicated = 0,
    endocrine_indicated = 0,
    digestive_indicated = 0,
    congenital_indicated = 0
  )
  return(tmp)	
}

# Method that adds a 1 to an ID if a particular ltc_term is in the ltc_data
# OR a particular prescription is present at least once
indicated = function(ltc_data,pis_data,ltc_term,pis_term) {
  # Check for the LTC
  indicated_ltc = ltc_data %>% filter(LTC_TYPE == ltc_term) %>%
    transmute(id = id,ltc = 1)
  # Check the PIC
  if (!is.null(pis_term)) {
    indicated_pis = pis_data %>% filter(BNF_section == pis_term) %>%
      group_by(id) %>%
      filter(row_number() == 1) %>% # note - just take the first if there is more than one
      ungroup() %>%
      transmute(id=id,pis=1)
    
    # Combine the two
    indicated = indicated_ltc %>%
      full_join(indicated_pis, by = c("id")) %>%
      mutate(
        ltc = replace_na(ltc,0),
        pis = replace_na(pis,0)
      ) %>%
      transmute(id = id, indicated = ltc + pis)
  } else {
    indicated=indicated_ltc %>%
      mutate(
        ltc = replace_na(ltc,0),
      ) %>%
      transmute(id = id, indicated = ltc)
  }
  return(indicated)
}

##################################################################
##        SMR01 and Systemwatch merge with flag creation        ##
##################################################################

merge_smr01 = function(list_of_data_tables,time_cutoff) {  #### ARGUMENT CHANGED BY JL ON NSH: ADDED time_cutoff
  # Merges SMR01 and SystemWatch to create a 
  # new SMR01
  # Also adds emergency drug and alcohol flags
  list_of_data_tables$SMR01 = list_of_data_tables$SMR01 %>%
    group_by(id,CIS_MARKER) %>%
    mutate(LENGTH_OF_STAY = sum(LENGTH_OF_STAY)) %>%
    filter(admission_type != 18) %>%
    ungroup()
  
  # Add a LOCATION_CODE column to align with SystemWatch
  list_of_data_tables$SMR01 = list_of_data_tables$SMR01 %>%
    #mutate(LOCATION_CODE = EPISODE_LOCATION) #### CHANGED IN NSH BY JL TO SUBSEQUENT LINE
    mutate(LOCATION_CODE = LOCATION)
  
  # Replace any NAs in SystemWatch time_discharge column with the time cutoff
  list_of_data_tables$SystemWatch = list_of_data_tables$SystemWatch %>%
    replace_na(list(time_discharge = as.POSIXct(time_cutoff,origin = "1970-01-01")))
  
  # Add a LENGTH_OF_STAY column to SystemWatch
  list_of_data_tables$SystemWatch$LENGTH_OF_STAY = as.numeric(
    list_of_data_tables$SystemWatch$time_discharge - list_of_data_tables$SystemWatch$time
  )/(24*3600)
  
  # Add an INPATIENT_DAYCASE_IDENTIFIER column to SW
  # This avoids NAs after rowbinding
  list_of_data_tables$SystemWatch$INPATIENT_DAYCASE_IDENTIFIER = "SW"
  
  ### ADDED BY JL ON NSH: column ADMISSION_TRANSFER_FROM and MANAGEMENT_OF_PATIENT needs to be character type in SystemWatch to merge
  list_of_data_tables$SystemWatch$ADMISSION_TRANSFER_FROM = as.character(list_of_data_tables$SystemWatch$ADMISSION_TRANSFER_FROM)
  list_of_data_tables$SystemWatch$MANAGEMENT_OF_PATIENT = as.character(list_of_data_tables$SystemWatch$MANAGEMENT_OF_PATIENT)
  
  # Add systemwatch records to SMR01
  # Systemwatch records with (id, time, LOCATION_CODE) that
  # don't exist in SMR01 should be added.
  # The anti join finds them
  sw_rows_to_add = anti_join(
    list_of_data_tables$SystemWatch,
    list_of_data_tables$SMR01,
    by = c("id","time","LOCATION_CODE")
  )
  # Add the rows found by anti_join
  list_of_data_tables$SMR01 = list_of_data_tables$SMR01 %>%
    bind_rows(sw_rows_to_add)
  
  # Add flags
  icd10groupings = get_icd10_grouping_drug_alcohol_selfharm()
  selfharm_grouping=
    output = list_of_data_tables$SMR01 %>%
    mutate(
      emergency_admin=ifelse(
        ((admission_type %in% 20:22 | admission_type %in% 30:36 | admission_type==39) |
            (admission_type==38 & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))),1,0),
      elective_admin = ifelse(
        ((admission_type %in% c(10:12,19)) & (INPATIENT_DAYCASE_IDENTIFIER=="IP")) |
          ((admission_type %in% c(10:12,19)) & (source_table=="SystemWatch")),1,0),
      dc_admin = ifelse(
        INPATIENT_DAYCASE_IDENTIFIER=="DC",1,0
      ),
      alcohol_admin = ifelse((main_condition %in% icd10groupings$alcohol|OTHER_CONDITION_1 %in% icd10groupings$alcohol|
          OTHER_CONDITION_2 %in% icd10groupings$alcohol|OTHER_CONDITION_3 %in% icd10groupings$alcohol|
          OTHER_CONDITION_4 %in% icd10groupings$alcohol|OTHER_CONDITION_5 %in% icd10groupings$alcohol),1,0),
      drug_admin = ifelse(((main_condition %in% icd10groupings$drug|OTHER_CONDITION_1 %in% icd10groupings$drug|
          OTHER_CONDITION_2 %in% icd10groupings$drug|OTHER_CONDITION_3 %in% icd10groupings$drug|
          OTHER_CONDITION_4 %in% icd10groupings$drug|OTHER_CONDITION_5 %in% icd10groupings$drug)),1,0),
      selfharm_admin = ifelse(((main_condition %in% icd10groupings$selfharm|OTHER_CONDITION_1 %in% icd10groupings$selfharm|
          OTHER_CONDITION_2 %in% icd10groupings$selfharm|OTHER_CONDITION_3 %in% icd10groupings$selfharm|
          OTHER_CONDITION_4 %in% icd10groupings$selfharm|OTHER_CONDITION_5 %in% icd10groupings$selfharm)),1,0),
    ) %>%
    mutate(
      emergency_alcohol_substance=ifelse(
        emergency_admin==1 & (alcohol_admin==1|drug_admin==1),1,0
      ),
      emergency_alcohol = ifelse(
        emergency_admin==1 & alcohol_admin==1,1,0
      ),
      emergency_selfharm = ifelse(
        emergency_admin==1 & selfharm_admin==1,1,0
      )
    )
  
  return(output)
}

##################################################################
##                   General outpatient                         ##
##################################################################

# Note: REFERRAL_TYPE < 3

general_outpatient = function(SMR00) {
  codes_to_include = c(
    "A1", "A2", "A3", "A6", "A7", "A8", "A81", "A82", 
    "A9", "AB", "AD", "H2", "AF", "CA", "AG", "AH", "AM",
    "AP", "AQ", "AR", "C1", "C11", "C12", "C13", "C3",
    "C31", "C4", "C41", "C42", "C5", "C51", "C6", "C7", "C8",
    "C9", "CB", "D3", "D4", "D6", "E12", "F2", "G1", "G1A",
    "G2", "G21", "G22", "G3", "G4", "G5", "G6", "H1", "J3",
    "J4", "J5"
  )
  num_op_appointments = SMR00 %>%
    filter(
      SPECIALTY %in% codes_to_include,
      REFERRAL_TYPE < 3
    ) %>%
    group_by(id) %>%
    summarise(num_outpatient_appointment_general = n())
  
  tmp= list()
  tmp$matrix = num_op_appointments
  tmp$missing_default = list(num_outpatient_appointment_general = 0)
  
  return(tmp)
}



##################################################################
##               Psychiatric outpatient                         ##
##################################################################

# Note: REFERRAL_TYPE < 3

psych_outpatient = function(SMR00) {
  codes_to_include = c(
    "G1", "G1A", "G2", "G21", "G22", "G3", "G4", "G5", "G6"
  )
  num_psych_op_appointments = SMR00 %>%
    filter(
      SPECIALTY %in% codes_to_include,
      REFERRAL_TYPE < 3
    ) %>%
    group_by(id) %>%
    summarise(num_outpatient_appointment_mental = n())
  
  tmp= list()
  tmp$matrix = num_psych_op_appointments
  tmp$missing_default = list(num_outpatient_appointment_mental = 0)
  
  return(tmp)
}







##################################################################
##                         LTC features                         ##
##################################################################

general_ltc_names=c("FIRST_ARTHRITIS_EPISODE", 
  "FIRST_ASTHMA_EPISODE",
  "FIRST_ATRIAL_FIBRILLATION_EPISODE", 
  "FIRST_COPD_EPISODE", 
  "FIRST_CANCER_EPISODE", 
  "FIRST_CEREBROVASCULAR_DISEASE_EPISODE", 
  "FIRST_CHRONIC_LIVER_DISEASE_EPISODE",
  "FIRST_DEMENTIA_EPISODE", 
  "FIRST_DIABETES_EPISODE",
  "FIRST_EPILEPSY_EPISODE", 
  "FIRST_HEART_DISEASE_EPISODE", 
  "FIRST_HEART_FAILURE_EPISODE", 
  "FIRST_MULTIPLE_SCLEROSIS_EPISODE", 
  "FIRST_PARKINSON_DISEASE_EPISODE", 
  "FIRST_RENAL_FAILURE_EPISODE"
)

# Additional yed codes
yed_ltc_names = c(
  "FIRST_CONGENITAL_PROB_EPISODE",
  "FIRST_DIS_BLOOD_EPISODE",
  "FIRST_ENDOCRINE_MET_EPISODE",
  "FIRST_OTHER_DIGESTIVE_EPISODE"
)

ltc_features = function(ltc_table) {
  # ltc_table should be the sparse-style table
  ltc_sub_table =ltc_table %>%
    select(id, LTC_TYPE, time)
  
  wide_table = ltc_sub_table %>%
    pivot_wider(names_from = "LTC_TYPE",values_from="time")
  
  general_ltc_names_sub = c()
  # filter names as some might not be present
  for (name in general_ltc_names) {
    if (name %in% names(wide_table)) {
      general_ltc_names_sub = c(general_ltc_names_sub,name)
    }
  }
  yed_ltc_names_sub=c()
  # filter names as some might not be present
  for (name in yed_ltc_names) {
    if (name %in% names(wide_table)) {
      yed_ltc_names_sub = c(yed_ltc_names_sub,name)
    }
  }
  
  wide_table$num_general_ltc = rowSums(!is.na(wide_table[general_ltc_names_sub]))
  wide_table$num_yed_ltc = wide_table$num_general_ltc #+ rowSums(!is.na(wide_table[yed_ltc_names_sub])) ### CHANGED BY JL: we have none of the YED/U16 specific LTCs
  
  wide_table = wide_table %>%
    select(id,num_general_ltc,num_yed_ltc)
  
  tmp = list()
  tmp$matrix = wide_table
  tmp$missing_default = list(
    num_general_ltc = 0,
    num_yed_ltc = 0
  )
  return(tmp)
}



##################################################################
##                      PIS definitions                         ##
##################################################################

bnf_groupings=list()
bnf_groupings$pis_gastro_int=c(
  "NUM_BNF_0101","NUM_BNF_0102","NUM_BNF_0103","NUM_BNF_0105",
  "NUM_BNF_0106","NUM_BNF_0107","NUM_BNF_0109"
)
bnf_groupings$pis_respiratory = c("NUM_BNF_0301", "NUM_BNF_0302", 
  "NUM_BNF_0303", "NUM_BNF_0304", "NUM_BNF_0305", "NUM_BNF_0306",
  "NUM_BNF_0307", "NUM_BNF_0308", "NUM_BNF_0309", "NUM_BNF_0310")
bnf_groupings$pis_cns = c("NUM_BNF_0401", "NUM_BNF_0402", "NUM_BNF_0403",
  "NUM_BNF_0404", "NUM_BNF_0405", "NUM_BNF_0406", "NUM_BNF_0407",
  "NUM_BNF_0408", "NUM_BNF_0409", "NUM_BNF_0410", "NUM_BNF_0411")
bnf_groupings$pis_infections = c("NUM_BNF_0501", "NUM_BNF_0502", "NUM_BNF_0503", "NUM_BNF_0504",
  "NUM_BNF_0505")
bnf_groupings$pis_endocrine = c( "NUM_BNF_0601", "NUM_BNF_0602", "NUM_BNF_0603",
  "NUM_BNF_0604", "NUM_BNF_0605", "NUM_BNF_0606", "NUM_BNF_0607")
bnf_groupings$pis_incontinence = c( "NUM_BNF_2201", "NUM_BNF_2202",
  "NUM_BNF_2205", "NUM_BNF_2210", "NUM_BNF_2215", "NUM_BNF_2220",
  "NUM_BNF_2230", "NUM_BNF_2240", "NUM_BNF_2250", "NUM_BNF_2260",
  "NUM_BNF_2270", "NUM_BNF_2280", "NUM_BNF_2285", "NUM_BNF_2290")
bnf_groupings$pis_stoma = c( "NUM_BNF_2305", "NUM_BNF_2310", "NUM_BNF_2315", "NUM_BNF_2320",
  "NUM_BNF_2325", "NUM_BNF_2330", "NUM_BNF_2335", "NUM_BNF_2340",
  "NUM_BNF_2345", "NUM_BNF_2346", "NUM_BNF_2350", "NUM_BNF_2355",
  "NUM_BNF_2360", "NUM_BNF_2365", "NUM_BNF_2370", "NUM_BNF_2375",
  "NUM_BNF_2380", "NUM_BNF_2385", "NUM_BNF_2390", "NUM_BNF_2392",
  "NUM_BNF_2393", "NUM_BNF_2394", "NUM_BNF_2396", "NUM_BNF_2398")
bnf_groupings$pis_nutrition = c(
  "NUM_BNF_0908","NUM_BNF_0909","NUM_BNF_0911","NUM_BNF_0912"
)
bnf_groupings$pis_skin = c(
  "NUM_BNF_1301" ,"NUM_BNF_1302" ,"NUM_BNF_1303" ,"NUM_BNF_1304" ,
  "NUM_BNF_1305" ,"NUM_BNF_1308" ,"NUM_BNF_1310" ,"NUM_BNF_1311" ,
  "NUM_BNF_1313" ,"NUM_BNF_1314"
)
bnf_groupings$pis_supplements = c("NUM_BNF_0911","NUM_BNF_0912")
bnf_groupings$pis_vitamins=c("NUM_BNF_0905","NUM_BNF_0906")
bnf_groupings$pis_bandages=c("NUM_BNF_2002","NUM_BNF_2003")

bnf_groupings$pis_gut_motility=c("NUM_BNF_0102")
bnf_groupings$pis_antisecretory=c("NUM_BNF_0103")
bnf_groupings$pis_intestinal=c("NUM_BNF_0109")
bnf_groupings$pis_heart=c("NUM_BNF_0208")
bnf_groupings$pis_anticoagulant=c("NUM_BNF_0208")
bnf_groupings$pis_antifibrinolytic=c("NUM_BNF_0211")
bnf_groupings$pis_respiratory_corticosteroids=c("NUM_BNF_0302")
bnf_groupings$pis_bronco=c("NUM_BNF_0301")
bnf_groupings$pis_cromo=c("NUM_BNF_0303")
bnf_groupings$pis_mucolytics=c("NUM_BNF_0307")
bnf_groupings$pis_antiepileptic_Drugs=c("NUM_BNF_0408")
bnf_groupings$pis_parkinsonism=c("NUM_BNF_0409")
bnf_groupings$pis_sub_depend=c("NUM_BNF_0410")
bnf_groupings$pis_dementia=c("NUM_BNF_0411")
bnf_groupings$pis_antibacterial=c("NUM_BNF_0501")
bnf_groupings$pis_diabetes=c("NUM_BNF_0601")
bnf_groupings$pis_corticosteroids=c("NUM_BNF_0603")
bnf_groupings$pis_fluids=c("NUM_BNF_0902")
bnf_groupings$pis_Minerals=c("NUM_BNF_0905")
bnf_groupings$pis_rheumatic=c("NUM_BNF_1001")
bnf_groupings$pis_neuromuscular=c("NUM_BNF_1002")
bnf_groupings$pis_mydriatics=c("NUM_BNF_1105")
bnf_groupings$pis_catheters=c("NUM_BNF_2102")
bnf_groupings$pis_inotropic = c("NUM_BNF_0201")
bnf_groupings$pis_diuretics = c("NUM_BNF_0202")
bnf_groupings$pis_antiarrhythmics = c("NUM_BNF_0203")
bnf_groupings$pis_betablockers = c("NUM_BNF_0204")
bnf_groupings$pis_hypertensive_heart_failure = c("NUM_BNF_0205")
bnf_groupings$pis_antianginal = c("NUM_BNF_0206")
bnf_groupings$pis_antiplatelets = c("NUM_BNF_0209")
bnf_groupings$pis_lipid = c("NUM_BNF_0212")
bnf_groupings$pis_genitourinary = c("NUM_BNF_0704")
bnf_groupings$pis_cytotoxics = c("NUM_BNF_0801")
bnf_groupings$pis_immune = c("NUM_BNF_0802")
bnf_groupings$pis_sex_hormone_antagonists = c("NUM_BNF_0803")
bnf_groupings$pis_antianaemics = c("NUM_BNF_0901")
bnf_groupings$pis_topical_pain_relief = c("NUM_BNF_1003")
bnf_groupings$pis_antibacterial_eyes = c("NUM_BNF_1103")
bnf_groupings$pis_antiinflammatory_corticosteroids = c("NUM_BNF_1104")
bnf_groupings$pis_glaucoma = c("NUM_BNF_1106")
bnf_groupings$pis_local_anaesthetics = c("NUM_BNF_1107")
bnf_groupings$pis_ophthalmic = c("NUM_BNF_1108")
bnf_groupings$pis_ear = c("NUM_BNF_1201")
bnf_groupings$pis_nose = c("NUM_BNF_1202")
bnf_groupings$pis_oropharynx = c("NUM_BNF_1203")
bnf_groupings$pis_hosiery = c("NUM_BNF_2107")
bnf_groupings$pis_metabolic = c("NUM_BNF_0908")
bnf_groupings$pis_food = c("NUM_BNF_0909")



# BNF sections for the total count
# has various fields: individaul = all individual codes in the table on p31 of specs
# groups = codes that should be considered as a group
# excluded = the codes that should be excluded for ltc and fe
bnf_section_ref=list()
bnf_section_ref$individual=c(
  "NUM_BNF_0101","NUM_BNF_0102","NUM_BNF_0103","NUM_BNF_0105",
  "NUM_BNF_0106","NUM_BNF_0107","NUM_BNF_0109",
  
  "NUM_BNF_0201" ,"NUM_BNF_0202" ,"NUM_BNF_0203" ,"NUM_BNF_0204" ,
  "NUM_BNF_0205" ,"NUM_BNF_0206" ,"NUM_BNF_0208" ,"NUM_BNF_0209" ,
  "NUM_BNF_0211" ,"NUM_BNF_0212",
  
  "NUM_BNF_0301" ,"NUM_BNF_0302" ,"NUM_BNF_0303" ,"NUM_BNF_0304" ,
  "NUM_BNF_0306" ,"NUM_BNF_0307" ,"NUM_BNF_0309" ,"NUM_BNF_0310",
  
  "NUM_BNF_0401" ,"NUM_BNF_0402" ,"NUM_BNF_0403" ,"NUM_BNF_0404" ,
  "NUM_BNF_0405" ,"NUM_BNF_0406" ,"NUM_BNF_0407" ,"NUM_BNF_0408" ,
  "NUM_BNF_0409" ,"NUM_BNF_0410" ,"NUM_BNF_0411",
  
  "NUM_BNF_0501" ,"NUM_BNF_0502" ,"NUM_BNF_0503" ,"NUM_BNF_0504" ,
  "NUM_BNF_0505",
  
  "NUM_BNF_0601" ,"NUM_BNF_0602" ,"NUM_BNF_0603" ,"NUM_BNF_0604" ,
  "NUM_BNF_0605" ,"NUM_BNF_0605" ,"NUM_BNF_0607",
  
  "NUM_BNF_0704",
  
  "NUM_BNF_0801" ,"NUM_BNF_0802" ,"NUM_BNF_0803",
  
  "NUM_BNF_0901" ,"NUM_BNF_0902" ,"NUM_BNF_0904" ,"NUM_BNF_0905" ,
  "NUM_BNF_0906",
  
  "NUM_BNF_1001" ,"NUM_BNF_1002" ,"NUM_BNF_1003",
  
  "NUM_BNF_1103" ,"NUM_BNF_1104" ,"NUM_BNF_1105" ,"NUM_BNF_1106" ,
  "NUM_BNF_1107" ,"NUM_BNF_1108",
  
  "NUM_BNF_1201" ,"NUM_BNF_1202" ,"NUM_BNF_1203",
  
  "NUM_BNF_2002" ,"NUM_BNF_2003",
  
  "NUM_BNF_2102" ,"NUM_BNF_2107"
)


bnf_section_ref$groups = list()
bnf_section_ref$groups$nutrition = c(
  "NUM_BNF_0908","NUM_BNF_0909","NUM_BNF_0911","NUM_BNF_0912"
)
bnf_section_ref$groups$skin = c(
  "NUM_BNF_1301" ,"NUM_BNF_1302" ,"NUM_BNF_1303" ,"NUM_BNF_1304" ,
  "NUM_BNF_1305" ,"NUM_BNF_1308" ,"NUM_BNF_1310" ,"NUM_BNF_1311" ,
  "NUM_BNF_1313" ,"NUM_BNF_1314"
)
bnf_section_ref$groups$combined_X = c(
  "NUM_BNF_2201" ,"NUM_BNF_2202" ,"NUM_BNF_2205" ,"NUM_BNF_2210" ,
  "NUM_BNF_2215" ,"NUM_BNF_2220" ,"NUM_BNF_2230" ,"NUM_BNF_2240" ,
  "NUM_BNF_2250" ,"NUM_BNF_2260" ,"NUM_BNF_2270" ,"NUM_BNF_2280" ,
  "NUM_BNF_2285" ,"NUM_BNF_2290"
)
bnf_section_ref$groups$combined_Y = c(
  "NUM_BNF_2305" ,"NUM_BNF_2310" ,"NUM_BNF_2315" ,"NUM_BNF_2320" ,
  "NUM_BNF_2325" ,"NUM_BNF_2330" ,"NUM_BNF_2335" ,"NUM_BNF_2340" ,
  "NUM_BNF_2345" ,"NUM_BNF_2346" ,"NUM_BNF_2350" ,"NUM_BNF_2355" ,
  "NUM_BNF_2360" ,"NUM_BNF_2365" ,"NUM_BNF_2370" ,"NUM_BNF_2375" ,
  "NUM_BNF_2380" ,"NUM_BNF_2385" ,"NUM_BNF_2390" ,"NUM_BNF_2392" ,
  "NUM_BNF_2393" ,"NUM_BNF_2394" ,"NUM_BNF_2396" ,"NUM_BNF_2398"
)

bnf_section_ref$excluded = list()
bnf_section_ref$excluded$excluded_ltc=c(
  "NUM_BNF_0309" ,"NUM_BNF_0310" ,"NUM_BNF_0503" ,"NUM_BNF_1103" ,
  "NUM_BNF_1104" ,"NUM_BNF_1105" ,"NUM_BNF_1107" ,"NUM_BNF_1201" ,
  "NUM_BNF_1202"
)
bnf_section_ref$excluded$excluded_fe = c(
  "NUM_BNF_0404" ,"NUM_BNF_0405" ,"NUM_BNF_0410" ,"NUM_BNF_0505" ,
  "NUM_BNF_0604" ,"NUM_BNF_0605" ,"NUM_BNF_0607" ,"NUM_BNF_0801" ,
  "NUM_BNF_0803" ,"NUM_BNF_1105" ,"NUM_BNF_1107" ,"NUM_BNF_1201" ,
  "NUM_BNF_1202" ,"NUM_BNF_1203"
)



##################################################################
##              ICD10 groupings                                 ##
##################################################################

get_icd10_grouping_drug_alcohol_selfharm=function() {
  
  list(
    alcohol = c(
      "E244", "E512", "F100", "F102", "F103", "F104", 
      "F105", "F106", "F107", "F108", "F109", "G312", "G621", "G721", 
      "I426", "K292", "K700", "K701", "K702", "K703", "K704", "K705", 
      "K706", "K707", "K708", "K709", "K860", "O354", "P043", "Q860", 
      "R780", "T510", "T511", "T519", "X450", "X451", "X452", "X453", 
      "X454", "X455", "X456", "X457", "X458", "X459", "X650", "X651", 
      "X652", "X653", "X654", "X655", "X656", "X657", "X658", "X659", 
      "Y150", "Y151", "Y152", "Y153", "Y154", "Y155", "Y156", "Y157", 
      "Y158", "Y159", "Y573", "Y900", "Y902", "Y903", "Y904", "Y905", 
      "Y906", "Y907", "Y908", "Y909", "Y910", "Y911", "Y912", "Y913", 
      "Y914", "Y915", "Y916", "Y917", "Y918", "Y919", "Z502", "Z714", 
      "Z721"), 
    drug = c(
      "F110", "F111", "F112", "F113", "F114", "F115", 
      "F116", "F117", "F118", "F119", "F120", "F121", "F122", "F123", 
      "F124", "F125", "F126", "F127", "F128", "F129", "F130", "F131", 
      "F132", "F133", "F134", "F135", "F136", "F137", "F138", "F139", 
      "F140", "F141", "F142", "F143", "F144", "F145", "F146", "F147", 
      "F148", "F149", "F150", "F151", "F152", "F153", "F154", "F155", 
      "F156", "F157", "F158", "F159", "F160", "F161", "F162", "F163", 
      "F164", "F165", "F166", "F167", "F168", "F169", "F180", "F181", 
      "F182", "F183", "F184", "F185", "F186", "F187", "F188", "F189", 
      "F190", "F191", "F192", "F193", "F194", "F195", "F196", "F197", 
      "F198", "F199"),
    selfharm = paste0("X",600:849)
  )
}





##################################################################
##              Patient info (age, SIMD etc)                    ##
##################################################################

#' transformer_patient_info
#' Returns basic information about each patient as a feature
#' That includes age, gender and SIMD_DECILE_2016_SCT
#' The latter is a categorised deprivation status
#' 
transformer_patient_info <- function(
  # Input data
  patients, episodes, list_of_data_tables,
  time_cutoff # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  
){
  
  
  # Create the out_matrix
  out_matrix = patients %>% 
    transmute(
      id = id,
      gender = gender,
      age = as.integer((as.integer(as_date(as_datetime(time_cutoff))) - as.integer(as_date(date_of_birth)))/365),
      SIMD_DECILE_2016_SCT = SIMD_DECILE_2016_SCT
    ) #%>%
  #mutate(age = ifelse(age<0, 0, age))
  
  out_missing_default = as.list(rep(NA, length=length(names(out_matrix)[-1])))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}


##################################################################
##              Turn matrix into factors                        ##
##################################################################


##' turn_data_into_sparrav3_factors() 
##' Turn the SPARRAv3 columns into factor levels (as used by SPARRAv3). Written by James with aid of v3 coefficients.
##'  Also renames variables to those used in spreadsheet. 
##' @param df input matrix of type used in SPARRAv4
##' @param mode sPARRAv3 class: factor levels differ between classes
##' @return matrix for SPARRAv3 computation
turn_data_into_sparrav3_factors <- function(df,mode="fe"){
  df0=df;
  if (mode=="fe") {
    df=df %>%
      mutate(
        
        # Demographics
        AGE = cut(age, breaks = c(-Inf, seq(75,89), Inf),labels=c(75:89,"90+")),
        # GENDER=factor(2-sexM),
        SIMD_QUINTILE=factor(floor((as.numeric(SIMD_DECILE_2016_SCT)+1)/2)),
        
        # Attendances
        ED_ATTENDANCES = cut(num_ae2_attendances, breaks = c(-Inf, 0:9, Inf),labels=c(0:9,"10+")),
        EMERGENCY_ADMISSIONS = cut(num_emergency_admissions, breaks = c(-Inf, 0:5, Inf),labels=c(0:5,"6+")), # This is the proper version, but it has 3-year lookback
        EMERGENCY_BEDDAYS = cut(emergency_bed_days, breaks = c(-Inf, 0:14, Inf),labels=c(0:14,"15+")),
        #ELECTIVE_ADMISSIONS = cut(num_elective_admissions, breaks = c(-Inf, 0:2, Inf),labels=c(0:2,"3+")),
        #ELECTIVE_BEDDAYS = cut(elective_bed_days, breaks = c(-Inf, 0:28, Inf),labels=c(0:28,"29+")),
        OP_APPOINTMENTS = cut(num_outpatient_appointment_general, breaks = c(-Inf, 0:4, Inf),labels=c(0:4,"5+")),
        #OP_PSYCHIATRIC_APPOINTMENTS = cut(num_outpatient_appointment_mental,c(-Inf,0,Inf),labels=c("0","1+")),
        #DC_ADMISSIONS=cut(num_dc_elective_admissions+num_dc_emergency_admissions,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        ELECTIVE_DC_ADMISSIONS=cut(num_elective_dc_admissions,c(-Inf,0:6,Inf),labels=c(0:6,"7+")), 
        #ALCOHOL_SUB_MISUSE_ADMISSIONS = cut(emergency_drugAndalcohol_admin,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        ALCOHOL_ADMISSIONS =  cut(num_alcohol_admissions, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #PSYCHIATRIC_ADMISSIONS = cut(num_psych_admissions,c(-Inf,0,Inf),labels=c("0","1+")),
        
        # LTCs
        LTC=cut(num_general_ltc, breaks=c(-Inf,0:5,Inf),labels=c(0:5,"6+")),
        PARKINSONS_DISEASE=cut(parkinsons_indicated,breaks=c(-Inf,0,Inf),labels=c(0:1)),
        #MS=cut(MS_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        #EPILEPSY=cut(epilepsy_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        #DEMENTIA=cut(dementia_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        #DEMENTIA 
        #ASTHMA 
        #DIABETES 
        #CANCER
        #ENDOCRINE 
        #DIGESTIVE 
        #BLOOD 
        #CONGENITAL 
        #SELF_HARM 
        
        # BNF
        BNF_SECTIONS = cut(num_bnf_sections_fe, breaks = c(-Inf, 0:19, Inf),labels=c(0:19,"20+")),
        BNF_RESPIRATORY = cut(pis_respiratory, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_CNS = cut(pis_cns, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_INFECTIONS = cut(pis_infections, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_ENDOCRINE = cut(pis_endocrine, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        #BNF_SUBSTANCE = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_DEMENTIA = cut(pis_Dementia, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_CORTICOSTEROIDS = cut(pis_Corticosteroids_Respiratory, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_FLUIDS = cut(pis_Fluids_And_Electrolytes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NUTRITION = cut(pis_Oral_Nutrition, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_VITAMINS = cut(pis_Vitamins, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_BANDAGES = cut(pis_Arm_Sling_Bandages, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_CATHETERS = cut(pis_Catheters, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_STOMA = cut(pis_CombinedStomaDevices, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_INCONTINENCE = cut(pis_incontinence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_GUT_MOTILITY = cut(pis_Antispasmod_Other_Drgs_Alt_Gut_Motility, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTISECRETORY = cut(pis_Antisecretory_Drugs_Mucosal_Protectants, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_INTESTINAL = cut(pis_Drugs_Affecting_Intestinal_Secretions, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTICOAGULENT = cut(pis_Anticoagulants_And_Protamine, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTIFIBRINOLYTIC = cut(`pis_Antifibrinolytic_Drugs _Haemostatics`, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_MUCOLYTICS = cut(pis_Mucolytics, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_DIABETES = cut(pis_Drugs_Used_In_Diabetes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_GASTRO_INT = cut(pis_XXX, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_HEART = cut(v3t$pis, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_BRONCHO = cut(, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_CROMO = cut(v3t$, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_SUB_DEPEND = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTIBACTERIAL = cut(v3t$pi, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_RHEUMATIC = cut(v3t$pis_, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Neuromuscular_Disorders, breaks = c(-Inf,0,Inf),labels=c("0","1+"))
        #BNF_MYDRIATICS = cut(, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
      )
    fecol=c("AGE", "ALCOHOL_ADMISSIONS", "BNF_CNS", "BNF_ENDOCRINE", "BNF_INCONTINENCE", 
      "BNF_INFECTIONS", "BNF_RESPIRATORY", "BNF_SECTIONS", "ED_ATTENDANCES", 
      "ELECTIVE_DC_ADMISSIONS", "EMERGENCY_ADMISSIONS", "EMERGENCY_BEDDAYS", 
      "LTC", "OP_APPOINTMENTS", "PARKINSONS_DISEASE", "SIMD_QUINTILE")
    df=df[,fecol]
  }
  if (mode=="ltc") {
    df=df %>%
      mutate(
        
        # Demographics
        # AGE = cut(age, breaks = c(-Inf, seq(0,89), Inf),labels=c(0:89,"90+")),
        #GENDER=factor(2-sexM),
        SIMD_QUINTILE=factor(floor((as.numeric(SIMD_DECILE_2016_SCT)+1)/2)),
        
        # Attendances
        ED_ATTENDANCES = cut(num_ae2_attendances, breaks = c(-Inf, 0:9, Inf),labels=c(0:9,"10+")),
        EMERGENCY_ADMISSIONS = cut(num_emergency_admissions, breaks = c(-Inf, 0:5, Inf),labels=c(0:5,"6+")), # This is the proper version, but it has 3-year lookback
        EMERGENCY_BEDDAYS = cut(emergency_bed_days, breaks = c(-Inf, 0:21, Inf),labels=c(0:21,"22+")),
        ELECTIVE_ADMISSIONS = cut(num_elective_admissions, breaks = c(-Inf, 0:2, Inf),labels=c(0:2,"3+")),
        ELECTIVE_BEDDAYS = cut(elective_bed_days, breaks = c(-Inf, 0:28, Inf),labels=c(0:28,"29+")),
        OP_APPOINTMENTS = cut(num_outpatient_appointment_general, breaks = c(-Inf, 0:6, Inf),labels=c(0:6,"7+")),
        # OP_PSYCHIATRIC_APPOINTMENTS = cut(num_outpatient_appointment_mental,c(-Inf,0,Inf),labels=c("0","1+")),
        DC_ADMISSIONS=cut(num_dc_admissions,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        # ELECTIVE_DC_ADMISSIONS=cut(num_dc_elective_admissions,c(-Inf,0:6,Inf),labels=c(0:6,"7+")), 
        ALCOHOL_SUB_MISUSE_ADMISSIONS = cut(alcohol_substance,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        # ALCOHOL_ADMISSIONS =  cut(alcohol_admin, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        PSYCHIATRIC_ADMISSIONS = cut(num_psych_admissions,c(-Inf,0,Inf),labels=c("0","1+")),
        
        # LTCs
        LTC=cut(num_general_ltc, breaks=c(-Inf,0:5,Inf),labels=c(0:5,"6+")),
        PARKINSONS_DISEASE=cut(parkinsons_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        MS=cut(MS_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        EPILEPSY=cut(epilepsy_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        # DEMENTIA=cut(dementia_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        #DEMENTIA 
        #ASTHMA 
        #DIABETES 
        #CANCER
        #ENDOCRINE 
        #DIGESTIVE 
        #BLOOD 
        #CONGENITAL 
        #SELF_HARM 
        
        # BNF
        BNF_SECTIONS = cut(num_bnf_sections_ltc, breaks = c(-Inf, 0:13, Inf),labels=c(0:13,"14+")),
        # BNF_RESPIRATORY = cut(pis_Respiratory, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_CNS = cut(pis_CentralNervousSystem, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_INFECTIONS = cut(pis_infections, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        #BNF_ENDOCRINE = cut(pis_EndocrineSystem, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_SUBSTANCE = cut(pis_sub_depend, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_DEMENTIA = cut(pis_dementia, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_CORTICOSTEROIDS = cut(pis_corticosteroids, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_FLUIDS = cut(pis_fluids, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_NUTRITION = cut(pis_nutrition, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_VITAMINS = cut(pis_vitamins_binary, breaks = c(-Inf,0:1,Inf),labels=c(0:1,"2+")),
        BNF_BANDAGES = cut(pis_bandages_binary, breaks = c(-Inf,0:1,Inf),labels=c(0:1,"2+")),
        BNF_CATHETERS = cut(pis_catheters, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_STOMA = cut(pis_stoma, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_INCONTINENCE = cut(pis_IncontinenceDevices, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_GUT_MOTILITY = cut(pis_Antispasmod_Other_Drgs_Alt_Gut_Motility, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_ANTISECRETORY = cut(pis_Antisecretory_Drugs_Mucosal_Protectants, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_INTESTINAL = cut(pis_Drugs_Affecting_Intestinal_Secretions, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_ANTICOAGULENT = cut(pis_Anticoagulants_And_Protamine, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_ANTIFIBRINOLYTIC = cut(`pis_Antifibrinolytic_Drugs _Haemostatics`, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_MUCOLYTICS = cut(pis_Mucolytics, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_DIABETES = cut(pis_Drugs_Used_In_Diabetes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_GASTRO_INT = cut(pis_XXX, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_HEART = cut(v3t$pis, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_BRONCHO = cut(, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_CROMO = cut(v3t$, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_SUB_DEPEND = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTIBACTERIAL = cut(v3t$pi, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_RHEUMATIC = cut(v3t$pis_, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Neuromuscular_Disorders, breaks = c(-Inf,0,Inf),labels=c("0","1+"))
        #BNF_MYDRIATICS = cut(, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
      )
    ltccol=c("ALCOHOL_SUB_MISUSE_ADMISSIONS", "BNF_BANDAGES", "BNF_CATHETERS", 
      "BNF_CORTICOSTEROIDS", "BNF_DEMENTIA", "BNF_FLUIDS", "BNF_INFECTIONS", 
      "BNF_NUTRITION", "BNF_SECTIONS", "BNF_STOMA", "BNF_SUBSTANCE", 
      "BNF_VITAMINS", "DC_ADMISSIONS", "ED_ATTENDANCES", "ELECTIVE_ADMISSIONS", 
      "ELECTIVE_BEDDAYS", "EMERGENCY_ADMISSIONS", "EMERGENCY_BEDDAYS", 
      "EPILEPSY", "LTC", "MS", "OP_APPOINTMENTS", "PARKINSONS_DISEASE", 
      "PSYCHIATRIC_ADMISSIONS", "SIMD_QUINTILE")
    df=df[,ltccol]
  }
  if (mode=="yed") {
    df=df %>%
      mutate(
        
        # Demographics
        #AGE = cut(age, breaks = c(-Inf, seq(0,89), Inf),labels=c(0:89,"90+")),
        #GENDER=factor(gender),
        # SIMD_QUINTILE=factor(floor((as.numeric(SIMD_DECILE_2016_SCT)+1)/2)),
        
        # Attendances
        ED_ATTENDANCES = cut(num_ae2_attendances, breaks = c(-Inf, 0:9, Inf),labels=c(0:9,"10+")),
        EMERGENCY_ADMISSIONS = cut(num_emergency_admissions, breaks = c(-Inf, 0:5, Inf),labels=c(0:5,"6+")), # This is the proper version, but it has 3-year lookback
        EMERGENCY_BEDDAYS = cut(emergency_bed_days, breaks = c(-Inf, 0:21, Inf),labels=c(0:21,"22+")),
        #ELECTIVE_ADMISSIONS = cut(num_elective_admissions, breaks = c(-Inf, 0:2, Inf),labels=c(0:2,"3+")),
        #ELECTIVE_BEDDAYS = cut(elective_bed_days, breaks = c(-Inf, 0:28, Inf),labels=c(0:28,"29+")),
        OP_APPOINTMENTS = cut(num_outpatient_appointment_general, breaks = c(-Inf, 0:4, Inf),labels=c(0:4,"5+")),
        OP_PSYCHIATRIC_APPOINTMENTS = cut(num_outpatient_appointment_mental,c(-Inf,0,Inf),labels=c("0","1+")),
        # DC_ADMISSIONS=cut(num_dc_elective_admissions+num_dc_emergency_admissions,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        ELECTIVE_DC_ADMISSIONS=cut(num_elective_dc_admissions,c(-Inf,0:6,Inf),labels=c(0:6,"7+")),
        ALCOHOL_SUB_MISUSE_ADMISSIONS = cut(alcohol_substance,c(-Inf,0:1,Inf),labels=c(0:1,"2+")),
        # ALCOHOL_ADMISSIONS =  cut(alcohol_admin, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        PSYCHIATRIC_ADMISSIONS = cut(num_psych_admissions,c(-Inf,0,Inf),labels=c("0","1+")),
        
        # LTCs
        LTC=cut(num_yed_ltc, breaks=c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        # PARKINSONS_DISEASE=cut(parkinsons_indicated,breaks=c(-Inf,0,Inf),labels=c(0:1)),
        # MS=cut(MS_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        # EPILEPSY=cut(epilepsy_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        DEMENTIA=cut(dementia_indicated,breaks=c(-Inf,0,Inf),labels=c(0:1)),
        #DEMENTIA
        #ASTHMA=0,
        #DIABETES
        #CANCER=0,
        #ENDOCRINE=0,
        #DIGESTIVE
        #BLOOD=0,
        #CONGENITAL=0,
        #SELF_HARM=0,
        
        # BNF
        BNF_SECTIONS = cut(num_bnf_sections_yed, breaks = c(-Inf, 0:13, Inf),labels=c(0:13,"14+")),
        # BNF_RESPIRATORY = cut(pis_Respiratory, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        BNF_CNS = cut(pis_cns, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_INFECTIONS = cut(pis_Infections, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_ENDOCRINE = cut(pis_EndocrineSystem, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_SUBSTANCE = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_DEMENTIA = cut(pis_Dementia, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_CORTICOSTEROIDS = cut(pis_corticosteroids, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_FLUIDS = cut(pis_fluids, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_NUTRITION = cut(pis_Oral_Nutrition, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_VITAMINS = cut(pis_vitamins, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_BANDAGES = cut(pis_Arm_Sling_Bandages, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_CATHETERS = cut(pis_Catheters, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_STOMA = cut(pis_stoma, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_INCONTINENCE = cut(pis_IncontinenceDevices, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_GUT_MOTILITY = cut(pis_gut_motility, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_ANTISECRETORY = cut(pis_antisecretory, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_INTESTINAL = cut(pis_intestinal, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_ANTICOAGULENT = cut(pis_anticoagulent, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_ANTIFIBRINOLYTIC = cut(pis_antifibrinolytic, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_MUCOLYTICS = cut(pis_mucolytics, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_DIABETES = cut(pis_diabetes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_GASTRO_INT = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_HEART = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_BRONCHO = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_CROMO = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_SUB_DEPEND = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTIBACTERIAL = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_RHEUMATIC = cut(v3t$pis_, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Neuromuscular_Disorders, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_MYDRIATICS = cut(0, breaks = c(-Inf,0,Inf),labels=c("0","1+"))
      )
    yedcol=c("ALCOHOL_SUB_MISUSE_ADMISSIONS", "BNF_ANTICOAGULENT", "BNF_ANTIFIBRINOLYTIC", 
      "BNF_ANTISECRETORY", "BNF_CNS", "BNF_CORTICOSTEROIDS", "BNF_DIABETES", 
      "BNF_FLUIDS", "BNF_GUT_MOTILITY", "BNF_INTESTINAL", "BNF_MUCOLYTICS", 
      "BNF_SECTIONS", "BNF_STOMA", "BNF_VITAMINS", "DEMENTIA", "ED_ATTENDANCES", 
      "ELECTIVE_DC_ADMISSIONS", "EMERGENCY_ADMISSIONS", "EMERGENCY_BEDDAYS", 
      "LTC", "OP_APPOINTMENTS", "OP_PSYCHIATRIC_APPOINTMENTS", "PSYCHIATRIC_ADMISSIONS")
    df=df[,yedcol]
    
  }
  if (mode=="u16") {
    df=df %>%
      mutate(
        
        # Demographics
        AGE = cut(age, breaks = c(-Inf, seq(0,14), Inf),labels=c(0:14,"15+")),
        GENDER=factor(2-sexM),
        # SIMD_QUINTILE=factor(floor((as.numeric(SIMD_DECILE_2016_SCT)+1)/2)),
        
        # Attendances
        ED_ATTENDANCES = cut(num_ae2_attendances, breaks = c(-Inf, 0:3, Inf),labels=c(0:3,"4+")),
        EMERGENCY_ADMISSIONS = cut(num_emergency_admissions, breaks = c(-Inf, 0:6, Inf),labels=c(0:6,"7+")), # This is the proper version, but it has 3-year lookback
        # EMERGENCY_BEDDAYS = cut(emergency_bed_days, breaks = c(-Inf, 0:21, Inf),labels=c(0:21,"22+")),
        ELECTIVE_ADMISSIONS = cut(num_elective_admissions, breaks = c(-Inf, 0:2, Inf),labels=c(0:2,"3+")),
        #ELECTIVE_BEDDAYS = cut(elective_bed_days, breaks = c(-Inf, 0:28, Inf),labels=c(0:28,"29+")),
        OP_APPOINTMENTS = cut(num_outpatient_appointment_general, breaks = c(-Inf, 0:3, Inf),labels=c(0:3,"4+")),
        # OP_PSYCHIATRIC_APPOINTMENTS = cut(num_outpatient_appointment_mental,c(-Inf,0,Inf),labels=c("0","1+")),
        # DC_ADMISSIONS=cut(num_dc_elective_admissions+num_dc_emergency_admissions,c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        # ELECTIVE_DC_ADMISSIONS=cut(num_dc_elective_admissions,c(-Inf,0:6,Inf),labels=c(0:6,"7+")),
        # ALCOHOL_SUB_MISUSE_ADMISSIONS = cut(emergency_drugAndalcohol_admin,c(-Inf,0:1,Inf),labels=c(0:1,"2+")),
        # ALCOHOL_ADMISSIONS =  cut(alcohol_admin, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # PSYCHIATRIC_ADMISSIONS = cut(num_psych_admissions,c(-Inf,0,Inf),labels=c("0","1+")),
        
        # LTCs
        # LTC=cut(numLTCs_resulting_in_admin, breaks=c(-Inf,0:2,Inf),labels=c(0:2,"3+")),
        # PARKINSONS_DISEASE=cut(parkinsons_indicated,breaks=c(-Inf,0,Inf),labels=c(0:1)),
        # MS=cut(MS_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        EPILEPSY=cut(epilepsy_indicated_u16,breaks=c(-Inf,0,Inf),labels=c(0:1)),
        # DEMENTIA=cut(dementia_indicated,breaks=c(-Inf,0:1,Inf),labels=c(0:2)),
        ASTHMA=cut(asthma_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        DIABETES=cut(diabetes_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        CANCER=cut(cancer_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        ENDOCRINE=cut(endocrine_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        DIGESTIVE=cut(digestive_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        BLOOD=cut(blood_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        CONGENITAL=cut(congenital_indicated,breaks=c(-Inf,0,Inf),labels=0:1),
        SELF_HARM=cut(selfharm,breaks=c(-Inf,0,Inf),labels=0:1),
        
        # BNF
        BNF_SECTIONS = cut(num_bnf_sections_u16, breaks = c(-Inf, 0:8, Inf),labels=c(0:8,"9+")),
        # BNF_RESPIRATORY = cut(pis_Respiratory, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        #BNF_CNS = cut(pis_CentralNervousSystem, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_INFECTIONS = cut(pis_Infections, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_ENDOCRINE = cut(pis_EndocrineSystem, breaks = c(-Inf, 0:7, Inf),labels=c(0:7,"8+")),
        # BNF_SUBSTANCE = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_DEMENTIA = cut(pis_Dementia, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_CORTICOSTEROIDS = cut(pis_corticosteroids, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_FLUIDS = cut(pis_Fluids_And_Electrolytes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_NUTRITION = cut(pis_Oral_Nutrition, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_VITAMINS = cut(pis_Vitamins, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_BANDAGES = cut(pis_Arm_Sling_Bandages, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_CATHETERS = cut(pis_Catheters, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_STOMA = cut(pis_CombinedStomaDevices, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        # BNF_INCONTINENCE = cut(pis_IncontinenceDevices, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_GUT_MOTILITY = cut(pis_Antispasmod_Other_Drgs_Alt_Gut_Motility, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTISECRETORY = cut(pis_Antisecretory_Drugs_Mucosal_Protectants, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_INTESTINAL = cut(pis_Drugs_Affecting_Intestinal_Secretions, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTICOAGULENT = cut(pis_Anticoagulants_And_Protamine, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_ANTIFIBRINOLYTIC = cut(`pis_Antifibrinolytic_Drugs _Haemostatics`, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_MUCOLYTICS = cut(pis_Mucolytics, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_DIABETES = cut(pis_Drugs_Used_In_Diabetes, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_GASTRO_INT = cut(pis_gastro_int, breaks = c(-Inf,0:9,Inf),labels=c(0:9,"10+")),
        BNF_HEART = cut(pis_heart, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_BRONCHO = cut(pis_bronco, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_CROMO = cut(pis_cromo, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_SUB_DEPEND = cut(pis_sub_depend, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_ANTIBACTERIAL = cut(pis_antibacterial, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_RHEUMATIC = cut(pis_rheumatic, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        #BNF_NEUROMUSCULAR = cut(pis_Drugs_Used_In_Substance_Dependence, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_NEUROMUSCULAR = cut(pis_neuromuscular, breaks = c(-Inf,0,Inf),labels=c("0","1+")),
        BNF_MYDRIATICS = cut(pis_mydriatics, breaks = c(-Inf,0,Inf),labels=c("0","1+"))
      )
    u16col=c("AGE", "ASTHMA", "BLOOD", "BNF_ANTIBACTERIAL", "BNF_BRONCHO", 
      "BNF_CORTICOSTEROIDS", "BNF_CROMO", "BNF_GASTRO_INT", "BNF_HEART", 
      "BNF_MYDRIATICS", "BNF_NEUROMUSCULAR", "BNF_RHEUMATIC", "BNF_SECTIONS", 
      "BNF_SUB_DEPEND", "CANCER", "CONGENITAL", "DIABETES", "DIGESTIVE", 
      "ED_ATTENDANCES", "ELECTIVE_ADMISSIONS", "EMERGENCY_ADMISSIONS", 
      "ENDOCRINE", "EPILEPSY", "GENDER", "OP_APPOINTMENTS", "SELF_HARM")
    df=df[,u16col]
    
  }
  
  df=df[,setdiff(colnames(df),colnames(df0))]
  return(df)
}




