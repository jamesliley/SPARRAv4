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
  #output_of_sparsemat <- primitive_all_codes_docterm_sparsematrix(
  #  patients, episodes, list_of_data_tables,
  #  all_unique_codes = colnames(posterior(topic_model_fit)$terms), # This gets all the unique codes of the topic model fit
  #  return_extra_codes = return_extra_codes,
  #  table_of_all_codes_long = table_of_all_codes_long,
  #  min_total_codes = 1,
  #  min_unique_codes = 1
  #)
  
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




#' Returns the number of ICD10 diagnosis codes for each patient
#' @param source_table_names (optional, default is all) Only look for ICD10 diagnosis codes in the given data tables
#' @param unique_vals (optional, default = FALSE) 
#' If False, Returns the total number of condition codes, 
#' if True, returns the number of unique condition codes
#'
#' @export
transformer_diagnosis_codes_count <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  source_table_names = c("SMR00", "SMR01", "SMR01E", "SMR04"),
  unique_vals = FALSE
){
  
  table_of_diag_codes <- primitive_diagnosis_codes_count(
    patients, episodes, list_of_data_tables, 
    source_table_names
    )
  
  if (unique_vals){
    out_matrix <- table_of_diag_codes %>% select(id, count_unique_condition_codes )
  } else {
    out_matrix <- table_of_diag_codes %>% select(id, count_total_condition_codes )
  }

  
  out_missing_default = as.list(rep(0, length=(ncol(out_matrix)-1)))
  names(out_missing_default) <- names(out_matrix)[-1]

  # Return in the required transformer output
  list(
    matrix = out_matrix,
    missing_default = out_missing_default
  )
  
}




#' Returns posterior assignment probabilities of a patient belonging to each of the diagnosis clusters (=topics)
#' @param topic_model_fit a fitted topicmodels::LDA() object, trained on ICD10 codes representing medical history summaries
#' @param return_extra_codes binary, if TRUE, also return a column with the number of total and number of unique codes that are unrecognised by the trained topic model
#' 
#' 
transformer_diagnosis_codes_topic <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  topic_model_fit,
  
  # Further optional parameters (these have defaults)
  return_extra_codes = TRUE,
  table_of_diag_codes_long = NULL # One may supply a pre-computed table_of_diag_codes in the long format (to save 25 mins), 
  #  this should NOT be done?? unless the table_of_diag_codes_long object was created with the same time filtering than the rest of the inputs?
){
  
  
  
  # Get a doc_term matrix to infer on
  output_of_sparsemat <- primitive_diagnosis_codes_docterm_sparsematrix(
    patients, episodes, list_of_data_tables,
    all_unique_codes = colnames(posterior(topic_model_fit)$terms), # This gets all the unique codes of the topic model fit
    return_extra_codes = return_extra_codes,
    table_of_diag_codes_long = table_of_diag_codes_long,
    min_total_codes = 1,
    min_unique_codes = 1
  )
  
  print("Diagnosis code matrix generated")
  
  if (return_extra_codes == TRUE){
    doc_term_sparse_matrix <- output_of_sparsemat$doc_term_sparse_matrix
  } else {
    doc_term_sparse_matrix <- output_of_sparsemat
  }
  
  # Do the inference given the fit and the new doc_term matrix
  lda_post <- posterior(
    topic_model_fit, 
    newdata = doc_term_sparse_matrix)
  
  print("Posterior calculated")
  
    
  # Create the inference matrix 
  out_matrix <- as_tibble(
      lda_post$topics,
      rownames = "id"
    ) %>% 
    rename_at(.vars = 2:ncol(.), ~ paste0("topic_", .)) %>%
    mutate(
      id = as.integer(id)
    )
  
  print("Inference matrix generated")
  
  # Join extra columns with some statistics of missing/extra codes
  if (return_extra_codes == TRUE){
    out_matrix <- out_matrix %>%
      full_join(
        # Create a summary of extra codes by id
        output_of_sparsemat$table_of_extra_codes %>% 
          group_by(id) %>%
          summarise(
            topics_missing_code_total = sum(value),
            topics_missing_code_unique = n()
          )
      )
  }
  
  out_missing_default = as.list(rep(0, length=(ncol(out_matrix)-1)))
  names(out_missing_default) <- names(out_matrix)[-1]
  
  print("Post-processing completed")
  
  # Return in the required transformer output
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
  source_table_names = c("SMR01", "SMR01E")
  
){
  
  print("Check (last emergency admission days) 1/3")
  
  # Filter the appropriate tables
  data_inpatient <- list_of_data_tables[c("SMR01", "SMR01E")]
  data_inpatient_emergency <- lapply(
    data_inpatient,
    function(df) df %>% 
      select(id, time, source_table, source_row, admission_type, ADMISSION_TRANSFER_FROM ) %>%
      filter(
      (admission_type %in% 20:22 | admission_type %in% 30:36 | admission_type==39) |
        (admission_type==38 & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
    )
    ) %>%
    bind_rows()
  
  print("Check (last emergency admission days) 2/3")
  
  # Call the helper function "transformer_time_since ..."
  out <- transformer_last_episode_days_ago(
    patients = NULL,
    episodes = episodes %>%
      select(id, time, source_table, source_row) %>%
      semi_join(data_inpatient_emergency),
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
  # filter_conditions_quosure=quo(TRUE), # for example for loosely emergency admissions filter_conditions_quosure = quo(!is.na(admission_type) & admission_type < 40 & admission_type >= 30)
  # filter_conditions_name_added = ""
){
  
  print("Check (last episode days ago) 1/3")
  
  # Get the days since last visit for each table in "source_table_names"
  out_matrix <- episodes %>%
    ungroup() %>%
    mutate(time = as.integer(time)) %>% 
    select(id, time, source_table) %>%
    filter(source_table %in% source_table_names) %>% # Get rid of unnecessary information (conditional filtering should be done prior to calling this helper funciton)
    mutate(source_table = paste0("days_since_last_", source_table)) %>%
    mutate(source_table = factor(source_table, levels = paste0("days_since_last_", source_table_names))) %>% # ensure that returns as many columns as entries in source_table_names
    group_by(id, source_table) %>%
    summarise(last_time = max(time)) %>% 
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
#'
#' @export
transformer_ltcs <- function(
  # Input data
  patients, episodes, list_of_data_tables,
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further required parameters
  output_type # %in% c("total_count", "years_since_diag", "binary", "rawdata_NUMBEROFLTCs")
  
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
  
  ltc_names <- c("FIRST_ARTHRITIS_EPISODE", "FIRST_ASTHMA_EPISODE", "FIRST_ATRIAL_FIBRILLATION_EPISODE",
                 "FIRST_CANCER_EPISODE", "FIRST_CHRONIC_LIVER_DISEASE_EPISODE", 
                 "FIRST_COPD_EPISODE", "FIRST_DEMENTIA_EPISODE", "FIRST_DIABETES_EPISODE", 
                 "FIRST_EPILEPSY_EPISODE", "FIRST_HEART_DISEASE_EPISODE", "FIRST_HEART_FAILURE_EPISODE", 
                 "FIRST_MULTIPLE_SCLEROSIS_EPISODE", "FIRST_PARKINSON_DISEASE_EPISODE",
                 "FIRST_RENAL_FAILURE_EPISODE", "FIRST_CEREBROVASCULAR_DISEASE_EPISODE"
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
      mutate(value = as.integer(as.integer(time)/(60*60*24*365))) %>% # Divide here to avoid overflow issues
      select(-time) %>%
      mutate(value = sapply(as.integer((time_cutoff/(60*60*24*365)) - value), function(x)max(x, as.integer(1)))) %>% # Subtract for time_cutoff (in years)
      spread(key = LTC_TYPE, value = value, fill=as.integer(0), drop=FALSE)
    
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




#' transformer_sparrav3_grouped_inpatient_admissions
#' This function counts the number of 
#'    a) emergency admissions with drug or alcohol diagnosis
#'    b1) either (sum_over_LTCs == TRUE ) the number of different LTCs that has resulted in an admission (question: emergency or normal is also ok, ambiguous)
#'    b2) or (sum_over_LTCs == FALSE ) the number of admissions each LTC has resulted in
#'    as seen in the SPARRAv3 report, Table 2 / "Hospital Inpatient Admissions"
#'    
#' The intended grouping of the drugs may be found in the uploaded excel file:
#'    DRUG AND ALCOHOL DIAGNOSIS
#'    "\\Farr-FS1\Study Data\1718-0370\Linked Data\Upload to SH 4\admissions and prescribing (for Gergo) RP.xlsx"
#'    on Sheet 1
#'    
#'    GROUPING OF DIAGNOSIS CODES INTO LTCS
#'    "\\Farr-FS1\Study Data\1718-0370\Linked Data\Upload to SH 4\Long Term Conditions ICD9 and ICD10 codes (for Gergo) RP.xlsx"
#'  
#'    
#' Inputs
#' @param patients A tibble of patients with related basic information
#' @param episodes A tibble of individual records with related basic information
#' @param list_of_data_tables A list of the individual data tables with all detailed columns from the raw data
#' 
#' # Extra inputs
#' @param sum_over_LTCs If TRUE, the number of different LTCs that led to emergency admissions are output. Otherwise for each LTC we get a binary feature whether or not it has led to emergency admission.
#' 
#' Output
#' @return A named list with two elements: \itemize{
#'    \item{matrix}{a PATIENTS (rows) x FEATURES (columns) matrix, one patient per row, the two columns represent 
#'                        "number of alcohol/drug related admissions" and "number of different LTCs which have resulted in admission" }
#'    \item{missing_default}{a FEATURES (length) list of default value (0 here) per feature to set for all patients that have that particular feature as missing}
#' }
#' 
#' Written by Sam Oduro
#' 
#' @export
transformer_sparrav3_grouped_inpatient_admissions <- function(
  patients, episodes, list_of_data_tables,
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Extra variables
  sum_over_LTCs = TRUE
){
  
  # One needs to filter the "SMR01" table to only keep episodes with relevant diagnosis codes (see the groupings in the above cited excel files)
  #   (for filtering examples see the existing transformer_sparra_v3 code where we seperate data into emergency/elective admissions or mental/non-mental ones)
  # The appropriate prescription sections or diagnosis codes for filtering may be found in the above cited excel file
  
  # Then the number of filtered episodes per person per category may be counted via the primitive transformer 
  # "transformer_count()". 
  # Note that for LTCs the diagnosis codes first need to be grouped to LTCs, and then transformer_count should be called unique_vals=TRUE
  
  # Rename the columns into meaningful features
  
  # Create the final matrix
  
  # Output a list containing the final matrix and the defaults (0 for both categories)
  
  # Source the manually created ltc lookup file
  # source("Gergo/Lookups_and_groupings/get_icd10_grouping_ltc.R")
  # source("Gergo/Lookups_and_groupings/get_icd10_grouping_drug_alcohol.R")
  
  ltc_groupings <- get_icd10_grouping_ltc()
  # Add "ltc_" prefix to names and "_admin" postfix
  names(ltc_groupings) <- sapply(names(ltc_groupings), function(x)str_c(str_c("ltc_", x), "_admin"))
  
  alcohol_drug_groupings <- get_icd10_grouping_drug_alcohol()
  
  
  
  
  data_inpatient <- list_of_data_tables[c("SMR01", "SMR01E")]
  #table(data_inpatient$SMR01E$admission_type)
  
  data_inpatient_emergency <- lapply(
    data_inpatient,
    function(df) {
      # Add new ltc columns
      for (ltc_name in names(ltc_groupings)){
        df <- df %>% mutate(
          !!ltc_name := ifelse((main_condition %in% ltc_groupings[[ltc_name]]|OTHER_CONDITION_1 %in% ltc_groupings[[ltc_name]]|
                               OTHER_CONDITION_2 %in% ltc_groupings[[ltc_name]]|OTHER_CONDITION_3 %in% ltc_groupings[[ltc_name]]|
                               OTHER_CONDITION_4 %in% ltc_groupings[[ltc_name]]|OTHER_CONDITION_5 %in% ltc_groupings[[ltc_name]]),1,0)
        )
      }
      
      # Add new alcohol, drug and emergency admission columns
      df <- df %>% 
        mutate(
          alcohol_admin=ifelse((main_condition %in% alcohol_drug_groupings$alcohol|OTHER_CONDITION_1 %in% alcohol_drug_groupings$alcohol|
                               OTHER_CONDITION_2 %in% alcohol_drug_groupings$alcohol|OTHER_CONDITION_3 %in% alcohol_drug_groupings$alcohol|
                               OTHER_CONDITION_4 %in% alcohol_drug_groupings$alcohol|OTHER_CONDITION_5 %in% alcohol_drug_groupings$alcohol),1,0),
          drug_admin=ifelse((main_condition %in% alcohol_drug_groupings$drug|OTHER_CONDITION_1 %in% alcohol_drug_groupings$drug|
                               OTHER_CONDITION_2 %in% alcohol_drug_groupings$drug|OTHER_CONDITION_3 %in% alcohol_drug_groupings$drug|
                               OTHER_CONDITION_4 %in% alcohol_drug_groupings$drug|OTHER_CONDITION_5 %in% alcohol_drug_groupings$drug),1,0),
          emergency_admin=ifelse(
            ((admission_type %in% 20:22 | admission_type %in% 30:36| admission_type==39) |
               (admission_type==38 & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))),1,0)
        ) %>% 
        # Mutate the emergency drug and alcohol admissions together
        mutate(emergency_drugAndalcohol_admin=ifelse(((alcohol_admin==1|drug_admin==1) & emergency_admin==1),1,0)) %>%
        select(-c(emergency_admin, alcohol_admin, drug_admin))
    }
  )
      
             
  
  colNames<- names(data_inpatient_emergency$SMR01[,
                                                  grepl("_admin",names(data_inpatient_emergency$SMR01))])
  
  # colNamesSMR01<-paste0("SMR01.",colNames)
  # colNamesSMR01E<-paste0("SMR01E.",colNames)
  
  # Emergency alcohol & drugs and LTC admission count
  out <- transformer_sum(
    patients, episodes,
    list_of_data_tables = data_inpatient_emergency,
    time_cutoff = NULL,
    column_names_flat = c(paste0("SMR01.",colNames),paste0("SMR01E.",colNames))
  )
  
  # Sum SMR01 and SMR01 within same categories
  patterns <- setdiff(unique(str_extract(names(out$matrix),"[^.]+$")),"id")
  out$matrix<-  out$matrix %>% replace(.,is.na(.),0)
  out$matrix<-as_tibble(cbind(id=out$matrix$id,
                              sapply(patterns,function(x) 
                                rowSums(out$matrix[,grep(x, names(out$matrix)),drop=F]))))
  

  # We either want to collect the number of ltcs that has resulted in admission (for a single feature - see default sparra) OR the 
  if (sum_over_LTCs){
    out$matrix <- out$matrix %>%
      mutate_at(
        vars(starts_with("ltc_")),
        function(x)x>0
      ) %>%
      transmute(
        id = id, 
        emergency_drugAndalcohol_admin = emergency_drugAndalcohol_admin,
        numLTCs_resulting_in_admin = rowSums(select(., starts_with("ltc_")))
      )
  }
  
  # Set missing defaults
  out$missing_default = as.list(rep(0, length=length(names(out$matrix)[-1])))
  names(out$missing_default) <- names(out$matrix)[-1]
  
  # Output a list containing the final matrix and the defaults (0 for each category)
  out
}# ------------------------------------------------------
# Function 1

#' transformer_sparrav3_grouped_prescriptions
#' This function counts the number of prescriptions for a specific drugs or groups of drugs
#'    as seen in the SPARRAv3 report, Table 2 / "Prescriptions for specific drugs or groups of drugs"
#'    
#' The intended grouping of the drugs may be found in the uploaded excel file:
#'    "\\Farr-FS1\Study Data\1718-0370\Linked Data\Upload to SH 4\admissions and prescribing (for Gergo) RP.xlsx"
#'    on Sheet 6
#'  
#'    
#' Inputs
#' @param patients A tibble of patients with related basic information
#' @param episodes A tibble of individual records with related basic information
#' @param list_of_data_tables A list of the individual data tables with all detailed columns from the raw data
#' 
#' @param only_sparrav3 If TRUE, keep only prescription information that was in SPARRAv3. If false we keep all features (named appropriately)
#' 
#' Output
#' @return A named list with two elements: \itemize{
#'    \item{matrix}{a PATIENTS (rows) x FEATURES (columns) matrix, one patient per row, the columns are the appropriate BNF groupings}
#'    \item{missing_default}{a FEATURES (length) list of default value (0 here) per feature to set for all patients that have that particular feature as missing}
#' }
#' 
#' Written by Sam Oduro
#' Update by Gergo Bohner 26/02/2019 due to different input format of PIS table (long instead of wide)
#' 
#' @export
transformer_sparrav3_grouped_prescriptions <- function(
  patients, episodes, list_of_data_tables,
  time_cutoff, # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Extra variables
  only_sparrav3 = TRUE
){
  
  
  
  #rm(list=ls())
  
  
  
  data_PIS<-list_of_data_tables$PIS
  # feasKeep<-c("id","source_table","source_row","time", "PAID_GIC_INCL_BB","TWELVE_MONTH_TIME_PERIOD","NUMBER_OF_PAID_ITEMS" )
  # colNames<- paste0("PIS.",setdiff(names(data_PIS),feasKeep))
  
  bnf_lookup<-get_bnf_lookup() #data.frame(readxl::read_excel("Sam/All_Available_BNF_Sections.xlsx"))
  
  
  
  
  # Set up a list of BNF sections to retain/combine
  combine_BNF_list <- list(
    Respiratory = c("NUM_BNF_0301", "NUM_BNF_0302", 
                    "NUM_BNF_0303", "NUM_BNF_0304", "NUM_BNF_0305", "NUM_BNF_0306",
                    "NUM_BNF_0307", "NUM_BNF_0308", "NUM_BNF_0309", "NUM_BNF_0310"),
    CentralNervousSystem = c("NUM_BNF_0401", "NUM_BNF_0402", "NUM_BNF_0403",
                            "NUM_BNF_0404", "NUM_BNF_0405", "NUM_BNF_0406", "NUM_BNF_0407",
                            "NUM_BNF_0408", "NUM_BNF_0409", "NUM_BNF_0410", "NUM_BNF_0411"),
    Infections = c("NUM_BNF_0501", "NUM_BNF_0502", "NUM_BNF_0503", "NUM_BNF_0504",
                       "NUM_BNF_0505"),
    EndocrineSystem = c( "NUM_BNF_0601", "NUM_BNF_0602", "NUM_BNF_0603",
                         "NUM_BNF_0604", "NUM_BNF_0605", "NUM_BNF_0606", "NUM_BNF_0607"),
    IncontinenceDevices = c( "NUM_BNF_2201", "NUM_BNF_2202",
                             "NUM_BNF_2205", "NUM_BNF_2210", "NUM_BNF_2215", "NUM_BNF_2220",
                             "NUM_BNF_2230", "NUM_BNF_2240", "NUM_BNF_2250", "NUM_BNF_2260",
                             "NUM_BNF_2270", "NUM_BNF_2280", "NUM_BNF_2285", "NUM_BNF_2290"),
    CombinedStomaDevices = c( "NUM_BNF_2305", "NUM_BNF_2310", "NUM_BNF_2315", "NUM_BNF_2320",
                             "NUM_BNF_2325", "NUM_BNF_2330", "NUM_BNF_2335", "NUM_BNF_2340",
                             "NUM_BNF_2345", "NUM_BNF_2346", "NUM_BNF_2350", "NUM_BNF_2355",
                             "NUM_BNF_2360", "NUM_BNF_2365", "NUM_BNF_2370", "NUM_BNF_2375",
                             "NUM_BNF_2380", "NUM_BNF_2385", "NUM_BNF_2390", "NUM_BNF_2392",
                             "NUM_BNF_2393", "NUM_BNF_2394", "NUM_BNF_2396", "NUM_BNF_2398")
  )
  
  BNF_RetainAsSingleSection = c('NUM_BNF_0102','NUM_BNF_0103','NUM_BNF_0109','NUM_BNF_0208','NUM_BNF_0211',
                            'NUM_BNF_0302','NUM_BNF_0307','NUM_BNF_0408','NUM_BNF_0409','NUM_BNF_0410',
                            'NUM_BNF_0411','NUM_BNF_0601','NUM_BNF_0902','NUM_BNF_0904','NUM_BNF_0905',
                            'NUM_BNF_0906','NUM_BNF_1002','NUM_BNF_2002','NUM_BNF_2102')
  
  rename_lookup <- bnf_lookup %>% 
    select(BNF_Section2, Description1) %>%
    rename(BNF_section:=BNF_Section2) %>%
    filter(BNF_section %in% BNF_RetainAsSingleSection) %>%
    mutate(
      BNF_section = factor(BNF_section, levels = levels(data_PIS$BNF_section)),
      Description1 = factor(Description1)
      )
  
  library(forcats) # for combining factor levels
  
  # base information tables
  d1 = data_PIS %>% 
    select(id, YEARMONTH, PAID_GIC_INCL_BB, NUMBER_OF_PAID_ITEMS) %>%
    distinct() %>%
    select(-YEARMONTH)
  dp=as.tibble(aggregate(d1[,2:3],by=list(id=d1$id),FUN=sum));
  
  d2 = data_PIS %>%
    select(id, BNF_section) 
  dp2=as.tibble(aggregate(d2[,2],by=list(id=d2$id),FUN=n_distinct))
  colnames(dp2)[2]="countBNFsections"
  
  # combined sections table
  d3=data_PIS %>% 
    filter(BNF_section %in% unlist(combine_BNF_list)) %>%
    mutate(
      Description1 = do.call(fct_collapse, c(list(.f=.$BNF_section), combine_BNF_list))
    ) %>%
    select(id, Description1, NUM_sold)
  dp3=as_tibble(aggregate(d3[,3],by=list(id=d3$id,Description1=d3$Description1),FUN=sum,na.rm=TRUE))
  dp3= dp3 %>% spread(key=Description1, value=NUM_sold)
  
  # single sections table
  d4=data_PIS %>% 
    filter(BNF_section %in% BNF_RetainAsSingleSection) %>% 
    left_join(rename_lookup, by=c("BNF_section")) %>%
    select(id, Description1, NUM_sold)
  dp4=as_tibble(aggregate(d4[,3],by=list(id=d4$id,Description1=d4$Description1),FUN=sum,na.rm=TRUE))	
  dp4= dp4 %>% spread(key=Description1, value=NUM_sold)
  
  
  out_matrix <- (
    as_tibble(dp) %>% full_join(dp2,by = c("id")) %>% 
      full_join(dp3,by = c("id")) %>% full_join(dp4,by = c("id"))
  )
  
  colnames(out_matrix) = c("id", sapply(colnames(out_matrix)[-1], function(x)paste0("pis_", x)))
  
  out = list(
    matrix = out_matrix,
    missing_default = as.list(rep(0, length=length(names(out_matrix)[-1])))
  )

  names(out$missing_default) <- names(out$matrix)[-1]

  gc()

  # Return
  out
  
  
  
  
  # Set up for parallel map processing
  # library(furrr)
  # plan(multiprocess)
  # 
  # tic()
  # out_matrix <- data_PIS %>% 
  #   # filter(id < 100000) %>% # Just here for testing
  #   group_by(id) %>%
  #   nest() %>%
  #   mutate(
  #     tmp_col = future_map(data, 
  #                   function(x){tibble(
  #                     allBNFSections= length(unique(x$BNF_section)),
  #                     Respiratory = sum(x %>% filter(BNF_section %in% c("NUM_BNF_0301", "NUM_BNF_0302", 
  #                                                                       "NUM_BNF_0303", "NUM_BNF_0304", "NUM_BNF_0305", "NUM_BNF_0306", 
  #                                                                       "NUM_BNF_0307", "NUM_BNF_0308", "NUM_BNF_0309", "NUM_BNF_0310")) %>% pull(NUM_sold), na.rm=TRUE),
  #                     CentralNervousSystem = sum(x %>% filter(BNF_section %in% c("NUM_BNF_0401", "NUM_BNF_0402", "NUM_BNF_0403", 
  #                                                                       "NUM_BNF_0404", "NUM_BNF_0405", "NUM_BNF_0406", "NUM_BNF_0407", 
  #                                                                       "NUM_BNF_0408", "NUM_BNF_0409", "NUM_BNF_0410", "NUM_BNF_0411")) %>% pull(NUM_sold), na.rm=TRUE),
  #                     Infections = sum(x %>% filter(BNF_section %in% c("NUM_BNF_0501", "NUM_BNF_0502", "NUM_BNF_0503", "NUM_BNF_0504", 
  #                                                                      "NUM_BNF_0505")) %>% pull(NUM_sold), na.rm=TRUE),
  #                     EndocrineSystem = sum(x %>% filter(BNF_section %in% c( "NUM_BNF_0601", "NUM_BNF_0602", "NUM_BNF_0603", 
  #                                                                        "NUM_BNF_0604", "NUM_BNF_0605", "NUM_BNF_0606", "NUM_BNF_0607")) %>% pull(NUM_sold), na.rm=TRUE),
  #                     IncontinenceDevices = sum(x %>% filter(BNF_section %in% c( "NUM_BNF_2201", "NUM_BNF_2202", 
  #                                                                                "NUM_BNF_2205", "NUM_BNF_2210", "NUM_BNF_2215", "NUM_BNF_2220", 
  #                                                                                "NUM_BNF_2230", "NUM_BNF_2240", "NUM_BNF_2250", "NUM_BNF_2260", 
  #                                                                                "NUM_BNF_2270", "NUM_BNF_2280", "NUM_BNF_2285", "NUM_BNF_2290")) %>% pull(NUM_sold), na.rm=TRUE),
  #                     CombinedStomaDevices = sum(x %>% filter(BNF_section %in% c( "NUM_BNF_2305", "NUM_BNF_2310", "NUM_BNF_2315", "NUM_BNF_2320", 
  #                                                                                "NUM_BNF_2325", "NUM_BNF_2330", "NUM_BNF_2335", "NUM_BNF_2340", 
  #                                                                                "NUM_BNF_2345", "NUM_BNF_2346", "NUM_BNF_2350", "NUM_BNF_2355", 
  #                                                                                "NUM_BNF_2360", "NUM_BNF_2365", "NUM_BNF_2370", "NUM_BNF_2375", 
  #                                                                                "NUM_BNF_2380", "NUM_BNF_2385", "NUM_BNF_2390", "NUM_BNF_2392", 
  #                                                                                "NUM_BNF_2393", "NUM_BNF_2394", "NUM_BNF_2396", "NUM_BNF_2398")) %>% pull(NUM_sold), na.rm=TRUE),
  #                                           
  #                     
  #                     ) %>% # Add the "interesting" single sections with appropriate names
  #                     bind_cols(
  #                       x %>% 
  #                         filter(
  #                         BNF_section %in% c('NUM_BNF_0102','NUM_BNF_0103','NUM_BNF_0109','NUM_BNF_0208','NUM_BNF_0211',
  #                                                       'NUM_BNF_0302','NUM_BNF_0307','NUM_BNF_0408','NUM_BNF_0409','NUM_BNF_0410',
  #                                                       'NUM_BNF_0411','NUM_BNF_0601','NUM_BNF_0902','NUM_BNF_0904','NUM_BNF_0905',
  #                                                       'NUM_BNF_0906','NUM_BNF_1002','NUM_BNF_2002','NUM_BNF_2102')
  #                         ) %>% 
  #                         left_join(
  #                           bnf_lookup %>% 
  #                             select(BNF_Section2, Description1) %>%
  #                             rename(BNF_section:=BNF_Section2)
  #                         ) %>% 
  #                         select(Description1, NUM_sold) %>%
  #                         group_by(Description1) %>% 
  #                         summarise(NUM_sold = sum(NUM_sold, na.rm=TRUE)) %>%
  #                         spread(key=Description1, value=NUM_sold) %>%
  #                         mutate_all(as.integer)
  #                     )
  #               })
  #   ) %>%
  #   select(-data) %>%
  #   unnest()
  # toc()
  # 
  # colnames(out_matrix) = c("id", sapply(colnames(out_matrix)[-1], function(x)paste0("pis_", x)))
  # 
  # out = list(
  #   matrix = out_matrix,
  #   missing_default = as.list(rep(0, length=length(names(out_matrix)[-1])))
  # )
  # 
  # names(out$missing_default) <- names(out$matrix)[-1]
  # 
  # gc()
  # 
  # # Return
  # out
  # 
}

#' transformer_sparrav3
#' This transformer attempts to create features similar to those used by SPARRAv3
#' 
#' Usage example
#' transformer_sparrav3(patients, episodes, list_of_data_tables, time_cutoff)
#' 
transformer_sparrav3 <- function(
  # Input data
  patients, episodes, list_of_data_tables,
  time_cutoff # time_cutoff is the latest allowed time for an episode as a unix integer time (seconds since 1970)
  
  # Further optional parameters (these have defaults)
){
  
  # Save actual input for debugging reasons
  #save("patients", "episodes", "list_of_data_tables", file="Gergo/Transformation/Outputs/TMP_INPUT_To_SPARRA_TRANSFORMER.RData")
  
  out <- NULL
  join_transformer_outputs <- function(out, new_inp){
    if (is.null(out)) {out <- new_inp
    } else {
      out$matrix <- out$matrix %>%
        full_join(new_inp$matrix, by="id")
      out$missing_default = c(out$missing_default, new_inp$missing_default)
    }
    out
  }
  
  print("Check (v3) 1/12")
  
  # Define the columns we need to use for each feature
  
  # -----------------------------------
  # HOSPITAL INPATIENT ADMISSIONS
  # -----------------------------------
  
  
  # Split to emergency and non-emergency admission
  # (based on more complex definition of emergency admission)
  data_inpatient <- list_of_data_tables[c("SMR01", "SMR01E")]
  data_inpatient_emergency <- lapply(
    data_inpatient,
    function(df) df %>% filter(
        (admission_type %in% 20:22 | admission_type %in% 30:36 | admission_type==39) |
           (admission_type==38 & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
    )
  )
  
  
  data_inpatient_elective <- lapply(
    data_inpatient,
    function(df) df %>% filter(
      !((admission_type %in% 20:22 | admission_type %in% 30:36 | admission_type==39) |
        (admission_type==38 & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
      )
    )
  )
  
  print("Check (v3) 2/12")
  
  
    
  # Emergency bed days
  tmp <- transformer_sum(
    patients, episodes,
    list_of_data_tables = data_inpatient_emergency,
    time_cutoff = NULL,
    column_names_flat = c("SMR01.LENGTH_OF_STAY", "SMR01E.LENGTH_OF_STAY")
  )
  tmp$matrix <- tmp$matrix %>% 
    replace( is.na(.), 0 ) %>% # This is needed as the + in the following line doesn't ignore NA!
    transmute(id = id, emergency_bed_days = SMR01.LENGTH_OF_STAY + SMR01E.LENGTH_OF_STAY)
  tmp$missing_default <- list(emergency_bed_days = 0)
  out <- join_transformer_outputs(out, tmp)
  
  # Emergency admission count
  tmp <- transformer_count(
    patients, episodes,
    list_of_data_tables = data_inpatient_emergency,
    time_cutoff = NULL,
    column_names_flat = c("SMR01.LENGTH_OF_STAY", "SMR01E.LENGTH_OF_STAY")
  )
  tmp$matrix <- tmp$matrix %>% 
    replace( is.na(.), 0 ) %>% # This is needed as the + in the following line doesn't ignore NA!
    transmute(id = id, num_emergency_admissions = SMR01.LENGTH_OF_STAY + SMR01E.LENGTH_OF_STAY)
  tmp$missing_default <- list(num_emergency_admissions = 0)
  out <- join_transformer_outputs(out, tmp)
  
  print("Check (v3) 3/12")
  
  
    
  # elective bed days
  tmp <- transformer_sum(
    patients, episodes,
    list_of_data_tables = data_inpatient_elective,
    time_cutoff = NULL,
    column_names_flat = c("SMR01.LENGTH_OF_STAY", "SMR01E.LENGTH_OF_STAY")
  )
  tmp$matrix <- tmp$matrix %>% 
    replace( is.na(.), 0 ) %>% # This is needed as the + in the following line doesn't ignore NA!
    transmute(id = id, elective_bed_days = SMR01.LENGTH_OF_STAY + SMR01E.LENGTH_OF_STAY)
  tmp$missing_default <- list(elective_bed_days = 0)
  out <- join_transformer_outputs(out, tmp)

  print("Check (v3) 4/12")
  
  
  # elective admission count
  tmp <- transformer_count(
    patients, episodes,
    list_of_data_tables = data_inpatient_elective,
    time_cutoff = NULL,
    column_names_flat = c("SMR01.LENGTH_OF_STAY", "SMR01E.LENGTH_OF_STAY")
  )
  tmp$matrix <- tmp$matrix %>% 
    replace( is.na(.), 0 ) %>% # This is needed as the + in the following line doesn't ignore NA!
    transmute(id = id, num_elective_admissions = SMR01.LENGTH_OF_STAY + SMR01E.LENGTH_OF_STAY)
  tmp$missing_default <- list(num_elective_admissions = 0)
  out <- join_transformer_outputs(out, tmp)
  
  
  # Add emeregency admissions with drug or alcohol diagnosis and number of different LTCs which have resulted in admission
  #  source("Gergo/Transformation/transformer_functions/transformer_sparrav3_grouped_inpatient_admissions.R")
  tmp <- transformer_sparrav3_grouped_inpatient_admissions(patients, episodes, list_of_data_tables, time_cutoff = NULL, sum_over_LTCs = TRUE)
  out <- join_transformer_outputs(out, tmp)
  
  print("Check (v3) 5/12")
  
  
  # -----------------------------------
  # PRESCRIBING INFORMATION
  # -----------------------------------
  
  # Add all prescribing related features
  #  source("James/Transformation/transformer_functions/transformer_sparrav3_grouped_prescriptions.R")
  tmp <- transformer_sparrav3_grouped_prescriptions(patients, episodes, list_of_data_tables, time_cutoff = NULL, only_sparrav3 = TRUE)
  out <- join_transformer_outputs(out, tmp)

  print("Check (v3) 6/12")
  
  
  # -----------------------------------
  # AE2 ATTENDANCES & PSYCHATRIC ADMISSIONS
  # -----------------------------------
  tmp <- transformer_count(
    patients, episodes,
    list_of_data_tables = list_of_data_tables,
    time_cutoff = NULL,
    column_names_flat = c("AE2.time", "SMR04.time")
  )
  tmp$matrix <- tmp$matrix %>% transmute(id = id, num_psych_admissions = SMR04.time, num_ae2_attendances = AE2.time)
  tmp$missing_default <- list(num_psych_admissions = 0, num_ae2_attendances = 0)
  out <- join_transformer_outputs(out, tmp)
  
  print("Check (v3) 7/12")
  
  
  # -----------------------------------
  # OUTPATIENT APPOINTMENTS
  # -----------------------------------
  data_outpatient <- list_of_data_tables[c("SMR00")]
  specialty_codes_mental <- c("AV", "AH", "G1", "G2", "G21", "G22", "G3", "G4", "G5", "G6")
  
  data_outpatient_nonmental <- lapply(
    data_outpatient,
    function(df) df %>% filter(!(SPECIALTY %in% specialty_codes_mental))
  )
  data_outpatient_mental <- lapply(
    data_outpatient,
    function(df) df %>% filter(SPECIALTY %in% specialty_codes_mental)
  )

  print("Check (v3) 8/12")
  
  
  
  # Outpatient appointment count (non-mental)
  tmp <- transformer_count(
    patients, episodes,
    list_of_data_tables = data_outpatient_nonmental,
    time_cutoff = NULL,
    column_names_flat = c("SMR00.time")
  )
  tmp$matrix <- tmp$matrix %>% transmute(id = id, num_outpatient_appointment_general = SMR00.time)
  tmp$missing_default <- list(num_outpatient_appointment_general = 0)
  out <- join_transformer_outputs(out, tmp)
  
  # Outpatient appointment count (mental)
  tmp <- transformer_count(
    patients, episodes,
    list_of_data_tables = data_outpatient_mental,
    time_cutoff = NULL,
    column_names_flat = c("SMR00.time")
  )
  tmp$matrix <- tmp$matrix %>% transmute(id = id, num_outpatient_appointment_mental = SMR00.time)
  tmp$missing_default <- list(num_outpatient_appointment_mental = 0)
  out <- join_transformer_outputs(out, tmp)
  
  print("Check (v3) 9/12")
  
  
  
  # -----------------------------------
  # COMBINED DATA SOURCES
  # -----------------------------------
  
  # Add age and SIMD
  tmp <- transformer_patient_info(patients, episodes, list_of_data_tables, time_cutoff = time_cutoff)
  tmp$matrix <- tmp$matrix %>% select(-c("gender"))
  tmp$missing_default  <- tmp$missing_default[!names(tmp$missing_default) %in% "gender"]
  out <- join_transformer_outputs(out, tmp)
  
  print("Check (v3) 10/12")
  
  # Add indication of Parkinson's, MS, epilepsy or Alzheimer's/dementia
  tmp <- transformer_sparrav3_grouped_inpatient_admissions(patients, episodes, list_of_data_tables, time_cutoff = NULL, sum_over_LTCs = FALSE)
  
  tmp$matrix <- tmp$matrix %>%
    select(id, ltc_parkinsons_disease_admin, ltc_multiple_Sclerosis_admin, ltc_epilepsy_admin, ltc_dementia_admin) %>%
    full_join(
      list_of_data_tables$PIS %>%
        filter(BNF_section %in% c("NUM_BNF_0409", "NUM_BNF_1002", "NUM_BNF_0408", "NUM_BNF_0411")) %>%
        replace(.,is.na(.),0) %>%
        select(
          id, YEARMONTH, BNF_section, NUM_sold
        ) %>% 
        spread(key=BNF_section, value=NUM_sold) %>%
        select(-YEARMONTH) %>%
        group_by(id) %>%
        summarise_all(~sum(., na.rm=TRUE))
    ) %>%
    replace(.,is.na(.),0) %>%
    transmute(
      id = id,
      parkinsons_indicated = ifelse((ltc_parkinsons_disease_admin>0 | NUM_BNF_0409>0),1,0),
      MS_indicated = ifelse((ltc_multiple_Sclerosis_admin>0| NUM_BNF_1002>0),1,0),
      epilepsy_indicated =ifelse( (ltc_epilepsy_admin>0| NUM_BNF_0408>0),1,0),
      dementia_indicated =ifelse( (ltc_dementia_admin>0|NUM_BNF_0411>0),1,0)
  )

  print("Check (v3) 11/12")
  
  
  # Old version before PIS data format changes on 26 Feb 2019
  # tmp$matrix <- tmp$matrix %>%
  #   select(id, ltc_parkinsons_disease_admin, ltc_multiple_Sclerosis_admin, ltc_epilepsy_admin, ltc_dementia_admin) %>%
  #   full_join(
  #     list_of_data_tables$PIS %>% 
  #       replace(.,is.na(.),0) %>%
  #       select(
  #         id, NUM_BNF_0409, NUM_BNF_1002, NUM_BNF_0408, NUM_BNF_0411
  #       ) %>% group_by(id) %>% summarise_all(sum)
  #   ) %>% 
  #   replace(.,is.na(.),0) %>%
  #   transmute(
  #     id = id,
  #     parkinsons_indicated = ifelse((ltc_parkinsons_disease_admin>0 | NUM_BNF_0409>0),1,0),
  #     MS_indicated = ifelse((ltc_multiple_Sclerosis_admin>0| NUM_BNF_1002>0),1,0),
  #     epilepsy_indicated =ifelse( (ltc_epilepsy_admin>0| NUM_BNF_0408>0),1,0),
  #     dementia_indicated =ifelse( (ltc_dementia_admin>0|NUM_BNF_0411>0),1,0)
    # )
    
  
  tmp$missing_default = as.list(rep(0, length=length(names(tmp$matrix)[-1])))
  names(tmp$missing_default) <- names(tmp$matrix)[-1]
  
  
  out <- join_transformer_outputs(out, tmp)

  print("Check (v3) 12/12")
  
  
  # Return
  out
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


  # Define the function to apply
  # function_to_apply <- function(x){
  #   if (!is.numeric(x) & !is.null(x)) {
  #     stop("transformer_sum is attempting to sum a non-numeric column")
  #   } else {
  #     out = if_else(length(x)>0,
  #            sum(x),
  #            0)
  #   }
  #   out
  # }
  
  function_to_apply <- quo(sum(., na.rm=TRUE))

  
  # function_to_apply <- quo(
  #   if (!is.numeric(.) & !is.null(.)) {
  #     stop("transformer_sum is attempting to sum a non-numeric column")
  #   } else {
  #     out = if_else(length(.)>0,
  #                   sum(.),
  #                   0)
  #   }
  # )
  # 
  # # Doesn't work
  # function_to_apply <- quo(
  #   function(var_name){
  #     var_name <- enquo(var_name)
  #     if (!is.numeric(var_name) & !is.null(var_name)) {
  #         stop("transformer_sum is attempting to sum a non-numeric column")
  #       } else {
  #         out = if_else(length(var_name)>0,
  #                sum(var_name),
  #                0)
  #       }
  #     
  #   }
  # )
  
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




