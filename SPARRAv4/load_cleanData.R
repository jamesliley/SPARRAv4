# Update on Gergos original script

# Packages and scripts
library(tidyverse)
library(lubridate)
library(fst)

# source("James/Util/basic_operations.R")

#' Load the cleaned dataset
#' @param partition If "train", loads only the training data (a fixed 80%  of patients, with episodes in the last 12 months censored)
#'                  if "test", loads only the test patients (with no time censoring)
#'                  if "all" loads all data with neither patient nor time exclusions
#' @param subset_frac Between 0 and 1, randomly subsets the patients in the data to reduce in-memory data
#' @param subset_seed Random seed of data subsetting
#' @param load_to_global if TRUE, loads the outputs into the global environment, otherwise returns them as a list
#' @param load_sparrav3_scores if TRUE, loads all SPARRAv3 scores (huge data, we generally do not want it)
#' @param dir_cleanData Defines where to look for the clean data (in case it moves)
#' @param load_tables By default all data tables are loaded, but the user can decide to load a limited set only
#' 
#' Output
#' @return Loads 3 objects (patients, episodes, list_of_data_tables), either into the global environment, or outputs them as a named list
#' @export
load_cleanData <- function(
  partition="train", 
  subset_frac=1, subset_seed=1234,
  load_to_global = TRUE,
  load_sparrav3_scores = FALSE,
  dir_cleanData = "../Data/Data_clean", # Define the clean data folder
  load_tables = c("AE2", "deaths", "PIS", "SMR00", "SMR01", "SMR01E", "SMR04", 
                  "SPARRALTC", "SystemWatch")
){
  # Call the actual function
  out = load_cleanData_internal(partition=partition, 
                                subset_frac=subset_frac, subset_seed=subset_seed,
                                load_to_global = load_to_global,
                                load_sparrav3_scores =load_sparrav3_scores,
                                dir_cleanData = dir_cleanData,
                                load_tables = load_tables)
  
  # Run garbage collection to remove ~14 GB of unnecessarily allocated memory
  gc()
  
  # Return the output
  out
  
}


  
  
load_cleanData_internal <- function(
    partition="train", 
    subset_frac=1, subset_seed=1234,
    load_to_global = TRUE,
    load_sparrav3_scores = FALSE,
    dir_cleanData = "//FARR-FS1/Study Data/1718-0370/Research/Data/Data_clean", # Define the clean data folder
    load_tables = c("AE2", "deaths", "PIS", "SMR00", "SMR01", "SMR01E", "SMR04", 
                         "SPARRALTC", "SystemWatch")
  ){
 
  
  # If load_sparrav3_scores, add it to load_tables
  if (load_sparrav3_scores){
    load_tables <- c(load_tables, "SPARRA")
  }
  
  
  # --------------------------------------------------------------------------------------
  # Load patients and episodes, and filter appropriately based on train/test and subsetting (for now only on patient level)
  patients <- as_tibble(read.fst(path = file.path(dir_cleanData, "patients.fst")))
  
  episodes_metadata <- metadata_fst(
    path = file.path(
    dir_cleanData,
    "episodes.fst"))
  
  if (file.exists(file.path(
    dir_cleanData,
    "episodes_lookup_by_datatable.fst"))){
    episodes_lookup_by_datatable <- as_tibble(read_fst(file.path(
      dir_cleanData,
      "episodes_lookup_by_datatable.fst")))
    
    episodes_lookup_by_datatable <- episodes_lookup_by_datatable %>% 
      filter(source_table %in% load_tables)
    
    # Load the appropriate rows for each data table
    episodes <- bind_rows(
      episodes_lookup_by_datatable %>% 
        group_by(source_table) %>% 
        nest() %>% 
        mutate(cur_table = map(data, 
                               ~ as_tibble(
                                 read_fst(path = file.path(dir_cleanData, "episodes.fst"),
                                          from = .$row_min,
                                          to = .$row_max)
                               )
          )) %>%
        pull(cur_table)
                                 
    )
    
    
  } else {
  
    
    if (load_sparrav3_scores){
      episodes <- as_tibble(read.fst(path = file.path(dir_cleanData, "episodes.fst")))  
    } else {
      episodes <- as_tibble(read.fst(path = file.path(dir_cleanData, "episodes_no_sparrav3.fst")))
    }
  
  }
  
  # Create the train-test split for use within this loading function
  
  # Select 20% of patient ids (that will constitute the test set)
  set.seed(2713)
  patient_ids_unique <- (patients %>% select(id) %>% distinct())$id
  patient_test_sample <- sample(patient_ids_unique, ceiling(.20 * length(patient_ids_unique)), FALSE )
  
  # Use the sampled test patient set (alongside the type of partition) to choose a part of episodes we will use
  if (partition=="test"){
    # Take only patients within the test sample
    patients <- patients %>% 
      filter(id %in% patient_test_sample)
    
    # Take all times (no removal of last year of data)
    episodes <- episodes %>% 
      filter(id %in% patient_test_sample)
  } else if (partition == "train") {
    # Take only patients NOT within the test sample
    patients <- patients %>% 
      filter(!(id %in% patient_test_sample))
    
    #  Remove episodes within the last year of data
    episodes <- episodes %>% 
      filter(!(id %in% patient_test_sample)) %>%
      add_mutate_timeUNIX() %>%
      filter(timeUNIX < as.integer(as_datetime("2017-04-29"))) %>% # Last time in the dataset is "2018-04-30 23:59:00 UTC" (except in SPARRALTC!)
      remove_mutate_timeUNIX()
  } else if (partition=="all") {
    # placeholder, do nothing
  }
  
  rm(patient_ids_unique)
  rm(patient_test_sample)
  
  
  # If subset_frac is not 1, load only part of the data
  if (subset_frac < 1){
    set.seed(subset_seed)
    patient_ids_unique <- (patients %>% select(id) %>% distinct())$id
    patient_subset_sample <- sample(patient_ids_unique, ceiling(subset_frac * length(patient_ids_unique)), FALSE )
    # Take only patients within the subset sample
    patients <- patients %>% 
      filter(id %in% patient_subset_sample)
    
    # Take only patients within the subset sample
    episodes <- episodes %>% 
      filter(id %in% patient_subset_sample)
    
    rm(patient_ids_unique)
    rm(patient_subset_sample)
  }
  
  
  # Episodes and patients done

  
  # ---------------------------------------------------
  # Load the information from the detailed data tables, corresponding to the chosen patients / episodes
  
  # Define the table files to be used
  
  cleanData_table_filenames = c("Final_AE2_extract_incl_UniqueStudyID.fst", 
                          "Final_deaths_extract_incl_UniqueStudyID_v2.fst", 
                          "Final_PIS_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR00_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR01_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR01E_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR04_extract_incl_UniqueStudyID.fst", 
                          "Final_SPARRA_extract_incl_UniqueStudyID.fst", 
                          "final_LTCs_using_new_ICD10codes_incl_uniquestudyid.fst", 
                          "Final_SystemWatch_extract_incl_UniqueStudyID_v2.fst")
  
  table_names = c("AE2",
                  "deaths",
                  "PIS",
                  "SMR00",
                  "SMR01",
                  "SMR01E",
                  "SMR04",
                  "SPARRA",
                  "SPARRALTC",
                  "SystemWatch"
  )
  
  names(cleanData_table_filenames) <- table_names
  
  data_filenames <- cleanData_table_filenames[names(cleanData_table_filenames) %in% load_tables]
  
  
  
  # Load each file, filter the necessary rows and collect them into a list
  list_of_data_tables = list()
  
  for (name in names(data_filenames)){
    table_current <- as_tibble(read_fst(
      path = file.path(
        dir_cleanData,
        data_filenames[name]
      )
    ))
    
    # Keep only the rows that are part of the filtered "episodes"
    table_current <- table_current[episodes %>% filter(source_table==name) %>% pull(source_row), ]
    
    # Add the filtered table to our list of tables
    list_of_data_tables[[name]] = table_current
    
    rm(table_current)
  }
  
  
  # # TEsting random access speed (not worth row-by-row, may be worth chunk-by-chunk but probably not, only if subset_frac << 1)
  # tic(); tmp <- read.fst(path = file.path(dir_cleanData, "episodes.fst"), from = 1000000, to = 2000000); toc();
  # 
  # tmp <- tmp[1,]
  # 
  # tic()
  # for (i in 1000000:1002000){
  #   tmp <- tmp %>% rbind(read.fst(path = file.path(dir_cleanData, "episodes.fst"), from = i, to = i))
  # }
  # toc()
  
  
  if (load_to_global){
    assign("patients", patients, envir=.GlobalEnv) 
    assign("episodes", episodes, envir=.GlobalEnv) 
    assign("list_of_data_tables", list_of_data_tables, envir=.GlobalEnv) 
    out = NULL
  } else {
    out = list(
      patients = patients,
      episodes = episodes,
      list_of_data_tables = list_of_data_tables
    )
  }
  
  
  
  # Return value (either NULL or list of objects)
  out
}





# Get a minibatch from loaded data ( this subsets patients, episodes and list_of_data_tables in .GlobalEnv based on patient id !!!)
extract_patient_minibatch_from_data <- function(
  cur_batchnum = 1,
  num_batches = 10
){
  
  patient_ids_unique <- (patients %>% select(id) %>% distinct() )$id
  patient_ids_unique <- sort(patient_ids_unique)
  
  # Cut the patient ids into num_batches chunks and select the cur_batchnum chunk
  cur_minibatch_patient_ids <- split(patient_ids_unique, cut(seq_along(patient_ids_unique), num_batches, labels=FALSE))[[cur_batchnum]]
  
  # Replace the objects within the global environment
  assign("patients", patients %>% filter(id %in% cur_minibatch_patient_ids), envir=.GlobalEnv) 
  assign("episodes", episodes %>% filter(id %in% cur_minibatch_patient_ids), envir=.GlobalEnv) 
  assign("list_of_data_tables", lapply(list_of_data_tables, 
                                       FUN = function(cur_data_table){
                                         cur_data_table %>% filter(id %in% cur_minibatch_patient_ids)
                                       }), envir=.GlobalEnv)
  
  # Run garbage collection to remove unnecessarily allocated memory
  gc()
}



load_cleanData_minibatch <- function(
  partition="train", 
  subset_frac=1., subset_seed=1234,
  cur_batchnum = 1,
  num_batches = 10,
  load_to_global = TRUE,
  load_sparrav3_scores = FALSE
){
  # Call the actual function
  out = load_cleanData_internal(partition=partition, 
                                subset_frac=subset_frac, subset_seed=subset_seed,
                                load_to_global = load_to_global,
                                load_sparrav3_scores =load_sparrav3_scores)
  
  # Run garbage collection to remove ~14 GB of unnecessarily allocated memory
  gc()
  
  # Filter the output data objects based on people's ids to create minibatches
  if (load_to_global == FALSE){
    stop("Minibatch data loading not implemented for \"load_to_global==FALSE\" in load_cleanData_minibatch()")
  }
  
  # Subset the loaded data in place
  extract_patient_minibatch_from_data(cur_batchnum = cur_batchnum, num_batches = num_batches)
  
  # Return NULL
  NULL
  
}