# Reads raw data from database files, and generates usable objects
# Assigns cross-validation folds at the patient level.
# Merges SMR01, SystemWatch and SMR01E into a table SMR01M
# Gergo Bohner and James Liley
# 11/19

library(Matrix)
library(tidyverse)
library(lubridate)
library(fst)
library(haven)



##' Clean raw data and generate tidy versions of data tables. Will only update episodes table if force_redo = TRUE.
##' 
##' @param force_redo Redo if clean data are not already there
##' @param dir_rawData directory containing raw data tables (from EHRs)
##' @param dir_cleanData directory to which clean data tables will be written
##' @returns NULL
clean_rawData <- function(force_redo = TRUE, 
                          dir_rawData = "//Farr-FS1/Study Data/1718-0370/Linked Data",
                          dir_cleanData = "//FARR-FS1/Study Data/1718-0370/Research/Data/Data_clean"){
  
  # Raw data files
  data_filenames = c("Final_AE2_extract_incl_UniqueStudyID.zsav", 
                     "Final_deaths_extract_incl_UniqueStudyID_v2.zsav", 
                     "Prescribing_pivoted_201305_to_201404.csv.gz",
                     "Prescribing_pivoted_201405_to_201504.csv.gz",
                     "Prescribing_pivoted_201505_to_201604.csv.gz",
                     "Prescribing_pivoted_201605_to_201704.csv.gz",
                     "Prescribing_pivoted_201705_to_201804.csv.gz",
                     "Final_SMR00_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR01_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR01E_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR04_extract_v2_for_ATI_incl_prev_psych_care.sav",
                     "FINAL_v4_LTC.zsav",
                     "Final_SystemWatch_extract_incl_UniqueStudyID_v2.zsav",
                     "LOOKUP_SIMD_sex_dob_for_ATI.zsav")
  
  table_names = c("AE2",
                  "deaths",
                  "PIS_1314",
                  "PIS_1415",
                  "PIS_1516",
                  "PIS_1617",
                  "PIS_1718",
                  "SMR00",
                  "SMR01",
                  "SMR01E",
                  "SMR04",
                  "SPARRALTC",  
                  "SystemWatch",
                  "Patient_lookup"
                  )
  
  names(data_filenames) <- table_names
  
  source_names = c("PIS","SMR01M",grep("PIS",table_names,val=T,invert=T))
  
  
  # Helper function to turn "numerise-able" columns into real numeric columns
  numerise_columns <- function(.data) {
    .data %>% 
      mutate_if(
        ~(is.character(.) | is.factor(.)),
        function(x){
          if (all(check.numeric(x))){ 
            as.numeric(x)
          } else
            x
        } 
      )
  }
  
  cleanData_filenames = c("episodes.fst",
                          "episodes_lookup.fst",
                          "AE2.fst",
                          "deaths.fst",
                          "PIS.fst",
                          "SMR00.fst",
                          "SMR01.fst",
                          "SMR01E.fst",
                          "SMR04.fst",
                          "SPARRALTC.fst",
                          "SystemWatch.fst",
                          "SMR01M.RDS",
                          "patients.fst")
  

  if (force_redo==FALSE && all(cleanData_filenames %in% list.files(dir_cleanData))){
    cat("Data is already cleaned and all clean data files exist. Set force_redo==TRUE to redo data cleaning")
    
  } else {

    episodes = tibble(id = as.integer(1), 
                      time = lubridate::as_datetime("1990/10/21"),
                      time_discharge = lubridate::as_datetime("1994/01/13"),
                      source_table = factor("SMR01", levels = source_names),
                      source_row = as.integer(c(1)),
                      admission_type = as.integer(c(1)),
                      main_condition = as.character(c("G405H")), # ICD10 codes
                      CIS_MARKER = as.integer(c(1)) # Continuous inpatient stay number (from "SMR01" and "SMR04")
    )[-1,]
  
    update_episodes <- function(episodes, new_table){
      episodes %>%
        bind_rows(
          new_table %>%
            select(colnames(episodes)[colnames(episodes) %in% colnames(.)])
        )
    }
  
  
    #**********************************************************#  
    # TIDYING INDIVIDUAL DATA TABLES -------------------------
    #**********************************************************#  
    
    cat("Tidying individual tables.\n \n")
      
    
    #**********************************************************#  
    # AE2 --------------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"AE2.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      
      table_AE2 = haven::read_sav(file.path(dir_rawData, data_filenames[["AE2"]]))
      
      glimpse(table_AE2)
      
      # Parse the times
      table_AE2 <- table_AE2 %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("AE2", levels = source_names),
          source_row = as.integer(row_number()),
          time = lubridate::parse_date_time(
            str_c(ADMISSION_DATE, ADMISSION_TIME, sep=" "),
            orders = c("d m Y H M")
          ),
          time_discharge = lubridate::parse_date_time(
            str_c(TRANSFER_DISCHARGE_DATE, TRANSFER_DISCHARGE_TIME, sep=" "),
            orders = c("d m Y H M")
          )
        ) %>%
        select(-c(UNIQUE_STUDY_ID, ADMISSION_DATE, ADMISSION_TIME, TRANSFER_DISCHARGE_DATE, TRANSFER_DISCHARGE_TIME)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_AE2)
      
      # Save the parsed table to compressed RData
      write.fst(table_AE2,
           path = filepath_clean, compress = 100
           )
      
      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_AE2)
      
      # Remove from memory
      rm(table_AE2)
      
    }
    
    
    
    
    #**********************************************************#  
    # SMR00 ------------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"SMR00.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SMR00 = haven::read_sav(file.path(dir_rawData, data_filenames[["SMR00"]])) 
      
      glimpse(table_SMR00)
      
      # Parse the useful information
      table_SMR00 <- table_SMR00 %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("SMR00", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_birth = as_date(DOB), 
          gender = labelled(as.integer(SEX), c(Male =1, Female=2)),
          time = as_datetime(CLINIC_DATE),
          main_condition = as.character(MAIN_CONDITION)
        ) %>%
        select(-c(UNIQUE_STUDY_ID, DOB, CLINIC_DATE, SEX, MAIN_CONDITION)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_SMR00)
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR00,
           path = filepath_clean, compress = 100
      )
      
      
      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SMR00)
      
      # Remove from memory
      rm(table_SMR00)
      gc()
      
    }
    
    
    
    
    #**********************************************************#  
    # SMR01 ------------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"SMR01.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SMR01 = haven::read_sav(file.path(dir_rawData, data_filenames[["SMR01"]]))
      
      glimpse(table_SMR01)
      
      # Parse the useful information
      table_SMR01 <- table_SMR01 %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("SMR01", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_birth = as_date(DOB),
          gender = labelled(as.integer(SEX), c(Male =1, Female=2)),
          time = as_datetime(ADMISSION_DATE),
          time_discharge = as_datetime(DISCHARGE_DATE),
          admission_type = as.integer(ADMISSION_TYPE),
          main_condition = as.character(MAIN_CONDITION)
        
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID, DOB, ADMISSION_DATE, DISCHARGE_DATE, SEX, ADMISSION_TYPE, MAIN_CONDITION)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_SMR01)
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR01,
           path = filepath_clean, compress = 100
      )

      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SMR01)
      
      # Remove from memory
      rm(table_SMR01)
      gc()
      
    }
    
    
    
    
    
    
    #**********************************************************#  
    # SMR01E -----------------------------
    #**********************************************************#  
    
    
    
    filepath_clean = file.path(dir_cleanData,"SMR01E.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SMR01E = haven::read_sav(file.path(dir_rawData, data_filenames[["SMR01E"]]))
      
      glimpse(table_SMR01E)
      
      # Parse the useful information
      table_SMR01E <- table_SMR01E %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("SMR01E", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_birth = as_date(DOB),
          gender = labelled(as.integer(SEX), c(Male =1, Female=2)),
          time = as_datetime(ADMISSION_DATE),
          time_discharge = as_datetime(DISCHARGE_DATE),
          admission_type = as.integer(ADMISSION_TYPE),
          main_condition = as.character(MAIN_CONDITION)
          
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID, DOB, ADMISSION_DATE, DISCHARGE_DATE, SEX, ADMISSION_TYPE, MAIN_CONDITION)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR01E,
           path = filepath_clean, compress = 100
      )

      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SMR01E)
      
      # Remove from memory
      rm(table_SMR01E)
      
    }
    
    

    
    
    
    #**********************************************************#  
    # SystemWatch ------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"SystemWatch.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SystemWatch = haven::read_sav(file.path(dir_rawData, data_filenames[["SystemWatch"]]))
      
      glimpse(table_SystemWatch)
      
      # Parse the useful information
      table_SystemWatch <- table_SystemWatch %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("SystemWatch", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_birth = as_date(parse_date_time(DATE_OF_BIRTH, orders = "d m y")),
          gender = labelled(as.integer(SEX), c(Male =1, Female=2)),
          time = as_datetime(parse_date_time(DATE_OF_ADMISSION, orders = "d m y")),
          time_discharge = as_datetime(parse_date_time(DATE_OF_DISCHARGE, orders = "d m y")),
          admission_type = as.integer(ADMISSION_TYPE),
          main_condition = ifelse(as.character(DIAGNOSIS_1) %in% c(""), NA, as.character(DIAGNOSIS_1))
          
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID, DATE_OF_BIRTH, DATE_OF_ADMISSION, DATE_OF_DISCHARGE, SEX, ADMISSION_TYPE, DIAGNOSIS_1)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SystemWatch,
                path = filepath_clean, compress = 100
      )

      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SystemWatch)
      
      
      # Remove from memory
      rm(table_SystemWatch)
      
    }
    
    
    
    #**********************************************************#  
    # Merge SMR01,SMR01E, SystemWatch  ---------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"SMR01M.RDS")
    cat(paste0("... merging: ", filepath_clean, "\n"))

    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      
    # List of data tables and convert names to lowercase
    l1=list(
      SMR01=as_tibble(read_fst(path = file.path(dir_cleanData,"SMR01.fst"))),
      SMR01E=as_tibble(read_fst(path = file.path(dir_cleanData,"SMR01E.fst"))),
      SystemWatch=as_tibble(read_fst(path = file.path(dir_cleanData,"SystemWatch.fst")))
    )
    for (i in 1:length(l1)) colnames(l1[[i]])=tolower(colnames(l1[[i]]))
    
    SMR01M=smr01_systemwatch_merge(
      l1,
      lubridate::ymd("2018-05-01"),
      use_sw_diagnosis_codes = TRUE,
      include_smr01e= TRUE)
    
    # Get rid of length_of_stay, so we don't accidentally use it
    SMR01M = SMR01M %>% select(-length_of_stay)
    
    # Convert date_of_birth to POSIXct
    SMR01M$date_of_birth= as.POSIXct(SMR01M$date_of_birth)
    
    # Glimpse
    glimpse(SMR01M)
    
    # Save the parsed table to RData. Can't compress lists so can't use FST here.
    saveRDS(SMR01M,file=filepath_clean)
    
    # Add information to the episodes table for tables SMR01, SMR01E, SystemWatch
    if  (force_redo)  {
      episodes = episodes %>% update_episodes(SMR01M)
    
      # Now add information to the episodes table for SMR01M
      SMR01M$source_table=factor("SMR01M",levels=source_names)
      SMR01M$source_row=1:dim(SMR01M)[1]
      episodes = episodes %>% update_episodes(SMR01M)
      
    }        
    # Remove from memory
    rm(SMR01M)
    rm(l1)

  }
  
    
        
    
    #**********************************************************#  
    # SMR04 ------------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"SMR04.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SMR04 = haven::read_sav(file.path(dir_rawData, data_filenames[["SMR04"]]))
      
      glimpse(table_SMR04)
      
      # Parse the useful information
      table_SMR04 <- table_SMR04 %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("SMR04", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_birth = as_date(DOB),
          gender = labelled(as.integer(SEX), c(Male =1, Female=2)),
          time = as_datetime(ADMISSION_DATE),
          time_discharge = as_datetime(DISCHARGE_DATE),
          admission_type = as.integer(ADMISSION_TYPE),
          main_condition = as.character(MAIN_CONDITION)
          
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID, DOB, ADMISSION_DATE, DISCHARGE_DATE, SEX, ADMISSION_TYPE, MAIN_CONDITION)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR04,
           path = filepath_clean, compress = 100
      )
      
      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SMR04)
      
      # Remove from memory
      rm(table_SMR04)
      
    }
    
    
    
    
    
    
    #**********************************************************#  
    # PIS --------------------------------
    #**********************************************************#  
    
    # Combine all 5 PIS files
    # data_filenames[["PIS"]] = "Final_PIS_extract_incl_UniqueStudyID.zsav"
    
    filepath_clean = file.path(dir_cleanData,"PIS.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      tmp_first_load = TRUE
      for (fname in data_filenames[grep("PIS_", names(data_filenames))]){
        
        table_PIS_tmp = read.csv(file.path(dir_rawData, fname))
        cat(paste(fname, " loaded, processing..."))
        glimpse(table_PIS_tmp)
        
        # Set "id" column
        table_PIS_tmp <- table_PIS_tmp %>%
          mutate(id = as.integer(unique_study_id)) %>%
          mutate(time=as_datetime(lubridate::ym(month_time_period)))

        
        if (tmp_first_load){
          table_PIS = table_PIS_tmp
          tmp_first_load = FALSE
        } else {
          table_PIS <- table_PIS %>%
            bind_rows(
              table_PIS_tmp
            )
        }
        
        rm(table_PIS_tmp)
        gc()
      }
  
      # Add table key for the collected data
      table_PIS <- table_PIS %>%
        mutate(
          source_table = factor("PIS", levels = source_names),
          source_row = as.integer(row_number()),
          bnf_section = as.factor(toupper(bnf_section))
        ) 
        
      
      glimpse(table_PIS)
      
      
      write.fst(table_PIS,
           path = filepath_clean, compress = 100
      )
      
      
      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_PIS)
      
      # Remove from memory
      rm(table_PIS)
      gc()
      
    }
    
    
   

    
    #**********************************************************#  
    # SPARRALTCs -------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData, "SPARRALTC.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SPARRALTC = haven::read_sav(file.path(dir_rawData, data_filenames[["SPARRALTC"]]))
      
      glimpse(table_SPARRALTC)
      
      # Parse the useful information
      table_SPARRALTC <- table_SPARRALTC %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          NUMBEROFLTCs = as.integer(NUMBER_OF_LTCS)
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID)) %>%
        select(-c(NUMBER_OF_LTCS)) %>%
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        mutate_if(is.Date, ~as_date(format(.,"%d-%m-%Y"), format="%d-%m-%Y")) %>% 
        tidyr::gather(LTC_TYPE, time, 
                      FIRST_ARTHRITIS_EPISODE:FIRST_OTHER_ENDOCRINE_METABOLIC_DISEASES_EPISODE,
                      factor_key=TRUE) %>% 
        filter(!is.na(time)) %>% 
        mutate(
          source_table = factor("SPARRALTC", levels = source_names),
          source_row = as.integer(row_number())
        ) %>%
        mutate(
          time = as_datetime(time)
        ) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SPARRALTC,
           path = filepath_clean, compress = 100
      )
      
      
      # Add the information to the episodes table
      if  (force_redo)  episodes = episodes %>%
        update_episodes(table_SPARRALTC)
      
      # Remove from memory
      rm(table_SPARRALTC)
      
    }
    
    

        
    #**********************************************************#  
    # Deaths --------------------------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"deaths.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_deaths = haven::read_sav(file.path(dir_rawData, data_filenames[["deaths"]]))
      
      # Fix incorrectly formatted times
      d1=lubridate::as_date(table_deaths$DATE_OF_DEATH, format="%d-%b-%y")
      w=which(is.na(d1))
      d2=lubridate::as_date(table_deaths$DATE_OF_DEATH[w])
      d1[w]=d2
      table_deaths$DATE_OF_DEATH=d1
      
      
      
      
      # Parse the times
      table_deaths <- table_deaths %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("deaths", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_death = lubridate::as_date(DATE_OF_DEATH, format="%d-%b-%y"),
          time = as_datetime(date_of_death)
        ) %>%
        
        select(-c(UNIQUE_STUDY_ID, DATE_OF_DEATH)) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_deaths)
      
      # Save the parsed table to compressed RData
      write.fst(table_deaths,
                path = filepath_clean, compress = 100
      )
      
      # Add the information to the episodes table
      if  (force_redo==TRUE)  episodes = episodes %>% update_episodes(table_deaths)
      
      # Remove from memory
      #rm(table_deaths)
      
    } else {
      table_deaths <- as_tibble(read_fst(path = filepath_clean))
    }
    
    
    
    
    #**********************************************************#  
    # Patient_lookup table ---------------
    #**********************************************************#  
    
    filepath_clean = file.path(dir_cleanData,"patients.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      
      # Retain deaths (that is not in the lookup_table)
      table_deaths <- as_tibble(read_fst(path = file.path(dir_cleanData,"deaths.fst")))
      patient_deaths <- table_deaths %>% select(id, date_of_death)
      rm(table_deaths)
      
      
      table_patients = haven::read_sav(file.path(dir_rawData, data_filenames[["Patient_lookup"]]))
      
      glimpse(table_patients)
      
      table_patients <- table_patients %>% 
        rename(id:=UNIQUE_STUDY_ID ) %>% 
        mutate_if(is.character, ~ifelse(nchar(.)==0, NA, .)) %>%
        mutate(
          id=as.integer(id),
          date_of_birth = lubridate::as_date(lubridate::fast_strptime(DOB, format="%d-%m-%Y")),
          gender=as.factor(GENDER),
          SIMD_DECILE_2016_SCT = as.integer(SIMD_DECILE_2016_SCT)
        ) %>% 
        select(-c(SIMD_QUINTILE_2016_SCT))
      
      patients <- table_patients %>%
        select(id, gender, date_of_birth, SIMD_DECILE_2016_SCT) %>%
        left_join(
          patient_deaths, by="id"
        )
    }
    
    
    
    #**********************************************************#  
    # Save the resulting patients and episodes files
    #**********************************************************#  
    
    if (force_redo) {
    
    # Ensure correct data types
    patients <- patients %>% 
      transmute(
        id = as.integer(id),
        gender = as.factor(gender),
        date_of_birth = as_date(date_of_birth),
        date_of_death = as_date(date_of_death),
        SIMD_DECILE_2016_SCT = as.integer(SIMD_DECILE_2016_SCT)
      )
      
    
    episodes <- episodes %>%
      transmute(
        id = as.integer(id),
        time = as_datetime(time),
        time_discharge = as_datetime(time_discharge),
        source_table = factor(source_table, levels = source_names),
        source_row = as.integer(source_row),
        admission_type = as.integer(admission_type),
        main_condition = as.character(main_condition), # ICD10 codes
        CIS_MARKER = as.integer(CIS_MARKER)
      )

    
    # Generate cross-validation fold assignments by patient
    k=3
    set.seed(212)
    nt=nrow(patients)
    cvs=rep(1:k,1+floor(nt/k))[1:nt][order(runif(nt))] # fold assignments; equal size folds, randomised
    patients$cv=cvs
        
    write.fst(patients,
         path = file.path(
           dir_cleanData,
           "patients.fst"))
    
    write.fst(episodes,
         path = file.path(
           dir_cleanData,
           "episodes.fst"), 
         compress = 100)
    
    episodes_lookup_by_datatable <- episodes %>% 
      mutate(rownum = as.integer(row_number())) %>% 
      group_by(source_table) %>% 
      summarise(row_min = min(rownum), row_max = max(rownum)) %>%
      arrange(row_max) %>%
      mutate_at(vars(row_min, row_max), as.integer)
    
    write.fst(episodes_lookup_by_datatable ,
              path = file.path(dir_cleanData,"episodes_lookup.fst"))
    
    cat(paste0("\n\nDATA CLEANING FINISHED\n----------------------\n"))
    
    }
  }
  
}

