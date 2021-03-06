# Update to Gergo's original script.
# Reads raw data from database files, and generates usable objects
# Update also assigns cross-validation folds at the patient level.
# James Liley
# 11/19

library(Matrix)
library(tidyverse)
library(lubridate)
#library(varhandle)
library(fst)
library(haven)



##' Clean raw data and generate tidy versions of data tables
##' 
##' @param force_redo Redo if clean data are not already there
##' @param dir_rawData directory containing raw data tables (from EHRs)
##' @param dir_cleanData directory to which clean data tables will be written
##' @returns NULL
clean_rawData <- function(force_redo = FALSE, dir_rawData = "//Farr-FS1/Study Data/1718-0370/Linked Data",
  dir_cleanData = "//FARR-FS1/Study Data/1718-0370/Research/Data/Data_clean"){
  
  # Vector generated by dput(list.files(dir_rawData))
  data_filenames = c("Final_AE2_extract_incl_UniqueStudyID.zsav", 
                     "Final_deaths_extract_incl_UniqueStudyID_v2.zsav", 
                     "PIS_201305_to_201404_for_ATI_monthly_ID.zsav", 
                     "PIS_201405_to_201504_for_ATI_monthly_ID.zsav", 
                     "PIS_201505_to_201604_for_ATI_monthly_ID.zsav", 
                     "PIS_201605_to_201704_for_ATI_monthly_ID.zsav", 
                     "PIS_201705_to_201804_for_ATI_monthly_ID.zsav", 
                     "Final_SMR00_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR01_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR01E_extract_incl_UniqueStudyID.zsav", 
                     "Final_SMR04_extract_v2_for_ATI_incl_prev_psych_care.sav",
                     "SPARRAscores_1July13_1Apr16_ATI.sav",  # Old SPARRA v3 scores
                     "Final_SPARRA_extract_incl_UniqueStudyID.zsav",  # Recent SPARRA v3 scores
 #                    "Final_SPARRALTC_extract_incl_UniqueStudyID.zsav",  # Old SPARRA long-term conditions
                     "final_LTCs_using_new_ICD10codes_incl_uniquestudyid.zsav", # Recent SPARRA long-term conditions as of 11/19
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
                  "SPARRA_old",
                  "SPARRA_recent",
 #                 "SPARRALTC_old",
                  "SPARRALTC",  # New as of 11/19
                  "SystemWatch",
                  "Patient_lookup"
                  )
  
  
  table_names
  
  names(data_filenames) <- table_names
  
  source_names = c(table_names[1:2], list("PIS"), table_names[8:11], list("SPARRA"), table_names[14:17])
  
  
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
  
  
  cleanData_filenames = c("episodes.fst", "episodes_lookup_by_datatable.fst", "Final_AE2_extract_incl_UniqueStudyID.fst", 
                          "Final_deaths_extract_incl_UniqueStudyID_v2.fst", "Final_PIS_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR00_extract_incl_UniqueStudyID.fst", "Final_SMR01_extract_incl_UniqueStudyID.fst", 
                          "Final_SMR01E_extract_incl_UniqueStudyID.fst", "Final_SMR04_extract_incl_UniqueStudyID.fst", 
                          "Final_SPARRA_extract_incl_UniqueStudyID.fst", "final_LTCs_using_new_ICD10codes_incl_uniquestudyid.fst", 
                          "Final_SystemWatch_extract_incl_UniqueStudyID_v2.fst", "patients.fst")
  
  if (force_redo==FALSE && all(cleanData_filenames %in% list.files(dir_cleanData))){
    cat("Data is already cleaned, all clean data files exist. Set force_redo==TRUE to redo data cleaning")
    
  } else {
    # Create organised and small tables that can refer to the raw data
    patients = tibble(id = as.integer(1),
                      #sex = factor("male", levels = c("male", "female")),
                      gender = labelled(as.integer(c(1)), c(Female=1, Male=2)),
                      date_of_birth = lubridate::as_date("1990/10/21"),
                      date_of_death = lubridate::as_date("2190/10/21"),
                      source_table = factor("SMR01", levels = source_names),
                      source_row = as.integer(c(1))
    )[-1,]
    
    episodes = tibble(id = as.integer(1), 
                      time = lubridate::as_datetime("1990/10/21"),
                      time_discharge = lubridate::as_datetime("1994/01/13"),
                      source_table = factor("SMR01", levels = source_names),
                      source_row = as.integer(c(1)),
                      admission_type = as.integer(c(1)),
                      main_condition = as.character(c("G405H")), # ICD10 codes
                      CIS_MARKER = as.integer(c(1)) # Continuous inpatient stay number (from "SMR01" and "SMR04")
    )[-1,]
  
  
    # Create update functions for patients and episodes
    update_patients <- function(patients, new_table){
      patients %>%
        bind_rows(
          new_table %>%
            select(colnames(patients)[colnames(patients) %in% colnames(.)]) %>%
            distinct()
        )
    }
    
    update_episodes <- function(episodes, new_table){
      episodes %>%
        bind_rows(
          new_table %>%
            select(colnames(episodes)[colnames(episodes) %in% colnames(.)])
        )
    }
  
  
  
  
  
  
  
  
    # ---------------------------------------------------------------------------
    # TIDYING INDIVIDUAL DATA TABLES --------------------------------------------
    # ---------------------------------------------------------------------------
  
    cat("Tidying individual tables...\n \n")
      
    # ------------------------------------
    # AE2 --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["AE2"]], 1, -5), "fst"))
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_AE2)
      
      # Save the parsed table to compressed RData
      write.fst(table_AE2,
           path = filepath_clean, compress = 100
           )
    } else {
      table_AE2 <- as_tibble(read_fst(path = filepath_clean))
    }
    
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_AE2)
    
    episodes <- episodes %>%
      update_episodes(table_AE2)
    
    # Remove from memory
    rm(table_AE2)
    
    
    
    # ------------------------------------
    # Deaths --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["deaths"]], 1, -5), "fst"))
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_deaths = haven::read_sav(file.path(dir_rawData, data_filenames[["deaths"]]))
      
      # Parse the times
      table_deaths <- table_deaths %>% 
        mutate(
          id = as.integer(UNIQUE_STUDY_ID),
          source_table = factor("deaths", levels = source_names),
          source_row = as.integer(row_number()),
          date_of_death = lubridate::as_date(DATE_OF_DEATH, tz="", format="%d-%b-%y"),
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
    } else {
      table_deaths <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables 
    patients <- patients %>%
      update_patients(table_deaths)
    
    episodes <- episodes %>%
      update_episodes(table_deaths)
    
    # Remove from memory
    rm(table_deaths)
    
    
    
    # ------------------------------------
    # SMR00 --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["SMR00"]], 1, -5), "fst"))
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_SMR00 = haven::read_sav(file.path(dir_rawData, data_filenames[["SMR00"]])) # This table already has times as "DATE" type, so parsed correctly
      
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_SMR00)
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR00,
           path = filepath_clean, compress = 100
      )
    } else {
      table_SMR00 <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_SMR00)
    
    episodes <- episodes %>%
      update_episodes(table_SMR00)
    
    # Remove from memory
    rm(table_SMR00)
    gc()
    
    
    # ------------------------------------
    # SMR01 --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["SMR01"]], 1, -5), "fst"))
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      glimpse(table_SMR01)
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR01,
           path = filepath_clean, compress = 100
      )
    } else {
      table_SMR01 <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables
    patients <- patients %>%
      #filter(!(source_table %in% c("SMR01"))) %>%
      update_patients(table_SMR01)
    
    episodes <- episodes %>%
      #filter(!(source_table %in% c("SMR01"))) %>%
      update_episodes(table_SMR01)
    
    # Remove from memory
    rm(table_SMR01)
    gc()
    
    # ------------------------------------
    # SMR01E --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["SMR01E"]], 1, -5), "fst"))
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR01E,
           path = filepath_clean, compress = 100
      )
    } else {
      table_SMR01E <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables
    patients <- patients %>%
      #filter(!(source_table %in% c("SMR01E"))) %>%
      update_patients(table_SMR01E)
    
    episodes <- episodes %>%
      #filter(!(source_table %in% c("SMR01E"))) %>%
      update_episodes(table_SMR01E)
    
    # Remove from memory
    rm(table_SMR01E)
    
    
    
    # ------------------------------------
    # SMR04 --------------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      "Final_SMR04_extract_incl_UniqueStudyID.fst")
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        numerise_columns() %>% 
        filter(!is.na(id)) # No point keeping records with no ID. Typically very few.
      
      # Save the parsed table to compressed RData
      write.fst(table_SMR04,
           path = filepath_clean, compress = 100
      )
    } else {
      table_SMR04 <- as_tibble(read_fst(path = filepath_clean))
    }
    
    
    # Add the information to the small tables
    patients <- patients %>%
      #filter(!(source_table %in% c("SMR04"))) %>%
      update_patients(table_SMR04)
    
    episodes <- episodes %>%
      #filter(!(source_table %in% c("SMR04"))) %>%
      update_episodes(table_SMR04)
    
    # Remove from memory
    rm(table_SMR04)
    
    
    
    
    # ------------------------------------
    # PIS --------------------------------
    # ------------------------------------
    
    # Combine all 5 PIS files
    data_filenames[["PIS"]] = "Final_PIS_extract_incl_UniqueStudyID.zsav"
    
    
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["PIS"]], 1, -5), "fst"))
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      tmp_first_load = TRUE
      for (fname in data_filenames[grep("PIS_", names(data_filenames))]){
        
        table_PIS_tmp = haven::read_sav(file.path(dir_rawData, fname))
        
        cat(paste(fname, " loaded, processing..."))
        
        glimpse(table_PIS_tmp)
        
        # Set "id" column
        table_PIS_tmp <- table_PIS_tmp %>%
          mutate(id = as.integer(UNIQUE_STUDY_ID)) 
        
        # Retain the pieces of data that we need to NOT sparsify
        table_PIS_tmp_NUM_BNF_colnames <- colnames(table_PIS_tmp %>% select(starts_with("NUM_BNF")))
        table_PIS_tmp_nonsparse_data <- table_PIS_tmp %>% 
          select(PAID_GIC_INCL_BB, NUMBER_OF_PAID_ITEMS, id, YEARMONTH ) %>%
          mutate(NUMBER_OF_PAID_ITEMS = as.integer(NUMBER_OF_PAID_ITEMS))
        
        # This has ~140 GB of overhead..... (WHY?!)
        table_PIS_tmp <- table_PIS_tmp %>%
          select(starts_with("NUM_BNF"), id, YEARMONTH) %>% 
          mutate_at(vars(starts_with("NUM_BNF")), as.integer) %>%
          mutate(YEARMONTH_id = str_c(YEARMONTH, as.character(id))) %>%
          select(-id, -YEARMONTH)
        
        gc()
        
        table_PIS_tmp <- table_PIS_tmp %>%
          data.table::as.data.table() 
        
        gc()
        
        table_PIS_tmp <- table_PIS_tmp %>%
          as.matrix(rownames="YEARMONTH_id") %>% 
          as("sparseMatrix") 
        
        gc()
        
        table_PIS_tmp <- table_PIS_tmp %>%
          broom::tidy() %>%
          as_tibble %>%
          rename(
            YEARMONTH_id:=row, 
            BNF_section:=column,
            NUM_sold:=value
            ) %>%
          transmute(
            id = as.integer(str_sub(YEARMONTH_id, start = 6)),
            YEARMONTH = str_sub(YEARMONTH_id, end=5),
            BNF_section = factor(BNF_section, 
                                 levels = table_PIS_tmp_NUM_BNF_colnames),
            NUM_sold = as.integer(NUM_sold)
            ) %>%
          left_join(
            table_PIS_tmp_nonsparse_data
          )
        
        gc()
        
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
      
      # Create recoding table to turn "YEARMONTH" strings into actual times
      yearmonth_codes <- c()
      for (i11 in 1:5){
        for (j11 in 1:12){
          yearmonth_codes <- c(yearmonth_codes, stringr::str_interp("Y$[01i]{i11}M$[02i]{j11}"))
        } 
      }
      
      yearmonth_times = seq(as_datetime(ymd_hm("2013-05-01 00:00")), as_datetime(ymd_hm("2018-04-01 00:00")), by = "months")
      
      # Add table key for the collected data
      table_PIS <- table_PIS %>%
        left_join( # Recode the YEARMONTH column into actual times via this left_join
          tibble(
            YEARMONTH=yearmonth_codes,
            time = yearmonth_times)
        ) %>%
        arrange(id, time, BNF_section) %>% 
        mutate(
          source_table = factor("PIS", levels = source_names),
          source_row = as.integer(row_number())
        ) 
        
      
      glimpse(table_PIS)
      
      
      write.fst(table_PIS,
           path = filepath_clean, compress = 100
      )
      
    
    }  else {
      table_PIS <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # this leads to other issues down the line, AVOID! (for now)
    # table_PIS_for_update <- table_PIS %>%
    #   group_by(id, time) %>%
    #   summarise(source_row = min(source_row)) %>%
    #   ungroup() %>%
    #   mutate(source_table = factor("PIS", levels = source_names))
    
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_PIS)
    
    episodes <- episodes %>%
      update_episodes(table_PIS)
    
    # Remove from memory
    rm(table_PIS)
    #rm(table_PIS_for_update)
    gc()
    
    
    
    
    # ------------------------------------
    # SystemWatch ------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["SystemWatch"]], 1, -5), "fst"))
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
    } else {
      table_SystemWatch <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_SystemWatch)
    
    episodes <- episodes %>%
      update_episodes(table_SystemWatch)
    
    # Remove from memory
    rm(table_SystemWatch)
    
   
    
    
    # ------------------------------------
    # # SPARRALTC_old -----------------------------
    # # ------------------------------------
    # filepath_clean = file.path(
    #   dir_cleanData,
    #   str_c(str_sub(data_filenames[["SPARRALTC_old"]], 1, -5), "fst"))
    # cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    # 
    # if ((!file.exists(filepath_clean)) | force_redo==TRUE){
    #   table_SPARRALTC_old = haven::read_sav(file.path(dir_rawData, data_filenames[["SPARRALTC_old"]]))
    #   
    #   glimpse(table_SPARRALTC_old)
    #   
    #   # Parse the useful information
    #   table_SPARRALTC_old <- table_SPARRALTC_old %>% 
    #     mutate(
    #       id = as.integer(UNIQUE_STUDY_ID),
    #       NUMBEROFLTCs = as.integer(NUMBEROFLTCs)
    #     ) %>%
    #     
    #     select(-c(UNIQUE_STUDY_ID)) %>%
    #     mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
    #     mutate_if(is.character, funs(as_date(., tz="", format="%d-%m-%Y"))) %>% 
    #     tidyr::gather(LTC_TYPE, time, ARTHRITIS:CEREBROVASCULARDISEASE, factor_key=TRUE) %>% 
    #     filter(!is.na(time)) %>% 
    #     mutate(
    #       source_table = factor("SPARRALTC", levels = source_names),
    #       source_row = as.integer(row_number())
    #     ) %>%
    #     mutate(
    #       time = as_datetime(time)
    #     ) %>%
    #     numerise_columns()
    #   
    #   # Save the parsed table to compressed RData
    #   write.fst(table_SPARRALTC_old,
    #        path = filepath_clean, compress = 100
    #   )
    # } else {
    #   table_SPARRALTC_old <- as_tibble(read_fst(path = filepath_clean))
    # }
    # 
    # 
    # # Add the information to the small tables
    # patients <- patients %>%
    #   update_patients(table_SPARRALTC_old)
    # 
    # episodes <- episodes %>%
    #   update_episodes(table_SPARRALTC_old)
    # 
    # # Remove from memory
    # rm(table_SPARRALTC_old)
    # 
    # 
    
    
    # # ------------------------------------
    # SPARRALTC_recent ---------------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      str_c(str_sub(data_filenames[["SPARRALTC"]], 1, -5), "fst"))
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
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        mutate_if(is.Date, funs(as_date(format(.,"%d-%m-%Y"), tz="", format="%d-%m-%Y"))) %>% 
        tidyr::gather(LTC_TYPE, time, FIRST_ARTHRITIS_EPISODE:FIRST_RENAL_FAILURE_EPISODE, factor_key=TRUE) %>% 
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
    } else {
      table_SPARRALTC <- as_tibble(read_fst(path = filepath_clean))
    }
    
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_SPARRALTC)
    
    episodes <- episodes %>%
      update_episodes(table_SPARRALTC)
    
    # Remove from memory
    rm(table_SPARRALTC)
    
    
    
    # ------------------------------------
    # # SPARRA Scores ----------------------
    # # ------------------------------------
    # filepath_clean = file.path(
    #   dir_cleanData,
    #   str_c(str_sub(data_filenames[["SPARRA"]], 1, -5), "fst"))
    # cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    # 
    # if ((!file.exists(filepath_clean)) | force_redo==TRUE){
    #   table_SPARRA = haven::read_sav(file.path(dir_rawData, data_filenames[["SPARRA"]]))
    #   
    #   glimpse(table_SPARRA)
    #   
    #   
    #   table_SPARRA <- table_SPARRA %>% 
    #     tidyr::gather(key=period,value=SPARRAv3_score,
    #                   SPARRA_RISK_SCORE_01_04_2016:SPARRA_RISK_SCORE_01_05_2018,
    #                   na.rm=F, factor_key = F)%>% 
    #     tidyr::separate(period,c("first","second","third","day","month","year"),sep="_") %>% 
    #     mutate(time= lubridate::dmy(paste0(day,month,year,sep=""))) %>% 
    #     select(UNIQUE_STUDY_ID,time,SPARRAv3_score) %>% 
    #     rename(
    #       id = UNIQUE_STUDY_ID
    #       ) %>% 
    #     mutate(
    #       id = as.integer(id),
    #       time = as_datetime(time),
    #       SPARRAv3_score = as.integer(SPARRAv3_score)
    #     ) %>% mutate(
    #       source_table = factor("SPARRA", levels = source_names),
    #       source_row = as.integer(row_number())
    #     )
    #   
    #   # Save the parsed table to compressed RData
    #   write.fst(table_SPARRA,
    #             path = filepath_clean, compress = 100
    #   )
    # } else {
    #   table_SPARRA <- as_tibble(read_fst(path = filepath_clean))
    # }
    
    
    # ------------------------------------
    # SPARRA Scores ----------------------
    # ------------------------------------
    filepath_clean = file.path(
      dir_cleanData,
      "Final_SPARRA_extract_incl_UniqueStudyID.fst"
      )
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      # Load both SPARRAv3 score files, and merge them
      table_SPARRA_old = haven::read_sav(file.path(dir_rawData, data_filenames[["SPARRA_old"]]))
      table_SPARRA_recent = haven::read_sav(file.path(dir_rawData, data_filenames[["SPARRA_recent"]]))
      
      table_SPARRA_old <- table_SPARRA_old %>% 
        tidyr::gather(key=period,value=SPARRAv3_score,
                      SPARRA_RISK_SCORE_01_07_2013:SPARRA_RISK_SCORE_01_03_2016,
                      na.rm=F, factor_key = F) %>%
        select(-SPARRA_RISK_SCORE_01_04_2016)
      
      table_SPARRA_recent <- table_SPARRA_recent %>%
        tidyr::gather(key=period,value=SPARRAv3_score,
                      SPARRA_RISK_SCORE_01_04_2016:SPARRA_RISK_SCORE_01_05_2018,
                      na.rm=F, factor_key = F)
      
      table_SPARRA <- bind_rows(
        table_SPARRA_old,
        table_SPARRA_recent
      )
      
      
      glimpse(table_SPARRA)
      
      rm(table_SPARRA_old)
      rm(table_SPARRA_recent)
      
      table_SPARRA <- table_SPARRA %>% 
        tidyr::separate(period,c("first","second","third","day","month","year"),sep="_") %>% 
        mutate(time= lubridate::dmy(paste0(day,month,year,sep=""))) %>% 
        select(UNIQUE_STUDY_ID,time,SPARRAv3_score) %>% 
        rename(
          id = UNIQUE_STUDY_ID
        ) %>% 
        mutate(
          id = as.integer(id),
          time = as_datetime(time),
          SPARRAv3_score = as.integer(SPARRAv3_score)
        ) %>% mutate(
          source_table = factor("SPARRA", levels = source_names),
          source_row = as.integer(row_number())
        )
      
      # Save the parsed table to compressed RData
      write.fst(table_SPARRA,
                path = filepath_clean, compress = 100
      )
    } else {
      table_SPARRA <- as_tibble(read_fst(path = filepath_clean))
    }
    
    # Add the information to the small tables
    patients <- patients %>%
      update_patients(table_SPARRA)
    
    episodes <- episodes %>%
      update_episodes(table_SPARRA)
    
    # Remove from memory
    rm(table_SPARRA)
    
    
    
    # ------------------------------------
    # Patient_lookup table (new!) -------
    # ------------------------------------

    
    # # Save the "gathered" patients table, then remove source_row (to avoid conflicts with previous version)
    # write.fst(patients, 
    #           path = file.path(
    #             dir_cleanData,
    #             "patients_gathered.fst"),
    #           compress = 100)
    
    # Retain deaths (that is not in the lookup_table)
    patient_deaths <- patients %>% filter(source_table=="deaths") %>% select(id, date_of_death)
    
    rm(patients)
    
    
    # Load patients lookup table and compare to the "gathered" patients table
    # Give warnings about discrepencies, but ultimately the lookup table is accepted as the ground truth.
    
    filepath_clean = file.path(
      dir_cleanData,
      "patients.fst")
    cat(paste0("... cleaning file: ", filepath_clean, "\n"))
    
    if ((!file.exists(filepath_clean)) | force_redo==TRUE){
      table_patients = haven::read_sav(file.path(dir_rawData, data_filenames[["Patient_lookup"]]))
      
      glimpse(table_patients)
      
      table_patients <- table_patients %>% 
        rename(id:=UNIQUE_STUDY_ID ) %>% 
        mutate_if(is.character, funs(ifelse(nchar(.)==0, NA, .))) %>%
        mutate(
          id=as.integer(id),
          date_of_birth = lubridate::as_date(lubridate::fast_strptime(DOB, tz="", format="%d-%m-%Y")),
          gender=as.factor(GENDER),
          SIMD_DECILE_2016_SCT = as.integer(SIMD_DECILE_2016_SCT)
        ) %>% 
        select(-c(SIMD_QUINTILE_2016_SCT))
      
      # Investigate discrepancies
      # tmp <- patients %>% 
      #   group_by(id) %>% 
      #   # Choose smallest/earliest value in case of discrepancies
      #   summarise(
      #     sex= min(sex, na.rm = TRUE),
      #     date_of_birth = min(date_of_birth, na.rm = TRUE),
      #     date_of_death = min(date_of_death, na.rm = TRUE)
      #     ) %>%
      #   mutate(sex = 
      #          recode_factor(sex,
      #            `1` = "Male",
      #            `2` = "Female",
      #            .default = NA_character_),
      #          date_of_birth = as_date(as.integer(date_of_birth)), # Make sure NAs are true NAs
      #          date_of_death = as_date(as.integer(date_of_death))  # Make sure NAs are true NAs
      #     )
      # 
      # tmp2 <- tmp %>% 
      #   full_join(table_patients)
      # 
      # rm(tmp)
      # rm(table_patients)
      # 
      # glimpse(tmp2)
        
    
      # # Find discrepancies and report them
      # cat("PATIENTS LOOKUP TABLE vs information in individual data tables")
      # cat("--------------------------------------------------------------")
      # cat(sprintf("GENDER\n agree: %d, disagree: %d, new: %d",
      #             sum(tmp2$GENDER==tmp2$sex, na.rm=TRUE),
      #             sum(!tmp2$GENDER==tmp2$sex, na.rm=TRUE),
      #             sum(is.na(tmp2$sex) & !is.na(tmp2$GENDER), na.rm=TRUE)
      # ))
      # 
      # cat(sprintf("DATE OF BIRTH\n agree: %d, disagree: %d, new: %d",
      #             sum(tmp2$DOB==tmp2$date_of_birth, na.rm=TRUE),
      #             sum(!tmp2$DOB==tmp2$date_of_birth, na.rm=TRUE),
      #             sum(is.na(tmp2$date_of_birth) & !is.na(tmp2$DOB), na.rm=TRUE)
      # ))
      # 
      #
      # rm(tmp2)
      # # Select the appropriate columns and add back death times
      patients <- table_patients %>%
        select(id, gender, date_of_birth, SIMD_DECILE_2016_SCT) %>%
        left_join(
          patient_deaths, by="id"
        )
      
    } else {
      patients <- as_tibble(read_fst(path = filepath_clean))
    }
    
    
    
    # ----------------------------------------------------------------------
    # Save the resulting patients and episodes files
    # ----------------------------------------------------------------------
    
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

    
    ######################################################################         
    ## Generate cross-validation fold assignments by sample             ##
    ######################################################################         
    
    k=3
    set.seed(212)
    nt=nrow(patients)
    cvs=rep(1:k,1+floor(nt/k))[1:nt][order(runif(nt))] # fold assignments; equal size folds, randomised
    patients$cv=cvs
    
    ######################################################################         
    
    
        
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
              path = file.path(
                dir_cleanData,
                "episodes_lookup_by_datatable.fst"))
    
    cat(paste0("\n\nDATA CLEANING FINISHED\n----------------------\n"))
    
    # Create a version of episodes without SPARRAv3 scores as well, it is a lot of data and we generally won't use it as features
    
    # episodes_no_sparrav3 <- episodes %>% filter(! (source_table %in% c("SPARRA")))
    # write.fst(episodes_no_sparrav3,
    #           path = file.path(
    #             dir_cleanData,
    #             "episodes_no_sparrav3.fst"))
    
    
  
  }
  
}

