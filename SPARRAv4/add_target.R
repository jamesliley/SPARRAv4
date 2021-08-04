# Add the target column to a long format tibble with "id" and "time" columns

library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)

#source("James/Transformation/get_episodes_within_time.R")
#source("James/Util/basic_operations.R")

#' add_target
#' This adds the target column to a long format tibble with "id" and "time" columns, based on information in an "episodes" table evaluated by an "event_func"
#' 
#' Input:
#' @param data a tibble that have at least "id" <integer> and "time" <dttm> columns. The "time" column should be the same for all ids
#' @param episodes a tibble that is used as input to event_func
#' @param list_of_data_tables The list of data tables used as input to event_func
#' @param event_func A function handle that takes the "episodes" table and the "list_of_data_tables" as input and and returns a tibble with "id", "time", and "event" <binary>
#' @param prediction_period = 365, # days, length of the outcome aggregation (also: censors this much at the end of the dataset)
#' @param fill_times = FALSE, Whether or not have all unique times in data filled for each patient (even if some of them missed some rows)
#' @param target_colname Name of the added target column
#' 
#' Output:
#' @return The data tibble with an extra binary column named "target" <binary>. The value in "target":
#'       It is TRUE if for a given id and time, the same id had a true "event" (event_func returns positive) within "prediction_period"
#'       It is NA if the given time is within "prediction_period" to the last supplied time in "episodes", indicating we do not have data to determine the target,
#'       It is FALSE otherwise
#' 
#' The function also saves this table to the disk within the intermediate data folder
add_target <- function(
  data,
  episodes,
  list_of_data_tables,
  event_func = event_func_v4,
  prediction_period = 365, # days, length of the outcome aggregation (also: censor this much at the end of the dataset)
  fill_times = FALSE,
  target_colname="target"
){
  
  # Create emergency admission event_func if NULL
  if (is.null(event_func)){
    event_func=event_func_v4
    # Complicated event func, including "urgent" admissions and taking "ADMISSION_TRANSFER" into account
    # See complicated definition of emergency admission on sheet 1 of "\\Farr-FS1\Study Data\1718-0370\Linked Data\Upload to SH 4\admissions and prescribing (for Gergo) RP.xlsx"
    
    # # Simple event_func (according to online data dictionary)
    # event_func <- function(ep){
    #   ep %>% transmute(
    #     id = id,
    #     time = time,
    #     event = (!is.na(admission_type) & admission_type >= 30 & admission_type < 40) 
    #   )
    # }
  }
  
  
  # ------------------------------
  # Get the true events to the episodes table
  # As we'll add evaluation times, we only need events that are TRUE to propagate back
  # This tends to be memory-intensive, so do separately for bits of 'episodes'
  
  dim_e=dim(episodes)[1]; n_e=floor(dim_e/1e7)
  
  true_event_times=c()
  for (i_ep in 1:n_e) {
  true_event_timesx <- episodes[(1+(1e7)*(i_ep-1)):min((1e7)*i_ep,dim_e),] %>%
    ungroup() %>%
    filter(id %in% data$id) %>%
    event_func(., list_of_data_tables) %>%
    filter(event) %>%
    distinct() # As actual admission tables and SystemWatch sometimes contains exact duplicates based on just "id" and "time", we want to keep only distinct records here 
  true_event_times=rbind(true_event_times,true_event_timesx)
  gc()
  print(paste0("Processed event times ",i_ep," of ",n_e))
  }
  true_event_times<- true_event_times %>% distinct() 
  true_event_times$time[which(is.na(true_event_times$time))]=as_datetime("2000-01-01") # If admission or death date is unknown, set to an early timepoint so as to remove from dataset (set target to NA).
  
  print("Determined true event times")
  
  # Get death dates
  death_dates=c()
  for (i_ep in 1:n_e) {
  death_datesx=episodes %>%
    ungroup() %>%
    filter(id %in% data$id) %>%
    mutate(death=(source_table=="deaths")) %>%
    filter(death) %>%
    distinct()
  death_dates=rbind(death_dates,death_datesx)
  print(paste0("Processed death dates ",i_ep," of ",n_e))
  }
  death_dates=death_dates %>% distinct()

  print("Determined death dates")
    
  # Get all the times we need 
  all_times <- data %>%
    select(time) %>%
    distinct()
  
  print("Attained times")
  
  # If we are only predicting at a few unique times, or want to "fill_times", use expand and filter, otherwise propagate individually
  if ((nrow(all_times) < 50) || fill_times) {
  
    # Compute the target at all times for all people with true events
    target_per_patient_per_time <- true_event_times %>%
      select(id, time, event) %>%
      bind_rows(
        tidyr::expand(
          tibble(time = all_times$time, source_table="TMP_TIMES"),
          id = .$id,
          nesting(time, source_table)
        )
      ) %>% 
      mutate(event = ifelse(is.na(event), FALSE, event)) %>% # Make sure we don't have NA events
      add_mutate_timeUNIX() %>%
      arrange(id, timeUNIX) %>%
      group_by(id) %>%
      mutate(
        #   find all times that are associated with an event within a time window,
        #   and add target = 1. (If a time is not within such a time window, set target= 0 )
        target = get_episodes_within_time(timeUNIX, event, max_days = prediction_period)
      ) %>% 
      ungroup() %>%
      select(id, time, source_table, target) %>%
      filter(source_table=="TMP_TIMES") %>% 
      select(-c(source_table)) %>% 
      # Add back all the people who had no events whatsoever
      bind_rows(
        tidyr::expand(
          tibble(time = all_times$time, target=FALSE),
          id = unique(data$id)[!unique(data$id) %in% unique(.$id)],
          nesting(time, target)
        )
      )
    print("Generated targets per patient per time")
    
    if (!fill_times){
      out <- data %>%
        left_join(target_per_patient_per_time)
    } else {
      out <- data %>%
        full_join(target_per_patient_per_time)
    }
    
    
  } else { # We have lots of unique times in the "data" and thus "all_times", and we want to avoid storing predictions at each unique time, even temporarily.
  
    # Compute the target for all people with true events, but only at their own times
    target_per_patient <- true_event_times %>%
      full_join(
        data %>% select(id, time) %>% filter(id %in% unique(true_event_times$id)) %>% distinct(),
        by = c("id", "time")
      ) %>%
      mutate(event = ifelse(is.na(event), FALSE, event)) %>% # Make sure we don't have NA events
      add_mutate_timeUNIX() %>%
      arrange(id, timeUNIX) %>%
      group_by(id) %>%
      mutate(
        #   find all times that are associated with an event within a time window,
        #   and add target = 1. (If a time is not within such a time window, set target= 0 )
        target = get_episodes_within_time(timeUNIX, event, max_days = prediction_period)
      ) %>% 
      ungroup() %>%
      select(id, time, target)
      
    print("Generated target per patient")  
    
    
    # For all the people and times who did not have any true events, just set the target to FALSE at all times
    out <- data %>%
      left_join(target_per_patient) %>%
      mutate(target = ifelse(is.na(target), FALSE, target)) # Make sure we don't have NA targets
    
    # Check what is the time of the last episode in our dataset, and set targets to NA that are within the prediction period to the last episode
    last_episode <- max(episodes$time)
    out <- out %>% 
      mutate_cond( # This function is from Gergo/Util/basic_operations.R
        condition = (as.integer(time) > (as.integer(max(episodes$time)) - prediction_period*24*60*60 )), # prediction period, convert from days to seconds
        target = NA
        ) 
  }
  
  
  print("Finalising results")
  
  # Rename column
  colnames(out)[which(colnames(out)=="target")]=target_colname
  
  # Remove individuals who died before the time cutoff
  death_date_by_id=death_dates$time[match(out$id,death_dates$id)] # mostly NA
  out$target[which(out$time>death_date_by_id)]=NA
  
  # Return the augmented tibble
  out
  
}


##' event_func_v4()
##'
##' Determines events as per SPARRAv4 definition
##'
##' @param ep episodes tibble; output from loadCleanData
##' @param list_of_data_tables list of data tables; output from loadCleanData
##' @return tibble with "id", "time", and "event" <binary>
event_func_v4 <- function(ep, list_of_data_tables){
  ep_tmp <- ep %>% 
    mutate(
      event = ( (!is.na(admission_type) & admission_type %in% c(20,21,22,30,31,32,33,34,35,36,38,39)) | # get simple event definition first
                (source_table=="deaths")) # death in outcome period counts as an admission
    )
  
  # Fix edge case events (admission type == 38),  we need to check ADMISSION_TRANSFER_FROM in SMR01 and SMR01E
  ep_tmp_edgecases <- ep_tmp %>%
    # We only need the edge-cases to fix (admission type is 38)
    filter(admission_type == 38) %>%
    select(id, source_table, source_row, event) %>%
    
    # First sort out SMR01 admission transfer
    left_join(
      list_of_data_tables$SMR01 %>% select(id, source_row, source_table, ADMISSION_TRANSFER_FROM),
      by = c("id", "source_row", "source_table")
    ) %>% 
    mutate(
      event = (event & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
    ) %>%
    select(-c("ADMISSION_TRANSFER_FROM")) %>%
    
    # Now sort out SMR01E admission transfer
    left_join(
      list_of_data_tables$SMR01E %>% select(id, source_row, source_table, ADMISSION_TRANSFER_FROM),
      by = c("id", "source_row", "source_table")
    ) %>% 
    mutate(
      event = (event & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
    ) %>%
    select(-c("ADMISSION_TRANSFER_FROM"))
  
  ep_tmp <- ep_tmp %>% 
    left_join(ep_tmp_edgecases, by = c("id", "source_table", "source_row")) %>% 
    transmute(
      id = id,
      time = time,
      event= ifelse(!is.na(event.y), event.y, event.x)
    )
  
  rm(ep_tmp_edgecases)
  
  # Return
  ep_tmp
}




# Currently unused parameters
#' @param dir_results Directory to store the output tibble
#' @param file_results File name to store the output tibble
#' @param force_redo Binary flag, indicates whether the process should be carried out anew, or one may just load the latest stored intermediate version
#' 

#' get_event_func_icd
#' ICD-code specific event function
#' 
get_event_func_icd = function(icd_start="S") {

icdchar=nchar(icd_start[1])

out_function <- function(ep, list_of_data_tables){
      ep_tmp <- ep %>% 
        mutate(
          event = (!is.na(admission_type) & admission_type %in% c(20,21,22,30,31,32,33,34,35,36,38,39)) & # get simple event definition first
                  (!is.na(main_condition) & substring(main_condition,1,icdchar) %in% icd_start) 
        )
      
      # Fix edge case events (admission type == 38),  we need to check ADMISSION_TRANSFER_FROM in SMR01 and SMR01E
      ep_tmp_edgecases <- ep_tmp %>%
        # We only need the edge-cases to fix (admission type is 38)
        filter(admission_type == 38) %>%
        select(id, source_table, source_row, event) %>%
        
        # First sort out SMR01 admission transfer
        left_join(
          list_of_data_tables$SMR01 %>% select(id, source_row, source_table, ADMISSION_TRANSFER_FROM),
          by = c("id", "source_row", "source_table")
        ) %>% 
        mutate(
          event = (event & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
          ) %>%
        select(-c("ADMISSION_TRANSFER_FROM")) %>%
        
        # Now sort out SMR01E admission transfer
        left_join(
          list_of_data_tables$SMR01E %>% select(id, source_row, source_table, ADMISSION_TRANSFER_FROM),
          by = c("id", "source_row", "source_table")
        ) %>% 
        mutate(
          event = (event & !str_sub(ADMISSION_TRANSFER_FROM,1,1) %in% c("4","5"))
        ) %>%
        select(-c("ADMISSION_TRANSFER_FROM"))
      
      ep_tmp <- ep_tmp %>% 
        left_join(ep_tmp_edgecases, by = c("id", "source_table", "source_row")) %>% 
        transmute(
          id = id,
          time = time,
          event= ifelse(!is.na(event.y), event.y, event.x)
          )
      
      rm(ep_tmp_edgecases)
      
      # Return
      ep_tmp
    }
   out_function
}

 