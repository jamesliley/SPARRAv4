#' Returns all diagnosis and prescription codes for all patients as a list column
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_all_codes_docterm_sparsematrix <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  min_total_codes = 1,
  min_unique_codes = 1,
  return_extra_codes = FALSE,
  all_unique_codes = NULL, # Supply a pre-existing set of condition codes to use as column names 
  table_of_all_codes_long = NULL # One may supply a pre-computed table of diagnosis codes
  
){
  
  if (is.null(table_of_all_codes_long)){
    
    # This function is provided in "James/Transformation/primitive_function/primitive_diagnosis_codes_count.R" folder
    table_of_diag_codes <- primitive_diagnosis_codes_count(
      patients, episodes, list_of_data_tables
    ) 
    gc()
    
    print("Diagnosis code count complete")
    
    # This function is provided in "James/Transformation/primitive_function/primitive_prescription_count.R" folder
    table_of_prescription_codes <- primitive_prescriptions_count(
      patients, episodes, list_of_data_tables
    ) 
    gc()
    
    print("Prescription code count complete")
    
    # Combine tables
    table_of_diag_codes <- table_of_diag_codes %>% 
      rename(all_codes=all_condition_codes,count_codes=count_condition_codes,count_unique_codes=count_unique_condition_codes)
    table_of_prescription_codes <- table_of_prescription_codes %>% 
      rename(all_codes=all_prescription_codes,count_codes=count_prescription_codes,count_unique_codes=count_unique_prescription_codes)
    table_of_all_codes <- rbind(table_of_prescription_codes,table_of_diag_codes)
    
    # Clean up
    rm(table_of_prescription_codes)
    rm(table_of_diag_codes)
    gc()
    
    print("Table of all codes generated")
    
    # Turn into a long table format
    table_of_all_codes_long <- primitive_base_R_longtable(table_of_all_codes)
    
  }
  
  print("Long table generated")
  
  # If all_unique_codes are not given, compute from current data, 
  # otherwise create matrix from the supplied one (presumably that of training data)
  if (is.null(all_unique_codes)){
    all_unique_codes <- unique(unlist(c(table_of_all_codes_long$all_codes)))
  }
  
  print("All unique codes generated")
  
  # Create table of extra codes
  if (return_extra_codes==TRUE){
    table_of_extra_codes <- table_of_all_codes_long %>%
      ungroup() %>% 
      filter(
        !(all_codes %in% all_unique_codes)
      )
  }
  
  print("Table of extra codes generated")
  
  # Filter table (first remove rows unknown condition codes, then make sure the remaining is still in accordance with given filtering)
  tmp <- table_of_all_codes_long %>%
    ungroup() %>%
    filter(
      (all_codes %in% all_unique_codes)
    ) %>%
    filter(
      count_codes >= min_total_codes & count_unique_codes >= min_unique_codes
    ) %>% select(
      -count_codes, -count_unique_codes
    )
  
  print("Table filtered")
  
  # Get the remaining unique ids to use as row names later
  all_ids = unique(tmp$id)
  
  tmp <- tmp %>%
    # First convert id and condition codes to factors, 
    # then pass factor levels to sparseMatrix constructor
    mutate(
      id = factor(id, levels = all_ids),
      all_codes = factor(
        all_codes, levels = all_unique_codes)
    ) %>% filter(
      !is.na(all_codes)
    )
  
  
  print("Ready to produce sparse matrix")
  
  library(Matrix)
  doc_term_sparse_matrix <- sparseMatrix(i = as.integer(tmp$id), 
                                         j = as.integer(tmp$all_codes), 
                                         x = as.numeric(tmp$value),
                                         dims = c(max(as.integer(tmp$id)), length(all_unique_codes))
  )
  
  rownames(doc_term_sparse_matrix) <- all_ids[1:nrow(doc_term_sparse_matrix)]
  colnames(doc_term_sparse_matrix) <- all_unique_codes[1:ncol(doc_term_sparse_matrix)]
  
  print("Produced sparse matrix")
  
  if (return_extra_codes){
    # Return the sparse matrix + extra codes
    out <- list(
      doc_term_sparse_matrix = doc_term_sparse_matrix,
      table_of_extra_codes = table_of_extra_codes
    )
  } else {
    
    # Return the sparse doc_term_matrix
    out <- doc_term_sparse_matrix
  }
  
  
  # Return
  out
}




# This function is a rewrite in base R of the dplyr command
#
#  doc_term_long_table_training <- table_of_all_codes %>%
#    unnest(cols=all_codes) %>%
#    group_by(id, count_codes, count_unique_codes, all_codes) %>%
#    summarise(value = as.integer(n()))
#
# If verbose is set, print progress every 10k rows.

primitive_base_R_longtable=function(tab,verbose=T) {
  
  
  xall_codes=unlist(unlist(tab$all_codes))
  xid=rep(tab$id,tab$count_codes)
  xcount_codes=rep(tab$count_codes,tab$count_codes)
  xcount_unique_codes=rep(tab$count_unique_codes,tab$count_codes)
  
  newtab=as_tibble(data.frame(id=xid,all_codes=xall_codes,
                              count_codes=xcount_codes,
                              count_unique_codes=xcount_unique_codes,
                              stringsAsFactors=FALSE))
  
  rm(xall_codes)
  rm(xid)
  rm(xcount_codes)
  rm(xcount_unique_codes)
  gc()
  
  if (verbose) print("Unnested")
  
  newtab <- newtab %>% group_by(id, count_codes, count_unique_codes, all_codes) 
  newtab=newtab[order(newtab$id),]
  
  if (verbose) print("Grouped")
  
  nx=newtab %>% summarise(value = as.integer(n()))
  
  # Add back in any samples with no codes
  w0=tab$id[which(tab$count_codes==0)]
  w1=tab$id[which(tab$count_codes>0)]
  if (length(w0)>0) {
    wx=setdiff(w0,w1)
    tx=nx[1:length(wx),]
    tx[[1]]=wx
    tx[[2]]=0
    tx[[3]]=0
    tx[[4]]=""
    tx[[5]]=0
    
    nx=rbind(nx,tx)
  } 
  
  nx=nx[order(nx$id),]
  
  
  return(nx)
  
}


# Test times
#xx=10*c(1000,2000,3000,4000)
#xt=c()
#
#for (i in 1:length(xx)) {
#    tabx=table_of_all_codes[sample(dim(table_of_all_codes)[1],xx[i]),]
#    
#xt=rbind(xt,system.time(dt<- expand_tab(tabx))[1:3])
#print(i)    
#}#' Returns a sparse matrix (to used for topic modelling) given a code table:
#' 
#' with columns "id", "code", "value"
#' 
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_codes_docterm_sparsematrix <- function(
  # Input data
  table_of_codes_long,
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  min_total_codes = 1,
  min_unique_codes = 1,
  max_unique_codes = Inf,
  return_extra_codes = FALSE,
  all_unique_codes = NULL # Supply a pre-existing set of condition codes to use as column names 
){
  
  
  
  # If all_unique_codes are not given, compute from current data, 
  # otherwise create matrix from the supplied one (presumably that of training data)
  if (is.null(all_unique_codes)){
    all_unique_codes <- unique(unlist(c(table_of_codes_long$all_codes)))
  }
  
  
  
  # Create table of extra codes
  if (return_extra_codes==TRUE){
    table_of_extra_codes <- table_of_codes_long %>%
      ungroup() %>% 
      filter(
        !(all_codes %in% all_unique_codes)
      )
  }
  
  # Create smaller table that computes the total and unique code requirements
  tmp_summary <- table_of_codes_long %>% 
    group_by(id) %>%
    summarise(total_codes = sum(value), num_unique_codes = n())
  
  # Filter table (first remove rows unknown condition codes, then make sure the remaining is still in accordance with given filtering)
  tmp <- table_of_codes_long %>%
    ungroup() %>%
    left_join(tmp_summary) %>%
    filter(
      (all_codes %in% all_unique_codes)
    ) %>%
    filter(
      count_codes >= min_total_codes & count_unique_codes >= min_unique_codes & count_unique_codes <= max_unique_codes 
    ) %>% select(
      -count_codes, -count_unique_codes
    )
  
  # Get the remaining unique ids to use as row names later
  all_ids = unique(tmp$id)
  
  tmp <- tmp %>%
    # First convert id and condition codes to factors, 
    # then pass factor levels to sparseMatrix constructor
    mutate(
      id = factor(id, levels = all_ids),
      code = factor(
        all_codes, levels = all_unique_codes)
    ) %>% filter(
      !is.na(code)
    )
  
  
  library(Matrix)
  doc_term_sparse_matrix <- sparseMatrix(i = as.integer(tmp$id), 
                                         j = as.integer(tmp$code), 
                                         x = as.numeric(tmp$value),
                                         dims = c(max(as.integer(tmp$id)), length(all_unique_codes))
  )
  
  rownames(doc_term_sparse_matrix) <- all_ids[1:nrow(doc_term_sparse_matrix)]
  colnames(doc_term_sparse_matrix) <- all_unique_codes[1:ncol(doc_term_sparse_matrix)]
  
  if (return_extra_codes){
    # Return the sparse matrix + extra codes
    out <- list(
      doc_term_sparse_matrix = doc_term_sparse_matrix,
      table_of_extra_codes = table_of_extra_codes
    )
  } else {
    
    # Return the sparse doc_term_matrix
    out <- doc_term_sparse_matrix
  }
  
  
  # Return
  out
}



#' Returns all diagnosis codes for all patients as a list column
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_diagnosis_codes_count_simple <- function(
  # Input data
  list_of_data_tables,
  
  time_cutoff=dmy_hm("2-5-2100 00:00"), # Filter to only times before this
  time_min=dmy_hm("2-5-1900 00:00") # Filter to only times after this
) {
  
  # Extract
  SMR01M = list_of_data_tables$SMR01M
  
  # Filter
  SMR01M = SMR01M[which(SMR01M$time < time_cutoff &
                          SMR01M$time >= time_min), ]
  
  # List of ICD10 codes
  out = SMR01M[, c("id", "all_stay_code_list", "time", "time_discharge","cis_marker","source_table")]
  
  # Remove NAs
  out$all_stay_code_list = lapply(out$all_stay_code_list, function(x)
    x[which(!is.na(x) & length(x)>0)])
  
  # Fix ICD10 codes with spaces
  out$all_stay_code_list = lapply(out$all_stay_code_list, function(x) 
    gsub(" ","_",x))
  
  # Rename lists of condition codes
  names(out$all_stay_code_list) = rep(NULL, dim(out)[1])
  
  print("Extracted all_stay_code_list")
  
  # Fix records generated after time cutoff
  w=which(is.na(out$time_discharge) | (out$time_discharge >= time_cutoff))
  sn=c("main_condition",paste0("OTHER_CONDITION_",1:5)) # Pertinent column names in SMR01
  sen=c("main_condition",paste0("OTHER_CONDITION_",1:5)) # Pertinent column names in SMR01E
  swn=c("main_condition",paste0("DIAGNOSIS_",2:6)) # Pertinent column names in SystemWatch
  SMR01 = list_of_data_tables$SMR01
  SMR01E = list_of_data_tables$SMR01E
  SystemWatch = list_of_data_tables$SystemWatch
  SMR01=SMR01[which(SMR01$id %in% out$id[w]),]
  SMR01E=SMR01E[which(SMR01E$id %in% out$id[w]),]
  SystemWatch=SystemWatch[which(SystemWatch$id %in% out$id[w]),]
  for (i in 1:length(w)) {
    if ((i %% 500)==0) print(paste0("Completed ",i," of ",length(w)))
    s1=SMR01[which(
      SMR01$id==out$id[w[i]] & 
        SMR01$CIS_MARKER==out$cis_marker[w[i]] &
        SMR01$time < time_cutoff &
        out$source_table[w[i]] =="SMR01"),sn]
    s1e=SMR01E[which(
      (SMR01E$id==out$id[w[i]] &
         SMR01E$CIS_MARKER==out$cis_marker[w[i]] &
         SMR01E$time < time_cutoff &
         out$source_table[w[i]] =="SMR01E") |
        (SMR01E$id==out$id[w[i]] &
           SMR01E$time==out$time[w[i]] &
           SMR01E$time < time_cutoff) &
        out$source_table[w[i]] =="SMR01E"),sen]
    sw1=SystemWatch[which(
      SystemWatch$id==out$id[w[i]] & 
        SystemWatch$time==out$time[w[i]] & 
        SystemWatch$time < time_cutoff &
        out$source_table[w[i]] =="SystemWatch"),swn]
    xlist=c(unlist(s1),unlist(s1e),unlist(sw1))
    if (length(xlist)>0) {
      out$all_stay_code_list[w[i]][[1]]=unique(xlist[which(!is.na(xlist))]) 
    } else out$all_stay_code_list[w[i]][[1]]=c("")
  }
  out=out[which(!unlist(lapply(out$all_stay_code_list,function(x) x[1]==""))),]
  out=out[,1:2]
  
  # Remove NAs
  out$all_stay_code_list = lapply(out$all_stay_code_list, function(x)
    x[which(x!= "")])
  
  # Combine by ID 
  out = out %>% 
    group_by(id) %>% 
    unnest(all_stay_code_list) %>% 
    summarise(all_condition_codes=list(all_stay_code_list))
  
  # Count
  out$count_condition_codes = unlist(lapply(out$all_condition_codes, function(x)
    length(x)))
  out$count_unique_condition_codes = unlist(lapply(out$all_condition_codes, function(x)
    length(unique(x))))
  
  
  # Trim to non-zero
  out = out[which(out$count_condition_codes > 0), ]
  
  return(out)
}


#' Returns all diagnosis codes for all patients as a list column
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_diagnosis_codes_count <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  source_table_names = c("SMR00", "SMR01", "SMR01E", "SMR04"),
  
  time_cutoff=dmy_hm("2-5-2100 00:00"), # Filter to only times before this
  time_min=dmy_hm("2-5-1900 00:00") # Filter to only times after this
){
  
  # Collect relevant column names that contain condition (ICD10 code)
  column_names_flat <- c()
  
  for (cur_table in c("SMR01M", "SMR00","SMR04")) { #, "SMR01", "SMR01E", "SMR04")){
    if (cur_table %in% source_table_names){
      column_names_flat <- c(
        column_names_flat,
        paste0(cur_table, ".", # dput(names(list_of_data_tables$SMR00 %>% select(contains("CONDITION"))))
               names(list_of_data_tables[[cur_table]] %>% select(contains("CONDITION"))))
      ) 
    }
    
  }
  
  # Change ""'s to NAs
  for (cur_table in c("SMR01M", "SMR00","SMR04")) { #, "SMR01", "SMR01E", "SMR04")){
    if (cur_table %in% source_table_names){
      cx=names(list_of_data_tables[[cur_table]] %>% select(contains("CONDITION")))
      for (i in 1:length(cx)) {
        w=which(list_of_data_tables[[cur_table]][,cx[i]]=="")
        list_of_data_tables[[cur_table]][w,cx[i]]=NA
      }
    }
    
  }
  
  # This is all condition codes for all people as lists in columns over many columns
  lists_elementwise <- collect_and_full_join_from_list_of_tables(
    list_of_data_tables,
    column_names_flat,
    table_key = "id",
    time_cutoff=time_cutoff,
    time_min=time_min
  )
  
  # Gather all codes into a single character list
  #lists_elementwise
  table_of_diag_codes <- lists_elementwise %>% 
    gather(key="key", value = "value", 2:ncol(.)) %>% 
    group_by(id) %>% 
    summarise(all_condition_codes = list(unlist(value)[!is.na(unlist(value))]))
  
  # Now we have an ID and a vector of all diagnosis codes for the patient
  
  table_of_diag_codes <- table_of_diag_codes %>%
    mutate(
      count_condition_codes =  unlist(map(.x = all_condition_codes, ~ length(.))),
      count_unique_condition_codes = unlist(map(.x = all_condition_codes, ~ length(unique(.))))
    )
  
  # We get to do various statistics on the diagnosis codes
  #tmp2 <- tmp %>% mutate(allcodes2 = map(.x = allcodes, ~ .[[1]][length(.[[1]])]))
  
  # Return the diagnosis codes table
  return(table_of_diag_codes)
  
}




#' Returns all diagnosis codes for all patients as a list column
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_diagnosis_codes_docterm_sparsematrix <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  min_total_codes = 1,
  min_unique_codes = 1,
  return_extra_codes = FALSE,
  all_unique_codes = NULL, # Supply a pre-existing set of condition codes to use as column names 
  table_of_diag_codes_long = NULL # One may supply a pre-computed table of diagnosis codes
  
){
  
  if (is.null(table_of_diag_codes_long)){
    table_of_diag_codes <- primitive_diagnosis_codes_count(
      patients, episodes, list_of_data_tables
    ) # 20 min 
    
    # Convert it to a sparse matrix
    #table_of_diag_codes_long <- table_of_diag_codes %>% 
    #  unnest() %>%
    #  group_by(id, count_condition_codes, count_unique_condition_codes, all_condition_codes) %>%
    #  summarise(value = as.integer(n())) # 5 min
    
    # dplyr's 'unnest' is slow in this release (known problem) so using the following workaround for the moment.
    
    tab=table_of_diag_codes
    
    xall_codes=unlist(unlist(tab$all_condition_codes))
    xid=rep(tab$id,tab$count_condition_codes)
    xcount_codes=rep(tab$count_condition_codes,tab$count_condition_codes)
    xcount_unique_codes=rep(tab$count_unique_condition_codes,tab$count_condition_codes)
    
    newtab=as_tibble(data.frame(id=xid,all_condition_codes=xall_codes,
                                count_condition_codes=xcount_codes,
                                count_unique_condition_codes=xcount_unique_codes,
                                stringsAsFactors=FALSE))
    
    table_of_diag_codes_long <- newtab %>% 
      group_by(id, count_condition_codes, count_unique_condition_codes, all_condition_codes) %>%
      summarise(value = as.integer(n()))
    
  }
  
  print("Long table generated")
  
  # If all_unique_codes are not given, compute from current data, 
  # otherwise create matrix from the supplied one (presumably that of training data)
  if (is.null(all_unique_codes)){
    all_unique_codes <- unique(unlist(c(table_of_diag_codes_long$all_condition_codes)))
  }
  
  print("All unique codes generated")
  
  # Create table of extra codes
  if (return_extra_codes==TRUE){
    table_of_extra_codes <- table_of_diag_codes_long %>%
      ungroup() %>% 
      filter(
        !(all_condition_codes %in% all_unique_codes)
      )
  }
  
  print("Table of extra codes generated")
  
  # Filter table (first remove rows unknown condition codes, then make sure the remaining is still in accordance with given filtering)
  tmp <- table_of_diag_codes_long %>%
    ungroup() %>%
    filter(
      (all_condition_codes %in% all_unique_codes)
    ) %>%
    filter(
      count_condition_codes >= min_total_codes & count_unique_condition_codes >= min_unique_codes
    ) %>% select(
      -count_condition_codes, -count_unique_condition_codes
    )
  
  print("Table filtered")
  
  # Get the remaining unique ids to use as row names later
  all_ids = unique(tmp$id)
  
  tmp <- tmp %>%
    # First convert id and condition codes to factors, 
    # then pass factor levels to sparseMatrix constructor
    mutate(
      id = factor(id, levels = all_ids),
      all_condition_codes = factor(
        all_condition_codes, levels = all_unique_codes)
    ) %>% filter(
      !is.na(all_condition_codes)
    )
  
  
  print("Ready to produce sparse matrix")
  
  library(Matrix)
  doc_term_sparse_matrix <- sparseMatrix(i = as.integer(tmp$id), 
                                         j = as.integer(tmp$all_condition_codes), 
                                         x = as.numeric(tmp$value),
                                         dims = c(max(as.integer(tmp$id)), length(all_unique_codes))
  )
  
  rownames(doc_term_sparse_matrix) <- all_ids[1:nrow(doc_term_sparse_matrix)]
  colnames(doc_term_sparse_matrix) <- all_unique_codes[1:ncol(doc_term_sparse_matrix)]
  
  print("Produced sparse matrix")
  
  if (return_extra_codes){
    # Return the sparse matrix + extra codes
    out <- list(
      doc_term_sparse_matrix = doc_term_sparse_matrix,
      table_of_extra_codes = table_of_extra_codes
    )
  } else {
    
    # Return the sparse doc_term_matrix
    out <- doc_term_sparse_matrix
  }
  
  
  # Return
  return(out)
}




#' Returns all prescriptions for all patients as a long table
#' furthermore returns the number of condition codes and number of unique cond. codes for each patient
#' 
primitive_prescriptions_count <- function(
  # Input data
  patients, episodes, list_of_data_tables, 
  
  # Further required parameters
  
  # Further optional parameters (these have defaults)
  keep_time = FALSE,
  time_cutoff=dmy_hm("2-5-2100 00:00"), # Filter to only times before this
  time_min=dmy_hm("2-5-1900 00:00") # Filter to only times after this
){
  
  
  bnf_char=as.character(list_of_data_tables$PIS$bnf_section)[
    which(list_of_data_tables$PIS$time<time_cutoff & 
            list_of_data_tables$PIS$time >= time_min)] # pre-emptively do factor->character conversion, using R's stored list of factor labels
  table_of_prescription_codes= list_of_data_tables$PIS %>% 
    filter(time<time_cutoff) %>%
    filter(time>= time_min) %>%
    mutate(bnf_char=bnf_char) %>%
    group_by(id) %>%
    summarise(all_prescription_codes = list(bnf_char[!is.na(bnf_char)]))
  
  # above takes 1 min
  
  table_of_prescription_codes <- table_of_prescription_codes %>%
    mutate(
      count_prescription_codes =  unlist(map(.x = all_prescription_codes, ~ length(.))),
      count_unique_prescription_codes = unlist(map(.x = all_prescription_codes, ~ length(unique(.))))
    )
  
  # Return the prescription codes table
  return(table_of_prescription_codes)
  
}





