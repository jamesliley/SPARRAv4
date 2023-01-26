######################################################
## Function to fit topic models to long-term        ##
##  table of data                                   ##
######################################################
##
## James Liley, 2020. Largely copied from original script by Gergo Bohner, 2018-19
##




######################################################
## Create table for topic model fitting             ##
######################################################


##' Function to produce table to which topic model can be fitted, given objects patients, episodes, and list_of_data_tables from load_cleandata.R.
##' @param patients patients table, output from load_cleanData
##' @param episodes episodes table, output from load_cleanData
##' @param list_of_data_tables, output from load_cleanData
##' @param time_cutoff optional time cutoff(s); only use records from before this date. Default 2100AD.
##' @param time_min optional time cutoff(s); only use records from after this date. Default 1900AD.
##' @param save_file set this to save the episode as ${save_file}.rds
##' @return a matrix of prescription and diagnosis codes
topic_training_matrix=function(patients,episodes,list_of_data_tables,
                               time_cutoff=dmy_hm("2-5-2100 00:00"),
                               time_min=dmy_hm("2-5-1900 00:00"),
                               save_file=NULL) {
  
  
  # Initialise
  doc_term_sparse_matrix=c()    
  
  for (time_i in 1:length(time_cutoff)) {
    
    # Extract diagnostic codes
    table_of_diag_codes = primitive_diagnosis_codes_count_simple(
      list_of_data_tables, 
      time_cutoff=time_cutoff[time_i],
      time_min=time_min[time_i])
    gc()
    print("Generated table of diagnosis codes")
    
    # Extract prescription codes
    table_of_prescription_codes = primitive_prescriptions_count(
      patients, episodes, list_of_data_tables,
      time_cutoff=time_cutoff[time_i]-(30*24*3600), 
      time_min=time_min[time_i]
    ) 
    gc()
    print("Generated table of prescription codes")
    
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
    
    # Turn into a long table format
    # dplyr is slow here - formerly, this read 
    #doc_term_long_table_training <- table_of_all_codes %>%
    #  unnest(cols=all_codes) %>%
    #  group_by(id, count_codes, count_unique_codes, all_codes) %>%
    #  summarise(value = as.integer(n()))
    #
    doc_term_long_table_training <- primitive_base_R_longtable(table_of_all_codes)
    
    
    # Get a matrix to train on
    # This function is provided in "James/Transformation/primitive_function/primitive_codes_docterm_sparsematrix.R" folder
    doc_term_sparse_matrix0 = primitive_codes_docterm_sparsematrix(doc_term_long_table_training,
                                                                   min_total_codes = -1,
                                                                   min_unique_codes = -1,
                                                                   max_unique_codes = Inf,
                                                                   return_extra_codes = FALSE,
                                                                   all_unique_codes = NULL # Supply a pre-existing set of condition codes to use as column names 
    )
    
    # Trim to patients in database
    doc_term_sparse_matrix0=doc_term_sparse_matrix0[which(rownames(doc_term_sparse_matrix0) %in% patients$id),]
    rownames(doc_term_sparse_matrix0)=paste0(rownames(doc_term_sparse_matrix0),"_",as.numeric(time_cutoff[time_i]))
    
    # Concatenate by time cutoff. Need to do rbind.fill manually.
    
    if (time_i==1) doc_term_sparse_matrix=doc_term_sparse_matrix0 else {
      d1=colnames(doc_term_sparse_matrix)
      d0=colnames(doc_term_sparse_matrix0)
      m1x=as(matrix(0,dim(doc_term_sparse_matrix)[1],length(setdiff(d0,d1))),"dgCMatrix")
      m0x=as(matrix(0,dim(doc_term_sparse_matrix0)[1],length(setdiff(d1,d0))),"dgCMatrix")
      m1=cbind(doc_term_sparse_matrix,m1x)
      m0=cbind(doc_term_sparse_matrix0,m0x)
      colnames(m1)=c(colnames(doc_term_sparse_matrix),setdiff(d0,d1))
      colnames(m0)=c(colnames(doc_term_sparse_matrix0),setdiff(d1,d0))
      dx=union(d0,d1)
      m0=m0[,dx]; m1=m1[,dx]
      doc_term_sparse_matrix=rbind(m1,m0)
      rm(list=c("m0","m1","doc_term_sparse_matrix0"))
    }
    
  }
  
  if (!is.null(save_file)) saveRDS(doc_term_sparse_matrix,file=paste0(save_file,".rds"))
  
  return(doc_term_sparse_matrix)
  
}


# eg
#topic_model_fit <- LDA(x = doc_term_sparse_matrix, k = 30)
#saveRDS(topic_model_fit,file=paste0(output_dir,"topic_model_fit.rds"))

