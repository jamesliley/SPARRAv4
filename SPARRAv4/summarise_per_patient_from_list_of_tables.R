#' summarise_per_patient_from_list_of_tables
#' Collecting values from multiple columns from multiple data tables that are stored in a list
#' @param list_of_data_tables A named list of data_tables (tibbles)
#' @param column_names_flat A vector of column names to join, each in a format <TABLE_NAME.COLUMN_NAME>
#' @param summary_func A function to apply to the collected observations per patient per column
#' @param table_key A vector of column names shared by each data_table in the list. 
#'    The individual tables will be grouped_by these columns, and other values will be collected
#' @return A table with column names = c(table_key, column_names_flat) where each entry 
#'    - for table_key columns is part of the joint unique key that is shared across all tables 
#'    - for column_names_flat columns is a list of the values collected from the respective table for that unique key
summarise_per_patient_from_list_of_tables <- function(
  list_of_data_tables, 
  column_names_flat,
  summary_func = NULL, # Importantly passing summary_func as a "quosure" [using quo() in the function definition] enables hybrid evaluation
  table_key=c("id")){
  
  # Check if summary_func is provided
  if (is.null(summary_func)) stop("You need to provide a summary function as input to summarise_per_patient_from_list_of_tables")
  
  column_names_matrix <- str_split(column_names_flat, pattern = '\\.', n=2, simplify = TRUE)
  column_names <- as_tibble(as.data.frame(column_names_matrix)) %>% rename(source_table = V1, column_name = V2)
  
  if ("quosure" %in% class(summary_func)) {
    summary_func_general=eval(parse(text=paste0("~",quo_name(summary_func))))
  } else summary_func_general=summary_func
  out = NULL
  for (cur_table_name in unique(column_names$source_table)){
    
    tmp2 <- list_of_data_tables[[cur_table_name]] %>%
      select(table_key,
             as.character(
               (column_names %>%
                  filter(source_table %in% cur_table_name)
               )$column_name)
      ) %>%
      group_by_at(vars(one_of(table_key))) %>%
      summarise_all(summary_func_general) %>%
      rename_at(vars(-one_of(table_key)), ~paste0(cur_table_name,".",.))
    #funs(paste0(cur_table_name, ".", .)))
    
    # # Parallel version
    # tmp2 <- list_of_data_tables[[cur_table_name]] %>%
    #   select(table_key,
    #          as.character(
    #            (column_names %>%
    #               filter(source_table %in% cur_table_name)
    #            )$column_name)
    #   )# %>%
    #   
    #   tmp2 <- tmp1 %>% arrange(id) %>% mutate(TMP_MULTIPROCESS_GROUP = ntile(id, detectCores())) 
    #   tmp22 <- tmp2 %>% group_by(id) %>% mutate(TMP_MULTIPROCESS_GROUP=min(TMP_MULTIPROCESS_GROUP))
    #   tmp3 <- tmp2 %>% select_if(is.numeric) %>% group_by(TMP_MULTIPROCESS_GROUP) %>% nest()
    #   
    #   tic()
    #   tmp4 <- tmp3 %>% 
    #     transmute( summs = map(data, function(x)x %>% group_by(id) %>% summarise_all(sum))) %>% 
    #     unnest(summs)
    #   toc()
    #   
    #   library(furrr)
    #   options(future.globals.maxSize = 20000*(1024^2))
    #   plan(multiprocess)
    #   tic()
    #   tmp4 <- tmp3 %>% 
    #     transmute( summs = future_map(data, function(x)x %>% group_by(id) %>% summarise_all(sum))) %>% 
    #     unnest(summs)
    #   toc()
    #   
    #   group_by_at(vars(one_of(table_key))) %>%
    #   summarise_all(summary_func) %>%
    #   rename_at(vars(-one_of(table_key)), funs(paste0(cur_table_name, ".", .)))
    
    if (is.null(out)) {out <- tmp2
    } else {
      out <- out %>%
        full_join(tmp2, by=table_key)
    }
    
  }
  
  out
  
}