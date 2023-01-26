library(tidyverse)
library(haven)
library(lubridate)
library(purrr)
#library(varhandle)



# CONDITIONAL MUTATE (of certain rows)

mutate_cond <- function(.data, condition, ..., envir = parent.frame()){
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


# KEEP CLASSES INTACT

# Build an "as.*()" type call automatically using eval
convert_to_same_class_as <- function(x, y_to_convert_to){
  # Only works for converting the class itself, NOT copying factor levels!
  eval(call( paste0("as.", class(y_to_convert_to)), x ))
}

convert_to_same_factor_as <- function(x, y_to_convert_to){
  y_attrs <- attributes(y_to_convert_to)
  y_attrs$class <- NULL
  do.call(factor, args = c(list(x=x), y_attrs))
}

# Or using structure(., class=class(*))
sapply_and_infer_class <- function(x_vec, cur_fun, keep_orig_class=FALSE){
  # This function is necessary! 
  # As soon as the first element in x_vec is NULL, 
  # the output vector's class would otherwise become "NULL",
  # and it will "gain" more and more classes until everything is coercible 
  # (eg NULL will become logical if you add NA, and settles on integer as a compromise with dttm... ) 
  if (keep_orig_class){
    out <- lapply(x_vec, function(y)structure(cur_fun(y), class=class(y)))
    #out <- sapply(x_vec, function(y)structure(cur_fun(y), class=class(y))) # surpirisingly sapply doesnt work here (???)
  } else {
    out <- lapply(x_vec, cur_fun)
    #out <- sapply(x_vec, cur_fun)
  }
  tmp_classes <- unique(sapply(out, class))
  # Find first non "NULL" class
  non_null_class <- tmp_classes[sapply(tmp_classes, function(x)x[1]!="NULL")][[1]]
  out <- as.vector(out)
  class(out) <- non_null_class
  out
}


# IMPUTATION

# This function can be used to impute missing values (to have no effect on prediction!) when predicting via glm models 
impute_zero_or_first_factor_level <- function(df){
  df %>%
    mutate_if(
      is.numeric, 
      function(x) replace(x, is.na(x), convert_to_same_class_as(0, x))
      ) %>% # Replace numbers with 0
    mutate_if(is.factor, function(x) replace(x, is.na(x), convert_to_same_factor_as(levels(x)[1], x)))
}


# ------------------------------------
# TRANSFORMER HELPER FUNCTIONS

#' Collecting values from multiple columns from multiple data tables that are stored in a list
#' @param list_of_data_tables A named list of data_tables (tibbles)
#' @param column_names_flat A vector of column names to join, each in a format <TABLE_NAME.COLUMN_NAME>
#' @param table_key A vector of column names shared by each data_table in the list. 
#'    The individual tables will be grouped_by these columns, and other values will be collected
#' @param time_cutoff filter to records made before this time. Default 2100.
#' @param time_min filter to records made after this time. Default 1900.
#' @return A table with column names = c(table_key, column_names_flat) where each entry 
#'    - for table_key columns is part of the joint unique key that is shared across all tables 
#'    - for column_names_flat columns is a list of the values collected from the respective table for that unique key
collect_and_full_join_from_list_of_tables <- function(
  list_of_data_tables, 
  column_names_flat, 
  table_key=c("id"),
  time_cutoff=dmy_hm("2-5-2100 00:00"),
  time_min=dmy_hm("2-5-1900 00:00")){
  
  column_names_matrix <- str_split(column_names_flat, pattern = '\\.', n=2, simplify = TRUE)
  column_names <- as_tibble(as.data.frame(column_names_matrix)) %>% rename(source_table = V1, column_name = V2)
  
  out = NULL
  for (cur_table_name in unique(column_names$source_table)){
    
    tmp2 <- list_of_data_tables[[cur_table_name]] %>% 
      filter(time< time_cutoff) %>%
      filter(time>= time_min) %>%
      select(table_key, 
        as.character(
          (column_names %>% 
          filter(source_table %in% cur_table_name)
         )$column_name)
        ) %>% 
      group_by_at(vars(one_of(table_key))) %>%
      summarise_all(~list(.)) %>% #funs(list(.))) %>%
      rename_at(vars(-one_of(table_key)), ~ paste0(cur_table_name, ".", .)) #funs(paste0(cur_table_name, ".", .)))
    
    if (is.null(out)) {out <- tmp2
    } else {
      out <- out %>%
        full_join(tmp2, by=table_key)
    }
    
  }
  
  out
  
}





apply_function_elementwise <- function(.data, function_to_apply, table_key=c("id"), keep_orig_class=FALSE){
  .data %>% 
    mutate_at(
      vars(-one_of(table_key)), 
      funs(
        sapply_and_infer_class(
          ., function_to_apply, keep_orig_class = keep_orig_class
        )
      )
    )
}



# Work with UNIX time (faster, it is an integer describing the number of seconds that have passed since 1970 Jan 1)

add_mutate_timeUNIX <- function(.data){
  .data %>% 
    mutate(timeUNIX=as.integer(time))
}

remove_mutate_timeUNIX <- function(.data){
  .data %>%
    select(-c("timeUNIX"))
}

get_days_as_seconds <- function(time_in_days){
  time_in_days*24*60*60
}
