### VERBATIM COPY OF maintranformer.R in Louis' directory. Here to keep all files in one place.

##' transformer()
##' 
##' General function to transform data tables into a training matrix.
##' 
##' @param patients a data frame with patient data; see function loadCleanData
##' @param episodes a data frame of healthcare contacts; see function loadCleanData
##' @param list_of_data_tables list of data tables containing all SPARRA information; see function loadCleanData
##' @param variable_list a list of transformed variables to request, together with parameters to call the transformation function with. Format is variable_list=list(name=list(option1=value1,option2=value2...)) where {name} corresponds to a function transformer_{name} which takes arguments {option1},{option2}...
##' @param time_cutoffs a vector of time cutoffs: training matrix uses data before these cutoffs
##' @param hard_max_lookback maximum lookback in days; training matrix is agnostic to any data before time_cutoffs-hard_max_lookback and after time_cutoffs
##' @returns a tibble containing the training matrix  with column time indicating the time cutoff
transformer <- function(patients, episodes, list_of_data_tables,
                        variable_list, time_cutoffs, hard_max_lookback) {
  # Error checking of inputs
  if(!("data.frame" %in% class(patients))) {
    stop("patients must be a data frame")
  }
  ### Table patients may have a column 'cv' affixed to the end; this is not used here  ### Added by James, 19/03/20
  if (ncol(patients)==6) patients=patients[,1:5]
  ### end addition
  if(!all.equal(sapply(patients, class),
                c(id = "integer", gender = "factor", date_of_birth = "Date",
                  date_of_death = "Date", SIMD_DECILE_2016_SCT = "integer"))) {
    stop("column name or type mismatch in patients data frame")
  }
  
  if(!("data.frame" %in% class(episodes))) {
    stop("episodes must be a data frame")
  }
  if(!all.equal(sapply(episodes, class),
                list(id = "integer", time = c("POSIXct", "POSIXt"),
                     time_discharge = c("POSIXct", "POSIXt"),
                     source_table = "factor", source_row = "integer",
                     admission_type = "integer",
                     main_condition = "character",
                     CIS_MARKER = "integer"))) {
    stop("Column name or type mismatch in episodes data frame")
  }
  
  if(!is.list(list_of_data_tables)) {
    stop("list_of_data_tables should be a list")
  }
  if(!all(sapply(lapply(list_of_data_tables, class), function(x) { "data.frame" %in% x }))) {
    stop("the list elements of list_of_data_tables should be data frames")
  }
  if(!(setequal(names(list_of_data_tables),
                c("AE2", "deaths", "PIS", "SMR00", "SMR01", "SMR01E",
                  "SMR04", "SPARRALTC", "SystemWatch")))) {
    stop("list_of_data_tables does not contain the correct list of data frames")
  }
  if(!all(sapply(sapply(list_of_data_tables, names), function(x) { "time" %in% x }))) {
    stop("all the data frames in list_of_data_tables must contain a time variable")
  }
  if(!all(sapply(sapply(list_of_data_tables, names), function(x) { "id" %in% x }))) {
    stop("all the data frames in list_of_data_tables must contain an id variable")
  }
  # targets <- get_targets_regularly(episodes) # GERGO_20190430 - NOTE THAT THIS FUNCTION HAS MOVED AND NOT SOURCED BY DEFAULT ANYMORE!!!
  # if(any(!(time_cutoffs %in% as.integer(unique(targets$time))))) {
  #   stop("some time_cutoffs are not present in targets")
  # }
  # ADD CHECKS of other args
  
  
  
  # Transform all to data tables for performance when filtering time cutoffs
  episodes <- episodes %>% mutate(unixtime = as.integer(time))
  for(i in 1:length(list_of_data_tables)) {
    list_of_data_tables[[i]] <- list_of_data_tables[[i]] %>% mutate(unixtime = as.integer(time))
  }
  hard_max_lookback <- hard_max_lookback * (3600*24) # Convert days to seconds (to use with unix time)
  
  
  # Collect transformed variables separately
  variable_list_names <- names(variable_list)
  var.res <- list()
  all.var.names <- NULL
  for(t in time_cutoffs) {
    # Subset time
    episodes2 <- episodes %>% filter(unixtime < t, unixtime > t-hard_max_lookback) %>% select(-unixtime)
    list_of_data_tables2 <- list_of_data_tables
    for(i in 1:length(list_of_data_tables2)) {
      list_of_data_tables2[[i]] <- list_of_data_tables2[[i]] %>% filter(unixtime < t, unixtime > t-hard_max_lookback) %>% select(-unixtime)
      if(nrow(list_of_data_tables2[[i]]) == 0) {
        warning(glue("the cutoff and lookback choices result in no observations in data table '{names(list_of_data_tables2)[i]}'"))
      }
    }
    patients2 <- patients %>% filter(id %in% episodes2$id)
    
    # Fetch variables
    for(i in 1:length(variable_list)) {
      # Retrieve transformed variable
      res <- do.call(glue("transformer_{variable_list_names[i]}"),
                     c(list(patients = patients2,
                            episodes = episodes2,
                            list_of_data_tables = list_of_data_tables2,
                            time_cutoff = t),
                       variable_list[[i]]))
      gc()  ########## Added by James 
      # Check return value
      if(!is.list(res) ||
         !setequal(names(res), c("matrix", "missing_default")) ||
         !("data.frame" %in% class(res$matrix))) {
        stop(glue("transformer {variable_list_names[i]} is not returning a correctly formed list(matrix=..., missing_default=...) with matrix a data frame"))
      }
      # Append the cutoff time
      res$matrix <- add_column(res$matrix, time_cutoff = t, .before = 1)
      # Store it
      if(is.null(var.res[i][[1]])) {
        # This is the first time we've seen this variable
        var.res[[i]] <- res
        
        # Since first time seeing, we need to check the variable name is not taken
        # and modify if so
        j <- which(names(res$matrix) %in% names(res$missing_default))
        k <- 1
        new.names <- names(res$missing_default)
        while(any(new.names %in% all.var.names)) {
          new.names <- glue("{names(res$missing_default)}.{k}")
          k <- k+1
        }
        names(var.res[[i]]$missing_default) <- new.names
        names(var.res[[i]]$matrix)[j] <- new.names
        
        # Add to list of seen names
        all.var.names <- c(all.var.names, names(var.res[[i]]$missing_default))
      } else {
        # The names might have changed, so first job is to get the names from what
        # we've seen already and convert the newly generated
        new.names <- names(var.res[[i]]$missing_default)
        j <- which(names(res$matrix) %in% names(res$missing_default))
        names(res$missing_default) <- new.names
        names(res$matrix)[j] <- new.names
        
        # Now we proceed to check defaults have not changed and add new rows
        if(!all.equal(res$missing_default,
                      var.res[[i]]$missing_default)) {
          stop(glue("transformer {variable_list_names[i]} is returning inconsistent missing_default object"))
        }
        var.res[[i]]$matrix <-
          bind_rows(var.res[[i]]$matrix, res$matrix)
      }
    }
  }
  
  
  
  # Perform full joins
  res <- var.res[[1]]$matrix
  if(length(var.res) > 1) {
    for(i in 2:length(var.res)) {
      res <- res %>%
        full_join(var.res[[i]]$matrix)
    }
  }
  
  
  
  # Handle possible NAs introduced by joins
  for(i in 1:length(var.res)) {
    for(v in names(var.res[[i]]$missing_default)) {
      res <- res %>%
        mutate_at(v, function(x) { ifelse(is.na(x), var.res[[i]]$missing_default[[v]], x) })
    }
  }
  
  res
}

# targets gives us id, time, target
# my function needs to check that the time_cutoffs given match at least one time in that column
# then join the target by id and time_cutoff
