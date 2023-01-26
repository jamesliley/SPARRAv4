req_packages = list(
  # h2o
  "h2o",
  
  # Package building and documentation
  "devtools",
  "roxygen2",
  
  # Tidyverse recommended
  "tidyverse",
  "styler", # # Style fixing
  "lintr", # Style checing
  "roxygen2", # Documentation
  "testthat", # Testing
  "stringr",
  "purrr",
  "forcats",
  
  # Data table manipulation
  "haven",
  "data.table",
  "dplyr",
  "reshape2",
  "tidyr",
  "varhandle",
  "ggalluvial",
  "labelled", # Deals with labelled data formats
  "fuzzyjoin",  # Enables partial joins (used for aggregated code descriptions)
  
  # Parallel processing
  "doParallel",
  "plyr",
  "future",
  "furrr",
  "compiler",
  "fst", # Parallel read-write
  "import", # (caret needs this)
  "foreach", # (caret needs this)
  
  # Time manipulation
  "lubridate",
  
  # Plotting and statistics
  "plotly",
  "ggplot2",
  "corrplot",
  "matrixStats",
  "Matrix",
  "AUC",
  "ROCR",
  "pROC",
  "PRROC",
  "cvAUC",
  "corrplot",
  "fBasics",
  "car",
  "latex2exp",
  
  # Machine learning
  "caret", "e1071",
  "mlr",
  "kernlab",
  "randomForest",
  "topicmodels",
  "maptpx", # alternative topic models
  "nlme", # AIC.loglik()
  "xgboost",
  "Metrics",
  "groupdata2",
  "speedglm",
  "glmnet",
  "ranger", # Fast random forest implementation
  "PReMiuM", # installed by Sebastian or Louis?
  
  
  # Example data
  "nycflights13",
  
  # Benchmarking
  "benchmarkme",
  
  # Reading from excel
  "readxl",
  
  # Exporting to excel
  # "xlsx", # doesn't work, never use
  
  # Categorical association
  "lsr",
  
  # Profiling
  "profvis",
  "tictoc",
  "pryr", #Memory profiling
  
  #Reporting
  "gridExtra", 
  "rmarkdown", # Not found
  "knitr" # Not found
  
)



##' Install packages
##' 
##' @name install_required_packages()
##' @returns list of packages which failed to install
install_required_packages=function() {
fail=c()

# Function to test if package is installed
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])


for (mypkg in req_packages){
  if (is.installed(mypkg) == FALSE){
    install.packages(mypkg)
    if (is.installed(mypkg)){
      cat(sprintf("Successfully installed package \"%s\" \n", mypkg))  
    } else {
      cat(sprintf("\n FAILED TO INSTALL package \"%s\" \n \n", mypkg))
      fail=c(fail,mypkg)
    }
    
  } else {
    cat(sprintf("Package \"%s\" is already installed \n", mypkg))
  }
}

return(fail)
}

##' Load required packages, as many as possible. On NSH, windows machine cannot load h2o, linux cannot load topicmodels.
##' 
##' @name load_packages()
##' @returns list of packages which failed to load
load_packages=function() {
  ip=installed.packages()
  to_load=intersect(ip[,1],unlist(req_packages))
  for (i in 1:length(to_load)) library(to_load[i],character.only=TRUE)
}


# Install newer version of "haven", as we need it to work with .zsav files
# Doesn't work
# eDRISInstall_custom("haven", force_R_version = "3.5", force_reinstall = TRUE)
