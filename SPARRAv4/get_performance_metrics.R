# Compute various metrics of performance, similar to the ones in the SPARRAv3 report and the data study group report

#' get_performance_metrics()
#' Generates a named list of performance metrics and tables
#' Requires Metrics package installed
#' Input:
#' @param data a tibble that have at least "id" <integer> and "time" <dttm>, "risk_score" <numeric> and "target" <binary> columns
#'     Time is expected to be regularly sampled, and each time that has a non-NA target, may be used as a prediction time. 
#'     Across the multiple prediction times, there will be variations in the performance metrics, so we can get temporal error bars around them.
#'     
#'     Later on we might also want to estimate error at a particular prediction time via k-fold cross-validation on patient ids.
#'     This information may be captured by an addition column in data, "crossval_group". The following scenarios are possible:
#'        - If this is not present or all the same, then we are doing within-sample prediction, no error bars
#'        - If this is present, but there is only one group + NAs, we are doing out of sample prediction, but no error bars
#'        - If this is present and multiple groups, compute the across-patients error and propagate accordingly
#'     
#'     Then both temporal and across-patient variation will contribute to the actual error bars?
#' @param threshold a single value, or vector of values, $\in$ [0, 1] specifying a threshold to split risk scores into a binary prediction
#' Default value of 0.5 
#' 
#' Output:
#' @return A list of named tibbles of performance metrics, each list item reflecting a different threshold:
#' \itemize{
#'   \item{"MMCE"}{Mean misclassification error: the proportion of observations where the prediction, for a given threshold, is incorrect. Range: 0-1; Optimum: 0}
#'   \item{"PPV"}{Positive predictive value: of all observations where TRUE (readmission) is predicted, what proportion were correct. Range: 0-1; Optimum: 1}
#'   \item{"TNR"}{True negative rate: also 'specificity', of all observations who do not readmit (FALSE), what proportion were correctly predicted Range: 0-1; Optimum: 1}
#'   \item{"TPR"}{True positive rate: also 'sensitivity' of all observations who do readmit (TRUE), what proportion were correctly predicted/ Range: 0-1; Optimum: 1}
#'   \item{"AUC"}{Area under the receiver operating characteristic curve (ROC). The ROC shows the trade-off between sensitivity and specificity of a model for differing threshold values. 
#'   The area under the ROC reflects the probability that the model will assign an observation who did readmit with a higher probability of readmitting than one who did not. 
#'   Range: 0-1; Optimum: 1. A value <0.5 indicates worse performance than random chance}
#'   \item{"Brier"}{The Brier score: the mean squared difference between the predicted risk score and the true outcome. Range: 0-1; Optimum: 0}
#'   \item{"F1"}{The F1 score: The harmonic mean of PPV and TPR. Range: 0-1; Optimum: 1}
#'   \item{"logLoss"}{Minus the average log difference between the predicted risk and the true outcome. Range: 0-1; Optimum: 0}
#'}
#' 
#' The function also saves this table to the disk within the intermediate data folder
get_performance_metrics <- function(data, threshold = 0.5){
  # Return a list of performance metrics
  metrics <- data %>%
    group_by(crossval_group) %>%
    summarise(
      F1 = Metrics::f1(target, risk_score),
      AUC = fast_AUC(risk_score, target),
      Brier = mean((risk_score - target) ^ 2),
      logLoss = Metrics::logLoss(target, risk_score)
    ) %>%
    gather(measure, value, F1:logLoss) %>%
    group_by(measure) %>%
    summarise(min = min(value),
              mean = mean(value),
              max = max(value))
  metrics_list <- lapply(threshold, function(x) {
    data %>%
    mutate(prediction = as.numeric(risk_score >= x)) %>%
    group_by(crossval_group) %>%
    summarise(
      MMCE = mean(prediction != target),
      TPR = sum(prediction * target) / sum(target),
      TNR = sum((1 - prediction) * (1 - target)) / sum(1 - target),
      PPV = sum(prediction * target) / sum(prediction)
      ) %>%
    gather(measure, value, MMCE:PPV) %>%
    group_by(measure) %>%
    summarise(min = min(value),
              mean = mean(value),
              max = max(value)) %>%
      bind_rows(metrics)
  }
  )
  names(metrics_list) <- sprintf("threshold=%.2f", threshold)
  return(metrics_list)
}

#' get_calibration()
#' Generates a binned calibration table
#' @param data a tibble that have at least "id" <integer> and "time" <dttm>, "risk_score" <numeric> in [0, 1] and "target" <binary> columns
#'     Time is expected to be regularly sampled, and each time that has a non-NA target, may be used as a prediction time. 
#'     Across the multiple prediction times, there will be variations in the performance metrics, so we can get temporal error bars around them.
#'     
#'     Later on we might also want to estimate error at a particular prediction time via k-fold cross-validation on patient ids.
#'     This information may be captured by an addition column in data, "crossval_group". The following scenarios are possible:
#'        - If this is not present or all the same, then we are doing within-sample prediction, no error bars
#'        - If this is present, but there is only one group + NAs, we are doing out of sample prediction, but no error bars
#'        - If this is present and multiple groups, compute the across-patients error and propagate accordingly
#'     
#'     Then both temporal and across-patient variation will contribute to the actual error bars?
#' @param nbins an <integer> specifying the number of bins to split `risk_score`
#' @return A tibble with mean admitted (observed and predicted) for each of `crossval_group` 
#' @export 
get_calibration <- function(data, nbins = 10) {
  out <- data %>%
    mutate(bin = cut(risk_score, # Treating cut as factor enables "complete" for empty bins
                                  seq(0, 1, length.out = nbins + 1),
                                  include.lowest = TRUE)
           ) %>%
    # mutate(bin = as.character(cut(risk_score,
    #                               seq(0, 1, length.out = nbins + 1),
    #                               include.lowest = TRUE)),
    #        bin = gsub("[^0-9,\\.]", "", bin)) %>%
    group_by(crossval_group, bin) %>%
    summarise(actual = sum(target) / n(), n = n()) %>%
    group_by(bin) %>%
    complete(bin) %>%
    #filter(!is.na(bin)) %>%
    summarise(min      = min(actual * 100),
              mean     = mean(actual * 100),
              max      = max(actual * 100),
              min_n    = min(n),
              se       = sqrt(mean / 100 * (1 - mean / 100) / mean(n)) * 100, 
              lower_95 = max(mean - 1.96 * se, 0),
              upper_95 = min(mean + 1.96 * se, 100)) %>%
    mutate(expected = (seq(0, 100 * (1 - 1/ nbins), length.out = nbins) +
                         seq(100 / nbins, 100, length.out = nbins)) / 2)
  return(out)
}


#' fast_auc()
#' Based on github/traversc/fastAUC.R
#' Generates the area under the receiver operating characteristic curve
#' @param probs A vector of predicted probabilities
#' @param class The true probabilities
#' @return The AUC, a numeric value
fast_AUC <- function(probs, class) {
  probs[probs > 1] <- 1
  probs[probs< 0] <- 0
  x1 <- probs[class]; n1 <- length(x1)
  x2 <- probs[!class]; n2 <- length(x2)
  r <- frank(c(x1, x2))
  auc <- (sum(r[1:n1]) - n1 * (n1 + 1) / 2) / n1 / n2
  return(auc)
}

#' get_auc()
#' Generates the mean and variance of the AUC according to the Cortes paper
#' I don't believe this is relevant, except in very small cases
#' and even then...
#' @param k A fixed number of misclassifications
#' @param m The number of positive cases
#' @param n The number of negative cases
#' @return A list of the  
#' @export 
get_auc <- function(k, m, n) {
  mean_auc <- 1 - 
      (k / (m + n)) - 
      ((n - m) ^ 2 * (m + n + 1) / 
         (4 * m * n)) *
      (k / (m + n) - 
      sum(choose(m + n, 0:(k - 1))) / 
      sum(choose(m + n + 1, 0:k)))
  T  <- 3 * ((m - n) ^ 2 + m + n) + 2
  Q0 <- (m + n + 1) * T * k ^ 2 +
    ((- 3 * n ^ 2 + 3 * m * n + 3 * m + 1) * T
     - 12 * (3 * m * n + m + n) - 8) * k + 
    (- 3 * m ^ 2 + 7 * m + 10 * n + 3 * n * m + 10) * T - 
    4 * (3 * m * n + m + n + 1)
  
  Q1 <- T * k ^ 3 + 3 * (m - 1) * T * k ^ 2 + 
    ((- 3 * n ^ 2 + 3 * m * n - 3 * m + 8) * T - 
    6 * (6 * m * n + m + n)) * k + (- 3 * m ^ 2 + 7 * (m + n) +3 * m * n) * T - 
  2 * (6 * m * n + m + n)
  
  sigma_auc <- ((m + n + 1) * (m + n) * (m + n - 1) * 
                  T * ((m + n - 2) * Zi(k, m, n, 4) - (2 * m - n + 3 * k - 10) * Zi(k, m, n, 3))) / 
                  (72 * m ^ 2 * n ^ 2) + 
    ((m + n + 1) * (m + n) * T * (m ^ 2 - n * m + 3 * k * m - 5 * m + 2 * k ^ 2 - n * k + 12 - 9 * k) * Zi(k, m, n, 2)) / 
    (48 * m ^ 2 * n ^ 2) -
    ((m + n + 1) ^ 2 * (m - n) ^ 4 * Zi(k, m, n, 1) ^ 2) / 
    (16 * m ^ 2 * n ^ 2) - 
    ((m + n + 1) * Q1 * Zi(k, m, n, 1)) / (72 * m ^ 2 * n ^ 2) +
    k * Q0 / (144 * m ^ 2 * n ^ 2)
  return(list(mean_auc, sigma_auc))
  
}

Zi <- function(k, m, n, i) {
  sum(choose(m + n + 1 - i, 0:(k - i))) / 
    sum(choose(m + n + 1, 0:(k)))
    
}
get_auc(500, 5000, 2000)



#' save_performance_metrics()
#' Runs all performance metrics and saves to file
#' Input:
#' @param pred a tibble that have at least "ID" <integer> "risk_score" <numeric> and "target" <binary> and "crossval_group" columns
#' @param dir The root directory to store output
#' @param suffix A string appended to the end of the saved output to identify
#' 
save_performance_metrics <- function(pred, dir, suffix) {
  pred <- drop_na(pred)
  
  performance_metrics <- get_performance_metrics(pred, threshold = c(0.5, 0.4, 0.3))
  saveRDS(performance_metrics,
          paste(dir, "metrics_", suffix, ".RDS", sep = ""))
  
  plot_metrics(performance_metrics, what = "TPR")
  ggsave(paste(dir, "TPR_", suffix, ".png", sep = ""), width = 10, height = 6, dpi = 300)
  
  # Getting a confidence interval on AUC via DeLong's method
  roc <- roc(pred$target, pred$risk_score, ci = TRUE, of = "auc")
  saveRDS(list(`95% CI` = c(roc$ci[1], roc$ci[3]) , AUC = roc$ci[2]),
          paste(dir, "ROC_", suffix, ".RDS", sep = ""))
  
  pred %>%
    plot_roc
  ggsave(paste(dir, "ROC_", suffix, ".png", sep = ""),
         width = 7, height = 7, dpi = 300)
  
  
  calibration_table <- pred %>%
    get_calibration() 
  plot_calibration_se(calibration_table)
  ggsave(paste(dir, "calibration_se_", suffix, ".png", sep = ""),
         width = 7, height = 7, dpi = 300)
  
  plot_calibration_crossval(calibration_table)
  ggsave(paste(dir, "calibration_cv_", suffix, ".png", sep = ""),
         width = 7, height = 7, dpi = 300)
}
