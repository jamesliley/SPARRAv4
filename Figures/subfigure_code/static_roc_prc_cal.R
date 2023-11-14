# List of libraries to check and install
libs <- c("tidyverse", "patchwork", "ggdist", "dplyr", "tidyr",
          "forcats", "ggrepel", "purrr","lubridate")

# Create an empty list to store the missing libraries
missing_libs <- list()

# Loop through each library and check if it's installed
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    missing_libs[[lib]] <- TRUE
  }
}

# Install and load the missing libraries
if (length(missing_libs) > 0) {
  install.packages(names(missing_libs))
}

# Load all libraries
lapply(libs, require, character.only = TRUE)

source("Figures/util.R") # for import_sparra_expr()


## ROC ####
plot_dir = "Analysis/time_attenuation/"
eval(import_sparra_expr(paste0(plot_dir, "cohort/drift/drift_roc.txt")))

## old code
test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)
# Number of time cutoffs
ntx=5

## ROC curves for predictor at various time cutoffs
pdf(paste0(plot_dir,"cohort/drift/drift_roc.pdf"),width=3,height=3.5)
roc_2panel(list(xy1,xy2,xy3,xy4,xy5),labels=as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

# use ggplot2 to create ROC curves for a predictor at various time cutoffs

# The following code uses the default colours in ggplot2.
# It also creates the lower panel with the differences wrt to the curve
# corresponding to the first time cutoff

# ROC curves for each time cutoff
rocs <- list(xy1, xy2, xy3, xy4, xy5)

# Generate colors
colors <- scales::hue_pal()(length(rocs))

# Combine the ROC curves into a single data frame
df <- data.frame(
  Time = rep(as.character(test_times), each = length(xy1$sens)),
  Sensitivity = c(xy1$sens, xy2$sens, xy3$sens, xy4$sens, xy5$sens),
  Specificity = c(1 - xy1$spec, 1 - xy2$spec, 1 - xy3$spec, 1 - xy4$spec,
                  1 - xy5$spec)
)

# Convert test_times to character
test_times <- as.character(test_times)

# Create the ROC plot
roc_plot <- ggplot(df, aes(x = Specificity, y = Sensitivity, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")

# extract spec column from each roc object
spec_list <- lapply(rocs, function(x) x$spec)

# subtract the first roc's spec from each other roc's spec
diffs_list <- lapply(spec_list[-1], function(x) x - rocs[[1]]$spec)

# combine the diffs and add a Time column
Time <- rep.int(rep(test_times[2:5], each = 100), length(diffs_list))
diffs <- data.frame(
              Delta = unlist(diffs_list),
              Specificity = rep(rocs[[1]]$spec, length(diffs_list)),
              Time = Time)

# Create the difference plot
diff_plot <- ggplot(diffs, aes(x = 1 - Specificity, y = Delta, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors[2:5]) +
  labs(x = "1 - Specificity", y = expression(paste(Delta, " Sensitivity"))) +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "gray")

# Combine the two plots vertically - alternative heights = c(3, 1)
print(roc_plot / diff_plot) + plot_layout(heights = c(2/3, 1/3))

# Save the plot to a file
ggsave("Figures/pdfs/Unsorted/static_score_roc.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)



### PRC ####

# old code
plot_dir = "Analysis/time_attenuation/"
eval(import_sparra_expr(paste0(plot_dir, "cohort/drift/drift_prc.txt")))

test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)
# Number of time cutoffs
ntx=5
pdf(paste0(plot_dir,"cohort/drift/drift_prc.pdf"),width=3,height=3.5)
prc_2panel(list(px1,px2,px3,px4,px5),labels=as.character(test_times),
           col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

## new code
# PRC curves for each time cutoff
prc <- list(px1, px2, px3, px4, px5)

# Generate colors
colors <- scales::hue_pal()(length(prc))

# Combine the PRC curves into a single data frame
df <- data.frame(
  Time = rep(as.character(test_times), each = length(px1$sens)),
  sens = c(px1$sens, px2$sens, px3$sens, px4$sens, px5$sens),
  ppv = c(px1$ppv, px2$ppv, px3$ppv, px4$ppv, px5$ppv)
)

# Convert test_times to character
test_times <- as.character(test_times)

# Create the PRC plot
prc_plot <- ggplot(df, aes(x = sens, y = ppv, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(x = "Recall", y = "Precision") +
  theme_minimal(base_size = 8)

# extract sens column (precision) from each PRC object
sens_list <- lapply(prc, function(x) x$sens)

# subtract the first prc's sens from each other prc's sens
diffs_list <- lapply(sens_list[-1], function(x) x - prc[[1]]$sens) #prc[[1]]$sens - x

# combine the diffs and add a Time column
Time <- rep.int(rep(test_times[2:5], each = 100), length(diffs_list))
diffs <- data.frame(
  Delta = unlist(diffs_list),
  Precision = rep(prc[[1]]$sens, length(diffs_list)),
  Time = Time)

# Create the difference plot
diff_plot <- ggplot(diffs, aes(x = Precision, y = Delta, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors[2:5]) +
  labs(x = "Recall", y = expression(paste(Delta, " Precision"))) +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "gray")

# Combine the two plots vertically
print(prc_plot / diff_plot) + plot_layout(heights = c(2/3, 1/3))

# Save the plot to a file
ggsave("Figures/pdfs/Unsorted/static_score_prc.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)




### CAL ####

# old code
plot_dir = "Analysis/time_attenuation/"
eval(import_sparra_expr(paste0(plot_dir, "cohort/drift/drift_cal.txt")))

test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)

## code that generates old plot
# pdf(paste0(plot_dir,"recalculate/CAL/cal.",xname,".pdf"),width=3,height=3.5)
# cci=rgb((ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),(ntx-1)/(1.3*ntx),alpha=0.5) # colour for confidence envelope
# cal_2panel(list(cx1,cx2,cx3,cx4,cx5),as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)),
#            ci_col=c(NA,NA,NA,NA,cci))
# dev.off()

## code to generate the new plot
ntx <- length(test_times)

# Calibration curves for each time cutoff
cals <- list(cx1, cx2, cx3, cx4, cx5)

# Generate colors
colors <- scales::hue_pal()(length(cals))

# Combine the PRC curves into a single data frame
df <- data.frame(
  Time = rep(as.character(test_times), each = length(cx1$x)),
  x = c(cx1$x, cx2$x, cx3$x, cx4$x, cx5$x),
  y = c(cx1$y, cx2$y, cx3$y, cx4$y, cx5$y)
)

# Convert test_times to character
test_times <- as.character(test_times)

# Create the calibration plot
cal_plot <- ggplot(df, aes(x = x, y = y, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(x = "Predicted", y = "Observed") +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")

# Differences in lower panel
# Define the reference time cutoff
ref_time <- 1

# Calculate the differences between each calibration curve and the reference curve
cal_diffs <- lapply(cals[-ref_time], function(x) {
  data.frame(
    x = x$x,
    y = x$y - cals[[ref_time]]$y[x$x == cals[[ref_time]]$x]
  )
})

# Combine the differences into a single data frame
df <- data.frame(
  Time = rep(as.character(test_times[-ref_time]), each = length(cals[[ref_time]]$x)),
  x = rep(cals[[ref_time]]$x, length(test_times[-ref_time])),
  y = do.call("c", lapply(cal_diffs, `[[`, "y"))
)

# Create the plot
cal_diff_plot <- ggplot(df, aes(x = x, y = y, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  labs(x = "Predicted", y = expression(paste(Delta," Calibration"))) +
  scale_color_manual(values = colors[2:5]) +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "gray")


# Combine the two plots vertically
print(cal_plot / cal_diff_plot) + plot_layout(heights = c(2/3, 1/3))

# Save the plot to a file
ggsave("Figures/pdfs/Unsorted/static_score_cal.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)

