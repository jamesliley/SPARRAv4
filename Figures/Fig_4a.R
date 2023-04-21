# List of libraries to check and install
libs <- c("tidyverse", "patchwork", "ggdist", "dplyr", "tidyr",
          "forcats", "ggrepel", "purrr")

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
# Number of time points
ntx=5

## ROC curves for predictor at various time points
pdf(paste0(plot_dir,"cohort/drift/drift_roc.pdf"),width=3,height=3.5)
roc_2panel(list(xy1,xy2,xy3,xy4,xy5),labels=as.character(test_times),col=gray(((1:ntx)-1)/(1.3*ntx)))
dev.off()

# use ggplot2 to create ROC curves for a predictor at various time points

# The following code uses the default colours in ggplot2.
# It also creates the lower panel with the differences wrt to the curve
# corresponding to the first time point

# ROC curves for each time point
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
roc_plot <- ggplot(df, aes(x = 1 - Specificity, y = Sensitivity, color = Time, linetype = Time)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 8) +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "gray")

# extract spec column from each roc object
spec_list <- lapply(rocs, function(x) x$spec)

# subtract the first roc's spec from each other roc's spec
diffs_list <- lapply(spec_list[-1], function(x) rocs[[1]]$spec - x)

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
ggsave("Figures/pdfs/Fig_4a.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)
