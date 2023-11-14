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
  labs(x = "Predicted", y = "Difference in Calibration") +
  scale_color_manual(values = colors[2:5]) +
  theme_minimal(base_size = 8) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "gray")


# Combine the two plots vertically
print(cal_plot / cal_diff_plot) + plot_layout(heights = c(2/3, 1/3))

# Save the plot to a file
ggsave("Figures/pdfs/Fig_4c.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)
