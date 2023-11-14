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
ggsave("Figures/pdfs/Fig_4b.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)
