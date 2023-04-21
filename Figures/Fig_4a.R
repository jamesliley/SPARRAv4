# List of libraries to check and install
libs <- c("tidyverse", "patchwork", "viridis", "ggdist")

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

# code converted to use ggplot2 for creating ROC curves for a predictor at various time points:

# Combine ROC data into a single data frame
df <- data.frame(
  Model = rep(as.character(test_times), each = length(xy1$sens)),
  Sensitivity = c(xy1$sens, xy2$sens, xy3$sens, xy4$sens, xy5$sens),
  Specificity = 1 - c(xy1$spec, xy2$spec, xy3$spec, xy4$spec, xy5$spec)
)

# Convert test_times to character
test_times <- as.character(test_times)

# Define custom colours
# viridis colour palette is known to be colour-blind friendly
# Define a custom gradient color scale centered on blue
colors <- viridis(5, begin = 0.3, end = 0.7, direction = -1)

# Create ggplot object
p <- ggplot(df, aes(x = 1 - Specificity, y = Sensitivity, color = Model, linetype = Model)) +
  geom_line(linetype = "solid") +
  scale_color_manual(values = colors) + # manually specify colours
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal() +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "gray") # Add dashed gray line

# Print the ggplot object
print(p)

# Save the plot to a file
ggsave("Figures/pdfs/Fig_4a.pdf",
        width = 7.5, height = 7.25, units = "cm",
        device = cairo_pdf)



