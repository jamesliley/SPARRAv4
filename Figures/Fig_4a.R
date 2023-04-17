# Load the required libraries
library("tidyverse")
library("patchwork")

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

library(ggplot2)

# Combine ROC data into a single data frame
df <- data.frame(
  Model = rep(as.character(test_times), each = length(xy1$sens)),
  Sensitivity = c(xy1$sens, xy2$sens, xy3$sens, xy4$sens, xy5$sens),
  Specificity = 1 - c(xy1$spec, xy2$spec, xy3$spec, xy4$spec, xy5$spec)
)

# Convert test_times to character
test_times <- as.character(test_times)


# Define custom colors
custom_colors <- c("red", "blue", "green", "purple", "orange")

# Create ggplot object
p <- ggplot(df, aes(x = 1 - Specificity, y = Sensitivity, color = Model, linetype = Model)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = custom_colors) +  # Set custom colors
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()

# Print the ggplot object
print(p)

# Save the plot to a file
ggsave("Figures/pdfs/Fig_4a.pdf",
        width = 7.5, height = 7.25, units = "cm",
        device = cairo_pdf)



