# Load the required libraries
library("tidyverse")
library("patchwork")
library("egg")


## Stratification by age ####
#setwd("~/SPARRAv4") # CAV: If you are using a project, this is not needed
source("Figures/util.R") # for import_sparra_expr()
#source("SPARRAv4/auxiliary.R") # for roc_2panel()

# conversion of figure 3a performance_by_age
plot_dir="Analysis/full_model/"
eval(import_sparra_expr(paste0(plot_dir, "Analytics/performance_by_age.txt")))

# Define the age group labels
age_labels <- c("(0, 20]", "(20, 30]", "(30, 40]", "(40, 50]",
                "(50, 60]", "(60, 70]", "(70, 80]", "(80, 120]")

# Create a data frame for plotting
plot_data <- data.frame(
  age = factor(rep(seq_len(ncol(perf)), times = 2)),  # x-axis values for positioning each version
  version = rep(c("V4", "V3"), each = ncol(perf)),  # Specify the versions
  auroc = c(perf[1, ], perf[2, ]),  # Repeat the AUROC values for each version
  sd = c(perf[3,], perf[4,]) # AUROC sd
)

# Plot per-group AUCs for both models
p3a1 <-
    ggplot(plot_data, aes(x = age, y = auroc, group = version)) +
    geom_point(shape = 18, size = 2, aes(col = version)) +  # Use versions as points
    geom_errorbar(aes(ymin = auroc - 3*sd, ymax = auroc + 3*sd, col = version),
                  width = 0.2, linewidth = 0.4) + # position = position_dodge(width = 1)
    xlab("Age group") +
    ylab("AUROC") +
    scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
    scale_x_discrete(labels = age_labels) +
    #scale_x_continuous(minor_breaks = NULL, breaks = 1:(length(age_labels)), labels = age_labels) +  # Set custom x-axis labels +  # Set custom x-axis labels
    # Set y-axis to start from 0
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.68, 0.81)) +
    guides(col = guide_legend(title = NULL)) +
    coord_fixed(ratio = 60) + # Set aspect ratio to 1:1
    theme_bw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8),
      #axis.text.y = element_text(angle = 90, hjust = 1, size = 7),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank() # Remove vertical grid lines
    )

# Sub-panels

plot_data2 <- data.frame(
  age = factor(seq_len(ncol(perf))),
  difference = plot_data$auroc[plot_data$version == "V4"] -
    plot_data$auroc[plot_data$version == "V3"],
  percentage_increase = ((plot_data$auroc[plot_data$version == "V4"] / plot_data$auroc[plot_data$version == "V3"]) - 1) * 100,
  prop_events = xrate
)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3a2 <-
  ggplot(plot_data2, aes(x = age, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Age group") +
  ylab("% increase") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 0.5) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 4.4)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3a3 <-
  ggplot(plot_data2, aes(x = age, y = prop_events)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Age group") +
  ylab("% events") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 7) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 0.33)

p3a <- ggarrange( p3a1 + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank(),
                       plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.2, unit = "cm")),
           p3a2 + theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
           p3a3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
           nrow = 3)

ggsave("Figures/pdfs/Unsorted/strat_by_age.pdf", p3a,
       width = 10, height = 20, units = "cm",
       device = cairo_pdf)






## Stratification by SIMD ####

#Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
eval(import_sparra_expr("Analysis/full_model/Analytics/performance_by_simd.txt"))

# old code
plot_dir="Analysis/full_model/"

#### new version of the figure
library(ggplot2)

labs <- 1:10

xsc <- 25
swidth <- 0.1
m0 <- dim(perf)[1] / 2
n0 <- dim(perf)[2]
perfmin <- min(perf[1:m0, ]) - 0.02

# Create a data frame for plotting
plot_data <- data.frame(
  x = factor(rep((m0 + 2) * (1:n0) - floor(m0 / 2) - 1, each = m0)),
  y = c(perf[1:m0, ], perf[1:m0, ]),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)

sd_values <- c(0.0006950962, 0.0007206066, 0.0007468640, 0.0007848581, 0.0007990873, 0.0008351305, 0.0008600128, 0.0008998324, 0.0009158358, 0.0009865907,
               0.0007128243, 0.0007392061, 0.0007657252, 0.0008038775, 0.0008194903, 0.0008559886, 0.0008828196, 0.0009234900, 0.0009403307, 0.0010172065)
plot_data$sd <- sd_values
# Set up the plot
p3b1 <-
  ggplot(plot_data, aes(x = x, y = y, group = version)) +
  geom_point(shape = 18, size = 2, aes(col = version)) +
  geom_errorbar(aes(ymin = y - 3*sd, ymax = y + 3*sd, col = version), width = 0.2, size = 0.4) + # position = position_dodge(width = 1)
  scale_x_discrete(labels = labs) + # breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1,
  #scale_y_continuous(expand = c(0, 0), limits = c(min(plot_data$y)-0.1, 0.9)) +
  scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.68, 0.81)) +
  labs(x = "SIMD", y = "AUROC") +
  guides(col = guide_legend(title = NULL)) +
  coord_fixed(ratio = 60) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

plot_data2 <- data.frame(
  simd = factor(seq_len(ncol(perf))),
  difference = plot_data$y[plot_data$version == "V4"] -
    plot_data$y[plot_data$version == "V3"],
  percentage_increase = ((plot_data$y[plot_data$version == "V4"] / plot_data$y[plot_data$version == "V3"]) - 1) * 100,
  prop_events = xrate
)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3b2 <-
  ggplot(plot_data2, aes(x = simd, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("SIMD") +
  ylab("% increase") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = labs) +
  coord_fixed(ratio = 0.5) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 4.4)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3b3 <-
  ggplot(plot_data2, aes(x = simd, y = prop_events)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("SIMD") +
  ylab("% events") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = labs) +
  coord_fixed(ratio = 7) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 0.33)

p3b <- ggarrange( p3b1 + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank(),
                               plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                  p3b2 + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank(),
                               plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                  p3b3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
                  nrow = 3)

ggsave("Figures/pdfs/Unsorted/strat_by_simd.pdf", p3b,
       width = 35, height = 30, units = "cm",
       device = cairo_pdf)




### Stratification by v3 cohort ####

# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Analytics/performance_by_v3_cohort.txt")))

### new version
labs <- c("FEC", "LTC", "YED", "U16")

xcol <- c("black", "red")
xsc <- 20
swidth <- 0.1
m0 <- dim(perf)[1] / 2
n0 <- dim(perf)[2]
perfmin <- min(perf[1:m0, ]) - 0.02

# Create a data frame for plotting
plot_data <- data.frame(
  x = factor(rep((m0 + 2) * (1:n0) - floor(m0 / 2) - 1, each = m0)),
  y = c(perf[1:m0, ], perf[1:m0, ]),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)
se <- c(0.0005187255, 0.0003386057, 0.0004669980, 0.0008822474,
        0.0005326287, 0.0003475607, 0.0004784836, 0.0008933571)
plot_data$se <- se
# Set up the plot
p3c1 <-
  ggplot(plot_data, aes(x = x, y = y, group = version)) +
  geom_point(aes(color = version), shape = 18, size = 2) +
  geom_errorbar(aes(ymin = y - 3*se, ymax = y + 3*se, col = version), width = 0.2, size = 0.4) + # position = position_dodge(width = 1)
  
  scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(
    #breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1,
    labels = labs,
    #expand = c(0, 0), limits = c(0, 15),  # Adjust the expand argument to extend the x-axis limits
    position = "bottom"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.68, 0.81)) +
  labs(x = "Cohort", y = "AUROC") +
  guides(col = guide_legend(title = NULL)) +
  coord_fixed(ratio = 60) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )



plot_data2 <- data.frame(
  x = factor(seq_len(ncol(perf))),
  difference = plot_data$y[plot_data$version == "V4"] -
    plot_data$y[plot_data$version == "V3"],
  percentage_increase = ((plot_data$y[plot_data$version == "V4"] / plot_data$y[plot_data$version == "V3"]) - 1) * 100,
  prop_events = xrate
)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3c2 <-
  ggplot(plot_data2, aes(x = x, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Cohort") +
  ylab("% increase") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = labs) +
  coord_fixed(ratio = 0.5) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 4.4)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3c3 <-
  ggplot(plot_data2, aes(x = x, y = prop_events)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Cohort") +
  ylab("% events") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = labs) +
  coord_fixed(ratio = 7) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 0.33)

p3c <- ggarrange( p3c1 + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank(),
                               plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                  p3c2 + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank(),
                               plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                  p3c3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
                  nrow = 3)

ggsave("Figures/pdfs/Fig_3c.pdf", p3c,
       width = 10, height = 20, units = "cm",
       device = cairo_pdf)

# Set the height ratios
heights <- c(2, 1, 1)  # Specify the heights of each row

p3 <- ggarrange( p3a1 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3b1 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3c1 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3a2 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3b2 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3c2 + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
                 p3a3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
                 p3b3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
                 p3c3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
                 nrow = 3, heights=heights)


ggsave("Figures/pdfs/Unsorted/strat_by_v3.pdf", p3,
       width = 35, height = 30, units = "cm",
       device = cairo_pdf)
