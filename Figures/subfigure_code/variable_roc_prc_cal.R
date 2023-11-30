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



## old code
test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)
#timecol=palette.colors(palette = "R4")[1:length(test_times)]
timecol=as.vector(palette.colors(palette = "Okabe-Ito")[c(1,2,4,5,8)])

# Number of time cutoffs
ntx=5




## ROC ####
plot_dir = "Analysis/time_attenuation/"
## ROC ####
eval(import_sparra_expr("Analysis/time_attenuation/recalculate/ROC/roc.super.txt"))

## ROC curves for predictor at various time cutoffs
rx=roc_2panel_gg(
  list(xy1,xy2,xy3,xy4,xy5),
  labels=as.character(test_times),
  col=timecol,
  legend_title="Time")

# Save the plot to a file
sc=1
ggsave("Figures/pdfs/Unsorted/variable_score_roc.pdf",rx,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)



### PRC ####

plot_dir = "Analysis/time_attenuation/"
eval(import_sparra_expr("Analysis/time_attenuation/recalculate/PRC/prc.super.txt"))

px=prc_2panel_gg(
  list(px1,px2,px3,px4,px5),
  labels=as.character(test_times),
  col=timecol)

ggsave("Figures/pdfs/Unsorted/variable_score_prc.pdf",px,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)




### CAL ####

plot_dir = "Analysis/time_attenuation/"
eval(import_sparra_expr("Analysis/time_attenuation/recalculate/CAL/cal.super.txt"))


ccix=col2rgb(timecol[5])/256
cci=c(rep(NA,4),rgb(ccix[1],ccix[2],ccix[3],alpha=0.5))
cx=cal_2panel_gg(
  list(cx1,cx2,cx3,cx4,cx5),
  as.character(test_times),
  col=timecol,
  ci_col=cci)

ggsave("Figures/pdfs/Unsorted/variable_score_cal.pdf",cx,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)




