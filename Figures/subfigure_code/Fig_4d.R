
# Pipelines/time_attenuation.R
plot_dir="Analysis/time_attenuation/"

source("Figures/util.R") # for import_sparra_expr()

eval(import_sparra_expr(paste0(plot_dir,"cohort/change/cohort_change_over_time.txt")))

# old code
pdf(paste0(plot_dir,"cohort/change/cohort_change_over_time.pdf"),width=4,height=4)