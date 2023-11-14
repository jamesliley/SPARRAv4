# Generate figures
source("Figures/subfigure_code/variable_roc_prc_cal.R")
source("Figures/subfigure_code/time_comparison_density.R")
source("Figures/subfigure_code/change_over_time.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_time_att/"
file.copy("Figures/pdfs/Unsorted/time1_vs_time2.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/variable_score_roc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/variable_score_prc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/variable_score_cal.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/cohort_change_over_time.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/individual_change_over_time.pdf",ff,overwrite=TRUE)



