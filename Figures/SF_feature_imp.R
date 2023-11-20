# Generate figures
source("Figures/subfigure_code/topic_roc_prc_cal.R")
source("Figures/subfigure_code/shapley_value_plots.R")
source("Figures/subfigure_code/n_adm_effect.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/SF_feature_imp/"
file.copy("Figures/pdfs/Unsorted/topic_roc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/topic_prc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/topic_cal.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/shapley_ltc_FIRST_COPD_EPISODE_yearssincediag_small_range.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/shapley_days_since_last_emergency_admission_large_range.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/n_adm_effect.pdf",ff,overwrite=TRUE)

