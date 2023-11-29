# Generate figures
source("Figures/subfigure_code/static_roc_prc_cal.R")
source("Figures/subfigure_code/auroc_comparison.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_staticscore/"
file.copy("Figures/pdfs/Unsorted/static_score_roc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/static_score_prc.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/static_score_cal.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/auroc_comparison_time.pdf",ff,overwrite=TRUE)



