# Generate figures
source("Figures/subfigure_code/v3v4_roc_prc_cal.R")
source("Figures/subfigure_code/calibration_differential.R")
source("Figures/subfigure_code/v3v4_difference_topN.R")
source("Figures/subfigure_code/v3v4_nnt.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_v3v4/"
file.copy("Figures/pdfs/Unsorted/roc_v3v4.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/prc_v3v4.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/cal_v3v4.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/calibration_differential.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/v3v4_difference_topN.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/v3v4_nnt.pdf",ff,overwrite=TRUE)
