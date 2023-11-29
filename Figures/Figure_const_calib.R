# Generate figures
source("Figures/subfigure_code/calibration_constituents.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_const_calib/"
fx=list.files("Figures/pdfs/Unsorted/",pattern="cal_*",full.names=TRUE)
for (i in 1:length(fx)) file.copy(fx[i],ff,overwrite=TRUE)
