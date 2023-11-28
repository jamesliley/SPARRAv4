# Generate figures
source("Figures/subfigure_code/v34_density.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_v34_density/"
file.copy("Figures/pdfs/Unsorted/v3_v4_density.pdf",ff,overwrite=TRUE)


