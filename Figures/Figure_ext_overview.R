# Generate figures
source("Figures/subfigure_code/age_simd_event_exclusion_distributions.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_ext_overview//"
file.copy("Figures/pdfs/Unsorted/age_simd_distributions.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/event_type.pdf",ff,overwrite=TRUE)


