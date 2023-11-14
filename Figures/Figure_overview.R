# Generate figures
source("Figures/subfigure_code/age_simd_event_exclusion_distributions.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_overview/"
file.copy("Diagrams/SPARRA_scriberia.jpeg",ff,overwrite=TRUE)
file.copy("Diagrams/flow.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/age_simd_density.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/exclusions_venn.pdf",ff,overwrite=TRUE)
####### NEED TO ADD UP TO DATE FLOW DIAGRAM 



