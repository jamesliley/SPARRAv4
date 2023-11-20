# Generate figures
source("Figures/subfigure_code/stratified_v3v4.R")
source("Figures/subfigure_code/admission_type_violin.R")
source("Figures/subfigure_code/time_to_first_event.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_stratified_perf/"
file.copy("Figures/pdfs/Unsorted/strat_by_age.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/strat_by_simd.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/strat_by_v3.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/disease_class_violin.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/time_to_first_event.pdf",ff,overwrite=TRUE)


