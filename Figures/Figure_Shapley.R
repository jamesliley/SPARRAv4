# Generate figures
source("Figures/subfigure_code/shapley_age_simd.R")
source("Figures/subfigure_code/age_simd.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/Figure_Shapley/"
file.copy("Figures/pdfs/Unsorted/age_effect.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/simd_effect.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/equivalent_years_older_shapley.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/age_simd_equivalent.pdf",ff,overwrite=TRUE)



