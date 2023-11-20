# Generate figures
source("Figures/subfigure_code/max_roc_cal_prc.R")

# Move ROC/PRC/CAL files
ff="Figures/pdfs/SF_updating/"
file.copy("Diagrams/causality_diagram_panel1.pdf",ff,overwrite=TRUE)
file.copy("Diagrams/causality_diagram_panel2.pdf",ff,overwrite=TRUE)
file.copy("Diagrams/causality_diagram_panel3.pdf",ff,overwrite=TRUE)
file.copy("Diagrams/causality_diagram_panel4.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/roc_max.pdf",ff,overwrite=TRUE)
file.copy("Figures/pdfs/Unsorted/cal_max.pdf",ff,overwrite=TRUE)

