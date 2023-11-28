library("tidyverse")
library("patchwork")
library("glue")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v3_vs_v4/v3_vs_v4_calibration_differential.txt"))



# Old Figure 2(d)
del=0.1
cic=c(rgb(1,0,0,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(cv3h3,cv3h4,cv4h3,cv4h4),labels=c("v3, v3>v4","v3, v4>v3","v4, v3>v4", "v4, v4>v3"),
           title=paste0("|v3-v4| > ",del),xy_col="gray",xy_lty=2,col=c("red","red","black","black"),
           lty=c(1,2,1,2),ci_col=cic)

# New figure
cic2=c(rgb(0,0,1,alpha=0.5),rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(1,0,0,alpha=0.5))
cx=cal_2panel_gg(list(cv3h3,cv3h4,cv4h3,cv4h4),
              labels=c("v3, v3>v4","v3, v4>v3","v4, v3>v4", "v4, v4>v3"),
           legend_title=" ", #paste0("|v3-v4| > ",del),
           xy_col="gray",
           col=c("blue","blue","red","red"),
           ci_col=cic)

sc=1.5
ggsave("Figures/pdfs/Unsorted/calibration_differential.pdf",cx,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)

