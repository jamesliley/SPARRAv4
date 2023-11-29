library("tidyverse")
library("patchwork")

source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()

eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/roc_max.txt"))

# ROC ####
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2);
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"))


# New Figure 2(a)
groc=roc_2panel_gg(list(xroc3b,xroc4b),labels=labs,
           col=c("blue","red"),xy_col="gray",
           legend_title="AUROC (SE)")


sc=1.3
ggsave("Figures/pdfs/Unsorted/roc_v3v4.pdf",plot=groc,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)






####### PRC #####
eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/prc_max.txt"))



# Old Figure 2(b)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2);
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"))



# New Figure 2(b)
gprc=prc_2panel_gg(list(xprc3b,xprc4b),labels=labs,col=c("blue","red"),
                legend_title="AUPRC (SE)")


ggsave("Figures/pdfs/Unsorted/prc_v3v4.pdf",plot = gprc,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)





### CAL #####

eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/cal_max.txt"))



# Old Figure 2(c)
labs=c("v3", "v4")
cic=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5))


# New Figure 2(c)
gcal=cal_2panel_gg(list(xcal3,xcal4),labels=labs,col=c("blue","red"),
                   ci_col=cic,xy_col=NA)


ggsave("Figures/pdfs/Unsorted/cal_v3v4.pdf",plot=gcal,
       width = sc*7.5, height = sc*7.25, units = "cm",
       device = cairo_pdf)
