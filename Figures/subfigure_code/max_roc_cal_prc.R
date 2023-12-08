source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/roc_max.txt"))

pdf("Figures/pdfs/Unsorted//roc_max.pdf",width=3,height=3.5)
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2); 
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2); 
aucm=signif(xrocm$auc[4],digits=3); sem=signif(xrocm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
#roc_2panel_gg(list(xroc3b,xroc4b,xrocmb),labels=labs,
#           col=c("blue","black","red"),xy_col="gray",
#           legend_title="AUROC (SE)")
#roc_2panel_gg(list(xrocmb),labels=labs[3],
#           col=c("blue","black","red")[2],xy_col="gray",
#           legend_title="AUROC (SE)")
roc_2panel_gg(list(xroc4b,xrocmb),labels=labs[2:3],
              col=c("red","black"),xy_col="gray",
              legend_title="AUROC (SE)")
dev.off()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/prc_max.txt"))

pdf("Figures/pdfs/Unsorted/prc_max.pdf",width=3,height=3.5)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2); 
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2); 
aucm=signif(xprcm$auc[4],digits=3); sem=signif(xprcm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
#prc_2panel_gg(list(xprc3b,xprc4b,xprcmb),labels=labs,col=c("blue","red","black"),
#           legend_title="AUPRC (SE)")
#prc_2panel_gg(list(xprcmb),labels=labs[3],col=c("black"),
#           legend_title="AUPRC (SE)")
prc_2panel_gg(list(xprc4b,xprcmb),labels=labs[2:3],col=c("red","black"),
           legend_title="AUPRC (SE)")
dev.off()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/cal_max.txt"))

pdf("Figures/pdfs/Unsorted/cal_max.pdf",width=3,height=3.5)
labs=c("v3", "v4","Max")
cic=c(rgb(0,0,1,alpha=0.5),rgb(1,0,1,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
#cal_2panel_gg(list(xcal3,xcal4,xcalm),labels=labs,col=c("blue","red","black"),
#             ci_col=cic,xy_col="gray")
#cal_2panel_gg(list(xcalm),labels=labs[3],col=c("blue","red","black")[3],
#             ci_col=cic[3],xy_col="gray")
cal_2panel_gg(list(xcal4,xcalm),labels=labs[2:3],col=c("red","black"),
             ci_col=cic[2:3],xy_col="gray")
dev.off()
