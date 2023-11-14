source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/roc_max.txt"))

pdf("Figures/pdfs/Unsorted//roc_max.pdf",width=3,height=3.5)
auc3=signif(xroc3$auc[4],digits=3); se3=signif(xroc3$se[4],digits=2); 
auc4=signif(xroc4$auc[4],digits=3); se4=signif(xroc4$se[4],digits=2); 
aucm=signif(xrocm$auc[4],digits=3); sem=signif(xrocm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
roc_2panel(list(xroc3b,xroc4b,xrocmb),labels=labs,
           col=c("blue","black","red"),xy_lty=2,xy_col="gray",
           text.col=c("blue","red","black"),title="AUROC (SE)",
           title.col="black")
dev.off()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/prc_max.txt"))

pdf("Figures/pdfs/Unsorted/prc_max.pdf",width=3,height=3.5)
auc3=signif(xprc3$auc[4],digits=3); se3=signif(xprc3$se[4],digits=2); 
auc4=signif(xprc4$auc[4],digits=3); se4=signif(xprc4$se[4],digits=2); 
aucm=signif(xprcm$auc[4],digits=3); sem=signif(xprcm$se[4],digits=2);
labs=c(paste0("v3: ",auc3," (",se3,")"),paste0("v4: ",auc4," (",se4,")"),
       paste0("Max: ",aucm," (",sem,")"))
prc_2panel(list(xprc3b,xprc4b,xprcmb),labels=labs,col=c("blue","red","black"),
           title="AUPRC (SE)",title.col="black",text.col=c("blue","red","black"))
dev.off()



eval(import_sparra_expr("Analysis/full_model/Performance/v4_final/cal_max.txt"))

pdf("Figures/pdfs/Unsorted/cal_max.pdf",width=3,height=3.5)
labs=c("v3", "v4","Max")
cic=c(rgb(0,0,1,alpha=0.5),rgb(1,0,1,alpha=0.5),rgb(0.5,0.5,0.5,alpha=0.5))
cal_2panel(list(xcal3,xcal4,xcalm),labels=labs,col=c("blue","red","black"),
           text.col=c("blue","red","black"),ci_col=cic,xy_col="gray",xy_lty=2)
dev.off()
