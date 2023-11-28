

# Scripts
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel() etc

# Read in data
eval(import_sparra_expr("Analysis/full_model/Topics/topic_roc.txt"))
eval(import_sparra_expr("Analysis/full_model/Topics/topic_prc.txt"))
eval(import_sparra_expr("Analysis/full_model/Topics/topic_cal.txt"))

# Draw ROC
pdf("Figures/pdfs/Unsorted/topic_roc.pdf",width=3,height=3.5)
aucT=signif(xrocT$auc,digits=4); seT=signif(xrocT$se,digits=2);
aucNT=signif(xrocNT$auc,digits=4); seNT=signif(xrocNT$se,digits=2);
labs=c(paste0("Topics: ",aucT," (",seT,")"),
       paste0("No topics: ",aucNT," (",seNT,")"))
xcol=c("red","black")
roc_2panel_gg(list(xrocT,xrocNT),
              labels=labs,
              col=xcol,
              legend_title="AUROC (SE)",
              xy_col="blue")
dev.off()


### PR curves
pdf("Figures/pdfs/Unsorted/topic_prc.pdf",width=3,height=3.5)
aucT=signif(xprcT$auc,digits=4); seT=signif(xprcT$se,digits=2);
aucNT=signif(xprcNT$auc,digits=4); seNT=signif(xprcNT$se,digits=2);
labs=c(paste0("Topics: ",aucT," (",seT,")"),
       paste0("No topics: ",aucNT," (",seNT,")"))
xcol=c("red","black")
prc_2panel_gg(list(xprcT,xprcNT),
              labels=labs,
              col=xcol,
              legend_title="AUPRC (SE)")
dev.off()

## Calibration curves
pdf("Figures/pdfs/Unsorted/topic_cal.pdf",width=3,height=3.5)
labs=c("Topics","No topics")
xcol=c("red","black")
cci=rgb(0.5,0.5,0.5,alpha=0.5) # colour for confidence envelope
cal_2panel_gg(list(xcalT,xcalNT),
           labels=labs,
           col=xcol,
           ci_col=c(NA,cci),
           xy_col="blue")
dev.off()
