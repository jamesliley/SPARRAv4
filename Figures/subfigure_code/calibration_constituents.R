# Scripts
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel() etc


px=list.files("Analysis/full_model/Performance/constituent_predictors/",pattern="roc*")
px=px[grep("pdf",px)]
px=setdiff(gsub("roc_","",gsub(".pdf","",px)),"all")



# Separate
for (i in 1:length(px)) {
  
  eval(import_sparra_expr(paste0("Analysis/full_model/Performance/constituent_predictors/cal_",px[i],".txt")))
  
  # Calibration
  pdf(paste0("Figures/pdfs/Unsorted/cal_",px[i],".pdf"),width=3,height=3.5)
  cal_2panel_gg(list(xcal1,xcal2,xcal3),labels=1:3,col=c("black","red","blue"),legend_title ="CV fold" )
  dev.off()
  
}




