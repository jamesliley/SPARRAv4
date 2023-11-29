
# Need case_when
library(dplyr)

# Scripts
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel() etc

# Read data
lxs=list.files("Analysis/full_model/Shapley_values/All/",pattern="*small_range.txt",full.names=TRUE)
lxr=list.files("Analysis/full_model/Shapley_values/All/",pattern="*large_range.txt",full.names=TRUE)
lx=c(lxs,lxr)

# Ranges
small_range=c(-0.02,0.06)
large_range=c(-0.02,0.15)

# Draw figures
for (i in 1:length(lx)) {
  eval(import_sparra_expr(lx[i]))
  
  if (lx[i] %in% lxr) {
    yvar=large_range 
    suffix="large_range"
  } else {
    yvar=small_range  
    suffix="small_range"
  }
  
  xname=gsub("_small_range.txt","",gsub("_large_range.txt","",gsub("shapley_","",basename(lx[i]))))
  
  
  # Forgot to export labels when exporting Shapley value plots. Fill in manually for a few.
  lab=x0
  if (xname=="days_since_last_emergency_admission") {
    lab=c(
      "1 - 110",
      "110 - 210",
      "210 - 320",
      "320 - 430",
      "430 - 540",
      "540 - 640",
      "640 - 760",
      "760 - 870",
      "870 - 980",
      "980 - 1100"
    )
  }
  if (xname=="ltc_FIRST_COPD_EPISODE_yearssincediag") {
    lab=c(
      "0 - 2",
      "3 - 5",
      "6 - 7",
      "8 - 10",
      "11 - 12",
      "13 - 15",
      "16 - 18"
    )
  }
  
  
  pdf(paste0("Figures/pdfs/Unsorted/",gsub("txt","pdf",basename(lx[i]))),width=7,height=6)
  par(mar=c(7,4.1,4.1,2.1))
  plot(0,type="n",xlim=range(x0),ylim=yvar,main=longvarnames(xname),xaxt="n",
       xlab="",ylab="Shapley value") # note fixed ylim
  lines(x0,mx,lty=2,col="black",lwd=2)
  segments(x0,mx-sdx,x0,mx+sdx,col="black",lwd=2)
  points(x0,mx,pch=16,col="black")
  axis(1,at=x0,labels=lab,las=2)
  abline(h=0,lty=2)
  dev.off()
}
