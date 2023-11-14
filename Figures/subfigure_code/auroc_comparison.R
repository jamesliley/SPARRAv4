library(matrixStats)
library(lubridate)

source("Figures/util.R") # for import_sparra_expr()

eval(import_sparra_expr("Analysis/time_attenuation/cohort/auroc_comparison.txt"))


# Recent vs earlier predictions
pdf("Figures/pdfs/Unsorted/auroc_comparison_time.pdf",width=5,height=5)
ntx=5
cols=gray(((1:ntx)-1)/(1.3*ntx))
ltyx=c(1,2,3,4,5)

plot(0,type="n",xlim=c(1,ntx),ylim=c(0.755,0.795),xlab="Time point",ylab=expression(paste("AUROC")))
for (i in 1:ntx) {
  if (i<ntx) lines(i:ntx,dx[i,i:ntx],col=cols[i],lty=ltyx[i]) else points(ntx,dx[ntx,ntx],col=cols[i],pch=16)
} 
legend("bottomleft",title="Predictor for:", legend=1:5, #as.character(test_times),
       col=cols,lty=c(1:4,NA),pch=c(rep(NA,ntx-1),16))
dev.off()