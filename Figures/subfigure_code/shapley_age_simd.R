# Scripts
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel() etc

## Age ####

eval(import_sparra_expr("Analysis/full_model/Shapley_values/age_effect.txt"))

pdf("Figures/pdfs/Unsorted/age_effect.pdf",width=6,heigh=5)
xrange=c(0,85); yrange=c(-0.02,0.06)
icol=gray(1-(1:100)/200)

plot(0,type="n",xlim=xrange,ylim=yrange,bty="n",
     xlab="Age",ylab="Shapley value: Age") #,main="Importance of age as predictor")
image(scx,add=T,col=icol) 

lines(xtrue,ytrue,lwd=2)
lines(xtrue,ytrue + sqrt(yvar),lty=2,lwd=2)
lines(xtrue,ytrue - sqrt(yvar),lty=2,lwd=2)

legend("topleft",c("Mean","Mean +/- SD","Density"),lty=c(1,2,NA),lwd=2,
       pch=c(NA,NA,16),pt.cex=c(NA,NA,2),col=c("black","black","gray"),bty="n")

dev.off()




### SIMD ####


eval(import_sparra_expr("Analysis/full_model/Shapley_values/simd_effect.txt"))

# Forgot to export dxx; just copy plot
file.copy("Analysis/full_model/Shapley_values/simd_effect.pdf",
          "Figures/pdfs/Unsorted/",overwrite=TRUE)
