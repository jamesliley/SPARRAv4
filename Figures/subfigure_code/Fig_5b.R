## SIMD vs Shapley value for SIMD
##
#Pipelines/main_pipeline.R

source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Shapley_values/simd_effect.txt")))

# old code
pdf(paste0(plot_dir,"Shapley_values/simd_effect.pdf"),width=6,heigh=5)
yrange=mean(mx) + rcex*(range(mx)-mean(mx))
plot(0,type="n",xlim=c(1,11),ylim=yrange,bty="n",
     xlab="SIMD",ylab="Shapley value: SIMD") #,main="Importance of SIMD as predictor")

for (i in 1:10) polygon(i+dxy[i,],dxx[i,],col="gray",border=NA)
points(1:10,mx,pch=16)

legend("topright",c("Mean","Density"),pch=c(16,16),col=c("black","gray"),bty="n")

dev.off()