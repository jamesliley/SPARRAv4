#Pipelines/main_pipeline.R

source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Shapley_values/age_simd_equivalent.txt")))

# old code
pdf(paste0(plot_dir,"Shapley_values/age_simd_equivalent.pdf"),width=5,height=5)
plot(0,type="n",xlab="Chron. age",ylab="Effective age",xlim=c(15,max(xx)),ylim=c(20,max(xx)))
for (i in 1:10) lines(xx,yeq[,i],col=xcol[i])
abline(0,1,col="black",lty=2)
legend("bottomright",c("1","...","10"),title="SIMD",lty=1,col=c(xcol[1],xcol[5],xcol[10]))
dev.off()