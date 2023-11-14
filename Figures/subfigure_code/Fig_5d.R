#Pipelines/main_pipeline.R

source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Shapley_values/age_simd_direct_effect.txt")))

# old code

pdf(paste0(plot_dir,"Shapley_values/age_simd_direct_effect.pdf"),width=5,height=5)
plot(0,type="n",xlim=range(xx),xlab="Age",ylab="EA freq.",col="red",ylim=range(ysimd))
for (i in 1:10) lines(xx,ysimd[,i],col=gray(i/20))
lines(xx,yy,col="red")
#lines(xx,yyf,col="red")
legend("topleft",c("All","1","...","10"),lty=c(1,1,NA,1),col=c("red",gray(1/20),NA,gray(1/2)))

dev.off()