## Age vs Shapley value for age
##
#Pipelines/main_pipeline.R

source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Shapley_values/age_effect.txt")))

# old code
pdf(paste0(plot_dir,"Shapley_values/age_effect.pdf"),width=6,heigh=5)
xrange=c(0,85)
icol=gray(1-(1:100)/200)

plot(0,type="n",xlim=xrange,ylim=yrange,bty="n",
     xlab="Age",ylab="Shapley value: Age") #,main="Importance of age as predictor")
image(scx,add=T,col=icol)
#plot(xx,yy,col="gray",cex=0.5) # leave this out for the moment for privacy reasons

lines(xtrue,ytrue)
lines(xtrue,ytrue + sqrt(yvar),lty=2)
lines(xtrue,ytrue - sqrt(yvar),lty=2)

# split by sex: uninteresting
#lines(xtruem,ytruem,col="red")
#lines(xtruef,ytruef,col="blue")

legend("topleft",c("Mean","SD","Density"),lty=c(1,2,NA),
       pch=c(NA,NA,16),pt.cex=c(NA,NA,2),col=c("black","black","gray"),bty="n")

dev.off()