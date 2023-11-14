# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Analytics/disease_class.txt")))

# old code

pdf(paste0(plot_dir,"Analytics/disease_class.pdf"),width=8,height=5)
par(mar=c(8.1,4.1,4.1,2.1))
xmax=8; mh=0
gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("blue",gcol) # Changed from 'red' in response to reviewer suggestion

# Truncate xfreq for plotting simplicity
xfreq1=xfreq
xfreq1[,14]=rowSums(xfreq[,14:dim(xfreq)[2]])
xfreq1=xfreq1[,1:14]

xlabs=colnames(xfreq1)

m0=dim(xfreq1)[1]; n0=dim(xfreq1)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq1)),bty="n",xaxt="n",
     xlab="",ylab="% patients")
for (i in 1:dim(xfreq1)[2]) {
  xcol2=xcol;
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq1[,i],col=xcol,lty=1,lwd=1)
}
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=xlabs,las=2,cex.axis=0.5)

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"blue"),
       c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
       lwd=2)

dev.off()