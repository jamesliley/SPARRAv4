# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Analytics/n_events.txt")))

# old code

pdf(paste0(plot_dir,"Analytics/n_events.pdf"),width=5,height=3)
xmax=8; mh=0
gcol=gray((1:xmax)/(1.5*xmax))
xcol=c("red",gcol)

# Truncate for ease of viewing
xfreq1=xfreq
xfreq1[,6]=rowSums(xfreq[,6:dim(xfreq)[2]])
xfreq1=xfreq1[,1:6]

m0=dim(xfreq1)[1]; n0=dim(xfreq1)[2]
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(0,max(xfreq1)),bty="n",xaxt="n",
     xlab="N admissions",ylab="% patients")
for (i in 1:dim(xfreq1)[2])
  segments((i-1)*(m0+2) + 1:m0,rep(0,m0),(i-1)*(m0+2) + 1:m0,xfreq1[,i],col=xcol,lty=1,lwd=1)
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=c(1:(n0-1),paste0(">",n0-1)))

legend("topright",lty=1,col=c(gcol[1],"white",gcol[xmax],"red"),
       c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"),
       lwd=2)

dev.off()