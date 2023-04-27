# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"
## If running from loaded data, first run
lscore=c(0,0); lscore0=c(0,0)
eval(import_sparra_expr(paste0(plot_dir,"Analytics/disease_class_violin.txt")))

# old code
pdf(paste0(plot_dir,"Analytics/disease_class_violin.pdf"),width=6,height=5)

par(mar=c(8.1,4.1,4.1,2.1))
dmax=0
for (i in 1:length(xmed)) {
  dz=get(paste0("d",i))
  dmax=c(dmax,max(dz$y))
}
sc=0.3/max(dmax)

plot(0,type="n",xlim=c(0,2+length(xl0)),ylim=c(0,2),bty="n",xaxt="n",yaxt="n",
     xlab="",ylab="Score (%), log scale")

polygon(c(sc*d0$y,rev(-sc*d0$y))+1/3,c(d0$x,rev(d0$x)),col="red",border=NA)
points(1/3,xmed0,pch=16)

polygon(c(sc*dx$y,rev(-sc*dx$y))+5/3,c(dx$x,rev(dx$x)),col="red",border=NA)
points(5/3,xmedx,pch=16)

ox=order(-xmed)
for (i0 in 1:length(xl0)) {
  i=ox[i0]
  dz=get(paste0("d",i))
  polygon(c(sc*dz$y,rev(-sc*dz$y))+i0+2,c(dz$x,rev(dz$x)),col="gray",border=NA)
  points(i0+2,xmed[i],pch=16)
}
axis(1,at=c(1/3,5/3,3:(2+length(xl0))),labels=c("Not admitted","All admitted",xl0[ox]),las=2,cex.axis=0.5)
yp=pretty(c(0,2),n=5)
axis(2,at=yp,labels=round((10^yp)),las=2,cex.axis=0.5)
dev.off()