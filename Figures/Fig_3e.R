# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="Analysis/full_model/"
## If running from loaded data, first run
lscore=c(0,0); lscore0=c(0,0); lscoreeb=c(0,0); lscored=c(0,0)
eval(import_sparra_expr(paste0(plot_dir,"Analytics/disease_class_violin.txt")))

# old code
pdf(paste0(plot_dir,"Analytics/disease_class_violin.pdf"),width=6,height=5)

col1="blue" # Colour for categories

par(mar=c(8.1,4.1,4.1,2.1))
dmax=0
for (i in c(1:length(xmed))) {
  dz=get(paste0("d",i))
  dmax=c(dmax,max(dz$y))
}
space1=4/3
space2=1

sc=0.3/max(dmax)

plot(0,type="n",xlim=c(0,4*space1 + length(xmed)*space2),ylim=c(0,2),bty="n",xaxt="n",yaxt="n",yaxs="i",
     xlab="",ylab="Score (%), log scale")

locs=c()
for (i in 1:4) {
  suff = c("0","x","eb","d")[i]
  
  d=get(paste0("d",suff))
  m=get(paste0("xmed",suff))
  polygon(c(sc*d$y,rev(-sc*d$y))+(i-0.5)*space1,c(d$x,rev(d$x)),col=col1,border=NA)
  points((i-0.5)*space1,m,pch=16)
  locs=c(locs,(i-0.5)*space1)
}
ofs=(i+0.5)*space1

ox=order(-xmed)
for (i0 in 1:length(xl0)) {
  i=ox[i0]
  dz=get(paste0("d",i))
  polygon(c(sc*dz$y,rev(-sc*dz$y))+ofs + (i0-0.5)*space2,c(dz$x,rev(dz$x)),col="gray",border=NA)
  points(ofs + (i0-0.5)*space2,xmed[i],pch=16)
  locs=c(locs,ofs + (i0-0.5)*space2)
}
yp=pretty(c(0,2),n=5)
axis(2,at=yp,labels=round((10^yp)),las=2,cex.axis=0.5)
axis(1,at=locs[1:4],
     labels=c("No admission/death","All admitted/died","Admitted","Died"),
     col=col1, # suppress line
     col.ticks=col1,
     col.axis=col1,
     las=2,cex.axis=0.5)
axis(1,at=locs[5:length(locs)],
     labels=xl0[ox],
     col="black", # suppress line
     col.ticks="black",
     col.axis="black",
     las=2,cex.axis=0.5)


dev.off()