# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="Analysis/full_model/"
# To recreate: first set
x=c(1,1); tmax=1

eval(import_sparra_expr(paste0(plot_dir,"Analytics/time_to_first_event.txt")))

# old code
# We mostly predict imminent events
pdf("Figures/pdfs/Unsorted/time_to_first_event.pdf",width=5,height=5)


xmax=8
mh=0; for (i in 1:xmax) mh=max(mh,get(paste0("d",i))$y)

plot(0,type="n",xlim=c(0,365),xlab="Days after time cutoff",ylim=c(0,mh),
     ylab="Density", #main="Density of time-to-first-admission",
     yaxt="n")

gcol=gray((1:xmax)/(1.5*xmax))
sub=c(1,3,5,7)
for (i in sub) lines(get(paste0("d",i)),col=gcol[i],lty=1,lwd=2)
lines(d0,col="black",lwd=4)

legend("topright",lty=1,lwd=c(2,2,2,4),col=c(gcol[1],"white",gcol[xmax],"black"),,
       c("v4 > 0.1","...","v4 > 0.8","All"))

dev.off()