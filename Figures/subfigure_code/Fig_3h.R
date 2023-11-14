# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"
# To recreate: first set
x=c(1,1); tmax=1
eval(import_sparra_expr(paste0(plot_dir,"Analytics/time_to_only_event.txt")))

# old code

pdf(paste0(plot_dir,"Analytics/time_to_only_event.pdf"),width=5,height=5)

xmax=8
mh=0; for (i in 1:xmax) mh=max(mh,get(paste0("d",i))$y)


plot(0,type="n",xlim=c(0,365),xlab="Days after time cutoff",ylim=c(0,mh),
     ylab="Density",#main="Density of time-to-only-admission",
     yaxt="n")

gcol=gray((1:xmax)/(1.5*xmax))
for (i in 1:xmax) lines(get(paste0("d",i)),col=gcol[i],lty=1,lwd=2)
lines(d0,col="red",lwd=2)

legend("topright",lty=1,lwd=2,col=c(gcol[1],"white",gcol[xmax],"red"),
       c(expression(paste(hat(P),"(EA)>0.1")),"...",bquote(paste(hat(P),"(EA)>",.(xmax/10))),"All"))

dev.off()