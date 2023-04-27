# Pipelines/time_attenuation.R
plot_dir="~/SPARRAv4/Analysis/time_attenuation/"

source("Figures/util.R") # for import_sparra_expr()

eval(import_sparra_expr(paste0(plot_dir,"cohort/change/individual_change_over_time.txt")))

# old code
pdf(paste0(plot_dir,"cohort/change/individual_change_over_time.pdf"),width=4,height=4)
par(mar=c(4,4,0.5,0.5))
plot(0,type="n",xlim=c(1,ntx),ylim=c(0,0.05+max(qmat)),
     xaxt="n",xlab="Time cutoff",
     ylab="P(EA)")
for (i in 0:10) abline(h=i/10,col="lightgray",lty=2)
axis(1,labels=c(expression("t"[1]),expression("t"[2]),expression("t"[3]),expression("t"[4]),expression("t"[5])),
     at=1:length(test_times))
for (i in 1:nquant) lines(1:ntx,qmat[i,],col="lightblue",lty=1)
for (i in 1:10) lines(1:ntx,qmat[round(nquant*i/10),],col="blue",lty=1)
dev.off()