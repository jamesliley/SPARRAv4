# Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
plot_dir="~/SPARRAv4/Analysis/full_model/"

eval(import_sparra_expr(paste0(plot_dir,"Analytics/performance_by_v3_cohort.txt")))

# old code

pdf(paste0(plot_dir,"Analytics/performance_by_v3_cohort.pdf"),width=6,height=5)
labs=c("FEC","LTC","YED","U16") #,"LTC,YED")

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=20; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
     xlab="Cohort",ylab="AUROC") #,main="ROC and EA frequency in v3 cohorts")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}
lines((m0+2)*(-0.5 + 1:n0),perfmin+xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=1)
axis(4,at=perfmin + seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomleft",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()

### new version
labs <- c("FEC", "LTC", "YED", "U16")

xcol <- c("black", "red")
xsc <- 20
swidth <- 0.1
m0 <- dim(perf)[1] / 2
n0 <- dim(perf)[2]
perfmin <- min(perf[1:m0, ]) - 0.02

# Create a data frame for plotting
plot_data <- data.frame(
  x = rep((m0 + 2) * (1:n0) - floor(m0 / 2) - 1, each = m0),
  y = c(t(perf[1:m0, ])),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)

# Set up the plot
ggplot(plot_data, aes(x = x, y = y, group = group)) +
  geom_point(aes(color = version), shape = 18, size = 3) +
  scale_x_continuous(
    breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1,
    labels = labs,
    expand = c(0, 0), limits = c(0, 15),  # Adjust the expand argument to extend the x-axis limits
    position = "bottom"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +  # Extend the y-axis to 1
  labs(x = "Cohort", y = "AUROC") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
