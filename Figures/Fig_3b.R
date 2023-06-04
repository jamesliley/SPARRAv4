#Pipelines/main_pipeline.R
source("Figures/util.R") # for import_sparra_expr()
eval(import_sparra_expr("~/SPARRAv4/Analysis/full_model/Analytics/performance_by_simd.txt"))

# old code
plot_dir="~/SPARRAv4/Analysis/full_model/"

pdf(paste0(plot_dir,"Analytics/performance_by_simd.pdf"),width=6,height=5)
labs=1:10

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=25; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
     xlab="SIMD",ylab="AUROC") #,main="ROC and EA frequency by SIMD")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}
lines((m0+2)*(-0.5 + 1:n0),perfmin + xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=1)
axis(4,at=perfmin + seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomright",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()

#### new version of the figure
library(ggplot2)

labs <- 1:10

xsc <- 25
swidth <- 0.1
m0 <- dim(perf)[1] / 2
n0 <- dim(perf)[2]
perfmin <- min(perf[1:m0, ]) - 0.02

# Create a data frame for plotting
plot_data <- data.frame(
  x = factor(rep((m0 + 2) * (1:n0) - floor(m0 / 2) - 1, each = m0)),
  y = c(t(perf[1:m0, ])),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)

# Set up the plot
p3b <-
  ggplot(plot_data, aes(x = x, y = y, group = group, color = version)) +
  geom_point(shape = 18, size = 3) +
  scale_x_discrete(breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1, labels = labs) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "SIMD", y = "AUROC") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )


ggsave("Figures/pdfs/Fig_3b.pdf", p3b,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)