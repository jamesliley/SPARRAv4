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
  y = c(perf[1:m0, ], perf[1:m0, ]),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)

sd_values <- c(0.0006950962, 0.0007206066, 0.0007468640, 0.0007848581, 0.0007990873, 0.0008351305, 0.0008600128, 0.0008998324, 0.0009158358, 0.0009865907,
  0.0007128243, 0.0007392061, 0.0007657252, 0.0008038775, 0.0008194903, 0.0008559886, 0.0008828196, 0.0009234900, 0.0009403307, 0.0010172065)
plot_data$sd <- sd_values
# Set up the plot
p3b1_min <-
  ggplot(plot_data, aes(x = x, y = y, group = version)) +
  geom_point(shape = 18, size = 2, aes(col = version)) +
  geom_errorbar(aes(ymin = y - 3*sd, ymax = y + 3*sd, col = version), width = 0.2, size = 0.4) + # position = position_dodge(width = 1)
  scale_x_discrete(breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1, labels = labs) +
  #scale_y_continuous(expand = c(0, 0), limits = c(min(plot_data$y)-0.1, 0.9)) +
  scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "SIMD", y = "AUROC") +
  guides(col = guide_legend(title = NULL)) +
  coord_fixed(ratio = 200) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )


ggsave("Figures/pdfs/Fig_3b1.pdf", p3b1_min,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)

p3b1_zero <-
  ggplot(plot_data, aes(x = x, y = y, group = version)) +
  geom_point(shape = 18, size = 3, aes(col = version)) +
  scale_x_discrete(breaks = (m0 + 2) * (1:n0) - floor(m0 / 2) - 1, labels = labs) +
  #scale_y_continuous(expand = c(0, 0), limits = c(min(plot_data$y)-0.1, 0.9)) +
  scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x = "SIMD", y = "AUROC") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# ggsave("Figures/pdfs/Fig_3b1.pdf", p3b1_zero,
#        width = 8.5, height = 9, units = "cm",
#        device = cairo_pdf)


# difference in AUC
plot_data <- data.frame(
  x = factor(rep((m0 + 2) * (1:n0) - floor(m0 / 2) - 1, each = 1)),
  y = c(perf[1:m0, ], perf[1:m0, ]),
  group = factor(rep(1:n0, each = m0)),
  version = rep(c("V4", "V3"), times = n0)
)

plot_data$difference <- plot_data$y[plot_data$version == "V4"] -
  plot_data$y[plot_data$version == "V3"]

# Set up the plot
p3b2 <-
  ggplot(plot_data, aes(x= x, y = difference)) +
  geom_point(shape = 18, size = 3) +
  scale_x_discrete(labels = labs) +
  xlab("SIMD") +
  ylab("Difference in AUC") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Figures/pdfs/Fig_3b2.pdf", p3b2,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)



# Percentage increase
plot_data$percentage_increase <- ((plot_data$y[plot_data$version == "V4"] / plot_data$y[plot_data$version == "V3"]) - 1) * 100


p3b3 <-
  ggplot(plot_data, aes(x = x, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("SIMD") +
  ylab("Percentage Increase") +
  scale_x_discrete(labels = labs) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7)
  )


ggsave("Figures/pdfs/Fig_3b3.pdf", p3b3,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)


diffs <- data.frame(
  x = (m0 + 2)*(-0.5 + 1:n0) - 0.5*(m0 + 2)/n0,
  y = round(xrate, 3),
  group = factor(rep(1, n0))
)

freq <- ggplot(diffs, aes(x = x, y = y, group = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  #geom_bar(stat = "identity", color = "black", fill = "red") +
  scale_x_continuous(limits = c(0, n0*(m0 + 2)),
                     breaks = (m0 + 2)*(1:n0) - floor(m0/2) - 1 ,
                     labels = labs) +
  scale_y_continuous(limits = c(0, max(diffs$y) + 0.01),
                     breaks = seq(0, max(diffs$y), length.out = 5),
                     labels = function(x) round(x, 2),
                     expand = expansion(mult = c(0.0000005, 0.00000005))) +
  labs(x = "Cohort", y = "EA Frequency") +
  coord_fixed(ratio = 300) +  # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7)
  )


ggsave("Figures/pdfs/Fig_3b4.pdf", freq,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)

