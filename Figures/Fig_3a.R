# Load the required libraries
library("tidyverse")
library("patchwork")
library("egg")

#setwd("~/SPARRAv4") # CAV: If you are using a project, this is not needed
source("Figures/util.R") # for import_sparra_expr()
#source("SPARRAv4/auxiliary.R") # for roc_2panel()

# conversion of figure 3a performance_by_age
plot_dir="Analysis/full_model/"
eval(import_sparra_expr(paste0(plot_dir, "Analytics/performance_by_age.txt")))

# Define the age group labels
age_labels <- c("(0, 20]", "(20, 30]", "(30, 40]", "(40, 50]",
                "(50, 60]", "(60, 70]", "(70, 80]", "(80, 120]")

# Create a data frame for plotting
plot_data <- data.frame(
  age = factor(rep(seq_len(ncol(perf)), times = 2)),  # x-axis values for positioning each version
  version = rep(c("V4", "V3"), each = ncol(perf)),  # Specify the versions
  auroc = c(perf[1, ], perf[2, ]),  # Repeat the AUROC values for each version
  sd = c(perf[3,], perf[4,]) # AUROC sd
)

# Plot per-group AUCs for both models
p3a1 <-
    ggplot(plot_data, aes(x = age, y = auroc, group = version)) +
    geom_point(shape = 18, size = 2, aes(col = version)) +  # Use versions as points
    geom_errorbar(aes(ymin = auroc - 3*sd, ymax = auroc + 3*sd, col = version),
                  width = 0.2, linewidth = 0.4) + # position = position_dodge(width = 1)
    xlab("Age group") +
    ylab("AUROC") +
    scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
    scale_x_discrete(labels = age_labels) +
    #scale_x_continuous(minor_breaks = NULL, breaks = 1:(length(age_labels)), labels = age_labels) +  # Set custom x-axis labels +  # Set custom x-axis labels
    # Set y-axis to start from 0
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.68, 0.81)) +
    guides(col = guide_legend(title = NULL)) +
    coord_fixed(ratio = 60) + # Set aspect ratio to 1:1
    theme_bw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8),
      #axis.text.y = element_text(angle = 90, hjust = 1, size = 7),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank() # Remove vertical grid lines
    )

# Sub-panels

plot_data2 <- data.frame(
  age = factor(seq_len(ncol(perf))),
  difference = plot_data$auroc[plot_data$version == "V4"] -
    plot_data$auroc[plot_data$version == "V3"],
  percentage_increase = ((plot_data$auroc[plot_data$version == "V4"] / plot_data$auroc[plot_data$version == "V3"]) - 1) * 100,
  prop_events = xrate
)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3a2 <-
  ggplot(plot_data2, aes(x = age, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Age group") +
  ylab("% increase") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 0.5) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 4.4)

# create a plot that shows the AUC differences as a percentage increase (v4 vs v3)
p3a3 <-
  ggplot(plot_data2, aes(x = age, y = prop_events)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Age group") +
  ylab("% events") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 7) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8)
  ) +
  ylim(0, 0.33)

p3a <- ggarrange( p3a1 + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank(),
                       plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.2, unit = "cm")),
           p3a2 + theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm")),
           p3a3 + theme(plot.margin = margin(t = 0, r = 0.2, l = 0.2, unit = "cm")),
           nrow = 3)

ggsave("Figures/pdfs/Fig_3a.pdf", p3a,
       width = 10, height = 20, units = "cm",
       device = cairo_pdf)

# old code
pdf(paste0(plot_dir,"Analytics/performance_by_age-new.pdf"),width=6,height=5)
labs=c();
age_split=c(0,20,30,40,50,60,70,80,120)
for (i in 1:(length(age_split)-1)) labs=c(labs,paste0("(",age_split[i],",",age_split[i+1],"]"))

par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=12; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
perfmin=min(perf[1:m0,])-0.02
plot(0,type="n",xlim=c(0,n0*(m0+2)),ylim=c(perfmin,max(perf)),xaxt="n",
     xlab="",ylab="AUROC") #,main="ROC and EA frequency by age band")
for (i in 1:n0) {
  segments((i-1)*(m0+2) + 1:m0+ 0.5,rep(0,m0),(i-1)*(m0+2) + 1:m0 + 0.5,perf[1:m0,i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]+perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]+perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
  segments((i-1)*(m0+2) + 1:m0+ 0.5-swidth,perf[1:m0,i]-perf[m0+(1:m0),i],
           (i-1)*(m0+2) + 1:m0 + 0.5+swidth,perf[1:m0,i]-perf[m0+(1:m0),i],col=xcol,lty=1,lwd=1)
}
lines((m0+2)*(-0.5 + 1:n0),perfmin + xrate/(xsc*max(xrate)),lty=2)

axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=2)
axis(4,at=perfmin+seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomright",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()
