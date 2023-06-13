# Load the required libraries
library("tidyverse")
library("patchwork")

#setwd("~/SPARRAv4") # CAV: If you are using a project, this is not needed
source("Figures/util.R") # for import_sparra_expr()
#source("SPARRAv4/auxiliary.R") # for roc_2panel()

# conversion of figure 3a performance_by_age
plot_dir="Analysis/full_model/"
eval(import_sparra_expr(paste0(plot_dir, "Analytics/performance_by_age.txt")))

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

#### new version of the figure
library(ggplot2)

# Define the age group labels
age_labels <- c("(0, 20]", "(20, 30]", "(30, 40]", "(40, 50]", "(50, 60]", "(60, 70]", "(70, 80]", "(80, 120]")

# Create a data frame for plotting
plot_data <- data.frame(
  age_group = factor(rep(1:(length(age_split)-1), each = 1)),
  version = rep(c("V4", "V3"), times = length(age_split)-1),  # Specify the versions
  auroc = c(perf[1:m0, ], perf[1:m0, ]),  # Repeat the AUROC values for each version
  x = rep(1:(length(age_split)-1), each = 2)  # x-axis values for positioning each version
)

sd_values <- c(0.0008363494, 0.001016852, 0.0009855209, 0.0008637745, 0.0007508718, 0.0006955779, 0.0006557762, 0.0006895856,
               0.0008486927, 0.001045860, 0.0010092724, 0.0008837348, 0.0007658613, 0.0007116656, 0.0006741307, 0.0007079138)

# Update the plot_data dataframe
plot_data$sd <- sd_values

# # Calculate the maximum label length of the three plots which will be combined
# max_label_length1 <- max(nchar(as.character(plot_data$auroc)))
# max_label_length2 <- max(nchar(as.character(plot_data$percentage_increase)))
# max_label_length3 <- max(nchar(as.character(diffs$y)))
#
# # Set the maximum label length (in inches) for both plots
# max_label_length <- max(max_label_length1, max_label_length2, max_label_length3) * 0.25  # Adjust the multiplier as needed
#
# # Set up the main plot
# # The errorbars overlapped, so use position_dodge to move them horizontally
# pd <- position_dodge(0.1) # move them .05 to the left and right

p3aa_min <-
    ggplot(plot_data, aes(x = x, y = auroc, group = version)) +
    geom_point(shape = 18, size = 2, aes(col = version)) +  # Use versions as points
  #geom_linerange(aes(ymin = auroc - 3*sd, ymax = auroc + 3*sd), position = position_dodge(width = 0.5), col = "black") +
  #geom_errorbarh(aes(xmin = auroc - 3*sd, xmax = auroc + 3*sd), height = 0.2, col = "black") +  #geom_errorbar(aes(ymin = auroc - sd, ymax = auroc + sd), width = 0.2, col = "black") +
    geom_errorbar(aes(ymin = auroc - 3*sd, ymax = auroc + 3*sd, col = version), width = 0.2, size = 0.4) + # position = position_dodge(width = 1)
    #geom_errorbarh(aes(xmin = auroc - 3*sd, xmax = auroc + 3*sd), height = 0.2) +
    xlab("Age Group") +
    ylab("AUROC") +
    scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
    scale_x_continuous(minor_breaks = NULL, breaks = 1:(length(age_labels)), labels = age_labels) +  # Set custom x-axis labels +  # Set custom x-axis labels
    # Set y-axis to start from 0
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.68, 0.79)) +
    guides(col = guide_legend(title = NULL)) +
    coord_fixed(ratio = 100) + # Set aspect ratio to 1:1
    theme_bw() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 1, size = 7),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank() # Remove vertical grid lines
    )

# Set up subpanel
plot_data$difference <- plot_data$auroc[plot_data$version == "V4"] -
  plot_data$auroc[plot_data$version == "V3"]

# Create the plot
p3ab <-
  ggplot(plot_data, aes(x = age_group, y = difference)) +
  geom_point(shape = 18, size = 3) +
  xlab("Age Group") +
  ylab("Difference in AUC") +
  #scale_color_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 500) +
  #coord_equal()+
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7)
  )

###

# create a plot that shows the AUC differences as a percentage increase
plot_data$percentage_increase <- ((plot_data$auroc[plot_data$version == "V4"] / plot_data$auroc[plot_data$version == "V3"]) - 1) * 100

# Create the plot
p3ac <-
  ggplot(plot_data, aes(x = age_group, y = percentage_increase)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  xlab("Age Group") +
  ylab("Percentage Increase") +
  scale_fill_manual(values = c("V4" = "black", "V3" = "red")) +
  scale_x_discrete(labels = age_labels) +
  coord_fixed(ratio = 1.5) + # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7)
  )

ggsave("Figures/pdfs/Fig_3a1.pdf", p3aa_min,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)

ggsave("Figures/pdfs/Fig_3a2.pdf", p3ab,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)

ggsave("Figures/pdfs/Fig_3a3.pdf", p3ac,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)


# past attempts
labs = c()
age_split = c(0, 20, 30, 40, 50, 60, 70, 80, 120)

for (i in 1:(length(age_split) - 1)) {
  labs = c(labs, paste0("(", age_split[i], ",", age_split[i+1], "]"))
}

xcol = c("black", "red")
xsc = 12
swidth = 0.1

m0 = dim(perf)[1]/2
n0 = dim(perf)[2]
perfmin = min(perf[1:m0, ]) - 0.02

df <- data.frame(
  x = rep((m0 + 2)*(-0.5 + 1:n0), each = m0),
  y = c(perf[1:m0, ] + perf[m0 + (1:m0), ]),
  group = factor(rep(1:n0, each = m0))
)

df2 <- data.frame(
  x = rep((m0 + 2)*(-0.5 + 1:n0), each = m0),
  y = c(perf[1:m0, ] - perf[m0 + (1:m0), ]),
  group = factor(rep(1:n0, each = m0))
)

p <- ggplot() +
  geom_line(data = df, aes(x = x, y = y, group = group), color = "black") +
  geom_line(data = df2, aes(x = x, y = y, group = group), color = "black") +
  geom_line(aes(x = (m0 + 2)*(-0.5 + 1:n0), y = perfmin + xrate/(xsc*max(xrate))), linetype = "dashed") +
  scale_x_continuous(limits = c(0, n0*(m0 + 2)),
                     breaks = (m0 + 2)*(1:n0) - floor(m0/2) - 1,
                     labels = labs) +
  scale_y_continuous(limits = c(perfmin, max(perf)),
                     breaks = perfmin + seq(0, 1/xsc, length = 4),
                     labels = signif(seq(0, max(xrate), length = 4), digits = 2)) +
  labs(x = "", y = "AUROC", title = "ROC and EA frequency by age band") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottomright") +
  guides(linetype = FALSE, color = FALSE, fill = FALSE) +
  annotate(geom = "text", x = n0*(m0 + 2) - floor(m0/2) - 1, y = perfmin,
           label = "EA frequency", hjust = 1, vjust = -1.5, size = 5)

p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# the difference and frequency in two different sub-figures below the main one,
# with the latter showing the AUC of each version in different colours

####
labs = c()
age_split = c(0, 20, 30, 40, 50, 60, 70, 80, 120)

for (i in 1:(length(age_split) - 1)) {
  labs = c(labs, paste0("(", age_split[i], ",", age_split[i+1], "]"))
}

xcol = c("black", "red")
xsc = 12
swidth = 0.1

m0 = dim(perf)[1]/2
n0 = dim(perf)[2]
perfmin = min(perf[1:m0, ]) - 0.02

df <- data.frame(
  x = rep((m0 + 2)*(-0.5 + 1:n0), each = m0),
  y = c(perf[1:m0, ] + perf[m0 + (1:m0), ]),
  group = factor(rep(1:n0, each = m0))
)

df2 <- data.frame(
  x = rep((m0 + 2)*(-0.5 + 1:n0), each = m0),
  y = c(perf[1:m0, ] - perf[m0 + (1:m0), ]),
  group = factor(rep(1:n0, each = m0))
)

diffs <- data.frame(
  x = (m0 + 2)*(-0.5 + 1:n0) - 0.5*(m0 + 2)/n0,
  y = round(xrate, 3),
  group = factor(rep(1, n0))
)

# attempt to reproduce the combined plot
p <- ggplot() +
  geom_line(data = df, aes(x = x, y = y, group = group), color = "black") +
  geom_line(data = df2, aes(x = x, y = y, group = group), color = "black") +
  geom_line(aes(x = (m0 + 2)*(-0.5 + 1:n0), y = perfmin + xrate/(xsc*max(xrate))), linetype = "dashed") +
  scale_x_continuous(limits = c(0, n0*(m0 + 2)),
                     breaks = (m0 + 2)*(1:n0) - floor(m0/2) - 1,
                     labels = labs) +
  scale_y_continuous(limits = c(0, max(perf)),
                     breaks = seq(0, max(perf), length = 4)) + # start y-axis at zero
  labs(x = "", y = "AUROC", title = "ROC and EA frequency by age band") +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottomright") +
  guides(linetype = FALSE, color = FALSE, fill = FALSE) +
  # annotate(geom = "text", x = n0*(m0 + 2) - floor(m0/2) - 1, y = perfmin,
  #          label = "EA frequency", hjust = 1, vjust = -0.5, size = 5) +
  geom_col(data = diffs, aes(x = x, y = y, group = group), fill = "red") + # add bottom sub-panel
  labs(subtitle = "Difference in AUROC")  # label the sub-panel

p + theme(axis.text.x = element_text(angle = 45, hjust = 1))



# plot 1: AUC plot
p1 <- ggplot() +
  geom_line(data = df, aes(x = x, y = y, group = group), color = "black") +
  geom_line(data = df2, aes(x = x, y = y, group = group), color = "black") +
  scale_x_continuous(limits = c(0, n0*(m0 + 2)),
                     breaks = (m0 + 2)*(1:n0) - floor(m0/2) - 1,
                     labels = labs) +
  scale_y_continuous(limits = c(0, max(perf)),
                     breaks = seq(0, max(perf), length = 4)) + # start y-axis at zero
  labs(x = "", y = "AUROC", title = "ROC and EA frequency by \nage band") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.position = "bottomright") +
  labs(subtitle = "Difference in AUROC")

# plot 2: frequency plot
p2 <- ggplot(diffs, aes(x = x, y = y, group = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  #geom_col(color = "black", fill = "red") +
  scale_x_continuous(limits = c(0, n0*(m0 + 2)),
                     breaks = (m0 + 2)*(1:n0) - floor(m0/2) - 1 ,
                     labels = labs) +
  scale_y_continuous(limits = c(0, max(diffs$y) + 0.01),
                     breaks = seq(0, max(diffs$y), length.out = 5),
                     labels = function(x) round(x, 2),
                     expand = expansion(mult = c(0.0000005, 0.00000005))) +
  labs(x = "Age Group", y = "EA Frequency") +
  coord_fixed(ratio = 100) +  # Set aspect ratio to 1:1
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7)
  )


ggsave("Figures/pdfs/Fig_3a4.pdf", p2,
       width = 8.5, height = 9, units = "cm",
       device = cairo_pdf)




# combine the plots vertically
library(gridExtra)
install.packages("egg")
library(egg)

figure_width <- 6
figure_height <- 6

library(patchwork)

# Enforce consistent figure size
p3aa_min <- p3aa_min + theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = figure_width / figure_height)
p2 <- p2 + theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = figure_width / figure_height)
p3ac <- p3ac + theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = figure_width / figure_height)

p3a <- grid.arrange(p3aa_min, p3ac, p2, ncol = 1)
plot_layout(p3a, ncol = 1, heights = c(10, 10, 10))

# Adjust the width and height as needed
#set_panel_size(width = 2, height = 2) # returns an error
#grid.arrange(p3aa_min, p3ac, p2, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1))
# unit(1, "npc"), heights = unit(1, "npc")) , widths = c(4, 4, 4),  heights = c(4, 4, 4)
ggsave("Figures/pdfs/Fig_3_comb.pdf", p3a,
       #width = 8.5, height = 9, units = "cm",
       width = 10.5, height = 10, units = "cm",
       device = cairo_pdf)

