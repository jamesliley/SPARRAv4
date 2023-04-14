library("tidyverse")
library("patchwork")

# Load the ggplot2 library
library(ggplot2)

setwd("~/SPARRAv4")
source("Figures/util.R") # for import_sparra_expr()
source("SPARRAv4/auxiliary.R") # for roc_2panel()

# conversion of this figure is pending
eval(import_sparra_expr("Analysis/full_model/Shapley_values/age_by_simd_effect.txt"))


# conversion of figure 3a performance_by_age
plot_dir="Analysis/full_model/"
eval(import_sparra_expr(paste0(plot_dir, "Analytics/performance_by_age.txt")))

# old code
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

# new code

# Create data for the plot
age_split <- c(0, 20, 30, 40, 50, 60, 70, 80, 120)
xcol <- c("black", "red")
xsc <- 12
swidth <- 0.1
m0 <- dim(perf)[1]/2
n0 <- dim(perf)[2]
perfmin <- min(perf[1:m0,]) - 0.02
labs <- c()
for (i in 1:(length(age_split)-1)) {
  labs <- c(labs, paste0("(", age_split[i], ",", age_split[i+1], "]"))
}
data <- data.frame(x = rep((m0+2)*(-0.5 + 1:n0), each = m0),
                   y = c(perf[1:m0,]) + xrate/(xsc*max(xrate)),
                   col = rep(xcol, each = m0))

# Create the plot using ggplot2
ggplot(data, aes(x = x, y = y, col = col)) +
  geom_line() +
  geom_segment(aes(xend = x, yend = y - perf[m0 + (1:m0),], col = col),
               linetype = "dashed", size = 1) +
  geom_segment(aes(xend = x, yend = y + perf[m0 + (1:m0),], col = col),
               linetype = "dashed", size = 1) +
  scale_color_manual(values = xcol) +
  labs(title = "ROC and EA Frequency by Age Band", x = "X Axis Label", y = "AUROC") +
  theme_minimal() +
  theme(legend.position = "bottomright") +
  scale_x_continuous(breaks = (m0+2)*(1:n0) - floor(m0/2) - 1, labels = labs) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "EA Frequency",
                                         breaks = perfmin + seq(0, 1/xsc, length = 4),
                                         labels = signif(seq(0, max(xrate), length = 4), digits = 2))) +
  guides(col = guide_legend(title = "Legend Title"))



###
# New Figure 'age by simd effect'

df <- rbind(
  cbind(Model = "v3", data.frame(sens = xroc3b$sens[1,],
                                 spec = xroc3b$spec[1,])),
  cbind(Model = "v4", data.frame(sens = xroc4b$sens[1,],
                                 spec = xroc4b$spec[1,])),
  cbind(Model = "Max", data.frame(sens = xrocmb$sens[1,],
                                  spec = xrocmb$spec[1,]))
) |> mutate(Model = fct_relevel(Model, "v3", "v4", "Max"))
df2 <- build_diff(df, spec, sens)

p1 <- ggplot(df |>
               mutate(spec = 1-spec) |>
               filter(Model != "Max")) +
  geom_line(aes(x = spec, y = sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("") + ylab("Sensitivity") +
  theme_minimal(base_size = 8) + theme(legend.justification = c(1,0),
                                       legend.position = c(1,0),
                                       legend.spacing = unit(0, "npc"),
                                       legend.margin = unit(0, "npc"),
                                       legend.background = element_rect(fill = "white", size = 0, colour = "white"))

p2 <- ggplot(df2 |>
               mutate(spec = 1-spec,
                      v4 = v4 - v3,
                      Max = Max - v3,
                      v3 = v3 - v3) |>
               pivot_longer(v3:Max, names_to = "Model", values_to = "delta_sens") |>
               mutate(Model = fct_relevel(Model, "v3", "v4", "Max")) |>
               filter(Model != "Max")) +
  geom_line(aes(x = spec, y = delta_sens, col = Model), linewidth = 0.4) +
  xlim(0, 1) +
  xlab("Recall") + ylab("Δ Sensitivity") +
  theme_minimal(base_size = 8) + theme(legend.position = "none")

p <- p1 / p2 + plot_layout(heights = c(3,1))

print(p)

ggsave("Figures/pdfs/Fig_2a.pdf",
       width = 7.5, height = 7.25, units = "cm",
       device = cairo_pdf)




# Old Figure 'age by simd effect'

plot(0,type="n",xlim=range(x0),ylim=yvar,main=longvarnames(colnames(shapley)[ii]),xaxt="n",
     xlab="",ylab="Shapley value") # note fixed ylim
lines(x0,mx,lty=2,col="red")
segments(x0,mx-sdx,x0,mx+sdx,col="red")
points(x0,mx,pch=16,col="red")
axis(1,at=x0,labels=lab,las=2)
abline(h=0,lty=2)
dev.off()

