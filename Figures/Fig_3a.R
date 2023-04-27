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
par(mar=c(5.1,4.1,4.1,4.1))

xcol=c("black","red"); xsc=12; swidth=0.1
m0=dim(perf)[1]/2; n0=dim(perf)[2]
#perfmin=min(perf[1:m0,])-0.02
perfmin = 0
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

# CAV: The next line fails because of "labs". However, this can be fixed
# by manually generating values for it (you can take it from the paper)
# Note that it would be best to change the name to e.g. "age.labs" as
# labs() is a function within ggplot2
axis(1,at=(m0+2)*(1:n0)- floor(m0/2)-1,labels=labs,cex.axis=1,las=2)
axis(4,at=perfmin+seq(0,1/xsc,length=4),labels=signif(seq(0,max(xrate),length=4),digits=2))
mtext("EA frequency",side=4,line=3)

legend("bottomright",lty=c(1,1,2),col=c(xcol,"black"),c("V4","V3","Freq."),bg="white")

dev.off()

# CAV: Louis and I discussed that the y-axis can be a bit misleading
# as it visually suggest a large difference, when it's not if we consider
# the scale starting from zero. One potential solution would be to start
# the y-axis at zero, and then generate a bottom sub-panel that shows the
# difference in AUROC (similar to what was done in Figure 2).
# Perhaps you could work on that?

# new code

# Define the age split values
age_split <- c(0, 20, 30, 40, 50, 60, 70, 80, 120)

# Loop through age split values
for (i in 1:(length(age_split)-1)) {
  # Subset data based on age split values
  sub <- which(all_pred$age > age_split[i] & all_pred$age <= age_split[i+1] & is.finite(all_pred$super+all_pred$v3))

  # Compute ROC metrics
  psp <- getroc(all_pred$target[sub], all_pred$super[sub])
  p3 <- getroc(all_pred$target[sub], all_pred$v3[sub])

  # Store performance values
  perf <- cbind(perf, c(psp$auc, p3$auc, psp$se, p3$se))

  # Compute EA frequency
  xrate <- c(xrate, sum(all_pred$target[sub]) / length(sub))
}

# Create a data frame for plotting
df <- data.frame(perf)
df$xrate <- xrate

# Define labels for x-axis
age.labs <- c()
for (i in 1:(length(age_split)-1)) {
  age.labs <- c(age.labs, paste0("(", age_split[i], ",", age_split[i+1], "]"))
}

# Plot the data using ggplot2
ggplot(df, aes(x = 1:nrow(df))) +
  geom_segment(aes(y = perf[,1], yend = perf[,2], color = "V4"), size = 1) +
  geom_segment(aes(y = perf[,1] - perf[,3], yend = perf[,1] + perf[,3], color = "V4"), size = 1) +
  geom_segment(aes(y = perf[,2] - perf[,4], yend = perf[,2] + perf[,4], color = "V4"), size = 1) +
  geom_segment(aes(y = xrate / max(xrate), yend = xrate / max(xrate), color = "Freq."), linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = 1:nrow(df), labels = age.labs) +
  scale_y_continuous(limits = c(min(perf) - 0.02, max(perf)), expand = c(0, 0)) +
  scale_color_manual(values = c("V4" = "black", "Freq." = "black", "V3" = "red")) +
  labs(x = "", y = "AUROC", title = "ROC and EA frequency by age band") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

