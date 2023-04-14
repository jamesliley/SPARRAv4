# Load the required libraries
library("tidyverse")
library("patchwork")
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
perf=c(); xrate=c()
for (i in 1:(length(age_split)-1)) {
  sub=which(all_pred$age>age_split[i] & all_pred$age <= age_split[i+1] & is.finite(all_pred$super+all_pred$v3))
  psp=getroc(all_pred$target[sub],all_pred$super[sub])
  p3=getroc(all_pred$target[sub],all_pred$v3[sub])
  perf=cbind(perf,c(psp$auc,p3$auc,psp$se,p3$se))
  xrate=c(xrate,sum(all_pred$target[sub])/length(sub))
}

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

# Define the age split values
age_split <- c(0, 20, 30, 40, 50, 60, 70, 80, 120)

# Initialize empty vectors to store performance values
perf <- c()
xrate <- c()

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
labs <- c()
for (i in 1:(length(age_split)-1)) {
  labs <- c(labs, paste0("(", age_split[i], ",", age_split[i+1], "]"))
}

# Plot the data using ggplot2
ggplot(df, aes(x = 1:nrow(df))) +
  geom_segment(aes(y = perf[,1], yend = perf[,2], color = "V4"), size = 1) +
  geom_segment(aes(y = perf[,1] - perf[,3], yend = perf[,1] + perf[,3], color = "V4"), size = 1) +
  geom_segment(aes(y = perf[,2] - perf[,4], yend = perf[,2] + perf[,4], color = "V4"), size = 1) +
  geom_segment(aes(y = xrate / max(xrate), yend = xrate / max(xrate), color = "Freq."), linetype = "dashed", size = 1) +
  scale_x_continuous(breaks = 1:nrow(df), labels = labs) +
  scale_y_continuous(limits = c(min(perf) - 0.02, max(perf)), expand = c(0, 0)) +
  scale_color_manual(values = c("V4" = "black", "Freq." = "black", "V3" = "red")) +
  labs(x = "", y = "AUROC", title = "ROC and EA frequency by age band") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###

