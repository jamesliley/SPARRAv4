library(matrixStats)
library(lubridate)

source("Figures/util.R") # for import_sparra_expr()

plot_dir = "Analysis/time_attenuation/"


# Kernel density estimates normalised to unit diagonal
xcol=colorRampPalette(c("white","darkgray","red","yellow"))(100)
xx=1; nn=100; sp=7; scmax=100
qplot= ((1:1000)/1000)^(1/4) # Add q-q plots at these quantiles. >500 samples in each bin.


test_times=c(
  dmy_hm("1-5-2015 00:00"),
  dmy_hm("1-12-2015 00:00"),
  dmy_hm("1-5-2016 00:00"),
  dmy_hm("1-12-2016 00:00"),
  dmy_hm("1-5-2017 00:00")
)

for (i in 1:length(test_times)) {
  for (j in 1:length(test_times)) {
    if (i<j) {
      
      pi=c(0,0); pj=c(0,0)
      eval(import_sparra_expr(paste0(plot_dir, "cohort/crosstime/time",i,"_vs_time",j,".txt")))
      
      pdf(paste0("Figures/pdfs/Other/time",i,"_vs_time",j,".pdf"),width=7,height=7)
      
      
      xout=xout0
      
      # Normalise columns, then rows
      for (ii in 1:dim(xout)[2]) xout[ii,]=xout[ii,]/sum(xout[ii,],na.rm=T);
      for (ii in 1:dim(xout)[1]) xout[,ii]=xout[,ii]/sum(xout[,ii],na.rm=T);  
      xmax=min(colMaxs(xout[,10:(scmax-15)],na.rm=T));
      xmin=min(xout[which(xout>0)],na.rm=T); ncol0=5
      
      image((1:scmax)/100,(1:scmax)/100,xout,xlab=paste0("Prediction for time cutoff ",i,", ",as.character(test_times[i])),
            ylab=paste0("Prediction for time cutoff ",j,", ",as.character(test_times[j])),
            col=xcol,breaks=c(0,rep(xmin/2,ncol0),seq(xmin,4*xmax/3,length=length(xcol)-ncol0-1),1+max(xout,na.rm=T)),
            xaxs="i",yaxs="i",xlim=c(-1/sp,scmax/100),ylim=c(-1/sp,scmax/100),bty="n")
      
      polygon(d1$x,d1$y/(sp*max(d1$y))-(1/sp),col="gray"); 
      polygon(d2$y/(sp*max(d2$y))-(1/sp),d2$x,col="gray")
      
      abline(0,1,col="black",lty=2,lwd=2)
      #abline(0,1,col="white",lty=1,lwd=2)
      
      #lines(q0,q1,col="black",lwd=4)
      lines(q0,q1,col="blue",lwd=2)
      
      #lines(q0a,q1a,col="black",lwd=4)
      lines(q0a,q1a,col="darkgreen",lwd=2)
      
      legend("bottomright",c("Density (low)","", "Density (high)", "Marginal", "X-Y","Q-Q","Median (col)","Median (row)"),
             pch=c(16,16,16,NA,NA,NA,NA,NA),lty=c(NA,NA,NA,1,2,1,1,2),lwd=c(NA,NA,NA,1,2,2,2,2),
             col=c("darkgray","red","yellow","black","black","blue","darkgreen","darkgreen"))
      dev.off()
    }
  }
}

