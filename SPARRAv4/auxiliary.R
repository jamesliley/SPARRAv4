# Various useful auxiliary functions

##' Function for evaluating hyperparameter performance.
##' @name hyp_evaluate
##' @param method function to fit model from h2o
##' @param parameters set of parameters to evaluate, named according to arguments of h2o.deeplearning
##' @param train training matrix
##' @param test test matrix
##' @param metric metric to evaluate performance by
##' @param ... additional parameters passed to h2o
hyp_evaluate=function(method,parameters,train,test,metric,...) {
	modfit=do.call(method,modifyList(list(training_frame=train,...),parameters))
	#  print(modfit@parameters)
	print(parameters)
	#h2o.rm(train)
	fitx=h2o.performance(modfit, test)
	#h2o.rm(test)
	fitx@metrics[metric]
}



##' plotcal() 
##' Produces a set of points for a calibration plot, and optionally plots them.
##' Uses either a binning method or a kernel method to determine height of points. 
##' In both method, divides unit interval into subintervals [0,1/n], [1/n,2/n), ... [(n-1)/n,1). 
##' For bin \[a,b\)
##'    x co-ordinate is (a+b)/2
##'   For binning method
##'    y co_ordinate is mean({y:predicted value for y is in [a,b)})
##'   For kernel method
##'    y co-ordinate is weighted mean of all y values with the weight of value yi given by dnorm(y-yi,sd=kernel_sd)
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param n number of subintervals/points
##' @param kernel set to TRUE to use kernel method
##' @param kernel_sd kernel width for kernel method; see above
##' @param conf include a confidence interval; alpha, c0, c2 are only relevant if conf=TRUE
##' @param alpha return a pointwise confidence envolope for conservative 1-alpha confidence interval
##' @param c0 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @param c2 for computing maximum bias; assume true covariance function is of the form a0+ a1x + a2x^2, with |a0|<c0, |a2|<c2 (c1 does not matter)
##' @param plot set to FALSE to suppress plot
##' @param ... further parameters passed to plot()
##' @return n x 2 matrix containing co-ordinates of points on the curve.
plotcal=function(y,ypred,n=10,kernel=F,kernel_sd=0.05,alpha=0.05,c0=0,c2=0.1,plot=TRUE,conf=TRUE,...)  {
  if (!kernel) {
   ycal=rep(0,n); xcal=ycal
   xup=rep(0,n); xdown=rep(0,n)
	 for (i in 1:n) {
		  sub=which((ypred> (i-1)/n) & (ypred< i/n))
 	  	ycal[i]=mean(y[sub])
 	  	xcal[i]=mean(ypred[sub])
 	  	if (conf) {
 	  	  xse=sqrt(ycal[i]*(1-ycal[i])/length(sub))
 	  	  xup[i]=ycal[i] -qnorm(alpha/2)*xse
 	  	  xdown[i]=ycal[i] + qnorm(alpha/2)*xse
 	  	}
	 }
	 if (plot) plot(xcal,ycal,...)
	 return(list(x=xcal,y=ycal,upper=xup,lower=xdown))
  } else {  # use kernel method with given sd
   ypredsub=seq(0,1,length=n)
   kern=function(x,y) dnorm(x-y,sd=kernel_sd)
   wt=outer(ypredsub,ypred,kern)
   x1=(wt %*% ypred)
   y1=(wt %*% y)
   csub=ypredsub*y1/x1
   csub=pmax(pmin(csub,1),0)

   if (plot) plot(ypredsub,csub,...)
   
   if (conf) {
   # Confidence intervals
   wts=ypredsub*wt/as.vector(x1) # now csub= wts %*% y
   yvar=(wts^2 %*% as.vector(ypred*(1-ypred)))
   
   # Max bias
   bias=rowSums(outer(ypredsub,ypred,function(x,y) kern(x,y)*(c0*(y-x) + c2*(x^2 * y - y^2 * x))))/x1
   
   upper=csub - qnorm(alpha/2)*sqrt(yvar) + bias
   lower=csub + qnorm(alpha/2)*sqrt(yvar) - bias
   
   return(list(x=ypredsub,y=csub,lower=lower,upper=upper))
   } else {
   return(cbind(ypredsub,csub))
   }
  }
}

##' fastroc() 
##' Quick plotting function for receiver-operator characteristic curve.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param res resolution; computes sensitivity and specificity for this many uniformly-spaced points in the range of ypred
##' @param cuts cutoffs, between 0 and 1 (scaled). If NULL use res to generate uniform cuts
##' @param plot set to FALSE to suppress plot
##' @param ... further parameters passed to plot()
##' @return n x 2 matrix containing co-ordinates of points on the curve.
fastroc=function(y,ypred,res=100,cuts=NULL,plot=T,...) {
if (is.null(cuts)) cuts=seq(min(ypred),max(ypred),length.out=res) else cuts=min(ypred) + cuts*(max(ypred)-min(ypred))
sens=0*cuts; spec=0*cuts
yt=sum(y)
for (i in 1:res) {
	sens[i]=sum(y* (ypred >= cuts[i]) )/yt
	spec[i]=sum((1-y)*(ypred<cuts[i]))/(length(y)-yt)
}
if (plot) plot(1-spec,sens,...)
return(cbind(sens,spec))
}

##' getroc() 
##' Comprehensive plotting function for receiver-operator characteristic curve. Also calculates AUROC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' SE of AUROC with no CV structure is from Hanley and McNeil 1982. SE of AUROC with CV folds is from LeDell et al 2012
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @return list containing: spec, specificity for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
getroc=function(y,ypred,cv=NULL,res=100,addauc=FALSE,cols=NULL) {
if (is.null(cv)) cv=rep(1,length(y))
if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")

sens=c(); spec=c(); auc=c(); se=c()
for (i in 1:max(cv)) {
y0=y[which(cv==i)]; 
ypred0=ypred[which(cv==i)]

yt=sum(y0); yl=length(y0)
opred=order(ypred0)
#ipred=order(opred) # can use ipred to reorder in the order of original ypred

sy=y0[opred]; sp=ypred0[opred]

sens0=1- (cumsum(sy)/yt)
spec0=cumsum(1-sy)/(yl-yt) 
auc0=integral(sens0,spec0)
se0=aucse(as.numeric(yt),as.numeric(yl-yt),auc0)

if (!is.null(res)) {
  ds=cumsum(sqrt((spec0[1:(yl-1)]-spec0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
  ds=ds/ds[yl-1]
  lsp=(1:(yl-1))/yl
  sub=round(yl*approx(ds,lsp,n=res)$y)
  sens0=sens0[sub]
  spec0=spec0[sub]
}

auc=c(auc,auc0)
se=c(se,se0)
spec=rbind(spec,spec0)
sens=rbind(sens,sens0)
}

if (length(auc)>1) {
  auc=c(auc,mean(auc))
  se=c(se,ci.cvAUC(ypred,y,folds=cv)$se)
}

out=list(sens=sens,spec=spec,auc=auc,se=se)
class(out)="sparraROC"
return(out)
}

# Internal function to compute SE of AUC
aucse=function(n1,n2,auc) {
  q1=auc/(2-auc); q2=2*(auc^2)/(1+auc)
  num=auc*(1-auc) + (n1-1)*(q1- (auc^2)) + (n2-1)*(q2-(auc^2))
  return(sqrt(num/(n1*n2)))
}

##' Plot function for class above
##' @param out output from getroc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
plot.sparraROC=function(out,addauc=FALSE,cols=rep("black",dim(out$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="1-Spec.",ylab="Sens.",...)
  ncv=dim(out$spec)[1]
  for (i in 1:ncv) lines(1-out$spec[i,],out$sens[i,],col=cols[i])
  abline(0,1,col="red",lty=2)
  auc=out$auc[length(out$auc)]
  se=out$se[length(out$se)]
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}


##' getprc() 
##' Comprehensive plotting function for precision-recall curve. Also calculates AUPRC and standard error. 
##' 
##' Rather than returning points corresponding to every cutoff, only returns a representative sample of equally-spaced points along the curve.
##'
##' Does not plot anything. Object can be plotted in a default way.
##'
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param cv cross-validation fold assignments, if relevant. Changes estimate of standard error.
##' @param res resolution. Returns this many equally-spaced points along the curve. Set res to null to return all points.
##' @return list containing: ppv, ppv for res points in every cv fold; sens, sensitivity for res points in every cv fold; auc, areas under the curve for each fold and average (note length is 1 greater than number of CV folds); se, standard error for AUC in each fold and standard error for average auc (note length is 1 greater than number of CV folds)
getprc=function(y,ypred,cv=NULL,res=100,addauc=FALSE,cols=NULL) {
  if (is.null(cv)) cv=rep(1,length(y))
  if (!(length(y)==length(ypred))) stop("Parameters y and ypred should have the same length")
  
  sens=c(); ppv=c(); auc=c(); se=c()
  for (i in 1:max(cv)) {
    y0=y[which(cv==i)]; 
    ypred0=ypred[which(cv==i)]
    
    yt=sum(y0); yl=length(y0)
    opred=order(ypred0)
    #ipred=order(opred) # can use ipred to reorder in the order of original ypred
    sy=y0[opred]; sp=ypred0[opred]

    ppv0=rev(cumsum(rev(sy))/(1:length(sy)))
    sens0=1- (cumsum(sy)/yt)

    auc0=integral(sens0,ppv0)
    
    
    se0=sqrt(auc0*(1-auc0)/sum(y0))

    if (!is.null(res)) {
      ds=cumsum(sqrt((ppv0[1:(yl-1)]-ppv0[2:yl])^2 + (sens0[1:(yl-1)]-sens0[2:yl])^2))
      ds=ds/ds[yl-1]
      lsp=(1:(yl-1))/yl
      sub=suppressWarnings(round(yl*approx(ds,lsp,n=res)$y))
      sens0=sens0[sub]
      ppv0=ppv0[sub]
    }
    
    auc=c(auc,auc0)
    se=c(se,se0)
    ppv=rbind(ppv,ppv0)
    sens=rbind(sens,sens0)
  }

  if (length(auc)>1) {
    auc=c(auc,mean(auc))
    se=c(se,mean(se)/sqrt(3))
  }
  
    
  out=list(sens=sens,ppv=ppv,auc=auc,se=se)
  class(out)="sparraPRC"
  return(out)
}


##' Plot function for class above
##' @param out output from getprc()
##' @param addauc set to TRUE to add text to the plot showing the (mean) AUC and SE.
##' @param cols colour to draw lines
plot.sparraPRC=function(out,addauc=FALSE,cols=rep("black",dim(out$sens)[1]),...) {
  plot(0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",...)
  ncv=dim(out$sens)[1]
  for (i in 1:ncv) lines(out$sens[i,],out$ppv[i,],col=cols[i])
  auc=mean(out$auc)
  se=mean(out$se)/sqrt(3)
  txx=paste0(signif(auc,digits=2),"+/-",signif(se,digits=2))
  if (addauc) text(0.6,0.4,txx)
}









##' integral() 
##' Quick form for trapezoidal integration over range of x
##'
##' @param x x co-ordinates, or nx2 matrix of points 
##' @param y y co-ordinates
##' @return trapezoidal estimate of integral of y[x] over range of x.
integral=function(x,y=NULL) {
	if (is.null(y)) {
		y=x[,2]; x=x[,1]
	}
	ox=order(x); xs=x[ox]; ys=y[ox]
	sum((xs[-1]-xs[-length(xs)])*(ys[-1]+ys[-length(ys)]))/2
}

##' ab() 
##' Shorthand to draw a red x-y line
ab=function(...) abline(0,1,col="red",...)

##' logit() 
##' Logistic function
##' @param x parameter
logit=function(x) 1/(1+exp(-x))

##' ilogit() 
##' Inverse ogistic function
##' @param x parameter between 0 and 1
ilogit=function(x) -log((1/x)-1)





##' roc_2panel
##' Draws a ROC curve (with legend) with a second panel underneath showing sensitivity difference.
##' 
##' @param rocs list of sparraROC objects (if one object, plots folds separately)
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to red
##' @param ... other parameters passed to legend()
roc_2panel=function(rocs,labels,col,lty=rep(1,length(col)),xy_lty=2,xy_col="red",...) {
  
  if ("sens" %in% names(rocs)) {
    nfold=length(rocs$auc)-1; np=length(rocs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=rocs$sens[i,],
         spec=rocs$spec[i,])
    rocs=r0
  }
  
  # Set up plot parameters
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Sensitivity",type="n")
  abline(0,1,col=xy_col,lty=xy_lty)
  
  # x-values to draw sensitivity at
  xspec=seq(0,1,length=100)[2:99]; xr=c()
  
  # Draw ROC curves on top panel
  for (i in 1:length(rocs)) {
    xy=rocs[[i]]
    lines(1-xy$spec,xy$sens,col=col[i],lty=lty[i])
    xsens=suppressWarnings(approx(1-xy$spec,xy$sens,xspec)$y)
    xr=rbind(xr,xsens)
  }
  
  # Add legend
  legend("bottomright",legend=labels,col=col,lty=lty,...)
  
  # Bottom panel setup
  par(mar=c(4,4,0.1,0.1))
  plot(0,xlim=c(0,1),ylim=range(t(xr)-xr[1,]),type="n",
    xlab="1-Specificity",ylab=expression(paste(Delta,"(sens.)")),yaxt="n")
  axis(2,at=pretty(range(t(xr)-xr[1,]),n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(rocs)) lines(xspec,xr[i,]-xr[1,],col=col[i],lty=lty[i])
  #abline(h=0,col=xy_col,lty=xy_lty)
  
}  




##' prc_2panel
##' Draws a PRC curve (with legend) with a second panel underneath showing precision difference.
##' 
##' @param prcs list of sparraPRC objects. If of length 1, splits into folds
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param ... other parameters passed to legend()
prc_2panel=function(prcs,labels,col,lty=rep(1,length(col)),...) {
  
  if ("ppv" %in% names(prcs)) {
    nfold=length(prcs$auc)-1; np=length(prcs$sens)/nfold
    r0=list()
    for (i in 1:nfold) 
      r0[[i]]=list(sens=prcs$sens[i,],
        ppv=prcs$ppv[i,])
    prcs=r0
  }
  
  
  
  # Set up plot 
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Precision",type="n")
  xr=c(); 
  
  # X-values at which difference in precision will be plotted
  xsens=seq(0,1,length=100)[2:99]
  
  # Draw curves on top panel
  for (i in 1:length(prcs)) {
    px=prcs[[i]]
    lines(px$sens,px$ppv,col=col[i],lty=lty[i])
    xppv=suppressWarnings(approx(px$sens,px$ppv,xsens)$y)
    xr=rbind(xr,xppv)
  }
  
  # Add legend
  legend("topright",legend=labels,col=col,lty=lty,...)
  
  # Initialise bottom panel
  par(mar=c(4,4,0.1,0.1))
  plot(0,xlim=c(0,1),ylim=range(t(xr)-xr[1,]),type="n",
    xlab="Recall",ylab=expression(paste(Delta,"(Prec.)")),
    yaxt="n")
  axis(2,at=pretty(range(t(xr)-xr[1,]),n=2))
  
  # Draw lines on bottom panel
  for (i in 1:length(prcs)) lines(xsens,xr[i,]-xr[1,],col=col[i],lty=lty[i])
}



##' cal_2panel
##' Draws calibration curves (with legend) with a second panel underneath showing predicted differences.
##' 
##' @param cals list of calibration objects, output from plotcal(). 
##' @param labels labels to use in legend
##' @param col line colours
##' @param lty line type, defaults to 1
##' @param xy_lty line type for x-y line, defaults to 2 (dashed)
##' @param xy_col line colour for x-y line, defaults to red
##' @param ci_col colours to draw confidence intervals on lower panel; NA to not draw. 
##' @param 
##' @param ... other parameters passed to legend()
cal_2panel=function(cals,labels,col,lty=rep(1,length(col)),xy_lty=2,xy_col="red",ci_col=rep(NA,length(col)),...) {

    
  # Setup plot
  par(mar=c(1,4,0.1,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))
  
  # Initialise plot
  plot(0,xlim=c(0,1),ylim=c(0,1),xaxt="n",ylab="Observed",type="n")
  
  # Draw lines on top panel
  xr=c()
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    lines(cx$x,cx$y,col=col[i],lty=lty[i])
    xr=rbind(xr,cx$y-cx$x)
  }
  abline(0,1,col=xy_col,lty=xy_lty,lwd=1)
  
  # Draw legend
  legend("topleft",legend=labels,col=col,lty=lty,...)
  
  
  # Initialise bottom panel
  par(mar=c(4,4,0.1,0.1))
  xpred=seq(0,1,length=100)[2:99]
  plot(0,xlim=c(0,1),ylim=range(xr),type="n",
    xlab="Predicted",ylab=expression(paste(Delta,"(cal.)")),
    yaxt="n")
  axis(2,at=pretty(range(xr),n=3))
  
  # Draw lines on bottom panel
  for (i in 1:length(cals)) {
    cx=cals[[i]]
    if (!is.na(ci_col[i]))
      polygon(c(cx$x,rev(cx$x)),c(cx$lower,rev(cx$upper))-c(cx$x,rev(cx$x)),
        col=ci_col[i],border=NA)
    lines(cx$x,cx$y-cx$x,col=col[i],lty=lty[i])
  }
  abline(h=0,col=xy_col,lty=xy_lty)
  
}  






##' forcecal() 
##' Attempt to find a transformation t(ypred) which is better-calibrated than ypred.
##' 
##' Denote C(p) as the (kernel) calibration of the model at {estimated P(Y=1)}=p,  C:[0,1] -> [0,1]. This is the function computed by plotcal(...,kernel=T).
##' 
##' We assess overall mis-calibration as the integral of C(p)-p over 0<p<1. We seek to find a function t(p) which minimises
##' 
##' F{t}=integral of C(t(p))-p over 0<p<1
##' 
##' We consider basis functions f1,f2,.. fn and model
##'  
##' t(p)=a + b1 f1(p) + b2 f2(p) ... + bn fn(p)
##'
##' maximising over parameters a,b1,..bn.
##' 
##' This function does NOT automatically enforce monotonicity of t, but this can be enforced by only using basis functions with nonnegative derivative and setting parameter poscoef to TRUE
##' 
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param nspace estimate integral using this many points; this is passed as 'n' to plotcal
##' @param kernel_sd passed to plotcal
##' @param flist list of basis functions. Defaults to monomials of increasing degree.
##' @param par_start values from which to start search for parameters b1, b2... bn. Does NOT include starting value for constant term a, which starts at 0.
##' @param poscoef set to TRUE to force positive values for a,b1,b2,... bn
##' @param coef_max maximum absolute value of coefficients a, b1,b2...bn
##' @param ... further parameters passed to R's 'optim' function. Uses method L-BFGS-B. Use control=list(trace=6) to track progress.
##' @return vector of optimal coefficients.
forcecal_approximate=function(y,ypred,nspace=20,kernel_sd=0.05,
  flist=list({function(x) x},{function(x) x^2},{function(x) x^3},{function(x) x^4}),
  par_start=rep(1,length(flist)),poscoef=FALSE,coef_max=10,...) {

# objective function (of parameters) to be maximised  
obj=function(par) {
  a=par[1]; b=par[2:length(par)]
  
  # transform function
  transf=function(x) pmax(pmin(a+ t(do.call(rbind,lapply(flist,function(f) f(x)))) %*% b,1),0)
  
  # calibration assessment (kernel)
  transvec=as.vector(transf(ypred)); transvec=pmax(transvec,0); transvec=pmin(transvec,1)
  out=plotcal(y,transvec,plot=F,n=nspace,kernel_sd=kernel_sd,kernel=TRUE,conf=FALSE)
  
  xout=sum((out[,2]-out[,1])^2)
  
  if (is.na(xout)) xout=nspace
  # integral between curve and x-y line, approx.
  return(xout)
}

if (poscoef) lst=c(-coef_max,rep(0,length(par_start))) else lst=c(-coef_max,rep(0,length(par_start)))

fit=optim(c(0,par_start),obj,
  lower=lst,upper=rep(coef_max,1+length(par_start)),
  method="L-BFGS-B",...)

out=fit
f0=function(x) 
  as.vector(
    pmax(pmin(
      fit$par[1] + 
        t(do.call(rbind,lapply(flist,function(f) f(x)))) %*%
        fit$par[2:length(fit$par)],
      1),0))
f1=function(x) (1- (1e-5))*f0(x) + (1e-5)*x # ensure function is strictly increasing
out[[length(fit)+1]]=f1
names(out)=c(names(fit),"transform")

return(out)
}




##' forcecal() 
##' Finds a function t(ypred) which is better-calibrated than ypred.
##' 
##' Unlike forcecal_approximate, does this by simply inverting the empirical calibration function.
##' 
##' Denote C(p) as the (kernel) calibration of the model at {estimated P(Y=1)}=p,  C:[0,1] -> [0,1]. This is the function computed by plotcal(...,kernel=T).
##' 
##' Unlike forcecal_general, this function does not include a constant factor and automatically enforces monotonicity of t, by constraining its derivative to be positive.
##' 
##' @param y class labels, 0/1 or logical
##' @param ypred predictions Pr(Y=1), numeric vector
##' @param nspace estimate integral using this many points; this is passed as 'n' to plotcal
##' @return function t
forcecal=function(y,ypred,nspace=50,...) {

  ab=mmt(ypred,y)
  a=ab[1]; b=ab[2]
  ypred=a*ypred + b
  
  # Compute empirical calibration function (this needs to be non-kernel based)
  px=plotcal(y,ypred,n=nspace,kernel=F,plot=F)
  xc=px$x; yc=px$y

  # Remove NA
  w=which(!is.na(xc+yc))
  xc=xc[w]; yc=yc[w]
  
  # Ensure px$y is increasing and include 0,1 values
  yc=cummax(yc)
  y_e=c(0,yc,1)
  x_e=c(0,(xc-b)/a,1)
  
  y_e=y_e*(1-(1e-5)) + x_e*(1e-5)
  
  # Calibrating transform
  transf=function(x) suppressWarnings(approx(x_e,y_e,x)$y)
  
  out=list(x_e=x_e,y_e=y_e,transform=restrictEnvironment(transf,c("x_e","y_e")))
  
  return(out)
}


##' General utility function to restrict environment of function
##' @param f function to restrict
##' @param vars variables to put in environment
restrictEnvironment=function(f,vars) {
  oldenv=environment(f)
  newenv=new.env(parent=parent.env(oldenv))
  for (v in vars) {
    assign(v,get(v,envir=oldenv),envir=newenv)
  }
  environment(f)=newenv
  return(f)
}



##' Auxiliary function for forcecal
##' Given x, y finds a transform which correctly calibrates mean(x) and mode(x), and 
##' ensures 0 <= x <= 1
##' @param x predictor
##' @param y target
##' @return pair of coefficients a,b where transform is x-> ax + b
mmt=function(x,y) {
  med1=mean(x,na.rm=T)
  t1=table(round(100*x))
  mode1=as.numeric(names(t1)[which.max(t1)])/100
  
  wmed=which(abs(x-med1)<0.01)
  pmed=mean(x[wmed])
  tmed=mean(y[wmed])
  
  wmode=which(abs(x-mode1)<0.01)
  pmode=mean(x[wmode])
  tmode=mean(y[wmode])
  
  a=(tmode-tmed)/(pmode-pmed)
  b=tmode - a*pmode
  
  xt=a*x + b
  
  if (min(xt,na.rm=T)<0) trs=min(xt,na.rm=T) else trs=0
  if (max(xt,na.rm=T)>1) sc=max(xt,na.rm=T) else sc=1
  
  
  return(c(a/sc, (b-trs)/sc))
}








##' turn_data_into_sparrav3_factors() 
##' Turn the SPARRAv3 columns into factor levels (as used by SPARRAv3). Written by  Gergo. 
##' @param df input matrix of type used in SPARRAv4
##' @return matrix for SPARRAv3 computation
# Redacted for privacy -  not needed for pipeline



##' combine_training_matrices()
##' 
##' Function to tidy up and combine training matrices
##' @param v3_matrix output from eg transformer(...,list(sparrav3 = list())). This can be either a matrix or a filename.
##' @param topic_matrix output from eg transformer(...,list(all_codes_topic = list(topic_model_fit = topic_model_fit))). This can be either a matrix or a filename.
##' @param time_matrix output from eg transformer(...,list(last_episode_days_ago = list(source_table_names = c("AE2", "SMR00", "SMR01", "SMR01E", "SMR04")),last_emergency_admission_days_ago = list(),ltcs = list(output_type = "years_since_diag")))). This can be either a matrix or a filename.
##' @param ltc_matrix output from eg transformer(...,list(ltcs = list(output_type = "rawdata_NUMBEROFLTCs"))). This can be either a matrix or a filename.
##' @param patients,episodes,list_of_data_tables from loadCleanData. Reloaded if null.
##' @param keep_id_and_time set to TRUE to retain ID and time columns
##' @returns a training matrix with the above matrices left-joined
combine_training_matrices=function(v3_matrix,topic_matrix,time_matrix,ltc_matrix,
  patients=NULL,episodes=NULL,list_of_data_tables=NULL,keep_id_and_time=TRUE) {

  ##### Begin with v3 matrix
  if (is.character(v3_matrix)) train_matrix=readRDS(v3_matrix) else train_matrix=v3_matrix
  print("Loaded v3 matrix")
  w=which(colnames(train_matrix)=="time"); colnames(train_matrix)[w]="time_cutoff"
  train_matrix$time_cutoff=as.numeric(train_matrix$time_cutoff)
  
  
  
  
  
  ##### Join with time matrix
  if (is.character(time_matrix)) time_matrix=readRDS(time_matrix)
  w=which(colnames(time_matrix)=="time"); colnames(time_matrix)[w]="time_cutoff"
  time_matrix$time_cutoff=as.numeric(time_matrix$time_cutoff)
  train_matrix = train_matrix  %>% 
    left_join(time_matrix,by=c("id","time_cutoff"))

  # Clean up
  print("Combined with time matrix")
  rm(time_matrix)
  gc()
  
  
  
  
  
  ##### Join with LTC matrix
  if (is.character(ltc_matrix)) ltc_matrix=readRDS(ltc_matrix)
  w=which(colnames(ltc_matrix)=="time"); colnames(ltc_matrix)[w]="time_cutoff"
  ltc_matrix$ltc_cutoff=as.numeric(ltc_matrix$time_cutoff)
  train_matrix=train_matrix %>% 
    left_join(ltc_matrix,by=c("id","time_cutoff")) %>% 
    replace_na(list(ltc_rawdata_NUMBEROFLTCS=0))
  
  print("Combined with LTC matrix")
  rm(ltc_matrix)
  gc()


  
  ##### Join with topic matrix
  if (is.character(topic_matrix)) topic_matrix=readRDS(topic_matrix)
  w=which(colnames(topic_matrix)=="time"); colnames(topic_matrix)[w]="time_cutoff"
  topic_matrix$topic_cutoff=as.numeric(topic_matrix$time_cutoff)
  
  # Remove some unnecessary variables
  if ("topics_missing_code_total" %in% colnames(topic_matrix)) topic_matrix = topic_matrix %>% select(-c("topics_missing_code_total","topics_missing_code_unique"))
  
  train_matrix=train_matrix %>%
    left_join(topic_matrix,by=c("id","time_cutoff"))
  rm(topic_matrix)
  gc()
  
  print("Combined matrices, adding target")


  ## Reload data
  if (is.null(patients)) {
    load_cleanData(partition = "all",
    subset_frac=1, subset_seed = 1234,
    load_to_global=TRUE,
    load_sparrav3_scores = FALSE)
  }

  print("Loaded data")
  
  
  # Fix column names and add target
  colnames(train_matrix) = make.names(colnames(train_matrix))
  
  train_matrix = 	train_matrix %>% 
    mutate(time=as_datetime(time_cutoff)) %>%
    select(-c("time_cutoff")) %>% 
    mutate(time = as_datetime(time)) %>%
    add_target(episodes, list_of_data_tables, fill_times = FALSE)
  
  if (!keep_id_and_time) train_matrix=train_matrix %>% select(-c("id","time")) 
  
  return(train_matrix)  
}  


##' icd2sick()
##' 
##' Transforms a vector of ICD-10 codes into rough description of disease class
##' @param x vector of strings of ICD10 codes
##' @returns rough description of disease class
icd2sick=function(x0) {
  x1=substr(x0,1,1)
  x2=substr(x0,1,2)
  case_when(
    x1 %in% c("A","B") ~ 'Infectious disease',
    x1 == "C" ~ 'Neoplasm',
    x2 %in% c("D1","D2","D3","D4") ~ 'Neoplasm',
    x2 %in% c("D5","D6","D7","D8","D9") ~ 'Blood',
    x1 == 'E' ~ 'Endocrine/metabolic',
    x1 == 'F' ~ 'Mental/behavioural',
    x1 == 'G' ~ 'Nervous system',
    x2 %in% c('H1','H2','H3','H4','H5') ~ 'Eye',
    x2 %in% c('H6','H7','H8','H9') ~ 'Ear',
    x1 == 'I' ~ 'Circulatory',
    x1 == 'J' ~ 'Respiratory',
    x1 == 'K' ~ 'Digestive',
    x1 == 'L' ~ 'Skin',
    x1 == 'M' ~ 'Musculoskeletal',
    x1 == 'N' ~ 'Genitourinary',
    x1 == 'O' ~ 'Obstetric/puerperium',
    x1 == 'P' ~ 'Perinatal',
    x1 == 'Q' ~ 'Congenital',
    x1 == 'R' ~ 'Abnormality NEC',
    x1 %in% c('S',"T","V","X","Y") ~ 'External',
    x1 %in% c("U","Z") ~ 'Other')
}



##' adcode2lit()
##' 
##' Transforms a vector of admission codes to a literal description of the admission type
##' @param x vector of strings of admission codes
##' @returns rough description of admission types
adcode2lit=function(x) {
  case_when(
    x == 36 ~ 'Patient non-injury',
    x == 30 ~ 'Em. adm., unknown',
    x == 38 ~ 'Other incl transfer',
    x == 33 ~ 'Patient home injury',
    x == 34 ~ 'Patient work injury',
    x == 35 ~ 'Patient other injury',
    x == 20 ~ 'Urg. adm., unknown',
    x == 31 ~ 'Patient self-inflicted',
    x == 39 ~ 'Em. adm., unknown',
    x == 32 ~ 'Road accident',
    x == 22 ~ 'Urg. adm., hospital delay',
    x == 21 ~ 'Urg. adm., patient delay')
}







## General function
##' get_shap()
##' 
##' General-purpose Shapley value estimator
##' @param Xt matrix to generate Shapley values for, assumed to  include target
##' @param fitx function taking a matrix like Xt and returning predictions
##' @param ni number of iterations. About ni/e samples per Shapley value
##' @param verbose set to TRUE to print progress
##' @return matrix of dimension Xt of Shapley values
get_shap=function(Xt,fitx,ni=500,verbose=T) {
  
  shap_tst=0*Xt; shap_n=shap_tst[1,]
  nvar=dim(Xt)[2] # number of variables
  nsamp = dim(Xt)[1] # number of samples
  
  # random subset
  rsub=function(s) s[which(runif(length(s))>0.5)]
  
  for (i in 1:ni) {
    rx=rsub(1:nvar)
    rc=setdiff(1:nvar,rx)
    
    ox=order(runif(nsamp)) # random permutation of samples
    Xv=Xt
    Xv[,rx]=Xv[ox,rx] # replace observations in ox with observations from other samples
    predx=fitx(Xv)
    for (j in 1:length(rc)) {
      rxj=c(rx,rc[j])
      Xv=Xt
      Xv[,rxj]=Xv[ox,rxj] # replace observations in ox with observations from other samples
      predxj=fitx(Xv)		
      shap_tst[,rc[j]]=shap_tst[,rc[j]] + predx-predxj
      shap_n[rc[j]]=1+shap_n[rc[j]]
      if (verbose) print(paste0(j," of ",length(rc)))
    }
    if (verbose) print(paste0("Completed iteration ", i))
  }
  
  shap_a=shap_tst/outer(rep(1,100),as.numeric(shap_n))
  colnames(shap_a)=colnames(Xt)
  
  return(shap_a)
}




## Function to compute Shapley values on split machines.
##' get_shap.split()
##' 
##' General-purpose Shapley value estimator for split machines. Run on linux first, will save approximately ni*nrow(Xt)/2 files, Then run on linux with same parameters.
##' @param Xt matrix to generate Shapley values for, assumed to  include target
##' @param model_id Model ID (file name) as used in SPARRAv4fit.split
##' @param proc_id Unique ID to this procedure which will appear in all filenames
##' @param ni number of iterations. About ni/e samples per Shapley value
##' @param seed random seed
##' @param verbose set to TRUE to print progress
##' @param linux set to TRUE if on linux machine
##' @param remove_when_done set to TRUE to delete auxiliary files after calculating
##' @return matrix of dimension Xt of Shapley values
get_shap.split=function(Xt,model_id,proc_id=NULL,seed=220,ni=500,verbose=T,linux=FALSE,remove_when_done=FALSE) {
  
  # random subset function
  rsub=function(s) s[which(runif(length(s))>0.5)]
  
  if (is.null(proc_id))  proc_id="proc"
  
  model_dir=paste0(dirname(model_id),"/",proc_id)
  model_base=basename(model_id)
  
  if (!file.exists(model_dir)) dir.create(model_dir)
  
  existing=list.files(model_dir,pattern=paste0(model_base,".shapley.*"))
  
  if (!linux & length(existing)==0) stop("Run on linux machine first")
  
  shap_tst=0*Xt; 
  for (i in 1:dim(shap_tst)[2]) shap_tst[which(is.na(shap_tst[,i])),i]=0
  shap_n=shap_tst[1,]
  nvar=dim(Xt)[2] # number of variables
  nsamp = dim(Xt)[1] # number of samples
  
  # Load everything first to speed things up
  xlist=readRDS(paste0(model_id,".pred.RDS"))
  nlist=names(xlist)
  models=gsub("function.","",grep("function.",nlist,value=TRUE))
  models_linux=grep("h2o",models,value=TRUE); models_windows=setdiff(models,models_linux)  
  if (linux) {
    h2o.init(bind_to_localhost = FALSE, max_mem_size = "125G")
    h2o.removeAll()
    for (i in 1:length(models_linux)) xlist[[i+length(nlist)]]=h2o.loadModel(paste0(model_id,".",models_linux[i]))
    names(xlist)=c(nlist,models_linux)
  } else {
    for (i in 1:length(models_windows)) xlist[[i+length(nlist)]]=readRDS(paste0(model_id,".",models_windows[i]))
    names(xlist)=c(nlist,models_windows)
  }
  if (verbose) print("Loaded models and predictor functions")
  
  for (i in 1:ni) {
    iseed=seed+ nvar*i
    set.seed(iseed)
    rx=rsub(1:nvar)
    rc=setdiff(1:nvar,rx)
    
    ox=order(runif(nsamp)) # random permutation of samples
    Xv=Xt
    Xv[,rx]=Xv[ox,rx] # replace observations in ox with observations from other samples
    i_name=paste0(model_dir,"/",model_base,".shapley.",iseed,".",proc_id)
    #if (file.exists(paste0(i_name,".RDS"))) {
    #  predx.saved=readRDS(paste0(i_name,".RDS"))
    #  predx=predx.saved$final.full
    #} else {
      predx=predict.sparra(Xv,model_id,save_id=i_name,verbose=verbose,linux=linux,use_loaded=xlist)
    #}

    Xv=c()
    for (j in 1:length(rc)) {
      jseed=seed+ i*nvar + j
      set.seed(jseed)
      rxj=c(rx,rc[j])
      cmat=Xt
      cmat[,rxj]=cmat[ox,rxj] # replace observations in ox with observations from other samples
      Xv=rbind(Xv,cmat)
    }
    if (verbose) print("Concatenated matrices")

    predj_name=paste0(i_name,".comp")
    #if (file.exists(paste0(predj_name,".RDS"))) {
    #  predj.saved=readRDS(paste0(predj_name,".RDS"))
    #  predj=predj.saved$final.full
    #} else {
      predj=predict.sparra(Xv,model_id,save_id=paste0(i_name,".comp"),verbose=verbose,linux=linux,use_loaded=xlist)
    #}
    
    if (!linux) {
      for (j in 1:length(rc)) {
        predxj=predj[(1+ (j-1)*dim(Xt)[1]):(j*dim(Xt)[1])]
        shap_tst[,rc[j]]=shap_tst[,rc[j]] + predx-predxj
        shap_n[rc[j]]=1+shap_n[rc[j]]
      }
      if (remove_when_done) {
        file.remove(paste0(model_dir,"/",model_id,".",i_name,".RDS"))
        file.remove(paste0(model_dir,"/",model_id,".",i_name,".comp.RDS"))
      }
    }
    
    if (verbose) print(paste0("Completed iteration ", i))
  }
  
  if (!linux) {
    shap_a=shap_tst/outer(rep(1,dim(shap_tst)[1]),as.numeric(shap_n))
    colnames(shap_a)=colnames(Xt)
    if (remove_when_done) dir.remove(model_dir)
    return(list(values=shap_a,counts=shap_n))
  } else return(NULL)
}

## Function to convert machine variable names into long variable names.
##' longvarnames()
##' 
##' Converts machine variable names to long names. Does not do anything clever; add new variables as necessary..
##' @param x variable names
##' @return vector of long form variable names
longvarnames=function(x)
case_when(
x=="emergency_bed_days"~"Number of emergency bed days",
x=="num_emergency_admissions"~"Number of emergency admissions",
x=="elective_bed_days"~"Number of elective bed days",
x=="num_elective_admissions"~"Number of elective admissions",
x=="emergency_drugAndalcohol_admin"~"Number of emergency drug and alcohol related admissions",
x=="numLTCs_resulting_in_admin"~"Number of long-term conditions resulting in admission",
x=="pis_PAID_GIC_INCL_BB"~"Total amount paid in prescription costs",
x=="pis_NUMBER_OF_PAID_ITEMS"~"Total number of paid prescription items",
x=="pis_countBNFsections"~"Number of BNF sections from which a prescription was filled",
x=="pis_Respiratory"~"Number of respiratory related prescriptions",
x=="pis_CentralNervousSystem"~"Number of central nervous system related prescriptions",
x=="pis_Infections"~"Number of infection related prescriptions",
x=="pis_EndocrineSystem"~"Number of endocrine related prescriptions",
x=="pis_IncontinenceDevices"~"Number of incontinence device presciptions",
x=="pis_CombinedStomaDevices"~"Number of stoma device prescriptions",
x=="pis_Anticoagulants_And_Protamine"~"Number of anticoagulant and protamine prescriptions",
x=="pis_Antiepileptic_Drugs"~"Number of antiepileptic prescriptions",
x=="pis_Antifibrinolytic_Drugs._Haemostatics"~"Number of antifibrinolytic and haemostatic prescriptions",
x=="pis_Antisecretory_Drugs_Mucosal_Protectants" ~"Number of antisecretory and mucosal protectant prescriptions",
x=="pis_Antispasmod_Other_Drgs_Alt_Gut_Motility"~"Number of antispasmodic and gut-motility altering prescriptions",
x=="pis_Arm_Sling_Bandages"~"Number of arm sling bandage prescriptions",
x=="pis_Catheters"~"Number of catheter prescriptions",
x=="pis_Corticosteroids_Respiratory"~"Number of respiratory corticosteroid prescriptions",
x=="pis_Dementia"~"Number of dementia related prescriptions",
x=="pis_Drugs_Affecting_Intestinal_Secretions"~"Number of prescriptions for drugs affecting intestinal secretions",
x=="pis_Drugs_Used_In_Diabetes"~"Number of prescriptions for drugs used in diabetes mellitus",
x=="pis_Drugs_Used_In_Neuromuscular_Disorders"~"Number of prescriptions for drugs used in neuromuscular disorders",
x=="pis_Drugs_Used_In_Park.ism_Related_Disorders"~"Number of prescriptiosn for drugs used in Parkinsonism and related disorders",
x=="pis_Drugs_Used_In_Substance_Dependence"~"Number of prescriptions for drugs used in substance dependence",
x=="pis_Fluids_And_Electrolytes"~"Number of fluid and electrolyte prescriptions",
x=="pis_Minerals"~"Number of mineral prescriptions",
x=="pis_Mucolytics"~"Number of mucolytic prescriptions",
x=="pis_Oral_Nutrition"~"Number of oral nutrition prescriptions",
x=="pis_Vitamins"~"Number of vitamin prescriptions",
x=="num_psych_admissions"~"Number of previous psychiatric admissions",
x=="num_ae2_attendances"~"Number of previous A&E attendances",
x=="num_outpatient_appointment_general"~"Number of previous outpatient appointments",
x=="num_outpatient_appointment_mental"~"Number of previous outpatient mental health appointments",
x=="age"~"Age at time cutoff",
x=="sexM"~"Sex",  
x=="SIMD_DECILE_2016_SCT"~"SIMD decile",
x=="parkinsons_indicated"~"Parkinsons disease",
x=="MS_indicated"~"Multiple sclerosis",
x=="epilepsy_indicated"~"Epilepsy",
x=="dementia_indicated"~"Dementia",
x=="days_since_last_AE2"~"Days since last A&E attendance",
x=="days_since_last_SMR00"~"Days since last outpatient attendance",
x=="days_since_last_SMR01"~"Days since last acute/inpatient or day case attendance",
x=="days_since_last_SMR01E"~"Days since last geriatric long stay attendance",
x=="days_since_last_SMR04"~"Days since last mental health and day case attedance",
x=="days_since_last_SMR01_emergency_only"~"Days since last emergency inpatient or day case attendance",
x=="days_since_last_SMR01E_emergency_only"~"Days since last emergency geriatric long stay attendance",
x=="ltc_FIRST_ARTHRITIS_EPISODE_yearssincediag"~"Years since first arthritis diagnosis",
x=="ltc_FIRST_ASTHMA_EPISODE_yearssincediag"~"Years since first asthma diagnosis",
x=="ltc_FIRST_ATRIAL_FIBRILLATION_EPISODE_yearssincediag"~"Years since first atrial fibrillation diagnosis",
x=="ltc_FIRST_CANCER_EPISODE_yearssincediag"~"Years since first cancer diagnosis",
x=="ltc_FIRST_CHRONIC_LIVER_DISEASE_EPISODE_yearssincediag"~"Years since first chronic liver disease diagnosis",
x=="ltc_FIRST_COPD_EPISODE_yearssincediag"~"Years since first chronic obstructive pulmonary disease diagnosis",
x=="ltc_FIRST_DEMENTIA_EPISODE_yearssincediag"~"Years since first dementia diagnosis",
x=="ltc_FIRST_DIABETES_EPISODE_yearssincediag"~"Years since first diabetes mellitus diagnosis",
x=="ltc_FIRST_EPILEPSY_EPISODE_yearssincediag"~"Years since first epilepsy diagnosis",
x=="ltc_FIRST_HEART_DISEASE_EPISODE_yearssincediag"~"Years since first heart disease diagnosis",
x=="ltc_FIRST_HEART_FAILURE_EPISODE_yearssincediag"~"Years since first heart failure diagnosis",
x=="ltc_FIRST_MULTIPLE_SCLEROSIS_EPISODE_yearssincediag"~"Years since first multiple sclerosis diagnosis",
x=="ltc_FIRST_PARKINSON_DISEASE_EPISODE_yearssincediag"~"Years since first Parkinsons disease diagnosis",
x=="ltc_FIRST_RENAL_FAILURE_EPISODE_yearssincediag"~"Years since first renal failure diagnosis",
x=="ltc_FIRST_CEREBROVASCULAR_DISEASE_EPISODE_yearssincediag"~"Years since first cerebrovascular disease diagnosis",
x=="ltc_rawdata_NUMBEROFLTCS"~"Number of recorded long term conditions",
x=="target"~"Emergency admission in year following cutoff date",
)




##' Moving kernel mean
##' 
##' @name mkernm
##' @param xgrid x-values at which to calculate
##' @param xx x input values
##' @param yy y input values
##' @param kern kernel width for Gaussian smoothing (sd)
##' @return kernel-smoothed values at xx, yy and variances
mkernm=function(xgrid,xx,yy,kern=0.5) {
  ytrue=rep(0,res); xtrue=rep(0,res); yvar=ytrue
  for (i in 1:res) {
    wt=dnorm(xgrid[i]-xx,sd=kern)
    s1=sum(wt); s2=sum(wt^2)
    ytrue[i]=sum(wt*yy) / s1
    xtrue[i]=sum(wt*xx) / s1
    yvar[i]= (s1/(s1^2 - s2))*sum(wt*((yy-ytrue[i])^2))
    if (sum(wt/s1)^2 < 1/5) { # THIS IS IMPORTANT: no density estimate is equivalent to a mean of <5 elements
      xtrue[i]=NA
      ytrue[i]=NA
    }
  }
  return(cbind(xtrue,ytrue,yvar))
}

