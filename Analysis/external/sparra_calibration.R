## Quick simulation to check consistency, bias and variance of kernel-based calibration estimator
## James Liley, 1 May 2020

###################################################################
## Exact bias - cubic function                                   ##
###################################################################

# General parameters.
# N=number of observations (large); delta=kernel width (large for high bias)
# g=parameter of sampling distribution for p(Xi), xsamp=values p(Xi)
# kernel=kernel function, assumed symmetric and ergodic so ie K_delta(x,y)=kernel(x-y,width=delta)
N=c(100000); delta=0.2; g=1; 
xsamp=(1/g)-(1/g)*sqrt(runif(N))
kernel=function(x,width=delta[1]) dnorm(x,sd=width)

# PDF of distribution of p(Xi)
fz= function(z) 2*(1/g -z)*g^2 

# Link function: P(Y|Xi)=P(Y|cp(Xi))=Bernoulli(cp(p(Xi))), and derivative of link function. We are trying to estimate this.
cp=function(z) (z + z^2 + z^3)/3
dcp=function(z) (1 + 2*z + 3* z^2)/3

# Sample y
ysamp=rbinom(N,1,prob=cp(xsamp))

# Estimate c_p at these values
xp=seq(0,1,length=100)

# Estimate c_p
ox= outer(xp,xsamp,function(x,y) kernel(x-y,width=delta[1]))
x1= (ox %*% xsamp) # x1[i]= sum_i(K(xp,xsamp[i])*xsamp[i])
y1= (ox %*% ysamp) # y1[i]= sum_i(K(xp,xsamp[i])*ysamp[i])
cpe=xp*(y1/x1) # now y2[i]=hat(c_p)(xp[i])

# True c_p
cpt=cp(xp)

# Bias
bias_num=rowSums(outer(xp, xsamp, 
  function(x,y) kernel(x-y,width=delta[1])*(y*cp(x)-x*cp(y))))
bias=bias_num/x1


# Plot
plot(0,xlim=c(0,1),ylim=c(-0.1,1),type="n",xlab="p(Xi)",ylab="c_p",
  main=expression(paste("Bias if c"[p],"(x)=(x+x"^2,"+x"^3,")/3")))
lines(xp, cpt,lwd=2,col="gray")
lines(xp, cpe,col="red")
lines(xp,cpe+bias,col="black",lty=2)
lines(xp,bias)
legend("topleft",c("True c_p","Est c_p","Bias","Bias-corrected"),
  lwd=c(2,1,1,1),lty=c(1,1,1,2),col=c("gray","red","black","black"),bty="n")

# lines(x1/rowSums(ox),y1/rowSums(ox),col="blue") # Uncomment to also draw parametric version of plot




###################################################################
## Exact bias - linear function                                  ##
###################################################################

# General parameters.
# N=number of observations (large); delta=kernel width (large for high bias)
# g=parameter of sampling distribution for p(Xi), xsamp=values p(Xi)
# kernel=kernel function, assumed symmetric and ergodic so ie K_delta(x,y)=kernel(x-y,width=delta)
N=c(100000); delta=0.2; g=1; 
xsamp=(1/g)-(1/g)*sqrt(runif(N))
kernel=function(x,width=delta[1]) dnorm(x,sd=width)

# PDF of distribution of p(Xi)
fz= function(z) 2*(1/g -z)*g^2 

# Link function: P(Y|Xi)=P(Y|cp(Xi))=Bernoulli(cp(p(Xi))), and derivative of link function. We are trying to estimate this.
cp=function(z) z/2 
dcp=function(z) 1/2

# Sample y
ysamp=rbinom(N,1,prob=cp(xsamp))

# Estimate c_p at these values
xp=seq(0,1,length=100)

# Estimate c_p
ox= outer(xp,xsamp,function(x,y) kernel(x-y,width=delta[1]))
x1= (ox %*% xsamp) # x1[i]= sum_i(K(xp,xsamp[i])*xsamp[i])
y1= (ox %*% ysamp) # y1[i]= sum_i(K(xp,xsamp[i])*ysamp[i])
cpe=xp*(y1/x1) # now y2[i]=hat(c_p)(xp[i])

# True c_p
cpt=cp(xp)

# Bias
bias_num=rowSums(outer(xp, xsamp, 
  function(x,y) kernel(x-y,width=delta[1])*(y*cp(x)-x*cp(y))))
bias=bias_num/x1


# Unbiased if cp(x)= kx
plot(0,xlim=c(0,1),ylim=c(-0.1,1),type="n",xlab="p(Xi)",ylab="c_p",
  main=expression(paste("Bias if c"[p],"(x)=x/2")))
lines(xp, cpt,lwd=2,col="gray")
lines(xp, cpe,col="red")
lines(xp,cpe+bias,col="black",lty=2)
lines(xp,bias)
legend("topleft",c("True c_p","Est c_p","Bias","Bias-corrected"),
  lwd=c(2,1,1,1),lty=c(1,1,1,2),col=c("gray","red","black","black"),bty="n")





###################################################################
## Consistency                                                   ##
###################################################################

# General parameters.
# N=number of observations (large); delta=vector of kernel widths
# g=parameter of sampling distribution for p(Xi), xsamp=values p(Xi)
# kernel=kernel function, assumed symmetric and ergodic so ie K_delta(x,y)=kernel(x-y,width=delta)
N=c(1000000); delta=sqrt(1/seq(10,100,length=10)); g=1; 
xsamp=(1/g)-(1/g)*sqrt(runif(N))
kernel=function(x,width=delta[1]) dnorm(x,sd=width)

# PDF of distribution of p(Xi)
fz= function(z) 2*(1/g -z)*g^2 

# Link function: P(Y|Xi)=P(Y|cp(Xi))=Bernoulli(cp(p(Xi))), and derivative of link function. We are trying to estimate this.
cp=function(z) (z + z^2 + z^3)/3
dcp=function(z) (1 + 2*z + 3* z^2)/3

# Sample y
ysamp=rbinom(N,1,prob=cp(xsamp))


xp=c(0.1,0.3,0.5,0.8) # we'll estimate f(x) at x in xp


estcps=c()
for (i in 1:length(delta)) {
  ox= outer(xp,xsamp,function(x,y) kernel(x-y,width=delta[i]))
  x1= (ox %*% xsamp) # x1[i]= sum_i(K(xp,xsamp[i])*xsamp[i])
  y1= (ox %*% ysamp) # y1[i]= sum_i(K(xp,xsamp[i])*ysamp[i])
  y2=xp*(y1/x1) # now y2[i]=hat(c_p)(xp[i])
  estcps=cbind(estcps,as.vector(y2)) # assemble
}


yt=cp(xp)
estcps=abs(estcps-yt)

# Estimate is consistent as delta->0 and bias is O(1/delta^2). 
plot(0,xlim=range(1/delta^2),ylim=range(1/estcps),type="n",
  xlab=expression(paste("1/delta"^2)),ylab="1/bias",main="Consistency")
for (i in 1:dim(estcps)[1]) lines(1/delta^2,1/estcps[i,],col=i)
legend("topleft",legend=xp,lty=1,col=1:4)







###################################################################
## Standard error and confidence intervals                       ##
###################################################################

# General parameters
# N=number of observations; delta=kernel width, 
# g=parameter of sampling distribution for p(Xi), xsamp=values p(Xi)
# kernel=kernel function, assumed symmetric and ergodic so ie K_delta(x,y)=kernel(x-y,width=delta)
N=1000; delta=0.05; g=1; 
xsamp=(1/g)-(1/g)*sqrt(runif(N))
kernel=function(x,width=delta) dnorm(x,sd=width)

# PDF of distribution of p(Xi)
fz= function(z) 2*(1/g -z)*g^2 

# Link function: P(Y|Xi)=P(Y|cp(Xi))=Bernoulli(cp(p(Xi))), and derivative of link function. We are trying to estimate this.
cp=function(z) (z + z^2 + z^3)/3
dcp=function(z) (1 + 2*z + 3* z^2)/3

y2z=c()

# We compute confidence intervals conditioning on X so we fix xsamp
xp=seq(0,1,length.out=100) # we'll estimate f(x) at x in xp
ox= outer(xp,xsamp,function(x,y) kernel(x-y,width=delta))
x1= (ox %*% xsamp) # x1[i]= sum_i(K(xp,xsamp[i])*xsamp[i])

# Resample Y 100 times to empirically estimate standard error conditioning on X
for (i in 1:100) {
  ysamp=rbinom(N,1,prob=cp(xsamp))
  y1= (ox %*% ysamp) # y1[i]= sum_i(K(xp,xsamp[i])*ysamp[i])
  y2=xp*(y1/x1) # now y2[i]=hat(c_p)(xp[i])
  y2z=rbind(y2z,as.vector(y2)) # assemble
}



# Bias assuming c_p(x) = kx + k2 x^2, k2<1/2
k2=1/2
bias_num=k2*rowSums(outer(xp, xsamp, 
  function(x,y) kernel(x-y,width=delta[1])*(x^2* y  - x* y^2)))
bias=abs(bias_num/x1)


# Now do a single resample of Y, and we'll estimate the confidence envelope
ysamp=rbinom(N,1,prob=f(xsamp))
y1= (ox %*% ysamp) #/rowSums(ox)
y2=xp*(y1/x1)

xx=xsamp*(1-xsamp)
xs2=((ox^2) %*% xx) # Estimated standard error 

yse=(xp^2)*xs2/((ox %*% xsamp)^2)

plot(0,xlim=c(0,1),ylim=c(0,1),type="n",xlab="p(Xi)",ylab="c_p")
for (i in 1:dim(y2z)[1]) lines(xp,y2z[i,],col="gray")
lines(xp, f(xp),lwd=2)
lines(xp, y2,col="red")
lines(xp,y2-bias-sqrt(yse),col="red",lty=2)
lines(xp,y2+bias+sqrt(yse),col="red",lty=2)
legend("topleft",legend=c("True c_p","100 Est c_ps","Est c_p","+/-SE"),
  col=c("black","gray","red","red"),lty=c(1,1,1,2),bty="n")






