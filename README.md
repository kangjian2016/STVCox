# R package STVCox

## Overview 
This package implements the parameter estimation and statistical inferences on sparse time-varying effects in survival models based on a soft-thresholding approach.  

## Install from Github

Install package `STVCox` from Github with 

```r
devtools::install_github("kangjian2016/STVCox")
library(STVCox)
```

## Examples:
```r
# Simulate data
dat <- sim_data_multi_fun_unstra()
# Setup the time grid points
test.time <- seq(from = 0, to = 3,length.out = 500)
# Model fitting
res <- STVCox_fit(time=dat$time,
                  delta=dat$delta,
                  z=dat$z,
                  time_grid = test.time,
                  qn = 8,rho = 0.5)                
# Show results                  
par(mfcol=c(1,3),mar=c(4,2,2,1))
xpl <- c(test.time,test.time[length(test.time):1])
for(k in 1:3){
    ypl <- c(res$CI[[k]][,1],res$CI[[k]][length(test.time):1,2])
    plot(test.time,res$beta_grid[,k],type="n",col="blue",
    ylim=range(res$CI[[k]]),main=paste("beta",k),xlab="Time")
    polygon(xpl,ypl,col="gray",lty=0)
    lines(test.time,dat$fun.list[[k]](test.time),col="red",lwd=3)
    lines(test.time,res$beta_grid[,k],col="black",lwd=1)
    if(k==2){ 
       legend("top",c("True Effect","Estimate","95% Sparse CI"),lwd=c(3,1,6),col=c("red","black","gray"),lty=1)
    }
}
```
