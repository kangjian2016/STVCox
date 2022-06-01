### Common parameters
# Ti: survival outcome
# d: event indicator
# z: covariates
# alph: thresholding parameter, previously was lambda.
# rho: penalization coefficient
# B: b spline basis
# gma: coefficient, a vector of p*q
# thta: inner functions
# bta: varying coefficients
# p: number of covariates
# qn: number of basis



# CI for one variable
# sparse confidence interval
STV_CI <- function(theta_grid, sd_grid, alpha){
  L = theta_grid-1.96*sd_grid
  R = theta_grid+1.96*sd_grid
  L1 = L
  R1 = R
  for(i in 1:length(theta_grid)){
    if(L[i]>=alpha | R[i] <= -alpha | (L[i]<=-alpha & R[i] >= alpha)){
      L1[i] = thre_beta_t(theta_grid[i], alph = alpha)-1.96*sd_grid[i]
      R1[i] = thre_beta_t(theta_grid[i], alph = alpha)+1.96*sd_grid[i]
    }
    if(abs(R[i])< alpha & L[i] < - alpha){
      A = (alpha - theta_grid[i])/sd_grid[i]
      B = qnorm(0.05-pnorm(-A),lower.tail = FALSE)
      L1[i] = thre_beta_t(theta_grid[i], alph = alpha)-B*sd_grid[i]
      R1[i] = 0
    }
    if(abs(L[i]) < alpha & R[i]>alpha){
      B = (alpha + theta_grid[i])/sd_grid[i] # 
      A = qnorm(0.05-pnorm(-B),lower.tail = FALSE)
      L1[i] = 0
      R1[i] = thre_beta_t(theta_grid[i], alph = alpha)+A*sd_grid[i]
    }
    if(abs(L[i])<=alpha & abs(R[i]) <= alpha){
      L1[i] = 0
      R1[i] = 0
    }
  }
  return(cbind(L1,R1))
}


#'@importFrom MASS ginv
### Calculate the covariance matrix for STV estimation
cov.coxstv <-function(Ti, d, z, B, B_grid, gma, alph){
  p = ncol(z)
  n = length(Ti)
  qn = ncol(B)
  
  d = d[order(Ti)]
  z = z[order(Ti),]
  B = B[order(Ti),]
  Ti = Ti[order(Ti)]
  
  Gma = matrix(gma,ncol=p)
  thta = matrix(NA, nr = n, nc = p)
  bta = matrix(NA, nr = n, nc = p)
  zbta = rep(0,n)
  for (j in 1:p) {
    thta[,j] =  B%*%Gma[,j]
    bta[,j] = thre.beta(thta[,j],alph[j])
  }
  
  S0 = vapply(1:n, function(x){ sum(exp(z[x:n,]%*%bta[x,]))},FUN.VALUE=numeric(1))
  S1 = t(vapply(1:n, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%bta[x,])), ncol=p)),FUN.VALUE=numeric(p)))
  
  S2 = array(NA, dim = c(p,p,n))
  i2=0
  repeat{
    i2=i2+1
    S2[,i2,]=(vapply(1:n, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%bta[x,]))*z[x:n,i2], ncol=p)),FUN.VALUE=numeric(p)))
    if (i2==p) break
  }
  
  S0_2 = S0[d==1]
  S1_2 = S1[d==1,]
  S2_2 = S2[,,d==1]
  B2 = B[d==1,]
  GVG_pre=vapply(1:sum(d), function(x)kronecker((S2_2[,,x]/S0_2[x] - matrix(S1_2[x,],ncol=1)%*% matrix(S1_2[x,],ncol=p)/(S0_2[x] *S0_2[x])),B2[x,]%*%t(B2[x,]), FUN = "*") ,outer(1:(p*qn),1:(p*qn)))
  GVG=matrix(0, nrow=p*qn,ncol=p*qn)
  j=0
  repeat{
    j=j+1
    GVG[,j]=rowSums(GVG_pre[,j,])
    if (j==p*qn) break
  }
  
  GVG_inv = MASS::ginv(GVG, tol = 1e-5)
  #eigGVG = eigen(GVG_inv)
  #range(eigen(GVG)$values)
  n_grid = nrow(B_grid)
  sd_all = matrix(0,n_grid,p)
  for(j in 1:p){
    ej = rep(0,p)
    ej[j] = 1
    sd_all[,j] = sqrt(sapply(1:n_grid, function(x) t(kronecker(ej,B_grid[x,]))%*%GVG_inv%*%kronecker(ej,B_grid[x,])))
  }
  #range(sd_all)
  return(list(GVG_inv = GVG_inv,sd_all=sd_all))
}


cov.coxreg <-function(Ti, d, z, B, B_grid, gma){
  p = ncol(z)
  n = length(Ti)
  qn = ncol(B)
  
  d = d[order(Ti)]
  z = z[order(Ti),]
  B = B[order(Ti),]
  Ti = Ti[order(Ti)]
  
  Gma = matrix(gma,ncol=p)
  bta = matrix(NA, nr = n, nc = p)
  zbta = rep(0,n)
  for (j in 1:p) {
    bta[,j] =  B%*%Gma[,j]
  }
  
  S0 = vapply(1:n, function(x){ sum(exp(z[x:n,]%*%bta[x,]))},FUN.VALUE=numeric(1))
  S1 = t(vapply(1:n, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%bta[x,])), ncol=p)),FUN.VALUE=numeric(p)))
  
  S2 = array(NA, dim = c(p,p,n))
  i2=0
  repeat{
    i2=i2+1
    S2[,i2,]=(vapply(1:n, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%bta[x,]))*z[x:n,i2], ncol=p)),FUN.VALUE=numeric(p)))
    if (i2==p) break
  }
  
  S0_2 = S0[d==1]
  S1_2 = S1[d==1,]
  S2_2 = S2[,,d==1]
  B2 = B[d==1,]
  GVG_pre=vapply(1:sum(d), function(x)kronecker((S2_2[,,x]/S0_2[x] - matrix(S1_2[x,],ncol=1)%*% matrix(S1_2[x,],ncol=p)/(S0_2[x] *S0_2[x])),B2[x,]%*%t(B2[x,]), FUN = "*") ,outer(1:(p*qn),1:(p*qn)))
  GVG=matrix(0, nrow=p*qn,ncol=p*qn)
  j=0
  repeat{
    j=j+1
    GVG[,j]=rowSums(GVG_pre[,j,])
    if (j==p*qn) break
  }
  
  GVG_inv = MASS::ginv(GVG, tol = 1e-5)
  #eigGVG = eigen(GVG_inv)
  #range(eigen(GVG)$values)
  n_grid = nrow(B_grid)
  sd_all = matrix(0,n_grid,p)
  for(j in 1:p){
    ej = rep(0,p)
    ej[j] = 1
    sd_all[,j] = sqrt(sapply(1:n_grid, function(x) t(kronecker(ej,B_grid[x,]))%*%GVG_inv%*%kronecker(ej,B_grid[x,])))
  }
  #range(sd_all)
  return(list(GVG_inv = GVG_inv,sd_all=sd_all))
}


stv.beta <-function(B, gma, alph){
  qn = ncol(B)
  nc = length(gma)/qn
  h_t = NULL
  for(j in 1:nc){
    h_t = cbind(h_t,thre.beta(B%*%gma[(1+(j-1)*qn):(j*qn)], alph[j]))
  }
  return(h_t)
}

#'@title Soft-threshold time-varying effects in survival model.
#'@param time numeric vector of n observations for event time.
#'@param delta binary (0,1) vector of observations for event indicator.
#'@param z design matrix (n x p) for the covariates .
#'@param time_grid time grid points for estimating the effects.
#'@param qn integer to specify the number of knots in the B-spline.
#'@param maxit integer to specify the maximum number of iterations for optimization.
#'@param bounds vector of two values for the range for knots of the B-spline.
#'@param rho positive numeric value for the tuning parameter for the penalty term
#'@return list object of model fitting results
#'\describe{
#'\item{gma}{basis coeffients}
#'\item{beta_grid}{estimated sparse time varying effects (after soft-thresholding) on the grid points.}
#'\item{theta_grid}{estimated unthresholded time varying effects on the grid points.}
#'\item{sd_grid}{estiamted standard deviation of the unthresholded time varying effect estimates on the grid points.}
#'\item{CI}{list object of 95 percent sparse confidence intervals of each covariate, number of grid points by two matrix for lower CI and upper CI.}
#'\item{Bs}{Basis function.}
#'\item{ini.para}{initial values.}
#'\item{alph}{threhold values.}
#'\item{BFGS.fit}{output of the function optim for maximizing the partial likelihood function.}
#'\item{elapsed}{running time in seconds.}
#'}
#'@author Jian Kang <jiankang@umich.edu> Yuan Yang and Yi Li
#'@examples
#'dat <- sim_data_multi_fun_unstra()
#'test.time <- seq(from = 0, to = 3,length.out = 500)
#'res <- STVCox_fit(time=dat$time,
#'delta=dat$delta,
#'z=dat$z,
#'time_grid = test.time,
#'qn = 8,rho = 0.5)
#'par(mfcol=c(1,3),mar=c(4,2,2,1))
#'xpl <- c(test.time,test.time[length(test.time):1])
#'for(k in 1:3){
#'ypl <- c(res$CI[[k]][,1],res$CI[[k]][length(test.time):1,2])
#'plot(test.time,res$beta_grid[,k],type="n",col="blue",
#'ylim=range(res$CI[[k]]),main=paste("beta",k),xlab="Time")
#'polygon(xpl,ypl,col="gray",lty=0)
#'lines(test.time,dat$fun.list[[k]](test.time),col="red",lwd=3)
#'lines(test.time,res$beta_grid[,k],col="black",lwd=1)
#'legend("top",c("True Effect","Estimate","95% Sparse CI"),lwd=c(3,1,6),col=c("red","black","gray"),lty=1)
#'}
#'@export
STVCox_fit <- function(time, delta, z, time_grid, qn=8,maxit=200,bounds=NULL,
                   rho=0.01){
  if(is.null(bounds)){
    bounds <- range(time)
  }
  Bs = gen_basis_cubic(time,qn,delta,bounds)
  #initialization
  ini.para = ini_para1(time,delta,z,qn)
  theta = ini.para$theta
  alph = ini.para$alph
  #alph[alph<0.5]=1
  #llk.true=likelihood_true_fun(delta,z, time=time,fun.list)
  
  elapsed=proc.time()[3]
  BFGS.fit2=optim(par = theta, lplk_Rcpp_penalty,alph= alph,time=time,delta=delta,z=z,Bs=Bs,rho=rho,
                  method ="BFGS",control = list(maxit=maxit))
  elapsed=proc.time()[3] - elapsed

  #evaluate performance
  time2 = time[delta==1]
  n_grid = length(time_grid)
  knot_set = quantile(time2,prob=seq(1:(qn-4))/(qn-3))
  test.bs = matrix(splines::bs(time_grid,knot=knot_set, intercept=TRUE, degree=3, Boundary.knots = c(0,max(time))), nrow = n_grid)
  sd_result = cov.coxstv(Ti=time, d=delta, z=z, B=Bs, B_grid=test.bs, gma=BFGS.fit2$par, alph=alph)
  STV_gma <- matrix(BFGS.fit2$par,nc=ncol(z))
  STV_theta_grid <- test.bs%*%STV_gma
  CI <- list()
  for(j in 1:ncol(z)){
    CI[[j]] <- STV_CI(STV_theta_grid[,j],sd_result$sd_all[,j],alph[j])
  }
  
  res <- list(gma =  STV_gma,
              beta_grid = stv.beta(test.bs, STV_gma, alph=alph),
              theta_grid = STV_theta_grid,
              sd_grid = sd_result$sd_all,
              CI = CI,
              Bs = Bs,
              ini.para = ini.para,
              alph = alph,
              BFGS.fit = BFGS.fit2,
              elapsed=elapsed)
  return(res)
}



STV_cp_all <- function(theta_grid, sd_grid, true_beta, alph){
  cp_mat = NULL
  zero_mat = NULL
  LR_mat = NULL
  for(j in 1:length(alph)){
    LR = STV_CI(STV_theta_grid[,j], STV_sd_grid[,j], alph[j])
    cp_temp = (LR[,1] <= true_beta[,j]) &  (LR[,2] >= true_beta[,j])
    cp_mat = cbind(cp_mat,as.numeric(cp_temp))
    zero_temp = !((LR[,1] <= 0) &  (LR[,2] >= 0))
    zero_mat = cbind(zero_mat, as.numeric(zero_temp))
    LR_mat = cbind(LR_mat,LR)
  }
  return(list(cp_mat=cp_mat,zero_mat=zero_mat,LR_mat=LR_mat))
}
# CI for one variable
# sparse confidence interval
STV_CI <- function(theta_grid, sd_grid, alph){
  L = theta_grid-1.96*sd_grid
  R = theta_grid+1.96*sd_grid
  L1 = L
  R1 = R
  for(i in 1:length(theta_grid)){
    if(L[i]>=alph | R[i] <= -alph | (L[i]<=-alph & R[i] >= alph)){
      L1[i] = thre.beta(theta_grid[i], alph = alph)-1.96*sd_grid[i]
      R1[i] = thre.beta(theta_grid[i], alph = alph)+1.96*sd_grid[i]
    }
    if(abs(R[i])< alph & L[i] < - alph){
      A = (alph - theta_grid[i])/sd_grid[i]
      B = qnorm(0.05-pnorm(-A),lower.tail = FALSE)
      L1[i] = thre.beta(theta_grid[i], alph = alph)-B*sd_grid[i]
      R1[i] = 0
    }
    if(abs(L[i]) < alph & R[i]>alph){
      B = (alph + theta_grid[i])/sd_grid[i] # 
      A = qnorm(0.05-pnorm(-B),lower.tail = FALSE)
      L1[i] = 0
      R1[i] = thre.beta(theta_grid[i], alph = alph)+A*sd_grid[i]
    }
    if(abs(L[i])<=alph & abs(R[i]) <= alph){
      L1[i] = 0
      R1[i] = 0
    }
  }
  return(cbind(L1,R1))
}



STV_TFPN_all<-function(beta,beta0){
  rslt = matrix(0,ncol = ncol(beta),nrow = 4)
  for(j in 1:ncol(beta)){
    rslt[,j] = TFPN(beta[,j],beta0[,j])
  }
  return(rslt)
}

TFPN <-function(estimate,truth){
  n.zero = sum(truth == 0)
  n.nonzero = sum(truth !=0)
  p.TP = sum(truth !=0 & estimate != 0)/n.nonzero
  p.TN = sum(truth ==0 & estimate == 0)/n.zero
  p.FP = sum(truth ==0 & estimate != 0)/n.zero
  p.FN = sum(truth !=0 & estimate == 0)/n.nonzero
  return(c(p.TP,p.TN,p.FP,p.FN))
}


thre.beta <- function(thta, alph){
  bta = rep(0,length(thta))
  bta[abs(thta)>alph] = (sign(thta)*(abs(thta)-alph))[abs(thta)>alph]
  bta
}
thre.beta = Vectorize(thre.beta)



#'@importFrom survival coxph 
#initialize parameters
ini_para1<-function(time,delta,z,qn){
  
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  alph=abs(constant_beta/2)
  p=ncol(z)
  theta=rep(0,qn*p)
  for(i in 1:p){
    if(constant_beta[i]>0){
      theta[(1:qn)+(i-1)*qn]=matrix(rep(constant_beta[i],qn), nrow=qn, byrow=FALSE)+alph[i]
    }else{
      theta[(1:qn)+(i-1)*qn]=matrix(rep(constant_beta[i],qn), nrow=qn, byrow=FALSE)-alph[i]
    }
  }
  return(list(theta=theta,alph=alph))
}

ini_para2<-function(time,delta,z,qn){
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  alph=rep(max(abs(constant_beta)),ncol(z))
  p=ncol(z)
  theta=rep(0,qn*p)
  for(i in 1:p){
    theta[(1:qn)+(i-1)*qn]=matrix(rep(constant_beta[i],qn), nrow=qn, byrow=FALSE)
  }
  return(list(theta=theta,alph=alph))
}

ini_para3<-function(time,delta,z,qn){
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  alph=rep(max(abs(constant_beta)),ncol(z))
  p=ncol(z)
  theta=rep(0,qn*p)
  for(i in 1:p){
    theta[(1:qn)+(i-1)*qn]=rnorm(qn,mean=constant_beta[i],sd=abs(constant_beta[i])/2)
  }
  return(list(theta=theta,alph=alph))
}

ini_para4<-function(time,delta,z,qn){
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  p=ncol(z)
  alph=rep(10,p)
  theta=rep(0,qn*p)
  for(i in 1:p){
    theta[(1:qn)+(i-1)*qn]= rnorm(qn,mean = constant_beta[i]+alph[i],sd=0.5)
  }
  return(list(theta=theta,alph=alph))
}

#initialize parameters
ini_para5<-function(time,delta,z,qn,alph0 = 2){
  p=ncol(z)
  vfit <- survival::coxph(Surv(time,delta) ~z)
  constant_beta<-vfit$coef
  alph = rep(alph0,p)
  theta = rep(0,qn*p)
  for(i in 1:p){
    if(constant_beta[i]>0){
      theta[(1:qn)+(i-1)*qn]=matrix(rep(constant_beta[i],qn), nrow=qn, byrow=FALSE)+alph[i]
    }else{
      theta[(1:qn)+(i-1)*qn]=matrix(rep(constant_beta[i],qn), nrow=qn, byrow=FALSE)-alph[i]
    }
  }
  return(list(theta = theta,alph=alph))
}


#'@importFrom splines bs
#generate cubic spline basis
gen_basis_cubic<-function(time,qn,delta,w_bounds){
  
  N =  length(delta)
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(qn-4))/(qn-3))
  B0 = splines::bs(time, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = w_bounds)
  Bs =matrix(B0, nrow=N)
  return(Bs)
}

#generate cubic spline basis
gen_basis_cubic2 <-function(time,qn,w_bounds){
  N =  length(time)
  knot_set=quantile(time,prob=seq(1:(qn-4))/(qn-3))
  B0 = splines::bs(time, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = w_bounds)
  Bs = matrix(B0, nrow=N)
  return(Bs)
}

#generate cubic spline basis
gen_basis_linear<-function(w, qn, w_bounds = c(0,max(w))){
  N =  length(w)
  knot_set = quantile(seq(w_bounds[1],w_bounds[2],length.out = 100),(1:(qn-4))/(qn-3))
  extra = min(diff(knot_set),0)
  boundary = c(w_bounds[1]-extra,w_bounds[2]+extra)
  B0 = splines::bs(w, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = boundary)
  B = matrix(B0, nrow=N)
  return(B)
}

# generate covariates with different cov structure

#'@importFrom mvtnorm rmvnorm
covariate_sim <- function(n,p,z.type="Ind",z.para){
  z = matrix(NA,n,p)
  if(z.type=="Ind"){
    for(i in 1:p){
      z[,i] = runif(n,0,1)
    }
  }
  if(z.type == "AR1"){
    z_mat = diag(p) 
    z_mat = z.para^abs(row(z_mat)-col(z_mat)) 
    z = mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = z_mat)
  }
  if(z.type == "CS"){
    z_mat = matrix(z.para, p,p)
    diag(z_mat) = rep(1,p)
    z = mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = z_mat)
  }
  colnames(z) = paste("z",1:p,sep='')
  return(z)
}



#'@title Simulate unstratified survival data with time-varying effects
#'@param N sample size
#'@param seed random seed
#'@param censor.max censoring time
#'@param end.time the study end time
#'@param fun.list a list of time-varying functions for covariates
#'@param gamma0 the baseline hazard (a constant)
#'@param z.type a character indicating correlation structure of covariate variables. 
#'the possible options include: "Ind" for "independent",
#' "AR1" for "order 1 auto regressive" and "CS" for "compound symmetry".
#'@param z.para a parameter to specify the covariate correlation according to "z.type".
#'@return a list object of four elements
#'\describe{
#'\item{delta}{event or consoring indicator}
#'\item{z}{covariate matrix}
#'\item{time}{event time}
#'\item{event.rate}{the event rate}
#'\item{fun.list}{true functions}
#'}
#'@author Jian Kang <jiankang@umich.edu> Yuan Yang, Yi Li
#'@examples 
#'dat <- sim_data_multi_fun_unstra()
#'fit <- with(dat,survival::coxph(Surv(time,delta)~z))
#'print(fit)
#'@export
sim_data_multi_fun_unstra<-function(N = 1000,
                                    seed = 123,
                                    censor.max = 10,
                                    end.time = 3,
                                    fun.list = NULL,
                                    gamma0 = 0.5,
                                    z.type = "Ind", 
                                    z.para = 0.5){
  # censoring ~ U(0,censor.max)
  # end.time: study end time
  set.seed(seed)
  if(is.null(fun.list)){
    f1<-function(x){
      2*(-x^2+3)*(x<=sqrt(3))
    }
    f1 =  Vectorize(f1)
    
    f2<-function(x){
      6*log(x+0.01)*(x>=1)
    }
    f2 = Vectorize(f2)
    
    f3<-function(x){
      -2*(-x^2+3)*(x<=sqrt(3))
    }
    f3 =  Vectorize(f3)
    
    fun.list = list(f1,f2,f3)
  }
  p=length(fun.list)
  ############ generate data ########################
  # Sigma_z1=diag(p)
  # z= rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  z = covariate_sim(N,p, z.type = z.type, z.para = z.para)
  U=runif(N, 0,1)
  time=rep(0,N)
  kerf<-function(x,i){
    ker=0
    for(j in 1:p){
      ker=ker+fun.list[[j]](x)*z[i,j]
    }
    ker
  }
  for (i in 1:N) {
    f=function(t) {
      integrand <- function(x) {gamma0*exp(kerf(x,i))}
      alph=integrate(integrand, lower = 0, upper = t)$value
      alph+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = end.time+1), silent=TRUE))
    if (class(r1) == "try-error"){    
      time[i]= end.time+1 # if there is no root in (0,4), then censored by 4.
    }
    else time[i]=uniroot(f,  lower = 0, upper = end.time+1)$root
  }
  censoring=runif(N,1,censor.max)
  censoring=censoring*(censoring<end.time)+end.time*(censoring>=end.time) # about 1/4 censored at 4
  tcens = (censoring<time) # censoring indicator
  delta = 1-tcens
  time = time*(delta==1)+censoring*(delta==0)
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  return(list(delta = delta, z = z, time = time, 
              event.rate = mean(delta),fun.list=fun.list))
}


#simulate one variable
simulate_data_Kevin1<-function(n.F=100,alph=50,true_fun,seed=12,censor.max=30){
  set.seed(seed)
  #library(timereg)
  #n.F=100
  n_f = rpois(n.F, alph = alph)  #sample size: a lot of 00
  N=sum(n_f)
  p=1
  #gamma = exp(rnorm(F, mean=0, sd=0.4))
  gamma=rep(1,n.F)
  range(gamma)
  gamma_subject=rep(gamma,n_f)
  
  ############generate data########################
  Sigma_z1=diag(p)
  z= mvtnorm::rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  #z[,4]=rbinom(n=N, size=1, prob=0.5)
  
  U=runif(N, 0,1)
  F_pre=1:n.F
  facility=rep(F_pre, n_f)
  
  time=rep(0,N)
  
  for (i in 1:N) {
    f=function(t) {
      integrand <- function(x) {0.3*exp(true_fun(x)*z[i,1])}#3*sin(3*pi*x/4)*(x<3)*z[i,1]
      alph=integrate(integrand, lower = 0, upper = t)$value
      alph+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      time[i]=4
    }
    else time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  
  censoring=runif(N,0,censor.max)
  censoring=censoring*(censoring<3)+3*(censoring>=3)
  tcens=(censoring<time) # censoring indicator
  # rat4=sum(time==4)/N
  # if(rat4<censor.rate){
  #   censor.rate=censor.rate-rat4
  # }else{
  #   censor.rate=0
  # }
  # tcens=(time==4)
  # tcens[time!=4]=rbinom(N-sum(time==4),size=1,censor.rate)
  delta=1-tcens
  # censoring=runif(N,1,5)
  # censoring=censoring*(censoring<4)+4*(censoring>=4)
  time=time*(delta==1)+censoring*(delta==0)
  
  #mean(delta)  # 0.8
  
  delta = delta[order(time)]
  z = z[order(time),]
  facility=facility[order(time)]
  time = time[order(time)]
  return(list(delta=delta,z=z,time=time,facility=facility))
}
#simulate multi-variables
simulate_data_multivar_fun<-function(n.F=100,alph=50,seed=12,censor.max=30,fun.list=fun.list,gamma0=0.5){
  set.seed(seed)
  #library(timereg)
  #n.F=100
  n_f = rpois(n.F, alph = alph)  #sample size: a lot of 00
  N=sum(n_f)
  p=length(fun.list)
  #gamma = exp(rnorm(F, mean=0, sd=0.4))
  # gamma=rep(1,n.F)
  # range(gamma)
  # gamma_subject=rep(gamma,n_f)
  
  ############generate data########################
  Sigma_z1=diag(p)
  z= mvtnorm::rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  #z[,4]=rbinom(n=N, size=1, prob=0.5)
  
  U=runif(N, 0,1)
  F_pre=1:n.F
  facility=rep(F_pre, n_f)
  
  time=rep(0,N)
  
  kerf<-function(x,i){
    ker=0
    for(j in 1:p){
      ker=ker+fun.list[[j]](x)*z[i,j]
    }
    ker
  }
  
  for (i in 1:N) {
    f=function(t) {
      integrand <- function(x) {gamma0*exp(kerf(x,i))}#3*sin(3*pi*x/4)*(x<3)*z[i,1]
      alph=integrate(integrand, lower = 0, upper = t)$value
      alph+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      time[i]=4
    }
    else time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  censoring=runif(N,1,censor.max)
  censoring=censoring*(censoring<3)+3*(censoring>=3)
  tcens=(censoring<time) # censoring indicator
  delta=1-tcens
  time=time*(delta==1)+censoring*(delta==0)
  
  #mean(delta)  # 0.8
  
  delta = delta[order(time)]
  z = z[order(time),]
  facility=facility[order(time)]
  time = time[order(time)]
  return(list(delta=delta,z=z,time=time,facility=facility))
}

simulate_data_multivar_fun2<-function(n.F=100,alph=50,seed=12,censor.max=30,fun.list=fun.list,gamma0=0.5){
  set.seed(seed)
  #library(timereg)
  #n.F=100
  n_f = rpois(n.F, alph = alph)  #sample size: a lot of 00
  N=sum(n_f)
  p=length(fun.list)
  #gamma = exp(rnorm(F, mean=0, sd=0.4))
  # gamma=rep(1,n.F)
  # range(gamma)
  # gamma_subject=rep(gamma,n_f)
  
  ############generate data########################
  Sigma_z1=diag(p)
  z= mvtnorm::rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  #z[,4]=rbinom(n=N, size=1, prob=0.5)
  
  U=runif(N, 0,1)
  F_pre=1:n.F
  facility=rep(F_pre, n_f)
  
  time=rep(0,N)
  
  kerf<-function(x,i){
    ker=0
    for(j in 1:p){
      ker=ker+fun.list[[j]](x)*z[i,j]
    }
    ker
  }
  
  for (i in 1:N) {
    f=function(t) {
      integrand <- function(x) {gamma0*exp(fun.list[[1]](x)*z[i,1]+fun.list[[2]](x)*z[i,2])}#3*sin(3*pi*x/4)*(x<3)*z[i,1]
      alph=integrate(integrand, lower = 0, upper = t)$value
      alph+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      time[i]=4
    }
    else time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  censoring=runif(N,1,censor.max)
  censoring=censoring*(censoring<3)+3*(censoring>=3)
  tcens=(censoring<time) # censoring indicator
  delta=1-tcens
  time=time*(delta==1)+censoring*(delta==0)
  
  #mean(delta)  # 0.8
  
  delta = delta[order(time)]
  z = z[order(time),]
  facility=facility[order(time)]
  time = time[order(time)]
  return(list(delta=delta,z=z,time=time,facility=facility))
}


#simulate multi-variable functions by function list
simulate_data_multivar_fun3<-function(n.F=100,alph=50,seed=12,censor.max=30,fun.list=fun.list,gamma0=0.5){
  set.seed(seed)
  #library(timereg)
  #n.F=100
  n_f = rpois(n.F, alph = alph)  #sample size: a lot of 00
  N=sum(n_f)
  p=length(fun.list)
  
  ############generate data########################
  Sigma_z1=diag(p)
  z= mvtnorm::rmvnorm(N, mean=rep(0,p), sigma=Sigma_z1)
  #z[,4]=rbinom(n=N, size=1, prob=0.5)
  
  U=runif(N, 0,1)
  F_pre=1:n.F
  facility=rep(F_pre, n_f)
  
  time=rep(0,N)
  
  kerf<-function(x,i){
    ker=0
    for(j in 1:p){
      ker=ker+fun.list[[j]](x)*z[i,j]
    }
    ker
  }
  
  for (i in 1:N) {
    f=function(t) {
      integrand <- function(x) {gamma0*exp(kerf(x,i))}#3*sin(3*pi*x/4)*(x<3)*z[i,1]
      alph=integrate(integrand, lower = 0, upper = t)$value
      alph+log(1-U[i])
    }
    r1 <- suppressWarnings(try(uniroot(f,  lower = 0, upper = 4), silent=TRUE))
    if (class(r1) == "try-error"){    
      time[i]=4
    }
    else time[i]=uniroot(f,  lower = 0, upper = 4)$root
  }
  
  censoring=runif(N,1,censor.max)
  censoring=censoring*(censoring<3)+3*(censoring>=3)
  tcens=(censoring<time) # censoring indicator
  delta=1-tcens
  time=time*(delta==1)+censoring*(delta==0)
  
  #mean(delta)  # 0.8
  
  delta = delta[order(time)]
  z = z[order(time),]
  facility=facility[order(time)]
  time = time[order(time)]
  return(list(delta=delta,z=z,time=time,facility=facility))
}


#partial likelihood with multi-variables
lplk.NM3<-function(TV_data,par){
  theta=par
  alph=TV_data$alph
  time=TV_data$time
  delta=TV_data$delta
  z=TV_data$z # z should be a matrix
  z=as.matrix(z)
  Bs=TV_data$Bs
  # delta=delta[order(time)]
  # z=z[order(time),]
  # Bs=Bs[order(time),]
  # time=time[order(time)]
  n=length(time)
  p=ncol(z)
  if(is.null(p)){p=1}
  beta_t=NULL
  for(i in 1:p){
    beta_t=cbind(beta_t,h_beta_t0(tau=0.01,Bs%*%theta[(1:ncol(Bs))+(i-1)*ncol(Bs)],alph =alph[i]))
  }
  beta_t=as.matrix(beta_t)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(as.matrix(z[x:n,])%*%beta_t[x,])),FUN.VALUE=numeric(1))
  # partial_likelihood=sum(delta*((z%*%theta%*%t(B_spline))-log(S0)))
  psum.list=vapply(list, function(x) delta[x]*(z[x,]%*%beta_t[x,]-log(w0[x])),FUN.VALUE=numeric(1))
  partial_likelihood=sum(psum.list)
  return(-partial_likelihood)
}

#partial likelihood with multi-variables for Kevin's method
lplk.KE<-function(TV_data,KE.theta){
  Bs=TV_data$Bs
  z=TV_data$z
  z=as.matrix(z)
  time=TV_data$time
  delta=TV_data$delta
  beta_t=Bs%*%KE.theta
  n=length(time)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(as.matrix(z[x:n,])%*%beta_t[x,])),FUN.VALUE=numeric(1))
  KE.llk=sum(vapply(list, function(x) delta[x]*(z[x,]%*%beta_t[x,]-log(w0[x])),FUN.VALUE=numeric(1)))
  return(KE.llk) 
}


#neagative log partial likelihood for Nelder-Mead algorithm
lplk.NM<-function(NM_data,par){
  theta=matrix(par[1:(length(par)-1)],ncol = 1)
  alph=par[length(par)]
  time=NM_data$time
  delta=NM_data$delta
  z=NM_data$z # z is a vector. 
  n=length(time)
  Bs=NM_data$Bs
  beta_t=g_beta_t(Bs%*%theta,alph =alph)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(z[x:n]*beta_t[x])),FUN.VALUE=numeric(1))
  likelihood=sum(delta*(z*beta_t-log(w0)))
  return(-likelihood)
}

#neagative log partial likelihood for Nelder-Mead algorithm
#use approximation function
lplk.NM1<-function(NM_data,par){
  theta=matrix(par[1:(length(par)-1)],ncol = 1)
  alph=par[length(par)]
  time=NM_data$time
  delta=NM_data$delta
  z=NM_data$z # z is a vector. 
  n=length(time)
  Bs=NM_data$Bs
  beta_t=h_beta_t0(tau=0.01,Bs%*%theta,alph =alph)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(z[x:n]*beta_t[x])),FUN.VALUE=numeric(1))
  likelihood=sum(delta*(z*beta_t-log(w0)))
  return(-likelihood)
}

#neagative log partial likelihood for Nelder-Mead algorithm
#use approximation function
#do not update alph
lplk.NM2<-function(NM_data,par){
  theta=par
  alph=NM_data$alph
  time=NM_data$time
  delta=NM_data$delta
  z=NM_data$z # z is a vector. 
  n=length(time)
  Bs=NM_data$Bs
  beta_t=h_beta_t0(tau=0.01,Bs%*%theta,alph =alph)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(z[x:n]*beta_t[x])),FUN.VALUE=numeric(1))
  likelihood=sum(delta*(z*beta_t-log(w0)))
  return(-likelihood)
}


#true likelihood function
likelihood_true_fun<-function(delta,z, time=time,fun.list){
  z=as.matrix(z)
  n=length(time)
  p=length(fun.list)
  beta_t=NULL
  for(j in 1:p){
    beta_t=cbind(beta_t,fun.list[[j]](time))
  }
  beta_t=as.matrix(beta_t)
  list=1:n
  w0<-vapply(list, function(x) sum(exp(as.matrix(z[x:n,])%*%beta_t[x,])),FUN.VALUE=numeric(1))
  # partial_likelihood=sum(delta*((z%*%theta%*%t(B_spline))-log(S0)))
  partial_likelihood=sum(vapply(list, function(x) delta[x]*(z[x,]%*%beta_t[x,]-log(w0[x])),FUN.VALUE=numeric(1)))
  return(partial_likelihood)
}

#soft thresholding function
g_beta_t0<-function(bt,alph){
  if(bt>alph){ 
    bt-alph
  }else if(bt<(-alph)){
    bt+alph
  }else{
    0
  }
}
g_beta_t=Vectorize(g_beta_t0)

#continuous approximation of soft thresholding function
h_beta_t0<-function(tau, beta_t,alph){
  0.5*(1+2/pi*atan((beta_t-alph)/tau))*(beta_t-alph)+
    0.5*(1-2/pi*atan((beta_t+alph)/tau))*(beta_t+alph)
}
h_beta_t=Vectorize(h_beta_t0)

#Kevin's method, stratified NR
stra_NR<-function(time,delta,z=z,facility=NULL,qn=qn,ini.theta=NULL,max_ite=1000,tol=10e-4){
  #qn: dimension of B spline, should be greater than 4
  #initial values
  library(survival)
  library(splines)
  N=length(time)
  if(is.null(facility)){
    facility = rep(1,N)
  }
  z=as.matrix(z)
  p=ncol(z)
  if(is.null(ini.theta)){
    vfit <- coxph(Surv(time,delta) ~z)
    constant_beta<-vfit$coef
    theta_NR_stratify=matrix(rep(constant_beta, qn), nrow=p, byrow=FALSE)
  }else{
    theta_NR_stratify=ini.theta
  }
  loglik_cox <-vfit$loglik[2]
  # to build B spline
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(qn-4))/(qn-3))
  bs7=bs(time,knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = c(0,max(time)))
  Bs=matrix(bs7, nrow=N) 
  
  llk_all=NULL
  key=0
  repeat{
    key=key+1
    temp=ddloglik_schoenfeld_bs_t_stratify(n=N,delta,z=z,B_spline=Bs, knot=qn, theta=theta_NR_stratify, facility=facility)
    llk_all=c(llk_all,temp$partial_likelihood)
    dist=matrix(solve(temp$GVG)%*%matrix(temp$GR_test,ncol=1), nrow=p, byrow=TRUE)
    theta_NR_stratify=theta_NR_stratify+dist
    if(max(abs(dist))<tol|key>=max_ite) break
  }
  
  test_NR_stratify=rep(0,p+1)
  test_NR_stratify[p+1]=1-pchisq((2*temp$partial_likelihood-2*loglik_cox), p*qn)
  constrain=-diff(diag(qn*1),differences=1)
  j=0
  repeat{
    j=j+1
    theta_NR_j= theta_NR_stratify[j,]
    L2=solve(temp$GVG)[((j-1)*qn+1):((j-1)*qn+qn),((j-1)*qn+1):((j-1)*qn+qn)]
    test_contrast=t(constrain%*%theta_NR_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_NR_j)
    test_NR_stratify[j]=1-pchisq(test_contrast, (qn-1))
    if(j==p) break
  }
  return(list(pvalue=test_NR_stratify,H=temp$GVG,Bs=Bs,theta=theta_NR_stratify,qn=qn,llk_all=llk_all))
}

#used in stra_NR function
ddloglik_schoenfeld_bs_t_stratify = function(n,delta,z, B_spline=Bs, knot=knot,theta=theta, facility=facility){
  z=as.matrix(z)
  p = ncol(z)
  L = matrix(0,nrow=p,ncol=p)
  number_facility=length(unique(facility))
  V=rep(0, p*p*n)
  dim(V)=c(p,p,n)
  schoenfeld=matrix(0,nrow=n,ncol=p)
  partial_likelihood=0
  
  GVG=matrix(0,nrow=p*knot,ncol=p*knot)
  #GV=matrix(0,nrow=p*knot,ncol=p*knot)
  
  GR_test=rep(0,p*knot)
  GR=rep(0,p*knot)
  
  S0=rep(0,n)
  S1=matrix(0,nrow=n,ncol=p)
  S2= rep(0,p*p*n)
  dim(S2)=c(p,p,n)
  F_key=0
  while (F_key<number_facility){
    F_key=F_key+1
    size_temp=length(delta[facility==F_key])
    list=1:size_temp
    z_temp=as.matrix(z[facility==F_key,])
    B_spline_temp=B_spline[facility==F_key,]
    S0[facility==F_key]=vapply(list, function(x) sum(exp(z_temp[x:size_temp,]%*%theta%*%(B_spline_temp[x,]))),FUN.VALUE=numeric(1))
    S1[facility==F_key,]<-t(vapply(list, function(x) colSums(matrix(z_temp[x:size_temp,]*as.vector(exp(z_temp[x:size_temp,]%*%theta%*%(B_spline_temp[x,]))), ncol=p)),FUN.VALUE=numeric(p)))
    
    S2_temp= rep(0,p*p*size_temp)
    dim(S2_temp)=c(p,p,size_temp)
    
    i2=0
    repeat{
      i2=i2+1
      S2_temp[,i2,]=(vapply(list, function(x) colSums(matrix(z_temp[x:size_temp,]*as.vector(exp(z_temp[x:size_temp,]%*%theta%*%(B_spline_temp[x,])))*z_temp[x:size_temp,i2], ncol=p)),FUN.VALUE=numeric(p)))
      if (i2==p) break
    }
    S2[,,facility==F_key]= S2_temp
  }
  
  schoenfeld=delta*(z - S1/S0)
  # partial_likelihood=sum(delta*((z%*%theta%*%t(B_spline))-log(S0)))
  partial_likelihood=sum(vapply(list, function(x) delta[x]*(z[x,]%*%theta%*%B_spline[x,]-log(S0[x])),FUN.VALUE=numeric(1)))
  
  
  #V=vapply(list, function(x)(S2[,,x]/S0[x] - matrix(S1[x,],ncol=1)%*% matrix(S1[x,],ncol=p)/(S0[x] *S0[x])),outer(1:p,1:p))
  S0_2=as.matrix(S0[delta==1])
  S1_2=as.matrix(S1[delta==1,])
  S2_2=S2[,,delta==1]
  dim(S2_2)=c(p,p,sum(delta))
  list2=1:sum(delta)
  B_spline2=B_spline[delta==1,]
  schoenfeld2=as.matrix(schoenfeld[delta==1,])
  
  
  
  
  GVG_pre=vapply(list2, function(x)kronecker((S2_2[,,x]/S0_2[x] - matrix(S1_2[x,],ncol=1)%*% matrix(S1_2[x,],ncol=p)/(S0_2[x] *S0_2[x])),B_spline2[x,]%*%t(B_spline2[x,]), FUN = "*") ,outer(1:(p*knot),1:(p*knot)))
  
  GVG=matrix(0, nrow=p*knot,ncol=p*knot)
  j=0
  repeat{
    j=j+1
    GVG[,j]=rowSums(GVG_pre[,j,])
    if (j==p*knot) break
  }
  
  GR_test=rowSums(vapply(list2, function(x) kronecker(schoenfeld2[x,],B_spline2[x,], FUN = "*"),numeric(p*knot)))
  
  list(schoenfeld=schoenfeld2, GVG=GVG, GR_test=GR_test, partial_likelihood=partial_likelihood) 
} 


#Kevin's unstratified Newton Raphson method
unstra_NR<-function(time,delta,z=z,qn=qn,ini.theta=NULL,max_ite=1000,tol=10e-4){
  #qn: dimension of B spline, should be greater than 4
  #initial values
  library(survival)
  library(splines)
  N=length(time)
  z=as.matrix(z)
  p=ncol(z)
  if(is.null(ini.theta)){
    vfit <- coxph(Surv(time,delta) ~z)
    constant_beta<-vfit$coef
    theta_NR=matrix(rep(constant_beta, qn), nrow=p, byrow=FALSE)
  }else{
    theta_NR=ini.theta
  }
  # to build B spline
  w_bounds = range(time)
  time2=time[delta==1]
  knot_set=quantile(time2,prob=seq(1:(qn-4))/(qn-3))
  bs7=bs(time, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = w_bounds)
  Bs=matrix(bs7, nrow=N) 
  
  key=0
  repeat{
    key=key+1
    temp=ddloglik_schoenfeld_bs_t(n=N,delta,z=z,B_spline=Bs, knot=qn, theta=theta_NR)
    dist=matrix(solve(temp$GVG)%*%matrix(temp$GR_test,ncol=1), nrow=p, byrow=TRUE)
    theta_NR=theta_NR+dist
    if(max(abs(dist))<tol | key>=max_ite) break
  }
  test_NR=rep(0,p+1)
  #test_NR[p+1]=1-pchisq((2*temp$partial_likelihood-2*loglik_cox), p*knot)
  constrain=-diff(diag(qn*1),differences=1)
  j=0
  repeat{
    j=j+1
    theta_NR_j= theta_NR[j,]
    L2=solve(temp$GVG)[((j-1)*qn+1):((j-1)*qn+qn),((j-1)*qn+1):((j-1)*qn+qn)]
    test_contrast=t(constrain%*%theta_NR_j)%*%solve(constrain%*%L2%*%t(constrain))%*%(constrain%*%theta_NR_j)
    test_NR[j]=1-pchisq(test_contrast, (qn-1))
    if(j==p) break
  }
  return(list(pvalue=test_NR,H=temp$GVG,Bs=Bs,theta=theta_NR,knot=knot))
  
}

#used in unstra_NR
ddloglik_schoenfeld_bs_t = function(n,delta,z, B_spline=Bs, knot=knot,theta=theta){
  p = ncol(z)
  L = matrix(0,nrow=p,ncol=p)
  V=rep(0, p*p*n)
  dim(V)=c(p,p,n)
  schoenfeld=matrix(0,nrow=n,ncol=p)
  partial_likelihood=0
  
  GVG=matrix(0,nrow=p*knot,ncol=p*knot)
  #GV=matrix(0,nrow=p*knot,ncol=p*knot)
  
  GR_test=rep(0,p*knot)
  GR=rep(0,p*knot)
  list=1:n
  S0=rep(0,n)
  S1=matrix(0,nrow=n,ncol=p)
  S2= rep(0,p*p*n)
  S0<-vapply(list, function(x) sum(exp(z[x:n,]%*%theta%*%(B_spline[x,]))),FUN.VALUE=numeric(1))
  S1<-matrix(vapply(list, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%theta%*%(B_spline[x,]))), ncol=p)),FUN.VALUE=numeric(p)),nrow=n,ncol=p)#for p=1, when p>1, should recheck this part
  #t(vapply(list, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%theta%*%(B_spline[x,]))), ncol=p)),FUN.VALUE=numeric(p)))
  schoenfeld=delta*(z - S1/S0)
  # partial_likelihood=sum(delta*((z%*%theta%*%t(B_spline))-log(S0)))
  partial_likelihood=sum(vapply(list, function(x) delta[x]*(z[x,]%*%theta%*%B_spline[x,]-log(S0[x])),FUN.VALUE=numeric(1)))
  
  S2= rep(0,p*p*n)
  dim(S2)=c(p,p,n)
  
  i2=0
  repeat{
    i2=i2+1
    S2[,i2,]=(vapply(list, function(x) colSums(matrix(z[x:n,]*as.vector(exp(z[x:n,]%*%theta%*%(B_spline[x,])))*z[x:n,i2], ncol=p)),FUN.VALUE=numeric(p)))
    if (i2==p) break
  }
  
  #V=vapply(list, function(x)(S2[,,x]/S0[x] - matrix(S1[x,],ncol=1)%*% matrix(S1[x,],ncol=p)/(S0[x] *S0[x])),outer(1:p,1:p))
  # S0_2=S0[delta==1]
  # S1_2=S1[delta==1,]
  S0_2=as.matrix(S0[delta==1])
  S1_2=as.matrix(S1[delta==1,])
  S2_2=S2[,,delta==1]
  dim(S2_2)=c(p,p,sum(delta))
  list2=1:sum(delta)
  B_spline2=B_spline[delta==1,]
  #schoenfeld2=schoenfeld[delta==1,]
  schoenfeld2=as.matrix(schoenfeld[delta==1,])
  GVG_pre=vapply(list2, function(x)kronecker((S2_2[,,x]/S0_2[x] - matrix(S1_2[x,],ncol=1)%*% matrix(S1_2[x,],ncol=p)/(S0_2[x] *S0_2[x])),B_spline2[x,]%*%t(B_spline2[x,]), FUN = "*") ,outer(1:(p*knot),1:(p*knot)))
  
  GVG=matrix(0, nrow=p*knot,ncol=p*knot)
  j=0
  repeat{
    j=j+1
    GVG[,j]=rowSums(GVG_pre[,j,])
    if (j==p*knot) break
  }
  
  GR_test=rowSums(vapply(list2, function(x) kronecker(schoenfeld2[x,],B_spline2[x,], FUN = "*"),numeric(p*knot)))
  
  list(schoenfeld=schoenfeld2, GVG=GVG, GR_test=GR_test, partial_likelihood=partial_likelihood) 
} 



#### For data summarization ####
# return a Bool vector D of the same length of B. If A[i] in B[j], then D[j]=T.
find.id <- function(A,B){
  D = rep(F,length(B))
  for (i in A) {
    D[which(B==i)] = T
  }
  return(D)
}

# remove outliers that are out of 5*sd range.
rm.outlier <- function(dat){
  dat = as.matrix(dat)
  bar = mean(dat)+5*sd(dat)
  nc = ncol(dat)
  outliers = rep(F,nrow(dat))
  for (j in 1:nc) {
    outliers = (outliers|dat[,j]>bar)
  }
  return(list(dat = dat[!outliers,],outliers = outliers))
}



# plot the simulation results
plot_md <-function(w_grid,f0,betas,main = NULL,ylab = "Coefficient function",xlab = "w",cex_value=1.5,byrow=T,ylim=NULL){
  par(mar = c(4.0,4.6,2,2))
  if(is.null(ylim)){
    ylim = range(f0,betas)+c(-0.25,0.25)
  }
  plot(w_grid,f0,main = main,ylim = ylim,type = "l",ylab = ylab,xlab = xlab,lty = 1,lwd=3,col="grey",las=1,cex.main=cex_value,cex.axis=cex_value,cex.lab=cex_value) 
  if(byrow){
    md =  apply(betas, 2, median)
    for (j  in 1:nrow(betas)) {
      lines(w_grid,betas[j,],lty=1,lwd=0.5,col="grey80")
    }
  }else{
    md =  apply(betas, 1, median)
    for (j  in 1:ncol(betas)) {
      lines(w_grid,betas[,j],lty=1,lwd=0.5,col="grey80")
    }
  }
  lines(w_grid,f0,lty = 1,lwd=2,col='red')
  lines(w_grid,md,lty = 1,lwd=2)
}

# create mean(sd) table for latex
mean_sd_table <- function(means,sds,digits=3,multiply = 1){
  means = as.matrix(means)
  sds = as.matrix(sds)
  cmeans <- c(means)
  csds <- c(sds)
  msd <- paste(round(cmeans,digits)*multiply," (",round(csds,digits)*multiply,")",sep="") #mean and standard deviation
  tab <- matrix(msd,nc = ncol(means),nr=nrow(means),dimnames=list(rownames(means),colnames(means)))
  return(tab)
}



