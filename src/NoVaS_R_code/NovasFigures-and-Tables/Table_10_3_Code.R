rm(list=ls(all=TRUE))


######### Table3(a) The L1-comparison of predictors for GARCH. 

####Input the required data
dataset<-as.matrix(read.csv("3timeseries.csv"))

###Source the required code
source("NoVaS.R")

estimate.quick.garch <- function(x)
{
  Nobs <- length(x)
  mu <- mean(x)
  e2 <- (x-mu)^2
  ares <- arima(e2,order=c(1,0,1),include.mean=T,method="CSS")
  alpha <- -ares$coef[2]
  theta <- ares$coef[1]+ares$coef[2]
  omega <- mean(e2)*(1-alpha-theta)
  garch.param <- as.vector(c(mu,omega,alpha,theta))
  h0 <- omega/(1-alpha-theta) 
  h2 <- filter(omega+theta*lag(e2),alpha,method="recursive",init=h0)
  hf <- as.numeric(omega + alpha*h2[Nobs] + theta*e2[Nobs])
  list(est=garch.param,vol=h2,volf1 =hf, volf2= hf ) 
}


garch.compare.alpha <- function(data, data.l, data.r,  firstinv){
  
  ### Calculate Benchmark forecasts 
  upn <- data.l*data.r
  MAD.b<-rep(0,upn)
  MAD.b<-abs(novas.recmoments(data,mean.rm=T)[(firstinv + 1):length(data),2]^2-data[(firstinv + 1):length(data)]^2)
  
  
  ###Set up required dataframe to save values (MAD)
  MAD.m <- matrix(ncol= 2,nrow=upn)
  MAD.v <- rep(0,2)
  
  
  xnplus1<-rep(0,upn) ### res  for garch(1,1) with normal errors
  xnplus2<-rep(0,upn) ### res for grach(1,1) with t-errors
  
  for ( i in 1:data.r){
    
    dat <- data[1:(firstinv + i*data.l)]
    
    for (j in 1:data.l){
      xnplus1[(i-1)*data.l+j]=estimate.quick.garch(data[1:(firstinv + (i-1)*data.l+j)])$volf1
      xnplus2[(i-1)*data.l+j]=estimate.quick.garch(data[1:(firstinv + (i-1)*data.l+j)])$volf2
    }
    
  }
  N <- length(data)
  zz<- rnorm(N, 0, 1 )
  zz<-zz^2
  mu <- mean(zz, na.rm = T)
  xnplus1 <- xnplus1*mu
  tt<- rt(N,df = 5)
  tt<-tt^2
  mut<- mean(tt, na.rm = T)
  xnplus2 <- xnplus2*mut
  
  MAD.m[,1]<-abs(xnplus1-data[(firstinv + 1):length(data)]^2)
  MAD.v[1]<-mean(MAD.m[,1][which(MAD.m[,1]!="NA")])
  MAD.m[,2]<-abs(xnplus2-data[(firstinv + 1):length(data)]^2)
  MAD.v[2]<-mean(MAD.m[,2][which(MAD.m[,2]!="NA")])
  
  return(list(MAD=MAD.v, MAD.b=mean(MAD.b)))
}


scores1 <- garch.compare.alpha(dataset[,1], ceiling(length(dataset[,1])/10), data.r = 9, firstinv = 359)
scores2 <- garch.compare.alpha(dataset[1:2000,2], ceiling(length(dataset[1:2000,2])/7), data.r = 6, firstinv = 284 )
scores3 <- garch.compare.alpha(dataset[1:2000,3], ceiling(length(dataset[1:2000,3])/10), data.r = 9, firstinv = 200)


resultmad<-matrix(ncol=4, nrow=4)
resultmad[,2]<-c(1.013, 0.989, scores1$MAD/scores1$MAD.b)
resultmad[,3]<-c(1.187, 0.943, scores2$MAD/scores2$MAD.b)
resultmad[,4]<-c(1.132, 0.983, scores3$MAD/scores3$MAD.b)
resultmad<- round(resultmad, 3)
resultmad[,1]<- c("Exponential Smoothing with CV","Eq.(10.14)-AR fit with AIC", "Equation(10.3)-Garch(1,1) with normal errors", "Equation(10.3)-Garch(1,1) with t-errors")
colnames(resultmad)<-c("Predictor type", "YEN/Dollar","SP500","IBM")
xtable(resultmad, caption="Entries give the empirical mean absolute(MAD) of prediction of squared returns relative to benchmark;note that the MAD of prediction achieved by the benchmark was 6.27e-005, 1.56e-004, 2.33e-004 in the three cases Yen/Dollar, S&P500, IBM respectively.", digits=3)
##############
