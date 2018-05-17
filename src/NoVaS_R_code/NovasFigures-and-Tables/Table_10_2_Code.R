rm(list=ls(all=TRUE))


######### Table10.2 The L2-comparison of predictors for GARCH.  

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
  MSE.b<-rep(0,upn)
  MSE.b<-(novas.recmoments(data,mean.rm=T)[(firstinv + 1):length(data),2]^2-data[(firstinv + 1):length(data)]^2)^2
  
  
  ###Set up required dataframe to save values (MSE)
  MSE.m <- matrix(ncol= 2,nrow=upn)
  MSE.v <- rep(0,2)
  
  
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
  
  MSE.m[,1]<-(xnplus1-data[(firstinv + 1):length(data)]^2)^2
  MSE.v[1]<-mean(MSE.m[,1][which(MSE.m[,1]!="NA")])
  MSE.m[,2]<-(xnplus2-data[(firstinv + 1):length(data)]^2)^2
  MSE.v[2]<-mean(MSE.m[,2][which(MSE.m[,2]!="NA")])
  
  return(list(MSE=MSE.v, MSE.b=mean(MSE.b)))
}


scores1 <- garch.compare.alpha(dataset[,1], ceiling(length(dataset[,1])/10), data.r = 9, firstinv = 359)
scores2 <- garch.compare.alpha(dataset[1:2000,2], ceiling(length(dataset[1:2000,2])/7), data.r = 6, firstinv = 284 )
scores3 <- garch.compare.alpha(dataset[1:2000,3], ceiling(length(dataset[1:2000,3])/10), data.r = 9, firstinv = 200)


resultMSE<-matrix(ncol=4, nrow=4)
resultMSE[,2]<-c(1.065, 1.013, scores1$MSE/scores1$MSE.b)
resultMSE[,3]<-c(1.151, 1.002, scores2$MSE/scores2$MSE.b)
resultMSE[,4]<-c(1.198, 1.034, scores3$MSE/scores3$MSE.b)
resultMSE<- round(resultMSE, 3)
resultMSE[,1]<- c("Exponential Smoothing with CV","Eq.(10.14)-AR fit with AIC", "Equation(10.3)-Garch(1,1) with normal errors", "Equation(10.3)-Garch(1,1) with t-errors")
colnames(resultMSE)<-c("Predictor type", "YEN/Dollar","SP500","IBM")
xtable(resultMSE, caption="Entries give the empirical Mean Squared Error (MSE) of prediction of squared returns relative to benchmark; note that the MSE of prediction achieved by the benchmark was 2.96e-008, 1.70e-006, 1.78e-006 in the three cases Yen/Dollar, S&P500, IBM respectively.", digits=3)
##############
