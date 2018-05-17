##Clean UP
rm(list=ls(all=TRUE))

####Input the required data
dataset<-as.matrix(read.csv("3timeseries.csv"))

###Source the required code

source("NoVaS.R")



simple.novas.compare.alpha <- function(data, data.l, data.r,  firstinv){
  
  ### Calculate Benchmark forecasts 
  upn <- data.l*data.r
  MAD.b<-rep(0,upn)
  MAD.b<-abs(novas.recmoments(data,mean.rm=T)[(firstinv + 1):length(data),2]^2-data[(firstinv + 1):length(data)]^2)
  
  
  ###Set up required dataframe to save values (MAD and for parameter value c)
  alpha <-seq(0,0.7,0.1)
  MAD.m <- matrix(ncol=length(alpha),nrow=upn)
  MAD.v <- rep(0,length(alpha))
  c.v  <- rep(0,length(alpha))
  for (k in 1 :length(alpha)){
    
    c.v[k]<-P.OPTsimple(data, alp=alpha[k],  plot = F)
    
    xnplus1<-rep(0,upn)
    #firstinv = 359
    
    for ( i in 1:data.r){
      
      dat <- data[1:(firstinv + i*data.l)]
      novas1<-P.OPTsimple(dat, alp=alpha[k],  plot = F)
      for (j in 1:data.l){
        xnplus1[(i-1)*data.l+j]=VARSTABdp_pred(data[1:(firstinv + (i-1)*data.l+j)], novas1, alp = alpha[k])
      }
      
    }
    MAD.m[,k]<-abs(xnplus1-data[(firstinv + 1):length(data)]^2)
    MAD.v[k]<-mean(MAD.m[,k][which(MAD.m[,k]!="NA")])
  }
  return(list(c=c.v, MAD=MAD.v, MAD.b=mean(MAD.b)))
}


scores1 <- simple.novas.compare.alpha(dataset[,1], ceiling(length(dataset[,1])/10), data.r = 9, firstinv = 359)
scores2 <- simple.novas.compare.alpha(dataset[1:2000,2], ceiling(length(dataset[1:2000,2])/7), data.r = 6, firstinv = 284 )
scores3 <- simple.novas.compare.alpha(dataset[1:2000,3], ceiling(length(dataset[1:2000,3])/10), data.r = 9, firstinv = 200)
alpha <-seq(0,0.7,0.1)

resultc<-matrix(ncol=4, nrow=length(alpha))
resultc[,2]<-scores1$c
resultc[,3]<-scores2$c
resultc[,4]<-scores3$c

resultc[,1]<-alpha
colnames(resultc)<-c("alpha", "YEN/Dollar","SP500","IBM")
xtable(resultc, caption="Entries give the optimal p (that is p+1 in paper) from kurtosis matching with parameter alpha", digits=1)
resultmad<-matrix(ncol=4, nrow=length(alpha))
resultmad[,2]<-scores1$MAD/scores1$MAD.b
resultmad[,3]<-scores2$MAD/scores2$MAD.b
resultmad[,4]<-scores3$MAD/scores3$MAD.b
resultmad[,1]<-alpha
colnames(resultmad)<-c("alpha", "YEN/Dollar","SP500","IBM")
xtable(resultmad, caption="Entries give the MAD of prediction of squared returns relative to benchmark with parameter alpha", digits=3)