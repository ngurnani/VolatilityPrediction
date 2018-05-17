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
  alpha <-0
  MAD.m <- matrix(ncol=2,nrow=upn)
  MAD.v <- rep(0,2)
    
    xnplus1<-rep(0,upn)
    #firstinv = 359
    
    for ( i in 1:data.r){
      
      dat <- data[1:(firstinv + i*data.l)]
      novas1<-P.OPTsimple(dat, alp=0,  plot = F)
      for (j in 1:data.l){
        xnplus1[(i-1)*data.l+j]=VARSTABdp_pred(data[1:(firstinv + (i-1)*data.l+j)], novas1, alp = 0)
      }
      
    }
    MAD.m[,1]<-abs(xnplus1-data[(firstinv + 1):length(data)]^2)
    MAD.v[1]<-mean(MAD.m[,1][which(MAD.m[,1]!="NA")])
  
    xnplus2<-rep(0,upn)
    for ( i in 1:data.r){
    
       dat <- data[1:(firstinv + i*data.l)]
       novas1<-novas.estimate(dat,target="normal",alpha=0)
       for (j in 1:data.l){
          xnplus2[(i-1)*data.l+j]=novas.forecasting(data[1:(firstinv + (i-1)*data.l+j)],novas1$zret[1:(firstinv + (i-1)*data.l+j)], alpha = 0 ,novas1$weights,0)$vf1
       }
    
     }
   MAD.m[,2]<-abs(xnplus2-data[(firstinv + 1):length(data)]^2)
   MAD.v[2]<-mean(MAD.m[,2][which(MAD.m[,2]!="NA")])
  
  return(list( MAD=MAD.v, MAD.b=mean(MAD.b)))
}


scores1 <- simple.novas.compare.alpha(dataset[,1], ceiling(length(dataset[,1])/10), data.r = 9, firstinv = 359)
scores2 <- simple.novas.compare.alpha(dataset[1:2000,2], ceiling(length(dataset[1:2000,2])/7), data.r = 6, firstinv = 284 )
scores3 <- simple.novas.compare.alpha(dataset[1:2000,3], ceiling(length(dataset[1:2000,3])/10), data.r = 9, firstinv = 200)



resultmad<-matrix(ncol=4, nrow=2)
resultmad[,2]<-scores1$MAD/scores1$MAD.b
resultmad[,3]<-scores2$MAD/scores2$MAD.b
resultmad[,4]<-scores3$MAD/scores3$MAD.b
resultmad <- round(resultmad, 4)
resultmad[,1]<-c("Simple NoVaS", "Exponential NoVaS")
colnames(resultmad)<-c("Predictor type", "YEN/Dollar","SP500","IBM")

xtable(resultmad, caption="Entries give the empicial mean absolute deviation (MAD) of prediction of squared returns relative to benchmark", digits = 3)







