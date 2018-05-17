################# Simulation in 10.7 Time-varing NoVaS: robustness against structural breaks ######

# 1. CP-GARCH(1,1)
# 2. TV-GARCH(1,1)
# 3. Methods condisidered were TV-GARCH based on windowed Gaussian MLE, and TV-NoVaS (simple or exponential)

print(Sys.time())

options(digits = 8)

### 1.1 Simulate a CP-GARCH(1,1) model with Gaussian#
#### based on equation (10.5) C = 10^(-5), t<=n/2, A = 0.10 and B = 0.73. 
### for t>n/2, A = 0.05 and B = 0.93. n = 1001. 

simulate.garch.breaks <- function(n,pbs,ctrl.param,ival)
{ # ival is the initial value of Y. # pbs is the time of  break point 
  Nobs <- n
  Y <- rep(0,Nobs)
  h <- rep(0,Nobs)
  z <- rnorm(Nobs)
  omega <- ctrl.param[1,] # VALUE OF C
  alpha <- ctrl.param[2,] # VALUE OF B
  theta <- ctrl.param[3,] # Value OF A
  
  # initial values  
  h[1] <- omega[1]/(1-alpha[1]-theta[1])
  Y[1] <- ival 
  
  j <- 1
  
  for (i in 2:Nobs ) 
  { 
    h[i] <- omega[j] + alpha[j]*h[i-1] + theta[j]*(Y[i-1]^2)
    Y[i] <- sqrt(h[i])*z[i]
   
    if (i == pbs) # the second half starts here when i>=501
    {
      j <- j + 1
    }
  }
  list(level=Y, vol = h)
}


pbs <- 501
ctrl.param1 <- matrix(NA,3,2)
ctrl.param1[1,] <- c(0.00001, 0.00001)

ctrl.param1[2,] <- c(0.73, 0.93)

ctrl.param1[3,] <- c(0.10, 0.05)

# set initial value of Y is close to zero e.g Y1 = 0.0001.

cpg <-simulate.garch.breaks(1001,pbs = 501, ctrl.param1, ival = 0.0001 ) 

cpg<- cpg$level
## cpg stands for change point garch model precess 

#################################################################################
####### #1.2 Simulate a time-varying GARCH(1,1) model with Gaussian
# T = 1001, omega is constant 0.00001

simulate.tv.garch <- function( T,omega, alpha, theta,ival)
{
  Nobs <- T
  Y <- rep(0,Nobs)
  h <- rep(0,Nobs)
  z <- rnorm(Nobs)
  
  h[1] <- omega[1]/(1-alpha[1]-theta[1])
  Y[1] <- ival
  for (i in seq(2,Nobs,1)) 
  {
    h[i] <- omega[i] + alpha[i]*h[i-1] + theta[i]*(Y[i-1]^2)
    Y[i] <- sqrt(h[i])*z[i]
  }
  list(level=Y, vol = h)
}

### theta = A 
theta2 <- rep(0, 1001)
theta2[1]<- 0.10

for (i in 1:1000){
   theta2[i+1]<- theta2[i] - 5e-05
 }

### alpha = B
alpha2 <- rep(0, 1001)
alpha2[1] <- 0.73

for ( i in 1:1000){
  alpha2[i+1] <- alpha2[i] + 2e-04
}

omega2 <- rep(0.00001, 1001)
##### simulate the time varing garch series
tvg <- simulate.tv.garch( T = 1001,omega2, alpha2, theta2, ival = 0.0001 )

tvg <- tvg$level

#tvg stands for time varying garch process
#################### plot the simulated series
par(mfrow = c( 2, 1))

plot(cpg, type = 'l', xlab = "time", main = " Simulated series of CP-GARCH")
plot(tvg, type = 'l', xlab = "time", main = " Simulated series of TV-GARCH")

####################################################
################# under model CP-GARCH 

#################  Quick estimation of a GARCH(1,1) model via ARMA approximation

estimate.quick.garch_with_mean <- function(x)
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
  hf <-0.45*hf
  hf## (h_t+1)^2 = E(Y_t+1^2|F_t)
  
}


######### time_varying garch(1,1) prediction. 

tv_garch_one_step_prediction <- function(xx, b){
  pred_garch <- rep(0, 16)
  j=1 
  for ( i in seq(250, 1000, 50)){
     yy <- xx[(i-b+1):i]
     pred_garch[j] <- estimate.quick.garch_with_mean(yy) #### h_t+1
     j <- j+1
    } 
  pred_garch
}

#######cp-garch-dataset

##### b = 125
gar_res1_125 = tv_garch_one_step_prediction (cpg, b = 125)

##### b = 250
gar_res1_250 = tv_garch_one_step_prediction (cpg, b = 250)


########### tv-garch-data

##### b = 125
gar_res2_125 = tv_garch_one_step_prediction (tvg, b = 125)

##### b = 250
gar_res2_250 = tv_garch_one_step_prediction (tvg, b = 250)

####################################

KURTOSIS<-function(x)
{ # regular definition 
  mm <- mean(x, na.rm = T)
  s2 <- mean((x - mm)^2, na.rm = T)
  kk <- mean((x - mm)^4, na.rm = T)
  kk <- kk/s2^2
  kk
}

##simple NoVaS Algorithm

#########################
VARSTABdp<-function(xx, p)
{ # xx is financial returns series (percentage returns)
  # the output is the transformed series using simple NoVaS
  # note that first p-1 entries of the output are unreliable
  # Here we have is average of p values (not p+1 as in the 2007 paper)
  # i.e., if the algorithm below indicates p=10, then it is a_0 and p=9 more
  n <- length(xx)
  yy <- c(1:(2 * p + n)) * 0
  yy[(2 * p + 1):(2 * p + n)] <- xx
  yy <- as.ts(yy)
  kk <- 0 * yy
  for(i in 0:(p - 1)) {
    kk <- kk + lag(yy^2,  - i)
  }
  zz <- yy/sqrt(kk/p)    # zz_as.vector(zz)
  #   zz_zz[(2*p+1):(2*p+n)]
  zz <- na.omit(zz)
  zz[1] <- sign(zz[1])
  zz
}


########## seclect p
P.OPTsimple<- function(x,  plot = F)
{ # x is financial returns series (percentage returns)
  # search for optimum P in Simple NOVAS 
  zz <- c(1:50)
  c4 <- c(1:50)
  ccc <- c(1:50)
  
  for(i in 1:50) {
    zz[i] <- KURTOSIS(VARSTABdp(x, i))
  }
 
  if(plot == F) {
    # aa=.9999 means no optimum was found because of refined search--redo with refine=F
    aa <- 0.9999
    if(min(zz) < 3 && max(zz) > 3) {
      zmin <- min(abs(zz - 3), na.rm = T)
      zminus3 <- abs(zz - 3)
      aa <- ccc[c4[zminus3 == zmin]]
      aa <- na.omit(aa)
    }
  }
  if(plot == T) {
    aa <- cbind(ccc, zz)
  }
  
  aa 
}


VARSTABdp_pred<-function(xx, p)
{ # xx is financial returns series (percentage returns)
  # the output is the transformed series using simple NoVaS
  # note that first p-1 entries of the output are unreliable
  # Here we have is average of p values (not p+1 as in the 2007 paper)
  # i.e., if the algorithm below indicates p=10, then it is a_0 and p=9 more
  
  n <- length(xx)
  ss <- c(1:(2 * p + n)) * 0
  ss[(2 * p + 1):(2 * p + n)] <- xx
  ss <- as.ts(ss)
  kk <- 0 * ss
  
  for(i in 0:(p - 1)) {
    kk <- kk + lag(ss^2,  - i)
  }
  ww <- ss/sqrt(kk/p)    
  
  ww <- na.omit(ww)
  ww[1] <- sign(ww[1])
  
  kk <- as.vector(kk)
  kk <-kk/p    
  N <- length(kk)
  An2 <-kk[N]  
  
  a0 <- 1/p
  ww2 <- ww^2/(1-a0*ww*ww)
  m2 <- median(ww2, na.rm = T)
  
  res <- m2*An2  
  res
}

########### simlpe tv-novas prediction


snovas_one_step_prediction <- function(xx, b){
  pred_snovas <- rep(0, 16)
  j=1 
  for ( i in seq(250, 1000, 50)){
    yy <- xx[(i-b+1):i]
    
    p <- P.OPTsimple(yy, plot = F)## here p is for a0, ....ap. 
    
    pred_snovas[j] <- VARSTABdp_pred(yy, p)  
    j <- j+1
  } 
  pred_snovas
}

###### under model cp-garch
snovas_res1_125 <- snovas_one_step_prediction(cpg, b = 125)

snovas_res1_250 <- snovas_one_step_prediction(cpg, b = 250)

################### exponential TV- NoVaS
#######get the transformed series using exponential NoVaS.

VARSTABexpoCCcorrected<-function(xx, d, cc = 0, epsilon = 0.01)
{  # xx is financial returns series (percentage returns)
  # the output is the transformed series using exponential NoVaS 
  # with constant cc*var in denominator, and exponent d(it is c in paper)
  ### Corrected on Aug. 3, 2009 
  
  n <- length(xx)
  p <- floor(n/4)
  a <- c(1:(p - 1)) * 0
  for(i in 1:(p - 1)) {
    a[i] <- exp( - d * i) }
  A <- sum(a)
  a0 <- (1 - cc)/(1 + A)
  a <- (a * (1 - cc))/(1 + A)
  a <- a[a > epsilon]
  asum <- a0 + sum(a)
  a0 <- (a0 * (1 - cc))/asum
  a <- (a * (1 - cc))/asum
  p <- length(a) + 1
  yy <- c(1:(2 * p + n)) * 0
  yy[(2 * p + 1):(2 * p + n)] <- xx
  kk <- yy * 0
  aa <- c(1:p) * 0
  aa[1] <- a0
  aa[2:p] <- a
  kk <- filter(yy^2, aa, sides = 1)
  kk <- as.vector(kk)
  V <- var(xx[1:(n-1)],na.rm=T)
  zz <- yy/sqrt(cc * V + kk  )  # Aug 2009 correction 
  zz <- zz[(2 * p + 1):(2 * p + n)] 
  
  zz}

######################### search optimum exponent in expo NOVAS################################
P.OPTexpoCCcorrected.robust<- function(x, CC = 0, refine = T, plot = F)
{ # x is financial returns series (percentage returns)
  # search for optimum exponent in expo NOVAS WITH CONSTANT
  # CCC is the proposed constant (CCC*var in NoVaS denominator,i.e,alpha)
  zz <- c(1:100)
  c4 <- c(1:100)
  ccc <- c(1:100)/100
  ccc <- ccc + 0.01
  if(refine == T) ccc <- ccc/2 + 0.01   # refine==T denotes a refined search
  for(i in 1:100) {
    zz[i] <- KURTOSIS(VARSTABexpoCCcorrected(x, d = ccc[i], cc = CC))
  }
  # return(cbind(ccc, zz))} # for test only
  if(plot == F) {
    # aa=.9999 means no optimum was found because of refined search--redo with refine=F
    aa <- 0.9999
    if(min(zz) < 3 && max(zz) > 3) {
      zmin <- min(abs(zz - 3), na.rm = T)
      zminus3 <- abs(zz - 3)
      aa <- ccc[c4[zminus3 == zmin]]
      aa <- na.omit(aa)
    }
  }
  if(plot == T) {
    aa <- cbind(ccc, zz)
  }
  
  aa 
}


VARSTABexpoCCcorrectedPREDVOL<- function(xx, d, cc = 0, epsilon = 0.01)
{  # xx is financial returns series (percentage returns)
  # as VARSTABexpoCCcorrectedVOL but returns a single volatility for
  # prediction one-step-ahead (not the whole series)
  
  n <- length(xx)
  p <- floor(n/3)
  a <- c(1:(p - 1)) * 0
  for(i in 1:(p - 1)) {
    a[i] <- exp( - d * i)
  }
  A <- sum(a)
  a0 <- (1 - cc)/(1 + A)
  a <- (a * (1 - cc))/(1 + A)
  a <- a[a > epsilon]
  asum <- a0 + sum(a)
  a0 <- (a0 * (1 - cc))/asum
  a <- (a * (1 - cc))/asum
  p <- length(a) + 1
  yy <- c(1:(2 * p + n)) * 0
  yy[(2 * p + 1):(2 * p + n)] <- xx
  kk <- yy * 0    #aa <- c(1:p) * 0
  #aa[1] <- a0
  #aa[2:p] <- a
  #kk <- filter(yy^2, aa, sides = 1)
  kk <- filter(yy^2, a, sides = 1)
  kk <- as.vector(kk)
  V <- var(xx[1:(n-1)],na.rm=T)  
  zz <- cc * V + (kk * (1 - 0)) # Corrected Aug 2009
  An2 <- zz[(2 * p + n)]
  ww <- yy/sqrt(cc * V + kk  )  # Aug 2009 correction 
  w2 <- ww^2/(1-a0*ww^2)
  m1 <- median(w2[(p+1):n], na.rm = T)
  one_step_pre <- m1*An2
  one_step_pre
}

####### exponential tv-novas one step prediction

enovas_one_step_prediction <- function(xx, b){
  pred_enovas <- rep(0, 16)
  j=1 
  for ( i in seq(250, 1000, 50)){
    ss <- xx[(i-b+1):i]
    
    d <- P.OPTexpoCCcorrected.robust(ss, CC = 0, refine = T, plot = F)
    pred_enovas[j]  <- VARSTABexpoCCcorrectedPREDVOL(ss, d, cc = 0, epsilon = 0.01)
    
    j <- j+1
  } 
  pred_enovas
}

##### under cp-garch model

enovas_res1_125 <- enovas_one_step_prediction(cpg, b = 125)

enovas_res1_250 <- enovas_one_step_prediction(cpg, b = 250)

#######under tv-garch model

enovas_res2_125 <- enovas_one_step_prediction(tvg, b = 125)

enovas_res2_250 <- enovas_one_step_prediction(tvg, b = 250)



##########################################################
##### cp-garch loop*500
### tv-garch MAD b=125
cpg_res_tg_125 = array(NA, c(16,500))
cpg_res_sn_125 = array(NA, c(16,500))
cpg_res_en_125 = array(NA, c(16,500))
cpg_res_tg_250 = array(NA, c(16,500))
cpg_res_sn_250 = array(NA, c(16,500))
cpg_res_en_250 = array(NA, c(16,500))


for (count in 1:500){
    cpg_i = simulate.garch.breaks(1001,pbs = 501, ctrl.param1, ival = 0.0001 ) 
    cpg_i <- cpg_i$level
    
    cpg_res_tg_125[, count] <- tv_garch_one_step_prediction (cpg_i, b = 125)
    cpg_res_sn_125[, count] <- snovas_one_step_prediction(cpg_i, b = 125)
    cpg_res_en_125[, count] <- enovas_one_step_prediction(cpg_i, b = 125)
    
    cpg_res_tg_250[, count] <- tv_garch_one_step_prediction (cpg_i, b = 250)
    cpg_res_sn_250[, count] <- snovas_one_step_prediction(cpg_i, b = 250)
    cpg_res_en_250[, count] <- enovas_one_step_prediction(cpg_i, b = 250)
    print (count)
}



#### compute the MAD of 16 points. 

mine_mad <- function(array_in){
  tmp<- rep(NA, 16)
  for (i in 1:16){
    mu <- median(array_in[i,][which(array_in[i,]>0 & array_in[i,]<0.01)], na.rm = T)
    tmp[i]<- mean(abs(array_in[i,][which(array_in[i,]>0 & array_in[i,]<0.01)]-mu))
  }
  tmp
}


gar_mad125_cpg = mine_mad(cpg_res_tg_125)
sn_mad125_cpg<- mine_mad(cpg_res_sn_125)
en_mad125_cpg<- mine_mad(cpg_res_en_125)

gar_mad250_cpg<- mine_mad(cpg_res_tg_250)
sn_mad250_cpg<- mine_mad(cpg_res_sn_250)
en_mad250_cpg<- mine_mad(cpg_res_en_250)




#### plot b=125 for cpg ##### the differences between simple tv-novas and exponent tv-novas are quite tiny
### so the red line and green line seem overlaped. 
setEPS()
postscript("MAD_CP_Garch.eps")
par(mfrow=c(1,2))
plot(seq(251, 1001, 50),gar_mad125_cpg, lty=1, ylim = c(0, max(gar_mad125_cpg)*1.1), xlab = "Updated Points", ylab = "MAD", type = "l", main = "MADs for CP-Garch, b=125")
par(new = T)
plot(seq(251, 1001, 50),sn_mad125_cpg, lty = 2, ylim = c(0, max(gar_mad125_cpg)*1.1),  xlab = "", ylab = " ", type = "l")
par(new = T)
plot(seq(251, 1001, 50),en_mad125_cpg, lty = 3, ylim = c(0, max(gar_mad125_cpg)*1.1),  xlab = "", ylab = " ", type = "l")
legend("topleft",lwd=1,lty=1:3, legend=c("TV-Garch","Simple TV-NoVaS","Exp. TV-NoVaS"), bty="n")


###################### 
#### plot b=250 for cpg
##### the differences between simple tv-novas and exponent tv-novas are quite tiny
### so the red line and green line seem overlaped. 

plot(seq(251, 1001, 50),gar_mad250_cpg, lty = 1, ylim = c(0, max(gar_mad250_cpg)*1.1), xlab = "Updated Points", ylab = "MAD", type = "l", main = "MADs for CP-Garch, b=250")
par(new = T)
plot(seq(251, 1001, 50),sn_mad250_cpg,  lty = 2, ylim = c(0, max(gar_mad250_cpg)*1.1),  xlab = "", ylab = " ", type = "l")
par(new = T)
plot(seq(251, 1001, 50),en_mad250_cpg,  lty = 3, ylim = c(0, max(gar_mad250_cpg)*1.1),  xlab = "", ylab = " ", type = "l")
legend("topleft",lwd=1,lty=1:3, legend=c("TV-Garch","Simple TV-NoVaS","Exp. TV-NoVaS"), bty="n")
dev.off()


#################################################################################

#################################################################################
############ based on datasets from tv-garch

tvg_res_tg_125 = array(NA, c(16,500))
tvg_res_sn_125 = array(NA, c(16,500))
tvg_res_en_125 = array(NA, c(16,500))
tvg_res_tg_250 = array(NA, c(16,500))
tvg_res_sn_250 = array(NA, c(16,500))
tvg_res_en_250 = array(NA, c(16,500))


for (count in 1:500){
  tvg_i = simulate.tv.garch( T = 1001,omega2, alpha2, theta2, ival = 0.0001 )
  tvg_i <- tvg_i$level
  
  tvg_res_tg_125[, count] <- tv_garch_one_step_prediction (tvg_i, b = 125)
  tvg_res_sn_125[, count] <- snovas_one_step_prediction(tvg_i, b = 125)
  tvg_res_en_125[, count] <- enovas_one_step_prediction(tvg_i, b = 125)
  
  tvg_res_tg_250[, count] <- tv_garch_one_step_prediction (tvg_i, b = 250)
  tvg_res_sn_250[, count] <- snovas_one_step_prediction(tvg_i, b = 250)
  tvg_res_en_250[, count] <- enovas_one_step_prediction(tvg_i, b = 250)
  print (count)
}


gar_mad125_tvg = mine_mad(tvg_res_tg_125)
sn_mad125_tvg<- mine_mad(tvg_res_sn_125)
en_mad125_tvg<- mine_mad(tvg_res_en_125)

gar_mad250_tvg<- mine_mad(tvg_res_tg_250)
sn_mad250_tvg<- mine_mad(tvg_res_sn_250)
en_mad250_tvg<- mine_mad(tvg_res_en_250)

############################################################
#### plot b=125 for tvg
setEPS()
postscript("MAD_TV_Garch.eps")
par(mfrow = c(1, 2))
plot(seq(251, 1001, 50),gar_mad125_tvg, lty=1, ylim = c(0, max(gar_mad125_tvg)*1.1), xlab = "Updated Points", ylab = "MAD", type = "l", main = "MADs for TV-Garch, b=125")
par(new = T)
plot(seq(251, 1001, 50),sn_mad125_tvg,  lty=2, ylim = c(0, max(gar_mad125_tvg)*1.1), xlab = "", ylab = "", type = "l", xaxt="n", yaxt="n")
par(new = T)
plot(seq(251, 1001, 50),en_mad125_tvg,  lty=3, ylim = c(0, max(gar_mad125_tvg)*1.1), xlab = "", ylab = "", type = "l", xaxt="n", yaxt="n")
legend("topleft",lty=1:3, legend=c("TV-Garch","Simple TV-Novas","Exp. TV-Novas"), bty="n")
#### plot b=250 for tvg
plot(seq(251, 1001, 50),gar_mad250_tvg, lty=1, ylim = c(0, max(gar_mad250_tvg)*1.1), xlab = "Updated Points", ylab = "MAD", type = "l", main = "MADs for TV-Garch, b=250")
par(new = T)
plot(seq(251, 1001, 50),sn_mad250_tvg,  lty=2, ylim = c(0, max(gar_mad250_tvg)*1.1), xlab = "", ylab = "", type = "l", xaxt="n", yaxt="n")
par(new = T)
plot(seq(251, 1001, 50),en_mad250_tvg,  lty=3, ylim = c(0, max(gar_mad250_tvg)*1.1), xlab = "", ylab = "", type = "l", xaxt="n", yaxt="n")
legend("topleft",lty=1:3, legend=c("TV-Garch","Simple TV-Novas","Exp. TV-Novas"), bty="n")
dev.off()

######################## END ###############################

