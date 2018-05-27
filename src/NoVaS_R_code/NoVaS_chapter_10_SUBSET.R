# Contains subset of useful functions from NoVaS_chapter_10.R

########### one step ahead prediction by NoVaS simple ########
VARSTABdp_pred<-function(xx, p, alp)
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
    V <- var(xx[1:(n-1)])
    ww <- ss/sqrt(alp*V+kk/p)    
    
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