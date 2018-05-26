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