# 
# NoVaS Related Functions
#
# All functions written by : Dimitrios D. Thomakos
# 
# updated 6/23/2015
#

# Include required libraries
library(boot)
library(MASS)
library(nortest)
library(quantmod)

#
# Simulation of NoVaS innovations, normal target with squares
#
sim.n1dist<-function(n,a0)
{
	n3 <- floor(1.3 * n)  # corrected 
	Z <- rnorm(n3)
	c <- 1/sqrt(a0)
	W <- Z[abs(Z) <= c]
	U <- W/sqrt(1 - a0 * W^2)
	U <- U[1:n]
	return(U)
}

#
# Simulation of NoVaS innovations, normal target with absolute values
#
sim.n2dist<-function(n,a0)
{
	n3 <- floor(1.3 * n)  # corrected 
	Z <- rnorm(n3)
	c <- 1/a0
	W <- Z[abs(Z) <= c]
	U <- W/(1 - a0 *abs(W))
	U <- U[1:n]
	return(U)
}

#
# Simulation of NoVaS innovations, uniform target with squares
#
sim.n3dist<-function(n,a0)
{
	n3 <- floor(1.3 * n)  # corrected 
	W <- runif(n3,min=-1/sqrt(a0),max=1/sqrt(a0))
	U <- W/sqrt(1 - a0 * W^2)
	U <- U[1:n]
	return(U)
}

#
# Simulation of NoVaS innovations, uniform target with absolute values
#
sim.n4dist<-function(n,a0)
{
	n3 <- floor(1.3 * n)  # corrected 
	W <- runif(n3,min=-1/a0,max=1/a0)
	U <- W/(1 - a0 *abs(W))
	U <- U[1:n]
	return(U)
}

#
# Quick estimation and printing of descriptive statistics
#
novas.describe<-function(x,doprint=FALSE)
{
	N  <- length(x)
	mx <- mean(x,na.rm=T)
	vx <- mean((x-mx)^2,na.rm=T)
	sx <- mean((x-mx)^3,na.rm=T)/(sqrt(vx)^3)
	kx <- mean((x-mx)^4,na.rm=T)/(vx^2)
	res<- list(mx,vx,sx,kx)
	if (doprint == T)
	{
	cat("Sample mean     = ",mx,"\n")
	cat("Sample std.dev. = ",sqrt(vx),"\n")
	cat("Sample skewness = ",sx,"\n")
	cat("Sample kurtosis = ",kx,"\n")
	cat("","\n")
	cat("Std. error of mean     = ",sqrt(vx/N),"\n")
	cat("Std. error of skewness = ",sqrt(6/N),"\n")
	cat("Std. error of kurtosis = ",sqrt(24/N),"\n")
	}
return(res)
}

# Similar
# Compute, and optionally print, descriptive statistics
#
novas.descriptives <- function(data,lags=20,doprint=F,...)
{
	# Convert data to matrix - works with data frames as well
	data <- as.matrix(data)
	# Dimensions of data
	N <- nrow(data)
	K <- ncol(data)
	# Row and column names
	rnames <- c("Obs.","Mean","Median","Std.Dev.","Skewness","(s.e.S.)","Kurtosis","(s.e.K.)","SW-test","(p-value1)",
			"Q-test","(p-value2)","Q^2-test","(p-value3)")
	cnames <- colnames(data) 
	names <- list(rnames,cnames)
	stats <- matrix(0,nrow=14,ncol=K,dimnames=names)
	# Run a loop through all the variables in data
	for (i in seq(1,K,1))
	{
		x  <- data[,i]
		mx <- mean(x,na.rm=T)
		dx <- median(x,na.rm=T)
		vx <- mean((x-mx)^2,na.rm=T)
		sx <- mean((x-mx)^3,na.rm=T)/(sqrt(vx)^3)
		sxs<- sqrt(6/N)
		kx <- mean((x-mx)^4,na.rm=T)/(vx^2)
		kxs<- sqrt(24/N)
		h  <- shapiro.test(x)
		htest <- h$statistic
		hpval <- h$p.value
		b1     <- Box.test(x,lag=lags,type="Ljung-Box")
		b1test <- b1$statistic
		b1pval <- b1$p.value
		b2     <- Box.test(x^2,lag=lags,type="Ljung-Box")
		b2test <- b2$statistic
		b2pval <- b2$p.value
		stats[,i] <- c(N,mx,dx,sqrt(vx),sx,sxs,kx,kxs,htest,hpval,b1test,b1pval,b2test,b2pval)
	}
	# Print if requested
	if (doprint == T) { print.default(stats,...) }
	# Return the matrix of statistics
	return(stats)
}

# Adapted from the respective package
#
mshapiro.test <- function(U,doboot=FALSE,...) 
{
	# Take transpose of matrix to fit into this code
	U <- t(as.matrix(U))
	
	# OK, the rest are the same
	if(!is.matrix(U))
  	stop("U[] is not a matrix with number of columns (sample size) between 3 and 5000")
	n     <- ncol(U)
	if(n < 3 || n > 5000)
	stop("sample size must be between 3 and 5000")
	rng <- range(U)
	rng <- rng[2] - rng[1]
	if(rng == 0)
	stop("all `U[]' are identical")

	# Remove the means
	Us <- apply(U,1,mean)
	R  <- U-Us

	M.1  <- solve(R%*%t(R),tol=1e-18)
	Rmax <- diag(t(R)%*%M.1%*%R)
	C    <- M.1%*%R[,which.max(Rmax)]
	Z    <- t(C)%*%U
	test <- shapiro.test(Z)
	
	# If bootstrap-based p-value is wanted then we have the following
	if (doboot==TRUE)
	{
		# Load what's needed
		# library(boot)
		# source("b.star.R")
		
		# Compute the optimal block length /note we transpose again!!!
		optbl <- b.star(t(U),round=T)[,1]

		# Call the bootstrap function /again we transpose!!!
		out <- tsboot(tseries=t(U),statistic=mshapiro.test.fun,sim="geom",l=optbl,...)
		
		# Compute the test-statistic and p-value
		test <- c(0,0)
		test[1] <- out$t0
		test[2] <- mean(out$t <= out$t0)
	}
	
	# Done, return
	return(test)
}

# The test function for use with bootstraping above!
mshapiro.test.fun <- function(xxx) { as.double(mshapiro.test(xxx)$statistic) }

# Estimate Mardia's multivariate kurtosis
mardia.kurtosis <- function(U)
{
	# Convert to matrix
	U <- as.matrix(U)
	
	# Get the sample size and number of variables
	Nobs <- NROW(U)
	Nvar <- NCOL(U)
	
	# Remove sample means
	Z <- U-matrix(apply(U,2,mean),Nobs,Nvar,byrow=T)
	
	# Compute the inverse covariance matrix
	iS <- solve(cov(Z))
	
	# Compute the test statistic
	test <- sum(diag((Z%*%iS%*%t(Z))^2))/Nobs - (Nvar*(Nvar+2))/sqrt(8*Nvar*(Nvar+2)/Nobs)
	
	# Compute the p-value
	pval <- 2*(1-pnorm(abs(test)))
	
	# Done, return
	return(c(test,pval))
}

#
# Fast line plot
#
plotline <- function(x,add=F,...) 
{ 
	if (add==F) { plot(x,type="l",...)  }
	else { lines(x,...) }
}
	
#
# Compute QQ-plot coordinates for the normal or uniform distribution
#
novas.qqplot <- function(x,target=c("normal","uniform"),minu,maxu,vname="",doplot=TRUE,...)
{
	x <- as.double(x)
	n <- length(x)
	p <- ppoints(n)
	if (target == "normal") { xx <- qnorm(p) } else { xx <- qunif(p,min=minu,max=maxu) }
	yy <- sort(x)
	if (doplot==TRUE)
	{
		qqplot(xx,yy,xlab="Quantiles",ylab="Order Statistics", 
	       	main = paste("Q-Q plot for ",target," distribution - ",vname,sep=""),...)
	       	bb <- (max(yy)-min(yy))/(max(xx)-min(xx))
	       	aa <- min(yy)-bb*min(xx)
		abline(aa,bb,col="red")	
	}
	return(list(xq=xx,yq=yy))
}

#
# Create a chi-square plot for assessing multivariate normality
#
novas.chi.square.plot <- function(X,...)
{
	# Preliminaries
	X <- as.matrix(X)
	Nobs <- NROW(X)
	Nvar <- NCOL(X)
	
	# Compute the distances
	SX <- cov(X)
	D2 <- mahalanobis(X, colMeans(X), SX)
	
	# Return the Q-Q plot and the distances
	qqplot(qchisq(ppoints(Nobs),df=Nvar),D2,xlab="Quantiles",ylab="Squared distances", 
	       main = expression("Q-Q plot of Mahalanobis" * ~D^2 * " vs. quantiles of" * ~ chi^2),...)
	abline(0, 1, col = "red")
	
	return(D2)
}

#
# Compute sample moments
#
novas.moments <- function(x,mean.rm=F)
{
	if (mean.rm == F) { mu <- 0 } else { mu <- mean(x,na.rm=T) }
	rmse <- sqrt(mean((x-mu)^2,na.rm=T))
	mad  <- mean(abs(x-mu),na.rm=T)
	skew <- mean((x-mu)^3,na.rm=T)/(rmse^3)
	kurt <- mean((x-mu)^4,na.rm=T)/(rmse^4)
	return(c(mu,rmse,mad,skew,kurt))
}

#
# Compute recursive moments and absolute moments - note that the function
# does not compute the full sample moments but only up to observation n-1
#
novas.recmoments <- function(x,...)
{
	n <- length(x)
	y <- matrix(0,n,5)
	for (i in seq(2,n,1))
	{		
		y[i,] <- novas.moments(x[1:i-1],...)
	}
	return(y)
}

#
# Plot original series along with recursive variance and kurtosis
#
novas.dataplot <- function(x,rmse,rkurt)
{
	# Make new window
	windows()
	# Prepare the screen
	split.screen(c(2,1))
	split.screen(c(1,2),2)
	screen(1)
	plot(x,type="l",main="Original Series")
	screen(3)
	plot(rmse,type="l",main="Recursive Std. Deviation",col="red",ylab="Std. Dev.")
	screen(4)
	plot(rkurt,type="l",main="Recursive Kurtosis",col="red",ylab="Kurtosis")
	close.screen(all=T)
}

# 
# Auxiliary function to extract the correct recursive moment based the 
# NoVaS transformation type
#
novas.extract.recmoments <- function(x,tfp,...)
{
	data <- novas.recmoments(x,...)
	if (tfp == 0) { recm <- data[,2] }
	if (tfp == 1) { recm <- data[,3] }
	return(recm)
}

#
# Compute simple, exponential and polynomial NoVaS weights
#
novas.weights <- function(p,tp,alpha,theta,thresh)
{
	# Equal weights
	if (tp == 0) 
	{
		g <- rep(1,p+1)/(p+1)	
		b <- g[g > thresh]
		b <- b/sum(b)
	}

	# Exponential weights
	if (tp == 1) 
	{
		# Make initial selection
		c <- theta[1]
		s <- seq(0,p)
		g <- exp(-c*s)/sum(exp(-c*s))
		# Trim the coefficients according to selected threshold
		b <- g[g > thresh]
		b <- b/sum(b)
	}

	# Polynomial weights
	if (tp == 2) 
	{
		# Preliminaries
		s  <- seq(0,p)
		b1 <- (p+1-s)
		b2 <- 2./((p+1)*(p+2))
		g <- b1*b2
		# Trim the coefficients according to selected threshold
		b <- g[g > thresh]
		b <- b/sum(b)
	}

	# Return
	return((1-alpha)*b)	
}

#
# Check function for bounds on NoVaS weights
#
novas.check.weights <- function(b,tfp)
{
	a0 <- b[1]
	v1 <- (tfp == 0) && (a0 > 0.11)
	v2 <- (tfp == 1) && (a0 > 0.33)
	vv <- v1 || v2
	#if (vv == T) { cat("NOTE --> NoVaS weight a0 = ",a0," exceeds recommended value","\n") } 
}

#
# For a given set of weights compute the NoVaS transformation; NOTE: the square root/abs.val. is returned always for volatility!!!
#
novas.transform <- function(x,p,tp,alpha,theta,thresh,recm,tfp)
{
	# Sample size
	n <- length(x)

	# Compute the NoVaS weights
	b <- novas.weights(p,tp,alpha,theta,thresh)
	novas.check.weights(b,tfp)
	M <- length(b)

	# Compute selector of observations
	Sobs <- seq(M,n,1)

	# First transformation based on squares
	if (tfp == 0)
	{
		x2 <- x^2
		g2 <- filter(x2,b,method="conv",sides=1)
		if (alpha > 0) 
		{ 	vol <- sqrt(alpha*(recm[Sobs]^2) + g2[Sobs])
			w2 <- x[Sobs]/vol 
		}
		else 
		{ 	vol <- sqrt(g2[Sobs])
			w2 <- x[Sobs]/vol 
		}
	}

	# Second transformation based on absolute values
	if (tfp == 1)
	{
		x2 <- abs(x)
		g2 <- filter(x2,b,method="conv",sides=1)
		if (alpha > 0) 
		{ 	vol <- (alpha*recm[Sobs] + g2[Sobs])
			w2 <- x[Sobs]/vol 
		}
		else 
		{ 	vol <- g2[Sobs]
			w2 <- x[Sobs]/vol 
		}
	}
	
	# Return
	list(weights=b,g=vol,w=w2)
}

# 
# For a given set of weights etc. evaluate a computed transformation on
# its proximity to the normal and the uniform distribution using various
# distributional measures
#
novas.evaluate <- function(w,target,b,tfp)
{
	# Evaluate the weights
	#novas.check.weights(b,tfp)

	# Number of observations
	n <- length(w)

	# Compute sample moments and extract mean, rmse and kurtosis
	moms <- novas.moments(w)
	mu   <- moms[1]
	rmse <- moms[2]
	kurt <- moms[5]

	# Compute the QQ correlation coefficient and the Kolmogorov-Smirnov
	# test objective function value
	if (target == "normal")
	{
		z <- (w-mu)/rmse
		Q <- sort(z)
		p <- (seq(1,n)-0.5)/n
		F <- qnorm(p)
		RQ <- 1 #- (cor(Q,F)^2)
		kstest <- ks.test(w,"pnorm")
		KS <- kstest$p.value
		EK <- abs(kurt - 3)
	}
	if (target == "uniform")
	{
		a0 <- b[1]
		if (tfp == 0) { r <- 1/sqrt(a0) }
		if (tfp == 1) { r <- 1/a0 }
		Q <- sort(w)
		p <- (seq(1,n)-0.5)/n
		F <- qunif(p,min=-r,max=+r)
		RQ <- 1 #- (cor(Q,F)^2)
		kstest <- ks.test(w,"punif",min=-r,max=+r)
		KS <- kstest$p.value
		EK <- abs(kurt - 1.85) # Fix this to 1.80
	}
	list(of1=EK,of2=RQ,of3=KS)
}

#
# Objective function used for optimizing the parameter in exponential NoVaS
#
novas.expw.objf <- function(par,data,order,wtp,alpha.parameter,
threshold,recmoments,transform,objft,trgt)
{
	# Check that the selected weights are the exponential
	if (wtp != 1) { stop("Caution! You can only optimize for exponential weights") }
 
	# Compute the NoVaS transform
	res <- novas.transform(data,order,wtp,alpha.parameter,par,threshold,recmoments,transform)

	# Get the transform and the weights
	w <- res$w
	b <- res$weights

	# Make the evalution
	eva <- novas.evaluate(w,target=trgt,b,transform)
	if (objft == 0) { obj <- eva$of1 }
	if (objft == 1) { obj <- eva$of2 }
	if (objft == 2) { obj <- -eva$of3 }

	# Return
	return(obj)
}

# 
# Optimizing the parameter in exponential NoVaS
#
novas.expw.opt <- function(x,p,tp,alpha,thresh,recm,tfp,otp,tgt)
{
	opt <- optimize(novas.expw.objf,interval=c(0,1),data=x,order=p,wtp=tp,
alpha.parameter=alpha,threshold=thresh,recmoments=recm,transform=tfp,objft=otp,trgt=tgt)
	return(opt)
}

#
# Objective function used for optimizing the lag in simple/polynomial NoVaS
#
novas.simplew.objf <- function(par,data,wtp,alpha.parameter,
threshold,recmoments,transform,objft,trgt)
{
	# Check that the selected weights are the simple or polynomial
	if (wtp == 1) { stop("Caution! Use simple or polynomial weights for this function") }
 
	# Compute the NoVaS transform
	res <- novas.transform(data,par,wtp,alpha.parameter,0,threshold,recmoments,transform)

	# Get the transform and the weights
	w <- res$w
	b <- res$weights

	# Make the evalution
	eva <- novas.evaluate(w,target=trgt,b,tfp)
	if (objft == 0) { obj <- eva$of1 }
	if (objft == 1) { obj <- eva$of2 }
	if (objft == 2) { obj <- eva$of3 }

	# Return
	return(obj)
}

#
# Optimizing the lags in simple/polynomial NoVaS
#
novas.simplew.opt <- function(x,pmax,tp,alpha,thresh,recm,tfp,otp,tgt)
{
	s <- seq(5,pmax,1)
	u <- lapply(s,FUN=novas.simplew.objf,data=x,wtp=tp,alpha.parameter=alpha,
threshold=thresh,recmoments=recm,transform=tfp,objft=otp,trgt=tgt)
	mini <- which.min(u)
	optp <- mini + 4
	objf <- as.double(u[mini])
	list(slag=optp,value=objf)
}

#
# Create a matrix of sequential lags
#
novas.seqlags <- function(x,p)
{
	n <- length(x)
	y <- matrix(0,n-p,p)
	for (i in seq(1,p,1))
	{
		z <- x[seq(1,n-i,1)]	
		y[,i] <- z[seq(p+1-i,length(z),1)]
	}
	return(y)
}

#
# Compute asymmetric exponential NoVaS weights
#
novas.expw.asym <- function(p,alpha,theta,thresh)
{
	b <- theta[1]
	c <- theta[2]
	s1 <- seq(0,p)
	s2 <- seq(1,p)
	a0 <- 1/( sum(exp(-b*s1)) + sum(exp(-c*s2)) )
	aj <- a0*exp(-b*s2)
	bk <- a0*exp(-c*s2)
	bk <- bk[aj > thresh]
	aj <- aj[aj > thresh]
	g  <- c(a0,aj,bk)
	b  <- g/sum(g)
	list(lagl=length(aj),wt=(1-alpha)*b)
}

#
# For a given set of weights compute the NoVaS transformation incorporating
# asymmetries in the variance
#
novas.transform.asym <- function(x,p,alpha,theta,thresh,recm,tfp)
{
	# Sample size
	n <- length(x)

	# Compute the NoVaS weights with asymmetries
	rw <- novas.expw.asym(p,alpha,theta,thresh)
	p  <- rw$lagl
	b  <- rw$wt
	#novas.check.weights(b,tfp)
	aj <- b[seq(1,p+1)]
	bk <- b[seq(p+2,length(b))]

	# Compute selector of observations
	Sobs <- seq(p+1,n,1)

	# First transformation based on squares
	if (tfp == 0)
	{
		x2 <- x^2
		g1 <- cbind(x2[seq(p+1,n,1)],novas.seqlags(x2,p))
		g2 <- (g1%*%aj) + (g1[,seq(2,p+1)]*(g1[,seq(2,p+1)] < 0))%*%bk
		if (alpha > 0) { w2 <- x[Sobs]/sqrt(alpha*(recm[Sobs]^2) + g2) }
		else { w2 <- x[Sobs]/sqrt(g2) }
	}

	# Second transformation based on absolute values
	if (tfp == 1)
	{
		x2 <- abs(x)
		g1 <- cbind(x2[seq(p+1,n,1)],novas.seqlags(x2,p))
		g2 <- (g1%*%aj) + (g1[,seq(2,p+1)]*(g1[,seq(2,p+1)] < 0))%*%bk
		if (alpha > 0) { w2 <- x[Sobs]/(alpha*recm[Sobs] + g2) }
		else { w2 <- x[Sobs]/g2 }
	}
	
	# Return
	list(slag=p,weights=b,g=g2,w=w2)
}

#
# Objective function used for optimizing the parameter in exponential NoVaS
# with asymmetries
#
# NOTE: use the nlminb or optim functions for optimization !!!
#
novas.expw.objf.asym <- function(par,data,order,alpha.parameter,
threshold,recmoments,transform,objft,trgt)
{
	# Compute the NoVaS transform with asymmetries
	res <- novas.transform.asym(x=data,p=order,alpha=alpha.parameter,theta=par,thresh=threshold,recm=recmoments,tfp=transform)

	# Get the transform and the weights
	w <- res$w
	b <- res$weights

	# Make the evalution
	eva <- novas.evaluate(w,target=trgt,b,tfp)
	if (objft == 0) { obj <- eva$of1 }
	if (objft == 1) { obj <- eva$of2 }
	if (objft == 2) { obj <- eva$of3 }

	# Return
	return(obj)
}

#
# Quick function for optimizing the asymmetric NoVaS method
#
quick.optimize.asym <- function(data,order,alpha.parameter,threshold,recmoments,transform,objft,trgt)
{
	# Create a sequence of values 
	s <- seq(0.01,1.0,0.01)
	l <- length(s)
	# A matrix to hold the results
	result <- matrix(0,nrow=l,ncol=l)
	# Now loop and store
	for (i in s)
	{
		for (j in s)
		{
			result[i,j] <- novas.expw.objf.asym(c(i,j),data,order,alpha.parameter,threshold,
			recmoments,transform,objft,trgt)
		}
	}
	# Locate position of minimum and return parameters
	indx <- which(result=min(result),arr.ind=T)
	b <- indx[1]
	c <- indx[2]
	return(c(b,c))

}

#
# NoVaS Forecasting for uncorrelated w
#
# NOTE: x is assumed to have zero mean
#
novas.forecasting <- function(x,w,alpha,b,tfp) 
{
	# Preparations
	#
	# Sample size
	n <- length(x)
	# Length of weights vector
	M <- length(b)
	# Lag length
	p <- M-1
	# Split the weights vector to separate the weight on current value - NOTE: reverse as well
	b1 <- rev(b[seq(2,M)])
	b0 <- b[1]
	# Compute unconditional variance
	s2 <- mean(x*x)
	# Only the last p x's are used in forming the A component
	xp <- x[seq(n-p+1,n)]
	
	# Using squared x's
	if (tfp == 0) 
	{
		An <- sqrt(alpha*s2 + sum(b1*(xp^2)))
		w1 <- w/sqrt(1-b0*(w*w))
		w2 <- ((w*w)/(1-b0*(w*w)))
		m1 <- median(w1)
		m2 <- median(w2)
		xf <- m1*An	
		x2 <- m2*(An^2)
		vf <- b0*x2 + (An^2)
	}
	
	# Using absolute x's
	if (tfp == 1)
	{
		An <- alpha*s2 + sum(b1*abs(xp))
		wt <- w/(1-b0*abs(w))
		mt <- median(wt)
		xf <- mt*An
		x2 <- median(abs(wt))*An
		vf <- b0*x2 + An
	}
	
	# Return
	list(rf=xf,vf1=x2,vf2=vf)
}	

#
# NoVaS Forecasting for uncorrelated w using asymmetries
#
# NOTE: x is assumed to have zero mean
#
novas.forecasting.asym <- function(x,w,alpha,b,p,tfp) 
{
	# Preparations
	#
	# Sample size
	n <- length(x)
	# Length of weights vector
	M <- length(b)
	# Split the weights vector to separate the weight on current value, past values and past values with
	# asymmetries - NOTE: reverse as well
	b0 <- b[1]
	b1 <- rev(b[seq(2,p+1)])
	c1 <- rev(b[seq(p+2,length(b))])
	# Compute unconditional variance
	s2 <- mean(x*x)
	# Only the last p x's are used in forming the A component
	xp <- x[seq(n-p+1,n)]
	
	# Using squared x's
	if (tfp == 0) 
	{
		An <- sqrt( alpha*s2 + sum(b1*(xp^2)) + sum(c1*(xp^2)*(xp < 0)) )
		w1 <- w/sqrt(1-b0*(w*w))
		w2 <- ((w*w)/(1-b0*(w*w)))
		m1 <- median(w1)
		m2 <- median(w2)
		xf <- m1*An	
		x2 <- m2*(An^2)
		vf <- b0*x2 + (An^2)
	}
	
	# Using absolute x's
	if (tfp == 1)
	{
		An <- alpha*s2 + sum(b1*abs(xp)) + sum(c1*abs(xp)*(xp < 0))
		wt <- w/(1-b0*abs(w))
		mt <- median(wt)
		xf <- mt*An
		x2 <- median(abs(wt))*An
		vf <- b0*x2 + An
	}
	
	# Return
	list(rf=xf,vf1=x2,vf2=vf)
}	

#
# Summary function for quick evaluation and forecasting using exponential NoVaS
#
novas.quick <- function(xi,pmax,tp,alpha,thresh,recm,tfp,otp,tgt)
{

	# Perform exponential NoVaS
	#
	opt <- novas.expw.opt(xi,pmax,tp,alpha,thresh,recm,tfp,otp,tgt)
	theta <- opt$minimum
	objf <- opt$objective

	# Use optimal parameters to compute the transformation 
	res <- novas.transform(xi,pmax,tp,alpha,theta,thresh,recm,tfp)
	b <- res$weights
	w <- res$w

	# Compute the NoVaS forecasts
	forc <- novas.forecasting(xi,w,alpha,b,tfp)
	
	# Return
	list(param=theta,vf=forc$vf2)
}

#
# Similar to the above function but for use with asymmetric NoVaS
#
novas.quick.asym <- function(xi,pmax,alpha,thresh,recm,tfp,otp,tgt,sv)
{

	# Perform exponential NoVaS with asymmetries
	opt <- optim(par=sv,fn=novas.expw.objf.asym,method="BFGS",data=xi,order=pmax,alpha.parameter=alpha,
	threshold=thresh,recmoments=recm,transform=tfp,objft=otp,trgt=tgt)
	theta <- opt$par

	# Use optimal parameters to compute the transformation 
	res <- novas.transform.asym(xi,pmax,alpha,theta,thresh,recm,tfp)
	b <- res$weights
	w <- res$w
	p <- res$slag

	# Compute the NoVaS forecasts with asymmetries
	forc <- novas.forecasting.asym(xi,w,alpha,b,p,tfp) 
	
	# Return
	list(param=theta,vf=forc$vf2)
}

#
# Concise function call to estimate exponential NoVaS with defaults
#
novas.estimate <- function(x,pmax=0,alpha=0,threshold=0.01,var.type=0,opt.type=0,target,doprint=TRUE)
{
	# Always do exponential NoVaS
	tp <- 1
	# Sample size
	n <- length(x)
	# Look at defaults and prepare
	if (pmax <= 0) { pmax <- as.integer(n/4) }
	# Extract the correct rec moments
	if (alpha > 0) { recm.x <- novas.extract.recmoments(x,var.type) }
	else { recm.x <- rep(0,n) }
	# Optimize and then transform
	opt.x <- novas.expw.opt(x,pmax,tp,alpha,threshold,recm.x,var.type,opt.type,target)	
	# Extract optimized values
	theta.x <- opt.x$minimum
	obj.x <- opt.x$objective
	# Use optimal parameters to compute the transformation and evaluate it
	res.x <- novas.transform(x,pmax,tp,alpha,theta.x,threshold,recm.x,var.type)
	b.x <- res.x$weights
	p.x <- length(b.x)-1
	w.x <- res.x$w
	v.x <- res.x$g
	eva <- novas.evaluate(w.x,target,b.x,var.type)
	# Do printing if requested
	if (doprint==TRUE)
	{
		cat("NoVaS results -- estimation assumes zero mean returns","\n")
		cat("","\n")
		cat("Details of estimation: ","\n")
		cat("% of unconditional variance = ",alpha,"\n")
		cat("Threshold value             = ",threshold,"\n")
		cat("Squared (0) or absolute (1)","\n")
		cat("returns used for variance   = ",var.type,"\n")
		cat("Optimizing on excess kurtosis (0)","\n")
		cat("              QQ-correlation  (1)","\n")
		cat("or            KS p-value      (2) = ",opt.type,"\n")
		cat("Target distribution is =",target,"\n")
		cat("","\n")
		cat("Estimate of exponent        = ",theta.x,"\n")
		cat("Value of objective function = ",obj.x,"\n")
		cat("a0 at the optimum           = ",b.x[1],"\n")
		cat("Lag length used             = ",floor(p.x),"\n")
		cat("Value of excess kurtosis    = ",eva$of1,"\n")
		cat("Value of 1-RQ               = ",eva$of2,"\n")
		cat("p-value of K/S statistic    = ",eva$of3,"\n")
	}
	# Done and return
	return(list(parameters=c(theta.x,obj.x,p.x),evaluation=c(eva$of1,eva$of2,eva$of3),weights=b.x,vol=v.x,zret=w.x))
}

#
# Model selection based on the above function
#
novas.model.selection <- function(x,doprint=FALSE,...)
{
	# Make vectors of possible combinations
	vtp <- c(0,1)
	otp <- c(0,1,2)
	tgt <- c("normal","uniform")
	# You will produce a 3 * 3 matrix, one row for each optimization type 
	# and the columns corresponding to the three measures
	results1 <- matrix(0,nrow=3,ncol=3)
	results2 <- matrix(0,nrow=3,ncol=3)
	results3 <- matrix(0,nrow=3,ncol=3)
	results4 <- matrix(0,nrow=3,ncol=3)	
	# OK, just do it brute force
	out11 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[1],target=tgt[1],doprint=FALSE,...)
	out21 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[2],target=tgt[1],doprint=FALSE,...)
	out31 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[3],target=tgt[1],doprint=FALSE,...)

	out12 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[1],target=tgt[1],doprint=FALSE,...)
	out22 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[2],target=tgt[1],doprint=FALSE,...)
	out32 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[3],target=tgt[1],doprint=FALSE,...)

	out13 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[1],target=tgt[2],doprint=FALSE,...)
	out23 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[2],target=tgt[2],doprint=FALSE,...)
	out33 <- novas.estimate(x,var.type=vtp[1],opt.type=otp[3],target=tgt[2],doprint=FALSE,...)

	out14 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[1],target=tgt[2],doprint=FALSE,...)
	out24 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[2],target=tgt[2],doprint=FALSE,...)
	out34 <- novas.estimate(x,var.type=vtp[2],opt.type=otp[3],target=tgt[2],doprint=FALSE,...)
	# Place them in the results matrix
	results1[1,] <- out11$evaluation
	results1[2,] <- out21$evaluation
	results1[3,] <- out31$evaluation
	
	results2[1,] <- out12$evaluation
	results2[2,] <- out22$evaluation
	results2[3,] <- out23$evaluation
	
	results3[1,] <- out13$evaluation
	results3[2,] <- out23$evaluation
	results3[3,] <- out33$evaluation
	
	results4[1,] <- out14$evaluation
	results4[2,] <- out24$evaluation
	results4[3,] <- out34$evaluation
	# GIven them row and column names
	colnames(results1) <- c("Exc.Kurt./NT","QQ-Corr/NT","KS-pvalue/NT")
	rownames(results1) <- c("Opt. on Exc.Kurt./SR","Opt. on QQ-Corr/SR","Opt. on KS-pvalue/SR")
	colnames(results2) <- c("Exc.Kurt./NT","QQ-Corr/NT","KS-pvalue/NT")
	rownames(results2) <- c("Opt. on Exc.Kurt./AR","Opt. on QQ-Corr/AR","Opt. on KS-pvalue/AR")
	#	
	colnames(results3) <- c("Exc.Kurt./UT","QQ-Corr/UT","KS-pvalue/UT")
	rownames(results3) <- c("Opt. on Exc.Kurt./SR","Opt. on QQ-Corr/SR","Opt. on KS-pvalue/SR")
	colnames(results4) <- c("Exc.Kurt./UT","QQ-Corr/UT","KS-pvalue/UT")
	rownames(results4) <- c("Opt. on Exc.Kurt./AR","Opt. on QQ-Corr/AR","Opt. on KS-pvalue/AR")
	# The comparison is based on their diagonals but only from Exc. Kurt and QQ-Corr!!!
	d1 <- diag(results1)
	d2 <- diag(results2)
	d3 <- diag(results3)
	d4 <- diag(results4)
	dd <- cbind(d1,d2,d3,d4)
	dd1<- dd[1:2,]
	dd2<- dd[3,]
	colnames(dd) <- c("ntsq","ntab","utsq","utab")
	rownames(dd) <- c("Exc. Kurtosis","QQ-Correlation","KS-pvalue")
	# Do the ranking
	rnk <- apply(dd1,1,which.min)
	# and return the optimized parameter combinations
	if (rnk[1]==1) { opc1 <- list(vtp=0,otp=0,tgt="normal") }
	if (rnk[1]==2) { opc1 <- list(vtp=1,otp=0,tgt="normal") }
	if (rnk[1]==3) { opc1 <- list(vtp=0,otp=0,tgt="uniform") }
	if (rnk[1]==4) { opc1 <- list(vtp=1,otp=0,tgt="uniform") }

	if (rnk[2]==1) { opc2 <- list(vtp=0,otp=1,tgt="normal") }
	if (rnk[2]==2) { opc2 <- list(vtp=1,otp=1,tgt="normal") }
	if (rnk[2]==3) { opc2 <- list(vtp=0,otp=1,tgt="uniform") }
	if (rnk[2]==4) { opc2 <- list(vtp=1,otp=1,tgt="uniform") }
	# Again,
	rnk <- which.max(dd2)
	if (rnk==1) { opc3 <- list(vtp=0,otp=2,tgt="normal") }
	if (rnk==2) { opc3 <- list(vtp=1,otp=2,tgt="normal") }
	if (rnk==3) { opc3 <- list(vtp=0,otp=2,tgt="uniform") }
	if (rnk==4) { opc3 <- list(vtp=1,otp=2,tgt="uniform") }
	
	# Done, return them all
	return(list(ntsq=results1,ntab=results2,utsq=results3,utab=results4,msc=dd,optm1=opc1,optm2=opc2,optm3=opc3))
}

# Supportive
novas.align <- function(x,y)
{
	nx <- length(x)
	ny <- length(y)
	if (nx >= ny) { xx <- x[(nx-ny+1):nx]; yy <- y }
	if (nx <  ny) { xx <- x; yy <- y[(ny-nx+1):ny] }
	return(cbind(xx,yy))
}

# 
# Simulate a GARCH(1,1) model with Gaussian or t-innovations
#
simulate.garch <- function(T,ctrl.param,ival,df)
{
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	mu <- ctrl.param[1]
	omega <- ctrl.param[2]
	alpha <- ctrl.param[3]
	theta <- ctrl.param[4]
	v[1] <- omega/(1-alpha-theta)
	X[1] <- ival[1]
	for (i in seq(2,Nobs,1)) 
	{
		v[i] <- omega + alpha*v[i-1] + theta*(X[i-1]^2)
		X[i] <- mu + sqrt(v[i])*w[i]
	}
	list(level=X,vol=v)
}

# 
# Simulate a GARCH(1,1) model with NoVaS innovations
#
simulate.garch.novas <- function(T,ctrl.param,ival,dt,a0)
{
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	if (dt == 1) { w <- sim.n1dist(Nobs,a0) }
	if (dt == 2) { w <- sim.n2dist(Nobs,a0) }
	if (dt == 3) { w <- sim.n3dist(Nobs,a0) }
	if (dt == 4) { w <- sim.n4dist(Nobs,a0) }
	mu <- ctrl.param[1]
	omega <- ctrl.param[2]
	alpha <- ctrl.param[3]
	theta <- ctrl.param[4]
	v[1] <- omega/(1-alpha-theta)
	X[1] <- ival[1]
	for (i in seq(2,Nobs,1)) 
	{
		v[i] <- omega + alpha*v[i-1] + theta*(X[i-1]^2)
		X[i] <- mu + sqrt(v[i])*w[i]
	}
	list(level=X,vol=v)
}

# 
# Simulate a GARCH(1,1) model with breaks with Gaussian or t-innovations
#

simulate.garch.breaks <- function(T,pb,ctrl.param,ival,df)
{
  Nobs <- T
  X <- rep(0,Nobs)
  v <- rep(0,Nobs)
  if (df == 0) w <- rnorm(Nobs)
  else w <- rt(Nobs,df)/sqrt(df/(df-2))
  nb <- ncol(ctrl.param)
  mu <- ctrl.param[1,]
  omega <- ctrl.param[2,]
  alpha <- ctrl.param[3,]
  theta <- ctrl.param[4,]
  v[1] <- omega[1]/(1-alpha[1]-theta[1])
  X[1] <- ival[1]
  j <- 1
  pbs <- pb[j]
  for (i in seq(2,Nobs,1)) 
  {
    v[i] <- omega[j] + alpha[j]*v[i-1] + theta[j]*(X[i-1]^2)
    X[i] <- mu[j] + sqrt(v[i])*w[i]
    if (i > pbs)
    {
      j <- j + 1
      pbs <- pb[j]
    }
  }
  list(level=X,vol=v)
}


#
# Create sinusoidal parameters for time-varying GARCH
#
tvpar.garch <- function(T,w1,w2,o1,o2,o3,o4) 
{
	t <- seq(T)
	p <- matrix(0,T,3)
	t1 <- 0
	t2 <- 0
	t3 <- 0
	k1 <- length(w1)
	k2 <- length(w2)
	for (i in seq(k1))
	{
		t1 <- t1 + sin(2*pi*w1[i]*t)
	}
	t1 <- (t1-min(t1)+o1)*o2
	for (i in seq(k2))
	{
		t2 <- t2 + sin(2*pi*w2[i]*t)
	}
	t2 <- (t2-min(t2)+o3)
	t2 <- (t2/max(t2))
	t3 <- 1-t2-runif(T,min=0,max=o4)
	p[,1] <- t1
	p[,2] <- t2
	p[,3] <- t3
	return(p)
}

#
# Simulate a GARCH(1,1) with a deterministic trend in the volatility
#
simulate.garch.time <- function(T,ctrl.param,ival,df)
{
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	mu <- ctrl.param[1]
	omega <- ctrl.param[2]
	alpha <- ctrl.param[3]
	theta <- ctrl.param[4]
	v[1] <- omega/(1-alpha-theta)
	X[1] <- ival[1]
	for (i in seq(2,Nobs,1)) 
	{
		v[i] <- omega + alpha*v[i-1] + theta*(X[i-1]^2)
		X[i] <- mu + (alpha+theta-(theta/alpha)*(i/Nobs))*sqrt(v[i])*w[i]
	}
	list(level=X,vol=v)
}

# 
# Simulate a time-varying GARCH(1,1) model with Gaussian or t-innovations
#
simulate.tv.garch <- function(T,tvpar,ival,df)
{
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	mu <- tvpar[,1]
	omega <- tvpar[,2]
	alpha <- tvpar[,3]
	theta <- tvpar[,4]
	v[1] <- omega[1]/(1-alpha[1]-theta[1])
	X[1] <- ival[1]
	for (i in seq(2,Nobs,1)) 
	{
		v[i] <- omega[i] + alpha[i]*v[i-1] + theta[i]*(X[i-1]^2)
		X[i] <- mu[i] + sqrt(v[i])*w[i]
	}
	list(level=X,vol=v)
}

#
# Simulate a Markov switching GARCH model with Gaussian or t-innovations
#
simulate.ms.garch <- function(T,ctrl.param,transm,ival,df)
{
	# Preliminaries
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	S <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	# Split state parameters
	mu1 <- ctrl.param[1,1]
	omega1 <- ctrl.param[2,1]
	alpha1 <- ctrl.param[3,1]
	theta1 <- ctrl.param[4,1]
	mu2 <- ctrl.param[1,2]
	omega2 <- ctrl.param[2,2]
	alpha2 <- ctrl.param[3,2]
	theta2 <- ctrl.param[4,2]	
	# Get elements of transition matrix
	p11 <- transm[1,1]
	p12 <- transm[1,2]
	p21 <- transm[2,1]
	p22 <- transm[2,2]
	# Check conformity of transition probabilities
	if ((p11+p12) != 1. | (p21+p22) != 1.) { stop("Transition probabilities do not sum to 1!") }
	# Get initial stationary distribution
	p1 <- p21/(p21+p12)
	p2 <- p12/(p21+p12)
	# Get initial state
	U0 <- runif(1)
	if (U0 <= p1) { S0 = 1. }
	else { S0 = 2. }
	# Get initial values
	X[1] <- ival[1]
	if (S0 == 1.) { v[1] <- omega1/(1-alpha1-theta1) }
	if (S0 == 2.) { v[1] <- omega2/(1-alpha2-theta2) }
	S[1] <- S0
	# Here is the main loop
	for (i in seq(2,Nobs,1))
	{
		U1 <- runif(1)
		C1 <- cumsum(transm[S0,])<=U1
		S1 <- sum(C1)+1
		if (S1 == 1.) 
		{ 
			v[i] <- omega1 + alpha1*v[i-1] + theta1*(X[i-1]^2) 
			X[i] <- mu1 + sqrt(v[i])*w[i]
		}
		if (S1 == 2.) 
		{ 
			v[i] <- omega2 + alpha2*v[i-1] + theta2*(X[i-1]^2) 
			X[i] <- mu2 + sqrt(v[i])*w[i]
		}
		S0 <- S1
		S[i] <- S0
	}
	list(level=X,vol=v,state=S)
}

#
# Simulate a smooth-transition GARCH with Gaussian or t-innovations
#
simulate.st.garch <- function(T,ctrl.param,g,ival,df)
{
	# Preliminaries
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	S <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	# Split state parameters
	mu1 <- ctrl.param[1,1]
	omega1 <- ctrl.param[2,1]
	alpha1 <- ctrl.param[3,1]
	theta1 <- ctrl.param[4,1]
	mu2 <- ctrl.param[1,2]
	omega2 <- ctrl.param[2,2]
	alpha2 <- ctrl.param[3,2]
	theta2 <- ctrl.param[4,2]	
	# Get initial values
	X[1] <- ival[1]
	L0 <- 1/(1+exp(-g*X[1]))
	S[1] <- L0
	v[1] <- (1-L0)*(omega1/(1-alpha1-theta1)) + L0*(omega2/(1-alpha2-theta2))
	# Here is the main loop
	for (i in seq(2,Nobs,1))
	{
		Li <- 1/(1+exp(-g*X[i-1]))
		h1 <- omega1 + alpha1*v[i-1] + theta1*(X[i-1]^2)
		h2 <- omega2 + alpha2*v[i-1] + theta2*(X[i-1]^2)
		v[i] <- (1-Li)*h1 + Li*h2
		X[i] <- (1-Li)*mu1 + Li*mu2 + sqrt(v[i])*w[i]	
		S[i] <- Li		
	}
	list(level=X,vol=v,state=S)
}

# 
# Simulate a Cox-Ingersoll-Ross (CIR) type diffusion
#
# The equation being simulated is dX(t) = k*(a - X(t))*dt + s*(X(t)^g)*W(t)
#
# where (k,a,s,g) are parameters and W(t) standard Brownian motion
#
# Inputs: T - scalar, number of years, months, weeks, days etc.
#     delta - scalar, sampling frequency per unit of T; must be less than one (e.g.
#             1/12 for monthly data, 1/250 for daily etc.)
#          M - scalar, step size factor; set to M >= 1
# ctrl.param - (4 * 1) vector of control parameters in order (k,a,s,g)
# 
# NOTE: to avoid negative values for X(t) the condition (2*k*a)/(s^2) >= 0
#       must be satisfied!
#
simulate.cir <- function(T,delta,M,ctrl.param)
{
	Nobs <- (M/delta)*T
	X <- rep(0,Nobs)
	w <- rnorm(Nobs)
	k <- ctrl.param[1]
	a <- ctrl.param[2]
	s <- ctrl.param[3]
	g <- ctrl.param[4]
	q <- (2*k*a)/(s^2) - 1
	X[1] <- rgamma(1,shape=q+1,scale=(s^2)/(2*k))
	for (i in seq(2,Nobs,1)) 
	{
		mu <- k*(a - X[i-1])
		sigma <- s*(X[i-1]^g)
		X[i] <- X[i-1] + mu*delta + sigma*sqrt(delta)*w[i]
	}
	if (M == 1) Sobs <- seq(1,Nobs,1)
	if (M > 1)  Sobs <- seq(1,Nobs,M)
	return(X[Sobs])
}

#
# Simulate a stochastic volatility model
#
simulate.sv <- function(T,ctrl.param,ival,df)
{
	Nobs <- T
	X <- rep(0,Nobs)
	v <- rep(0,Nobs)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))
	eta <- rnorm(Nobs)
	mu <- ctrl.param[1]
	alpha <- ctrl.param[2]
	theta <- ctrl.param[3]
	sigma <- ctrl.param[4]
	w <- w*sigma
	v[1] <- alpha/(1-theta)
	X[1] <- ival[1]
	for (i in seq(2,Nobs,1)) 
	{
		v[i] <- alpha + theta*v[i-1] + w[i]
		X[i] <- mu + sqrt(exp(v[i]))*eta[i]
	}
	list(level=X,vol=exp(v))
}

#
# Quick simulation of "realized" volatility series
#
simulate.rv <- function(T,delta,mu,sigma,pj,k)
{
	Nobs <- T
	P <- matrix(0,Nobs+k,4)
	v <- rep(0,Nobs+k)
	for (i in seq(1,Nobs+k,1)) 
	{
		ui <- runif(1)
		if (ui > pj) { ri <- rnorm(1/delta,mean=(mu-0.5*sigma[1]^2)*delta,sd=sigma[1]*sqrt(delta)) }
		else { ri <- rnorm(1/delta,mean=(mu-0.5*sigma[2]^2)*delta,sd=sigma[2]*sqrt(delta)) }
		v[i] <- sum((ri-mean(ri))^2)
		pi <- exp(cumsum(ri))
		P[i,] <- c(pi[1],max(pi),min(pi),pi[1/delta])
	}	
	v <- filter(v,rep(1,k),method="convolution",sides=1)
	v <- v[seq(k+1,Nobs+k,1)]
	P <- P[seq(k,Nobs+k,1),]
	X <- diff(log(P[,4]))
	list(level=X,vol=v,price=P)
}

#
# Combined simulation of GARCH, MS-GARCH, ST-GARCH, TV-GARCH and SV models to save computing time
#
# NOTE: a separate set of control parameters must be supplied for each model
#
simulate.combined.garch <- function(T,garch.param,ms.garch.param,transm,st.garch.param,g,tv.garch.param,sv.param,ival,df)
{
	# Preliminaries
	Nobs <- T
	X <- matrix(0,Nobs,5)
	v <- matrix(0,Nobs,5)
	S <- matrix(0,Nobs,2)
	if (df == 0) w <- rnorm(Nobs)
	else w <- rt(Nobs,df)/sqrt(df/(df-2))

	# GARCH control parameters
	mu <- garch.param[1]
	omega <- garch.param[2]
	alpha <- garch.param[3]
	theta <- garch.param[4]
	v[1,1] <- omega/(1-alpha-theta)
	X[1,1] <- ival[1]

	# MS GARCH control parameters
	mu1_ms <- ms.garch.param[1,1]
	omega1_ms <- ms.garch.param[2,1]
	alpha1_ms <- ms.garch.param[3,1]
	theta1_ms <- ms.garch.param[4,1]
	mu2_ms <- ms.garch.param[1,2]
	omega2_ms <- ms.garch.param[2,2]
	alpha2_ms <- ms.garch.param[3,2]
	theta2_ms <- ms.garch.param[4,2]	
	# Get elements of transition matrix
	p11 <- transm[1,1]
	p12 <- transm[1,2]
	p21 <- transm[2,1]
	p22 <- transm[2,2]
	# Check conformity of transition probabilities
	if ((p11+p12) != 1. | (p21+p22) != 1.) { stop("Transition probabilities do not sum to 1!") }
	# Get initial stationary distribution
	p1 <- p21/(p21+p12)
	p2 <- p12/(p21+p12)
	# Get initial state
	U0 <- runif(1)
	if (U0 <= p1) { S0 = 1. }
	else { S0 = 2. }
	# Get initial values
	X[1,2] <- ival
	if (S0 == 1.) { v[1,2] <- omega1_ms/(1-alpha1_ms-theta1_ms) }
	if (S0 == 2.) { v[1,2] <- omega2_ms/(1-alpha2_ms-theta2_ms) }
	S[1,1] <- S0

	# ST GARCH control parameters
	mu1_st <- st.garch.param[1,1]
	omega1_st <- st.garch.param[2,1]
	alpha1_st <- st.garch.param[3,1]
	theta1_st <- st.garch.param[4,1]
	mu2_st <- st.garch.param[1,2]
	omega2_st <- st.garch.param[2,2]
	alpha2_st <- st.garch.param[3,2]
	theta2_st <- st.garch.param[4,2]	
	# Get initial values
	X[1,3] <- ival
	L0 <- 1/(1+exp(-g*X[1,3]))
	S[1,2] <- L0
	v[1,3] <- (1-L0)*(omega1_st/(1-alpha1_st-theta1_st)) + L0*(omega2_st/(1-alpha2_st-theta2_st))

	# TV GARCH control parameters
	mu_tv <- tv.garch.param[,1]
	omega_tv <- tv.garch.param[,2]
	alpha_tv <- tv.garch.param[,3]
	theta_tv <- tv.garch.param[,4]
	# Get initial values
	v[1,4] <- omega_tv[1]/(1-alpha_tv[1]-theta_tv[1])
	X[1,4] <- ival[1]

	# SV control parameters
	eta <- rnorm(Nobs)
	mu_sv <- sv.param[1]
	alpha_sv <- sv.param[2]
	theta_sv <- sv.param[3]
	sigma_sv <- sv.param[4]
	w_sv <- w*sigma_sv
	v[1,5] <- alpha_sv/(1-theta_sv)
	X[1,5] <- ival[1]

	# Here is the main loop
	for (i in seq(2,Nobs,1)) 
	{
		# GARCH simulation
		v[i,1] <- omega + alpha*v[i-1,1] + theta*(X[i-1,1]^2)
		X[i,1] <- mu + sqrt(v[i,1])*w[i]

		# MS GARCH simulation
		U1 <- runif(1)
		C1 <- cumsum(transm[S0,])<=U1
		S1 <- sum(C1)+1
		if (S1 == 1.) 
		{ 
			v[i,2] <- omega1_ms + alpha1_ms*v[i-1,2] + theta1_ms*(X[i-1,2]^2) 
			X[i,2] <- mu1_ms + sqrt(v[i,2])*w[i]
		}
		if (S1 == 2.) 
		{ 
			v[i,2] <- omega2_ms + alpha2_ms*v[i-1,2] + theta2_ms*(X[i-1,2]^2) 
			X[i,2] <- mu2_ms + sqrt(v[i,2])*w[i]
		}
		S0 <- S1
		S[i,1] <- S0

		# ST GARCH simulation
		Li <- 1/(1+exp(-g*X[i-1,3]))
		h1 <- omega1_st + alpha1_st*v[i-1,3] + theta1_st*(X[i-1,3]^2)
		h2 <- omega2_st + alpha2_st*v[i-1,3] + theta2_st*(X[i-1,3]^2)
		v[i,3] <- (1-Li)*h1 + Li*h2
		X[i,3] <- (1-Li)*mu1_st + Li*mu2_st + sqrt(v[i,3])*w[i]	
		S[i,2] <- Li				

		# TV GARCH simulation
		v[i,4] <- omega_tv[i] + alpha_tv[i]*v[i-1,4] + theta_tv[i]*(X[i-1,4]^2)
		X[i,4] <- mu_tv[i] + sqrt(v[i,4])*w[i]

		# SV simulation
		v[i,5] <- alpha_sv + theta_sv*v[i-1,5] + w_sv[i]
		X[i,5] <- mu_sv + sqrt(exp(v[i,5]))*eta[i]

	}
	# Convert log volatility of SV model
	v[,5] <- exp(v[,5])
	# Return
	list(level=X,vol=v,state=S)
}

#
# GARCH(1,1) filter
#
filter.garch <- function(x,garch.param)
{
	mu <- garch.param[1]
	omega <- garch.param[2]
	alpha <- garch.param[3]
	theta <- garch.param[4]
	e2 <- (x-mu)^2
	h0 <- omega/(1-alpha-theta) 
	h2 <- filter(omega+theta*lag(e2),alpha,method="recursive",init=h0)
}

#
# Quick estimation of a GARCH(1,1) model via ARMA approximation
#
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
	list(est=garch.param,vol=h2,volf=c(hf,0.55*hf))
}

#
# Function to optimize exponent in computing conditional correlation
# based on the product of the two NoVaS series...
#
optimize.cor.ep <- function(ep,wbx,wby,wwx,wwy,threshold)
{
	# Length of weight series
	px <- length(wbx)-1
	py <- length(wby)-1

	# Align wx and wy - THIS IS DONE IN MAIN PROGRAM
	#if (px > py) { wwy <- wwy[seq(px-py+1,length(wwy),1)] }
	#if (px < py) { wwx <- wwx[seq(py-px+1,length(wwx),1)] }

	# Compute the conditional correlation...
	#
	# Adjust the weights
	if (px >= py) { bb <- (wbx[1:(py+1)]*wby)^ep }
	if (px < py) { bb <- (wby[1:(px+1)]*wbx)^ep }
	# Apply threshold
	bb <- bb[bb>threshold]
	# Sum-up to one as in univariate NoVaS
	bb <- bb/sum(bb)
	# OK!
	rxy <- filter(wwx*wwy,bb,method="conv",sides=1)
	# Return
	return((max(abs(rxy),na.rm=T)-0.999)^2)
}

#
# Another function to do what is done above but with 
# a freshly selected threshold...
#
optimize.cor.theta <- function(theta,p,wwx,wwy,threshold)
{
	s <- seq(0,p)
	g <- exp(-theta*s)/sum(exp(-theta*s))
	# Trim the coefficients according to selected threshold
	b <- g[g > threshold]
	# Again, the weights sum up to one - but careful in applications...
	b <- b/sum(b)
	rxy <- filter(wwx*wwy,b,method="conv",sides=1)		
	# Return
	return((max(abs(rxy),na.rm=T)-0.999)^2)
}

#
# Function to compute recursive values for Cholesky-based
# modeling of conditional variances and correlation...
#
bivariate.cholesky.filter <- function(a1,a2,tp,theta)
{
	# Preliminaries
	if (length(a1) != length(a2)) { stop("Length of vectors does not match") }
	T <- length(a1)
	g10 <- var(a1)
	c00 <- cov(a1,a2)
	q00 <- c00/g10
	g20 <- var(a2) - (c00^2)/g10
	if (g20 < 0) { g20 <- exp(g20) }
	
	# Initializations
	g11 <- rep(0,T)
	g11[1] <- g10
	g22 <- rep(0,T)
	g22[1] <- g20
	q21 <- rep(0,T)
	q21[1] <- q00 
	sT1 <- seq(1,T-1,1)
	
	# Compute the components...
	if (tp==0)
	{
		th1 <- theta[1]^2
		ss  <- sum(exp(theta[2:3]))+1
		th2 <- exp(theta[2])/ss
		th3 <- exp(theta[3])/ss
		x11 <- th1 + th2*(c(0,a1[sT1])^2)
		g11 <- filter(x11,th3,method="recursive",init=g10)
		x21 <- theta[4] + theta[5]*(c(0,a2[sT1])^2)
		q21 <- filter(x21,theta[6],method="recursive",init=q00)
		th7 <- theta[7]^2
		ss  <- sum(exp(theta[8:11]))+1
		th8 <- exp(theta[8])/ss
		th9 <- exp(theta[9])/ss
		th10<- exp(theta[10])/ss
		th11<- exp(theta[11])/ss
		zz21<- a2[sT1] - q21[sT1]*a1[sT1]
		x22 <- th7 + th8*(c(0,a1[sT1])^2) + th9*(c(0,zz21^2)) + th10*c(0,g11[sT1])
		g22 <- filter(x22,th11,method="recursive",init=g20)
		th  <- c(th1,th2,th3,theta[4:6],th7,th8,th9,th10,th11)
	}
	if (tp==1)
	{
		th1 <- theta[1]^2
		ss  <- sum(exp(theta[2:3]))+1
		th2 <- exp(theta[2])/ss
		th3 <- exp(theta[3])/ss
		x11 <- th1 + th2*(c(0,a1[sT1])^2)
		g11 <- filter(x11,th3,method="recursive",init=g10)
		x21 <- theta[4] + theta[5]*(c(0,a1[sT1])^2) + theta[6]*(c(0,a2[sT1]^2))
		q21 <- filter(x21,theta[7],method="recursive",init=q00)
		th8 <- theta[8]^2
		ss  <- sum(exp(theta[9:10]))+1
		th9 <- exp(theta[9])/ss
		th10<- exp(theta[10])/ss
		zz21<- a2[sT1] - q21[sT1]*a1[sT1]
		x22 <- th8 + th9*(c(0,zz21^2)) + th10*c(0,g11[sT1])
		g22 <- filter(x22,th10,method="recursive",init=g20)
		th  <- c(th1,th2,th3,theta[4:7],th8,th9,th10)
	}

	# OK, compute the likelihood function and the correlation
	rho <- (q21*sqrt(g11))/sqrt(g22 + (q21^2)*g11)
	b1 <- a1
	b2 <- a2 - q21*a1
	logl<- -0.5*( log(g11) + log(g22) + ((b1^2)/g11) + ((b2^2)/g22) )
	# Return
	return(list(var1=g11,var2=(g22 + (q21^2)*g11),corr=rho,loglik=logl,parms=th))
} 

#
# Extract the likelihood component from the above function
#
bivariate.choleskly.logl <- function(param,innov1,innov2,tp)
{
	# Just call the filtering function and return appropriatelly
	out <- bivariate.cholesky.filter(innov1,innov2,tp,param)
	logl <- out$loglik
	# Remember, sum and negate!
	return(sum(-logl))
	#return(as.vector(logl))
}

#
# Filter for GARCH(1,1) model with normal or t-innovations
#
garch11.filter <- function(innov,theta,lt=0,dof=3)
{
	T  <- length(innov)
	g0 <- var(innov)
	g1 <- rep(0,T)
	# Transform the coefficients to ensure stationarity
	# and positive variance
	th1 <- theta[1]^2
	ss  <- sum(exp(theta[2:3]))+1
	th2 <- exp(theta[2])/ss
	th3 <- exp(theta[3])/ss
	x1 <- th1 + th3*(c(0,innov[seq(1,T-1,1)])^2)
	g1 <- filter(x1,th2,method="recursive",init=g0)
	if (lt==0)	{ logl<- -0.5*( log(g1) + ((innov^2)/g1) ) }
	else if (lt==1)
	{
		logl <- -0.5*(dof+1)*log(1+(innov^2)/((dof-2)*g1)) - 0.5*log(g1)	
	}
	# Return
	return(list(garchvar=g1,loglik=logl,parms=c(th1,th2,th3)))
}

#
# The likelihood component of the above function
#
garch11.logl <- function(theta,lt,dof,innov)
{
	out <- garch11.filter(innov,theta,lt,dof)
	logl<- out$loglik
	return(sum(-logl))
}

#
# EWMA modeling for variances and correlations
#
bivariate.ewma.filter <- function(a1,a2,lambda)
{
	T  <- length(a1)
	sT1<- seq(1,T-1,1)
	S0 <- cov(cbind(a1,a2))
	x1 <- (1-lambda)*c(0,a1[sT1]^2)
	g1 <- filter(x1,lambda,method="recursive",init=S0[1,1])
	x21<- (1-lambda)*c(0,a1[sT1])*c(0,a2[sT1])
	g21<- filter(x21,lambda,method="recursive",init=S0[2,1])
	x2 <- (1-lambda)*c(0,a2[sT1]^2) 
	g2 <- filter(x2,lambda,method="recursive",init=S0[2,2])
	r12<- g21/sqrt(g1*g2)
	return(list(var1=g1,var2=g2,cov21=g21,cor21=r12))	
}

#
# Compute the likelihood component of the above function
#
bivariate.ewma.logl <- function(lambda,a1,a2)
{
	T <- length(a1)
	Q <- bivariate.ewma.filter(a1,a2,lambda)
	g1 <- Q$var1
	g2 <- Q$var2
	g21<- Q$cov21
	logl <- rep(0,T)
	for (i in seq(1,T,1))
	{
		Si <- matrix(0,nrow=2,ncol=2)
		Si[1,1] <- g1[i]
		Si[2,1] <- g21[i]
		Si[1,2] <- g21[i]
		Si[2,2] <- g2[i]
		xi <- matrix(c(a1[i],a2[i]),nrow=1,ncol=2)
		logl <- -0.5*log(det(Si)) - 0.5*xi%*%solve(Si)%*%t(xi)
	}
	return(sum(-logl))
}

#
# DCC modeling for variances and correlations
#
bivariate.dcc.filter <- function(innov1,innov2,theta,...)
{
	T  <- length(innov1)
	sT1<- seq(1,T-1,1)
	# Individual GARCH models
	f1 <- garch11.filter(innov1,theta[1:3],...)
	f2 <- garch11.filter(innov2,theta[4:6],...)
	# OK, standardize the innovations and form DCC filter
	h1 <- f1$garchvar
	h2 <- f2$garchvar
	a1 <- innov1/sqrt(h1)
	a2 <- innov2/sqrt(h2)
	ss <- sum(exp(theta[7:8]))+1
	l1 <- exp(theta[7])/ss
	l2 <- exp(theta[8])/ss
	S0 <- cov(cbind(a1,a2))
	x1 <- (1-l1-l2)*S0[1,1] + l1*c(0,a1[sT1]^2)
	g1 <- filter(x1,l2,method="recursive",init=S0[1,1])
	x21<- (1-l1-l2)*S0[2,1] + l1*c(0,a1[sT1])*c(0,a2[sT1])
	g21<- filter(x21,l2,method="recursive",init=S0[2,1])
	x2 <- (1-l1-l2)*S0[2,2] + l1*c(0,a2[sT1]^2) 
	g2 <- filter(x2,l2,method="recursive",init=S0[2,2])
	return(list(garchvar1=h1,garchvar2=h2,var1=g1,var2=g2,cov21=g21,i1=a1,i2=a2))	
}

#
# Compute the likelihood component of the above function
#
bivariate.dcc.logl <- function(theta,innov1,innov2,...)
{
	T <- length(innov1)
	Q <- bivariate.dcc.filter(innov1,innov2,theta,...)
	a1 <- Q$i1
	a2 <- Q$i2
	g1 <- Q$var1
	g2 <- Q$var2
	g21<- Q$cov21
	logl <- rep(0,T)
	for (i in seq(1,T,1))
	{
		Si <- matrix(0,nrow=2,ncol=2)
		Si[1,1] <- g1[i]
		Si[2,1] <- g21[i]
		Si[1,2] <- g21[i]
		Si[2,2] <- g2[i]
		Ji <- matrix(0,nrow=2,ncol=2)
		Ji[1,1] <- 1/Si[1,1]
		Ji[2,2] <- 1/Si[2,2]
		Ri <- Ji%*%Si%*%Ji
		xi <- matrix(c(a1[i],a2[i]),nrow=1,ncol=2)
		logl <- -0.5*log(det(Ri)) - 0.5*xi%*%solve(Ri)%*%t(xi)
	}
	return(sum(-logl))
}

#
# Modeling and order selection for the NoVaS correlations based on the 
# joint distribution of the product wx*wy and maximum likelihood estimation
# of the optimizing parameter...
#

# The error function - part of the density of the product, evaluated at single values
err.fun <- function(x,z,rho)
{
	A <- 1/(2*(1-rho^2))
	B <- x^2 - 2*rho*z + (z/x)^2
	C <- (1/(2*pi*sqrt(1-rho^2)))*exp(-A*B)/x
	return(C)
}

# Supportive for uniform distribution
uni.fun <- function(x) 
{ 
	C <- 1/abs(x) 
	return(C)
}

# Function for evaluating the density based on an exponentially weighted rolling mean
Fewmean.fun <- function(x,y,K,dtype,rtype,lambda)
{
	N    <- NROW(x)
	z    <- scale(x)*scale(y)
	rbar <- mean(z)
	sl   <- lambda^seq(0,K-1,1)
	Fmat <- matrix(0,N-K+rtype,2)
	sumw <- (1-lambda)*sum(sl)
	if (dtype=="normal")
	{
		for (i in seq(K+(1-rtype),N,1))
		{
			# OK, compute the correlation! NOTE that you reverse the order of the weights!!!
			R <- (1-lambda)*sum(rev(sl)*z[(i-K+rtype):(i+rtype-1)]) + (1-sumw)*rbar
			# R <- sum(dnorm((z-z[i+rtype-1])/lambda)*z[i])/N 
			# Squeeze rho? Use inverse Fisher transformation
			if (abs(R) > 1) { R <- (exp(2*R)-1)/(exp(2*R)+1) }
			#R <- (exp(2*R)-1)/(exp(2*R)+1)
			f1 <- integrate(err.fun,lower=1e-5,upper=Inf,stop.on.error=FALSE,z=z[i],rho=R)
			f2 <- integrate(err.fun,lower=-Inf,upper=1e-5,stop.on.error=FALSE,z=z[i],rho=R)
			Fmat[i-K+rtype,1] <- f1$value-f2$value
			Fmat[i-K+rtype,2] <- R
		}
	}
	if (dtype=="uniform")
	{
		for (i in seq(K+(1-rtype),N,1))
		{
			# OK, compute the correlation! NOTE that you reverse the order of the weights!!!
			R <- (1-lambda)*sum(rev(sl)*z[(i-K+rtype):(i+rtype-1)]) + (1-sumw)*rbar
			# R <- sum(dnorm((z-z[i+rtype-1])/lambda)*z[i])/N
			# Squeeze rho? Use inverse Fisher transformation
			if (abs(R) > 1) { R <- (exp(2*R)-1)/(exp(2*R)+1) }			
			#R <- (exp(2*R)-1)/(exp(2*R)+1)
			f1 <- 1/sqrt((1-R)*(1+R))
			limit <- sqrt(3)*(1+R)
			f2 <- integrate(uni.fun,lower=1e-5,upper=limit,stop.on.error=FALSE)
			Fmat[i-K+rtype,1] <- f1*(2*f2$value)
			Fmat[i-K+rtype,2] <- R
		}
	}	
	return(list(FF=Fmat,logl=sum(log(Fmat[,1]))-K,corr=Fmat[,2]))
}

# A function for evaluating the density based on the NoVaS weights
Fnovas.fun <- function(x,y,K,dtype,rtype,ep)
{
	N    <- NROW(x)
	z    <- scale(x)*scale(y)
	rbar <- mean(z)
	s    <- seq(0,K)
	b    <- exp(-ep*s)
	M    <- length(b)
	sumw <- sum(b/M)
	# Done, initialize and compute...
	Fmat <- matrix(0,N-M+rtype,2)
	#
	if (dtype=="normal")
	{
		for (i in seq(M+(1-rtype),N,1))
		{
			# OK, compute the correlation! NOTE that you reverse the order of the weights!!!
			R <- mean(rev(b)*z[(i-M+rtype):(i+rtype-1)]) + (1-sumw)*rbar
			# Squeeze rho? Use inverse Fisher transformation
			if (abs(R) > 1) { R <- (exp(2*R)-1)/(exp(2*R)+1) }			
			#R <- (exp(2*R)-1)/(exp(2*R)+1)
			f1 <- integrate(err.fun,lower=1e-5,upper=Inf,stop.on.error=FALSE,z=z[i],rho=R)
			f2 <- integrate(err.fun,lower=-Inf,upper=1e-5,stop.on.error=FALSE,z=z[i],rho=R)
			Fmat[i-M+rtype,1] <- f1$value-f2$value
			Fmat[i-M+rtype,2] <- R
		}
	}
	if (dtype=="uniform")
	{
		for (i in seq(M+(1-rtype),N,1))
		{
			# OK, compute the correlation! NOTE that you reverse the order of the weights!!!
			R <- mean(rev(b)*z[(i-M+rtype):(i+rtype-1)]) + (1-sumw)*rbar
			# Squeeze rho? Use inverse Fisher transformation
			if (abs(R) > 1) { R <- (exp(2*R)-1)/(exp(2*R)+1) }			
			#R <- (exp(2*R)-1)/(exp(2*R)+1)
			f1 <- 1/sqrt((1-R)*(1+R))
			limit <- sqrt(3)*(1+R)
			f2 <- integrate(uni.fun,lower=1e-5,upper=limit,stop.on.error=FALSE)
			Fmat[i-M+rtype,1] <- f1*(2*f2$value)
			Fmat[i-M+rtype,2] <- R
		}
	}	
	return(list(FF=Fmat,logl=sum(log(Fmat[,1]))-K,corr=Fmat[,2]))
}

# A function for evaluating the density when the matching is uniform

# Suportive function for optimizing the parameters for computing the correlations
optimize.Fnovas.fun <- function(param,.wx,.wy,.K,.dtype,.rtype,.Ftype)
{
	if (.Ftype=="ewmean")  { logl <- Fewmean.fun(.wx,.wy,.K,.dtype,.rtype,param)$logl }
	if (.Ftype=="novas") { logl <- Fnovas.fun(.wx,.wy,.K,.dtype,.rtype,param)$logl }
	return(-logl) 
}	

#
# NoVaS correlation estimation...
#
novas.correlation <- function(x,y,var.type,opt.type,target,mfact,whichd,whichr,whichF,...)
{
	xy <- as.matrix(cbind(x,y))
	x  <- xy[,1]
	y  <- xy[,2]
	
	# Split the parameters
	vtpx <- var.type[1]
	vtpy <- var.type[2]
	
	optx <- opt.type[1]
	opty <- opt.type[2]
	
	tgtx <- target[1]
	tgty <- target[2]
	
	# Re-estimate and extract what is needed
	outx <- novas.estimate(x,var.type=vtpx,opt.type=optx,target=tgtx,...)
	px <- outx$parameters[3]
	wx <- outx$zret
	vx <- outx$vol
	bx <- outx$weights
	
	outy <- novas.estimate(y,var.type=vtpy,opt.type=opty,target=tgty,...)
	py <- outy$parameters[3]
	wy <- outy$zret
	vy <- outy$vol
	by <- outy$weights
	
	# Alignment
	ww <- novas.align(wx,wy)
	wx <- ww[,1]
	wy <- ww[,2]
	
	bb <- novas.align(bx,by)
	bx <- bb[,1]
	by <- bb[,2]
	
	# Now for the correlation
	K <- ceiling(mfact*max(c(px,py)))
	outr <- optimize(optimize.Fnovas.fun,interval=c(0,1),.wx=wx,.wy=wy,.K=K,.dtype=whichd,.rtype=whichr,.Ftype=whichF)
	if (whichF=="ewmean") { outf <- Fewmean.fun(wx,wy,K,whichd,whichr,outr$minimum) }
	if (whichF=="novas") { outf <- Fnovas.fun(wx,wy,K,whichd,whichr,outr$minimum) }
	rxy <- outf$corr
	lgl <- outf$logl
	oep <- outr$minimum
	
	# Done, output everything...
	return(list(lambda=oep,zretx=wx,zrety=wy,wtx=bx,wty=by,volx=vx,voly=vy,corr=rxy,logl=lgl))
}

#
# Supportive functions
#
novas.vec <- function(x, byrow = FALSE) 
{
    if (is.vector(x)) 
        return(x)
    if (byrow) 
        x <- t(x)
    d <- ncol(x)
    vecx <- vector()
    for (j in 1:d) vecx <- c(vecx, x[, j])
    return(vecx)
}

novas.vech <- function(x) 
{
    if (is.vector(x)) {
        if (length(x) == 1) 
            return(x)
        else stop("vech undefined for vectors")
    }
    else if (is.matrix(x)) {
        d <- ncol(x)
        if (d != nrow(x)) 
            stop("vech only defined for square matrices")
        vechx <- vector()
        for (j in 1:d) vechx <- c(vechx, x[j:d, j])
        return(vechx)
    }
}

novas.invvec <- function(x, ncol, nrow, byrow = FALSE) 
{
    if (length(x) == 1) 
        return(x)
    d <- sqrt(length(x))
    if (missing(ncol) | missing(nrow)) {
        ncol <- d
        nrow <- d
        if (round(d) != d) 
            stop("Need to specify nrow and ncol for non-square matrices")
    }
    invvecx <- matrix(0, nrow = nrow, ncol = ncol)
    if (byrow) 
        for (j in 1:nrow) invvecx[j, ] <- x[c(1:ncol) + (j - 
            1) * ncol]
    else for (j in 1:ncol) invvecx[, j] <- x[c(1:nrow) + (j - 
        1) * nrow]
    return(invvecx)
}

novas.invvech <- function(x) 
{
    if (length(x) == 1) 
        return(x)
    d <- (-1 + sqrt(8 * length(x) + 1))/2
    if (round(d) != d) 
        stop("Number of elements in x will not form a square matrix.")
    invvechx <- matrix(0, nr = d, nc = d)
    for (j in 1:d) invvechx[j:d, j] <- x[1:(d - j + 1) + (j - 
        1) * (d - 1/2 * (j - 2))]
    invvechx <- invvechx + t(invvechx) - diag(diag(invvechx))
    return(invvechx)
}

#
# Rolling realized correlation via standardized returns...
#
novas.rolling.correlation <- function(zret,K)
{
	# Observations and variables
	Nobs <- NROW(zret)
	Nvar <- NCOL(zret)
	
	# Storage for the unique elements of the correlations
	corrs <- matrix(0,nrow=Nobs-K+1,ncol=0.5*Nvar*(Nvar-1))
	
	# OK, do the loop
	for (i in seq(1,Nobs-K+1,1))
	{
		# Rolling mean computation
		covi <- 0
		for (j in seq(1,K,1))
		{
			covi <- covi + matrix(zret[i-1+j,],Nvar,1)%*%matrix(zret[i-1+j,],1,Nvar)
		}
		covi <- covi/K
		corri <- cov2cor(covi)
		# Now save...
		#
		# Extract the unique elements...
		ucorri <- matrix(lower.tri(corri)*corri,Nvar*Nvar,1)
		ucorri <- subset(ucorri,abs(ucorri)>0)
		corrs[i,] <- ucorri
	}
	
	# Done, return it
	return(corrs)
}

#
# A function to simulate a "nice" multivariate GARCH type volatility model
#
simulate.mgarch.rcov <- function(.nobs,.nvar,.alpha,.beta,.gamma,.rho,.asym=0,.M)
{
	# Construct the unconditional covariance
	Sbar <- matrix(.rho,nrow=.nvar,ncol=.nvar)
	diag(Sbar) <- rep(1,.nvar)

	# Initializations
	St <- matrix(0,nrow=.nobs+100,ncol=.nvar^2)
	St[1,] <- novas.vec(Sbar)
	Qt <- matrix(0,nrow=.nobs+100,ncol=.nvar^2)
	Qt[1,] <- novas.vec(cov2cor(Sbar))
	Rt <- matrix(0,nrow=.nobs+100,ncol=.nvar)
	RC <- matrix(0,nrow=.nobs+100,ncol=.nvar^2)
	RR <- matrix(0,nrow=.nobs+100,ncol=.nvar^2)

	# Do the loop
	for (i in seq(2,.nobs+100,1))
	{
		# Get the intraday-innovations
		xmt <- matrix(rnorm(.M*.nvar,sd=(1/sqrt(.M))),nrow=.M,ncol=.nvar)

		# Get the daily innovations
		et <- apply(xmt,2,sum)

		# Compute the covariance matrix
		S0 <- novas.invvec(St[i-1,])
		if (.asym==0) { S1 <- .alpha*Sbar + .beta*S0 + .gamma*matrix(Rt[i-1,],.nvar,1)%*%matrix(Rt[i-1,],1,.nvar) }
		else
		{
			IA <- prod(as.double(Rt[i-1,] < 0))
			CR <- matrix(Rt[i-1,],.nvar,1)%*%matrix(Rt[i-1,],1,.nvar)
			S1 <- .alpha*Sbar + .beta*S0 + .gamma*CR + .asym*IA*CR
		}

		# Save the new covariance and correlation
		St[i,] <- novas.vec(S1)
		Qt[i,] <- novas.vec(cov2cor(S1))

		# Decomponse and get the return
		d1 <- eigen(S1,symmetric=TRUE)
		S12<- (d1$vectors)%*%diag(sqrt(d1$values))
		Rt[i,] <- S12%*%et

		# Get and save the realized covariance and correlation
		Sxmt <- t(xmt)%*%xmt
		RCV  <- S1%*%Sxmt
		RC[i,] <- novas.vec(RCV)
		RR[i,] <- novas.vec(cov2cor(RCV))
		# Squeeze the realized measure
	}

	# Extract the useful observations...
	St <- St[101:NROW(St),]
	Qt <- Qt[101:NROW(Qt),]
	Rt <- Rt[101:NROW(Rt),]
	RC <- RC[101:NROW(RC),]
	RR <- RR[101:NROW(RR),]	
	
	# Done, return...
	return(list(covariance=St,correlation=Qt,returns=Rt,realized.cov=RC,realized.cor=RR))
}

#
# Simple volatility forecasting/smoothing function
#
covariance.forecast <- function(Rt,.beta,.gamma,.k)
{
	# Preliminaries
	Nobs <- NROW(Rt)
	Nvar <- NCOL(Rt)
	Sbar <- cov(Rt)
	bg   <- (.beta+.gamma)
	
	# Initializations	
	St <- matrix(0,nrow=Nobs,ncol=Nvar^2)
	St[1,] <- novas.vec(Sbar)
	Qt <- matrix(0,nrow=Nobs,ncol=Nvar^2)
	Qt[1,] <- novas.vec(cov2cor(Sbar))
	
	# The loop
	for (i in seq(2,Nobs,1))
	{
		S0 <- novas.invvec(St[i-1,])
		S1 <- (1-.k)*Sbar + (.beta/bg)*.k*S0 + (.gamma/bg)*.k*matrix(Rt[i-1,],Nvar,1)%*%matrix(Rt[i-1,],1,Nvar)
		St[i,] <- novas.vec(S1)
		Qt[i,] <- novas.vec(cov2cor(S1))
	}
	
	# Return
	return(list(covariance=St,correlation=Qt))
}

#
# Repairing an indefinite correlation matrix
#
repair.correlation <- function(Q)
{
	# Do spectral decomposision
	E <- eigen(Q,symmetric = TRUE)   
	V <- E$vectors
	L <- E$values

	# Replace negative eigenvalues by zero
	L <- pmax(L,0)

	# Reconstruct correlation matrix
	QQ <- V%*%diag(L)%*%t(V)

	# Rescale correlation matrix
	S <- 1/sqrt(diag(QQ))
	SS<- outer(S,S)
	R <- QQ*SS
	
	# Done
	return(R)
}

#
# Optimal exponent selection for EWMA based on variance ratio
#
select.lambda <- function(x,prop)
{
	.rss <- function(lambda,.x)
	{
		out <- novas.ewma(.x,lambda)$fit
		varr<- sd(out)/sd(x)
		return(sum((varr-prop)^2))
	}
	out <- optimize(.rss,interval=c(0,1),.x=x)
	lambda <- out$minimum
	return(lambda)
}

#
# Fast, and correct, EWMA implementation (with forecasting)
#
novas.ewma <- function(x,lambda,iszoo=FALSE)
{
	Nobs <- NROW(x)
	startup <- floor(log(Nobs))
	out  <- as.double(filter((1-lambda)*as.double(x),lambda,method="recursive",init=0))
	fitt <- c(mean(as.double(x)[1:startup]),out[1:(Nobs-1)])
	forc <- out[Nobs]
	if (iszoo==TRUE) { fitt <- zoo(fitt,order.by=index(x)) }
	return(list(fit=fitt,forecast=forc))
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
VARSTABdp<-function(xx, p, alp)
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
  V = var(xx[1:(n-1)], na.rm = T)
  
  zz <- yy/sqrt(alp*V+kk/p)    # zz_as.vector(zz)
  zz<- zz[(2*p+1):(2*p+n)]
  zz <- na.omit(zz)
  zz[1] <- sign(zz[1])
  zz
}


########## seclect p
P.OPTsimple<- function(x, alp=0,  plot = F)
{ # x is financial returns series (percentage returns)
  # search for optimum P in Simple NOVAS 
  zz <- c(1:50)
  c4 <- c(1:50)
  ccc <- c(1:50)
  
  for(i in 1:50) {
    zz[i] <- KURTOSIS(VARSTABdp(x, i, alp))
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