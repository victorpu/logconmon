micma <- function(data,eps=10^-10,T1=200,T2=20,output='N',goalLL,goalTime){
  #  
  # Computation of a nonparametric maximum-likelihood estimator (MLE)
  # for a log-concave density function by minimising the second order 
  # taylor approximation to the log-likelihood-function by the 
  # iterative convex minorant algorithm
  #
  # Input:
  # - data       	: 	Consider observations X(1), X(2), ..., X(n) with unknown density
  #			function f. Suppose that we want to estimate this density under 
  #			the constraint that its logarithm is a concave function 
  #			phi: R --> [-infty,infty).
  #
  # - eps		: 	accuracy parameter in R+
  # - T1		: 	maximal number of Newton iterations (outer step) in N
  # - T2		: 	maximal number of inner step iterations in N
  # - output	:	Y: print output N: suppress output
  # 
  # Output: list containing two elements [[1]] and [[2]]:
  # - [[1]]  	: 	data frame containing two columns: x (inputed data) and
  #			f (estimated density at observation points)
  # - [[2]]  	: 	vector containing the values for the log-likelihood for every iteration
  # 
  # Auxiliary programs:
  # - auxiliary.r :	functions to calculate starting vector, reparametrisation, newton-step, 
  #			log-likelihood-function and its derivative and to display results,
  #	
  # Version: 	1.1
  # Date:		September 2004
  # Author:  	Kaspar Rufibach (Supported by Schweizer Nationalfonds)
  #
  # Source:
  # - Kaspar Rufibach (2004).
  #   Computing ML-Estimators of a log-concave density function.
  #   Preprint
  
  t1 <- proc.time()[1]
  #====== General settings ===============================
  x     <- sort(data)
  dx    <- c(0,diff(x))
  n     <- length(x)
  iter1 <- 0
  hist <- c()
  abbr <- 10
  dirder <- 2*eps
  
  #======= Start vectors ===============================
  #phi = log(2/sqrt(2*pi)) - x^2/2
  #eta = phieta.par1(x,phi)
  eta   <- starter.par1(x)
  etanew <- 1:n*0
  lik0 <- Lhat.par1(dx,eta)
  loglik <- lik0
  LL <- lik0
  time  <- proc.time() 
  
  #===============  start outer step ====================================
  #while (abs(dirder)>eps && iter1<T1){            			
  while (goalLL<=min(LL,na.rm=T) && proc.time()[1]-t1<=goalTime){
    iter1 <- iter1+1			
    b <- an.grad.par1(dx,eta)
    d <- -diag(an.hess.par1(dx,eta))
    y <- eta+b/d
    etanew[1] <- y[1]
#    etanew[2:n] <- -isoMean(-y[2:n],d[2:n])
    lambda = max(abs(d))*100
    mstar = 1
    while(abs(mstar)>1e-5){
      result <- -isoMean(c(0,-y[2:n]),c(lambda,d[2:n]))
      etanew[2:n] = result[2:n]
      mstar = result[1]
      lambda = 10*lambda
    }
    liknew <- Lhat.par1(dx,etanew)
    dirder <- as.numeric(t(b)%*%(etanew-eta))
    #---------------------- Robustification -------------------------
    iter2 <- 0
    while (liknew>loglik && iter2<T2){ 	
      iter2 <- iter2+1
      etanew <- (eta+etanew)/2
      liknew <- Lhat.par1(dx,etanew)
      dirder <- dirder/2}		
    #---------------------- Robustification -------------------------
    #---------------------- Hermite-Interpolation -------------------
    t0 <- (2-(liknew-loglik)/dirder)^(-1)
    if (t0<1){eta <- (1-t0)*eta+t0*etanew}
    else {eta <- etanew}
    #---------------------- Hermite-Interpolation -------------------
    
    loglik <- Lhat.par1(dx,eta)
    out <- data.frame("Iteration"=iter1,"LogLikelihood"=Lhat.par1(dx,eta))
    if (output=='Y'){print(out)}
    LL <- c(LL,out$LogLikelihood)
    time <- c(time,proc.time()[3])
  }
  
  # Generate output
  list(data.frame("x"=x,"f"=exp(etaphi.par1(x,eta))),LL,c(0,diff(time)))
}



