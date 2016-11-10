#=========================== pool(x,gap) ========================================================
# Pools x-values which are closer than gap   
pool <- function(x,gap=10^-4){
  n <- length(x)
  x2 <- sort(x)
  x3 <- c()
  w <- c()
  temp <- x2[1]
  for (i in 2:n){
    if ((x2[i]-mean(temp))<=gap){temp <- c(temp,x2[i])} else {
      if (length(temp)>1){
        x3 <- c(x3,mean(temp))
        w <- c(w,length(temp))
        temp <- x2[i]} else {
          temp <- x2[i]
          x3 <- c(x3,x2[i-1])
          w <- c(w,1)}
    }}
  if (length(temp)>0){
    x3 <- c(x3,mean(temp))
    w <- c(w,length(temp))}
  
  list(x3,w)}  

#################################
# Auxiliary function for icma #
#################################

isoMean <- function(y,w){
  
  # Pool-adjacent-violaters-algorithm for a weighted
  # mean.  
  #
  # Input:
  #   - y     : data points                                       
  #   - w     : corresonding weights                      
  
  n <- length(y)
  k <- 1:n*0
  gew <- 1:n*0
  ghat <- 1:n*0
  c <- 1
  k[c] <- 1
  gew[c] <- w[1]
  ghat[c] <- y[1]
  
  for (j in 2:n){     
    c <- c+1
    k[c] <- j
    gew[c] <- w[j]
    ghat[c] <- y[j]
    
    while (c>=2 && ghat[max(1,c-1)]>=ghat[c]){
      neu <- gew[c]+gew[c-1]
      ghat[c-1] <- ghat[c-1]+(gew[c]/neu)*(ghat[c]-ghat[c-1])
      gew[c-1] <- neu
      c <- c-1}}
  
  while (n>=1){
    for (j in k[c]:n){ghat[j] <- ghat[c]}
    n <- k[c]-1
    c <- c-1}
  ghat}


#================= starter.par1(x) ==========================================================
# start vector: fitted values for a model log(kern(x))~1+x+log(x) (corresponds to Gamma density)
starter.par1 <- function(x){
  X <- cbind(1,x,x^2)
  #X <- cbind(1,x,log(x))
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  phi <- H%*%(log(kern(x)))
  
  # it can happen that phi is convex instead of concave (e.g. when the true density is uniform)
  d <- diff(phi)/diff(x)
  if (d[2]<d[3]){phi <- -phi}
  
  eta <- phieta.par1(x,phi)
  eta[1] <- abs(eta[1])
  eta[-1] <- eta[-1] - eta[2] -10^-10
  eta}

#== Likelihood-function: ===========
Lhat.par1 <- function(dx,eta){
  n <- length(eta)
  L <- -n*eta[1]-sum((n-2:n+1)*alpha(dx,eta)[2:n])+n*exp(eta[1])*f.star(2:n,dx,eta)
  L}

#== Gradient of L: =================
an.grad.par1 <- function(dx,eta){
  n <- length(eta)                        
  dL2n <- 1:n*NA
  H0 <- n*exp(eta[1])
  dL1 <- n-H0*f.star(2:n,dx,eta)
  
  a <- alpha(dx,eta)
  b <- beta(dx,eta)
  k <- 2:n
  temp <- (dx*rev(cumsum(rev(c(0,1/eta[3:n]*(b[3:n]-b[2:(n-1)]),0)))))[k]    
  dL2n[k] <- (n-k+1)*dx[k]-H0*(-eta[k]^(-2)*(b[k]-b[k-1])+1/eta[k]*(dx[k]*b[k])+temp) 
  G <- c(dL1,dL2n[2:n])
  -G}

#== Diagonal of Hesse of L : =======
an.hess.par1 <- function(dx,eta){
  n <- length(eta)            
  d22n <- 1:n*NA
  H0 <- n*exp(eta[1])
  d21 <- -H0*f.star(2:n,dx,eta)
  
  a <- alpha(dx,eta)
  b <- beta(dx,eta)
  k <- 2:n
  temp <- (dx^2*rev(cumsum(rev(c(0,1/eta[3:n]*(b[3:n]-b[2:(n-1)]),0)))))[k]    
  d22n[k] <- -H0*(2*eta[k]^(-3)*(b[k]-b[k-1])-2*eta[k]^(-2)*(dx[k]*b[k])+dx[k]^2*b[k]/eta[k]+temp)    
  H <- c(d21,d22n[2:n])
  -diag(H)}

#=========================== kern(x) ============================================================
# kern(x) : calculates standard kernel density estimator at the observations x
kern <- function(x){
  x <- sort(x)
  n <- length(x)
  K <- function(x){1/sqrt(2*pi)*exp(-x^2/2)}
  h <- 1.06*sqrt(var(x))*length(x)^(-1/5)
  f <- NULL
  for (i in 1:n){f[i] <- sum(K((x[i]-x)/h))/(n*h)}
  f}

#================= phieta.par1(x,phi), etaphi.par1(x,eta) ===================================
# Changes between different parametrisations
phieta.par1 <- function(x,phi){ 
  n <- length(x)
  eta <- c(phi[1],diff(phi)/diff(x))
  eta}
etaphi.par1 <- function(x,eta){
  n <- length(x)
  dx <- c(0,diff(x))
  phi <- eta[1]+c(0,cumsum(dx[2:n]*eta[2:n]))
  phi}
#== divers auxiliary functions for derivatives =======
alpha <- function(dx,eta){c(0,(dx*eta)[-1])}
beta <- function(dx,eta){exp(cumsum(alpha(dx,eta)))}

#== 3rd summand in L: ==============
f.star <- function(j,dx,eta){sum(1/eta[j]*(beta(dx,eta)[j]-beta(dx,eta)[j-1]))}
