library(fdrtool)
n = 500
typ <- "normal"
random <- function(n){
  x <- sort(abs(rnorm(n,0,1)))
  x}

x <- random(n)
x <- pool(x,gap=10^-4)[[1]]
d1 <- micma(data=x,eps=10^-10,T1=200,T2=20,output='N',goalLL = 0,goalTime = 120)
plot(d1[[1]][,1],d1[[1]][,2],type = "l")
#LL = 802.1701

d2 <- icma(data=x,eps=10^-10,T1=200,T2=20,output='N',goalLL = 0,goalTime = 120)
plot(d2[[1]][,1],d2[[1]][,2],type = "l")
#LL = 805.2968
#log(6930.567)
# Initialization could be shaped-constraint density estimator.
# when apply monotone regression, ignore eta[1] at first, let continuity assign values to eta[1].

n = 50
typ <- "normal"
random <- function(n){
  x <- sort(abs(rnorm(n,0,1)))
  x}

x <- random(n)
x <- pool(x,gap=10^-4)[[1]]
d1 <- micma(data=x,eps=10^-10,T1=200,T2=20,output='Y',goalLL = 0,goalTime = 10)
plot(d1[[1]][,1],d1[[1]][,2],type = "l")
#LL = 87.5331

d2 <- icma(data=x,eps=10^-10,T1=200,T2=20,output='Y',goalLL = 0,goalTime = 10)
plot(d2[[1]][,1],d2[[1]][,2],type = "l")
#LL = 87.43504

e = ecdf(x)
d3 = grenander(e)
plot(d3) 