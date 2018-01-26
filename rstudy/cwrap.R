## Rscripts that uses the C++ made function to compile likelihood
## and uses R function to run Maximum likelihood algorithms

## This scripts attemtp different ML algorithm with sfs from simulated data
## Not giving clear results as of 10.12.17 (Eguillot)

setwd('/home/eguillot/Documents/unil/dfe/')
require(Rcpp)
require(stats4)
require(nloptr)
dyn.load('moran.so',now=T)
require(pso)

## from the output:   popsize1 popsize2 nbgen  f0 alpha beta mu


popsize2=200
nbgen=100
popsize1=100
f0=0
mu=0.0000001
alpha=0.5
beta=0.5
sfs="simulation_g200.sfs"
sfs="simu4.sfs"
sfs="slimS.sfs"

#try t ooptimize only tow parameters
LL <- function(alpha,beta)
{
    a=.Call("moranS",c(100,200,200,0.01,0.01,alpha,beta),sfs)
    return(-a)
}

mle(LL, start = list(alpha=0.5,beta=0.5), method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(0.99,0.99))


LL <- function(popsize2,nbgen,popsize1,f0,mu,alpha,beta)
{
    a=.Call("moranS",c(popsize2,nbgen,popsize1,f0,mu,alpha,beta),sfs)
    return(-a)
}

mle(LL, start = list(popsize2=200,nbgen=100,popsize1=100,f0=0.1,mu=0.0001,alpha=0.5,beta=0.5), method = "L-BFGS-B", lower = c(100,100,10,0.01,0.000001,0.001,0.001), upper = c(1000,10000,1000,0.5,0.01,0.99,0.99))



#try t ooptimize only tow parameters
LL <- function(x)
{
    a=.Call("moranS",c(100,200,2000,0,0.01,x[1],x[2]),sfs)
    return(-a)
}

mle(LL, start = list(alpha=0.5,beta=0.5), method = "L-BFGS-B", lower = c(0.001,0.001), upper = c(0.99,0.99))

bobyqa(c(0.5, 0.5), LL, lower = c(0.001,0.001), upper = c(0.99,0.99))

bobby=bobyqa(c(0.5,0.5), LL,upper=c(0.99,0.99),lower=c(0.001,0.001))#, lower = c(0.001,0.001), upper = c(0.9,0.9))

LL <- function(x)
{
    a=.Call("moranS",x,sfs)
    return(-a)
}

bobyqa(c(200,100,100,0.1,0.0001,0.001,0.001), LL, lower = c(100,100,10,0.01,0.000001,0.001,0.001), upper = c(1000,10000,1000,0.5,0.01,0.001,0.001))



LL0 <- function(x)
{
    a=.Call("moran0",c(x[1],x[2],x[3],x[4],.5,.5),sfs)
    return(-a)
}

bobby=bobyqa(c(400,100,100,0.1,0.0001), LL0, lower = c(100,100,10,0.01,0.000001), upper = c(1000,10000,1000,0.5,0.01))


bobby=nloptr(c(200,100,100,0.1,0.0001), LL0, lower = c(100,100,10,0.01,0.000001), upper = c(1000,10000,1000,0.5,0.01))


psoptim(c(200,100,100,0.1,0.0001), LL0, lower = c(100,100,10,0.01,0.000001), upper = c(1000,10000,1000,0.5,0.01))

