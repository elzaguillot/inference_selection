## R script
## created by eguillot december 2016

## attempt to reproduce equation 5 of paper by Keightley et al, 2007 because it was not clear why and how it would work
## (Weird reference in the paper)

require(matlab)

## fonction to fill a wright fisher matrix according to equation 1 and 2
Mjk <- function(j,k,n,s){
    q=j*1.0/(2*n)
    deltaq=(-s*q*(1.0-q))/(2*(1.0-s*q))
    return(choose(2*n,k)*(q+deltaq)^k*(1.0-q-deltaq)^(2*n-k))
}
## vectorize to use outer
Mjkv <- Vectorize(Mjk)

## function to create the matrix as function of N (2x pop size)
createM <- function(s,N)
{
    return(outer(seq(0, N), seq(0,N), Mjkv,n=N/2.0,s))
    
}

## function to computer the power of a matrix, e.g. M^3=M*M*M
## not optimized to keep numerical precision
matrixpower <- function(M,n)
{
    out=M
    for(i in 1:(n-1))        
        {out=out%*%M}
    return(out)
}

## Now test on an example
# compute equation 5 with Kemeny and Snell 
N=20
s=0.0001
myM=createM(s,N)
rowSums(myM) ## check the probability matrix is ok, all equal 1
I=eye(N-1) #identity matrix
Q=myM[2:N,2:N] # extract the submatrix
IQ=I-Q  # (i-q)
T=solve(IQ) # (i-q)^-1
Xeq=(T)[,1] # extract first column -> allele freq distribution at equilibrium


## using brute force to test: multiplying an inital freq by M for a large number of t

## What is the cumulative frequency  at equilibrium by assuming eq 4 
T1=10000 #big t
x0=rep(0,N+1);x0[2]=1 # initial vector with 1 mutant
xeq=rep(0,N+1) # x0 descrived in text before equation 3
for(t in 1:T1)
    xeq=xeq+x0%*%matrixpower(myM,t) ## eq 3+4 with big t
xeq=xeq/T1


T2=100000 #bigger t to test if stable
x0=rep(0,N+1);x0[2]=1 # initial vector with 1 mutant
xeq2=rep(0,N+1)
for(t in 1:T2)
{
    print(t)
    xeq2=xeq2+x0%*%matrixpower(myM,t) ## eq 3+4 with bg t
}
xeq2=xeq2/T2

## Compare

#png('keightley2.png')
plot(as.vector(xeq2[2:N]),type='l',col='green',xlab='nb of mutant alleles',ylab='frequency',main="SFS at equilibrium",lwd=2)
lines(as.vector(xeq2[2:N]/sum(xeq2[2:N])),col='blue',lwd=2,lty=2)
lines(as.vector(Xeq/sum(Xeq)),col='red',lwd=2)
legend('topright',c('Brute force T1','Brute force T2','Kemeny and Snell'),col=c('blue','green','red'),lty=1)
#dev.off()



createM2 <- function(s,N1,N2)
{
    return(outer(seq(0, N1), seq(0,N2), Mjkv,n=N2/2.0,s))
    
}

N2=100
M2=createM2(s,N,N2)
T1bis=100
myM2=createM(s,N2)

xeq3=(xeq)%*%(M2)%*%matrixpower(myM2,T1bis) ## eq 3+4 with big t
Xeq3=(Xeq)%*%(M2)%*%matrixpower(myM2,T1bis) ## eq 3+4 with big t
