## Script to reproduce dfe estimation by Keitghtley and Eyre Walker 2007
## Written by Elsa Guillot elsa.guillot@unil.ch
##
## First attempt to reproduce the equations in the paper
## However the matrices are too big and the computation too long for it to work
## Useful to somewhat understand the workflow but do not use it to run


require(expm)
require(stats)

## compute the transition matrix

## equation 1
Mjk <-function(j,k,n,s){ ## this is wrong because of computation approximations
    
    q=j*1.0/(n)
    out=choose(n,k)*(q+deltaq(q,s))^k*(1-q-deltaq(q,s))^(n-k)
    return(out)
}


## equation 2
deltaq <- function(q,s){
    return(-s*q*(1-q)*1.0/(2*(1-s*q)))
}

## equation 3
f_eq3<- function(t2,n,M){
    f0=c(1,rep(0,n-1))
    myf=f0 %*% (M %^% t2)
    return (myf)
}

## equation 4
x_eq4 <- function(tt,n,M)
    {
        t2 = seq(0,tt)
        myf = sapply(t2,f_eq3,n,M)
        return(rowSums(myf))
    }

## if s>1
xs1_eq4<- function(s,n)
    {
        myx <- rep(0,n)
        myx[1]=2.0/s
        return(myx)
    }


## equation 5
u_eq5<- function(M)
    {
        n=dim(M)[1]
        Q=M[1:(n-1),1:(n-1)]
        P=solve(diag(n-1)-Q)
        heatmap(M)
        return(P[,1])
    }

## between 5 and 6 described in paragraph
ws_eq<- function(M1,M2)
    {
        N1= dim(M1)[1]
        N2= dim(M2)[1]
        TrN1N2 <- matrix(0,ncol=N2,nrow=N1-1)
        TrN1N2[1:(N-1),1:(N1-1)]=matrix(1,nrow=N1-1,ncol=N1-1)
        return((u_eq5(M1)) %*% TrN1N2 %*% M2)
    }

## equation 6
vp_eq6<- function(M1,M2,t2,s)
    {
        N1= dim(M1)[1]
        N2= dim(M2)[1]
        if(s<=1){
            return( N1*ws_eq(M1,M2)+N2*x_eq4(t2,N2,M2))
        }
        else{
#            print(N1*ws_eq(M1,M2))
#            print(N2*xs1_eq4(s,N2))
            return( N1*ws_eq(M1,M2)+N2*xs1_eq4(s,N2))
        }
    }

## equation 7
vs <- function(N1,N2,s,t2,f0)
{
    M01 <- createM(0,N1)
    M02 <- createM(0,N2)
    M1 <- createM(s,N1)
    M2 <- createM(s,N2)
#    print(M2)
    denum=vp_eq6(M01,M02,t2,s)
    denum=1.0/sum(denum[1:(N2-1)])
    out=vp_eq6(M1,M2,t2,s)*denum
    out[1]=1-sum(out[2:N2])+f0
    out[2:N2]=out[2:N2]/(1.0-f0)
    return(out)
}

createM <- function(s,N)
{
    return(outer(seq(1, N), seq(1,N), Mjk,n=N,s))
    
}


probBinGam <- function(s,N2,t2,N1,f0,alpha,beta,ns)
    {
        for(i in 0:(ns-1))
            {
                binomvec <- pbinom(i,ns,seq(0,N2)/(2.0*N2))
                sumv=vs(N1,N2,s,t2,f0)*binomvec
                tointegrate=sumv*dgamma(s,alpha,beta)
            }
        return(tointegrate)                
    }

loglike <- function(N2,t2,N1,f0,alpha,beta,ns,sfsobs)
{
    p=0
    ll=integrate(probBinGam,0,Inf,N2,t2,N1,f0,alpha,beta,ns)
    p=p + sfs[i]*log(ll)
}
## equation 8

#?outer
#fille the matrix
##https://stat.ethz.ch/pipermail/r-help/2011-November/296893.html
N1=100
M1 <- outer(seq(1, N1), seq(1,N1), Mjk,n=N1,s=1)
un <- x(100,N1,M1) #compute
u(M1)

N2=500
M2 <- outer(seq(1, N2), seq(1,N2), Mjk,n=N2,s=1)

ws(M1,M2)

myvs <- vs(100,200,3,100)
plot(1:200,myvs)
myvs <- vs(100,200,0.3,100)
plot(1:200,myvs)
myvs <- vs(100,200,0,100) # must equal 1 according to paper
sum(myvs)## -> check yes
plot(1:200,myvs)
    
testy <- function(x,y,z){return(x+y+z)}
A3 <- outer(seq(1, N), seq(1,N), testy,z=1)


## now workon M to see if it's diagonalizable (so make the mulplication faster)

## check if matrix is diagonalizable
p <- eigen(M)$vectors # not real
d <- diag(eigen(M)$values) # not real
newM <- p%*%d%*%solve(p) # P D P^-1 should equal M
#M-newM #non zero, so not diagonalizable
M[1:10,1:10]
newM[1:10,1:10]
##solve(p) -> x such that px=I -> p inversed
## not working :(

## chekc the svd of M
Msvd=svd(M)
D=diag(Msvd$d)
newM <- Msvd$u %*% D  %*% t(Msvd$v)
#newM-M actually works tiny bits of difference due to approx
