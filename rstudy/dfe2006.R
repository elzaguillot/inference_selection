
## Script to reproduce dfe estimation by Eyre Walker 2006
## Written by Elsa Guillot elsa.guillot@unil.ch
## Unfinished as of 10.12.17


library(cubature) # 2d integration
## eq1
Pn_eq2 <- function(j,Ne,u,rj,Ln,s,psi)
    {
        out=4*Ne*u*rj*Ln
        integ <- intDHQ(Ne,psi,n,j)
        return(out*t1)
    }


intDHQ <- function(Ne,psi,n,j)
{
    integrate(function(y) { 
                  sapply(y, function(y) {
                             integrate(function(x) DHQ(x,y,Ne,psi,n,j), 0, 1)$value
                         })
              }, -Inf, Inf)$value
}

DHQ <- function(x,s,Ne,psi,n,j)
    {        
        out=D_eq5(psi,Ne,s)*H_eq3(Ne,s,x)*Q_eq4(n,j,x)
 #       print(D_eq5(psi,Ne,s))
        print(H_eq3(Ne,s,x))
#        print(out)
        return(out)
    }


## eq2
Pi_eq2 <- function(j,Ne,u,rj,Li)
    {
        return(4*Ne*u*rj*Li*(1.0/j+1.0/(n-j)))
    }

## eq3
H_eq3 <- function(Ne,s,x)
    {
        if((s*x*(1-x))!=0)
            {
                print(c(x))
                return(1-exp(4*Ne*s*(1-s)))*2.0/(x*(1-x)*(1-exp(4*Ne*s)))
            }       
        else
            return(1)
    }

## eq4
Q_eq4<- function(n,j,x)
    {
        if(j==(n*1.0/2))
            {
                return( choose(n,j)*x^j*(1-x)^(n-j) )
            }
        else
            {
                return( choose(n,j)* (x^j * (1-x)^(n-j) + x^(n-j) * (1-x)^j) )
            }
    }
##eq5
D_eq5 <- function(psi,Ne,s)
    {
        lambda = Ne*s*psi
#        print(lambda)
        return(psi^lambda*s^(lambda-1)*exp(-1.0*psi*s)/gamma(lambda))
    }


Nsample=10
Npop=1000
s=0.1
mu=0.000001
t=seq(0,10)
rj=1
nbsitesN=100
nbsitesS=50
psi=0.00001
for( i in t)
{
    print(i)
    Pn_eq2(i,Npop,mu,rj,nbsitesN,s,psi)
}
