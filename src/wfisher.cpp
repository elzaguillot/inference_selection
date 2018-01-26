#include "wfisher.h"
#include "utils.h"

WFmatrix::WFmatrix()
{
  popsize=0;
  popsize1=0;
  popsize2=0;
  s=0;
  u=0;
  v=0;
}

WFmatrix::WFmatrix(double s1, int N, double mmu)
{
  u=mmu;
  v=mmu;
  double t=time(0);
  s=s1;
  tmatrix mat(N+1,N+1);
  for(int i(0);i<N+1;++i)
    for(int j(0);j<N+1;++j)
      {
	mat(i,j)=Mjk(i,j,N,N,s1);
	//	std::cout<< i << " " << j << "  " << mat(i,j) << std::endl;
      }
  M=mat;
  popsize1=N+1;
  popsize2=N+1;
  popsize=popsize1;
  //  computeEG();
  //  std::cout<< "Matrix of size " << N+1<< " and s=" << s << " finished in " << time(0)-t << std::endl;
  
}

WFmatrix::WFmatrix(double s1, int N1, int N2,double mmu)
{
  u=mmu;
  v=mmu;
  double t=time(0);
  s=s1;
  tmatrix mat(N1+1,N2+1);
  for(int i(0);i<N1+1;++i)
    for(int j(0);j<N2+1;++j)
      mat(i,j)=Mjk(i,j,N1,N2,s1);
  M=mat;
  popsize1=N1+1;
  popsize2=N2+1;
  popsize=popsize2;
  //  computeEG();
  //  std::cout<< "Matrix of size " << N+1<< " and s=" << s << " finished in " << time(0)-t << std::endl;  
}

WFmatrix::~WFmatrix()
{}

WFmatrix::WFmatrix(const WFmatrix& a)
{
  M=a.M;
  P=a.P;
  P1=a.P1;
  popsize1=a.popsize1;
  popsize2=a.popsize2;
  s=a.s;
  u=a.u;
  v=a.v;
}

// double WFmatrix::Mjk(int j, int k, int popsize, int s)
// { // function to fill the moran transition matrix
//   double a=j*(popsize-j)*(1.0-s)/(popsize*popsize));
//   double b=j*(popsize-j)*(1.0)/(popsize*popsize);
//   //  std::cout << "s" << s << " a " << a << " b " << b << std::endl;
//   if(j==(k+1))
//     return(a);
//   if(j==(k-1))
//     return(b);
//   if(j==k)
//     return(1-b-a);
//   return(0);
// }

double WFmatrix::Mjk(int j, int k, int popsize1, int popsize2, double s)
{ // function to fill the moran transition matrix
  double myp=(1-s)*j/(1.0*(1-s)*j+popsize1-j);
  double x=nChoosek(popsize2,k);
  double a = pow(myp,k)*pow((1-myp),popsize2-k)*x;
  //  std::cout << j << ":"<<k  << " - " << myp << " myp " << a << " a " << x << " b " << std::endl;
  return(a);
}

double WFmatrix::computeEG()
{
  /*  if(s==0)
    {
      SelfAdjointEigenSolver<MatrixXd> es(M); // diagonalize matrix
      D=getM(es).real();
      P=getP(es).real();      
    }
    else*/
    {
      EigenSolver<MatrixXd> es(M); // diagonalize matrix
      D=getD(es).real();
      P=getP(es).real();
    }
  P1=P.inverse();
}

tmatrix WFmatrix::matrixpowerD(int n)
{
  int N = D.rows();
  double j(0);
  for(int i(0);i<popsize1;++i)
    {      
      D(i,i)=pow(D(i,i),n);
    }
  //  tmatrix A= P*D*P1;
  return(D);
}

tmatrix WFmatrix::matrixpowerM(int n)
{
  return(P*matrixpowerD(n)*P1);
}

