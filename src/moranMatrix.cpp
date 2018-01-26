#include "moranMatrix.h"

Mmatrix::Mmatrix()
{
  popsize=2;
  s=0;
  u=0;
  v=0;
}

Mmatrix::Mmatrix(double s1, int N,double mmu)
{
  u=mmu;
  v=mmu;
  double t=time(0);
  s=s1;
  std::vector<T> tripletList;
  tripletList.reserve((N+1)*3);
  SpMat mat(N+1,N+1);
  for(int i(0);i<N+1;++i)
    for(int j(std::max(i-1,0));j<std::min(i+2,N+1);++j)
      tripletList.push_back(T(i,j,Mjk(i,j,N,s)));
  //  tripletList.push_back(T(0,0,0.999));
  //  tripletList.push_back(T(N,N,0.999));
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  M=mat;
  popsize=N+1;
  //  computeEG();
  //  std::cout<< "Matrix of size " << N+1<< " and s=" << s << " finished in " << time(0)-t << std::endl;  
}

Mmatrix::~Mmatrix()
{}

Mmatrix::Mmatrix(const Mmatrix& a)
{
  M=a.M;
  P=a.P;
  P1=a.P1;
  popsize=a.popsize;
  s=a.s;
  u=a.u;
  v=a.v;
}

// double Mmatrix::Mjk(int j, int k, int popsize, int s)
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


double Mmatrix::Mjk(int j, int k, int popsize, double s)
{ // function to fill the moran transition matrix
  double myp=(double) j*1.0/popsize;
  double b=(1-myp)*((1-u)*myp+v*(1-myp));
  double a=((1-myp)*(1-v)*myp+u*myp*myp)*(1-s);
  //  std::cout << j << ":"<<k  << " - " << myp << " myp " << a << " a " << b << " b " << std::endl;
  if(j==(k+1))
    return(b);
  if(j==(k-1))
    return(a);
  if(j==k)
    if(j==0)
	return(1-b);
    else if(j==popsize)
      return(1-a);
    else
      return(1-b-a);
  return(0);
}

double Mmatrix::computeEG()
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

tmatrix Mmatrix::matrixpowerD(int n)
{
  int N = D.rows();
  double j(0);
  for(int i(0);i<N;++i)
    {      
      D(i,i)=pow(D(i,i),n);
    }
  //  tmatrix A= P*D*P1;
  return(D);
}

tmatrix Mmatrix::matrixpowerM(int n)
{
  return(P*matrixpowerD(n)*P1);
}
