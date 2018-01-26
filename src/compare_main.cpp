// program to compare Wright Fisher and Moran computations

#include "moran.h"
#include "wrightFisher.h"

int main(int argc, char* argv[])
{
  
  boost::random::mt19937 rng;
  //  tmatrix moran(10,10);
  std::string sfsfilename;
  //  std::cout << moran << std::endl;
  //  std::cout << matrixpower(moran,10) << std::endl;
  tparam myparam;
  if(argc>1)
    sfsfilename=argv[1];
  else
    {
      std::cerr<<"Error: input filename for sfs"<< std::endl;
    }
  tvector sfs=readSFS(sfsfilename);
  myparam.popsize2=200;
  myparam.nbgen=100;
  myparam.popsize1=100;
  myparam.f0=0.05;
  myparam.nsample=floor(sfs.rows());
  myparam.nbsite=sfs.sum();
  myparam.mu=0.1;
  myparam.alpha=0.5;
  myparam.beta=0.5;
  tvector sfsobs= sfs;
  double ll=0;
  Mmatrix M0N1(0,myparam.popsize1,myparam.mu);
  Mmatrix M0N2(0,myparam.popsize2,myparam.mu);
  ll=log(pb0(myparam,M0N1,M0N2,sfsobs)); //eq11
  std::cout <<   ll << std::endl;
  Mmatrix M0N1b(0,myparam.popsize1,myparam.mu);
  Mmatrix M0N2b(0,myparam.popsize2,myparam.mu);
  ll=log(pb0(myparam,M0N1b,M0N2b,sfsobs)); //eq11
  std::cout <<   ll << std::endl;

}