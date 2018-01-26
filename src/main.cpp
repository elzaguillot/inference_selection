
#include "moran.h"


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
  myparam.alpha=0.9;
  myparam.beta=0.1;
  
  //  for(double a(0);a<1;a+=0.01)
  //    for(double b(0);b<1;b+=0.01)
  //     {
  //	myparam.alpha=a;
  //	myparam.beta=b;
  //  double ll=integrate(&pb,myparam,0.1,0.3,3,sfsfilename);
  //	std::cout << a << " " << b << " " << loglike2(myparam,sfs) << std::endl;
	//      }

  std::cout <<   loglike(myparam,sfs) << std::endl;
  if(0)
    {
      myparam.popsize2=500;
      myparam.nbgen=100;
      myparam.popsize1=100;
      myparam.f0=0.05;
      myparam.nsample=floor(sfs.rows());
      myparam.nbsite=sfs.sum();
      myparam.mu=0.1;
      myparam.alpha=0.5;
      myparam.beta=0.2;
  
  //  for(double a(0);a<1;a+=0.01)
  //    for(double b(0);b<1;b+=0.01)
  //     {
  //	myparam.alpha=a;
  //	myparam.beta=b;
  //  double ll=integrate(&pb,myparam,0.1,0.3,3,sfsfilename);
  //	std::cout << a << " " << b << " " << loglike2(myparam,sfs) << std::endl;
	//      }
      std::cout <<   loglike(myparam,sfs) << std::endl;
         
      myparam.popsize2=500;
      myparam.nbgen=100;
      myparam.popsize1=100;
      myparam.f0=0.05;
      myparam.nsample=floor(sfs.rows());
      myparam.nbsite=sfs.sum();
      myparam.mu=0.1;
      myparam.alpha=0.3;
      myparam.beta=0.2;
  //  for(double a(0);a<1;a+=0.01)
  //    for(double b(0);b<1;b+=0.01)
  //     {
  //	myparam.alpha=a;
  //	myparam.beta=b;
  //  double ll=integrate(&pb,myparam,0.1,0.3,3,sfsfilename);
  //	std::cout << a << " " << b << " " << loglike2(myparam,sfs) << std::endl;
	//      }
      std::cout <<   loglike(myparam,sfs) << std::endl;
    }
  return(0);
}
