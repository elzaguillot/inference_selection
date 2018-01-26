// Moran.cpp
// adapt Eyre-Walker DFE to use moran matrix instead
// hopefully go faster
// by eguillot 

//#include</home/eguillot/softwarre/CBLAS/include/cblas.h>

// *** NOTES ***

// banded matrix of ublas perfect for moran

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//banded_matrix<T [, F, A]>   m(size1, size2, n_lower, n_upper); ublas
//typedef boost::numeric::ublas::matrix<double> tmatrix;
//typedef boost::numeric::ublas::vector<double> tvector;

#include "moran.h"
#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
 
 
RcppExport SEXP moranS(SEXP x, SEXP y){
  boost::random::mt19937 rng;
  //  tmatrix moran(10,10);
  std::string sfsfilename;
  //  std::cout << moran << std::endl;
  //  std::cout << matrixpower(moran,10) << std::endl;
  tparam myparam;
  Rcpp::NumericVector vector1(x);
  Rcpp::String str1(y);  
  sfsfilename=str1;
  tvector sfs=readSFS(sfsfilename);
  myparam.popsize2=vector1[0];
  myparam.nbgen=vector1[1];
  myparam.popsize1=vector1[2];
  myparam.f0=vector1[3];
  myparam.nsample=floor(sfs.rows());
  myparam.nbsite=sfs.sum();
  myparam.mu=vector1[4];
  myparam.alpha=vector1[5];
  myparam.beta=vector1[6];
  double ll =loglike2(myparam,sfs);
  if((ll<-99999)||(std::isnan(ll)))
    ll=-99999;
  std::cout<< vector1[0] << " " << vector1[1] << " " <<vector1[2] << " " <<vector1[3] << " " <<vector1[4] << " " <<vector1[5] << " " <<vector1[6] << " " <<" \n LL:" << ll << std::endl;
  return(Rcpp::wrap(ll));
}


RcppExport SEXP moran0(SEXP x, SEXP y){
  boost::random::mt19937 rng;
  //  tmatrix moran(10,10);
  std::string sfsfilename;
  //  std::cout << moran << std::endl;
  //  std::cout << matrixpower(moran,10) << std::endl;
  tparam myparam;
  Rcpp::NumericVector vector1(x);
  Rcpp::String str1(y);  
  sfsfilename=str1;
  tvector sfs=readSFS(sfsfilename);
  myparam.popsize2=vector1[0];
  myparam.nbgen=vector1[1];
  myparam.popsize1=vector1[2];
  myparam.f0=vector1[3];
  myparam.nsample=floor(sfs.rows());
  myparam.nbsite=sfs.sum();
  myparam.mu=vector1[4];
  myparam.alpha=vector1[5];
  myparam.beta=vector1[6];
  double ll =loglike(myparam,sfs);
  if((ll<-99999)||(std::isnan(ll)))
    ll=-99999;
  std::cout<< vector1[0] << " " << vector1[1] << " " <<vector1[2] << " " <<vector1[3] << " " <<vector1[4] << " " <<vector1[5] << " " <<vector1[6] << " " <<" \n LL:" << ll << std::endl;
  return(Rcpp::wrap(ll));
}
