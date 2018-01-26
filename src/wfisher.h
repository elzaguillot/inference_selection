#ifndef FISHERMATRIX_H
#define FISHERMATRIX_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "../Eigen/Sparse"
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include <boost/math/distributions/binomial.hpp> // for binomial_distribution
#include <boost/math/distributions/gamma.hpp> // gamma distribution
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>
#include "utils.h"

using boost::math::binomial_distribution;
using namespace Eigen;
typedef MatrixXd tmatrix;
typedef VectorXd tvector;
typedef Eigen::Triplet<double> T;


class WFmatrix
{
public:
  WFmatrix();
  WFmatrix(double s, int N,double mmu);
  WFmatrix(double s, int N1,int N2,double mmu);
  ~WFmatrix();
  WFmatrix(const WFmatrix &);
  double computeEG();
  //  MatrixXcd getP(EigenSolver<MatrixXd> & es);
  //  MatrixXcd getM(EigenSolver<MatrixXd> & es);
  //  MatrixXcd getP(SelfAdjointEigenSolver<MatrixXd> & es);
  //  MatrixXcd getM(SelfAdjointEigenSolver<MatrixXd> & es);
  
  //template <typename solver>
  MatrixXcd getD(EigenSolver<MatrixXd> & es) { // must return the Diagonal matrix of eigenvalues D in M=P.D.P^(-1)
    //  EigenvalueType egval=es.eigenvalues();
  double t=time(0);
  int nbval= es.eigenvalues().rows();
  MatrixXcd out(nbval,nbval);
  for(int i=0;i<nbval;++i)
    out(i,i)=es.eigenvalues()[i];
  return(out);
}
//template <typename solver>
MatrixXcd getP(EigenSolver<MatrixXd> & es ){
  // must return the matrix of eigenvectors P in M=P.D.P^(-1)
  double t=time(0);
  int nbval= es.eigenvectors().rows();
  MatrixXcd out;
  out=es.eigenvectors();
  assert(out.cols()==nbval);
  assert(out.rows()==nbval);
  //  std::cout << out.imag().sum() << std::endl;
  assert(out.imag().sum()<2*10^-16);
  return(out); 
}
 double Mjk(int,int,int,int,double);
  tmatrix matrixpowerD(int n);
  tmatrix matrixpowerM(int n);
  tmatrix M;
  tmatrix P;
  tmatrix D;
  tmatrix P1;
  int popsize;
  int popsize1; // store population size = N = matrix dimension
  int popsize2; // store population size = N = matrix dimension
  double s; // store selection factor
  double u; // mutation rate
  double v; // mutatiore rate
};

// idea : create sparse matrices after all the computation so that one ca store hem easily



#endif
