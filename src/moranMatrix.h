#ifndef MORANMATRIX_H
#define MORANMATRIX_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "Sparse"
#include "Dense"
#include "Eigenvalues"
#include <boost/math/distributions/binomial.hpp> // for binomial_distribution
#include <boost/math/distributions/gamma.hpp> // gamma distribution
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>


using boost::math::binomial_distribution;
using namespace Eigen;
typedef MatrixXd tmatrix;
typedef SparseMatrix<double> SpMat;
typedef VectorXd tvector;
typedef Eigen::Triplet<double> T;

// this class defines the transition matrix of a moran model
// In a population of constant size, at each time step, one individual dies and one new individual is binomial_distribution
// this transition matrix will be dependant on the mutation rate mmu, the population size N and the selection factor s


class Mmatrix
{
public:
  Mmatrix();
  Mmatrix(double s, int N,double mmu);
  ~Mmatrix();
  Mmatrix(const Mmatrix &);
  double computeEG(); // compute P D and P1 (eigenmatrix,eginvector and inverse of eigen matrix)
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
 double Mjk(int,int,int,double);
  tmatrix matrixpowerD(int n);
  tmatrix matrixpowerM(int n);
  SpMat M;
  tmatrix P;
  tmatrix D;
  tmatrix P1;
  int popsize; // store population size = N = matrix dimension
  double s; // store selection factor
  double u; // mutation rate
  double v; // mutatiore rate
};


#endif
