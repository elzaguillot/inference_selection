#ifndef UTILS_H
#define UTILS_H


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


struct tparam
{
  int popsize2;
  int nbgen;
  int popsize1;
  double f0;
  double alpha;
  double beta;
  int nsample;
  int nbsite;
  double mu;
};


Eigen::MatrixXd matrixpower(Eigen::MatrixXd &mat, int n);
unsigned nChoosek( unsigned n, unsigned k );

template<unsigned n, unsigned r>
struct Choose {
  enum {value = (n * Choose<n-1, r-1>::value) / r};
};
 
template<unsigned n>
struct Choose<n, 0> {
  enum {value = 1};
};

// from https://www.quora.com/What-are-some-cool-C++-tricks
#endif
