// define all the necessary step to compute the moran model of dfe

#ifndef MORAN_H
#define MORAN_H

#include "moranMatrix.h"
#include "utils.h"


// each function contains the number of the equation that it matches in the original paper from
tvector readSFS(std::string filename);
tvector f_eq3(int nbgen,  Mmatrix & MsN2);
tmatrix f_eq3_sped(int nbgen,  tmatrix &D);
tvector f_eq3_bis(int nbgen, Mmatrix &MsN2);
tvector x_eq4(int nbgen,Mmatrix &MsN2);
tvector x_eq4_sped(int nbgen,Mmatrix &MsN2);
tvector xs1_eq4(double &s,int popsize);
tvector u_eq5(Mmatrix &MsN1,int popsize);
tvector ws_eq(Mmatrix &MsN1,Mmatrix &MsN2,int nbgen);
tvector vp_eq6(Mmatrix &MsN1,Mmatrix &MsN2,int nbgen,double s);
tvector vs(int popsize1,int popsize2,double s,int nbgen,double f0,Mmatrix & M0N1,Mmatrix & M0N2,Mmatrix & MsN1,Mmatrix & MsN2);
tvector v0(int popsize1,int popsize2,int nbgen,double f0,Mmatrix & M0N1,Mmatrix & M0N2);
double fgamma(double s, double  alpha,double  beta);
double probBinGam(double s,int popsize2,int nbgen,int popsize1,double f0,double alpha,double beta,int nbsample,int nbsites, double mu,Mmatrix & M0N1,Mmatrix & M0N2,tvector & sfsobs);
double probBinGam0(int popsize2,int nbgen,int popsize1,double f0,double alpha,double beta,int nbsample,int nbsites, double mu,Mmatrix & M0N1,Mmatrix & M0N2,tvector & sfsobs);
double pb(double s,tparam allparam, Mmatrix& M0N1,Mmatrix& M0N2,tvector & sfsobs);
double pb0(tparam allparam, Mmatrix& M0N1,Mmatrix& M0N2,tvector & sfsobs);
double integrate(double (*fptr)(double,tparam,Mmatrix&,Mmatrix&,tvector&),tparam allparam,double a, double b,int len,Mmatrix& M0N1,Mmatrix& M0N2,tvector& sfsobs);
double loglike(tparam myparam,tvector sfsobs);
double loglike2(tparam myparam,tvector sfsobs);


template<typename Derived>
Eigen::Block<Derived>
topLeftCorner(MatrixBase<Derived>& m, int rows, int cols)
{
  return Eigen::Block<Derived>(m.derived(), 0, 0, rows, cols);
}

template<typename Derived>
const Eigen::Block<const Derived>
topLeftCorner(const MatrixBase<Derived>& m, int rows, int cols)
{
  return Eigen::Block<const Derived>(m.derived(), 0, 0, rows, cols);
}

#endif
