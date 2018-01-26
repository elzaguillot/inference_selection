// Moran.cpp
// adapt Eyre-Walker DFE to use moran matrix instead
// hopefully go faster
// by eguillot

//#include</home/eguillot/softwarre/CBLAS/include/cblas.h>

// *** NOTES ***


#include <iostream>
#include <fstream>
#include <iomanip>
#include "../Eigen/Dense"
#include <boost/math/distributions/binomial.hpp> // for binomial_distribution
#include <boost/math/distributions/gamma.hpp> // gamma distribution
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <mpi.h>


using boost::math::binomial_distribution;
using namespace Eigen;
typedef MatrixXd tmatrix;
typedef VectorXd tvector;




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
};


double pb(double,tparam,int);
double integrate(double (*)(double,tparam),tparam ,double , double ,int,int );
double loglike(tparam,tvector sfsobs);
tvector readSFS(std::string filename);

int main(int argc, char* argv[])
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int nb_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    // Get the rank of the process
    int id_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &id_proc);
    // Get the name of the processor
    if(id_proc==0)
      {
	boost::random::mt19937 rng;
	tmatrix moran(10,10);
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
	myparam.popsize2=10000;
	myparam.nbgen=100;
	myparam.popsize1=2000;
	myparam.f0=0.2;
	myparam.alpha=0.2;
	myparam.beta=0.5;
	myparam.nsample=floor(sfs.rows()/2);
	myparam.nbsite=sfs.sum();
	//  double ll=integrate(&pb,myparam,0.1,0.3,3,sfsfilename);
	std::cout << loglike(myparam,sfs) << std::endl;
	// Finalize the MPI environment.
      }
    MPI_Finalize();
    return(0);
}



// read a file with SFS stored with 1 number per line, send the sfs back
tvector readSFS(std::string filename)
{
  std::string line;
  std::ifstream myfile(filename.c_str(), std::ifstream::in);
  std::vector<double> sfs;
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  sfs.push_back(std::stod(line));
	}
      myfile.close();
    }
  else std::cout << "Unable to open file";
  int size=sfs.size();
  tvector sfst(size);
  for(int i(0);i<size;++i)
    sfst(i)=sfs[i];
  return sfst;
}


tmatrix matrixpowerD(tmatrix &D, int n)
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

double Mjk(int j, int k, int popsize, int s)
{ // function to fill the moran transition matrix
  double a=j*(popsize-j)*(1.0-s)/(pow(popsize,2));
  double b=j*(popsize-j)*(1.0)/pow(popsize,2);
  if(j==(k+1))
    return(a);
  if(j==(k-1))
    return(b);
  if(j==k)
    return(1-b-a);
  return(0);
}


tvector f_eq3(int nbgen, int popsize2, tmatrix &D,tmatrix  &P,tmatrix &P1)
{
  assert(popsize2==P.cols());
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(0)=1; // [1,0,0,0,0,...]
  tmatrix myf = P*matrixpowerD(D,nbgen)*P1; // f(t)=f0*M^t
  f0 =f0.transpose()*myf;
  return(f0);
}


tmatrix f_eq3_sped(int nbgen, int popsize2, tmatrix &D)
{
  return(matrixpowerD(D,nbgen)); // f(t)=f0*M^t
}


tvector f_eq3_bis(int nbgen, int popsize2,tmatrix &M)
{
  assert(popsize2==M.cols());
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(0)=1; // [1,0,0,0,0,...]
  tmatrix myf=matrixpower(M,nbgen); // f(t)=f0*M^t
  f0=f0.transpose()*myf;
  return(f0);
}

MatrixXcd getM(EigenSolver<MatrixXd> &es)
{ // must return the Diagonal matrix of eigenvalues D in M=P.D.P^(-1)
  //  EigenvalueType egval=es.eigenvalues();
  double t=time(0);
  int nbval= es.eigenvalues().rows();
  MatrixXcd out(nbval,nbval);
  for(int i=0;i<nbval;++i)
    out(i,i)=es.eigenvalues()[i];
  std::cout << "Time to compute getM " << time(0)-t << std::endl;
  return(out);
}


MatrixXcd getP(EigenSolver<MatrixXd> &es)
{
  // must return the matrix of eigenvectors P in M=P.D.P^(-1)
  double t=time(0);
  int nbval= es.eigenvectors().rows();
  MatrixXcd out;
  out=es.eigenvectors();
  assert(out.cols()==nbval);
  assert(out.rows()==nbval);
  assert(out.imag().sum()==0);
  std::cout << "Time to compute getP " << time(0)-t << std::endl;
  return(out);
}

tvector x_eq4(int nbgen,tmatrix &M)
{ // sum of the contributions of all generations
  tvector t2(nbgen);
  t2.fill(0);
  int N = M.rows();
  assert(N==M.cols());
  double t=time(0);
  EigenSolver<MatrixXd> es(M); // diagonalize matrix
  tmatrix D,P,P1;
  D=getM(es).real();
  P=getP(es).real();
  std::cout << "Time to compute eigen decomposition " << time(0)-t << std::endl;
  t=time(0);
  P1=P.inverse();
  std::cout <<"time to compute .inverse() " << time(0) -t << std::endl;
  //  std::cout << " the eigen values in D are " << D << std::endl;
  tmatrix M2(N,N); // matrix tor return
  for(int i(0);i<nbgen;++i)
    {
      //      M2.row(i)=f_eq3_bis(i,N,M); // add contribution of each time
      t=time(0);
      M2.row(i)=f_eq3(i,N,D,P,P1); // add contribution of each time
      // fill one per row
      std::cout << "time to compute f3 for i " << i << " is " << time(0)-t << std::endl;
    }
  tvector out=M2.colwise().sum();
  assert(out.cols()==1);
  assert(out.rows()==N);
  return(out); // sum for each time at the same site
}




tvector x_eq4_sped(int nbgen,tmatrix &M)
{ // sum of the contributions of all generations
  tvector t2(nbgen);
  t2.fill(0);
  int N = M.rows();
  assert(N==M.cols());
  double t=time(0);
  EigenSolver<MatrixXd> es(M); // diagonalize matrix
  tmatrix D,P,P1;
  D=getM(es).real();
  P=getP(es).real();
  std::cout << "Time to compute eigen decomposition " << time(0)-t << std::endl;
  t=time(0);
  P1=P.inverse();
  std::cout <<"time to compute .inverse() " << time(0) -t << std::endl;
  t=time(0);
  //  std::cout << " the eigen values in D are " << D << std::endl;
  tmatrix M2(N,N); // matrix tor return
  M2.fill(0);
  for(int i(0);i<nbgen;++i)
    {
      M2.noalias() += f_eq3_sped(i,N,D); // add contribution of each time
    }
  int popsize2=M.cols();
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(0)=1; // [1,0,0,0,0,...]
  tvector out=(f0.transpose()*P*M2*P1);
  assert(out.cols()==1);
  assert(out.rows()==N);
  std::cout <<" sume over all generations " << time(0) -t << std::endl;
  return(out); // sum for each time at the same site
}


tvector xs1_eq4(double &s,int popsize)
{ // in case of s>1 the computation becomes MUCH simplier
  tvector myx(popsize);
  myx.fill(0);
  myx(0)=2.0/s;
  return(myx);
}

tvector u_eq5(tmatrix &M,int popsize)
{
  popsize--;
  tmatrix Q = M.block(0,0,popsize,popsize);
  Q.noalias() = MatrixXd::Identity(popsize,popsize)-Q;
  //std::cout << Q.inverse() << std::endl;
  return(Q.inverse().col(0));
}


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

tvector ws_eq(tmatrix &M1,tmatrix &M2,int nbgen)
{
  double t =time(0);
  int popsize1 = M1.rows();
  int popsize2 = M2.rows();
  tmatrix TrN1N2(popsize2,popsize1-1);
  tmatrix littlemat(popsize1-1,popsize1-1);
  littlemat.fill(1);
  topLeftCorner(TrN1N2,popsize1-1,popsize1-1) = littlemat; // fill the transition matrix with 1
  tvector result5 = u_eq5(M1,popsize1);  // get the vector of equilibrium at nbgen ago
  result5.noalias() = (result5.transpose() * TrN1N2.transpose()); // apply the transition to bigger matrix
  //  result5 =result5;
  t=time(0);
  EigenSolver<MatrixXd> es(M2); // diagonalize matrix
  tmatrix D,P,P1;
  D=getM(es).real();
  P=getP(es).real();
  P1=P.inverse();
  std::cout << "Time to compute eigen decomposition " << time(0)-t << std::endl;
  result5.noalias() = result5.transpose()*P*matrixpowerD(D,nbgen)*P1;
  std::cout << " time to compute littlematrix " << time(0) -t<< std::endl;
  return(result5);
}


tvector vp_eq6(tmatrix &M1,tmatrix &M2,int nbgen,double s)
{
  int popsize1 = M1.rows();
  assert(popsize1==M1.cols());
  int popsize2 = M2.rows();
  assert(popsize2==M2.cols());
  double t=time(0);
  if(s<1)
    {
      tvector vp =popsize1*ws_eq(M1,M2,nbgen);
      assert(vp.rows()==popsize2);
      std::cout << "ws done in " << time(0) -t << std::endl;
      t=time(0);
      vp += popsize2*x_eq4_sped(nbgen,M2); // easier to debug in two lines
      std::cout << "eq4 done in " << time(0) -t << std::endl;
      return(vp);
    }
  else
    {
      tvector vp = popsize1*ws_eq(M1,M2,nbgen);
      assert(vp.rows()==popsize2);
      vp += popsize2*xs1_eq4(s,popsize2);
      return(vp);
    }
}


tmatrix createM(double s, int popsize)
{
  tmatrix M(popsize,popsize);
  for(int i(0);i<popsize;++i)
    for(int j(std::max(i-1,0));j<std::min(i+2,popsize);++j)
      M(i,j) = Mjk(i+1,j+1,popsize,s);
  M(0,0)=1;
  M(popsize-1,popsize-1)=1-M(popsize-1,popsize-2);
  return(M);
}

tvector vs(int popsize1,int popsize2,double s,int nbgen,double f0)
{
  double t=time(0);
  tmatrix M01 = createM(0,popsize1);
  tmatrix M02 = createM(0,popsize2);
  tmatrix M1 = createM(s,popsize1);
  tmatrix M2 = createM(s,popsize2);
  std::cout << "time to create matrices " << time(0)-t << std::endl;
  //std::cout << "created M01" << (M01) << std::endl;
  //std::cout << "created M02" << (M02) << std::endl;
  //  std::cout << "created M1" << (M1) << std::endl;
  //std::cout << "created M2" << (M2) << std::endl;
  tvector myvs=vp_eq6(M1,M2,nbgen,s);
  tvector myv0=vp_eq6(M01,M02,nbgen,0.0);
  //  denum(0)=0;
  double scaleFactor = myv0.sum();
  myvs /= scaleFactor;
  myvs(0)=1-myvs.sum();
  myvs*=1.0/(1-f0); // scale by f0
  myvs(0) += f0;
  //  double denum2 =1.0/denum.sum();
  //  tvector out = vp_eq6(M1,M2,t2,s).transpose()*denum;
  return (myvs);
}

double fgamma(double s, double  alpha,double  beta)
{
  // todo
  double mygamma;
  mygamma=pow(beta,alpha)*pow(s,(alpha-1))*exp(-beta*s)/tgamma(alpha);
  return(mygamma);
}


double probBinGam(double s,int popsize2,int nbgen,int popsize1,double f0,double alpha,double beta,int nbsample,int nbsites,int i)
{ // likelihood
  std::cout << "Measuring for S = " << s << std::endl;
  double t=time(0);
  tvector binomvec(2*popsize2);
  binomvec.fill(0);
  for(int j(0);j<(2*popsize2);++j)
    {
      boost::math::binomial_distribution<int> bn(2*nbsample, j/(2.0*popsize2));
      binomvec(j) =  boost::math::pdf(bn,i);
    }
  //  tmatrix binomvect = binomvec.transpose();
  tvector sumv=vs(2*popsize1,2*popsize2,s,nbgen,f0);
  double newsumv = binomvec.transpose()*sumv;
  newsumv = newsumv*fgamma(s,alpha,beta);
  //  double tot=newsumv.sum();
  //  tot=newsumv;
  std::cout << "Measuring for S = " << s << " took " << time(0) -t << std::endl;
  return(newsumv);
}

double pb(double s, tparam allparam,int i)
{ // compute the probBinGam with only two parameters, s that will vary and the rest that will stay fixed
  return(probBinGam(s,allparam.popsize2,allparam.nbgen,allparam.popsize1,allparam.f0,allparam.alpha,allparam.beta,allparam.nsample,allparam.nbsite,i));
}


double integrate(double (*fptr)(double,tparam,int),tparam allparam,double a, double b,int len,int i)
{ // see http://stackoverflow.com/questions/2982167/integration-math-in-c
  int nb_proc;
  int my_n;
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
  // Get the rank of the process
  int id_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &id_proc);
  double my_a;
  double my_b;
  int source;
  MPI::Status status;
  int tag;
  int target;
  double total;
  double x;
  if ( id_proc == 0 )
  {
//
//  We want N to be the total number of evaluations.
//  If necessary, we adjust N to be divisible by the number of processors.
//
    my_n = len / ( nb_proc - 1 );
    len = ( nb_proc - 1 ) * my_n;
  }
  double my_total;
  if ( id_proc == 0 )
  {
    for (int q = 1; q <= nb_proc - 1; q++ )
    {
      my_a = ( ( double ) ( nb_proc - q     ) * a
             + ( double ) (     q - 1 ) * b )
             / ( double ) ( nb_proc     - 1 );

      target = q;
      tag = 1;
      MPI::COMM_WORLD.Send ( &my_a, 1, MPI::DOUBLE, target, tag );

      my_b = ( ( double ) ( nb_proc - q - 1 ) * a
             + ( double ) (     q     ) * b )
             / ( double ) ( nb_proc     - 1 );

      target = q;
      tag = 2;
      MPI::COMM_WORLD.Send ( &my_b, 1, MPI::DOUBLE, target, tag );
    }
    total = 0.0;
    my_total = 0.0;
  }
//
//  Processes receive MY_A, MY_B, and compute their part of the integral.
//
  else
  {
    source = 0;
    tag = 1;
    MPI::COMM_WORLD.Recv ( &my_a, 1, MPI::DOUBLE, source, tag, status );
    source = 0;
    tag = 2;
    MPI::COMM_WORLD.Recv ( &my_b, 1, MPI::DOUBLE, source, tag, status );
    my_total = 0.0;
    for (int i = 1; i <= my_n; i++ )
    {
      x = (( double)(my_n - i)*my_a + (double)(i - 1)*my_b )/(double)(my_n - 1);
      my_total = my_total + (*fptr)( x,allparam,i );
    }
    my_total = (my_b - my_a)*my_total/( double ) ( my_n );
  }
//
//  Each process sends its value to the master process.
//
  MPI::COMM_WORLD.Reduce ( &my_total, &total, 1, MPI::DOUBLE, MPI::SUM, 0 );
  return(total);
}


double loglike(tparam myparam,tvector sfsobs)
{
  double ll=0;
  for(int i(0);i<myparam.nsample;++i)
    ll+=i*log(integrate(&pb,myparam,0,5,1000,i)); //eq11
  return(ll);
  // p=p + sfsobs[i]*log(ll[[1]]);
}
