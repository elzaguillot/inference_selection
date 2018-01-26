
#include "moran.h"
#include "wfisher.h"


// read a file with SFS stored with 1 number per line, send the sfs back
tvector readSFS(std::string filename)
{
  std::string line;
  std::ifstream myfile;
  myfile.open(filename.c_str(), std::ifstream::in);
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
  double tot=0;
  for(int i(0);i<size;++i)
    {
      sfst(i)=sfs[i];
      tot+=sfs[i];
    }
  return(sfst*1.0/tot);
}

tvector f_eq3(int nbgen,  Mmatrix & MsN2)
{
  int popsize2 =MsN2.popsize;
  //  assert(popsize2==P.cols());
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(1)=1; // [1,0,0,0,0,...]
  //  tmatrix myf = MsN2.P*MsN2.matrixpowerD(nbgen)*MsN2.P1; // f(t)=f0*M^t
  for(int i(0);i<nbgen;++i)
    {
      f0 =f0.transpose()*MsN2.M;
    }
  return(f0);
}

tmatrix f_eq3_sped(int nbgen,  tmatrix &D)
{
  return(matrixpower(D,nbgen)); // f(t)=f0*M^t
}

tvector f_eq3_bis(int nbgen, Mmatrix &MsN2)
{
  int popsize2 = MsN2.popsize;
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(1)=1; // [1,0,0,0,0,...]
  tmatrix myf=MsN2.matrixpowerD(nbgen); // f(t)=f0*M^t
  f0=f0.transpose()*myf;
  return(f0);
}

tvector x_eq4(int nbgen,Mmatrix &MsN2)
{ // sum of the contributions of all generations
  tvector t2(nbgen);
  t2.fill(0);
  int N = MsN2.popsize;
  assert(N==MsN2.M.cols());
  double t=time(0);
  tmatrix M2(nbgen,N); // matrix tor return
  int popsize2 =MsN2.popsize;
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(1)=1; // [1,0,0,0,0,...]
  tvector fi(popsize2); // initial Allele distribution in population
  fi=f0;
  for(int i(0);i<nbgen;++i)
    {
      //      M2.row(i)=f_eq3_bis(i,N,M); // add contribution of each time
      //  assert(popsize2==P.cols());
      //  tmatrix myf = MsN2.P*MsN2.matrixpowerD(nbgen)*MsN2.P1; // f(t)=f0*M^t
      fi =fi.transpose()*MsN2.M;
      M2.row(i)=fi;//f_eq3(i,MsN2); // add contribtdution of each time
      //      std::cout << M2.row(i)<< std::endl;
      // fill one per row
      //      std::cout << "time to compute f3 for i " << i << " is " << time(0)-t << std::endl;
    }
  tvector out=M2.colwise().sum(); // must devide by nbgen to obtain proportion
  out=out*1.0/nbgen;
  assert(out.cols()==1);
  assert(out.rows()==N);
  return(out); // sum for each time at the same site
}

tvector x_eq4_sped(int nbgen,Mmatrix &MsN2)
{ // sum of the contributions of all generations
  tvector t2(nbgen);
  t2.fill(0);
  int N = MsN2.popsize;
  assert(N==MsN2.M.cols());
  //  std::cout << " the eigen values in D are " << D << std::endl;
  double t=time(0);
  tmatrix M2(N,N); // matrix tor return
  M2.fill(0);
  for(int i(0);i<nbgen;++i)
    {
      M2.noalias() += f_eq3_sped(i,MsN2.D); // add contribution of each time
    }
  int popsize2=MsN2.popsize;
  tvector f0(popsize2); // initial Allele distribution in population
  f0.fill(0);
  f0(1)=1; // [1,0,0,0,0,...]
  tvector out=(f0.transpose()*MsN2.P*M2*MsN2.P1); // to check :-/
  assert(out.cols()==1);
  assert(out.rows()==N);
  //  std::cout <<" sume over all generations " << time(0) -t << std::endl;
  return(out); // sum for each time at the same site
}

tvector xs1_eq4(double &s,int popsize)
{ // in case of s>1 the computation becomes MUCH simplier
  tvector myx(popsize);
  myx.fill(0);
  myx(0)=2.0/s;
  return(myx);
}

tvector u_eq5(Mmatrix &MsN1,int popsize)
{
  /*  popsize--;
  tmatrix Q = MsN1.M.block(0,0,popsize+1,popsize+1);
  EigenSolver<MatrixXd> es(Q);
  /*std::cout << "M "  << MsN1.M << std::endl;
  std::cout << "Q "  << Q << std::endl;
  Q.noalias() = MatrixXd::Identity(popsize,popsize)-Q;
  std::cout << "I-Q "  << Q << " inverse " << Q.inverse() << std::endl;
  std::cout << Q.inverse() << std::endl; */
  //return(Q.inverse().col(0));;es.eigenvectors().col(0);es.eigenvectors().col(0)
  //    std::cout << "equi " << es.eigenvectors().real().col(popsize) << std::endl;
  //    std::cout << "equi " << es.eigenvalues().real() << std::endl;
  tvector x0(popsize);
  tvector u;
  x0.fill(1.0/popsize);
  //  tmatrix dmat;
  //  dmat = MatrixXd(MsN1.M);
  //  std::cout << x0.sum() << " " << dmat.rowwise().sum() << std::endl;
  for(int i(0);i<popsize*popsize*100;++i)
    {
      x0=x0.transpose()*MsN1.M;
      //      x0=x0*MsN1.M;
    }
  return(x0);
}

tvector ws_eq(Mmatrix &MsN1,Mmatrix &MsN2,int nbgen)
{
  double t =time(0);
  int popsize1 = MsN1.popsize;
  int popsize2 = MsN2.popsize;
  /*  tmatrix TrN1N2(popsize2,popsize1);
  TrN1N2.fill(0);
  tmatrix littlemat(popsize1,popsize1);
  littlemat.fill(0);
  for(int i(0);i<popsize1;++i)
    littlemat.diagonal()[i]=1;
    topLeftCorner(TrN1N2,popsize1,popsize1) = littlemat; // fill the transition matrix with 1*/
  WFmatrix Ttemp(MsN1.s,popsize1-1,popsize2-1,MsN1.u);
  tmatrix TrN1N2=Ttemp.M;
  tvector result5 = u_eq5(MsN1,popsize1);  // get the vector of equilibrium at nbgen ago
  //  std::cout << "result5 "  << result5 << std::endl;
  //  std::cout << "result5 "  << result5.sum() << std::endl;
  //  std::cout << "result5 "  << TrN1N2 << std::endl;
  result5.noalias() = (result5.transpose() * TrN1N2); // apply the transition to bigger matrix
  //  std::cout << "result5 AFTER 1* "  << result5 << " -- " << nbgen << std::endl;
  for(int i(0);i<nbgen;++i)
    {
      result5.noalias()=result5.transpose()*MsN2.M;
    }
  //  result5 =result5;
  //  result5.noalias() = (result5.transpose()*MsN2.P);
  //  std::cout << "result5 AFTER 2* "  << result5 << std::endl;
  //  std::cout << "result5 AFTER 2* "  << MsN2.M << std::endl;
  //  result5.noalias() = result5.transpose()*MsN2.matrixpowerD(nbgen)*MsN2.P1;
  return(result5);
}

tvector vp_eq6(Mmatrix &MsN1,Mmatrix &MsN2,int nbgen,double s)
{
  int popsize1 = MsN1.popsize;
  int popsize2 = MsN2.popsize;
  double t=time(0);
  if(s<1)
    {
      tvector vp =popsize1*ws_eq(MsN1,MsN2,nbgen);
      //      std::cout << "VP1\n" << vp << std::endl;
      //std::cout << "vp1=" << vp << " xeq4 " << x_eq4(nbgen,MsN2) << std::endl;
      assert(vp.rows()==popsize2);
      t=time(0);
      vp += popsize2*x_eq4(nbgen,MsN2); // easier to debug in two lines
      //      std::cout << "VP2\n" << x_eq4(nbgen,MsN2) << std::endl;
      //      std::cout << "VP2\n" << vp.sum() << " - " << popsize1+popsize2 << std::endl;
      return(vp);
    }
  else
    {
      tvector vp = popsize1*ws_eq(MsN1,MsN2,nbgen);
      assert(vp.rows()==popsize2);
      vp += popsize2*xs1_eq4(s,popsize2);
      return(vp);
    }
}

tvector vs(int popsize1,int popsize2,double s,int nbgen,double f0,Mmatrix & M0N1,Mmatrix & M0N2,Mmatrix & MsN1,Mmatrix & MsN2)
{
  tvector myvs=vp_eq6(MsN1,MsN2,nbgen,s);
  tvector myv0=vp_eq6(M0N1,M0N2,nbgen,0.0);
  double scaleFactor = myv0.sum();
  myvs /= scaleFactor;
  myvs(0) += f0;
  myvs*=1.0/myvs.sum(); // scale by f0
  //  std::cout << myvs.sum() << std::endl;
  return (myvs);
}

tvector v0(int popsize1,int popsize2,int nbgen,double f0,Mmatrix & M0N1,Mmatrix & M0N2)
{
  tvector myv0=vp_eq6(M0N1,M0N2,nbgen,0.0);
  return (myv0);
}


double fgamma(double s, double  alpha,double  beta)
{
  // todo
  double mygamma;
  mygamma=pow(beta,alpha)*pow(s,(alpha-1))*exp(-beta*s)/tgamma(alpha);
  //  std::cout << "Gamma " << mygamma << std::endl;
  return(mygamma);
}

double probBinGam(double s,int popsize2,int nbgen,int popsize1,double f0,double alpha,double beta,int nbsample,int nbsites, double mu,Mmatrix & M0N1,Mmatrix & M0N2,tvector & sfsobs)
{ // likelihood
  Mmatrix MsN1(s,popsize1,mu);
  //  MsN1.computeEG();
  Mmatrix MsN2(s,popsize2,mu);
  //  MsN2.computeEG();
  double t=time(0);
  // first independantly of i
  tvector sumv = vs(popsize1,popsize2,s,nbgen,f0,M0N1,M0N2,MsN1,MsN2);
  // then i speciific
  double tot=10;
  for(int i(0);i<nbsample;++i)
    {
      tvector binomvec(popsize2+1);
      binomvec.fill(0);
      for(int j(0);j<(popsize2);++j)
	{
	  boost::math::binomial_distribution<double> bn(2*nbsample, (double) j/(popsize2));
	  binomvec(j) =  boost::math::pdf(bn,i);
	}
      double newsumv = binomvec.dot(sumv);
      //      std::cout << tot << " " << newsumv << " sumv  " << sumv << " binomvec " << binomvec << " " << std::endl;
      newsumv = newsumv*fgamma(s,alpha,beta);
      //      std::cout << tot << " " << newsumv << " " << sfsobs[i] << std::endl;
      tot*=pow(newsumv,sfsobs[i]);
      //      std::cout << tot << " " << newsumv << " " << sfsobs[i] << std::endl;
    }
  return(tot);
}


double probBinGam0(int popsize2,int nbgen,int popsize1,double f0,double alpha,double beta,int nbsample,int nbsites, double mu,Mmatrix & M0N1,Mmatrix & M0N2,tvector & sfsobs)
{ // likelihood
  double t=time(0);
  // first independantly of i
  tvector sumv = v0(popsize1,popsize2,nbgen,f0,M0N1,M0N2);
  // then i speciific
  double tot=10;
  for(int i(0);i<nbsample;++i)
    {
      tvector binomvec(popsize2+1);
      binomvec.fill(0);
      for(int j(0);j<(popsize2);++j)
	{
	  boost::math::binomial_distribution<double> bn(2*nbsample, (double) j/(popsize2));
	  binomvec(j) =  boost::math::pdf(bn,i);
	}
      double newsumv = binomvec.dot(sumv);
      //      std::cout << tot << " " << newsumv << " sumv  " << sumv << " binomvec " << binomvec << " " << std::endl;
      newsumv = newsumv;
      std::cout << tot << " " << newsumv << " " << sfsobs[i] << std::endl;
      tot*=pow(newsumv,sfsobs[i]);

      //      std::cout << tot << " " << newsumv << " " << sfsobs[i] << std::endl;
    }
  return(tot);
}

double pb(double s,tparam allparam, Mmatrix& M0N1,Mmatrix& M0N2,tvector & sfsobs)
{ // compute the probBinGam with only two parameters, s that will vary and the rest that will stay fixed
  return(probBinGam(s,allparam.popsize2,allparam.nbgen,allparam.popsize1,allparam.f0,allparam.alpha,allparam.beta,allparam.nsample,allparam.nbsite,allparam.mu,M0N1,M0N2,sfsobs));
}

double pb0(tparam allparam, Mmatrix& M0N1,Mmatrix& M0N2,tvector & sfsobs)
{ // compute the probBinGam with only two parameters, s that will vary and the rest that will stay fixed
  return(probBinGam0(allparam.popsize2,allparam.nbgen,allparam.popsize1,allparam.f0,allparam.alpha,allparam.beta,allparam.nsample,allparam.nbsite,allparam.mu,M0N1,M0N2,sfsobs));
}


double integrate(double (*fptr)(double,tparam,Mmatrix&,Mmatrix&,tvector&),tparam allparam,double a, double b,int len,Mmatrix& M0N1,Mmatrix& M0N2,tvector& sfsobs)
{ // see http://stackoverflow.com/questions/2982167/integration-math-in-c
  // here we specifically integrate over a gamma distribition
  double tot=0;
  std::vector<double> x(len),y(len);
  double step=(b-a)*1.0/(len-1);
  for(int j=0;j<len;++j)
    {
      x[j]=a+step*j;
      //std::cout << "s - " << x[j] << std::endl;
    }
  for(int j=0;j<len;++j)
    {
      y[j]=(*fptr)(x[j],allparam,M0N1,M0N2,sfsobs);
      if(j!=0){
	tot+=(y[j]+y[j-1])/2.0*(x[j]-x[j-1]);
	//	std::cout << y[i] << " "  << y[i-1] << " x " << x[i] << " " << x[i-1] << std::endl;
      }
    }
  return(tot);
}

double loglike(tparam myparam, tvector sfsobs)
{
  double ll=0;
  Mmatrix M0N1(0,myparam.popsize1,myparam.mu);
  Mmatrix M0N2(0,myparam.popsize2,myparam.mu);
  ll=log(pb0(myparam,M0N1,M0N2,sfsobs)); //eq11
  //  std::cout<< "log (L) = "  << ll << std::endl;
  return(ll);
}


double loglike2(tparam myparam,tvector sfsobs)
// loglikelihood of the this set of parameters
{
  double ll=0;
  Mmatrix M0N1(0,myparam.popsize1,myparam.mu);
  Mmatrix M0N2(0,myparam.popsize2,myparam.mu);
  double t=time(0);
  ll=log(integrate(&pb,myparam,0.001,0.92,20,M0N1,M0N2,sfsobs)); //eq11
  //  std::cout<< "log (L) = "  << ll << std::endl;
  return(ll);
}
//
