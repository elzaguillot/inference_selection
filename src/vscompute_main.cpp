// script to compile vs value for different value of s
//g++ -O3 -std=c++11 precompile.cpp moran.cpp moranMatrix.cpp wfisher.cpp utils.cpp -o precompile

#include "moran.h"


int main(int argc, char** argv)
{
  if(argc>3)
    {
      double f0=0.2;
      double mu=0.00001;
      int popsize1 = atoi(argv[1]);
      int popsize2 = atoi(argv[2]);
      int nbgen = atoi(argv[3]);
      //      for(int popsize1(50);popsize1<10000; popsize1 += 50)
      //	{      
      //	  for(int popsize2(popsize1);popsize2<2000; popsize2 += 100)
      //	    {
      std::string out="s_values_"+std::to_string(popsize2)+".txt";
      std::ofstream fichier(out.c_str(),std::ios::out | std::ios::app);  // on ouvre le fichier en lecture
      Mmatrix M0N1(0,popsize1,mu);
      Mmatrix M0N2(0,popsize2,mu);
      //      for(double nbgen(1000);nbgen<100000;nbgen += 1000)
      //	{
      for(double s(0);s<0.2;s += 0.01)
	{      
	  Mmatrix MsN1(s,popsize1,mu);
	  Mmatrix MsN2(s,popsize2,mu);
	  tvector sumv = vs(popsize1,popsize2,s,nbgen,f0,M0N1,M0N2,MsN1,MsN2);
	  fichier << popsize1 << " " << popsize2 << " " << s << " " << nbgen << " ";
	  for(int i(0);i<popsize2+1;++i)
	    {
	      fichier << sumv(i) << " ";
	    }
	  fichier << std::endl;
	}
      fichier.close();
      //}
	      //	      std::cout << popsize1 << " " << popsize2 << std::endl;
      //    }
      //	}
    }
  return(0);
}
