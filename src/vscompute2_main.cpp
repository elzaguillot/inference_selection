// script to compile vs value for different value of s

//g++ -O3 -std=c++11 precompile2.cpp  moranMatrix.cpp wfisher.cpp utils.cpp wrightFisher.cpp -o precompile2

#include "wrightFisher.h"


int main(int argc, char** argv)
{
  if(argc>3)
    {
      int popsize1 = atoi(argv[1]);
      int popsize2 = atoi(argv[2]);
      int nbgen = atoi(argv[3]);
      double f0=0.2;
      double mu=0.00001;
      //      for(int popsize1(50);popsize1<10000; popsize1 += 50)
      //	{      
      //	  for(int popsize2(popsize1);popsize2<2000; popsize2 += 100)
      //	    {
      std::string out="s_valuesWF_"+std::to_string(popsize2)+".txt";
      std::ofstream fichier(out.c_str(),std::ios::out | std::ios::app);  // on ouvre le fichier en lecture	      
      WFmatrix M0N1(0,popsize1,mu);
      WFmatrix M0N2(0,popsize2,mu);	      
      //      for(double nbgen(100);nbgen<1000;nbgen += 10)
      //		{
      for(double s(0);s<0.2;s += 0.01)
	{      
	  WFmatrix MsN1(s,popsize1,mu);
	  WFmatrix MsN2(s,popsize2,mu);
	  tvector sumv = vs(popsize1,popsize2,s,nbgen,f0,M0N1,M0N2,MsN1,MsN2);
	  fichier << popsize1 << " " << popsize2 << " " << s << " " << nbgen << " ";
	  for(int i(0);i<popsize2+1;++i)
	    {
	      fichier << sumv(i) << " ";
	    }
	  fichier << std::endl;
	}
      fichier.close();
      //	}
      //	      std::cout << popsize1 << " " << popsize2 << std::endl;
      //	    }
    }
  return(0);
}
