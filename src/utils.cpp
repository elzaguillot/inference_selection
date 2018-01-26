#include "utils.h"



Eigen::MatrixXd matrixpower(Eigen::MatrixXd &mat, int n)
{
  int  tsize = mat.rows();
  //  std::cout << tsize << std::endl;
  Eigen::MatrixXd result(tsize,tsize);
  if(n==1)
    return(mat);
  while(n>0)
    {
      if (n % 2 != 0) {
	result *=  mat;
	n--;
      }
      mat=mat*mat;
      //mat = mat2;
      n /= 2;
    }
  return(result);
}

unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}
