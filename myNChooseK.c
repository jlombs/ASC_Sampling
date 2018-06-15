/*==========================================================
 * myNChooseK.c - MEX binomial coefficient calculator
 *
 * The calling syntax is:
 *
 *		out = myNChooseK(n, k)
 *
 *n and k must be double values. no input checking or output accuracy is done. out is a double.
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"

/* The computational routine */
void myNChooseK(double n[0], double k[0], double out[0])
{  
  double i, ai, bi;
  
  bi = n[0] - k[0];
  ai = bi + 1;

  for (i = 2; i <= k[0]; i++) {
     ai += (ai * bi) / i;
  }
  
  out[0] = ai;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *n,*k,*out;  
    
    /* get the value of the scalar input  */
    n = mxGetPr(prhs[0]);
    k = mxGetPr(prhs[1]);
    plhs[0]= mxCreateDoubleMatrix((mwSize)1,(mwSize)1,mxREAL);
    out = mxGetPr(plhs[0]);
  
    /* call the computational routine */
    myNChooseK(n,k,out);
}