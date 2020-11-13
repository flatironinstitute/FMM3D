#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "lfmm3d_c.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  cprin_skipline(2);


  int ns=2000;
  int ier=0;
  double *source = (double *)malloc(3*ns*sizeof(double));

  double *charge = (double *)malloc(ns*sizeof(double));

  double *pot = (double *)malloc(ns*sizeof(double));
  double *grad = (double *)malloc(3*ns*sizeof(double));

// initialize arrays
  for(int i=0;i<ns;i++)
  {
    source[3*i] = pow(rand01(),2);
    source[3*i+1] = pow(rand01(),2);
    source[3*i+2] = pow(rand01(),2);

    charge[i] = rand01(); 

  }

  double eps = 0.5e-6;
  cprin_message("this code is an example c driver.");
  cprin_message("on output, the code prints sample pot,grad");
  cprin_skipline(2);

// call the fmm routine
  lfmm3d_s_c_g_(&eps, &ns, source, charge, pot, grad, &ier);


  cprind("pot=",pot,12);
  cprind("grad=",grad,12);

  return 0;
}
