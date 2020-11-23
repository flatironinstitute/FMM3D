#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "hfmm3d_c.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  cprin_skipline(2);


  int ns=2000;
  int ier=0;
  double *source = (double *)malloc(3*ns*sizeof(double));

  CPX *charge = (CPX *)malloc(ns*sizeof(CPX));

  CPX *pot = (CPX *)malloc(ns*sizeof(CPX));
  CPX *grad = (CPX *)malloc(3*ns*sizeof(CPX));

// initialize arrays
  for(int i=0;i<ns;i++)
  {
    source[3*i] = pow(rand01(),2);
    source[3*i+1] = pow(rand01(),2);
    source[3*i+2] = pow(rand01(),2);

    charge[i] = rand01() + I*rand01();

  }

  double eps = 0.5e-6;
  CPX zk = 1.1 + 0.01*I;
  cprin_message("this code is an example c driver.");
  cprin_message("on output, the code prints sample pot,grad");
  cprin_skipline(2);

// call the fmm routine
  hfmm3d_s_c_g_(&eps, &zk, &ns, source, charge, pot, grad, &ier);



  cprinz("pot=",pot,12);
  cprinz("grad=",grad,12);

  return 0;
}
