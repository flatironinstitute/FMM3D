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
  int nt=1999;
  int nd=5;
  int ier=0;
  double *source = (double *)malloc(3*ns*sizeof(double));
  double *targ = (double *)malloc(3*nt*sizeof(double));

  double *charge = (double *)malloc(ns*nd*sizeof(double));
  double *dipvec = (double *)malloc(3*ns*nd*sizeof(double));

  double *pot = (double *)malloc(ns*nd*sizeof(double));
  double *pottarg = (double *)malloc(3*ns*nd*sizeof(double));

// initialize arrays
  for(int i=0;i<ns;i++)
  {
    source[3*i] = pow(rand01(),2);
    source[3*i+1] = pow(rand01(),2);
    source[3*i+2] = pow(rand01(),2);
  }
  for(int i=0;i<nd*ns;i++)
  {
    charge[i] = rand01(); 
    dipvec[3*i] = rand01();
    dipvec[3*i+1] = rand01();
    dipvec[3*i+2] = rand01();
  }

  for(int i=0;i<nt;i++)
  {
    targ[3*i] = rand01();
    targ[3*i+1] = rand01();
    targ[3*i+2] = rand01();
  }

  double eps = 0.5e-6;
  cprin_message("this code is an example c driver.");
  cprin_message("on output, the code prints sample pot,pottarg");
  cprin_skipline(2);

// call the fmm routine
  lfmm3d_st_cd_p_vec_(&nd, &eps, &ns, source, charge, dipvec, 
    pot, &nt, targ, pottarg, &ier);


  cprind("pot=",pot,12);
  cprind("grad=",pottarg,12);

  return 0;
}
