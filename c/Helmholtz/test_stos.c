#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "hfmm3dpartwrap.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  int ns=2000;
  int nt=2000;
  double *source = (double *)malloc(3*ns*sizeof(double));
  CPX *charge = (CPX *)malloc(ns*sizeof(CPX));
  CPX *dipstr = (CPX *)malloc(ns*sizeof(CPX));
  double *dipvec = (double *)malloc(3*ns*sizeof(double));

  CPX *pot = (CPX *)malloc(ns*sizeof(CPX));

  for(int i=0;i<ns;i++)
  {
    source[3*i] = pow(rand01(),2);
    source[3*i+1] = pow(rand01(),2);
    source[3*i+2] = pow(rand01(),2);

    charge[i] = rand01() + I*rand01();
    dipstr[i] = rand01() + I*rand01();

    dipvec[3*i] = rand01();
    dipvec[3*i+1] = rand01();
    dipvec[3*i+2] = rand01();

  }

  cprind("source=",source,24);
  cprinz("charge=",charge,24);
  cprinz("dipstr=",dipstr,24);
  cprind("dipvec=",dipvec,24);

  double eps = 0.5e-6;
  CPX zk = 1.1 + 0.01*I;

  hfmm3dpartstoscp_(&eps, &zk, &ns, source, charge, pot);

  cprinz("pot=",pot,24);

  return 0;
}  
