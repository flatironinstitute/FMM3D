#include "utils.h"
#include "math.h"

void dzero(int n, double *a)
{
  for(int i=0; i<n;i=i+1)
  {
    a[i] = 0;
  }
  return;
}


void czero(int n, CPX *a)
{
  for(int i=0; i<n;i=i+1)
  {
    a[i] = 0;
  }
  return;
}


void comp_err_helm(int n, int pg, int pgt, CPX *pot, CPX *potex, CPX *pottarg,
      CPX *pottargex, CPX *grad, CPX *gradex, CPX *gradtarg, CPX *gradtargex, double *err)
{
  double e = 0;
  double r = 0;
  if( pg >= 1)
  {
    for(int i=0; i<n; i++)
    {
      r = r + pow(cabs(potex[i]),2);
      e = e +  pow(cabs(potex[i]-pot[i]),2);
    }
  }

  if( pg == 2)
  {
    for(int i=0; i<3*n; i++)
    {
      r = r + pow(cabs(gradex[i]),2);
      e = e +  pow(cabs(gradex[i]-grad[i]),2);
    }
  }
  *err = sqrt(e/r);
  return;
}
