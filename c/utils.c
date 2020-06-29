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


void comp_err_lap(int n, int pg, int pgt, double *pot, double *potex, double *pottarg,
      double *pottargex, double *grad, double *gradex, double *gradtarg, double *gradtargex, 
      double *hess, double *hessex, double *hesstarg, double *hesstargex, 
      double *err)
{
  double e = 0;
  double r = 0;
  if( pg >= 1)
  {
    for(int i=0; i<n; i++)
    {
      r = r + pow(fabs(potex[i]),2);
      e = e +  pow(fabs(potex[i]-pot[i]),2);
    }
  }

  if( pg >= 2)
  {
    for(int i=0; i<3*n; i++)
    {
      r = r + pow(fabs(gradex[i]),2);
      e = e +  pow(fabs(gradex[i]-grad[i]),2);
    }
  }

  if ( pg >= 3)
  {
    for(int i=0; i<6*n; i++)
    {
      r = r + pow(fabs(hessex[i]),2);
      e = e +  pow(fabs(hessex[i]-hess[i]),2);
    }
  }

  if( pgt >= 1)
  {
    for(int i=0; i<n; i++)
    {
      r = r + pow(fabs(pottargex[i]),2);
      e = e +  pow(fabs(pottargex[i]-pottarg[i]),2);
    }
  }

  if( pgt >= 2)
  {
    for(int i=0; i<3*n; i++)
    {
      r = r + pow(fabs(gradtargex[i]),2);
      e = e +  pow(fabs(gradtargex[i]-gradtarg[i]),2);
    }
  }

  if ( pgt >= 3)
  {
    for(int i=0; i<6*n; i++)
    {
      r = r + pow(fabs(hesstargex[i]),2);
      e = e +  pow(fabs(hesstargex[i]-hesstarg[i]),2);
    }
  }

  *err = sqrt(e/r);
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

  if( pgt >= 1)
  {
    for(int i=0; i<n; i++)
    {
      r = r + pow(cabs(pottargex[i]),2);
      e = e +  pow(cabs(pottargex[i]-pottarg[i]),2);
    }
  }

  if( pgt == 2)
  {
    for(int i=0; i<3*n; i++)
    {
      r = r + pow(cabs(gradtargex[i]),2);
      e = e +  pow(cabs(gradtargex[i]-gradtarg[i]),2);
    }
  }
  *err = sqrt(e/r);
  return;
}
