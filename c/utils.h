#ifndef UTILS_H
#define UTILS_H
  
#include <stdint.h>
#include <complex.h>

typedef double complex CPX;
#define rand01() ((double)rand()/RAND_MAX)

void dzero(int n, double *a);
void czero(int n, CPX *a);

void comp_err_helm(int n,int pg, int pgt, CPX *pot, CPX *potex, CPX *pottarg, CPX *pottargex,
         CPX *grad, CPX *gradex, CPX *gradtarg, CPX *gradtargex, double *err);

void comp_err_lap(int n,int pg, int pgt, double *pot, double *potex, 
         double *pottarg, double *pottargex, double *grad, 
         double *gradex, double *gradtarg, double *gradtargex, 
         double *hess, double *hessex, double *hesstarg, 
         double *hesstargex, double *err);
#endif
