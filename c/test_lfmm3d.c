#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "lfmm3d_c.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  int ns=4000;
  int nt=3999;
  double *source = (double *)malloc(3*ns*sizeof(double));
  double *targ = (double *)malloc(3*nt*sizeof(double));

  double *charge = (double *)malloc(ns*sizeof(double));
  double *dipvec = (double *)malloc(3*ns*sizeof(double));

  int ntest = 10;

  int ntests = 54;
  int ipass[54];
  int ier = 0;

  for(int i=0;i<ntests;i++)
  {
    ipass[i] = 0;
  }


  double thresh = 1e-15;

  double err=0;

  double *pot = (double *)malloc(ns*sizeof(double));
  double *potex = (double *)malloc(ntest*sizeof(double));
  double *grad = (double *)malloc(3*ns*sizeof(double));
  double *gradex = (double *)malloc(3*ntest*sizeof(double));
  double *hess = (double *)malloc(6*ns*sizeof(double));
  double *hessex = (double *)malloc(6*ntest*sizeof(double));

  double *pottarg = (double *)malloc(nt*sizeof(double));
  double *pottargex = (double *)malloc(ntest*sizeof(double));
  double *gradtarg = (double *)malloc(3*nt*sizeof(double));
  double *gradtargex = (double *)malloc(3*ntest*sizeof(double));
  double *hesstarg = (double *)malloc(6*nt*sizeof(double));
  double *hesstargex = (double *)malloc(6*ntest*sizeof(double));

  int nd = 1;

  int pg = 0;
  int pgt = 0;



  for(int i=0;i<ns;i++)
  {
    source[3*i] = pow(rand01(),2);
    source[3*i+1] = pow(rand01(),2);
    source[3*i+2] = pow(rand01(),2);

    charge[i] = rand01() + I*rand01();

    dipvec[3*i] = rand01() +I*rand01();
    dipvec[3*i+1] = rand01() + I*rand01();
    dipvec[3*i+2] = rand01() + I*rand01();

  }

  for(int i=0; i<nt; i++)
  {
    targ[3*i] = rand01();
    targ[3*i+1] = rand01();
    targ[3*i+2] = rand01();

  }


  double eps = 0.51e-3;


  int itest = 0;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 0;

  lfmm3d_s_c_p_(&eps, &ns, source, charge, pot, &ier);

  dzero(ntest,potex);
  l3ddirectcp_(&nd, source, charge, &ns, source, &ntest, potex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 2;
  pgt = 0;
  lfmm3d_s_c_g_(&eps, &ns, source, charge, pot, grad, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  l3ddirectcg_(&nd, source, charge, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest += 1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  

  pg = 3;
  pgt = 0;
  lfmm3d_s_c_h_(&eps, &ns, source, charge, pot, grad, hess, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  l3ddirectch_(&nd, source, charge, &ns, source, &ntest, potex, gradex, 
    hessex,&thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 1;
  pgt = 0;
  lfmm3d_s_d_p_(&eps, &ns, source, dipvec, pot, &ier);

  dzero(ntest,potex);
  l3ddirectdp_(&nd, source, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  lfmm3d_s_d_g_(&eps, &ns, source, dipvec, pot, grad, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  l3ddirectdg_(&nd, source, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 0;
  lfmm3d_s_d_h_(&eps, &ns, source, dipvec, pot, grad, hess, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  l3ddirectdh_(&nd, source, dipvec, &ns, source, &ntest, potex, 
    gradex, hessex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  
  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 1;
  pgt = 0;
  lfmm3d_s_cd_p_(&eps, &ns, source, charge, dipvec, pot, &ier);

  dzero(ntest,potex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  lfmm3d_s_cd_g_(&eps, &ns, source, charge, dipvec, pot, grad, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 0;
  lfmm3d_s_cd_h_(&eps, &ns, source, charge, dipvec, pot, grad, hess, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, source, &ntest, 
    potex, gradex, hessex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 0;
  pgt = 1;

  lfmm3d_t_c_p_(&eps, &ns, source, charge, &nt, targ, pottarg, &ier);

  dzero(ntest,pottargex);
  l3ddirectcp_(&nd, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 0;
  pgt = 2;
  lfmm3d_t_c_g_(&eps, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectcg_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest += 1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  


  pg = 0;
  pgt = 3;
  lfmm3d_t_c_h_(&eps, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectch_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 0;
  pgt = 1;
  lfmm3d_t_d_p_(&eps, &ns, source, dipvec, &nt, targ, pottarg, &ier);

  dzero(ntest,pottargex);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 0;
  pgt = 2;
  lfmm3d_t_d_g_(&eps, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectdg_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 0;
  pgt = 3;
  lfmm3d_t_d_h_(&eps, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectdh_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");

  
  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 0;
  pgt = 1;
  lfmm3d_t_cd_p_(&eps, &ns, source, charge, dipvec, &nt, targ,
     pottarg, &ier);

  dzero(ntest,pottargex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 0;
  pgt = 2;
  lfmm3d_t_cd_g_(&eps, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 0;
  pgt = 3;
  lfmm3d_t_cd_h_(&eps, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 1;

  lfmm3d_st_c_p_(&eps, &ns, source, charge, pot, 
    &nt, targ, pottarg, &ier);

  dzero(ntest,potex);
  dzero(ntest,pottargex);
  l3ddirectcp_(&nd, source, charge, &ns, source, &ntest, 
       potex, &thresh);
  l3ddirectcp_(&nd, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 2;
  pgt = 2;
  lfmm3d_st_c_g_(&eps, &ns, source, charge, pot, grad, 
      &nt, targ, pottarg, gradtarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectcg_(&nd, source, charge, &ns, source, &ntest, potex, 
     gradex, &thresh);
  l3ddirectcg_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest += 1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  


  pg = 3;
  pgt = 3;
  lfmm3d_st_c_h_(&eps, &ns, source, charge, pot, grad, hess,
      &nt, targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectch_(&nd, source, charge, &ns, source, &ntest, potex, 
     gradex, hessex, &thresh);
  l3ddirectch_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 1;
  pgt = 1;
  lfmm3d_st_d_p_(&eps, &ns, source, dipvec, pot, 
     &nt, targ, pottarg, &ier);

  dzero(ntest,potex);
  dzero(ntest,pottargex);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 2;
  lfmm3d_st_d_g_(&eps, &ns, source, dipvec, pot, grad, 
     &nt, targ, pottarg, gradtarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectdg_(&nd, source, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  l3ddirectdg_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 3;
  lfmm3d_st_d_h_(&eps, &ns, source, dipvec, pot, grad, hess,
     &nt, targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectdh_(&nd, source, dipvec, &ns, source, &ntest, 
     potex, gradex, hessex, &thresh);
  l3ddirectdh_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  
  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 1;
  pgt = 1;
  lfmm3d_st_cd_p_(&eps, &ns, source, charge, dipvec, pot, 
     &nt, targ,  pottarg, &ier);

  dzero(ntest,potex);
  dzero(ntest,pottargex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 2;
  lfmm3d_st_cd_g_(&eps, &ns, source, charge, dipvec, pot, grad,
     &nt,targ, pottarg, gradtarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 3;
  lfmm3d_st_cd_h_(&eps, &ns, source, charge, dipvec, pot, grad, hess,
     &nt,targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(ntest,potex);
  dzero(3*ntest,gradex);
  dzero(6*ntest,hessex);
  dzero(ntest,pottargex);
  dzero(3*ntest,gradtargex);
  dzero(6*ntest,hesstargex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, hessex, &thresh);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  free(charge);
  free(dipvec);
  free(pot);
  free(potex);
  free(grad);
  free(gradex);
  free(pottarg);
  free(pottargex);
  free(gradtarg);
  free(gradtargex);
  free(hess);
  free(hessex);
  free(hesstarg);
  free(hesstargex);


  nd = 2;
  charge = (double *)malloc(ns*nd*sizeof(double));
  dipvec = (double *)malloc(3*ns*nd*sizeof(double));

  pot = (double *)malloc(ns*nd*sizeof(double));
  potex = (double *)malloc(ntest*nd*sizeof(double));
  grad = (double *)malloc(3*ns*nd*sizeof(double));
  gradex = (double *)malloc(3*ntest*nd*sizeof(double));
  hess = (double *)malloc(6*ns*nd*sizeof(double));
  hessex = (double *)malloc(6*ntest*nd*sizeof(double));

  pottarg = (double *)malloc(nt*nd*sizeof(double));
  pottargex = (double *)malloc(ntest*nd*sizeof(double));
  gradtarg = (double *)malloc(3*nt*nd*sizeof(double));
  gradtargex = (double *)malloc(3*ntest*nd*sizeof(double));
  hesstarg = (double *)malloc(6*nt*nd*sizeof(double));
  hesstargex = (double *)malloc(6*ntest*nd*sizeof(double));



  for(int i=0;i<nd*ns;i++)
  {

    charge[i] = rand01() + I*rand01();

    dipvec[3*i] = rand01() +I*rand01();
    dipvec[3*i+1] = rand01() + I*rand01();
    dipvec[3*i+2] = rand01() + I*rand01();

  }

  itest += 1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 0;

  lfmm3d_s_c_p_vec_(&nd, &eps, &ns, source, charge, pot, &ier);

  dzero(nd*ntest,potex);
  l3ddirectcp_(&nd, source, charge, &ns, source, &ntest, potex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 2;
  pgt = 0;
  lfmm3d_s_c_g_vec_(&nd, &eps, &ns, source, charge, pot, grad, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  l3ddirectcg_(&nd, source, charge, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest += 1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  

  pg = 3;
  pgt = 0;
  lfmm3d_s_c_h_vec_(&nd, &eps, &ns, source, charge, pot, grad, hess, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  l3ddirectch_(&nd, source, charge, &ns, source, &ntest, potex, gradex, 
    hessex,&thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 1;
  pgt = 0;
  lfmm3d_s_d_p_vec_(&nd, &eps, &ns, source, dipvec, pot, &ier);

  dzero(nd*ntest,potex);
  l3ddirectdp_(&nd, source, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  lfmm3d_s_d_g_vec_(&nd, &eps, &ns, source, dipvec, pot, grad, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  l3ddirectdg_(&nd, source, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 0;
  lfmm3d_s_d_h_vec_(&nd, &eps, &ns, source, dipvec, pot, grad, hess, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  l3ddirectdh_(&nd, source, dipvec, &ns, source, &ntest, potex, 
    gradex, hessex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  
  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 1;
  pgt = 0;
  lfmm3d_s_cd_p_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, &ier);

  dzero(nd*ntest,potex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  lfmm3d_s_cd_g_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, grad, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to source");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 0;
  lfmm3d_s_cd_h_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, grad, hess, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, source, &ntest, 
    potex, gradex, hessex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 0;
  pgt = 1;

  lfmm3d_t_c_p_vec_(&nd, &eps, &ns, source, charge, &nt, targ, pottarg, &ier);

  dzero(nd*ntest,pottargex);
  l3ddirectcp_(&nd, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 0;
  pgt = 2;
  lfmm3d_t_c_g_vec_(&nd, &eps, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectcg_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest += 1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  


  pg = 0;
  pgt = 3;
  lfmm3d_t_c_h_vec_(&nd, &eps, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectch_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 0;
  pgt = 1;
  lfmm3d_t_d_p_vec_(&nd, &eps, &ns, source, dipvec, &nt, targ, pottarg, &ier);

  dzero(nd*ntest,pottargex);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 0;
  pgt = 2;
  lfmm3d_t_d_g_vec_(&nd, &eps, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectdg_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 0;
  pgt = 3;
  lfmm3d_t_d_h_vec_(&nd, &eps, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectdh_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");

  
  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 0;
  pgt = 1;
  lfmm3d_t_cd_p_vec_(&nd, &eps, &ns, source, charge, dipvec, &nt, targ,
     pottarg, &ier);

  dzero(nd*ntest,pottargex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 0;
  pgt = 2;
  lfmm3d_t_cd_g_vec_(&nd, &eps, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 0;
  pgt = 3;
  lfmm3d_t_cd_h_vec_(&nd, &eps, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 1;

  lfmm3d_st_c_p_vec_(&nd, &eps, &ns, source, charge, pot, 
    &nt, targ, pottarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*ntest,pottargex);
  l3ddirectcp_(&nd, source, charge, &ns, source, &ntest, 
       potex, &thresh);
  l3ddirectcp_(&nd, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest += 1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: gradients");
  cprin_skipline(2);
  


  pg = 2;
  pgt = 2;
  lfmm3d_st_c_g_vec_(&nd, &eps, &ns, source, charge, pot, grad, 
      &nt, targ, pottarg, gradtarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectcg_(&nd, source, charge, &ns, source, &ntest, potex, 
     gradex, &thresh);
  l3ddirectcg_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest += 1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges");
  cprin_message("output: hessians");
  cprin_skipline(2);
  


  pg = 3;
  pgt = 3;
  lfmm3d_st_c_h_vec_(&nd, &eps, &ns, source, charge, pot, grad, hess,
      &nt, targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectch_(&nd, source, charge, &ns, source, &ntest, potex, 
     gradex, hessex, &thresh);
  l3ddirectch_(&nd, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);


  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 1;
  pgt = 1;
  lfmm3d_st_d_p_vec_(&nd, &eps, &ns, source, dipvec, pot, 
     &nt, targ, pottarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*ntest,pottargex);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  l3ddirectdp_(&nd, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 2;
  lfmm3d_st_d_g_vec_(&nd, &eps, &ns, source, dipvec, pot, grad, 
     &nt, targ, pottarg, gradtarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectdg_(&nd, source, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  l3ddirectdg_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 3;
  lfmm3d_st_d_h_vec_(&nd, &eps, &ns, source, dipvec, pot, grad, hess,
     &nt, targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectdh_(&nd, source, dipvec, &ns, source, &ntest, 
     potex, gradex, hessex, &thresh);
  l3ddirectdh_(&nd, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  
  itest +=1;
  cprin_message("testing source to target");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 1;
  pgt = 1;
  lfmm3d_st_cd_p_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, 
     &nt, targ,  pottarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*ntest,pottargex);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  l3ddirectcdp_(&nd, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 2;
  lfmm3d_st_cd_g_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, grad,
     &nt,targ, pottarg, gradtarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  l3ddirectcdg_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  itest +=1;
  cprin_message("testing source to source+target");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: hessians");
  cprin_skipline(2);


  pg = 3;
  pgt = 3;
  lfmm3d_st_cd_h_vec_(&nd, &eps, &ns, source, charge, dipvec, pot, grad, hess,
     &nt,targ, pottarg, gradtarg, hesstarg, &ier);

  dzero(nd*ntest,potex);
  dzero(nd*3*ntest,gradex);
  dzero(nd*6*ntest,hessex);
  dzero(nd*ntest,pottargex);
  dzero(nd*3*ntest,gradtargex);
  dzero(nd*6*ntest,hesstargex);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, hessex, &thresh);
  l3ddirectcdh_(&nd, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, hesstargex, &thresh);
  comp_err_lap(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,hess,hessex,hesstarg,hesstargex,&err);

  if(err<eps)
  {
    ipass[itest] = 1;
  }

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  free(charge);
  free(dipvec);
  free(pot);
  free(potex);
  free(grad);
  free(gradex);
  free(pottarg);
  free(pottargex);
  free(gradtarg);
  free(gradtargex);
  free(hess);
  free(hessex);
  free(hesstarg);
  free(hesstargex);



  int isum = 0;
  for(int i=0;i<ntests;i++)
  {
    isum = isum + ipass[i];
  }

  cprinf("Number of tests out of 54 passed =",&isum,1);

  return 0;
}  
