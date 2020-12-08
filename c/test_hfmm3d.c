#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "hfmm3d_c.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  int ns=4000;
  int nt=3999;
  double *source = (double *)malloc(3*ns*sizeof(double));
  double *targ = (double *)malloc(3*nt*sizeof(double));

  CPX *charge = (CPX *)malloc(ns*sizeof(CPX));
  CPX *dipvec = (CPX *)malloc(3*ns*sizeof(CPX));

  int ntest = 10;

  int ntests = 36;
  int ipass[36];

  for(int i=0;i<ntests;i++)
  {
    ipass[i] = 0;
  }


  double thresh = 1e-15;

  double err=0;
  int ier = 0;

  CPX *pot = (CPX *)malloc(ns*sizeof(CPX));
  CPX *potex = (CPX *)malloc(ntest*sizeof(CPX));
  CPX *grad = (CPX *)malloc(3*ns*sizeof(CPX));
  CPX *gradex = (CPX *)malloc(3*ntest*sizeof(CPX));

  CPX *pottarg = (CPX *)malloc(nt*sizeof(CPX));
  CPX *pottargex = (CPX *)malloc(ntest*sizeof(CPX));
  CPX *gradtarg = (CPX *)malloc(3*nt*sizeof(CPX));
  CPX *gradtargex = (CPX *)malloc(3*ntest*sizeof(CPX));

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
  CPX zk = 1.1 + 0.01*I;


  int itest = 0;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 0;

  hfmm3d_s_c_p_(&eps, &zk, &ns, source, charge, pot, &ier);

  czero(ntest,potex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, source, &ntest, potex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_c_g_(&eps, &zk, &ns, source, charge, pot, grad, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_s_d_p_(&eps, &zk, &ns, source, dipvec, pot, &ier);

  czero(ntest,potex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_d_g_(&eps, &zk, &ns, source, dipvec, pot, grad, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_cd_p_(&eps, &zk, &ns, source, charge, dipvec, pot, &ier);

  czero(ntest,potex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_cd_g_(&eps, &zk, &ns, source, charge, dipvec, pot, grad, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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

  hfmm3d_t_c_p_(&eps, &zk, &ns, source, charge, &nt, targ, pottarg, &ier);

  czero(ntest,pottargex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_c_g_(&eps, &zk, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, &ier);

  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_t_d_p_(&eps, &zk, &ns, source, dipvec, &nt, targ, pottarg, &ier);

  czero(ntest,pottargex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_d_g_(&eps, &zk, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, &ier);

  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_cd_p_(&eps, &zk, &ns, source, charge, dipvec, &nt, targ,
     pottarg, &ier);

  czero(ntest,pottargex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_cd_g_(&eps, &zk, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, &ier);

  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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

  hfmm3d_st_c_p_(&eps, &zk, &ns, source, charge, pot, 
    &nt, targ, pottarg, &ier);

  czero(ntest,potex);
  czero(ntest,pottargex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, source, &ntest, 
       potex, &thresh);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_c_g_(&eps, &zk, &ns, source, charge, pot, grad, 
      &nt, targ, pottarg, gradtarg, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, source, &ntest, potex, 
     gradex, &thresh);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_st_d_p_(&eps, &zk, &ns, source, dipvec, pot, 
     &nt, targ, pottarg, &ier);

  czero(ntest,potex);
  czero(ntest,pottargex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_d_g_(&eps, &zk, &ns, source, dipvec, pot, grad, 
     &nt, targ, pottarg, gradtarg, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_cd_p_(&eps, &zk, &ns, source, charge, dipvec, pot, 
     &nt, targ,  pottarg, &ier);

  czero(ntest,potex);
  czero(ntest,pottargex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_cd_g_(&eps, &zk, &ns, source, charge, dipvec, pot, grad,
     &nt,targ, pottarg, gradtarg, &ier);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  czero(ntest,pottargex);
  czero(3*ntest,gradtargex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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

  nd = 2;
  charge = (CPX *)malloc(ns*nd*sizeof(CPX));
  dipvec = (CPX *)malloc(3*ns*nd*sizeof(CPX));

  pot = (CPX *)malloc(ns*nd*sizeof(CPX));
  potex = (CPX *)malloc(ntest*nd*sizeof(CPX));
  grad = (CPX *)malloc(3*ns*nd*sizeof(CPX));
  gradex = (CPX *)malloc(3*ntest*nd*sizeof(CPX));

  pottarg = (CPX *)malloc(nt*nd*sizeof(CPX));
  pottargex = (CPX *)malloc(ntest*nd*sizeof(CPX));
  gradtarg = (CPX *)malloc(3*nt*nd*sizeof(CPX));
  gradtargex = (CPX *)malloc(3*ntest*nd*sizeof(CPX));



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

  hfmm3d_s_c_p_vec_(&nd, &eps, &zk, &ns, source, charge, pot, &ier);

  czero(nd*ntest,potex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, source, &ntest, potex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_c_g_vec_(&nd, &eps, &zk, &ns, source, charge, pot, grad, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_s_d_p_vec_(&nd, &eps, &zk, &ns, source, dipvec, pot, &ier);

  czero(nd*ntest,potex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_d_g_vec_(&nd, &eps, &zk, &ns, source, dipvec, pot, grad, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_cd_p_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, pot, &ier);

  czero(nd*ntest,potex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_s_cd_g_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, pot, grad, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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

  hfmm3d_t_c_p_vec_(&nd, &eps, &zk, &ns, source, charge, &nt, targ, pottarg, &ier);

  czero(nd*ntest,pottargex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_c_g_vec_(&nd, &eps, &zk, &ns, source, charge, &nt, targ, 
      pottarg, gradtarg, &ier);

  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_t_d_p_vec_(&nd, &eps, &zk, &ns, source, dipvec, &nt, targ, pottarg, &ier);

  czero(nd*ntest,pottargex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_d_g_vec_(&nd, &eps, &zk, &ns, source, dipvec, &nt, targ, 
     pottarg, gradtarg, &ier);

  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_cd_p_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, &nt, targ,
     pottarg, &ier);

  czero(nd*ntest,pottargex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_t_cd_g_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, &nt, 
     targ, pottarg, gradtarg, &ier);

  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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

  hfmm3d_st_c_p_vec_(&nd, &eps, &zk, &ns, source, charge, pot, 
    &nt, targ, pottarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*ntest,pottargex);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, source, &ntest, 
       potex, &thresh);
  h3ddirectcp_(&nd, &zk, source, charge, &ns, targ, &ntest, 
       pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_c_g_vec_(&nd, &eps, &zk, &ns, source, charge, pot, grad, 
      &nt, targ, pottarg, gradtarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, source, &ntest, potex, 
     gradex, &thresh);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, targ, &ntest, pottargex, 
     gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


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
  hfmm3d_st_d_p_vec_(&nd, &eps, &zk, &ns, source, dipvec, pot, 
     &nt, targ, pottarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*ntest,pottargex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_d_g_vec_(&nd, &eps, &zk, &ns, source, dipvec, pot, grad, 
     &nt, targ, pottarg, gradtarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_cd_p_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, pot, 
     &nt, targ,  pottarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*ntest,pottargex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     source, &ntest, potex, &thresh);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, 
     targ, &ntest, pottargex, &thresh);

  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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
  hfmm3d_st_cd_g_vec_(&nd, &eps, &zk, &ns, source, charge, dipvec, pot, grad,
     &nt,targ, pottarg, gradtarg, &ier);

  czero(nd*ntest,potex);
  czero(nd*3*ntest,gradex);
  czero(nd*ntest,pottargex);
  czero(nd*3*ntest,gradtargex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, 
     potex, gradex, &thresh);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, targ, &ntest, 
     pottargex, gradtargex, &thresh);
  comp_err_helm(nd*ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

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


  int isum = 0;
  for(int i=0;i<ntests;i++)
  {
    isum = isum + ipass[i];
  }

  cprinf("Number of tests out of 36 passed =",&isum,1);

  return 0;
}  
