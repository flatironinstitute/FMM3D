#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "hfmm3d_c.h"
#include "complex.h"
#include "cprini.h"

int main(int argc, char **argv)
{
  cprin_init("stdout", "fort.13");
  int ns=2000;
  int nt=1999;
  double *source = (double *)malloc(3*ns*sizeof(double));
  double *targ = (double *)malloc(3*nt*sizeof(double));

  CPX *charge = (CPX *)malloc(ns*sizeof(CPX));
  CPX *dipvec = (CPX *)malloc(3*ns*sizeof(CPX));

  int ntest = 10;

  int ntests = 18;
  int ipass[18];

  for(int i=0;i<ntests;i++)
  {
    ipass[i] = 0;
  }


  double thresh = 1e-15;

  double err=0;

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


  double eps = 0.5e-6;
  CPX zk = 1.1 + 0.01*I;


  int itest = 0;
  cprin_message("testing source to source");
  cprin_message("interaction: charges");
  cprin_message("output: potentials");
  cprin_skipline(2);


  pg = 1;
  pgt = 0;

  hfmm3dpartstoscp_(&eps, &zk, &ns, source, charge, pot);

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
  hfmm3dpartstoscg_(&eps, &zk, &ns, source, charge, pot, grad);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectcg_(&nd, &zk, source, charge, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");




  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);
  

  pg = 1;
  pgt = 0;
  hfmm3dpartstosdp_(&eps, &zk, &ns, source, dipvec, pot);

  czero(ntest,potex);
  h3ddirectdp_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  cprin_message("testing source to source");
  cprin_message("interaction: dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  hfmm3dpartstosdg_(&eps, &zk, &ns, source, dipvec, pot, grad);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectdg_(&nd, &zk, source, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");


  

  cprin_message("testing source to source");
  cprin_message("interaction: charges+dipoles");
  cprin_message("output: potentials");
  cprin_skipline(2);

  pg = 1;
  pgt = 0;
  hfmm3dpartstoscdp_(&eps, &zk, &ns, source, charge, dipvec, pot);

  czero(ntest,potex);
  h3ddirectcdp_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, &thresh);

  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);


  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");



  cprin_message("testing source to source");
  cprin_message("interaction: charges + dipoles");
  cprin_message("output: gradients");
  cprin_skipline(2);


  pg = 2;
  pgt = 0;
  hfmm3dpartstoscdg_(&eps, &zk, &ns, source, charge, dipvec, pot, grad);

  czero(ntest,potex);
  czero(3*ntest,gradex);
  h3ddirectcdg_(&nd, &zk, source, charge, dipvec, &ns, source, &ntest, potex, gradex, &thresh);
  comp_err_helm(ntest,pg,pgt,pot,potex,pottarg,pottargex,grad,gradex,
       gradtarg,gradtargex,&err);

  cprind("l2 rel error=",&err,1);

  cprin_skipline(2);
  cprin_message("========");












  return 0;
}  
