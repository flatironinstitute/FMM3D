  #include "utils.h"

  void lfmm3dpartstoscp_(double *eps, int *nsource,
                    double *source, double *charge, double *pot);


  void lfmm3dpartstoscg_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad);


  void lfmm3dpartstosdp_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot);


  void lfmm3dpartstosdg_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad);


  void lfmm3dpartstoscdp_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot);


  void lfmm3dpartstoscdg_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad);


  void lfmm3dpartstotcp_(double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstotcg_(double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstotdp_(double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstotdg_(double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstotcdp_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, int *nt,
                    double *targ, double *pottarg);


  void lfmm3dpartstotcdg_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg);



  void lfmm3dpartstostcp_(double *eps, int *nsource,
                    double *source, double *charge, double *pot,  
                    int *nt, double *targ, double *pottarg);


  void lfmm3dpartstostcg_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int *nt, double *targ, double *pottarg, double *gradtarg);


  void lfmm3dpartstostdp_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, 
                    int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstostdg_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int *nt, double *targ,
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstostcdp_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int *nt, double *targ, double *pottarg);


  void lfmm3dpartstostcdg_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *nt, 
                    double *targ, double *pottarg, double *gradtarg);


  void lfmm3dpartstoscp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot);


  void lfmm3dpartstoscg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad);


  void lfmm3dpartstosdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot);


  void lfmm3dpartstosdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad);


  void lfmm3dpartstoscdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot);


  void lfmm3dpartstoscdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad);



  void lfmm3dpartstotcp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstotcg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstotdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstotdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstotcdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, int *nt,
                    double *targ, double *pottarg);


  void lfmm3dpartstotcdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg);



  void lfmm3dpartstostcp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot,  
                    int *nt, double *targ, double *pottarg);


  void lfmm3dpartstostcg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int *nt, double *targ, double *pottarg, double *gradtarg);


  void lfmm3dpartstostdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, 
                    int *nt, double *targ, 
                    double *pottarg);


  void lfmm3dpartstostdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int *nt, double *targ,
                    double *pottarg, double *gradtarg);


  void lfmm3dpartstostcdp_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int *nt, double *targ, double *pottarg);


  void lfmm3dpartstostcdg_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *nt, 
                    double *targ, double *pottarg, double *gradtarg);


 void l3ddirectcp_(int *nd, double *source, double *charge, int *ns, 
          double *targ, int *nt, double *pot, double *thresh);


 void l3ddirectcg_(int *nd, double *source, double *charge, int *ns, 
          double *targ, int *nt, double *pot, double* grad, double *thresh);



 void l3ddirectdp_(int *nd, double *source, double *dipvec, int *ns, 
          double *targ, int *nt, double *pot, double *thresh);


 void l3ddirectdg_(int *nd, double *source, double *dipvec, int *ns, 
          double *targ, int *nt, double *pot, double* grad, double *thresh);



 void l3ddirectcdp_(int *nd, double *source, double *charge, 
          double *dipvec, int *ns, double *targ, int *nt, double *pot, 
          double *thresh);


 void l3ddirectcdg_(int *nd, double *source, double *charge,
          double *dipvec, int *ns, double *targ, int *nt, double *pot, 
          double* grad, double *thresh);


