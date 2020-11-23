  #include "utils.h"

  void lfmm3d_s_c_p_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, int *ier);


  void lfmm3d_s_c_g_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad, int *ier);

  void lfmm3d_s_c_h_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int* ier);


  void lfmm3d_s_d_p_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, int *ier);


  void lfmm3d_s_d_g_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, int *ier);

  void lfmm3d_s_d_h_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, double *hess, int* ier);


  void lfmm3d_s_cd_p_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot, int *ier);


  void lfmm3d_s_cd_g_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *ier);

  void lfmm3d_s_cd_h_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int* ier);


  void lfmm3d_t_c_p_(double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_t_c_g_(double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_c_h_(double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);


  void lfmm3d_t_d_p_(double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_t_d_g_(double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_d_h_(double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg, double *hesstarg, int *ier);


  void lfmm3d_t_cd_p_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, int *nt,
                    double *targ, double *pottarg, int *ier);


  void lfmm3d_t_cd_g_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_cd_h_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);


  void lfmm3d_st_c_p_(double *eps, int *nsource,
                    double *source, double *charge, double *pot,  
                    int *nt, double *targ, double *pottarg, int *ier);


  void lfmm3d_st_c_g_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int *nt, double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_c_h_(double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int *nt, double *targ, double *pottarg, 
                    double *gradtarg, double *hesstarg, int *ier);


  void lfmm3d_st_d_p_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, 
                    int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_st_d_g_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int *nt, double *targ,
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_d_h_(double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    double *hess, int *nt, double *targ,
                    double *pottarg, double *gradtarg, 
                    double *hesstarg, int *ier);


  void lfmm3d_st_cd_p_(double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int *nt, double *targ, double *pottarg, int *ier);


  void lfmm3d_st_cd_g_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *nt, 
                    double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_cd_h_(double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);


  void lfmm3d_s_c_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, int *ier);


  void lfmm3d_s_c_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad, int *ier);

  void lfmm3d_s_c_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int* ier);


  void lfmm3d_s_d_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, int *ier);


  void lfmm3d_s_d_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, int *ier);

  void lfmm3d_s_d_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, double *hess, int* ier);


  void lfmm3d_s_cd_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot, int *ier);


  void lfmm3d_s_cd_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *ier);

  void lfmm3d_s_cd_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int* ier);


  void lfmm3d_t_c_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_t_c_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_c_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, int *nt, double *targ, 
                    double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);


  void lfmm3d_t_d_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_t_d_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_d_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, int *nt, double *targ,
                    double *pottarg, double *gradtarg, double *hesstarg, int *ier);


  void lfmm3d_t_cd_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, int *nt,
                    double *targ, double *pottarg, int *ier);


  void lfmm3d_t_cd_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_t_cd_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, int *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);


  void lfmm3d_st_c_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot,  
                    int *nt, double *targ, double *pottarg, int *ier);


  void lfmm3d_st_c_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int *nt, double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_c_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int *nt, double *targ, double *pottarg, 
                    double *gradtarg, double *hesstarg, int *ier);


  void lfmm3d_st_d_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, 
                    int *nt, double *targ, 
                    double *pottarg, int *ier);


  void lfmm3d_st_d_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int *nt, double *targ,
                    double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_d_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    double *hess, int *nt, double *targ,
                    double *pottarg, double *gradtarg, 
                    double *hesstarg, int *ier);


  void lfmm3d_st_cd_p_vec_(int *nd, double *eps, int *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int *nt, double *targ, double *pottarg, int *ier);


  void lfmm3d_st_cd_g_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int *nt, 
                    double *targ, double *pottarg, double *gradtarg, int *ier);

  void lfmm3d_st_cd_h_vec_(int *nd, double *eps, int *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int *ier);





 void l3ddirectcp_(int *nd, double *source, double *charge, int *ns, 
          double *targ, int *nt, double *pot, double *thresh);


 void l3ddirectcg_(int *nd, double *source, double *charge, int *ns, 
          double *targ, int *nt, double *pot, double* grad, double *thresh);

 void l3ddirectch_(int *nd, double *source, double *charge, int *ns, 
          double *targ, int *nt, double *pot, double* grad, 
          double *hess, double *thresh);



 void l3ddirectdp_(int *nd, double *source, double *dipvec, int *ns, 
          double *targ, int *nt, double *pot, double *thresh);


 void l3ddirectdg_(int *nd, double *source, double *dipvec, int *ns, 
          double *targ, int *nt, double *pot, double* grad, double *thresh);

 void l3ddirectdh_(int *nd, double *source, double *dipvec, int *ns, 
          double *targ, int *nt, double *pot, double* grad, 
          double *hess, double *thresh);



 void l3ddirectcdp_(int *nd, double *source, double *charge, 
          double *dipvec, int *ns, double *targ, int *nt, double *pot, 
          double *thresh);


 void l3ddirectcdg_(int *nd, double *source, double *charge,
          double *dipvec, int *ns, double *targ, int *nt, double *pot, 
          double* grad, double *thresh);


 void l3ddirectcdh_(int *nd, double *source, double *charge,
          double *dipvec, int *ns, double *targ, int *nt, double *pot, 
          double* grad, double *hess, double *thresh);

