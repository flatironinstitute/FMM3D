  #include "utils.h"

  void lfmm3d_s_c_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, int64_t *ier);


  void lfmm3d_s_c_g_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad, int64_t *ier);

  void lfmm3d_s_c_h_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int64_t* ier);


  void lfmm3d_s_d_p_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, int64_t *ier);


  void lfmm3d_s_d_g_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, int64_t *ier);

  void lfmm3d_s_d_h_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t* ier);


  void lfmm3d_s_cd_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, double *pot, int64_t *ier);


  void lfmm3d_s_cd_g_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int64_t *ier);

  void lfmm3d_s_cd_h_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t* ier);


  void lfmm3d_t_c_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_t_c_g_(double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_c_h_(double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);


  void lfmm3d_t_d_p_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_t_d_g_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_d_h_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, double *hesstarg, int64_t *ier);


  void lfmm3d_t_cd_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, int64_t *nt,
                    double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_t_cd_g_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_cd_h_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);


  void lfmm3d_st_c_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot,  
                    int64_t *nt, double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_st_c_g_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int64_t *nt, double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_c_h_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int64_t *nt, double *targ, double *pottarg, 
                    double *gradtarg, double *hesstarg, int64_t *ier);


  void lfmm3d_st_d_p_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, 
                    int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_st_d_g_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_d_h_(double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    double *hess, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, 
                    double *hesstarg, int64_t *ier);


  void lfmm3d_st_cd_p_(double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int64_t *nt, double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_st_cd_g_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_cd_h_(double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);


  void lfmm3d_s_c_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, int64_t *ier);


  void lfmm3d_s_c_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad, int64_t *ier);

  void lfmm3d_s_c_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int64_t* ier);


  void lfmm3d_s_d_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, int64_t *ier);


  void lfmm3d_s_d_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, int64_t *ier);

  void lfmm3d_s_d_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t* ier);


  void lfmm3d_s_cd_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, double *pot, int64_t *ier);


  void lfmm3d_s_cd_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int64_t *ier);

  void lfmm3d_s_cd_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t* ier);


  void lfmm3d_t_c_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_t_c_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_c_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, int64_t *nt, double *targ, 
                    double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);


  void lfmm3d_t_d_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_t_d_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_d_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, double *hesstarg, int64_t *ier);


  void lfmm3d_t_cd_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, int64_t *nt,
                    double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_t_cd_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_t_cd_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);


  void lfmm3d_st_c_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot,  
                    int64_t *nt, double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_st_c_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    int64_t *nt, double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_c_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *pot, double *grad,
                    double *hess, int64_t *nt, double *targ, double *pottarg, 
                    double *gradtarg, double *hesstarg, int64_t *ier);


  void lfmm3d_st_d_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, 
                    int64_t *nt, double *targ, 
                    double *pottarg, int64_t *ier);


  void lfmm3d_st_d_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_d_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *dipvec, double *pot, double *grad,
                    double *hess, int64_t *nt, double *targ,
                    double *pottarg, double *gradtarg, 
                    double *hesstarg, int64_t *ier);


  void lfmm3d_st_cd_p_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double *charge, double *dipvec, double *pot,
                    int64_t *nt, double *targ, double *pottarg, int64_t *ier);


  void lfmm3d_st_cd_g_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg, int64_t *ier);

  void lfmm3d_st_cd_h_vec_(int64_t *nd, double *eps, int64_t *nsource,
                    double *source, double* charge, double *dipvec, double *pot,
                    double *grad, double *hess, int64_t *nt, 
                    double *targ, double *pottarg, double *gradtarg,
                    double *hesstarg, int64_t *ier);





 void l3ddirectcp_(int64_t *nd, double *source, double *charge, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double *thresh);


 void l3ddirectcg_(int64_t *nd, double *source, double *charge, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double* grad, double *thresh);

 void l3ddirectch_(int64_t *nd, double *source, double *charge, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double* grad, 
          double *hess, double *thresh);



 void l3ddirectdp_(int64_t *nd, double *source, double *dipvec, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double *thresh);


 void l3ddirectdg_(int64_t *nd, double *source, double *dipvec, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double* grad, double *thresh);

 void l3ddirectdh_(int64_t *nd, double *source, double *dipvec, int64_t *ns, 
          double *targ, int64_t *nt, double *pot, double* grad, 
          double *hess, double *thresh);



 void l3ddirectcdp_(int64_t *nd, double *source, double *charge, 
          double *dipvec, int64_t *ns, double *targ, int64_t *nt, double *pot, 
          double *thresh);


 void l3ddirectcdg_(int64_t *nd, double *source, double *charge,
          double *dipvec, int64_t *ns, double *targ, int64_t *nt, double *pot, 
          double* grad, double *thresh);


 void l3ddirectcdh_(int64_t *nd, double *source, double *charge,
          double *dipvec, int64_t *ns, double *targ, int64_t *nt, double *pot, 
          double* grad, double *hess, double *thresh);

