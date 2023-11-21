  #include "utils.h"

  void hfmm3d_s_c_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, int64_t *ier);


  void hfmm3d_s_c_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad, int64_t *ier);


  void hfmm3d_s_d_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, int64_t *ier);


  void hfmm3d_s_d_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *ier);


  void hfmm3d_s_cd_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot, int64_t *ier);


  void hfmm3d_s_cd_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *ier);


  void hfmm3d_t_c_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_t_c_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, int64_t *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_t_d_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_t_d_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, int64_t *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_t_cd_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, int64_t *nt,
                    double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_t_cd_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, int64_t *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);



  void hfmm3d_st_c_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int64_t *nt, double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_st_c_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int64_t *nt, double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_st_d_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_st_d_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int64_t *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_st_cd_p_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int64_t *nt, double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_st_cd_g_(double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_s_c_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, int64_t *ier);


  void hfmm3d_s_c_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad, int64_t *ier);


  void hfmm3d_s_d_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, int64_t *ier);


  void hfmm3d_s_d_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *ier);


  void hfmm3d_s_cd_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot, int64_t *ier);


  void hfmm3d_s_cd_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *ier);



  void hfmm3d_t_c_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_t_c_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, int64_t *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_t_d_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_t_d_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, int64_t *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_t_cd_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, int64_t *nt,
                    double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_t_cd_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, int64_t *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);



  void hfmm3d_st_c_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int64_t *nt, double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_st_c_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int64_t *nt, double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_st_d_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int64_t *nt, double *targ, 
                    CPX *pottarg, int64_t *ier);


  void hfmm3d_st_d_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int64_t *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int64_t *ier);


  void hfmm3d_st_cd_p_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int64_t *nt, double *targ, CPX *pottarg, int64_t *ier);


  void hfmm3d_st_cd_g_vec_(int64_t *nd, double *eps, CPX *zk, int64_t *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int64_t *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int64_t *ier);


 void h3ddirectcp_(int64_t *nd, CPX *zk, double *source, CPX *charge, int64_t *ns, 
          double *targ, int64_t *nt, CPX *pot, double *thresh);


 void h3ddirectcg_(int64_t *nd, CPX *zk, double *source, CPX *charge, int64_t *ns, 
          double *targ, int64_t *nt, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectdp_(int64_t *nd, CPX *zk, double *source, CPX *dipvec, int64_t *ns, 
          double *targ, int64_t *nt, CPX *pot, double *thresh);


 void h3ddirectdg_(int64_t *nd, CPX *zk, double *source, CPX *dipvec, int64_t *ns, 
          double *targ, int64_t *nt, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectcdp_(int64_t *nd, CPX *zk, double *source, CPX *charge, 
          CPX *dipvec, int64_t *ns, double *targ, int64_t *nt, CPX *pot, 
          double *thresh);


 void h3ddirectcdg_(int64_t *nd, CPX *zk, double *source, CPX *charge,
          CPX *dipvec, int64_t *ns, double *targ, int64_t *nt, CPX *pot, 
          CPX* grad, double *thresh);


