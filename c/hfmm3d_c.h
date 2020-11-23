  #include "utils.h"

  void hfmm3d_s_c_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, int *ier);


  void hfmm3d_s_c_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad, int *ier);


  void hfmm3d_s_d_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, int *ier);


  void hfmm3d_s_d_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad, int *ier);


  void hfmm3d_s_cd_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot, int *ier);


  void hfmm3d_s_cd_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *ier);


  void hfmm3d_t_c_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_t_c_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_t_d_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_t_d_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_t_cd_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, int *nt,
                    double *targ, CPX *pottarg, int *ier);


  void hfmm3d_t_cd_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int *ier);



  void hfmm3d_st_c_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int *nt, double *targ, CPX *pottarg, int *ier);


  void hfmm3d_st_c_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int *nt, double *targ, CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_st_d_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_st_d_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_st_cd_p_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int *nt, double *targ, CPX *pottarg, int *ier);


  void hfmm3d_st_cd_g_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_s_c_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, int *ier);


  void hfmm3d_s_c_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad, int *ier);


  void hfmm3d_s_d_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, int *ier);


  void hfmm3d_s_d_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad, int *ier);


  void hfmm3d_s_cd_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot, int *ier);


  void hfmm3d_s_cd_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *ier);



  void hfmm3d_t_c_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_t_c_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_t_d_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_t_d_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_t_cd_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, int *nt,
                    double *targ, CPX *pottarg, int *ier);


  void hfmm3d_t_cd_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int *ier);



  void hfmm3d_st_c_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int *nt, double *targ, CPX *pottarg, int *ier);


  void hfmm3d_st_c_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int *nt, double *targ, CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_st_d_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int *nt, double *targ, 
                    CPX *pottarg, int *ier);


  void hfmm3d_st_d_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg, int *ier);


  void hfmm3d_st_cd_p_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int *nt, double *targ, CPX *pottarg, int *ier);


  void hfmm3d_st_cd_g_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg, int *ier);


 void h3ddirectcp_(int *nd, CPX *zk, double *source, CPX *charge, int *ns, 
          double *targ, int *nt, CPX *pot, double *thresh);


 void h3ddirectcg_(int *nd, CPX *zk, double *source, CPX *charge, int *ns, 
          double *targ, int *nt, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectdp_(int *nd, CPX *zk, double *source, CPX *dipvec, int *ns, 
          double *targ, int *nt, CPX *pot, double *thresh);


 void h3ddirectdg_(int *nd, CPX *zk, double *source, CPX *dipvec, int *ns, 
          double *targ, int *nt, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectcdp_(int *nd, CPX *zk, double *source, CPX *charge, 
          CPX *dipvec, int *ns, double *targ, int *nt, CPX *pot, 
          double *thresh);


 void h3ddirectcdg_(int *nd, CPX *zk, double *source, CPX *charge,
          CPX *dipvec, int *ns, double *targ, int *nt, CPX *pot, 
          CPX* grad, double *thresh);


