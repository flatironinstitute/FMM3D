  #include "utils.h"

  void hfmm3dpartstoscp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot);


  void hfmm3dpartstoscg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad);


  void hfmm3dpartstosdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot);


  void hfmm3dpartstosdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad);


  void hfmm3dpartstoscdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot);


  void hfmm3dpartstoscdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad);



  void hfmm3dpartstotcp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstotcg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstotdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstotdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstotcdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, int *nt,
                    double *targ, CPX *pottarg);


  void hfmm3dpartstotcdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg);



  void hfmm3dpartstostcp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int *nt, double *targ, CPX *pottarg);


  void hfmm3dpartstostcg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int *nt, double *targ, CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstostdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstostdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstostcdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int *nt, double *targ, CPX *pottarg);


  void hfmm3dpartstostcdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstoscp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot);


  void hfmm3dpartstoscg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad);


  void hfmm3dpartstosdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot);


  void hfmm3dpartstosdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot,
                    CPX *grad);


  void hfmm3dpartstoscdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot);


  void hfmm3dpartstoscdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad);



  void hfmm3dpartstotcp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstotcg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, int *nt, double *targ, 
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstotdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstotdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstotcdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, int *nt,
                    double *targ, CPX *pottarg);


  void hfmm3dpartstotcdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg);



  void hfmm3dpartstostcp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot,  
                    int *nt, double *targ, CPX *pottarg);


  void hfmm3dpartstostcg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad,
                    int *nt, double *targ, CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstostdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, 
                    int *nt, double *targ, 
                    CPX *pottarg);


  void hfmm3dpartstostdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipvec, CPX *pot, CPX *grad,
                    int *nt, double *targ,
                    CPX *pottarg, CPX *gradtarg);


  void hfmm3dpartstostcdp_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipvec, CPX *pot,
                    int *nt, double *targ, CPX *pottarg);


  void hfmm3dpartstostcdg_vec_(int *nd, double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipvec, CPX *pot,
                    CPX *grad, int *nt, 
                    double *targ, CPX *pottarg, CPX *gradtarg);


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


