  #include "utils.h"

  void hfmm3dpartstoscp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot);


  void hfmm3dpartstoscg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *pot, CPX *grad);


  void hfmm3dpartstosdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipstr, double *dipvec, CPX *pot);


  void hfmm3dpartstosdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *dipstr, double *dipvec, CPX *pot,
                    CPX *grad);



  void hfmm3dpartstoscdp_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX *charge, CPX *dipstr, 
                    double *dipvec, CPX *pot);


  void hfmm3dpartstoscdg_(double *eps, CPX *zk, int *nsource,
                    double *source, CPX* charge, CPX *dipstr, 
                    double *dipvec, CPX *pot,
                    CPX *grad);


