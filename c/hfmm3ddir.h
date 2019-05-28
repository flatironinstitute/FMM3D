 #include "utils.h"

 void h3ddirectcp_(int *nd, CPX *zk, double *source, CPX *charge, int *ns, 
          double *targ, int *nd, CPX *pot, double *thresh);


 void h3ddirectcg_(int *nd, CPX *zk, double *source, CPX *charge, int *ns, 
          double *targ, int *nd, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectdp_(int *nd, CPX *zk, double *source, CPX *dipvec, int *ns, 
          double *targ, int *nd, CPX *pot, double *thresh);


 void h3ddirectdg_(int *nd, CPX *zk, double *source, CPX *dipvec, int *ns, 
          double *targ, int *nd, CPX *pot, CPX* grad, double *thresh);



 void h3ddirectcdp_(int *nd, CPX *zk, double *source, CPX *charge, 
          CPX *dipvec, int *ns, double *targ, int *nd, CPX *pot, 
          double *thresh);


 void h3ddirectcdg_(int *nd, CPX *zk, double *source, CPX *charge,
          CPX *dipvec, int *ns, double *targ, int *nd, CPX *pot, 
          CPX* grad, double *thresh);


