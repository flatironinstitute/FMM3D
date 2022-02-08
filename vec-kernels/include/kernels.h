#ifndef _KERNELS_H_
#define _KERNELS_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void helm3d_f_(const int64_t* nd, const float* zk, const float* sources, const float* charge, const int64_t* ns, const float* ztarg, const int64_t* nt, float* pot, const float* thresh);

void helm3d_vec_f_(const int64_t* nd, const float* zk, const float* sources, const float* charge, const int64_t* ns, const float* ztarg, const int64_t* nt, float* pot, const float* thresh);

void helm3d_d_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void helm3d_vec_d_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void h3ddirectcp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void h3ddirectcg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void h3ddirectch_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

void h3ddirectdp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void h3ddirectdg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void h3ddirectdh_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

void h3ddirectcdp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void h3ddirectcdg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void h3ddirectcdh_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

void l3ddirectcp_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void l3ddirectcg_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void l3ddirectch_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

void l3ddirectdp_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void l3ddirectdg_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void l3ddirectdh_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

void l3ddirectcdp_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh);

void l3ddirectcdg_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh);

void l3ddirectcdh_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh);

#ifdef __cplusplus
}
#endif

#endif //_KERNELS_H_
