#include <kernels.h>
#include <template-kernels.hpp>
#define VECDIM 4

#ifdef __AVX512F__
#undef VECDIM
#define VECDIM 8
#endif

#ifdef __cplusplus
extern "C" {
#endif

void helm3d_f_(const int64_t* nd, const float* zk, const float* sources, const float* charge, const int64_t* ns, const float* ztarg, const int64_t* nt, float* pot, const float* thresh) {
  h3ddirectcp_cpp<float>(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh);
}

void helm3d_vec_f_(const int64_t* nd, const float* zk, const float* sources, const float* charge, const int64_t* ns, const float* ztarg, const int64_t* nt, float* pot, const float* thresh) {
  h3ddirectcp_vec_cpp<float,VECDIM>(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh);
}

void helm3d_d_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh) {
  h3ddirectcp_cpp<double>(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh);
}

void helm3d_vec_d_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh) {
  h3ddirectcp_vec_cpp<double,VECDIM>(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh);
}

void h3ddirectcp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh) {
  h3ddirectcp_vec_cpp<double,VECDIM>(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh);
}

void h3ddirectcg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh) {
  h3ddirectcg_vec_cpp<double,VECDIM>(nd, zk, sources, charge, ns, ztarg, nt, pot, grad, thresh);
}

void h3ddirectch_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh) {
  h3ddirectch_vec_cpp<double,VECDIM>(nd, zk, sources, charge, ns, ztarg, nt, pot, grad, hess, thresh);
}

void h3ddirectdp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh) {
  h3ddirectdp_vec_cpp<double,VECDIM>(nd, zk, sources, dipvec, ns, ztarg, nt, pot, thresh);
}

void h3ddirectdg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh) {
  h3ddirectdg_vec_cpp<double,VECDIM>(nd, zk, sources, dipvec, ns, ztarg, nt, pot, grad, thresh);
}

void h3ddirectdh_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh) {
  h3ddirectdh_vec_cpp<double,VECDIM>(nd, zk, sources, dipvec, ns, ztarg, nt, pot, grad, hess, thresh);
}

void h3ddirectcdp_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh) {
  h3ddirectcdp_vec_cpp<double,VECDIM>(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, thresh);
}

void h3ddirectcdg_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh) {
  h3ddirectcdg_vec_cpp<double,VECDIM>(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, grad, thresh);
}

void h3ddirectcdh_cpp_(const int64_t* nd, const double* zk, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh) {
  h3ddirectcdh_vec_cpp<double,VECDIM>(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, grad, hess, thresh);
}

void l3ddirectcp_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh){
  l3ddirectcp_vec_cpp<double, VECDIM>(nd, sources, charge, ns, ztarg, nt, pot, thresh);
}

void l3ddirectcg_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh){
  l3ddirectcg_vec_cpp<double, VECDIM>(nd, sources, charge, ns, ztarg, nt, pot, grad, thresh);
}

void l3ddirectch_cpp_(const int64_t* nd, const double* sources, const double* charge, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh){
  l3ddirectch_vec_cpp<double, VECDIM>(nd, sources, charge, ns, ztarg, nt, pot, grad, hess, thresh);
}

void l3ddirectdp_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh){
  l3ddirectdp_vec_cpp<double, VECDIM>(nd, sources, dipvec, ns, ztarg, nt, pot, thresh);
}

void l3ddirectdg_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh){
  l3ddirectdg_vec_cpp<double, VECDIM>(nd, sources, dipvec, ns, ztarg, nt, pot, grad, thresh);
}

void l3ddirectdh_cpp_(const int64_t* nd, const double* sources, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh){
  l3ddirectdh_vec_cpp<double, VECDIM>(nd, sources, dipvec, ns, ztarg, nt, pot, grad, hess, thresh);
}

void l3ddirectcdp_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, const double* thresh){
  l3ddirectcdp_vec_cpp<double, VECDIM>(nd, sources, charge, dipvec, ns, ztarg, nt, pot, thresh);
}

void l3ddirectcdg_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, const double* thresh){
  l3ddirectcdg_vec_cpp<double, VECDIM>(nd, sources, charge, dipvec, ns, ztarg, nt, pot, grad, thresh);
}

void l3ddirectcdh_cpp_(const int64_t* nd, const double* sources, const double* charge, const double* dipvec, const int64_t* ns, const double* ztarg, const int64_t* nt, double* pot, double* grad, double* hess, const double* thresh){
  l3ddirectcdh_vec_cpp<double, VECDIM>(nd, sources, charge, dipvec, ns, ztarg, nt, pot, grad, hess, thresh);
}

#ifdef __cplusplus
}
#endif

