#include <sctl.hpp>
#include <kernels.h>

int main() {
  constexpr int32_t COORD_DIM = 3;

  omp_set_num_threads(1);
  std::cout<<"num of max omp threads: "<<omp_get_max_threads()<<"\n";

  double zk[2] = {1,1};
  double thresh = 1e-12;
  int32_t Ns = 1000, Nt = 1000, nd = 1;

  sctl::Vector<double> Xs(Ns*COORD_DIM), Xt(Nt*COORD_DIM), F(Ns*nd*2), U0(Nt*nd*2), U1(Nt*nd*2);
  for (auto& x : Xs) x = 10*M_PI*mydrand();
  for (auto& x : Xt) x = 10*M_PI*mydrand();
  for (auto& x : F) x = mydrand();
  U0 = 0;
  U1 = 0;

  { // Compute U0, U1
    sctl::Profile::Enable(true);

    sctl::Profile::Tic("Unvectorized");
    for (long i = 0; i < 100; i++) helm3d_d_(&nd, zk, &Xs[0], &F[0], &Ns, &Xt[0], &Nt, &U0[0], &thresh);
    sctl::Profile::Toc();

    sctl::Profile::Tic("Vectorized");
    for (long i = 0; i < 100; i++) helm3d_vec_d_(&nd, zk, &Xs[0], &F[0], &Ns, &Xt[0], &Nt, &U1[0], &thresh);
    sctl::Profile::Toc();

    sctl::Profile::print();
  }

  double max_err = 0, max_val = 0;
  for (long i = 0; i < U0.Dim(); i++) {
    max_err = std::max(max_err, fabs(U0[i]-U1[i]));
    max_val = std::max(max_val, fabs(U0[i]));
  }
  std::cout<<"Relative Error = "<<max_err<<"/"<<max_val<<" = "<<max_err/max_val<<'\n';

  return 0;
}

