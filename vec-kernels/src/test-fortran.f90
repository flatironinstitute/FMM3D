program main
  implicit none
  integer*4 :: nd, Ns, Nt

  nd = 1
  Ns = 1000
  Nt = 1000
  call test(nd, Ns, Nt)
end

subroutine test(nd, Ns, Nt)
  use iso_c_binding
  implicit none
  include 'kernels.f90'
  integer*4 :: i, nd, Ns, Nt
  real*8 :: Xs(Ns*3), Xt(Nt*3), thresh
  complex*16 :: zk, Vs(Ns*nd), Vt(Nt*nd), Vt_ref(Nt*nd)
  double precision :: omp_get_wtime, tt
  real*4 :: rand

  zk = (1.0,1.0)
  thresh = 1.0e-12

  call srand(0)
  do i=1, Ns*3
    Xs(i) = rand(0)
  enddo
  do i=1, Ns
    Vs(i) = rand(0)
  enddo
  do i=1, Nt*3
    Xt(i) = rand(0)
  enddo
  do i=1, Nt
    Vt(i) = (0,0)
    Vt_ref(i) = (0,0)
  enddo

  tt = -omp_get_wtime()
  do i=1,100
    call helm3d_d(nd, zk, Xs, Vs, ns, Xt, nt, Vt_ref, thresh)
  enddo
  tt = tt + omp_get_wtime()
  print*, "Unvectorized : ", tt

  tt = -omp_get_wtime()
  do i=1,100
    call helm3d_vec_d(nd, zk, Xs, Vs, ns, Xt, nt, Vt, thresh)
  enddo
  tt = tt + omp_get_wtime()
  print*, "Vectorized : ", tt

  print*, "Maximum relative error : ", maxval(abs(Vt_ref - Vt)) / maxval(abs(Vt_ref));
endsubroutine

