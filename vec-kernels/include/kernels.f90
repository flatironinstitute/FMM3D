interface

  subroutine helm3d_f(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="helm3d_f_")
    implicit none
    integer *4, intent(in) :: nd
    complex *8, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex *8, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex *8, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine helm3d_vec_f(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="helm3d_vec_f_")
    implicit none
    integer *4, intent(in) :: nd
    complex *8, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex *8, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex *8, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine helm3d_d(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="helm3d_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine helm3d_vec_d(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="helm3d_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectcp_vec_d(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectcp_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectcg_vec_d(nd, zk, sources, charge, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectcg_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    complex*16, intent(out):: grad(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectdp_vec_d(nd, zk, sources, dipvec, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectdp_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: dipvec(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectdg_vec_d(nd, zk, sources, dipvec, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectdg_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: dipvec(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    complex*16, intent(out):: grad(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectcdp_vec_d(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectcdp_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    complex*16, intent(in) :: dipvec(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    real    *8, intent(in) :: thresh
  end subroutine

  subroutine h3ddirectcdg_vec_d(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectcdg_vec_d_")
    implicit none
    integer *4, intent(in) :: nd
    complex*16, intent(in) :: zk
    real    *8, intent(in) :: sources(*)
    complex*16, intent(in) :: charge(*)
    complex*16, intent(in) :: dipvec(*)
    integer *4, intent(in) :: ns
    real    *8, intent(in) :: ztarg(*)
    integer *4, intent(in) :: nt
    complex*16, intent(out):: pot(*)
    complex*16, intent(out):: grad(*)
    real    *8, intent(in) :: thresh
  end subroutine

end interface
