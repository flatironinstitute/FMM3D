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

  subroutine h3ddirectcp(nd, zk, sources, charge, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectcp_")
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

  subroutine h3ddirectcg(nd, zk, sources, charge, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectcg_")
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

  subroutine h3ddirectdp(nd, zk, sources, dipvec, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectdp_")
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

  subroutine h3ddirectdg(nd, zk, sources, dipvec, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectdg_")
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

  subroutine h3ddirectcdp(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, thresh) !bind(C, name="h3ddirectcdp_")
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

  subroutine h3ddirectcdg(nd, zk, sources, charge, dipvec, ns, ztarg, nt, pot, grad, thresh) !bind(C, name="h3ddirectcdg_")
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
