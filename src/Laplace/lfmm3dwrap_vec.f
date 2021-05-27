c
c  This file contains the Laplace FMM wrappers
c
c  The laplace FMM evaluates the following potential
c  and its gradient
c    
c     u_{l}(x) = \sum_{j=1}^{N} c_{l,j}/|x-x_{j}| - 
c         v_{l,j}.\nabla (1/|x-x_{j}|)
c
c  Here x_{j} are the source locations, c_{l,j} are the charge strengths,
c  v_{l,j} are the dipole strengths. We refer to the collection of 
c  $x$ at which the potential and its gradient are evaluated as the
c  evaluation points
c
c  This file contains the wrappers to the vectorized routines
c
c  The subroutine names take the following form:
c    lfmm3d_<eval-pts>_<int-ker>_<out>_vec
c
c      <eval-pts>: evaluation points (sources/targets/sources+targets)
c        s: source locations
c        t: target locations
c        st: source and target locations
c
c      <int-ker>: kernel of interaction (charges/dipoles/both)
c        c: charges
c        d: dipoles
c        cd: charges+dipoles
c 
c      <out>: flag for potential/potential+gradient
c        p: potentials
c        g: potentials+gradients
c
c  The wrappers are:
c
c
c  -lfmm3d_s_c_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_c_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_s_c_h_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_s_d_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_d_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_s_d_h_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_s_cd_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_cd_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_s_cd_h_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_t_c_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_c_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_c_h_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_t_d_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_d_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_d_h_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_t_cd_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_cd_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_cd_h_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_st_c_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_c_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_c_h_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_st_d_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_d_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_d_h_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c  -lfmm3d_st_cd_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_cd_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_cd_h_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential, Gradient and Hessians
c-------------------------------------
c

      subroutine lfmm3d_s_c_p_vec(nd,eps,nsource,source,
     1    charge,pot,ier)
      implicit none
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c

      double precision eps
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource)

      double precision pottarg(nd,1)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)


      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_c_g_vec(nd,eps,nsource,source,
     1    charge,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      ier = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)


      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_s_c_h_vec(nd,eps,nsource,source,
     1    charge,pot,grad,hess,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient,
c  and its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision hess(nd,6,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hesstarg(nd,6)

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 0

      ntarg = 0

      ier = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)


      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_s_d_p_vec(nd,eps,nsource,source,
     1    dipvec,pot,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,1)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0


      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_d_g_vec(nd,eps,nsource,source,
     1    dipvec,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)


      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_d_h_vec(nd,eps,nsource,source,
     1    dipvec,pot,grad,hess,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient, and
c  its Hessians
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision hess(nd,6,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hesstarg(nd,6)

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      ntarg = 0

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)


      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_cd_p_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,1)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c
c
c
c
c

      subroutine lfmm3d_s_cd_g_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c
c

      subroutine lfmm3d_s_cd_h_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,grad,hess,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient, and 
c  its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision hess(nd,6,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      ntarg = 0

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_c_p_vec(nd,eps,nsource,source,
     1    charge,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_c_g_vec(nd,eps,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_c_h_vec(nd,eps,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient, and
c  its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)
      double precision hesstarg(nd,6,ntarg)

      double precision hess(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_d_p_vec(nd,eps,nsource,source,
     1    dipvec,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_d_g_vec(nd,eps,nsource,source,
     1    dipvec,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_d_h_vec(nd,eps,nsource,source,
     1    dipvec,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient,
c  and its hessian
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)
      double precision hesstarg(nd,6,ntarg)

      double precision hess(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_cd_p_vec(nd,eps,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_cd_g_vec(nd,eps,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_cd_h_vec(nd,eps,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient,
c  and its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)
      double precision hesstarg(nd,6,ntarg)

      double precision hess(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c

      subroutine lfmm3d_st_c_p_vec(nd,eps,nsource,source,
     1    charge,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_c_g_vec(nd,eps,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c
c

      subroutine lfmm3d_st_c_h_vec(nd,eps,nsource,source,
     1    charge,pot,grad,hess,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient,
c  and its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|}
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6,nsource),hesstarg(nd,6,ntarg)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_st_d_p_vec(nd,eps,nsource,source,
     1    dipvec,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_d_g_vec(nd,eps,nsource,source,
     1    dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_st_d_h_vec(nd,eps,nsource,source,
     1    dipvec,pot,grad,hess,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6,nsource),hesstarg(nd,6,ntarg)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_cd_p_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_cd_g_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_st_cd_h_vec(nd,eps,nsource,source,
     1    charge,dipvec,pot,grad,hess,ntarg,targ,pottarg,
     2    gradtarg,hesstarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad,hess
cf2py  intent(out) pottarg,gradtarg,hesstarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential, its gradient,
c  and its hessian
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j} \frac{1}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    nd: integer
c          number of densities
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double precision(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c
c    -    pot: double precision(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double precision(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    hess: double precision(nd,6,nsource)
c          Hessian at source locations, $\nabla^2 u_{\ell}(x_{j})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    pottarg: double precision(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double precision(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    hesstarg: double precision(nd,6,ntarg)
c          Hessian at target locations, $\nabla^2 u_{\ell}(t_{i})$
c          Hessian is ordered as
c          $u_{xx},u_{yy},u_{zz},u_{xy},u_{xz},u_{yz}$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6,nsource),hesstarg(nd,6,ntarg)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
