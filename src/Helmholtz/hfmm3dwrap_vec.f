c
c  This file contains the Helmholtz FMM wrappers
c
c  The Helmholtz FMM evaluates the following potential
c  and its gradient
c    
c     u_{l}(x) = \sum_{j=1}^{N} c_{l,j}e^{ik|x-x_{j}|}/|x-x_{j}| - 
c         v_{l,j}.\nabla (e^{ik|x-x_{j}|/|x-x_{j}|)
c
c  Here x_{j} are the source locations, c_{l,j} are the charge strengths,
c  v_{l,j} are the dipole strengths. We refer to the collection of 
c  $x$ at which the potential and its gradient are evaluated as the
c  evaluation points
c
c  This file contains the wrappers to the vectorized routines
c
c  The subroutine names take the following form:
c    hfmm3d_<eval-pts>_<int-ker>_<out>_vec
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
c  -hfmm3d_s_c_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_s_c_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_s_d_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_s_d_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_s_cd_p_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_s_cd_g_vec
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_t_c_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_t_c_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_t_d_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_t_d_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_t_cd_p_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_t_cd_g_vec
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_st_c_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_st_c_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_st_d_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_st_d_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -hfmm3d_st_cd_p_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -hfmm3d_st_cd_g_vec
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c
c

      subroutine hfmm3d_s_c_p_vec(nd,eps,zk,nsource,source,
     1    charge,pot,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_s_c_g_vec(nd,eps,zk,nsource,source,
     1    charge,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_s_d_p_vec(nd,eps,zk,nsource,source,
     1    dipvec,pot,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_s_d_g_vec(nd,eps,zk,nsource,source,
     1    dipvec,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_s_cd_p_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,pot,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_s_cd_g_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,pot,grad,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot,grad
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_t_c_p_vec(nd,eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_t_c_g_vec(nd,eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_t_d_p_vec(nd,eps,zk,nsource,source,
     1    dipvec,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_t_d_g_vec(nd,eps,zk,nsource,source,
     1    dipvec,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_t_cd_p_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_t_cd_g_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_st_c_p_vec(nd,eps,zk,nsource,source,
     1    charge,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_st_c_g_vec(nd,eps,zk,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3d_st_d_p_vec(nd,eps,zk,nsource,source,
     1    dipvec,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_st_d_g_vec(nd,eps,zk,nsource,source,
     1    dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = -\sum_{j=1}^{N} v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_st_cd_p_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,pot,ntarg,targ,pottarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c

      subroutine hfmm3d_st_cd_g_vec(nd,eps,zk,nsource,source,
     1    charge,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
cf2py  intent(out) ier
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u_{\ell}(x) = \sum_{j=1}^{N} c_{\ell,j}
c        \frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|} - 
c            v_{\ell,j} \cdot \nabla \left( 
c        \frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\right)
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
c    -    zk: double complex
c          Helmholtz parameter (k)
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double complex(nd,nsource)
c          Charge strengths, $c_{\ell,j}$
c    -    dipvec: double complex(nd,3,nsource)
c          Dipole strengths, $v_{\ell,j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double complex(nd,nsource)
c          Potential at source locations, $u_{\ell}(x_{j})$
c    -    grad: double complex(nd,3,nsource)
c          Gradient at source locations, $\nabla u_{\ell}(x_{j})$
cc    -    pottarg: double complex(nd,ntarg)
c          Potential at target locations, $u_{\ell}(t_{i})$
c    -    gradtarg: double complex(nd,3,ntarg)
c          Gradient at target locations, $\nabla u_{\ell}(t_{i})$
c    -    ier: integer
c          error flag
c           * ier = 0, for normal execution
c           * ier = 4/8, failed to allocate memory in fmm routine
c
c
c
c--------------------------------
c
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd,iper,ier
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(6),hesstarg(6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg,ier)

      return
      end
c
c
c
c
c
