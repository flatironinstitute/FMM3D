c  This file contains the Laplace FMM wrappers
c
c  The laplace FMM evaluates the following potential
c  and its gradient
c    
c     u(x) = \sum_{j=1}^{N} c_{j}/|x-x_{j}| - v_{j}.\nabla (1/|x-x_{j}|)
c
c  Here x_{j} are the source locations, c_{j} are the charge strengths,
c  v_{j} are the dipole strengths. We refer to the collection of 
c  $x$ at which the potential and its gradient are evaluated as the
c  evaluation points
c
c  The subroutine names take the following form:
c    lfmm3d_<eval-pts>_<int-ker>_<out>
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
c  -lfmm3d_s_c_p
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_c_g
c    - Evaluation points: Sources
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_s_d_p
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_d_g
c    - Evaluation points: Sources
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_s_cd_p
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_s_cd_g
c    - Evaluation points: Sources
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_c_p
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_c_g
c    - Evaluation points: Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_d_p
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_d_g
c    - Evaluation points: Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_t_cd_p
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_t_cd_g
c    - Evaluation points: Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_c_p
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_c_g
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_d_p
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_d_g
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c  -lfmm3d_st_cd_p
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential
c-------------------------------------
c  -lfmm3d_st_cd_g
c    - Evaluation points: Sources and Targets
c    - Interaction kernel: Charges and Dipoles
c    - Outputs requested: Potential and Gradient
c-------------------------------------
c
c

      subroutine lfmm3d_s_c_p(eps,nsource,source,
     1    charge,pot)
      implicit none
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot

c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c
c
c--------------------------------
c
      double precision eps
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nsource)
      
      double precision dipvec(3,1)

      double precision pot(nsource)

      double precision pottarg(1)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)


      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_c_g(eps,nsource,source,
     1    charge,pot,grad)
      implicit none
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
c
c
c--------------------------------
c
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nsource)
      
      double precision dipvec(3,1)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(1),gradtarg(3,1)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0


      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)


      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_s_d_p(eps,nsource,source,
     1    dipvec,pot)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(nsource)
      double precision pottarg(1)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0


      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_d_g(eps,nsource,source,
     1    dipvec,pot,grad)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot,grad
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(1),gradtarg(3,1)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)


      return
      end
c
c
c
c
c

      subroutine lfmm3d_s_cd_p(eps,nsource,source,
     1    charge,dipvec,pot)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(nsource)
      double precision pottarg(1)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

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

      subroutine lfmm3d_s_cd_g(eps,nsource,source,
     1    charge,dipvec,pot,grad)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot,grad
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source locations $x=x_{j}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(1),gradtarg(3,1)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_c_p(eps,nsource,source,
     1    charge,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,1)

      double precision pot(1)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_c_g(eps,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,1)

      double precision pot(1),grad(3,1)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_d_p(eps,nsource,source,
     1    dipvec,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(1)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_d_g(eps,nsource,source,
     1    dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(1),grad(3,1)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_t_cd_p(eps,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(1)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_t_cd_g(eps,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the target locations $x=t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(1),grad(3,1)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_c_p(eps,nsource,source,
     1    charge,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,1)

      double precision pot(nsource)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_c_g(eps,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|}
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,1)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine lfmm3d_st_d_p(eps,nsource,source,
     1    dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(nsource)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_d_g(eps,nsource,source,
     1    dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = -\sum_{j=1}^{N} v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(1)
      
      double precision dipvec(3,nsource)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_cd_p(eps,nsource,source,
     1    charge,dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
c-------------------------------------
c
c  This subroutine evaluates the potential 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(nsource)
      double precision pottarg(ntarg)

      double precision grad(3),gradtarg(3)
      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine lfmm3d_st_cd_g(eps,nsource,source,
     1    charge,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
c-------------------------------------
c
c  This subroutine evaluates the potential and its gradient 
c      u(x) = \sum_{j=1}^{N} c_{j} \frac{1}{\|x-x_{j}\|} - 
c            v_{j} \cdot \nabla \left( 
c            \frac{1}{\|x-x_{j}\|}\right)
c
c  at the source and target locations $x=x_{j},t_{i}$.
c  When $x=x_{j}$, the term corresponding to $x_{j}$ is 
c  dropped from the sum.
c
c  Input arguments:
c
c    -    eps: double precision
c          precision requested
c    -    nsource: integer
c          Number of sources
c    -    source: double precision(3,nsource)
c          Source locations, $x_{j}$
c    -    charge: double precision(nsource)
c          Charge strengths, $c_{j}$
c    -    dipvec: double precision(3,nsource)
c          Dipole strengths, $v_{j}$
c    -    ntarg: integer
c          Number of targets
c    -    targ: double precision(3,ntarg)
c          Target locations, $t_{i}$
c
c
c  Output arguments:
c
c    -    pot: double precision(nsource)
c          Potential at source locations, $u(x_{j})$
c    -    grad: double precision(3,nsource)
c          Gradient at source locations, $\nabla u(x_{j})$
cc    -    pottarg: double precision(ntarg)
c          Potential at target locations, $u(t_{i})$
c    -    gradtarg: double precision(3,ntarg)
c          Gradient at target locations, $\nabla u(t_{i})$
c
c
c--------------------------------
c
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nsource)
      double precision dipvec(3,nsource)

      double precision pot(nsource),grad(3,nsource)
      double precision pottarg(ntarg),gradtarg(3,ntarg)

      double precision hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
