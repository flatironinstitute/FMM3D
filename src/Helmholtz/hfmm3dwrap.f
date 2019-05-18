c     this file contains a collection a wrappers for the 
c     Helmholtz FMM. There is a separate wrapper for
c     each kind of interaction, and output requests
c
c     The Helmholtz FMM evaluates the following potential and it's
c     gradient at a collection of sources and targets
c
c     .. math::
c 
c          u(x) = \sum_{j=1}^{N} c_{j} G_{k}(x-x_{j})  - \nabla
c                    G_{k}(x-x_{j} \cdot v_{j} \, .
c     
c     Here $x_{j}$ are the source locations, $c_{j}$ are the 
c     charge strengths and $v_{j}$ are the dipole strengths.
c     
c
c     The subroutine names take the following form:
c       "hfmm3dpart<1><2><3>"
c       <1> = stos - Evaluate $u$ and it's gradient at the source
c                     locations $x_{i}$ 
c           = stot - Evaluate $u$ and it's gradient at $t_{i}$,
c                      a collection of target locations specified
c                      by the user
c           = stost - Evaluate $u$ and it's gradient at 
c                       both source and target locations
c                       $x_{i}$ and $t_{i}$.
c       <2> kind of interaction (charges/dipoles/both)
c           charge interactions = c_{j} G_{k}(x-x_{j})
c           dipole interactions = -\nabla G_{k}(x-x_{j}) \cdot v_{j}
c           = c - charges
c           = d - dipoles
c           = cd - charges + dipoles
c 
c       <3> = p - on output only $u$ is evaluated
c           = g - on output both $u$ and it's gradient are
c                 evaluated
c
c        These are all the single density routines. To get
c        a vectorized version of the routines use 
c        "<subroutine name>_vec"
c    
c    
c
c

      subroutine hfmm3dpartstoscp(eps,zk,nsource,source,
     1    charge,pot)
      implicit none
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot


      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(nsource)
      double complex pottarg(1)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)


      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0


      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstoscg(eps,zk,nsource,source,
     1    charge,pot,grad)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstosdp(eps,zk,nsource,source,
     1    dipvec,pot)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(1)
      
      double complex dipvec(3,nsource)

      double complex pot(nsource)
      double complex pottarg(1)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstosdg(eps,zk,nsource,source,
     1    dipvec,pot,grad)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(out) pot,grad
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(1)
      
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstoscdp(eps,zk,nsource,source,
     1    charge,dipvec,pot)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource)
      double complex pottarg(1)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstoscdg(eps,zk,nsource,source,
     1    charge,dipvec,pot,grad)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(out) pot,grad
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcp(eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(1)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcg(eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdp(eps,zk,nsource,source,
     1    dipvec,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      
      double complex dipvec(3,nsource)

      double complex pot(1)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdg(eps,zk,nsource,source,
     1    dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      
      double complex dipvec(3,nsource)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotcdp(eps,zk,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(1)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcdg(eps,zk,nsource,source,
     1    charge,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcp(eps,zk,nsource,source,
     1    charge,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(nsource)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcg(eps,zk,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,1)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstostdp(eps,zk,nsource,source,
     1    dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      
      double complex dipvec(3,nsource)

      double complex pot(nsource)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostdg(eps,zk,nsource,source,
     1    dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdp(eps,zk,nsource,source,
     1    charge,dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource)
      double complex pottarg(ntarg)

      double complex grad(3),gradtarg(3)
      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdg(eps,zk,nsource,source,
     1    charge,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) nsource,source,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
