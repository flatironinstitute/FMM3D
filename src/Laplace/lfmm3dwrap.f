c
c
c

      subroutine lfmm3d_s_c_p(eps,nsource,source,
     1    charge,pot)
      implicit none
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot

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
cf2py  intent(in) eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
      implicit none
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
