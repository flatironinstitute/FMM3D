c
c
c

      subroutine rfmm3dpartstoscp_vec(nd,eps,nsource,source,
     1    charge,pot)
      implicit none
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot

      double precision eps
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      double precision dipstr(nd,1)
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

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstoscg_vec(nd,eps,nsource,source,
     1    charge,pot,grad)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(out) pot,grad
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource)
      double precision dipstr(nd,1)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0


      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)


      return
      end
c
c
c
c
c
c

      subroutine rfmm3dpartstosdp_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,pot)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(out) pot
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
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


      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstosdg_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,pot,grad)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(out) pot,grad
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)


      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstoscdp_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,pot)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(out) pot
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource),dipstr(nd,nsource)
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

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
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

      subroutine rfmm3dpartstoscdg_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,pot,grad)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(out) pot,grad
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double precision charge(nd,nsource),dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,1),gradtarg(nd,3,1)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine rfmm3dpartstotcp_vec(nd,eps,nsource,source,
     1    charge,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,1)
      double precision dipvec(nd,3,1)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine rfmm3dpartstotcg_vec(nd,eps,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,1)
      double precision dipvec(nd,3,1)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstotdp_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstotdg_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstotcdp_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine rfmm3dpartstotcdg_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,1),grad(nd,3,1)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstostcp_vec(nd,eps,nsource,source,
     1    charge,pot,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,1)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstostcg_vec(nd,eps,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,1)
      double precision dipvec(nd,3,1)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine rfmm3dpartstostdp_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstostdg_vec(nd,eps,nsource,source,
     1    dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps
      
      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,1)
      double precision dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstostcdp_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,pot,ntarg,targ,pottarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot
cf2py  intent(out) pottarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource)
      double precision pottarg(nd,ntarg)

      double precision grad(nd,3),gradtarg(nd,3)
      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine rfmm3dpartstostcdg_vec(nd,eps,nsource,source,
     1    charge,dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source,dipstr,dipvec,charge
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,grad
cf2py  intent(out) pottarg,gradtarg
      implicit none
      double precision eps

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double precision charge(nd,nsource),dipstr(nd,nsource)
      double precision dipvec(nd,3,nsource)

      double precision pot(nd,nsource),grad(nd,3,nsource)
      double precision pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double precision hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm3dpart(nd,eps,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
