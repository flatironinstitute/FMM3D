c
c
c

      subroutine hfmm3dpartstoscp(eps,zk,nsource,source,
     1    charge,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource),dipstr(1)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
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
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource),dipstr(1)
      double complex dipvec(3,1)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstosdp(eps,zk,nsource,source,
     1    dipstr,dipvec,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(1)
      double complex dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstosdg(eps,zk,nsource,source,
     1    dipstr,dipvec,pot,grad)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(1)
      double complex dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstoscdp(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource),dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstoscdg(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,grad)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nsource),dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(1),gradtarg(3,1)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcp(eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(1)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcg(eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(1)
      double complex dipvec(3,1)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdp(eps,zk,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      double complex dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdg(eps,zk,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      double complex dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotcdp(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstotcdg(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(1),grad(3,1)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
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
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(1)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
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
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(1)
      double complex dipvec(3,1)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
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

      subroutine hfmm3dpartstostdp(eps,zk,nsource,source,
     1    dipstr,dipvec,pot,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      double complex dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostdg(eps,zk,nsource,source,
     1    dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(1)
      double complex dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdp(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(nsource)
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

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdg(eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nsource),dipstr(nsource)
      double complex dipvec(3,nsource)

      double complex pot(nsource),grad(3,nsource)
      double complex pottarg(ntarg),gradtarg(3,ntarg)

      double complex hess(6),hesstarg(6)

      nd = 1
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
