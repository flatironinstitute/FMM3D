c   These file contains the legacy FMM wrappers from the previous version
c     of fmmlib3d
c
c   
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Helmholtz potentials and fields due to
c     point charges and dipoles.  (FORTRAN 90 VERSION)
c
c     hfmm3dpart - Helmholtz FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     hfmm3dpartself - Helmholtz FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     hfmm3dparttarg - Helmholtz FMM in R^3: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     h3dpartdirect - Helmholtz interactions in R^3:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c

      subroutine hfmm3dpart(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
      implicit none

      integer ier,iprec,nsource
      integer ifcharge,ifdipole
      double complex zk
      double precision source(3,nsource)
      
      double complex charge(*),dipstr(*)
      double precision dipvec(3,*)

      integer ifpot,iffld
      double complex  pot(*),fld(3,*)

      integer nd,ifpgh,ifpghtarg
      integer ntarg
      double precision targ(3,1)
      double complex, allocatable :: dipvec_in(:,:)
      double complex, allocatable :: pottmp(:),gradtmp(:,:)

      double complex hess(6),hesstarg(6),pottarg,gradtarg(3)
      double precision eps
     
c     set fmm tolerance based on iprec flag.
c
      if( iprec .eq. -2 ) eps=.5d-0 
      if( iprec .eq. -1 ) eps=.5d-1
      if( iprec .eq. 0 ) eps=.5d-2
      if( iprec .eq. 1 ) eps=.5d-3
      if( iprec .eq. 2 ) eps=.5d-6
      if( iprec .eq. 3 ) eps.5d-9
      if( iprec .eq. 4 ) eps=.5d-12
      if( iprec .eq. 5 ) eps=.5d-15
      if( iprec .eq. 6 ) eps=0

      ntarg = 0
      ifpghtarg = 0


      if(ifdipole.eq.1) then
         allocate(dipvec_in(3,nsource))
         do i=1,nsource
           dipvec_in(1,i) = dipstr(i)*dipvec(1,i)
           dipvec_in(2,i) = dipstr(i)*dipvec(2,i)
           dipvec_in(3,i) = dipstr(i)*dipvec(3,i)
         enddo
      endif
      if(ifdipole.ne.1) allocate(dipvec_in(3))


      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpgh.eq.1) then
        allocate(pottmp(nsource),gradtmp(3))
        do i=1,nsource
          pottmp(i) = 0
        enddo
        gradtmp(1)= 0
        gradtmp(2)= 0
        gradtmp(3)= 0
      endif
      if(ifpgh.eq.2) then
        allocate(pottmp(nsource),gradtmp(3,nsource))
        do i=1,nsource
          pottmp(i) = 0
          gradtmp(1,i)= 0
          gradtmp(2,i)= 0
          gradtmp(3,i)= 0
        enddo
      endif
      nd = 1
      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
      ifdipole,dipvec_in,ifpgh,pottmp,gradtmp,hess,ntarg,
      targ,ifpghtarg,pottarg,gradtarg,hesstarg)

      if(ifpot.eq.1) then
        do i=1,nsource
          pot(i) = pottmp(i)
        enddo
      endif
      if(iffld.eq.1) then
        do i=1,nsource
          fld(1,i) = -grad(1,i)
          fld(2,i) = -grad(2,i)
          fld(3,i) = -grad(3,i)
        enddo
      endif

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartself(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
      implicit none

      integer ier,iprec,nsource
      integer ifcharge,ifdipole
      double complex zk
      double precision source(3,nsource)
      
      double complex charge(*),dipstr(*)
      double precision dipvec(3,*)

      integer ifpot,iffld
      double complex  pot(*),fld(3,*)

      integer nd,ifpgh,ifpghtarg
      integer ntarg
      double precision targ(3,1)
      double complex, allocatable :: dipvec_in(:,:)
      double complex, allocatable :: pottmp(:),gradtmp(:,:)

      double complex hess(6),hesstarg(6),pottarg,gradtarg(3)
      double precision eps
     
c     set fmm tolerance based on iprec flag.
c
      if( iprec .eq. -2 ) eps=.5d-0 
      if( iprec .eq. -1 ) eps=.5d-1
      if( iprec .eq. 0 ) eps=.5d-2
      if( iprec .eq. 1 ) eps=.5d-3
      if( iprec .eq. 2 ) eps=.5d-6
      if( iprec .eq. 3 ) eps.5d-9
      if( iprec .eq. 4 ) eps=.5d-12
      if( iprec .eq. 5 ) eps=.5d-15
      if( iprec .eq. 6 ) eps=0

      ntarg = 0
      ifpghtarg = 0


      if(ifdipole.eq.1) then
         allocate(dipvec_in(3,nsource))
         do i=1,nsource
           dipvec_in(1,i) = dipstr(i)*dipvec(1,i)
           dipvec_in(2,i) = dipstr(i)*dipvec(2,i)
           dipvec_in(3,i) = dipstr(i)*dipvec(3,i)
         enddo
      endif
      if(ifdipole.ne.1) allocate(dipvec_in(3))


      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpgh.eq.1) then
        allocate(pottmp(nsource),gradtmp(3))
        do i=1,nsource
          pottmp(i) = 0
        enddo
        gradtmp(1)= 0
        gradtmp(2)= 0
        gradtmp(3)= 0
      endif
      if(ifpgh.eq.2) then
        allocate(pottmp(nsource),gradtmp(3,nsource))
        do i=1,nsource
          pottmp(i) = 0
          gradtmp(1,i)= 0
          gradtmp(2,i)= 0
          gradtmp(3,i)= 0
        enddo
      endif
      nd = 1
      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
      ifdipole,dipvec_in,ifpgh,pottmp,gradtmp,hess,ntarg,
      targ,ifpghtarg,pottarg,gradtarg,hesstarg)

      if(ifpot.eq.1) then
        do i=1,nsource
          pot(i) = pottmp(i)
        enddo
      endif
      if(iffld.eq.1) then
        do i=1,nsource
          fld(1,i) = -grad(1,i)
          fld(2,i) = -grad(2,i)
          fld(3,i) = -grad(3,i)
        enddo
      endif

      return
      end
c
c
c
c
c

      subroutine hfmm3dparttarg(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,iffldtarg,
     $     fldtarg)
      implicit none

      integer ier,iprec,nsource
      integer ifcharge,ifdipole
      double complex zk
      double precision source(3,nsource)
      
      double complex charge(*),dipstr(*)
      double precision dipvec(3,*)

      integer ifpot,iffld
      double complex  pot(*),fld(3,*)
      double complex  pottarg(*),fldtarg(3,*)

      integer nd,ifpgh,ifpghtarg
      integer ntarg
      double precision targ(3,*)
      double complex, allocatable :: dipvec_in(:,:)
      double complex, allocatable :: pottmp(:),gradtmp(:,:)
      double complex, allocatable :: pottargtmp(:),gradtargtmp(:,:)

      double complex hess(6),hesstarg(6)
      double precision eps
     
c     set fmm tolerance based on iprec flag.
c
      if( iprec .eq. -2 ) eps=.5d-0 
      if( iprec .eq. -1 ) eps=.5d-1
      if( iprec .eq. 0 ) eps=.5d-2
      if( iprec .eq. 1 ) eps=.5d-3
      if( iprec .eq. 2 ) eps=.5d-6
      if( iprec .eq. 3 ) eps.5d-9
      if( iprec .eq. 4 ) eps=.5d-12
      if( iprec .eq. 5 ) eps=.5d-15
      if( iprec .eq. 6 ) eps=0


      if(ifdipole.eq.1) then
         allocate(dipvec_in(3,nsource))
         do i=1,nsource
           dipvec_in(1,i) = dipstr(i)*dipvec(1,i)
           dipvec_in(2,i) = dipstr(i)*dipvec(2,i)
           dipvec_in(3,i) = dipstr(i)*dipvec(3,i)
         enddo
      endif
      if(ifdipole.ne.1) allocate(dipvec_in(3))


      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpgh.eq.1) then
        allocate(pottmp(nsource),gradtmp(3))
        do i=1,nsource
          pottmp(i) = 0
        enddo
        gradtmp(1)= 0
        gradtmp(2)= 0
        gradtmp(3)= 0
      endif
      if(ifpgh.eq.2) then
        allocate(pottmp(nsource),gradtmp(3,nsource))
        do i=1,nsource
          pottmp(i) = 0
          gradtmp(1,i)= 0
          gradtmp(2,i)= 0
          gradtmp(3,i)= 0
        enddo
      endif

      if(ifpottarg.eq.1) ifpghtarg = 1
      if(iffldtarg.eq.1) ifpghtarg = 2

      if(ifpghtarg.eq.1) then
        allocate(pottargtmp(nsource),gradtargtmp(3))
        do i=1,ntarg
          pottargtmp(i) = 0
        enddo
        gradtargtmp(1)= 0
        gradtargtmp(2)= 0
        gradtargtmp(3)= 0
      endif
      if(ifpghtarg.eq.2) then
        allocate(pottargtmp(nsource),gradtargtmp(3,nsource))
        do i=1,ntarg
          pottargtmp(i) = 0
          gradtargtmp(1,i)= 0
          gradtargtmp(2,i)= 0
          gradtargtmp(3,i)= 0
        enddo
      endif

      nd = 1
      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
      ifdipole,dipvec_in,ifpgh,pottmp,gradtmp,hess,ntarg,
      targ,ifpghtarg,pottarg,gradtarg,hesstarg)

      if(ifpot.eq.1) then
        do i=1,nsource
          pot(i) = pottmp(i)
        enddo
      endif
      if(iffld.eq.1) then
        do i=1,nsource
          fld(1,i) = -grad(1,i)
          fld(2,i) = -grad(2,i)
          fld(3,i) = -grad(3,i)
        enddo
      endif

      if(ifpottarg.eq.1) then
        do i=1,ntarg
          pottarg(i) = pottargtmp(i)
        enddo
      endif
      if(iffld.eq.1) then
        do i=1,ntarg
          fldtarg(1,i) = -gradtarg(1,i)
          fldtarg(2,i) = -gradtarg(2,i)
          fldtarg(3,i) = -gradtarg(3,i)
        enddo
      endif


      return
      end
