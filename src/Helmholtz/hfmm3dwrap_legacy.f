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
c
c
c
        subroutine hfmm3dpart(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
        implicit real *8 (a-h,o-z)
c              
c              
c       Helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (exp(ikr)/r) for the Green's function, without the 
c       (1/4 pi) scaling. Self-interactions are not included.
c   
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
c
        dimension source(3,1)
        complex *16 charge(1),zk
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        dimension w(1)
c
        dimension targ(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarg=0
        ifpottarg=0
        iffldtarg=0
c
        call hfmm3dparttarg(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)
c
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
        implicit real *8 (a-h,o-z)
c              
c              
c       Helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (exp(ikr)/r) for the Green's function, without the 
c       (1/4 pi) scaling. Self-interactions are not included.
c   
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
c
        dimension source(3,1)
        complex *16 charge(1),zk
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        dimension w(1)
c
        dimension targ(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarg=0
        ifpottarg=0
        iffldtarg=0
c
        call hfmm3dparttarg(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarg,targ,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        return
        end
c
c
c
c
c

      subroutine hfmm3dparttarg(ier,iprec,zk,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,
     $     iffldtarg,fldtarg)
      implicit none
c       
c       Helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use (exp(ikr)/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine hfmm3dparttargmain.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c       dipvec: real *8 (3,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code (currently unused)
c
c       pot: complex *16 (nsource): potential at source locations
c       fld: complex *16 (3,nsource): field (-gradient) at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c
      integer ier,iprec,nsource,iper
      integer ifcharge,ifdipole
      double complex zk
      double precision source(3,nsource)
      
      double complex charge(nsource),dipstr(nsource)
      double precision dipvec(3,nsource)

      integer ifpot,iffld,ifpottarg,iffldtarg
      double complex  pot(nsource),fld(3,nsource)
      double complex  pottarg(ntarg),fldtarg(3,ntarg)

      integer nd,ifpgh,ifpghtarg
      integer ntarg
      double precision targ(3,ntarg)
      double complex, allocatable :: dipvec_in(:,:)
      double complex, allocatable :: pottmp(:),gradtmp(:,:)
      double complex, allocatable :: pottargtmp(:),gradtargtmp(:,:)

      double complex hess(6),hesstarg(6)
      double precision eps

      integer i
     
c     set fmm tolerance based on iprec flag.
c
      if( iprec .eq. -2 ) eps=.5d-0 
      if( iprec .eq. -1 ) eps=.5d-1
      if( iprec .eq. 0 ) eps=.5d-2
      if( iprec .eq. 1 ) eps=.5d-3
      if( iprec .eq. 2 ) eps=.5d-6
      if( iprec .eq. 3 ) eps=.5d-9
      if( iprec .eq. 4 ) eps=.5d-12
      if( iprec .eq. 5 ) eps=.5d-15
      if( iprec .eq. 6 ) eps=0


      if(ifdipole.eq.1) then
         allocate(dipvec_in(3,nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
         do i=1,nsource
           dipvec_in(1,i) = dipstr(i)*dipvec(1,i)
           dipvec_in(2,i) = dipstr(i)*dipvec(2,i)
           dipvec_in(3,i) = dipstr(i)*dipvec(3,i)
         enddo
C$OMP END PARALLEL DO        
      endif
      if(ifdipole.ne.1) allocate(dipvec_in(3,1))


      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpgh.eq.1) then
        allocate(pottmp(nsource),gradtmp(3,1))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,nsource
          pottmp(i) = 0
        enddo
C$OMP END PARALLEL DO        
        gradtmp(1,1)= 0
        gradtmp(2,1)= 0
        gradtmp(3,1)= 0
      endif
      if(ifpgh.eq.2) then
        allocate(pottmp(nsource),gradtmp(3,nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,nsource
          pottmp(i) = 0
          gradtmp(1,i)= 0
          gradtmp(2,i)= 0
          gradtmp(3,i)= 0
        enddo
C$OMP END PARALLEL DO        
      endif

      if(ifpottarg.eq.1) ifpghtarg = 1
      if(iffldtarg.eq.1) ifpghtarg = 2

      if(ifpghtarg.eq.1) then
        allocate(pottargtmp(nsource),gradtargtmp(3,1))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,ntarg
          pottargtmp(i) = 0
        enddo
C$OMP END PARALLEL DO        
        gradtargtmp(1,1)= 0
        gradtargtmp(2,1)= 0
        gradtargtmp(3,1)= 0
      endif
      if(ifpghtarg.eq.2) then
        allocate(pottargtmp(nsource),gradtargtmp(3,nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,ntarg
          pottargtmp(i) = 0
          gradtargtmp(1,i)= 0
          gradtargtmp(2,i)= 0
          gradtargtmp(3,i)= 0
        enddo
C$OMP END PARALLEL DO        
      endif

      nd = 1
      ier = 0

      call hfmm3d(nd,eps,zk,nsource,source,ifcharge,charge,
     1  ifdipole,dipvec_in,iper,ifpgh,pottmp,gradtmp,hess,ntarg,
     2  targ,ifpghtarg,pottargtmp,gradtargtmp,hesstarg,ier)

      if(ifpot.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,nsource
          pot(i) = pottmp(i)
        enddo
C$OMP END PARALLEL DO        
      endif
      if(iffld.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,nsource
          fld(1,i) = -gradtmp(1,i)
          fld(2,i) = -gradtmp(2,i)
          fld(3,i) = -gradtmp(3,i)
        enddo
C$OMP END PARALLEL DO        
      endif

      if(ifpottarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,ntarg
          pottarg(i) = pottargtmp(i)
        enddo
C$OMP END PARALLEL DO        
      endif
      if(iffldtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
        do i=1,ntarg
          fldtarg(1,i) = -gradtargtmp(1,i)
          fldtarg(2,i) = -gradtargtmp(2,i)
          fldtarg(3,i) = -gradtargtmp(3,i)
        enddo
C$OMP END PARALLEL DO        
      endif


      return
      end


      subroutine h3dpartdirect(zk,ns,
     $    source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,nt,
     $     targ,ifpottarg,pottarg,iffldtarg,fldtarg)

      implicit none
c
c       Helmholtz interactions in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use (exp(ikr)/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not-included.
c   
c       INPUT PARAMETERS:
c
c       zk: complex *16: Helmholtz parameter
c       ns: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c       dipvec: real *8 (3,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       nt: integer:  number of targets
c       targ: real *8 (3,nt):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at source locations
c       fld: complex *16 (3,nsource): field (-gradient) at source locations
c       pottarg: complex *16 (nt): potential at target locations 
c       fldtarg: complex *16 (3,nt): field (-gradient) at target locations 
c
      double complex zk
      integer ns,ifcharge,ifdipole,ifpot,iffld,nt
      integer ifpottarg,iffldtarg
      double precision source(3,*), targ(3,*)
      double complex charge(*),dipstr(*)
      double precision dipvec(3,*)
      double complex pot(*),fld(3,*),pottarg(*),fldtarg(3,*)

      double complex, allocatable :: pottmp(:),gradtmp(:,:)
      double complex, allocatable :: pottargtmp(:),gradtargtmp(:,:)

      double complex, allocatable :: dipvec_in(:,:)
      double precision thresh

      integer i,j,nd,n
      integer ifpgh,ifpghtarg

      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision bsize,btmp,sizex,sizey,sizez

      ifpgh = 0
      ifpghtarg = 0

      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpottarg.eq.1) ifpghtarg = 1
      if(iffldtarg.eq.1) ifpghtarg = 2

      if(ifpgh.eq.1) allocate(pottmp(ns),gradtmp(3,1))
      if(ifpgh.eq.2) allocate(pottmp(ns),gradtmp(3,ns))

      if(ifpghtarg.eq.1) allocate(pottargtmp(nt),gradtargtmp(3,1))
      if(ifpghtarg.eq.2) allocate(pottargtmp(nt),gradtargtmp(3,nt))


      if(ifdipole.eq.1) then
        allocate(dipvec_in(3,ns))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          dipvec_in(1,i) = dipstr(i)*dipvec(1,i)
          dipvec_in(2,i) = dipstr(i)*dipvec(2,i)
          dipvec_in(3,i) = dipstr(i)*dipvec(3,i)
        enddo
C$OMP END PARALLEL DO
      endif

      

      if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          pottmp(i) = 0
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          pottmp(i) = 0
          gradtmp(1,i) = 0
          gradtmp(2,i) = 0
          gradtmp(3,i) = 0
        enddo
C$OMP END PARALLEL DO
      endif


      if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
          pottargtmp(i) = 0
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
          pottargtmp(i) = 0
          gradtargtmp(1,i) = 0
          gradtargtmp(2,i) = 0
          gradtargtmp(3,i) = 0
        enddo
C$OMP END PARALLEL DO
      endif

c
c      compute threshold based on boxsize
c


      nd = 1

      xmin = source(1,1)
      xmax = source(1,1)
      ymin = source(2,1)
      ymax = source(2,1)
      zmin = source(3,1)
      zmax = source(3,1)

      do i=1,ns
        if(source(1,i).lt.xmin) xmin = source(1,i)
        if(source(1,i).gt.xmax) xmax = source(1,i)
        if(source(2,i).lt.ymin) ymin = source(2,i)
        if(source(2,i).gt.ymax) ymax = source(2,i)
        if(source(3,i).lt.zmin) zmin = source(3,i)
        if(source(3,i).gt.zmax) zmax = source(3,i)
      enddo

      do i=1,nt
        if(targ(1,i).lt.xmin) xmin = targ(1,i)
        if(targ(1,i).gt.xmax) xmax = targ(1,i)
        if(targ(2,i).lt.ymin) ymin = targ(2,i)
        if(targ(2,i).gt.ymax) ymax = targ(2,i)
        if(targ(3,i).lt.zmin) zmin = targ(3,i)
        if(targ(3,i).gt.zmax) zmax = targ(3,i)
      enddo

      sizex = xmax-xmin
      sizey = ymax-ymin
      sizez = zmax-zmin

      bsize = sizex
      if(sizey.gt.bsize) bsize = sizey
      if(sizez.gt.bsize) bsize = sizez

      btmp = sqrt(sizex**2+sizey**2+sizez**2)
      if(bsize/btmp.le.1.0d-16) then
        write(*,*) "Nothing to compute"
        write(*,*) "Sources and targets cannot be distinguished"
        write(*,*) "Exiting"
        return
      endif
      thresh = bsize*2.0d0**(-51)




c
c      compute potential at source
c 
      n = 1
      if(ifcharge.eq.1.and.ifdipole.eq.0) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectcp(nd,zk,source,charge,ns,
     1            source(1,i),n,pottmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectcg(nd,zk,source,charge,ns,
     1            source(1,i),n,pottmp(i),gradtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectcp(nd,zk,source,charge,ns,
     1            targ(1,i),n,pottargtmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectcg(nd,zk,source,charge,ns,
     1            targ(1,i),n,pottargtmp(i),
     2            gradtargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectdp(nd,zk,source,dipvec_in,ns,
     1            source(1,i),n,pottmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectdg(nd,zk,source,dipvec_in,ns,
     1            source(1,i),n,pottmp(i),gradtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectdp(nd,zk,source,dipvec_in,ns,
     1            targ(1,i),n,pottargtmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectdg(nd,zk,source,dipvec_in,ns,
     1            targ(1,i),n,pottargtmp(i),
     2            gradtargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectcdp(nd,zk,source,charge,dipvec_in,ns,
     1            source(1,i),n,pottmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call h3ddirectcdg(nd,zk,source,charge,dipvec_in,ns,
     1            source(1,i),n,pottmp(i),gradtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectcdp(nd,zk,source,charge,dipvec_in,ns,
     1            targ(1,i),n,pottargtmp(i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call h3ddirectcdg(nd,zk,source,charge,dipvec_in,ns,
     1            targ(1,i),n,pottargtmp(i),
     2            gradtargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

      endif

c
c       extract output arrays
c
      
      if(ifpot.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
           pot(i) = pottmp(i)
        enddo
C$OMP END PARALLEL DO
      endif
 
      if(iffld.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
           fld(1,i) = -gradtmp(1,i)
           fld(2,i) = -gradtmp(2,i)
           fld(3,i) = -gradtmp(3,i)
        enddo
C$OMP END PARALLEL DO
      endif
 
      
      if(ifpottarg.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
           pottarg(i) = pottargtmp(i)
        enddo
C$OMP END PARALLEL DO
      endif
 
      if(iffldtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
           fldtarg(1,i) = -gradtargtmp(1,i)
           fldtarg(2,i) = -gradtargtmp(2,i)
           fldtarg(3,i) = -gradtargtmp(3,i)
        enddo
C$OMP END PARALLEL DO
      endif
 
      return
      end
       
      
