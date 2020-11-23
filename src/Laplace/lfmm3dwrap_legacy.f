c   This file contains the legacy FMM wrappers from the previous version
c     of fmmlib3d
c
c   
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials and fields due to
c     point charges and dipoles.  (FORTRAN 90 VERSION)
c
c     lfmm3dpart - Laplace FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm3dpartself - Laplace FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm3dparttarg - Laplace FMM in R^3: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     l3dpartdirect - Laplace interactions in R^3:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c
c
c
        subroutine lfmm3dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
        implicit real *8 (a-h,o-z)
c              
c              
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (1/r) for the Green's function, without the 
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
        complex *16 charge(1)
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
        call lfmm3dparttarg(ier,iprec,nsource,source,
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
        subroutine lfmm3dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
        implicit real *8 (a-h,o-z)
c              
c              
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (1/r) for the Green's function, without the 
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
        complex *16 charge(1)
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
        call lfmm3dparttarg(ier,iprec,nsource,source,
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

      subroutine lfmm3dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarg,targ,ifpottarg,pottarg,
     $     iffldtarg,fldtarg)
      implicit none
c       
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dparttargmain.
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
      integer ier,iprec,nsource
      integer ifcharge,ifdipole,iper
      double precision source(3,nsource)
      
      double complex charge(*),dipstr(*)
      double precision dipvec(3,*)

      integer ifpot,iffld,ifpottarg,iffldtarg
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

      nd = 2
      ier = 0
      call lfmm3d(nd,eps,nsource,source,ifcharge,charge,
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


      subroutine l3dpartdirect(nsource,
     $    source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarg,
     $     targ,ifpottarg,pottarg,iffldtarg,fldtarg)

      implicit none
c
c       Laplace interactions in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not-included.
c   
c       INPUT PARAMETERS:
c
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
c       ntarg: integer:  number of targets
c       targ: real *8 (3,ntarg):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at source locations
c       fld: complex *16 (3,nsource): field (-gradient) at source locations
c       pottarg: complex *16 (ntarg): potential at target locations 
c       fldtarg: complex *16 (3,ntarg): field (-gradient) at target locations 
c
      integer nsource,ifcharge,ifdipole,ifpot,iffld,ntarg
      integer ifpottarg,iffldtarg
      integer nt,ns
      double precision source(3,*), targ(3,*)
      double complex charge(*),dipstr(*)
      double precision, allocatable :: charge_in(:,:)
      double precision dipvec(3,*)
      double complex pot(*),fld(3,*),pottarg(*),fldtarg(3,*)

      double precision, allocatable :: pottmp(:,:),gradtmp(:,:,:)
      double precision, allocatable :: pottargtmp(:,:),
     1   gradtargtmp(:,:,:)

      double precision, allocatable :: dipvec_in(:,:,:)
      double precision thresh

      integer i,j,nd
      integer ifpgh,ifpghtarg

      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision bsize,btmp,sizex,sizey,sizez

      double complex ima
      data ima/(0.0d0,1.0d0)/

      integer ntmp

      nt = ntarg
      ns = nsource

      ifpgh = 0
      ifpghtarg = 0

      if(ifpot.eq.1) ifpgh = 1
      if(iffld.eq.1) ifpgh = 2

      if(ifpottarg.eq.1) ifpghtarg = 1
      if(iffldtarg.eq.1) ifpghtarg = 2

      if(ifpgh.eq.1) allocate(pottmp(2,ns),gradtmp(2,3,1))
      if(ifpgh.eq.2) allocate(pottmp(2,ns),gradtmp(2,3,ns))

      if(ifpghtarg.eq.1) allocate(pottargtmp(2,nt),gradtargtmp(2,3,1))
      if(ifpghtarg.eq.2) allocate(pottargtmp(2,nt),gradtargtmp(2,3,nt))

      if(ifcharge.eq.1) then
        allocate(charge_in(2,ns))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          charge_in(1,i) = real(charge(i))
          charge_in(2,i) = imag(charge(i))
        enddo
C$OMP END PARALLEL DO
        if(ifdipole.ne.1) allocate(dipvec_in(2,3,1))
      endif


      if(ifdipole.eq.1) then
        allocate(dipvec_in(2,3,ns))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          dipvec_in(1,1,i) = real(dipstr(i))*dipvec(1,i)
          dipvec_in(2,1,i) = imag(dipstr(i))*dipvec(1,i)
          dipvec_in(1,2,i) = real(dipstr(i))*dipvec(2,i)
          dipvec_in(2,2,i) = imag(dipstr(i))*dipvec(2,i)
          dipvec_in(1,3,i) = real(dipstr(i))*dipvec(3,i)
          dipvec_in(2,3,i) = imag(dipstr(i))*dipvec(3,i)
        enddo
C$OMP END PARALLEL DO
        if(ifcharge.ne.1) allocate(charge_in(2,1))
      endif


      if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
          pottmp(1,i) = 0
          pottmp(2,i) = 0
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ns
          do j=1,2
            pottmp(j,i) = 0
            gradtmp(j,1,i) = 0
            gradtmp(j,2,i) = 0
            gradtmp(j,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


      if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
          pottargtmp(1,i) = 0
          pottargtmp(2,i) = 0
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,nt
          do j=1,2
            pottargtmp(j,i) = 0
            gradtargtmp(j,1,i) = 0
            gradtargtmp(j,2,i) = 0
            gradtargtmp(j,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

c
c      compute threshold based on boxsize
c


      nd = 2

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
      thresh = bsize*1.0d-16

      ntmp = 1

c
c      compute potential at source
c  
      if(ifcharge.eq.1.and.ifdipole.eq.0) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectcp(nd,source,charge_in,ns,
     1            source(1,i),1,pottmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectcg(nd,source,charge_in,ns,
     1            source(1,i),1,pottmp(1,i),gradtmp(1,1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectcp(nd,source,charge_in,ns,
     1            targ(1,i),1,pottargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectcg(nd,source,charge_in,ns,
     1            targ(1,i),1,pottargtmp(1,i),
     2            gradtargtmp(1,1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectdp(nd,source,dipvec_in,ns,
     1            source(1,i),1,pottmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectdg(nd,source,dipvec_in,ns,
     1            source(1,i),1,pottmp(1,i),gradtmp(1,1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectdp(nd,source,dipvec_in,ns,
     1            targ(1,i),1,pottargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectdg(nd,source,dipvec_in,ns,
     1            targ(1,i),1,pottargtmp(1,i),
     2            gradtargtmp(1,1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then

        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectcdp(nd,source,charge_in,dipvec_in,ns,
     1            source(1,i),1,pottmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,ns
            call l3ddirectcdg(nd,source,charge_in,dipvec_in,ns,
     1            source(1,i),1,pottmp(1,i),gradtmp(1,1,i),thresh)
          enddo
C$OMP END PARALLEL DO
        endif
        

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectcdp(nd,source,charge_in,dipvec_in,ns,
     1            targ(1,i),ntmp,pottargtmp(1,i),thresh)
          enddo
C$OMP END PARALLEL DO
          stop
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do i=1,nt
            call l3ddirectcdg(nd,source,charge_in,dipvec_in,ns,
     1            targ(1,i),1,pottargtmp(1,i),
     2            gradtargtmp(1,1,i),thresh)
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
           pot(i) = pottmp(1,i)+ima*(pottmp(2,i))
        enddo
C$OMP END PARALLEL DO
      endif
 
      if(iffld.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ns
           fld(1,i) = -(gradtmp(1,1,i)+ima*gradtmp(2,1,i))
           fld(2,i) = -(gradtmp(1,2,i)+ima*gradtmp(2,2,i))
           fld(3,i) = -(gradtmp(1,3,i)+ima*gradtmp(2,3,i))
        enddo
C$OMP END PARALLEL DO
      endif
 
      
      if(ifpottarg.eq.1) then 
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
           pottarg(i) = pottargtmp(1,i)+ima*pottargtmp(2,i)
        enddo
C$OMP END PARALLEL DO
      endif
 
      if(iffldtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nt
           fldtarg(1,i) = -(gradtargtmp(1,1,i)+ima*gradtargtmp(2,1,i))
           fldtarg(2,i) = -(gradtargtmp(1,2,i)+ima*gradtargtmp(2,2,i))
           fldtarg(3,i) = -(gradtargtmp(1,3,i)+ima*gradtargtmp(2,3,i))
        enddo
C$OMP END PARALLEL DO
      endif
 
      return
      end
       
      
