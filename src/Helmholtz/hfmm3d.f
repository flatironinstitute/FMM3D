c       
c       
c       Generalized helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets 
c
c       We use exp(ikr)/r for the Green's function., without
c       the 1/(4\pi ) scaling.
c
c   
c-----------------------------------------------------------
        subroutine hfmm3d(nd,eps,zk,nsource,source,ifcharge,
     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:    in: integer
c             number of densities
c   
c   eps:   in: double precision
c             requested precision
c
c   zk:    in: double complex
c               helmholtz parameter                
c
c   nsource in: integer  
c                number of sources
c
c   source  in: double precision (3,nsource)
c                source(k,j) is the kth component of the jth
c                source locations
c
c   ifcharge  in: integer  
c             charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c 
c   charge    in: double complex (nd,nsource) 
c              charge strengths
c
c   ifdipole   in: integer
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c
c   dipvec   in: double precision (nd,3,nsource) 
c              dipole orientation vectors
c
c   ifpgh   in: integer
c              flag for evaluating potential/gradient at the sources
c              ifpgh = 1, only potential is evaluated
c              ifpgh = 2, potential and gradients are evaluated
c
c
c   ntarg  in: integer  
c                 number of targs 
c
c   targ  in: double precision (3,ntarg)
c               targ(k,j) is the kth component of the jth
c               targ location
c
c   ifpghtarg   in: integer
c              flag for evaluating potential/gradient at the targs
c              ifpghtarg = 1, only potential is evaluated
c              ifpghtarg = 2, potential and gradient are evaluated
c
c
c     OUTPUT parameters:
c
c   pot:    out: double complex(nd,nsource) 
c               potential at the source locations
c
c   grad:   out: double complex(nd,3,nsource)
c               gradient at the source locations
c
c   hess    out: double complex(nd,6,nsource)
c               hessian at the source locations
c
c   pottarg:    out: double complex(nd,ntarg) 
c               potential at the targ locations
c
c   gradtarg:   out: double complex(nd,3,ntarg)
c               gradient at the targ locations
c
c   hesstarg    out: double complex(nd,6,ntarg)
c                hessian at the target locations
     
c------------------------------------------------------------------

      implicit none

      integer nd

      double complex zk
      double precision eps

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nsource,ntarg

      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource)

      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource),
     1     pottarg(nd,3,ntarg),
     1     gradtarg(nd,3,ntarg),hess(nd,6,*),hesstarg(nd,6,*)

c       Tree variables
      integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
      integer *8 ltree
      integer nlmax
      integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
      integer *8 ipointer(32)
      integer, allocatable :: itree(:)
      double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc      temporary sorted arrays
c
      double precision, allocatable :: sourcesort(:,:),targsort(:,:)
      double precision, allocatable :: radsrc(:)
      double complex, allocatable :: chargesort(:,:)
      double complex, allocatable :: dipvecsort(:,:,:)

      double complex, allocatable :: potsort(:,:),gradsort(:,:,:),
     1       hesssort(:,:,:)
      double complex, allocatable :: pottargsort(:,:),
     1    gradtargsort(:,:,:),hesstargsort(:,:,:)

c
cc       temporary fmm arrays
c
      double precision epsfmm
      integer, allocatable :: nterms(:)
      integer *8, allocatable :: iaddr(:,:)
      double precision, allocatable :: scales(:)
      double precision, allocatable :: rmlexp(:)

      integer lmptemp,nmax
      integer *8 lmptot
      double precision, allocatable :: mptemp(:),mptemp2(:)

c
cc       temporary variables not used in particle code
c
      double precision expc(3),scjsort(1),radexp
      double complex texpssort(100)
      double precision expcsort(3),radssort(1)
      integer ntj,nexpc,nadd

c
cc        other temporary variables
c
       integer i,iert,ifprint,ilev,idim,ier
       double precision time1,time2,omp_get_wtime,second

       
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c      
      ifprint=0



c
cc       figure out tree structure
c
   
c
cc        set criterion for box subdivision
c

       if(eps.ge.0.5d-0) then
         ndiv = 300
       else if(eps.ge.0.5d-1) then
         ndiv = 300
       else if(eps.ge.0.5d-2) then
         ndiv = 300
       else if(eps.ge.0.5d-3) then
         ndiv = 300
       else if(eps.ge.0.5d-6) then
         ndiv = 1000
       else if(eps.ge.0.5d-9) then
         ndiv = 1000
       else if(eps.ge.0.5d-12) then
         ndiv = 1000
       else if(eps.ge.0.5d-15) then
         ndiv = 1000
       else
         ndiv = nsource+ntarg
       endif






c
cc         set tree flags
c
       isep = 1
       nlmax = 200
       nlevels = 0
       nboxes = 0
       mhung = 0
       ltree = 0

       nexpc = 0
       nadd = 0
       ntj = 0

       idivflag = 0

       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       nbmax = 0

       allocate(radsrc(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
       do i=1,nsource
           radsrc(i) = 0
       enddo
C$OMP END PARALLEL DO   


       radexp = 0

c
cc      memory management code for constructing level restricted tree
        iert = 0


        call mklraptreemem(iert,source,nsource,radsrc,targ,ntarg,
     1        expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,
     2        nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3,
     3        mnlist4,mhung,ltree)

        if(ifprint.ge.1) print *, ltree/1.0d9


        if(iert.ne.0) then
           print *, "Error in allocating tree memory"
           stop
        endif


        allocate(itree(ltree))
        allocate(boxsize(0:nlevels))
        allocate(treecenters(3,nboxes))

c       Call tree code
        call mklraptree(source,nsource,radsrc,targ,ntarg,expc,
     1               nexpc,radexp,idivflag,ndiv,isep,mhung,mnbors,
     2               mnlist1,mnlist2,mnlist3,mnlist4,nlevels,
     2               nboxes,treecenters,boxsize,itree,ltree,ipointer)


c     Allocate sorted source and target arrays      

      allocate(sourcesort(3,nsource))
      allocate(targsort(3,ntarg))
      if(ifcharge.eq.1) allocate(chargesort(nd,nsource))
      if(ifdipole.eq.1) then
         allocate(dipvecsort(nd,3,nsource))
      endif

      if(ifpgh.eq.1) then 
        allocate(potsort(nd,nsource),gradsort(nd,3,1),hesssort(nd,6,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource),
     1       hesssort(nd,6,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource),
     1       hesssort(nd,6,nsource))
      else
        allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,6,1))
      endif

      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,1),
     1      hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg),
     1        hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg),
     1        hesstargsort(nd,6,ntarg))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,3,1),
     1     hesstargsort(nd,6,1))
      endif

      

c     scaling factor for multipole and local expansions at all levels
c
      allocate(scales(0:nlevels),nterms(0:nlevels))
      do ilev = 0,nlevels
       scales(ilev) = boxsize(ilev)*abs(zk)
       if(scales(ilev).gt.1) scales(ilev) = 1

cc       scales(ilev) = boxsize(ilev)
      enddo

c
cc      initialize potential and gradient at source
c       locations
c
      if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)

        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


      if(ifpgh.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
            hesssort(idim,4,i) = 0
            hesssort(idim,5,i) = 0
            hesssort(idim,6,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif



c
cc       initialize potential and gradient  at targ
c        locations
c
      if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
            hesstargsort(idim,4,i) = 0
            hesstargsort(idim,5,i) = 0
            hesstargsort(idim,6,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


c     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         call h3dterms(boxsize(i),zk,eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2*nd
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

c
cc       reorder sources
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) call dreorderf(2*nd,nsource,charge,chargesort,
     1                     itree(ipointer(5)))

      if(ifdipole.eq.1) then
         call dreorderf(6*nd,nsource,dipvec,dipvecsort,
     1       itree(ipointer(5)))
      endif

c
cc      reorder targs
c
      call dreorderf(3,ntarg,targ,targsort,itree(ipointer(6)))
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9

      allocate(rmlexp(lmptot),stat=iert)
      if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         stop
      endif


c     Memory allocation is complete. 
c     Call main fmm routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call hfmm3dmain(nd,eps,zk,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,radssort,
     $   iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $   scales,treecenters,itree(ipointer(1)),nterms,
     $   ifpgh,potsort,gradsort,hesssort,ifpghtarg,pottargsort,
     $   gradtargsort,hesstargsort,ntj,texpssort,scjsort)

      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)


      if(ifpgh.eq.1) then
        call dreorderi(2*nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
      endif
      if(ifpgh.eq.2) then 
        call dreorderi(2*nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(6*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
      endif

      if(ifpgh.eq.3) then 
        call dreorderi(2*nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(6*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
        call dreorderi(12*nd,nsource,hesssort,hess,
     1                 itree(ipointer(5)))
      endif


      if(ifpghtarg.eq.1) then
        call dreorderi(2*nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(2*nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(6*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(2*nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(6*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
        call dreorderi(12*nd,ntarg,hesstargsort,hesstarg,
     1     itree(ipointer(6)))
      endif


      return
      end
c
c
      subroutine hfmm3dmain(nd,eps,zk,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipvecsort,
     $     ntarg,targsort,nexpc,expcsort,radssort,
     $     iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,isep,ndiv,nlevels, 
     $     nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $     rscales,centers,laddr,nterms,ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     ntj,jsort,scjsort)
      implicit none

      integer nd
      double precision eps
      integer nsource,ntarg, nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      double complex zk,zk2

      double precision sourcesort(3,nsource)

      double complex chargesort(nd,*)
      double complex dipvecsort(nd,3,*)

      double precision targsort(3,ntarg)

      double complex pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      double complex pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

      integer ntj
      double precision expcsort(3,nexpc)
      double complex jsort(nd,0:ntj,-ntj:ntj,nexpc)


      integer lmptemp
      integer *8 iaddr(2,nboxes), lmptot
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)
       
      double precision timeinfo(10)
      double precision centers(3,nboxes)
c
cc      tree variables
c
      integer isep
      integer *8 ltree
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(32)
      integer itree(ltree)
      integer nboxes
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)
c
cc      pw stuff
c
      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

      integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
      double complex, allocatable :: rdminus2(:,:,:),zeyep(:)
      double complex, allocatable :: rdplus2(:,:,:)
      double precision, allocatable :: zmone(:)
      integer nn,nnn
  
      double complex, allocatable :: rlams(:),whts(:)

      double complex, allocatable :: rlsc(:,:,:)
      integer, allocatable :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:),yshift(:,:),zshift(:,:)

      double complex, allocatable :: fexp(:),fexpback(:)

      double complex, allocatable :: mexp(:,:,:,:)
      double complex, allocatable :: tmp(:,:,:),tmp2(:,:,:)
      double complex, allocatable :: mexpf1(:,:),mexpf2(:,:)
      double complex, allocatable :: mexpp1(:,:),mexpp2(:,:),
     1    mexppall(:,:,:)

      double precision, allocatable :: rsc(:)
      double precision r1

      double precision scjsort(nexpc),radssort(nexpc)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll,idim
      integer ibox,jbox,ilev,npts,npts0
      integer nchild,nlist1,nlist2,nlist3,nlist4

      integer istart,iend,istartt,iendt,istarte,iende
      integer istarts,iends
      integer jstart,jend

      integer ifprint,ifwrite

      integer ifhesstarg
      double precision d,time1,time2,omp_get_wtime

      double precision sourcetmp(3)
      double complex chargetmp(nd)

      integer ix,iy,iz
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)

      double precision thresh

      integer mnbors,mnlist1, mnlist2,mnlist3,mnlist4
      double complex eye, ztmp,zmult
      double precision alphaj
      integer ctr,ifinit2
      double precision, allocatable :: xnodes(:),wts(:)
      double precision radius
      integer nquad2
      integer maX_nodes
      double precision pi
      
      integer istart0,istart1,istartm1,nprin
      double precision rtmp1,rtmp2,rtmp3,rtmp4
      double precision ctmp(3)
      double complex ima

      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/

      integer nlfbox,ier


      ntmax = 1000
      allocate(nfourier(ntmax),nphysical(ntmax))
      allocate(rlams(ntmax),whts(ntmax))


      pi = 4.0d0*atan(1.0d0)

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo

      allocate(rsc(0:nmax))

c
cc     threshold for computing interactions,
c      interactions will be ignored
c      for all pairs of sources and targets
c      which satisfy |r| < thresh
c      where r is the disance between them

      thresh = 2.0d0**(-52)*boxsize(0)
      

      allocate(zeyep(-nmax:nmax),zmone(0:2*nmax))
      
      zeyep(0) = 1
      zmult = -ima
      do i=1,nmax
         zeyep(i) = zeyep(i-1)*zmult
         zeyep(-i) = zeyep(-i+1)/zmult
      enddo


      zmone(0) = 1
      do i=1,2*nmax
         zmone(i) = -zmone(i-1)
      enddo

c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=0
c
c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
      do i=1,nexpc
         do k=-ntj,ntj
           do j = 0,ntj
              do idim=1,nd
                jsort(idim,j,k,i)=0
              enddo
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO

c       
      do i=1,10
        timeinfo(i)=0
      enddo
        
      max_nodes = 10000
      allocate(xnodes(max_nodes))
      allocate(wts(max_nodes))


c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO          
       enddo


c
ccc       set scjsort
c
      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(3)+ibox-1)
            if(nchild.gt.0) then
               istart = itree(ipointer(16)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               do i=istart,iend
                  scjsort(i) = rscales(ilev)
                  radssort(i) = min(radssort(i),boxsize(ilev)/32*
     1                            sqrt(3.0d0))
               enddo
            endif
         enddo
C$OMP END PARALLEL DO
      enddo


c    initialize legendre function evaluation routines
      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)


c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions



      do ilev=2,nlevels
         if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = iend-istart+1

               nchild = itree(ipointer(3)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0) then
                  call h3dformmpc(nd,zk,rscales(ilev),
     1            sourcesort(1,istart),chargesort(1,istart),npts,
     2            centers(1,ibox),nterms(ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
C$OMP END PARALLEL DO            
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = iend-istart+1

               nchild = itree(ipointer(3)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0) then
                  call h3dformmpd(nd,zk,rscales(ilev),
     1            sourcesort(1,istart),
     2            dipvecsort(1,1,istart),npts,
     2            centers(1,ibox),nterms(ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
C$OMP END PARALLEL DO            
         endif

         if(ifdipole.eq.1.and.ifcharge.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)

               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = iend-istart+1

               nchild = itree(ipointer(3)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0) then
                  call h3dformmpcd(nd,zk,rscales(ilev),
     1            sourcesort(1,istart),chargesort(1,istart),
     2            dipvecsort(1,1,istart),npts,
     2            centers(1,ibox),nterms(ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
C$OMP END PARALLEL DO          
         endif
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1


      if(ifprint.ge.1)
     $   call prinf('=== STEP 2 (form lo) ===*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
      do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(26)+ibox-1)
            do i=1,nlist4
               jbox = itree(ipointer(27)+(ibox-1)*mnlist4+i-1)

c              Form local expansion for all boxes in list3
c              of the current box


               istart = itree(ipointer(10)+jbox-1)
               iend = itree(ipointer(11)+jbox-1)
               npts = iend-istart+1
               if(npts.gt.0) then
                  call h3dformtac(nd,zk,rscales(ilev),
     1             sourcesort(1,istart),chargesort(1,istart),npts,
     2             centers(1,ibox),nterms(ilev),
     3             rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo
      endif


      if(ifcharge.eq.0.and.ifdipole.eq.1) then
      do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(26)+ibox-1)
            do i=1,nlist4
               jbox = itree(ipointer(27)+(ibox-1)*mnlist4+i-1)

c              Form local expansion for all boxes in list3
c              of the current box


               istart = itree(ipointer(10)+jbox-1)
               iend = itree(ipointer(11)+jbox-1)
               npts = iend-istart+1
               if(npts.gt.0) then
                   call h3dformtad(nd,zk,rscales(ilev),
     1              sourcesort(1,istart),
     2              dipvecsort(1,1,istart),npts,centers(1,ibox),
     3              nterms(ilev),rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
      do ilev=2,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist4 = itree(ipointer(26)+ibox-1)
            do i=1,nlist4
               jbox = itree(ipointer(27)+(ibox-1)*mnlist4+i-1)

c              Form local expansion for all boxes in list3
c              of the current box


               istart = itree(ipointer(10)+jbox-1)
               iend = itree(ipointer(11)+jbox-1)
               npts = iend-istart+1
               if(npts.gt.0) then
                   call h3dformtacd(nd,zk,rscales(ilev),
     1              sourcesort(1,istart),chargesort(1,istart),
     2              dipvecsort(1,1,istart),npts,centers(1,ibox),
     3              nterms(ilev),rmlexp(iaddr(2,ibox)),wlege,nlege)
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      endif
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c


      do ilev=nlevels-1,0,-1
         nquad2 = nterms(ilev)*2.5
         nquad2 = max(6,nquad2)
         ifinit2 = 1
         call legewhts(nquad2,xnodes,wts,ifinit2)
         radius = boxsize(ilev)/2*sqrt(3.0d0)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = itree(ipointer(10)+jbox-1)
                  iend = itree(ipointer(11)+jbox-1)
                  npts = iend-istart+1

                  if(npts.gt.0) then
                     call h3dmpmp(nd,zk,rscales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),
     4               radius,xnodes,wts,nquad2)
                  endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO          
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1



      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc + mpeval) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels

c
cc       load the necessary quadrature for plane waves
c
      
         zk2 = zk*boxsize(ilev)
         if(real(zk2).le.16*pi.and.imag(zk2).le.12*pi) then
            ier = 0

c
c             get new pw quadrature
c
            
            call hwts3e(ier,eps,zk2,rlams,whts,nlams)
            call hnumfour(eps,zk2,nlams,nfourier)
            call hnumphys(eps,zk2,nlams,nphysical)

            
            nphmax = 0
            nthmax = 0
            nexptotp = 0
            nexptot = 0
            nn = 0
            do i=1,nlams
               nexptotp = nexptotp + nphysical(i)
               nexptot = nexptot + 2*nfourier(i)+1
               nn = nn + nfourier(i)*nphysical(i)
               if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
               if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
            enddo
            allocate(fexp(nn),fexpback(nn))

            allocate(xshift(-5:5,nexptotp))
            allocate(yshift(-5:5,nexptotp))
            allocate(zshift(5,nexptotp))
            allocate(rlsc(0:nterms(ilev),0:nterms(ilev),nlams))
            allocate(tmp(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))
            allocate(tmp2(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))
 
            allocate(mexpf1(nd,nexptot),mexpf2(nd,nexptot),
     1          mexpp1(nd,nexptotp))
            allocate(mexpp2(nd,nexptotp),mexppall(nd,nexptotp,16))


c
cc      NOTE: there can be some memory savings here
c
            bigint = 0
            bigint = nboxes
            bigint = bigint*6
            bigint = bigint*nexptotp*nd

            if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9


            allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
            if(iert.ne.0) then
              print *, "Cannot allocate pw expansion workspace"
              print *, "bigint=", bigint
              stop
            endif


            nn = nterms(ilev)
            allocate(carray(4*nn+1,4*nn+1))
            allocate(dc(0:4*nn,0:4*nn))
            allocate(rdplus(0:nn,0:nn,-nn:nn))
            allocate(rdminus(0:nn,0:nn,-nn:nn))
            allocate(rdsq3(0:nn,0:nn,-nn:nn))
            allocate(rdmsq3(0:nn,0:nn,-nn:nn))

c     generate rotation matrices and carray
            call getpwrotmat(nn,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


            call hrlscini(rlsc,nlams,rlams,rscales(ilev),zk2,
     1         nterms(ilev))
            call hmkexps(rlams,nlams,nphysical,nexptotp,zk2,xshift,
     1           yshift,zshift)
            
            call hmkfexp(nlams,nfourier,nphysical,fexp,fexpback)

c
cc      zero out mexp
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,k)
            do k=1,6
               do i=1,nboxes
                  do j=1,nexptotp
                     do idim=1,nd
                        mexp(idim,j,i,k) = 0.0d0
                     enddo
                  enddo
               enddo
            enddo
C$OMP END PARALLEL DO    



c
cc         compute powers of scaling parameter
c          for rescaling the multipole expansions
c
c          note: the scaling for helmholtz has been eliminated
c         since it is taken care in the scaling of the legendre
c         functions
c
          
cc           r1 = rscales(ilev)
           r1 = 1.0d0
           rsc(0) = 1.0d0
           do i=1,nterms(ilev)
             rsc(i) = rsc(i-1)*r1
           enddo

c
cc         create multipole to plane wave expansion for
c          all boxes at this level
c
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,tmp2)
            do ibox = laddr(1,ilev),laddr(2,ilev)
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = iend - istart+1
               if(npts.gt.0) then

c           rescale multipole expansion
                  call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1               rsc,tmp)
                
                  call hmpoletoexp(nd,tmp,nterms(ilev),
     1                  nlams,nfourier,nexptot,mexpf1,mexpf2,rlsc) 

                  call hftophys(nd,mexpf1,nlams,nfourier,nphysical,
     1                 mexp(1,1,ibox,1),fexp)           

                  call hftophys(nd,mexpf2,nlams,nfourier,nphysical,
     1                 mexp(1,1,ibox,2),fexp)


c             form mexpnorth, mexpsouth for current box

c             Rotate mpole for computing mexpnorth and
c             mexpsouth
                  call rotztoy(nd,nterms(ilev),tmp,
     1                           tmp2,rdminus)

                  call hmpoletoexp(nd,tmp2,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

                  call hftophys(nd,mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,3),fexp)           

                  call hftophys(nd,mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,4),fexp)   


c             Rotate mpole for computing mexpeast, mexpwest
                  call rotztox(nd,nterms(ilev),tmp,
     1                              tmp2,rdplus)
                  call hmpoletoexp(nd,tmp2,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

                  call hftophys(nd,mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,5),fexp)

                  call hftophys(nd,mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,6),fexp)           

               endif
            enddo
C$OMP END PARALLEL DO       
           


c
cc         loop over parent boxes and ship plane wave
c          expansions to the first child of parent 
c          boxes. 
c          The codes are now written from a gathering perspective
c
c          so the first child of the parent is the one
c          recieving all the local expansions
c          coming from all the lists
c
c          
c

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
C$OMP$PRIVATE(mexpf1,mexpf2,mexpp1,mexpp2,mexppall)
C$OMP$PRIVATE(nuall,uall,ndall,dall,nnall,nall,nsall,sall)
C$OMP$PRIVATE(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)
C$OMP$PRIVATE(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)
C$OMP$PRIVATE(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)
C$OMP$PRIVATE(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)
C$OMP$PRIVATE(nw2,w2,nw4,w4,nw6,w6,nw8,w8)
C$OMP$PRIVATE(npts0,nlist3,ctmp)
            do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
           
               npts = 0

               if(ifpghtarg.gt.0) then
                  istart = itree(ipointer(12)+ibox-1)
                  iend = itree(ipointer(13)+ibox-1)
                  npts = npts + iend-istart+1
               endif

               istart = itree(ipointer(14)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               npts = npts + iend-istart+1

               nchild = itree(ipointer(3)+ibox-1)

               if(ifpgh.gt.0) then
                  istart = itree(ipointer(10)+ibox-1)
                  iend = itree(ipointer(11)+ibox-1)
                  npts = npts + iend-istart+1
               endif


               if(npts.gt.0.and.nchild.gt.0) then

              
                  call getpwlistall(ibox,boxsize(ilev),nboxes,
     1            itree(ipointer(18)+ibox-1),itree(ipointer(19)+
     2            mnbors*(ibox-1)),nchild,itree(ipointer(4)),centers,
     3            isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,
     4            neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678,
     5            nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468,
     6            nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,
     7            e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7,
     8            nw2,w2,nw4,w4,nw6,w6,nw8,w8)


                  call hprocessudexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678,
     5            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     6            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     7            xshift,yshift,zshift,fexpback,rlsc)


                  call hprocessnsexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,
     5            ns3478,s3478,ns34,s34,ns78,s78,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8            mexppall(1,1,5),mexppall(1,1,6),mexppall(1,1,7),
     9            mexppall(1,1,8),rdplus,xshift,yshift,zshift,
     9            fexpback,rlsc)

                  call hprocessewexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,
     5            ne3,e3,ne5,e5,ne7,e7,nwall,wall,
     5            nw2468,w2468,nw24,w24,nw68,w68,
     5            nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8            mexppall(1,1,5),mexppall(1,1,6),
     8            mexppall(1,1,7),mexppall(1,1,8),mexppall(1,1,9),
     9            mexppall(1,1,10),mexppall(1,1,11),mexppall(1,1,12),
     9            mexppall(1,1,13),mexppall(1,1,14),mexppall(1,1,15),
     9            mexppall(1,1,16),rdminus,xshift,yshift,zshift,
     9            fexpback,rlsc)
               endif


c
c      handle mp eval
c

            nlist3 = itree(ipointer(24)+ibox-1)
            if(nlist3.gt.0.and.npts.gt.0) then
              call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1             nlist3,itree(ipointer(25)+(ibox-1)*mnlist3),isep,
     2             centers,nuall,uall,ndall,dall,nnall,nall,
     3                     nsall,sall,neall,eall,nwall,wall)
             
              ctmp(1) = centers(1,ibox) - boxsize(ilev)/2
              ctmp(2) = centers(2,ibox) - boxsize(ilev)/2
              ctmp(3) = centers(3,ibox) - boxsize(ilev)/2

              call hprocesslist3udexp(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),
     2             nexptotp,mexp,nuall,uall,ndall,dall,
     3             mexppall(1,1,1),mexppall(1,1,2),
     4             xshift,yshift,zshift)

              call hprocesslist3nsexp(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),
     2             nexptotp,mexp,nnall,nall,nsall,sall,
     3             mexppall(1,1,3),mexppall(1,1,4),
     4             xshift,yshift,zshift)

              call hprocesslist3ewexp(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),
     2             nexptotp,mexp,neall,eall,nwall,wall,
     3             mexppall(1,1,5),mexppall(1,1,6),
     4             xshift,yshift,zshift)

               if(ifpgh.eq.1) then
                 istart = itree(ipointer(10)+ibox-1)
                 iend = itree(ipointer(11)+ibox-1)
                 npts0 = iend-istart+1

                 if(npts0.gt.0) then
                   call hpw_ud_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,1),mexppall(1,1,2),
     3               pot(1,istart))

                   call hpw_ns_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,3),mexppall(1,1,4),
     3               pot(1,istart))

                   call hpw_ew_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,5),mexppall(1,1,6),
     3               pot(1,istart))
                 endif
               endif

               if(ifpgh.eq.2) then
                 istart = itree(ipointer(10)+ibox-1)
                 iend = itree(ipointer(11)+ibox-1)
                 npts0 = iend-istart+1
                 if(npts0.gt.0) then
                   call hpw_ud_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,1),mexppall(1,1,2),
     3               pot(1,istart),grad(1,1,istart))

                   call hpw_ns_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,3),mexppall(1,1,4),
     3               pot(1,istart),grad(1,1,istart))

                   call hpw_ew_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               sourcesort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,5),mexppall(1,1,6),
     3               pot(1,istart),grad(1,1,istart))
                 endif
               endif


               if(ifpghtarg.eq.1) then
                 istart = itree(ipointer(12)+ibox-1)
                 iend = itree(ipointer(13)+ibox-1)
                 npts0 = iend-istart+1
                 if(npts0.gt.0) then
                   call hpw_ud_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,1),mexppall(1,1,2),
     3               pottarg(1,istart))

                   call hpw_ns_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,3),mexppall(1,1,4),
     3               pottarg(1,istart))

                   call hpw_ew_eval_p(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,5),mexppall(1,1,6),
     3               pottarg(1,istart))
                 endif
               endif

               if(ifpghtarg.eq.2) then
                 istart = itree(ipointer(12)+ibox-1)
                 iend = itree(ipointer(13)+ibox-1)
                 npts0 = iend-istart+1
                 if(npts0.gt.0) then
                   call hpw_ud_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,1),mexppall(1,1,2),
     3               pottarg(1,istart),gradtarg(1,1,istart))

                   call hpw_ns_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,3),mexppall(1,1,4),
     3               pottarg(1,istart),gradtarg(1,1,istart))

                   call hpw_ew_eval_g(nd,zk2,ctmp,boxsize(ilev),npts0,
     1               targsort(1,istart),nlams,rlams,whts,nphysical,
     2               nexptotp,nphmax,mexppall(1,1,5),mexppall(1,1,6),
     3               pottarg(1,istart),gradtarg(1,1,istart))
                 endif
               endif
             endif
            enddo
C$OMP END PARALLEL DO        

            deallocate(xshift,yshift,zshift,rlsc,tmp,tmp2)
            deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)

            deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)
            deallocate(fexp,fexpback)

         else
            nquad2 = nterms(ilev)*2.2
            nquad2 = max(6,nquad2)

            ifinit2 = 1
            ier = 0

            call legewhts(nquad2,xnodes,wts,ifinit2)

            radius = boxsize(ilev)/2*sqrt(3.0d0)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nlist2,i,jbox)
            do ibox = laddr(1,ilev),laddr(2,ilev)

               npts = 0
               if(ifpghtarg.gt.0) then
                  istart = itree(ipointer(12)+ibox-1)
                  iend = itree(ipointer(13)+ibox-1)
                  npts = npts + iend - istart + 1
               endif

               istart = itree(ipointer(14)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               npts = npts + iend-istart+1

               if(ifpgh.gt.0) then
                  istart = itree(ipointer(10)+ibox-1)
                  iend = itree(ipointer(11)+ibox-1)
                  npts = npts + iend-istart+1
               endif


               nlist2 = itree(ipointer(22)+ibox-1)
               if(npts.gt.0) then
                  do i =1,nlist2
                     jbox = itree(ipointer(23)+mnlist2*(ibox-1)+i-1)

                     istart = itree(ipointer(10)+jbox-1)
                     iend = itree(ipointer(11)+jbox-1)
                     npts = iend-istart+1

                     if(npts.gt.0) then
                        call h3dmploc(nd,zk,rscales(ilev),
     1                  centers(1,jbox),
     1                  rmlexp(iaddr(1,jbox)),nterms(ilev),
     2                  rscales(ilev),centers(1,ibox),
     2                  rmlexp(iaddr(2,ibox)),nterms(ilev),
     3                  radius,xnodes,wts,nquad2)
                     endif
                  enddo
               endif
           enddo
C$OMP END PARALLEL DO     

c
c
c         handle list 3 interactions at this level
c


           if(ifpgh.eq.1) then         
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
             do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
               nlist3 = itree(ipointer(24)+ibox-1)
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)

               npts = iend-istart+1

               do i=1,nlist3
                 jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                 call h3dmpevalp(nd,zk,rscales(ilev),centers(1,jbox),
     1            rmlexp(iaddr(1,jbox)),nterms(ilev),
     2            sourcesort(1,istart),npts,pot(1,istart),wlege,nlege,
     3            thresh)
               enddo
             enddo
C$OMP END PARALLEL DO          
           endif

           if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
             do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
               nlist3 = itree(ipointer(24)+ibox-1)
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)

               npts = iend-istart+1

               do i=1,nlist3
                 jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                 call h3dmpevalg(nd,zk,rscales(ilev),centers(1,jbox),
     1             rmlexp(iaddr(1,jbox)),nterms(ilev),
     2             sourcesort(1,istart),npts,pot(1,istart),
     3             grad(1,1,istart),wlege,nlege,thresh)
               enddo
             enddo
C$OMP END PARALLEL DO          
           endif

           if(ifpghtarg.eq.1) then         
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
             do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
               nlist3 = itree(ipointer(24)+ibox-1)
               istart = itree(ipointer(12)+ibox-1)
               iend = itree(ipointer(13)+ibox-1)

               npts = iend-istart+1

               do i=1,nlist3
                 jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                 call h3dmpevalp(nd,zk,rscales(ilev),centers(1,jbox),
     1             rmlexp(iaddr(1,jbox)),nterms(ilev),
     2             targsort(1,istart),npts,pottarg(1,istart),
     3             wlege,nlege,thresh)
               enddo
             enddo
C$OMP END PARALLEL DO          
           endif

           if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
             do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
               nlist3 = itree(ipointer(24)+ibox-1)
               istart = itree(ipointer(12)+ibox-1)
               iend = itree(ipointer(13)+ibox-1)

               npts = iend-istart+1

               do i=1,nlist3
                 jbox = itree(ipointer(25)+(ibox-1)*mnlist3+i-1)
                 call h3dmpevalg(nd,zk,rscales(ilev),centers(1,jbox),
     1             rmlexp(iaddr(1,jbox)),nterms(ilev),
     2             targsort(1,istart),npts,pottarg(1,istart),
     3             gradtarg(1,1,istart),wlege,nlege,thresh)
               enddo
             enddo
C$OMP END PARALLEL DO
           endif
         endif
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

        nquad2 = nterms(ilev)*2
        nquad2 = max(6,nquad2)
        ifinit2 = 1
        call legewhts(nquad2,xnodes,wts,ifinit2)
        radius = boxsize(ilev+1)/2*sqrt(3.0d0)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0

            if(ifpghtarg.gt.0) then
               istart = itree(ipointer(12)+ibox-1)
               iend = itree(ipointer(13)+ibox-1)
               npts = npts + iend-istart+1
            endif

            istart = itree(ipointer(14)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)
            npts = npts + iend-istart+1

            if(ifpgh.gt.0) then
               istart = itree(ipointer(10)+ibox-1)
               iend = itree(ipointer(11)+ibox-1)
               npts = npts + iend-istart+1
            endif

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call h3dlocloc(nd,zk,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),rscales(ilev+1),centers(1,jbox),
     3                rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     4                radius,xnodes,wts,nquad2)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(5) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== step 6 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
c

      nquad2 = 2*ntj
      nquad2 = max(6,nquad2)
      ifinit2 = 1

      call legewhts(nquad2,xnodes,wts,ifinit2)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
C

c
cc       shift local expansion to local epxanion at expansion centers
c        (note: this part is not relevant for particle codes.
c        it is relevant only for qbx codes)

      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
C$OMP$SCHEDULE(DYNAMIC)      
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
               istart = itree(ipointer(16)+ibox-1)
               iend = itree(ipointer(17)+ibox-1)
               do i=istart,iend

                  call h3dlocloc(nd,zk,rscales(ilev),
     1             centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2             nterms(ilev),rscales(ilev),expcsort(1,i),
     3             jsort(1,0,-ntj,i),ntj,radssort(i),xnodes,wts,
     4             nquad2)
               enddo
            endif
         enddo
C$OMP END PARALLEL DO
      enddo

c
cc        evaluate local expansion at source and target
c         locations
c
      do ilev = 0,nlevels
        if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(10)+ibox-1)
              iend = itree(ipointer(11)+ibox-1)
              npts = iend-istart+1
              call h3dtaevalp(nd,zk,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO          
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(10)+ibox-1)
              iend = itree(ipointer(11)+ibox-1)
              npts = iend-istart+1
              call h3dtaevalg(nd,zk,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(12)+ibox-1)
              iend = itree(ipointer(13)+ibox-1)
              npts = iend-istart+1
              call h3dtaevalp(nd,zk,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(12)+ibox-1)
              iend = itree(ipointer(13)+ibox-1)
              npts = iend-istart+1

              call h3dtaevalg(nd,zk,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif
      enddo

    
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(6) = time2 - time1


      if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (direct) =====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()

c
cc       directly form local expansions for list1 sources
c        at expansion centers. 
c        (note: this part is not relevant for particle codes.
c         It is relevant only for qbx codes)


      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarte,iende,nlist1,i,jbox)
C$OMP$PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = itree(ipointer(16)+ibox-1)
            iende = itree(ipointer(17)+ibox-1)

            nlist1 = itree(ipointer(20)+ibox-1)
   
            do i =1,nlist1
               jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)


               jstart = itree(ipointer(10)+jbox-1)
               jend = itree(ipointer(11)+jbox-1)

               call hfmm3dexpc_direct(nd,zk,jstart,jend,istarte,
     1         iende,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipvecsort,expcsort,jsort,scjsort,ntj,
     3         wlege,nlege)
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo

c
cc        directly evaluate potential at sources and targets 
c         due to sources in list1

      do ilev=0,nlevels
c
cc           evaluate at the sources
c

        if(ifpgh.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)

              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcp(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO            
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectdp(nd,zk,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcdp(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif
        endif

        if(ifpgh.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)

              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcg(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),thresh)   
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectdg(nd,zk,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = itree(ipointer(10)+ibox-1)
              iends = itree(ipointer(11)+ibox-1)
              npts0 = iends-istarts+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcdg(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif
        endif

        if(ifpghtarg.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)

              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcp(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectdp(nd,zk,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcdp(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif
        endif

        if(ifpghtarg.eq.2) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)

              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcg(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             thresh)   
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectdg(nd,zk,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itree(ipointer(12)+ibox-1)
              iendt = itree(ipointer(13)+ibox-1)
              npts0 = iendt-istartt+1
              nlist1 = itree(ipointer(20)+ibox-1)
              do i=1,nlist1
                jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
                jstart = itree(ipointer(10)+jbox-1)
                jend = itree(ipointer(11)+jbox-1)
                npts = jend-jstart+1
                call h3ddirectcdg(nd,zk,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif
        endif
      enddo
 
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(7) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,7)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)


      return
      end


c------------------------------------------------------------------
      subroutine hfmm3dexpc_direct(nd,zk,istart,iend,jstart,
     $     jend,source,ifcharge,charge,ifdipole,
     $     dipvec,targ,texps,scj,ntj,wlege,nlege)
c---------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the existing local
c     expansions
c
c     INPUT arguments
c------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
c 
c     zk           in: double complex
c                  helmholtz parameter
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to compute the expansions
c 
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: double complex
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipvec      in: double complex(3,ns)
c                 Dipole orientation vector at the source locations
c
c     targ        in: double precision(3,nexpc)
c                 Expansion center locations
c
c     scj         in: double precision(nexpc)
c                 scaling parameters for local expansions
c
c     ntj         in: Integer
c                 Number of terms in expansion
c
c     wlege       in: double precision(0:nlege,0:nlege)
c                 precomputed array of recurrence relation
c                 coeffs for Ynm calculation.
c
c    nlege        in: integer
c                 dimension parameter for wlege
c------------------------------------------------------------
c     OUTPUT
c
c   Updated expansions at the targets
c   texps : coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j, nlege
        integer nd
        integer ifcharge,ifdipole,ier
        double complex zk
        double precision source(3,*)
        double precision wlege(0:nlege,0:nlege)
        double complex charge(nd,*)
        double complex dipvec(nd,3,*)
        double precision targ(3,*),scj(*)

        integer nlevels,ntj
c
        double complex texps(nd,0:ntj,-ntj:ntj,*)
        
c
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          do j=jstart,jend
            call h3dformtac(nd,zk,scj(j),
     1        source(1,istart),charge(1,istart),ns,
     2        targ(1,j),ntj,texps(1,0,-ntj,j),wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
          do j=jstart,jend
            call h3dformtad(nd,zk,scj(j),
     1        source(1,istart),
     2        dipvec(1,1,istart),ns,targ(1,j),ntj,texps(1,0,-ntj,j),
     3        wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.1.and.ifdipole.eq.1) then
          do j=jstart,jend
            call h3dformtacd(nd,zk,scj(j),
     1        source(1,istart),charge(1,istart),
     2        dipvec(1,1,istart),ns,targ(1,j),ntj,texps(1,0,-ntj,j),
     3        wlege,nlege)
           enddo
         endif

c
        return
        end
c------------------------------------------------------------------     
