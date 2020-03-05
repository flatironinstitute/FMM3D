c
       subroutine lfmm3d(nd,eps,nsource,source,ifcharge,
     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
c
c        Laplace FMM in R^{3}: evaluate all pairwise particle
c        interactions (ignoring self-interactions) and interactions
c        with targs.
c
c        We use (1/r) for the Green's function, without the
c        1/(4 \pi) scaling.
c
c
c        Input parameters:
c
c   nd:   number of densities
c
c   eps:  requested precision
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
c   charge    in: double precision (nsource) 
c              charge strengths
c
c   ifdipole   in: integer
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c
c
c   dipvec   in: double precision (3,nsource) 
c              dipole orientation vectors
c
c   ifpgh   in: integer
c              flag for evaluating potential/gradient at the sources
c              ifpgh = 1, only potential is evaluated
c              ifpgh = 2, potential and gradients are evaluated
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
c
c   pot:    out: double precision(nd,nsource) 
c               potential at the source locations
c
c   grad:   out: double precision(nd,3,nsource)
c               gradient at the source locations
c
c   hess    out: double precision(nd,6,nsource)
c               hessian at the source locations
c                 (currently not supported)
c
c   pottarg:    out: double precision(nd,ntarg) 
c               potential at the targ locations
c
c   gradtarg:   out: double precision(nd,3,ntarg)
c               gradient at the targ locations
c
c   hesstarg    out: double precision(nd,6,ntarg)
c                hessian at the target locations - currently not
c                supported
     
       implicit none

       integer nd

       double precision eps

       integer ifcharge,ifdipole
       integer ifpgh,ifpghtarg

       integer ntarg,nsource
       

       double precision source(3,*),targ(3,*)
       double precision charge(nd,*)
       double precision dipvec(nd,3,*)

       double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
       double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

c
cc       tree variables
c
       integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
       integer *8 ltree
       integer nlmax
       integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
       integer *8 ipointer(32)
       integer, allocatable :: itree(:)
       double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc       temporary sorted arrays
c
       double precision, allocatable :: sourcesort(:,:),targsort(:,:)
       double precision, allocatable :: radsrc(:)
       double precision, allocatable :: chargesort(:,:)
       double precision, allocatable :: dipvecsort(:,:,:)

       double precision, allocatable :: potsort(:,:),gradsort(:,:,:)
       double precision, allocatable :: hesssort(:,:,:)
       double precision, allocatable :: pottargsort(:,:)
       double precision, allocatable :: gradtargsort(:,:,:)
       double precision, allocatable :: hesstargsort(:,:,:)
c
cc        temporary fmm arrays
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
cc       temporary variables not main fmm routine but
c        not used in particle code
       double precision expc(3),scjsort(1),radexp
       double complex texpssort(100)
       double precision expcsort(3)
       integer ntj,nexpc,nadd,ifnear

c
cc         other temporary variables
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
cc        figure out tree structure
c
c
cc         set criterion for box subdivision
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
c       turn on computation of list 1
c
      ifnear = 1




c
cc      set tree flags
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
cc     memory management code for contructing level restricted tree
        iert = 0
        call mklraptreemem(iert,source,nsource,radsrc,targ,ntarg,
     1        expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,
     2        nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3,
     3        mnlist4,mhung,ltree)

        if(ifprint.ge.1) print *, ltree/1.0d9
        if(ifprint.ge.1) print *, "mnlist3 = ",mnlist3
        if(ifprint.ge.1) print *, "mnlist4 = ",mnlist4



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

c     Allocate sorted source and targ arrays      

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
        scales(ilev) = boxsize(ilev)
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
         call l3dterms(eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2*nd
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

c
cc     reorder sources 
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) call dreorderf(nd,nsource,charge,chargesort,
     1                     itree(ipointer(5)))

      if(ifdipole.eq.1) then
         call dreorderf(3*nd,nsource,dipvec,dipvecsort,
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


      allocate(rmlexp(lmptot),stat=ier)
      if(ier.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         stop
      endif

c     Memory allocation is complete. 
c     Call main fmm routine

      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call lfmm3dmain(nd,eps,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,
     $   epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $   scales,treecenters,itree(ipointer(1)),nterms,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,hesstargsort,ntj,
     $   texpssort,scjsort,ifnear)

      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)



      if(ifpgh.eq.1) then
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
      endif
      if(ifpgh.eq.2) then 
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(3*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
      endif

      if(ifpgh.eq.3) then 
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(3*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
        call dreorderi(6*nd,nsource,hesssort,hess,
     1                 itree(ipointer(5)))
      endif


      if(ifpghtarg.eq.1) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
        call dreorderi(6*nd,ntarg,hesstargsort,hesstarg,
     1     itree(ipointer(6)))
      endif

      return
      end

c       
c---------------------------------------------------------------
c
      subroutine lfmm3dmain(nd,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipvecsort,
     $     ntarg,targsort,nexpc,expcsort,
     $     epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,isep,ndiv,nlevels, 
     $     nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $     rscales,centers,laddr,nterms,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,ntj,
     $     tsort,scjsort,ifnear)
      implicit none

      integer nd
      double precision eps
      integer nsource,ntarg,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      double precision epsfmm

      double precision sourcesort(3,nsource)

      double precision chargesort(nd,*)
      double precision dipvecsort(nd,3,*)

      double precision targsort(3,ntarg)

      double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

      integer ntj
      integer ifnear
      double precision expcsort(3,nexpc)
      double complex tsort(nd,0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

      integer *8 iaddr(2,nboxes), lmptot
      integer lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)

      double precision thresh
       
      double precision timeinfo(10)
      double precision centers(3,nboxes)

      integer isep
      integer *8 ltree
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(32)
      integer itree(ltree)
      integer nboxes
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll,m,idim,igbox
      integer ibox,jbox,ilev,npts,npts0,kbox,dir
      integer nchild,nlist1,nlist2,nlist3,nlist4

      integer istart,iend,istarts,iends
      integer istartt,iendt,istarte,iende
      integer isstart,isend,jsstart,jsend
      integer jstart,jend

      integer ifprint

      double precision d,time1,time2,second,omp_get_wtime
      double precision pottmp,fldtmp(3),hesstmp(3)

c     PW variables
      integer nexpmax, nlams, nmax, nthmax, nphmax,nmax2
      integer lca
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
  
      double precision, allocatable :: rlams(:),whts(:)

      double precision, allocatable :: rlsc(:,:,:)
      integer, allocatable :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:)
      double complex, allocatable :: yshift(:,:)
      double precision, allocatable :: zshift(:,:)

      double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
      double complex, allocatable :: mexp(:,:,:,:)
      double complex, allocatable :: mexpf1(:,:),mexpf2(:,:)
      double complex, allocatable ::
     1       mexpp1(:,:),mexpp2(:,:),mexppall(:,:,:)

      double complex, allocatable :: tmp(:,:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)

      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      double precision, allocatable :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

c     list 3 variables
      double complex, allocatable :: iboxlexp(:,:)
      double precision iboxsubcenters(3,8)
      double precision, allocatable :: iboxpot(:,:)
      double precision, allocatable :: iboxgrad(:,:,:)
      double precision, allocatable :: iboxsrc(:,:)
      integer, allocatable :: iboxsrcind(:)
      integer iboxfl(2,8)
c     end of list 3 variables
c     list 4 variables
      integer cntlist4
      integer, allocatable :: list4(:),ilist4(:)
      double complex, allocatable :: gboxmexp(:,:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision  gboxsubcenters(3,8)
      double precision, allocatable :: gboxsort(:,:)
      integer, allocatable :: gboxind(:)
      integer gboxfl(2,8)
      double precision, allocatable :: gboxcgsort(:,:)
      double precision, allocatable :: gboxdpsort(:,:,:)
c     end of list 4 variables

      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/

      pi = 4.0d0*atan(1.0d0)

      thresh = 2.0d0**(-52)*boxsize(0)

c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, 
c     and other things if ifprint=2.
c       
      ifprint=0
      

c     Initialize routines for plane wave mp loc translation
 
      if(isep.eq.1) then
         if(eps.ge.0.5d-2) nlams = 7
         if(eps.lt.0.5d-2.and.eps.ge.0.5d-3) nlams = 12
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 20
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 29
         if(eps.lt.0.5d-9) nlams = 37
      endif
      if(isep.eq.2) then
         if(eps.ge.0.5d-3) nlams = 9
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 15
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 22
         if(eps.lt.0.5d-9) nlams = 29
      endif

      allocate(rlams(nlams),whts(nlams))
      allocate(nphysical(nlams),nfourier(nlams))

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo
      allocate(rscpow(0:nmax))
      allocate(carray(4*nmax+1,4*nmax+1))
      allocate(dc(0:4*nmax,0:4*nmax))
      allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rlsc(0:nmax,0:nmax,nlams))


c     generate rotation matrices and carray
      call getpwrotmat(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


c     generate rlams and weights (these are the nodes
c     and weights for the lambda integral)

      if(isep.eq.1) call vwts(rlams,whts,nlams)
      if(isep.eq.2) call lwtsexp3sep2(nlams,rlams,whts,errtmp)


c     generate the number of fourier modes required to represent the
c     moment function in fourier space

      if(isep.eq.1) call numthetahalf(nfourier,nlams)
      if(isep.eq.2) call numthetahalf_isep2(nfourier,nlams)
 
c     generate the number of fourier modes in physical space
c     required for the exponential representation
      if(isep.eq.1) call numthetafour(nphysical,nlams)
      if(isep.eq.2) call numthetasix(nphysical,nlams)

c     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)

c     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      nn = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
         nn = nn + nphysical(i)*nfourier(i)
      enddo

      allocate(fexpe(nn),fexpo(nn),fexpback(nn))
      allocate(tmp(nd,0:nmax,-nmax:nmax))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nd,nexptot),mexpf2(nd,nexptot),
     1   mexpp1(nd,nexptotp))
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

      allocate(list4(nboxes))
      allocate(ilist4(nboxes))
      do i=1,nboxes
        list4(i)=0
        ilist4(i)=0
      enddo
cccccc      allocate(gboxwexp(nd,nexptotp,8,6))
      cntlist4=0

c     Precompute table for shifting exponential coefficients in 
c     physical domain
      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

c     Precompute table of exponentials for mapping from
c     fourier to physical domain
      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      
c
cc    compute array of factorials

     
      nmax2 = 2*nmax
      allocate(fact(0:nmax2),cs(0:nmax,-nmax:nmax))
      
      d = 1
      fact(0) = d
      do i=1,nmax2
        d=d*sqrt(i+0.0d0)
        fact(i) = d
      enddo

      cs(0,0) = 1.0d0
      do l=1,nmax
        do m=0,l
          cs(l,m) = ((-1)**l)/(fact(l-m)*fact(l+m))
          cs(l,-m) = cs(l,m)
        enddo
      enddo


      
      if(ifprint.ge.1) 
     1   call prin2('end of generating plane wave info*',i,0)
c
c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
      do i=1,nexpc
        do k=-ntj,ntj
          do j = 0,ntj
            do idim=1,nd
              tsort(idim,j,k,i)=0
            enddo
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c       
      do i=1,10
        timeinfo(i)=0
      enddo

c
c       ... set all multipole and local expansions to zero
c

      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
          call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
        enddo
C$OMP END PARALLEL DO        
      enddo

c
c      set scjsort
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
c     count number of boxes are in list4
      lca = 4*nmax
      if(ifprint.ge.1)
     $   call prinf('=== STEP 0 list4===*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev=1,nlevels-1
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nlist3=itree(ipointer(24)+ibox-1)
            if(nlist3.gt.0) then
              cntlist4=cntlist4+1
              list4(ibox)=cntlist4
              ilist4(cntlist4)=ibox
            endif
         enddo
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
      allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
      allocate(gboxmexp(nd*(nterms(ilev)+1)*
     1                   (2*nterms(ilev)+1),8,cntlist4))
cccccc  bad code, note gboxmexp is an array not scalar
      pgboxwexp=0d0
      gboxmexp=0d0
c     form mexp for all list4 type box at first ghost box center
      do ilev=1,nlevels-1

         rscpow(0) = 1.0d0/boxsize(ilev+1)
         rtmp = rscales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i)
C$OMP$PRIVATE(gboxind,gboxsort,gboxfl,gboxsubcenters)
C$OMP$PRIVATE(gboxcgsort,gboxdpsort,gboxwexp)
C$OMP$PRIVATE(mexpf1,mexpf2,tmp,mptemp)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            if(list4(ibox).gt.0) then
              istart=itree(ipointer(10)+ibox-1)
              iend=itree(ipointer(11)+ibox-1)
              npts = iend-istart+1
              if(npts.gt.0) then
                allocate(gboxind(npts))
                allocate(gboxsort(3,npts))
                allocate(gboxwexp(nd,nexptotp,6,8))
                call subdividebox(sourcesort(1,istart),npts,
     1               centers(1,ibox),boxsize(ilev+1),
     2               gboxind,gboxfl,gboxsubcenters)
                call dreorderf(3,npts,sourcesort(1,istart),
     1               gboxsort,gboxind)
                if(ifcharge.eq.1) then
                  allocate(gboxcgsort(nd,npts))
                  call dreorderf(nd,npts,chargesort(1,istart),
     1                 gboxcgsort,gboxind)
                endif
                if(ifdipole.eq.1) then
                  allocate(gboxdpsort(nd,3,npts))
                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart),
     1                 gboxdpsort,gboxind)
                endif
                do i=1,8
                  if(gboxfl(1,i).gt.0) then
                    jstart=gboxfl(1,i)
                    jend=gboxfl(2,i)
                    npts0=jend-jstart+1
                    jbox=list4(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1),
     1                     gboxsort(1,jstart),
     2                     gboxcgsort(1,jstart),
     3                     npts0,gboxsubcenters(1,i),nterms(ilev+1),
     4                     gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1),
     1                     gboxsort(1,jstart),
     2                     gboxdpsort(1,1,jstart),
     3                     npts0,gboxsubcenters(1,i),nterms(ilev+1),
     4                     gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1),
     1                     gboxsort(1,jstart),
     2                     gboxcgsort(1,jstart),
     3                     gboxdpsort(1,1,jstart),
     4                     npts0,gboxsubcenters(1,i),nterms(ilev+1),
     5                     gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1),
     1                   gboxsubcenters(1,i),gboxmexp(1,i,jbox),
     2                   nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3                   rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),
     1                   rscpow,tmp)
c
cc                process up down for current box
c
                    call mpoletoexp(nd,tmp,nterms(ilev+1),nlams,
     1                   nfourier,nexptot,mexpf1,mexpf2,rlsc)

                    call ftophys(nd,mexpf1,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,1,i),fexpe,fexpo)

                    call ftophys(nd,mexpf2,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,2,i),fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i),
     1                   gboxwexp(1,1,2,i),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2),
     3                   xshift,yshift,zshift)
c
cc                process north-south for current box
c
                    call rotztoy(nd,nterms(ilev+1),tmp,mptemp,rdminus)
                    call mpoletoexp(nd,mptemp,nterms(ilev+1),nlams,
     1                   nfourier,nexptot,mexpf1,mexpf2,rlsc)

                    call ftophys(nd,mexpf1,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,3,i),fexpe,fexpo)

                    call ftophys(nd,mexpf2,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,4,i),fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i),
     1                   gboxwexp(1,1,4,i),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4),
     3                   xshift,yshift,zshift)
c
cc                process east-west for current box
c
                    call rotztox(nd,nterms(ilev+1),tmp,mptemp,rdplus)
                    call mpoletoexp(nd,mptemp,nterms(ilev+1),nlams,
     1                   nfourier,nexptot,mexpf1,mexpf2,rlsc)

                    call ftophys(nd,mexpf1,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,5,i),fexpe,fexpo)

                    call ftophys(nd,mexpf2,nlams,rlams,nfourier,
     1                   nphysical,nthmax,gboxwexp(1,1,6,i),fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i),
     1                   gboxwexp(1,1,6,i),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6),
     3                   xshift,yshift,zshift)
                  endif
                enddo
                deallocate(gboxind,gboxsort)
                if(ifcharge.eq.1) then
                  deallocate(gboxcgsort)
                endif
                if(ifdipole.eq.1) then
                  deallocate(gboxdpsort)
                endif
                deallocate(gboxwexp)
              endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
      timeinfo(2)=time2-time1
c     end of count number of boxes are in list4
c

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

               if(npts.gt.0.and.nchild.eq.0.and.list4(ibox).eq.0) then
                  call l3dformmpc(nd,rscales(ilev),
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

               if(npts.gt.0.and.nchild.eq.0.and.list4(ibox).eq.0) then
                  call l3dformmpd(nd,rscales(ilev),
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

               if(npts.gt.0.and.nchild.eq.0.and.list4(ibox).eq.0) then
                  call l3dformmpcd(nd,rscales(ilev),
     1            sourcesort(1,istart),chargesort(1,istart),
     2            dipvecsort(1,1,istart),npts,
     2            centers(1,ibox),nterms(ilev),
     3            rmlexp(iaddr(1,ibox)),wlege,nlege)          
               endif
            enddo
C$OMP END PARALLEL DO          
         endif
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"leaf:",cntlist4



      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1



      lca = 4*nmax


c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,0,-1
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
                     call l3dmpmp(nd,rscales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                  endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$        time1=omp_get_wtime()

c
cc     zero out mexp
c 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,idim)
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


      do ilev=2,nlevels

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,mptemp)
         do ibox=laddr(1,ilev),laddr(2,ilev)

            istart = itree(ipointer(10)+ibox-1)
            iend = itree(ipointer(11)+ibox-1)

            npts = iend-istart+1

            if(npts.gt.0) then
c            rescale the multipole expansion

                call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1                 rscpow,tmp)
c
cc                process up down for current box
c
                call mpoletoexp(nd,tmp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,1),fexpe,fexpo)

                call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,2),fexpe,fexpo)


c
cc                process north-south for current box
c
                call rotztoy(nd,nterms(ilev),tmp,mptemp,rdminus)
                call mpoletoexp(nd,mptemp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,3),fexpe,fexpo)

                call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,4),fexpe,fexpo)

c
cc                process east-west for current box

                call rotztox(nd,nterms(ilev),tmp,mptemp,rdplus)
                call mpoletoexp(nd,mptemp,nterms(ilev),nlams,nfourier,
     1              nexptot,mexpf1,mexpf2,rlsc)

                call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,5),fexpe,fexpo)


                call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical,
     1          nthmax,mexp(1,1,ibox,6),fexpe,fexpo)

            endif

         enddo
C$OMP END PARALLEL DO         
c
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
         rscpow(0) = 1.0d0
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
C$OMP$PRIVATE(mexpf1,mexpf2,mexpp1,mexpp2,mexppall)
C$OMP$PRIVATE(nuall,uall,ndall,dall,nnall,nall,nsall,sall)
C$OMP$PRIVATE(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)
C$OMP$PRIVATE(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)
C$OMP$PRIVATE(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)
C$OMP$PRIVATE(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)
C$OMP$PRIVATE(nw2,w2,nw4,w4,nw6,w6,nw8,w8)
C$OMP$PRIVATE(npts0,nlist3,ctmp,jstart,jend,i,iboxfl,iboxsubcenters)
C$OMP$PRIVATE(iboxpot,iboxgrad,iboxlexp,iboxsrc,iboxsrcind)
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
     1         itree(ipointer(18)+ibox-1),itree(ipointer(19)+
     2         mnbors*(ibox-1)),nchild,itree(ipointer(4)),centers,
     3         isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,neall,
     4         eall,nwall,wall,nu1234,u1234,nd5678,d5678,nn1256,n1256,
     5         ns3478,s3478,ne1357,e1357,nw2468,w2468,nn12,n12,nn56,n56,
     6         ns34,s34,ns78,s78,ne13,e13,ne57,e57,nw24,w24,nw68,w68,
     7         ne1,e1,ne3,e3,ne5,e5,ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)

               call processudexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678,
     5         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     6         mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),xshift,
     7         yshift,zshift,fexpback,rlsc,rscpow,
     8         pgboxwexp,cntlist4,list4,
     9         itree(ipointer(26)),itree(ipointer(27)),mnlist4)
               
               call processnsexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,
     5         ns3478,s3478,ns34,s34,ns78,s78,
     6         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7         mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8         mexppall(1,1,5),mexppall(1,1,6),mexppall(1,1,7),
     9         mexppall(1,1,8),rdplus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4,
     9         itree(ipointer(26)),itree(ipointer(27)),mnlist4)

               
               call processewexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(4)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,
     5         ne3,e3,ne5,e5,ne7,e7,nwall,wall,
     5         nw2468,w2468,nw24,w24,nw68,w68,
     5         nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6         mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7         mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8         mexppall(1,1,5),mexppall(1,1,6),
     8         mexppall(1,1,7),mexppall(1,1,8),mexppall(1,1,9),
     9         mexppall(1,1,10),mexppall(1,1,11),mexppall(1,1,12),
     9         mexppall(1,1,13),mexppall(1,1,14),mexppall(1,1,15),
     9         mexppall(1,1,16),rdminus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4,
     9         itree(ipointer(26)),itree(ipointer(27)),mnlist4)


            endif

            nlist3 = itree(ipointer(24)+ibox-1)
            if(nlist3.gt.0.and.npts.gt.0) then
              call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1             nlist3,itree(ipointer(25)+(ibox-1)*mnlist3),isep,
     2             centers,nuall,uall,ndall,dall,nnall,nall,
     3                     nsall,sall,neall,eall,nwall,wall)
              allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1                 (2*nterms(ilev)+1),8))
              iboxlexp=0
              call processlist3udexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2             nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,nuall,uall,ndall,dall,
     4             mexpf1,mexpf2,mexpp1,mexpp2,
     5             mexppall(1,1,1),mexppall(1,1,2),
     6             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3nsexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2             nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,nnall,nall,nsall,sall,
     4             mexpf1,mexpf2,mexpp1,mexpp2,
     5             mexppall(1,1,1),mexppall(1,1,2),rdplus,
     6             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3ewexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2             nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,neall,eall,nwall,wall,
     4             mexpf1,mexpf2,mexpp1,mexpp2,
     5             mexppall(1,1,1),mexppall(1,1,2),rdminus,
     6             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              if(ifpgh.eq.1) then
                istart = itree(ipointer(10)+ibox-1)
                iend = itree(ipointer(11)+ibox-1)
                npts = iend-istart+1
                if(npts.gt.0) then
                  allocate(iboxsrcind(npts))
                  allocate(iboxsrc(3,npts))
                  allocate(iboxpot(nd,npts))
                  call subdividebox(sourcesort(1,istart),npts,
     1                 centers(1,ibox),boxsize(ilev),
     2                 iboxsrcind,iboxfl,iboxsubcenters)
                  call dreorderf(3,npts,sourcesort(1,istart),
     1                 iboxsrc,iboxsrcind)
                  call dreorderf(nd,npts,pot(1,istart),
     1                 iboxpot,iboxsrcind)
                  do i=1,8
                    if(iboxfl(1,i).gt.0) then
                      jstart=iboxfl(1,i)
                      jend=iboxfl(2,i)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev),
     1                       iboxsubcenters(1,i),iboxlexp(1,i),
     2                       nterms(ilev),iboxsrc(1,jstart),npts0,
     3                       iboxpot(1,jstart),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot,pot(1,istart),
     1                 iboxsrcind)
                  deallocate(iboxsrcind,iboxsrc)
                  deallocate(iboxpot)
                endif
              endif

cccccc        todo
              if(ifpgh.eq.2) then
                istart = itree(ipointer(10)+ibox-1)
                iend = itree(ipointer(11)+ibox-1)
                npts = iend-istart+1
                if(npts.gt.0) then
                  allocate(iboxsrcind(npts))
                  allocate(iboxsrc(3,npts))
                  allocate(iboxpot(nd,npts))
                  allocate(iboxgrad(nd,3,npts))
                  call subdividebox(sourcesort(1,istart),npts,
     1                 centers(1,ibox),boxsize(ilev),
     2                 iboxsrcind,iboxfl,iboxsubcenters)
                  call dreorderf(3,npts,sourcesort(1,istart),
     1                 iboxsrc,iboxsrcind)
                  call dreorderf(nd,npts,pot(1,istart),
     1                 iboxpot,iboxsrcind)
                  call dreorderf(3*nd,npts,grad(1,1,istart),
     1                 iboxgrad,iboxsrcind)
                  do i=1,8
                    if(iboxfl(1,i).gt.0) then
                      jstart=iboxfl(1,i)
                      jend=iboxfl(2,i)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev),
     1                       iboxsubcenters(1,i),iboxlexp(1,i),
     2                       nterms(ilev),iboxsrc(1,jstart),npts0,
     3                       iboxpot(1,jstart),iboxgrad(1,1,jstart),
     4                       wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot,pot(1,istart),
     1                 iboxsrcind)
                  call dreorderi(3*nd,npts,iboxgrad,grad(1,1,istart),
     1                 iboxsrcind)
                  deallocate(iboxsrcind,iboxsrc)
                  deallocate(iboxpot,iboxgrad)
                endif
              endif

              if(ifpghtarg.eq.1) then
                istart = itree(ipointer(12)+ibox-1)
                iend = itree(ipointer(13)+ibox-1)
                npts = iend-istart+1
                if(npts.gt.0) then
                  allocate(iboxsrcind(npts))
                  allocate(iboxsrc(3,npts))
                  allocate(iboxpot(nd,npts))
                  call subdividebox(targsort(1,istart),npts,
     1                 centers(1,ibox),boxsize(ilev),
     2                 iboxsrcind,iboxfl,iboxsubcenters)
                  call dreorderf(3,npts,targsort(1,istart),
     1                 iboxsrc,iboxsrcind)
                  call dreorderf(nd,npts,pottarg(1,istart),
     1                 iboxpot,iboxsrcind)
                  do i=1,8
                    if(iboxfl(1,i).gt.0) then
                      jstart=iboxfl(1,i)
                      jend=iboxfl(2,i)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev),
     1                       iboxsubcenters(1,i),iboxlexp(1,i),
     2                       nterms(ilev),iboxsrc(1,jstart),npts0,
     3                       iboxpot(1,jstart),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot,pottarg(1,istart),
     1                 iboxsrcind)
                  deallocate(iboxsrcind,iboxsrc)
                  deallocate(iboxpot)
                endif
              endif

              if(ifpghtarg.eq.2) then
                istart = itree(ipointer(12)+ibox-1)
                iend = itree(ipointer(13)+ibox-1)
                npts = iend-istart+1
                if(npts.gt.0) then
                  allocate(iboxsrcind(npts))
                  allocate(iboxsrc(3,npts))
                  allocate(iboxpot(nd,npts))
                  allocate(iboxgrad(nd,3,npts))
                  call subdividebox(targsort(1,istart),npts,
     1                 centers(1,ibox),boxsize(ilev),
     2                 iboxsrcind,iboxfl,iboxsubcenters)
                  call dreorderf(3,npts,targsort(1,istart),
     1                 iboxsrc,iboxsrcind)
                  call dreorderf(nd,npts,pottarg(1,istart),
     1                 iboxpot,iboxsrcind)
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart),
     1                 iboxgrad,iboxsrcind)
                  do i=1,8
                    if(iboxfl(1,i).gt.0) then
                      jstart=iboxfl(1,i)
                      jend=iboxfl(2,i)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev),
     1                       iboxsubcenters(1,i),iboxlexp(1,i),
     2                       nterms(ilev),iboxsrc(1,jstart),npts0,
     3                       iboxpot(1,jstart),iboxgrad(1,1,jstart),
     4                       wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot,pottarg(1,istart),
     1                 iboxsrcind)
                  call dreorderi(3*nd,npts,iboxgrad,
     1                 gradtarg(1,1,istart),iboxsrcind)
                  deallocate(iboxsrcind,iboxsrc)
                  deallocate(iboxpot,iboxgrad)
                endif
              endif
ccccccc       end of todo
              deallocate(iboxlexp)
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(3) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

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
                     call l3dlocloc(nd,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),rscales(ilev+1),centers(1,jbox),
     3                rmlexp(iaddr(2,jbox)),nterms(ilev+1),dc,lca)
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval lo) ===*',i,0)

c     ... step 6, evaluate all local expansions
c

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

                  call l3dlocloc(nd,rscales(ilev),
     1             centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2             nterms(ilev),rscales(ilev),expcsort(1,i),
     3             tsort(1,0,-ntj,i),ntj,dc,lca)
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
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(10)+ibox-1)
              iend = itree(ipointer(11)+ibox-1)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(10)+ibox-1)
              iend = itree(ipointer(11)+ibox-1)
              npts = iend-istart+1
              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(12)+ibox-1)
              iend = itree(ipointer(13)+ibox-1)
              npts = iend-istart+1
              call l3dtaevalp(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(3)+ibox-1)
            if(nchild.eq.0) then 
              istart = itree(ipointer(12)+ibox-1)
              iend = itree(ipointer(13)+ibox-1)
              npts = iend-istart+1

              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif
      enddo

    
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(5) = time2 - time1


      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()

      if(ifnear.eq.0) goto 1000
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

               call lfmm3dexpc_direct(nd,jstart,jend,istarte,
     1         iende,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipvecsort,expcsort,tsort,scjsort,ntj,
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcp(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectdp(nd,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcdp(nd,sourcesort(1,jstart),
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcg(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),thresh)   
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectdg(nd,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcdg(nd,sourcesort(1,jstart),
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcp(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectdp(nd,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,nlist1,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcdp(nd,sourcesort(1,jstart),
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcg(nd,sourcesort(1,jstart),
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectdg(nd,sourcesort(1,jstart),
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
C$OMP$SCHEDULE(DYNAMIC)      
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
                call l3ddirectcdg(nd,sourcesort(1,jstart),
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
 1000 continue      
 
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
c------------------------------------------------
      subroutine lfmm3dexpc_direct(nd,istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,
     $     dipvec,expc,texps,scj,ntj,wlege,nlege)
c--------------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the expansion center array to the existing 
c     local expansions at the corresponding expansion centers.
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
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
c                  First index in the expansion center array at 
c                  which we  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in expansion center array at 
c                  which we wish to compute the expansions
c 
c     scjsort      in: double precision(*)
c                  Scale of expansions formed at the expansion centers
c
c     source       in: double precision(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: double precision
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipvec      in: double precision(3,ns)
c                 Dipole orientation vector at the source locations
c
c     expc        in: double precision(3,nexpc)
c                 Expansion center locations
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
c   Updated expansions at the targs
c   texps       out: double complex(0:ntj,-ntj:ntj,expc) 
c                 coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j, nlege
        integer ifcharge,ifdipole,ier,nd
        double precision source(3,*)
        double precision scj(*)
        double precision wlege(*)
        double precision charge(nd,*)
        double precision dipvec(nd,3,*)
        double precision expc(3,*)

        integer nlevels,ntj
c
        double complex texps(nd,0:ntj,-ntj:ntj,*)
        
c
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          do j=jstart,jend
            call l3dformtac(nd,scj(j),
     1        source(1,istart),charge(1,istart),ns,
     2        expc(1,j),ntj,texps(1,0,-ntj,j),wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.0.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtad(nd,scj(j),
     1        source(1,istart),
     2        dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j),
     3        wlege,nlege)
           enddo
         endif

         if(ifcharge.eq.1.and.ifdipole.eq.1) then
          do j=jstart,jend
            call l3dformtacd(nd,scj(j),
     1        source(1,istart),charge(1,istart),
     2        dipvec(1,1,istart),ns,expc(1,j),ntj,texps(1,0,-ntj,j),
     3        wlege,nlege)
           enddo
         endif

c
        return
        end
