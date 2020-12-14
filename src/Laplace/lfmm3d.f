c
       subroutine lfmm3d(nd,eps,nsource,source,ifcharge,
     $    charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,
     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
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
c   nsource in: integer(8)  
c                number of sources
c
c   source  in: double precision (3,nsource)
c                source(k,j) is the kth component of the jth
c                source locations
c
c   ifcharge  in: integer(8)  
c             charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c 
c   charge    in: double precision (nsource) 
c              charge strengths
c
c   ifdipole   in: integer(8)
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c
c
c   dipvec   in: double precision (3,nsource) 
c              dipole orientation vectors
c   iper    in: integer(8)
c             flag for periodic implmentations. Currently unused
c
c   ifpgh   in: integer(8)
c              flag for evaluating potential/gradient at the sources
c              ifpgh = 1, only potential is evaluated
c              ifpgh = 2, potential and gradients are evaluated
c              ifpgh = 3, potential, gradients, and hessian are
c                 evaluated
c
c   ntarg  in: integer(8)  
c                 number of targs 
c
c   targ  in: double precision (3,ntarg)
c               targ(k,j) is the kth component of the jth
c               targ location
c
c   ifpghtarg   in: integer(8)
c              flag for evaluating potential/gradient at the targs
c              ifpghtarg = 1, only potential is evaluated
c              ifpghtarg = 2, potential and gradient are evaluated
c              ifpghtarg = 3, potential, gradient, and hess are
c                  evaluated
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
c
c   pottarg:    out: double precision(nd,ntarg) 
c               potential at the targ locations
c
c   gradtarg:   out: double precision(nd,3,ntarg)
c               gradient at the targ locations
c
c   hesstarg    out: double precision(nd,6,ntarg)
c                hessian at the target locations 
c   ier         out: integer(8)
c                error flag
c                ier = 0, for successful execution
c                ier = 4, if failed to allocate workspace
c                      for multipole and local expansions
c                ier = 8, if failed to allocate workspace
c                      for plane wave expansions
c     
c     
       implicit none

       integer(8) nd,ier,iper,ndim

       double precision eps

       integer(8) ifcharge,ifdipole
       integer(8) ifpgh,ifpghtarg

       integer(8) ntarg,nsource
       

       double precision source(3,*),targ(3,*)
       double precision charge(nd,*)
       double precision dipvec(nd,3,*)

       double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
       double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

c
cc       tree variables
c
       integer(8) idivflag,ndiv,nboxes,nlevels
       integer(8) nlmax,nlmin,ifunif
       integer(8) ipointer(8),ltree
       integer(8), allocatable :: itree(:)
       integer(8), allocatable :: isrcse(:,:),itargse(:,:),isrc(:)
       integer(8), allocatable :: itarg(:)
       integer(8), allocatable :: iexpcse(:,:)
       integer(8) iexpc
       double precision, allocatable :: treecenters(:,:),boxsize(:)
       double precision b0,b0inv,b0inv2,b0inv3

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
       integer(8), allocatable :: nterms(:)
       integer(8), allocatable :: iaddr(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer(8) lmptemp,nmax
       integer(8) lmptot
       double precision, allocatable :: mptemp(:),mptemp2(:)

c
cc       temporary variables not main fmm routine but
c        not used in particle code
       double precision expc(3),scjsort(1),radexp
       double complex texpssort(100)
       double precision expcsort(3)
       integer(8) ntj,nexpc,nadd,ifnear

c
cc         other temporary variables
c
        integer(8) i,iert,ifprint,ilev,idim
        double precision time1,time2,omp_get_wtime,second

c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c      

      call cpu_time(time1)
C$     time1=omp_get_wtime()      
      ifprint=0

c
cc        figure out tree structure
c
c
cc         set criterion for box subdivision
c
      call lndiv(eps,nsource,ntarg,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag) 


c
c       turn on computation of list 1
c
      ifnear = 1


       nexpc = 0
       nadd = 0
       ntj = 0
       ndim = 3


c
cc      set tree flags
c 
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0
       iper = 0
       ifunif = 0

c
cc     memory management code for contructing level restricted tree
      call pts_tree_mem(source,nsource,targ,ntarg,idivflag,ndiv,nlmin,
     1  nlmax,iper,ifunif,nlevels,nboxes,ltree)
      

        if(ifprint.ge.1) print *, ltree/1.0d9

        allocate(itree(ltree))
        allocate(boxsize(0:nlevels))
        allocate(treecenters(3,nboxes))

c       Call tree code
      call pts_tree_build(source,nsource,targ,ntarg,idivflag,ndiv,
     1  nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,
     2  treecenters,boxsize)
      

      allocate(isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes))
      allocate(isrc(nsource),itarg(ntarg))

      call pts_tree_sort(nsource,source,itree,ltree,nboxes,nlevels,
     1   ipointer,treecenters,isrc,isrcse)
      
      call pts_tree_sort(ntarg,targ,itree,ltree,nboxes,nlevels,
     1   ipointer,treecenters,itarg,itargse)
      
      call pts_tree_sort(nexpc,expc,itree,ltree,nboxes,nlevels,
     1   ipointer,treecenters,iexpc,iexpcse)

c
c   End of tree build
c

c
c  Set rescaling parameters
c
      b0 = boxsize(0)
      b0inv = 1.0d0/b0
      b0inv2 = b0inv**2
      b0inv3 = b0inv2*b0inv

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

      allocate(nterms(0:nlevels))

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
      call dreorderf(ndim,nsource,source,sourcesort,isrc)

c
c       rescale sources to be contained in unit box
c
      call drescale(3*nsource,sourcesort,b0inv)

      if(ifcharge.eq.1) then
        call dreorderf(nd,nsource,charge,chargesort,
     1                     isrc)
        call drescale(nd*nsource,chargesort,b0inv)
      endif


      if(ifdipole.eq.1) then
         call dreorderf(3*nd,nsource,dipvec,dipvecsort,
     1       isrc)
         call drescale(3*nd*nsource,dipvecsort,b0inv2)
      endif

c
cc      reorder and rescale targs
c
      call dreorderf(ndim,ntarg,targ,targsort,itarg)
      call drescale(3*ntarg,targsort,b0inv)


c
c        update tree centers and boxsize
c
      call drescale(3*nboxes,treecenters,b0inv)
      call drescale(nlevels+1,boxsize,b0inv)

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
         ier = 4
         return
      endif

c     Memory allocation is complete. 
c     scaling factor for multipole and local expansions at all levels
c
      allocate(scales(0:nlevels))
      do ilev = 0,nlevels
        scales(ilev) = boxsize(ilev)
      enddo

      call cpu_time(time2)
C$     time2=omp_get_wtime()      

      if(ifprint.ge.1) 
     1   call prin2('time before fmm main=*',time2-time1,1)
c     Call main fmm routine

      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call lfmm3dmain(nd,eps,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,
     $   iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,ndiv,nlevels,
     $   nboxes,iper,boxsize,treecenters,isrcse,itargse,iexpcse,
     $   scales,itree(ipointer(1)),nterms,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,hesstargsort,ntj,
     $   texpssort,scjsort,ifnear,ier)
      if(ier.ne.0) return

      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)



      if(ifpgh.ge.1) then
        call dreorderi(nd,nsource,potsort,pot,isrc)
      endif
      if(ifpgh.ge.2) then 
        call dreorderi(3*nd,nsource,gradsort,grad,isrc)
        call drescale(nd*3*nsource,grad,b0inv)
      endif

      if(ifpgh.ge.3) then 
        call dreorderi(6*nd,nsource,hesssort,hess,isrc)
        call drescale(nd*6*nsource,hess,b0inv2)
      endif


      if(ifpghtarg.ge.1) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.ge.2) then
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,itarg)
        call drescale(nd*3*ntarg,gradtarg,b0inv)
      endif

      if(ifpghtarg.ge.3) then
        call dreorderi(6*nd,ntarg,hesstargsort,hesstarg,itarg)
        call drescale(nd*6*ntarg,hesstarg,b0inv2)
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
     $     iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,ndiv,nlevels, 
     $     nboxes,iper,boxsize,centers,isrcse,itargse,iexpcse,
     $     rscales,laddr,nterms,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,ntj,
     $     tsort,scjsort,ifnear,ier)
      implicit none

      integer(8) nd,ndim
      integer(8) ier
      double precision eps
      integer(8) nsource,ntarg,nexpc
      integer(8) ndiv,nlevels

      integer(8) ifcharge,ifdipole
      integer(8) ifpgh,ifpghtarg

      double precision sourcesort(3,nsource)

      double precision chargesort(nd,*)
      double precision dipvecsort(nd,3,*)

      double precision targsort(3,ntarg)

      double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
      double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

      integer(8) ntj
      integer(8) ifnear
      double precision expcsort(3,nexpc)
      double complex tsort(nd,0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

      integer(8) iaddr(2,nboxes), lmptot
      integer(8) lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)

      double precision thresh
       
      double precision timeinfo(10)
      double precision centers(3,nboxes)

      integer(8) isep,iper
      integer(8) laddr(2,0:nlevels)
      integer(8) nterms(0:nlevels)
      integer(8) ipointer(8),ltree
      integer(8) itree(ltree)
      integer(8) nboxes
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)
      integer(8) isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
      integer(8), allocatable :: nlist1(:),list1(:,:)
      integer(8), allocatable :: nlist2(:),list2(:,:)
      integer(8), allocatable :: nlist3(:),list3(:,:)
      integer(8), allocatable :: nlist4(:),list4(:,:)

      integer(8) nuall,ndall,nnall,nsall,neall,nwall
      integer(8) nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer(8) nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer(8) ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer(8), allocatable :: uall(:,:),dall(:,:),nall(:,:)
      integer(8), allocatable :: sall(:,:),eall(:,:),wall(:,:)
      integer(8), allocatable :: u1234(:,:),d5678(:,:)
      integer(8), allocatable :: n1256(:,:),s3478(:,:)
      integer(8), allocatable :: e1357(:,:),w2468(:,:)
      integer(8), allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      integer(8), allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      integer(8), allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      integer(8), allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

c     temp variables
      integer(8) i,j,k,l,ii,jj,kk,ll,m,idim,igbox
      integer(8) ibox,jbox,ilev,npts,npts0,kbox,dir
      integer(8) nchild

      integer(8) istart,iend,istarts,iends
      integer(8) istartt,iendt,istarte,iende
      integer(8) isstart,isend,jsstart,jsend
      integer(8) jstart,jend

      integer(8) ifprint

      double precision d,time1,time2,second,omp_get_wtime
      double precision pottmp,fldtmp(3),hesstmp(3)

c     PW variables
      integer(8) nexpmax, nlams, nmax, nthmax, nphmax,nmax2,nmaxt
      integer(8) lca
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
  
      double precision, allocatable :: rlams(:),whts(:)

      double precision, allocatable :: rlsc(:,:,:)
      integer(8), allocatable :: nfourier(:), nphysical(:)
      integer(8) nexptot, nexptotp
      double complex, allocatable :: xshift(:,:)
      double complex, allocatable :: yshift(:,:)
      double precision, allocatable :: zshift(:,:)

      double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
      double complex, allocatable :: mexp(:,:,:,:)
      double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, allocatable ::
     1       mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

      double complex, allocatable :: tmp(:,:,:,:)
      double precision, allocatable :: mptmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer(8) ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer(8) nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer(8) nterms_eval(4,0:nlevels)

      integer(8) mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer(8) ctr,nn,iptr1,iptr2
      double precision, allocatable :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

c     list 3 variables
      double complex, allocatable :: iboxlexp(:,:,:)
      double precision, allocatable :: iboxsubcenters(:,:,:)
      double precision, allocatable :: iboxpot(:,:,:)
      double precision, allocatable :: iboxgrad(:,:,:,:)
      double precision, allocatable :: iboxhess(:,:,:,:)
      double precision, allocatable :: iboxsrc(:,:,:)
      integer(8), allocatable :: iboxsrcind(:,:)
      integer(8), allocatable :: iboxfl(:,:,:)
c     end of list 3 variables
c     list 4 variables
      integer(8) cntlist4
      integer(8), allocatable :: list4ct(:),ilist4(:)
      double complex, allocatable :: gboxmexp(:,:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:,:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision, allocatable :: gboxsubcenters(:,:,:)
      double precision, allocatable :: gboxsort(:,:,:)
      integer(8), allocatable :: gboxind(:,:)
      integer(8), allocatable :: gboxfl(:,:,:)
      double precision, allocatable :: gboxcgsort(:,:,:)
      double precision, allocatable :: gboxdpsort(:,:,:,:)
c     end of list 4 variables

c
c   hessian variables
c
      double precision, allocatable :: scarray(:,:)

      integer(8) bigint
      integer(8) iert
      data ima/(0.0d0,1.0d0)/

      integer(8) nthd,ithd
      integer(8) omp_get_max_threads,omp_get_thread_num
      nthd = 1
C$    nthd=omp_get_max_threads()

      pi = 4.0d0*atan(1.0d0)

      thresh = 2.0d0**(-51)*boxsize(0)


c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, 
c     and other things if ifprint=2.
c       
      ifprint=0
      
c
c   initialize various tree lists
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      ndim = 3
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)
      

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

      call vwts(rlams,whts,nlams)


c     generate the number of fourier modes required to represent the
c     moment function in fourier space

      call numthetahalf(nfourier,nlams)
 
c     generate the number of fourier modes in physical space
c     required for the exponential representation
      call numthetafour(nphysical,nlams)

c     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)

c
c
c
      nn = 10*(nmax+2)**2
      allocate(scarray(nn,0:nlevels))

      do ilev=0,nlevels
        call l3dtaevalhessdini(nterms(ilev),scarray(1,ilev))
      enddo

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
      allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
      allocate(mptmp(lmptemp,nthd))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd),
     1   mexpp1(nd,nexptotp,nthd))
      allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))

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
        ier = 8
        return
      endif

      allocate(list4ct(nboxes))
      allocate(ilist4(nboxes))
      do i=1,nboxes
        list4ct(i)=0
        ilist4(i)=0
      enddo
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
            nchild = itree(ipointer(4)+ibox-1)
            if(nchild.gt.0) then
               istart = iexpcse(1,ibox)
               iend = iexpcse(2,ibox) 
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
            if(nlist3(ibox).gt.0) then
              cntlist4=cntlist4+1
              list4ct(ibox)=cntlist4
              ilist4(cntlist4)=ibox
            endif
         enddo
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
      allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
      allocate(gboxmexp(nd*(nterms(ilev)+1)*
     1                   (2*nterms(ilev)+1),8,cntlist4))



      allocate(gboxsubcenters(3,8,nthd))
      allocate(gboxfl(2,8,nthd))

      nmaxt = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
C$OMP$REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(list4ct(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
C$OMP END PARALLEL DO

      allocate(gboxind(nmaxt,nthd))
      allocate(gboxsort(3,nmaxt,nthd))
      allocate(gboxwexp(nd,nexptotp,6,8,nthd))
      allocate(gboxcgsort(nd,nmaxt,nthd))
      allocate(gboxdpsort(nd,3,nmaxt,nthd))

c   note gboxmexp is an array not scalar
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
C$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(sourcesort(1,istart),npts,
     1               centers(1,ibox),boxsize(ilev+1),
     2               gboxind(1,ithd),gboxfl(1,1,ithd),
     3               gboxsubcenters(1,1,ithd))
                call dreorderf(ndim,npts,sourcesort(1,istart),
     1               gboxsort(1,1,ithd),gboxind(1,ithd))
                if(ifcharge.eq.1) then
                  call dreorderf(nd,npts,chargesort(1,istart),
     1                 gboxcgsort(1,1,ithd),gboxind(1,ithd))
                endif
                if(ifdipole.eq.1) then
                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart),
     1                 gboxdpsort(1,1,1,ithd),gboxind(1,ithd))
                endif
                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1),
     1                   gboxsort(1,jstart,ithd),
     2                   gboxcgsort(1,jstart,ithd),
     3                   npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
     4                   gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1),
     1                   gboxsort(1,jstart,ithd),
     2                   gboxdpsort(1,1,jstart,ithd),
     3                   npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
     4                   gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1),
     1                   gboxsort(1,jstart,ithd),
     2                   gboxcgsort(1,jstart,ithd),
     3                   gboxdpsort(1,1,jstart,ithd),
     4                   npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
     5                   gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1),
     1                   gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox),
     2                   nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3                   rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),
     1                   rscpow,tmp(1,0,-nmax,ithd))
c
cc                process up down for current box
c
                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,1,i,ithd),
     3                   fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,2,i,ithd),
     3                   fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd),
     1                   gboxwexp(1,1,2,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2),
     3                   xshift,yshift,zshift)
c
cc                process north-south for current box
c
                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
     1                   mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,3,i,ithd),
     3                   fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,4,i,ithd),
     3                   fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd),
     1                   gboxwexp(1,1,4,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4),
     3                   xshift,yshift,zshift)
c
cc                process east-west for current box
c
                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
     1                   mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd),
     1                   nterms(ilev+1),nlams,
     2                   nfourier,nexptot,mexpf1(1,1,ithd),
     3                   mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,5,i,ithd),
     3                   fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
     1                   nlams,rlams,nfourier,
     2                   nphysical,nthmax,gboxwexp(1,1,6,i,ithd),
     3                   fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd),
     1                   gboxwexp(1,1,6,i,ithd),i,nexptotp,
     2                   pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6),
     3                   xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
      deallocate(gboxfl,gboxsubcenters,gboxwexp,gboxcgsort)
      deallocate(gboxdpsort,gboxind,gboxsort)

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

               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
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

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox) 
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
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

               istart = isrcse(1,ibox) 
               iend = isrcse(2,ibox)
               npts = iend-istart+1

               nchild = itree(ipointer(4)+ibox-1)

               if(npts.gt.0.and.nchild.eq.0.and.list4ct(ibox).eq.0) then
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
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = isrcse(1,jbox)
                  iend = isrcse(2,jbox)
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
      timeinfo(2)=timeinfo(2)+time2-time1

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

c     init uall,dall,...,etc arrays
      allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
      allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
      allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
      allocate(s3478(24,nthd))
      allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
      allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
      allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
      allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
      allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
      allocate(iboxsubcenters(3,8,nthd))
      allocate(iboxfl(2,8,nthd))
c
c  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot
c  and so on
c
      nmaxt = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
C$OMP$REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(nlist3(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts

          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = iend - istart + 1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
C$OMP END PARALLEL DO

      allocate(iboxsrcind(nmaxt,nthd))
      allocate(iboxsrc(3,nmaxt,nthd))
      allocate(iboxpot(nd,nmaxt,nthd))
      allocate(iboxgrad(nd,3,nmaxt,nthd))
      allocate(iboxhess(nd,6,nmaxt,nthd))

      do ilev=2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1           (2*nterms(ilev)+1),8,nthd))

         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts)
C$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox)

            npts = iend-istart+1

            if(npts.gt.0) then
c            rescale the multipole expansion

                call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1                 rscpow,tmp(1,0,-nmax,ithd))
c
cc                process up down for current box
c
                call mpoletoexp(nd,tmp(1,0,-nmax,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,1),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,2),fexpe,fexpo)


c
cc                process north-south for current box
c
                call rotztoy(nd,nterms(ilev),tmp(1,0,-nmax,ithd),
     1              mptmp(1,ithd),rdminus)
                call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),
     1              nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)

c
cc                process east-west for current box

                call rotztox(nd,nterms(ilev),tmp(1,0,-nmax,ithd),
     1              mptmp(1,ithd),rdplus)
                call mpoletoexp(nd,mptmp(1,ithd),
     1              nterms(ilev),nlams,nfourier,
     2              nexptot,mexpf1(1,1,ithd),
     3              mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,5),fexpe,fexpo)


                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,
     1          nphysical,nthmax,mexp(1,1,ibox,6),fexpe,fexpo)

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
C$OMP$PRIVATE(nuall,ndall,nnall,nsall)
C$OMP$PRIVATE(neall,nwall,nu1234,nd5678)
C$OMP$PRIVATE(nn1256,ns3478,ne1357,nw2468)
C$OMP$PRIVATE(nn12,nn56,ns34,ns78,ne13,ne57)
C$OMP$PRIVATE(nw24,nw68,ne1,ne3,ne5,ne7)
C$OMP$PRIVATE(nw2,nw4,nw6,nw8)
C$OMP$PRIVATE(npts0,ctmp,jstart,jend,i)
C$OMP$PRIVATE(ithd)
         do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
           ithd = 0
C$         ithd=omp_get_thread_num()
           ithd = ithd + 1
           npts = 0
           if(ifpghtarg.gt.0) then
             istart = itargse(1,ibox)
             iend = itargse(2,ibox) 
             npts = npts + iend-istart+1
           endif

           istart = iexpcse(1,ibox) 
           iend = iexpcse(2,ibox) 
           npts = npts + iend-istart+1

           nchild = itree(ipointer(4)+ibox-1)

           if(ifpgh.gt.0) then
             istart = isrcse(1,ibox) 
             iend = isrcse(2,ibox) 
             npts = npts + iend-istart+1
           endif


           if(npts.gt.0.and.nchild.gt.0) then

               call getpwlistall(ibox,boxsize(ilev),nboxes,
     1         itree(ipointer(6)+ibox-1),itree(ipointer(7)+
     2         mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers,
     3         isep,nuall,uall(1,ithd),ndall,dall(1,ithd),
     4         nnall,nall(1,ithd),nsall,sall(1,ithd),
     5         neall,eall(1,ithd),nwall,wall(1,ithd),
     6         nu1234,u1234(1,ithd),nd5678,d5678(1,ithd),
     7         nn1256,n1256(1,ithd),ns3478,s3478(1,ithd),
     8         ne1357,e1357(1,ithd),nw2468,w2468(1,ithd),
     9         nn12,n12(1,ithd),nn56,n56(1,ithd),ns34,s34(1,ithd),
     9         ns78,s78(1,ithd),ne13,e13(1,ithd),ne57,e57(1,ithd),
     9         nw24,w24(1,ithd),nw68,w68(1,ithd),ne1,e1(1,ithd),
     9         ne3,e3(1,ithd),ne5,e5(1,ithd),ne7,e7(1,ithd),
     9         nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9         nw8,w8(1,ithd))

               call processudexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nuall,uall(1,ithd),nu1234,u1234(1,ithd),
     5         ndall,dall(1,ithd),nd5678,d5678(1,ithd),
     6         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     7         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     8         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),xshift,
     8         yshift,zshift,fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,
     9         nlist4,list4,mnlist4)
               
               call processnsexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         nnall,nall(1,ithd),nn1256,n1256(1,ithd),
     5         nn12,n12(1,ithd),nn56,n56(1,ithd),nsall,sall(1,ithd),
     6         ns3478,s3478(1,ithd),ns34,s34(1,ithd),ns78,s78(1,ithd),
     7         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     8         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     9         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),
     9         mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9         mexppall(1,1,7,ithd),
     9         mexppall(1,1,8,ithd),rdplus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,
     9         nlist4,list4,mnlist4)

               
               call processewexp(nd,ibox,ilev,nboxes,centers,
     1         itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2         nterms(ilev),
     2         iaddr,rmlexp,rlams,whts,
     3         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4         neall,eall(1,ithd),ne1357,e1357(1,ithd),
     5         ne13,e13(1,ithd),ne57,e57(1,ithd),ne1,e1(1,ithd),
     6         ne3,e3(1,ithd),ne5,e5(1,ithd),
     7         ne7,e7(1,ithd),nwall,wall(1,ithd),
     8         nw2468,w2468(1,ithd),
     9         nw24,w24(1,ithd),nw68,w68(1,ithd),
     9         nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9         nw8,w8(1,ithd),
     9         mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     9         mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),
     9         mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9         mexppall(1,1,4,ithd),
     9         mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9         mexppall(1,1,7,ithd),mexppall(1,1,8,ithd),
     9         mexppall(1,1,9,ithd),
     9         mexppall(1,1,10,ithd),mexppall(1,1,11,ithd),
     9         mexppall(1,1,12,ithd),
     9         mexppall(1,1,13,ithd),mexppall(1,1,14,ithd),
     9         mexppall(1,1,15,ithd),
     9         mexppall(1,1,16,ithd),rdminus,xshift,yshift,zshift,
     9         fexpback,rlsc,rscpow,
     9         pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)


            endif

            if(nlist3(ibox).gt.0.and.npts.gt.0) then
              call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1             nlist3(ibox),list3(1,ibox),isep,
     2             centers,nuall,uall(1,ithd),ndall,dall(1,ithd),
     3             nnall,nall(1,ithd),
     4             nsall,sall(1,ithd),neall,eall(1,ithd),
     5             nwall,wall(1,ithd))
              do i=1,8
                call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
              enddo

              call processlist3udexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2             whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,nuall,uall(1,ithd),ndall,dall(1,ithd),
     4             mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5             mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6             mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),
     7             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3nsexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2             whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,nnall,nall(1,ithd),nsall,sall(1,ithd),
     4             mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5             mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6             mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdplus,
     7             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3ewexplong(nd,ibox,nboxes,centers,
     1             boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2             whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3             nexptotp,mexp,neall,eall(1,ithd),nwall,wall(1,ithd),
     4             mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     5             mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     6             mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdminus,
     7             xshift,yshift,zshift,fexpback,rlsc,rscpow)

              if(ifpgh.eq.1) then
                istart = isrcse(1,ibox) 
                iend = isrcse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts,
     1                    centers(1,ibox),boxsize(ilev),
     2                    iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                    iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,sourcesort(1,istart),
     1                    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pot(1,istart),
     1                    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1               pot(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpgh.eq.2) then
                istart = isrcse(1,ibox)
                iend = isrcse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts,
     1                    centers(1,ibox),boxsize(ilev),
     2                    iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                    iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,sourcesort(1,istart),
     1                    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pot(1,istart),
     1                    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart),
     1                   iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev),
     1                    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                    nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                    iboxpot(1,jstart,ithd),
     4                    iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1                pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd),
     1              grad(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif
c
c  continue from here
c
              

              if(ifpgh.eq.3) then
                istart = isrcse(1,ibox) 
                iend = isrcse(2,ibox)
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts,
     1                    centers(1,ibox),boxsize(ilev),
     2                    iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                    iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,sourcesort(1,istart),
     1                    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pot(1,istart),
     1                    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,grad(1,1,istart),
     1                   iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hess(1,1,istart),
     1                 iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
           
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),
     4                     iboxgrad(1,1,jstart,ithd),
     4                     iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1                pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd),
     1              grad(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd),
     1              hess(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif


              if(ifpghtarg.eq.1) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts,
     1                   centers(1,ibox),boxsize(ilev),
     2                   iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                   iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,targsort(1,istart),
     1                   iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart),
     1                   iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev),
     1                    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                    nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1              pottarg(1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.2) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts,
     1                   centers(1,ibox),boxsize(ilev),
     2                   iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                   iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,targsort(1,istart),
     1                   iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart),
     1                   iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart),
     1                    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalg(nd,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),
     4                     iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1               pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd),
     1                 gradtarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

              if(ifpghtarg.eq.3) then
                istart = itargse(1,ibox) 
                iend = itargse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(targsort(1,istart),npts,
     1                   centers(1,ibox),boxsize(ilev),
     2                   iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                   iboxsubcenters(1,1,ithd))
                  call dreorderf(ndim,npts,targsort(1,istart),
     1                   iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(nd,npts,pottarg(1,istart),
     1                   iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(3*nd,npts,gradtarg(1,1,istart),
     1                    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,hesstarg(1,1,istart),
     1                 iboxhess(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalh(nd,rscales(ilev),
     1                    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                    nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                    iboxpot(1,jstart,ithd),
     4                    iboxgrad(1,1,jstart,ithd),
     4                    iboxhess(1,1,jstart,ithd),scarray(1,ilev))
                      endif
                    endif
                  enddo
                  call dreorderi(nd,npts,iboxpot(1,1,ithd),
     1               pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(3*nd,npts,iboxgrad(1,1,1,ithd),
     1                 gradtarg(1,1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxhess(1,1,1,ithd),
     1                 hesstarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif

            endif
         enddo
C$OMP END PARALLEL DO        
        deallocate(iboxlexp)  
      enddo

      deallocate(iboxsrcind,iboxsrc,iboxpot,iboxgrad,iboxhess)
      deallocate(iboxsubcenters,iboxfl)
      deallocate(uall,dall,nall,sall,eall,wall)
      deallocate(u1234,d5678,n1256,s3478)
      deallocate(e1357,w2468,n12,n56,s34,s78)
      deallocate(e13,e57,w24,w68)
      deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
      deallocate(tmp,mptmp)
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
               istart = itargse(1,ibox)
               iend = itargse(2,ibox) 
               npts = npts + iend-istart+1
            endif

            istart = iexpcse(1,ibox) 
            iend = iexpcse(2,ibox) 
            npts = npts + iend-istart+1

            if(ifpgh.gt.0) then
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = npts + iend-istart+1
            endif

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
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
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
               istart = iexpcse(1,ibox) 
               iend = iexpcse(2,ibox) 
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
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
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
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox) 
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),grad(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif


        if(ifpgh.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = iend-istart+1
              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),sourcesort(1,istart),
     2         npts,pot(1,istart),grad(1,1,istart),hess(1,1,istart),
     3         scarray(1,ilev))
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
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
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalg(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
            endif
          enddo
C$OMP END PARALLEL DO         
        endif

        if(ifpghtarg.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = iend-istart+1

              call l3dtaevalh(nd,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),gradtarg(1,1,istart),
     3         hesstarg(1,1,istart),scarray(1,ilev))
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
C$OMP$PRIVATE(ibox,istarte,iende,i,jbox)
C$OMP$PRIVATE(jstart,jend)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarte = iexpcse(1,ibox) 
            iende = iexpcse(2,ibox) 

            

            do i =1,nlist1(ibox)
               jbox = list1(i,ibox)
               jstart = isrcse(1,jbox) 
               jend = isrcse(2,jbox) 

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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox) 
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox) 
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart =  isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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


        if(ifpgh.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),
     3             hess(1,1,istarts),thresh)   
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then

C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),
     3             hess(1,1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then

C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,sourcesort(1,istarts),
     2             npts0,pot(1,istarts),grad(1,1,istarts),
     3             hess(1,1,istarts),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif
        endif

        if(ifpghtarg.eq.1) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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

        if(ifpghtarg.eq.3) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectch(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             hesstarg(1,1,istartt),thresh)   
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectdh(nd,sourcesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             hesstarg(1,1,istartt),thresh)          
              enddo
            enddo
C$OMP END PARALLEL DO     
          endif

          if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)     
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
                npts = jend-jstart+1
                call l3ddirectcdh(nd,sourcesort(1,jstart),
     1             chargesort(1,jstart),
     2             dipvecsort(1,1,jstart),npts,targsort(1,istartt),
     2             npts0,pottarg(1,istartt),gradtarg(1,1,istartt),
     3             hesstarg(1,1,istartt),thresh)          
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
c     nd           in: integer(8)
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
c    nlege        in: integer(8)
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
        integer(8) istart,iend,jstart,jend,ns,j, nlege
        integer(8) ifcharge,ifdipole,ier,nd
        double precision source(3,*)
        double precision scj(*)
        double precision wlege(*)
        double precision charge(nd,*)
        double precision dipvec(nd,3,*)
        double precision expc(3,*)

        integer(8) nlevels,ntj
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
