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
     $    charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,
     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
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
c   iper    in: integer
c             flag for periodic implmentations. Currently unused
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
c
c   ier         out: integer
c                error flag
c                ier = 0, for successful execution
c                ier = 4, if failed to allocate workspace
c                      for multipole and local expansions
c                ier = 8, if failed to allocate workspace
c                      for plane wave expansions
c     
c------------------------------------------------------------------

      implicit none

      integer nd
      integer iper
      integer ier

      double complex zk
      double precision eps

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nsource,ntarg

      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,*)

      double complex dipvec(nd,3,*)

      double complex pot(nd,*),grad(nd,3,*),
     1     pottarg(nd,3,*),
     1     gradtarg(nd,3,*),hess(nd,6,*),hesstarg(nd,6,*)

c       Tree variables
      integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
      integer nlmax
      integer mnbors
      integer ifunif,nlmin
      integer *8 ipointer(8),ltree
      integer, allocatable :: itree(:)
      integer, allocatable :: isrcse(:,:),itargse(:,:),isrc(:)
      integer, allocatable :: itarg(:)
      integer, allocatable :: iexpcse(:,:)
      integer iexpc
      double precision, allocatable :: treecenters(:,:),boxsize(:)
      double precision b0,b0inv,b0inv2,b0inv3
      double complex zkfmm

c
cc      temporary sorted arrays
c
      double precision, allocatable :: sourcesort(:,:),targsort(:,:)
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
      integer ntj,nexpc,nadd,ifnear

c
cc        other temporary variables
c
      integer i,iert,ifprint,ilev,idim
      double precision time1,time2,omp_get_wtime,second

       
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c      
      ifprint=0

c
c       turn on computation of list 1
c
      ifnear = 1

c
cc        set criterion for box subdivision
c
      call hndiv(eps,nsource,ntarg,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag) 

      nexpc = 0
      nadd = 0
      ntj = 0
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

      allocate(nterms(0:nlevels)) 
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

c
cc       reorder sources
c
      call dreorderf(3,nsource,source,sourcesort,isrc)
      call drescale(3*nsource,sourcesort,b0inv)
      if(ifcharge.eq.1) then
        call dreorderf(2*nd,nsource,charge,chargesort,isrc)
        call drescale(2*nd*nsource,chargesort,b0inv)
      endif

      if(ifdipole.eq.1) then
         call dreorderf(6*nd,nsource,dipvec,dipvecsort,isrc)
         call drescale(6*nd*nsource,dipvecsort,b0inv2)
      endif

c
cc      reorder targs
c
      call dreorderf(3,ntarg,targ,targsort,itarg)
      call drescale(3*ntarg,targsort,b0inv)
c
c  update tree centers and boxsize
c
      call drescale(3*nboxes,treecenters,b0inv)
      call drescale(nlevels+1,boxsize,b0inv)

      zkfmm = zk*b0

c
c     scaling factor for multipole and local expansions at all levels
c
      allocate(scales(0:nlevels))
      do ilev = 0,nlevels
       scales(ilev) = boxsize(ilev)*abs(zkfmm)
       if(scales(ilev).gt.1) scales(ilev) = 1
      enddo


c     Compute length of expansions at each level      
      nmax = 0
      if(ifprint.ge.1) call prin2('boxsize=*',boxsize,nlevels+1)
      if(ifprint.ge.1) call prin2('zkfmm=*',zkfmm,2)
      do i=0,nlevels
         call h3dterms(boxsize(i),zkfmm,eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
      if(ifprint.ge.1) call prinf('nlevels=*',nlevels,1)
      if(ifprint.ge.1) call prinf('nterms=*',nterms,nlevels+1)
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
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9

      ier = 0
      allocate(rmlexp(lmptot),stat=iert)
      if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         ier = 4
         return
      endif


c     Memory allocation is complete.
c     Call main fmm routine
c
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      call hfmm3dmain(nd,eps,zkfmm,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,radssort,
     $   iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,ndiv,nlevels,
     $   nboxes,iper,boxsize,treecenters,isrcse,itargse,iexpcse,
     $   scales,itree(ipointer(1)),nterms,
     $   ifpgh,potsort,gradsort,hesssort,ifpghtarg,pottargsort,
     $   gradtargsort,hesstargsort,ntj,texpssort,scjsort,ifnear,ier)

      if(ier.ne.0) return


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)


      if(ifpgh.ge.1) then
        call dreorderi(2*nd,nsource,potsort,pot,isrc)
      endif
      if(ifpgh.ge.2) then 
        call dreorderi(6*nd,nsource,gradsort,grad,isrc)
        call drescale(6*nd*nsource,grad,b0inv)
      endif

      if(ifpgh.ge.3) then 
        call dreorderi(12*nd,nsource,hesssort,hess,isrc)
        call drescale(12*nd*nsource,hess,b0inv2)
      endif


      if(ifpghtarg.ge.1) then
        call dreorderi(2*nd,ntarg,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.ge.2) then
        call dreorderi(6*nd,ntarg,gradtargsort,gradtarg,itarg)
        call drescale(6*nd*ntarg,gradtarg,b0inv)
      endif

      if(ifpghtarg.ge.3) then
        call dreorderi(12*nd,ntarg,hesstargsort,hesstarg,itarg)
        call drescale(12*nd*ntarg,hesstarg,b0inv2)
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
     $     itree,ltree,ipointer,ndiv,nlevels, 
     $     nboxes,iper,boxsize,centers,isrcse,itargse,iexpcse,
     $     rscales,laddr,nterms,ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     ntj,jsort,scjsort,ifnear,ier)


      implicit none

      integer nd
      integer ier
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
      integer ifnear
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
      integer isep,iper
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(8),ltree
      integer itree(ltree)
      integer nboxes
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)
      integer isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)

c
cc      pw stuff
c
      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer, allocatable :: uall(:,:),dall(:,:),nall(:,:)
      integer, allocatable :: sall(:,:),eall(:,:),wall(:,:)
      integer, allocatable :: u1234(:,:),d5678(:,:)
      integer, allocatable :: n1256(:,:),s3478(:,:)
      integer, allocatable :: e1357(:,:),w2468(:,:)
      integer, allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      integer, allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      integer, allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      integer, allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

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
      double complex, allocatable :: tmp(:,:,:,:),tmp2(:,:,:,:)
      double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, allocatable :: mexpp1(:,:,:),mexpp2(:,:,:),
     1    mexppall(:,:,:,:)

      double precision, allocatable :: rsc(:)
      double precision r1

      double precision scjsort(nexpc),radssort(nexpc)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll,idim
      integer ibox,jbox,ilev,npts,npts0
      integer nchild

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
      double precision, allocatable :: wlege(:)

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

c     list 3 variables
      double complex, allocatable :: iboxlexp(:,:,:)
      double precision, allocatable :: iboxsubcenters(:,:,:)
      double complex, allocatable :: iboxpot(:,:,:)
      double complex, allocatable :: iboxgrad(:,:,:,:)
      double precision, allocatable :: iboxsrc(:,:,:)
      integer, allocatable :: iboxsrcind(:,:)
      integer, allocatable :: iboxfl(:,:,:)
c     end of list 3 variables
c     list 4 variables
      integer cntlist4
      integer, allocatable :: list4ct(:),ilist4(:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)


c     end of list 4 variables

      integer *8 bigint
      double precision zkiupbound,zi,zkrupbound,rz
      integer ilevcutoff

      integer iert
      data ima/(0.0d0,1.0d0)/

      integer nlfbox
      integer nthd,ithd
      integer omp_get_max_threads,omp_get_thread_num
      nthd = 1
C$    nthd=omp_get_max_threads()


      ntmax = 1000
      allocate(nfourier(ntmax),nphysical(ntmax))
      allocate(rlams(ntmax),whts(ntmax))

      pi = 4.0d0*atan(1.0d0)
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
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
      

c
c      If imaginary part is greater than 12*pi
c    don't do any multipole/local work
c

      zkiupbound = 12*pi
      zkrupbound = 16*pi
      zi = imag(zk)

      ilevcutoff = -1

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
         rz = exp(-zi*boxsize(i))/boxsize(i)
         if(rz.lt.eps) ilevcutoff = i
      enddo

      if(ifprint.ge.1) print *, "ilevcutoff=",ilevcutoff

      allocate(rsc(0:nmax))

c
cc     threshold for computing interactions,
c      interactions will be ignored
c      for all pairs of sources and targets
c      which satisfy |r| < thresh
c      where r is the disance between them

      thresh = 2.0d0**(-51)*boxsize(0)

      if(ifprint.ge.1) print *, "thresh=",thresh

      

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
            nchild = itree(ipointer(4)+ibox-1)
            if(nchild.gt.0) then
               istart = iexpcse(1,ibox)
               iend = iexpcse(2,ibox)
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
      nlege = nmax + 10
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege,wlege,lw7,lused7)

      allocate(list4ct(nboxes))
      allocate(ilist4(nboxes))
      do i=1,nboxes
        list4ct(i)=0
        ilist4(i)=0
      enddo
      cntlist4=0

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
        if(ilev.gt.ilevcutoff) then
          if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,npts,istart,iend,nchild)
            do ibox=laddr(1,ilev),laddr(2,ilev)
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = iend-istart+1

              nchild = itree(ipointer(4)+ibox-1)

              if(npts.gt.0.and.nchild.eq.0) then
                call h3dformmpc(nd,zk,rscales(ilev),
     1          sourcesort(1,istart),chargesort(1,istart),npts,
     2          centers(1,ibox),nterms(ilev),
     3          rmlexp(iaddr(1,ibox)),wlege,nlege)          
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

              if(npts.gt.0.and.nchild.eq.0) then
                call h3dformmpd(nd,zk,rscales(ilev),
     1          sourcesort(1,istart),
     2          dipvecsort(1,1,istart),npts,
     2          centers(1,ibox),nterms(ilev),
     3          rmlexp(iaddr(1,ibox)),wlege,nlege)          
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

              if(npts.gt.0.and.nchild.eq.0) then
                call h3dformmpcd(nd,zk,rscales(ilev),
     1          sourcesort(1,istart),chargesort(1,istart),
     2          dipvecsort(1,1,istart),npts,
     2          centers(1,ibox),nterms(ilev),
     3          rmlexp(iaddr(1,ibox)),wlege,nlege)          
              endif
            enddo
C$OMP END PARALLEL DO          
          endif
        endif
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c


      do ilev=nlevels-1,1,-1
        if(ilev.gt.ilevcutoff) then
          nquad2 = nterms(ilev)*2.5
          nquad2 = max(6,nquad2)
          ifinit2 = 1
          call legewhts(nquad2,xnodes,wts,ifinit2)
          radius = boxsize(ilev)/2*sqrt(3.0d0)

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
                  call h3dmpmp(nd,zk,rscales(ilev+1),
     1              centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2              nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3              rmlexp(iaddr(1,ibox)),nterms(ilev),
     4              radius,xnodes,wts,nquad2)
                endif
              endif
            enddo
          enddo
C$OMP END PARALLEL DO
        endif
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1



      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc+mpeval+formta) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions and big to small far and small to far big

      call cpu_time(time1)
C$    time1=omp_get_wtime()
c     init uall,dall,...,etc arrays
      allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
      allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
      allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
      allocate(s3478(24,nthd))
      allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
      allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
      allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
      allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd),w2(5,nthd))
      allocate(w4(5,nthd),w6(5,nthd),w8(5,nthd))


c
c  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot
c  and so on
c
      nmax = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
C$OMP$REDUCTION(max:nmax)
      do ibox=1,nboxes
        if(nlist3(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmax) nmax = npts

          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = iend - istart + 1
          if(npts.gt.nmax) nmax = npts
        endif
      enddo
C$OMP END PARALLEL DO

      allocate(iboxsrcind(nmax,nthd))
      allocate(iboxsrc(3,nmax,nthd))
      allocate(iboxpot(nd,nmax,nthd))
      allocate(iboxgrad(nd,3,nmax,nthd))
      allocate(iboxsubcenters(3,8,nthd))
      allocate(iboxfl(2,8,nthd))




      do ilev = 2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1           (2*nterms(ilev)+1),8,nthd))
        zk2 = zk*boxsize(ilev)
        if(real(zk2).le.zkrupbound.and.imag(zk2).lt.zkiupbound.and.
     1        ilev.gt.ilevcutoff) then
c             get new pw quadrature
            
          ier = 0
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
          allocate(tmp(nd,0:nterms(ilev),
     1            -nterms(ilev):nterms(ilev),nthd))
          allocate(tmp2(nd,0:nterms(ilev),
     1            -nterms(ilev):nterms(ilev),nthd))
 
          allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd),
     1          mexpp1(nd,nexptotp,nthd))
          allocate(mexpp2(nd,nexptotp,nthd),
     1          mexppall(nd,nexptotp,16,nthd))

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
c         note: the scaling for helmholtz has been eliminated
c         since it is taken care in the scaling of the legendre
c         functions
c
          r1 = 1.0d0
          rsc(0) = 1.0d0
          do i=1,nterms(ilev)
            rsc(i) = rsc(i-1)*r1
          enddo

c      generate ilev+1 list4 type box plane wave expansion
          cntlist4=0
          do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
            
            if(nlist3(ibox).gt.0) then
              cntlist4=cntlist4+1
              list4ct(ibox)=cntlist4
            endif
          enddo
          allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
          if(ifprint.ge.1) print *,"cntlist4:",cntlist4,"ilev:",ilev

          call h3dlist4pw(ilev-1,zk,nd,nexptotp,nexptot,nterms(ilev),
     1       nn,nlams,nlege,nlevels,ifcharge,ifdipole,list4ct,isrcse,
     2       laddr,nfourier,nphysical,
     3       rdminus,rdplus,rlsc,
     4       rscales(ilev),boxsize(ilev),xshift,yshift,zshift,
     5       sourcesort,chargesort,dipvecsort,centers,fexp,
     6       mexpf1,mexpf2,tmp,tmp2,wlege,rlams,rsc,pgboxwexp,
     7       cntlist4)

c
cc         create multipole to plane wave expansion for
c          all boxes at this level
c
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts)
C$OMP$PRIVATE(ithd)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            ithd = 0
C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend - istart+1
            if(npts.gt.0) then

c           rescale multipole expansion
              call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1               rsc,tmp(1,0,-nterms(ilev),ithd))
                
              call hmpoletoexp(nd,tmp(1,0,-nterms(ilev),ithd),
     1          nterms(ilev),nlams,nfourier,nexptot,mexpf1(1,1,ithd),
     2          mexpf2(1,1,ithd),rlsc) 

              call hftophys(nd,mexpf1(1,1,ithd),nlams,nfourier,
     1             nphysical,mexp(1,1,ibox,1),fexp)           

              call hftophys(nd,mexpf2(1,1,ithd),nlams,nfourier,
     1             nphysical,mexp(1,1,ibox,2),fexp)


c          form mexpnorth, mexpsouth for current box

c          Rotate mpole for computing mexpnorth and
c          mexpsouth
              call rotztoy(nd,nterms(ilev),tmp(1,0,-nterms(ilev),ithd),
     1                           tmp2(1,0,-nterms(ilev),ithd),rdminus)

              call hmpoletoexp(nd,tmp2(1,0,-nterms(ilev),ithd),
     1                  nterms(ilev),nlams,
     2                  nfourier,nexptot,mexpf1(1,1,ithd),
     3                  mexpf2(1,1,ithd),rlsc)

              call hftophys(nd,mexpf1(1,1,ithd),nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,3),fexp)           

              call hftophys(nd,mexpf2(1,1,ithd),nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,4),fexp)   


c         Rotate mpole for computing mexpeast, mexpwest
              call rotztox(nd,nterms(ilev),tmp(1,0,-nterms(ilev),ithd),
     1                            tmp2(1,0,-nterms(ilev),ithd),rdplus)
              call hmpoletoexp(nd,tmp2(1,0,-nterms(ilev),ithd),
     1                  nterms(ilev),nlams,
     2                  nfourier,nexptot,mexpf1(1,1,ithd),
     3                  mexpf2(1,1,ithd),rlsc)

              call hftophys(nd,mexpf1(1,1,ithd),nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,5),fexp)

              call hftophys(nd,mexpf2(1,1,ithd),nlams,nfourier,
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
C$          ithd=omp_get_thread_num()
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
     1           itree(ipointer(6)+ibox-1),itree(ipointer(7)+
     2           mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers,
     3           isep,nuall,uall(1,ithd),ndall,dall(1,ithd),
     4           nnall,nall(1,ithd),nsall,sall(1,ithd),
     5           neall,eall(1,ithd),nwall,wall(1,ithd),
     6           nu1234,u1234(1,ithd),nd5678,d5678(1,ithd),
     7           nn1256,n1256(1,ithd),ns3478,s3478(1,ithd),
     8           ne1357,e1357(1,ithd),nw2468,w2468(1,ithd),
     9           nn12,n12(1,ithd),nn56,n56(1,ithd),ns34,s34(1,ithd),
     9           ns78,s78(1,ithd),ne13,e13(1,ithd),ne57,e57(1,ithd),
     9           nw24,w24(1,ithd),nw68,w68(1,ithd),ne1,e1(1,ithd),
     9           ne3,e3(1,ithd),ne5,e5(1,ithd),ne7,e7(1,ithd),
     9           nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9           nw8,w8(1,ithd))


              call hprocessudexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nuall,uall(1,ithd),nu1234,u1234(1,ithd),
     5            ndall,dall(1,ithd),nd5678,d5678(1,ithd),
     6            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     7            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     8            mexppall(1,1,1,ithd),
     9            mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9            mexppall(1,1,4,ithd),
     9            xshift,yshift,zshift,fexpback,rlsc,
     9            pgboxwexp,cntlist4,list4ct,
     9            nlist4,list4,mnlist4)


              call hprocessnsexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nnall,nall(1,ithd),nn1256,n1256(1,ithd),
     5            nn12,n12(1,ithd),nn56,n56(1,ithd),nsall,sall(1,ithd),
     5            ns3478,s3478(1,ithd),ns34,s34(1,ithd),
     6            ns78,s78(1,ithd),
     7            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     8            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     9            mexppall(1,1,1,ithd),
     9            mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9            mexppall(1,1,4,ithd),
     9            mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9            mexppall(1,1,7,ithd),
     9            mexppall(1,1,8,ithd),rdplus,xshift,yshift,zshift,
     9            fexpback,rlsc,
     9            pgboxwexp,cntlist4,list4ct,
     9            nlist4,list4,mnlist4)

              call hprocessewexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(ipointer(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            neall,eall(1,ithd),ne1357,e1357(1,ithd),
     5            ne13,e13(1,ithd),ne57,e57(1,ithd),ne1,e1(1,ithd),
     5            ne3,e3(1,ithd),ne5,e5(1,ithd),
     6            ne7,e7(1,ithd),nwall,wall(1,ithd),
     7            nw2468,w2468(1,ithd),
     8            nw24,w24(1,ithd),nw68,w68(1,ithd),
     9            nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),
     9            nw8,w8(1,ithd),
     9            mexpf1(1,1,ithd),mexpf2(1,1,ithd),
     9            mexpp1(1,1,ithd),mexpp2(1,1,ithd),
     9            mexppall(1,1,1,ithd),
     9            mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),
     9            mexppall(1,1,4,ithd),
     9            mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),
     9            mexppall(1,1,7,ithd),mexppall(1,1,8,ithd),
     9            mexppall(1,1,9,ithd),
     9            mexppall(1,1,10,ithd),mexppall(1,1,11,ithd),
     9            mexppall(1,1,12,ithd),
     9            mexppall(1,1,13,ithd),mexppall(1,1,14,ithd),
     9            mexppall(1,1,15,ithd),
     9            mexppall(1,1,16,ithd),rdminus,xshift,yshift,zshift,
     9            fexpback,rlsc,pgboxwexp,cntlist4,list4ct,
     9            nlist4,list4,mnlist4)
            endif


c
c      handle mp eval
c

            if(nlist3(ibox).gt.0.and.npts.gt.0) then
              call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1              nlist3(ibox),list3(1,ibox),isep,
     2              centers,nuall,uall(1,ithd),ndall,dall(1,ithd),
     3              nnall,nall(1,ithd),
     4              nsall,sall(1,ithd),neall,eall(1,ithd),
     5              nwall,wall(1,ithd))
                 

              do i=1,8
                call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
              enddo

              call hprocesslist3udexplong(nd,zk2,ibox,nboxes,centers,
     1               boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),
     2               rlams,whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3               nexptotp,mexp,nuall,uall(1,ithd),
     4               ndall,dall(1,ithd),
     5               mexpf1(1,1,ithd),mexpf2(1,1,ithd),mexpp1(1,1,ithd),
     6               mexpp2(1,1,ithd),
     7               mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),
     8               xshift,yshift,zshift,fexpback,rlsc)

              call hprocesslist3nsexplong(nd,zk2,ibox,nboxes,centers,
     1           boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,
     2           whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3           nexptotp,mexp,nnall,nall(1,ithd),nsall,sall(1,ithd),
     4           mexpf1(1,1,ithd),mexpf2(1,1,ithd),mexpp1(1,1,ithd),
     5           mexpp2(1,1,ithd),
     6           mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdplus,
     7           xshift,yshift,zshift,fexpback,rlsc)

              call hprocesslist3ewexplong(nd,zk2,ibox,nboxes,centers,
     1              boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),
     2              rlams,whts,nlams,nfourier,nphysical,nthmax,nexptot,
     3              nexptotp,mexp,neall,eall(1,ithd),nwall,wall(1,ithd),
     4              mexpf1(1,1,ithd),mexpf2(1,1,ithd),mexpp1(1,1,ithd),
     6              mexpp2(1,1,ithd),
     7              mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdminus,
     8              xshift,yshift,zshift,fexpback,rlsc)

              if(ifpgh.eq.1) then
                istart = isrcse(1,ibox)
                iend = isrcse(2,ibox)
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts,
     1                    centers(1,ibox),boxsize(ilev),
     2                    iboxsrcind(1,ithd),iboxfl(1,1,ithd),
     3                    iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,sourcesort(1,istart),
     1                    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(2*nd,npts,pot(1,istart),
     1                    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call h3dtaevalp(nd,zk,rscales(ilev),
     1                    iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                    nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                    iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(2*nd,npts,iboxpot(1,1,ithd),
     1              pot(1,istart),iboxsrcind(1,ithd))
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
                  call dreorderf(3,npts,sourcesort(1,istart),
     1                    iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(2*nd,npts,pot(1,istart),
     1                    iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,grad(1,1,istart),
     1                   iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call h3dtaevalg(nd,zk,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),
     4                     iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(2*nd,npts,iboxpot(1,1,ithd),
     1               pot(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxgrad(1,1,1,ithd),
     1               grad(1,1,istart),iboxsrcind(1,ithd))
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
                  call dreorderf(3,npts,targsort(1,istart),
     1                   iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(2*nd,npts,pottarg(1,istart),
     1                   iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call h3dtaevalp(nd,zk,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(2*nd,npts,iboxpot(1,1,ithd),
     1                pottarg(1,istart),iboxsrcind(1,ithd))
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
                  call dreorderf(3,npts,targsort(1,istart),
     1                   iboxsrc(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(2*nd,npts,pottarg(1,istart),
     1                   iboxpot(1,1,ithd),iboxsrcind(1,ithd))
                  call dreorderf(6*nd,npts,gradtarg(1,1,istart),
     1                    iboxgrad(1,1,1,ithd),iboxsrcind(1,ithd))
                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call h3dtaevalg(nd,zk,rscales(ilev),
     1                     iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),
     2                     nterms(ilev),iboxsrc(1,jstart,ithd),npts0,
     3                     iboxpot(1,jstart,ithd),
     4                     iboxgrad(1,1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                  call dreorderi(2*nd,npts,iboxpot(1,1,ithd),
     1               pottarg(1,istart),iboxsrcind(1,ithd))
                  call dreorderi(6*nd,npts,iboxgrad(1,1,1,ithd),
     1                    gradtarg(1,1,istart),iboxsrcind(1,ithd))
                endif
              endif
            endif
          enddo
C$OMP END PARALLEL DO        

          deallocate(xshift,yshift,zshift,rlsc,tmp,tmp2)
          deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)
          deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)
          deallocate(fexp,fexpback)
          deallocate(pgboxwexp)


        else if((real(zk2).gt.zkrupbound.or.imag(zk2).gt.zkiupbound).
     1            and.ilev.gt.ilevcutoff) then
          nquad2 = nterms(ilev)*2.2
          if(ifprint.ge.1) print *, "In point and shoot regime"
          nquad2 = max(6,nquad2)

          ifinit2 = 1
          ier = 0

          call legewhts(nquad2,xnodes,wts,ifinit2)

          radius = boxsize(ilev)/2*sqrt(3.0d0)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,i,jbox)
          do ibox = laddr(1,ilev),laddr(2,ilev)
            npts = 0
            if(ifpghtarg.gt.0) then
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)
              npts = npts + iend - istart + 1
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
              do i=1,nlist2(ibox)
                jbox = list2(i,ibox)

                istart = isrcse(1,jbox)
                iend = isrcse(2,jbox)
                npts = iend-istart+1

                if(npts.gt.0) then
                  call h3dmploc(nd,zk,rscales(ilev),centers(1,jbox),
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
C$OMP$PRIVATE(ibox,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = iend-istart+1

              do i=1,nlist3(ibox)
                jbox = list3(i,ibox) 
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
C$OMP$PRIVATE(ibox,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
              
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)

              npts = iend-istart+1

              do i=1,nlist3(ibox)
                jbox = list3(i,ibox)
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
C$OMP$PRIVATE(ibox,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
              
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)

              npts = iend-istart+1

              do i=1,nlist3(ibox)
                jbox = list3(i,ibox)
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
C$OMP$PRIVATE(ibox,istart,iend,npts,i,jbox)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev-1),laddr(2,ilev-1)
              
              istart = itargse(1,ibox)
              iend = itargse(2,ibox)

              npts = iend-istart+1

              do i=1,nlist3(ibox)
                jbox = list3(i,ibox)
                call h3dmpevalg(nd,zk,rscales(ilev),centers(1,jbox),
     1             rmlexp(iaddr(1,jbox)),nterms(ilev),
     2             targsort(1,istart),npts,pottarg(1,istart),
     3             gradtarg(1,1,istart),wlege,nlege,thresh)
              enddo
            enddo
C$OMP END PARALLEL DO
          endif
        endif
        deallocate(iboxlexp)
      enddo

c
c    handle list 4 interactions not handled by plane waves
c    due to high frequency
c
    
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        do ilev=1,nlevels
          zk2 = zk*boxsize(ilev)
          if((real(zk2).gt.zkrupbound.or.imag(zk2).gt.zkiupbound).
     1            and.ilev.gt.ilevcutoff) then

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
              
              do i=1,nlist4(ibox)
                jbox = list4(i,ibox)

c              Form local expansion for all boxes in list3
c              of the current box


                istart = isrcse(1,jbox)
                iend = isrcse(2,jbox)
                npts = iend-istart+1
                if(npts.gt.0) then
                  call h3dformtac(nd,zk,rscales(ilev),
     1              sourcesort(1,istart),chargesort(1,istart),npts,
     2              centers(1,ibox),nterms(ilev),
     3              rmlexp(iaddr(2,ibox)),wlege,nlege)
                endif
              enddo
            enddo
C$OMP END PARALLEL DO
          endif
        enddo
      endif

    
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        do ilev=1,nlevels
          zk2 = zk*boxsize(ilev)
          if((real(zk2).gt.zkrupbound.or.imag(zk2).gt.zkiupbound).
     1            and.ilev.gt.ilevcutoff) then

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
              
              do i=1,nlist4(ibox)
                jbox = list4(i,ibox)

c              Form local expansion for all boxes in list3
c              of the current box


                istart = isrcse(1,jbox)
                iend = isrcse(2,jbox)
                npts = iend-istart+1
                if(npts.gt.0) then
                  call h3dformtad(nd,zk,rscales(ilev),
     1              sourcesort(1,istart),dipvecsort(1,1,istart),npts,
     2              centers(1,ibox),nterms(ilev),
     3              rmlexp(iaddr(2,ibox)),wlege,nlege)
                endif
              enddo
            enddo
C$OMP END PARALLEL DO
          endif
        enddo
      endif

    
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        do ilev=1,nlevels
          zk2 = zk*boxsize(ilev)
          if((real(zk2).gt.zkrupbound.or.imag(zk2).gt.zkiupbound).
     1            and.ilev.gt.ilevcutoff) then

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
              
              do i=1,nlist4(ibox)
                jbox = list4(i,ibox)

c              Form local expansion for all boxes in list3
c              of the current box


                istart = isrcse(1,jbox)
                iend = isrcse(2,jbox)
                npts = iend-istart+1
                if(npts.gt.0) then
                  call h3dformtacd(nd,zk,rscales(ilev),
     1              sourcesort(1,istart),chargesort(1,istart),
     2              dipvecsort(1,1,istart),npts,
     2              centers(1,ibox),nterms(ilev),
     3              rmlexp(iaddr(2,ibox)),wlege,nlege)
                endif
              enddo
            enddo
C$OMP END PARALLEL DO
          endif
        enddo
      endif

      deallocate(uall,dall,nall,sall,eall,wall,u1234,d5678,n1256)
      deallocate(s3478,e1357,w2468,n12,n56,s34,s78)
      deallocate(e13,e57,w24,w68,e1,e3,e5,e7,w2,w4,w6,w8)
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 1,nlevels-1
        if(ilev.gt.ilevcutoff) then
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
        endif
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval lo) ===*',i,0)

c     ... step 5, evaluate all local expansions
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
        if(ilev.gt.ilevcutoff) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
C$OMP$SCHEDULE(DYNAMIC)      
          do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              do i=istart,iend
                call h3dlocloc(nd,zk,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),rscales(ilev),expcsort(1,i),
     3          jsort(1,0,-ntj,i),ntj,radssort(i),xnodes,wts,nquad2)
              enddo
            endif
          enddo
C$OMP END PARALLEL DO
        endif
      enddo

c
cc        evaluate local expansion at source and target
c         locations
c
      do ilev = 0,nlevels
        if(zi*boxsize(ilev).lt.zkiupbound) then
          if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              nchild=itree(ipointer(4)+ibox-1)
              if(nchild.eq.0) then 
                istart = isrcse(1,ibox)
                iend = isrcse(2,ibox)
                npts = iend-istart+1
                call h3dtaevalp(nd,zk,rscales(ilev),centers(1,ibox),
     1           rmlexp(iaddr(2,ibox)),nterms(ilev),
     2           sourcesort(1,istart),npts,pot(1,istart),wlege,nlege)
              endif
            enddo
C$OMP END PARALLEL DO          
          endif

          if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)      
            do ibox = laddr(1,ilev),laddr(2,ilev)
              nchild=itree(ipointer(4)+ibox-1)
              if(nchild.eq.0) then 
                istart = isrcse(1,ibox)
                iend = isrcse(2,ibox)
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
              nchild=itree(ipointer(4)+ibox-1)
              if(nchild.eq.0) then 
                istart = itargse(1,ibox)
                iend = itargse(2,ibox)
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
              nchild=itree(ipointer(4)+ibox-1)
              if(nchild.eq.0) then 
                istart = itargse(1,ibox)
                iend = itargse(2,ibox)
                npts = iend-istart+1

                call h3dtaevalg(nd,zk,rscales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),targsort(1,istart),
     2         npts,pottarg(1,istart),gradtarg(1,1,istart),wlege,nlege)
              endif
            enddo
C$OMP END PARALLEL DO         
          endif
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istarts,iends,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istarts = isrcse(1,ibox)
              iends = isrcse(2,ibox)
              npts0 = iends-istarts+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              

              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
C$OMP$PRIVATE(ibox,istartt,iendt,npts0,i,jbox,jstart,jend,npts)
            do ibox = laddr(1,ilev),laddr(2,ilev)
              istartt = itargse(1,ibox)
              iendt = itargse(2,ibox)
              npts0 = iendt-istartt+1
              
              do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                jstart = isrcse(1,jbox)
                jend = isrcse(2,jbox)
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
