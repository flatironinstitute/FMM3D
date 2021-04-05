c
c   Estimate the memory used by Helmholtz FMM. Note that this is a 
c   rough estimate which may be off by up to 10% in either direction
c
c
c-----------------------------------------------------------
        subroutine hfmm3d_memest(nd,eps,zk,nsource,source,ifcharge,
     $    ifdipole,iper,ifpgh,ntarg,targ,ifpghtarg,rmem)
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
c   ifdipole   in: integer
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
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
c   rmem        out: real *8
c                 memory used in gb
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
      double precision rmem


c       Tree variables
      integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
      integer nlmax
      integer mnbors
      integer ifunif,nlmin
      integer *8 ipointer(8),ltree,lmem8,bigint,bigint0
      integer, allocatable :: itree(:)
      integer, allocatable :: isrcse(:,:),itargse(:,:),isrc(:)
      integer, allocatable :: itarg(:)
      integer, allocatable :: iexpcse(:,:)
      integer iexpc
      double precision, allocatable :: treecenters(:,:),boxsize(:)
      double precision b0,b0inv,b0inv2,b0inv3
      integer mnlist1,mnlist2,mnlist3,mnlist4
      


c
cc       temporary fmm arrays
c
      integer ilevcutoff
      double complex zkfmm,zk2
      double precision epsfmm
      integer, allocatable :: nterms(:)
      integer *8, allocatable :: iaddr(:,:)
      double precision rz,zi,zkiupbound,zkrupbound

      integer lmptemp,nmax
      integer *8 lmptot
c
c    plane wave arrays
c
      integer nlams,nphmax,ntmax,nthmax
      double complex, allocatable :: rlams(:),whts(:)
      integer, allocatable :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp, nn
      integer nexptot_max, nexptotp_max, nn_max


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
      integer nthread
      double precision time1,time2,omp_get_wtime,second
      integer omp_get_max_threads,ntread
      double precision done,pi



       
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
      

c
c   End of tree build
c
c
c   estimate mnlist1,mnlist2,mnlist3,mnlist4
c
      mnbors = 27
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      isep = 1

      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  treecenters,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      

c
c  Set rescaling parameters
c
      b0 = boxsize(0)
      b0inv = 1.0d0/b0
      b0inv2 = b0inv**2

      allocate(nterms(0:nlevels)) 

      zkfmm = zk*b0



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
      call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9

      lmem8 = 100000
c
c    note that things are doubled due to two copies of arrays
c    needed
c
c    source array
c
      bigint = 7
      bigint = bigint*nsource
      lmem8 = lmem8 + bigint 
c
c    charge and dipole arrays
c
      bigint0 = nsource
      bigint0 = bigint0*nd
      if(ifcharge.eq.1) then
        bigint = 4*bigint0
        lmem8 = lmem8 + bigint
      endif
      if(ifdipole.eq.1) then
        bigint = 12*bigint0
        lmem8 = lmem8 + bigint 
      endif

c
c    target array
c
      bigint = 7
      bigint = 7*ntarg
      lmem8 = lmem8 + bigint 
c
c
c     pot, grad, hess arrays
c
      bigint0 = nsource
      bigint0 = bigint0*nd
      if(ifpgh.eq.1) then
        bigint = bigint0*4
        lmem8 = lmem8 + bigint 
      endif
      if(ifpgh.eq.2) then
        bigint = bigint0*12
        lmem8 = lmem8 + bigint 
      endif
      if(ifpgh.eq.3) then
        bigint = bigint0*24
        lmem8 = lmem8 + bigint
      endif
      
      bigint0 = ntarg*nd
      if(ifpghtarg.eq.1) then
        bigint = bigint0*4
        lmem8 = lmem8 + bigint 
      endif
      if(ifpghtarg.eq.1) then
        bigint = bigint0*12
        lmem8 = lmem8 + bigint 
      endif
      if(ifpghtarg.eq.1) then
        bigint = bigint0*24
        lmem8 = lmem8 + bigint 
      endif

      nthread = 1
C$      nthread = omp_get_max_threads()

c
c
c    tree variables
c  
      lmem8 = lmem8 + ltree
      lmem8 = lmem8 + 7*nboxes

c
c   buffer for temporary arrays proportional to nboxes
c
      bigint = 20
      bigint = bigint*nboxes
      bigint = bigint*nthread
      lmem8 = lmem8 + bigint 
c
c   rmlexp array
c
      lmem8 = lmem8 + lmptot

c
c   buffer for mptemp arrays
c
      bigint = 4
      bigint = bigint*lmptemp
      bigint = bigint*nthread
      lmem8 = lmem8 + bigint 
c
c    list1-4
c
      bigint = mnlist1 + mnlist2 + mnlist3 + mnlist4 + 4
      bigint = bigint*nboxes
      lmem8 = lmem8 + bigint 
c
c    pwlist 
c
      lmem8 = lmem8 + 1000*nthread

c
c    temp list3-4 arrays
c
      bigint = 48
      bigint = bigint*(nsource+ntarg)
      bigint = bigint*nthread
      lmem8 = lmem8 + bigint

      bigint = 32
      bigint = bigint*nthread
      bigint = bigint*lmptemp
      lmem8 = lmem8 + bigint 
c
c    plane wave arrays
c
      ntmax = 1000
      allocate(nfourier(ntmax),nphysical(ntmax))
      allocate(rlams(ntmax),whts(ntmax))

      nexptot_max = 0
      nexptotp_max = 0
      nn_max = 0

      done = 1
      pi = atan(done)*4

      zkiupbound = 12*pi
      zkrupbound = 16*pi
      zi = imag(zkfmm)

      ilevcutoff = -1

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
         rz = exp(-zi*boxsize(i))/boxsize(i)
         if(rz.lt.eps) ilevcutoff = i
      enddo


      do ilev=2,nlevels
        zk2 = zkfmm*boxsize(ilev)
        if(real(zk2).le.zkrupbound.and.imag(zk2).lt.zkiupbound.and.
     1     ilev.gt.ilevcutoff) then

          ier = 0
          call hwts3e(ier,eps,zk2,rlams,whts,nlams)
          call hnumfour(eps,zk2,nlams,nfourier)
          call hnumphys(eps,zk2,nlams,nphysical)
          nphmax = 0
          nthmax = 0
          nexptotp = 0
          nexptot = 0
          nn = nterms(ilev)
          do i=1,nlams
            nexptotp = nexptotp + nphysical(i)
            nexptot = nexptot + 2*nfourier(i)+1
            if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
            if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
          enddo
          if(nn.gt.nn_max) nn_max = nn
          if(nexptotp.gt.nexptotp_max) nexptotp_max = nexptotp
          if(nexptot.gt.nexptot_max) nexptot_max = nexptot
        endif
      enddo
      
      bigint = 54
      bigint = bigint*nexptotp_max
      lmem8 = lmem8 + bigint

      bigint = 4*nd
      bigint = bigint*nthread
      bigint = bigint*nexptot_max
      lmem8 = lmem8 + bigint
      bigint = 36*nd
      bigint = bigint*nexptotp_max
      bigint = bigint*nthread
      lmem8 = lmem8 + bigint
c
c  adding buffer of 6*nboxes*nexptotp_max*nd to account for list4
c  processing through plane waves

      bigint = 18
      bigint = bigint*nd
      bigint = bigint*nexptotp_max
      bigint = bigint*nboxes
      lmem8 = lmem8 + bigint 


      lmem8 = lmem8 + 100*(nn_max+1)*(nn_max+1)*(2*nn_max+1)

      rmem = (lmem8+0.0d0)/1024/1024/1024*8
      if(ifprint.ge.1) call prinf('nthread=*',nthread,1)
      if(ifprint.ge.1) call prinf('nexptot_max=*',nexptot_max,1)
      if(ifprint.ge.1) call prinf('nn_max=*',nn_max,1)
      if(ifprint.ge.1) call prinf('nboxes=*',nboxes,1)
      if(ifprint.ge.1) call prinf('nexptotp_max=*',nexptotp_max,1)
      if(ifprint.ge.1) print *, "lmem8=",lmem8
      if(ifprint.ge.1) call prin2('mem required in GB=*',rmem,1)
       

      return
      end
c
