c       
c       Planning interfaces for the 
c       generalized helmholtz FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets 
c
c       We use exp(ikr)/r for the Green's function., without
c       the 1/(4\pi ) scaling.
c
c   
c-----------------------------------------------------------
        subroutine hfmm3d_plan(eps,zk,nsource,source,
     $    ifpgh,ntarg,targ,ifpghtarg,iptype,ltree,nboxes,nlevels,
     $    nwork)
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
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
c   iptype    in: integer
c             plan type
c             iptype = 1, store only direct interactions
c             iptype = 2, store direct interactions + 
c                         local evaluation
c             iptype = 3, store direct interactions + 
c                          local evaluation +
c                          multipole evaluation
c
c
c   OUTPUT parameters:
c   ltree     out: integer *8
c             length of tree to be allocated
c
c   nboxes    out: integer
c             number of boxes
c
c   nlevels   out: integer
c             number of levels in tree hierarchy
c
c   nwork     out: integer *8 
c             length of work array to be allocated
c
c
     
c------------------------------------------------------------------

      implicit none

      integer nd

      double complex zk
      double precision eps

      integer ifpgh,ifpghtarg

      integer nsource,ntarg

      integer iptype

      double precision source(3,nsource),targ(3,ntarg)
      integer *8 nwork, ltree
      integer nboxes,nlevels

c       Tree variables
      integer mhung,idivflag,ndiv,isep,nbmax
      integer nlmax
      integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
      integer *8 ipointer(32)
      integer, allocatable :: itree(:)
      double precision, allocatable :: treecenters(:,:),boxsize(:)
      double precision, allocatable :: radsrc(:)

c
cc      temporary sorted arrays
c
      double precision, allocatable :: sourcesort(:,:),targsort(:,:)
      integer isfac,itfac
      

c
cc       temporary fmm arrays
c
      integer, allocatable :: nterms(:)

c
cc       temporary variables not used in particle code
c
      double precision expc(3),radexp

c
cc        other temporary variables
c
       integer i,iert,ifprint,ilev,idim,ier
       integer nptss,nptst
       integer istarts,istartt,iends,iendt
       integer npts,nn,nmax,nlist1,nlist3,nexpc
       integer ibox,jbox,jstart,jend,nchild

c
cc       figure out tree structure
c
   
c
cc        set criterion for box subdivision
c
       ndiv = 100
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

        if(iert.ne.0) then
           print *, "Error in allocating tree memory in plan routine"
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



c     Compute length of expansions at each level     
      allocate(nterms(0:nlevels))
      nmax = 0
      do i=0,nlevels
         call h3dterms(boxsize(i),zk,eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

c
cc       reorder sources
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
c
cc      reorder targs
c
      call dreorderf(3,ntarg,targ,targsort,itree(ipointer(6)))

      print *, iptype
      if(iptype.ne.1.and.iptype.ne.2.and.iptype.ne.3) then
        print *, "Plan type must be 1,2 or 3, returning without plan"
        print *, "setting iptype = 0, to continue execution"
        iptype = 0
        return
      endif

      isfac = 0
      if(ifpgh.eq.1) isfac = 1
      if(ifpgh.eq.2) isfac = 4

      itfac = 0
      if(ifpghtarg.eq.1) itfac = 1
      if(ifpghtarg.eq.2) itfac = 4

      nwork = 0
      

      do ilev=0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istarts,iends,nptss,nlist1,i,jbox)
C$OMP$PRIVATE(istartt,iendt,nptst)
C$OMP$PRIVATE(jstart,jend)
C$OMP$REDUCTION(+:nwork)
        do ibox=itree(2*ilev+1),itree(2*ilev+2)
          istarts = itree(ipointer(10)+ibox-1)
          iends = itree(ipointer(11)+ibox-1)
          nptss = iends-istarts+1
          
          istartt = itree(ipointer(12)+ibox-1)
          iendt = itree(ipointer(13)+ibox-1)
          nptst = iendt-istartt+1
          
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
            jstart = itree(ipointer(10)+jbox-1)
            jend = itree(ipointer(11)+jbox-1)
            npts = jend-jstart+1

            nwork = nwork + 2*npts*(isfac*nptss+itfac*nptst)
          enddo
        enddo
C$OMP END PARALLEL DO
      enddo

      print *, "nwork (after direct)= ", nwork/1.0d9



      if(iptype.ge.2) then
        do ilev=0,nlevels
          nn = (2*nterms(ilev)+1)*2*(nterms(ilev)+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,istarts,iends,nptss)
C$OMP$PRIVATE(istartt,iendt,nptst)
C$OMP$REDUCTION(+:nwork)
          do ibox=itree(2*ilev+1),itree(2*ilev+2)
             nchild = itree(ipointer(3)+ibox-1)
             if(nchild.eq.0) then
               istarts = itree(ipointer(10)+ibox-1)
               iends = itree(ipointer(11)+ibox-1)
               nptss = iends-istarts+1
          
               istartt = itree(ipointer(12)+ibox-1)
               iendt = itree(ipointer(13)+ibox-1)
               nptst = iendt-istartt+1

               nwork = nwork + nn*(nptss*isfac + nptst*itfac)
             endif
          enddo
C$OMP END PARALLEL DO          
        enddo
      endif

      print *, "nwork (after taeval)= ", nwork/1.0d9


      if(iptype.eq.3) then
        do ilev=0,nlevels
          nn = (2*nterms(ilev)+1)*2*(nterms(ilev)+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istarts,iends,nptss)
C$OMP$PRIVATE(istartt,iendt,nptst)
C$OMP$REDUCTION(+:nwork)
          do ibox=itree(2*ilev+1),itree(2*ilev+2)
             nchild = itree(ipointer(3)+ibox-1)
             if(nchild.eq.0) then
               istarts = itree(ipointer(10)+ibox-1)
               iends = itree(ipointer(11)+ibox-1)
               nptss = iends-istarts+1
          
               istartt = itree(ipointer(12)+ibox-1)
               iendt = itree(ipointer(13)+ibox-1)
               nptst = iendt-istartt+1

               nlist3 = itree(ipointer(24)+ibox-1)

               nwork = nwork + nn*(nptss*isfac + nptst*itfac)*nlist3
             endif
          enddo
C$OMP END PARALLEL DO          
        enddo
      endif

      print *, "nwork (after mpeval)= ", nwork/1.0d9

      return
      end
c
c
