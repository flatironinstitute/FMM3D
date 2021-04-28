c
c
c
c
c   generate level restricted tree based on resolving points
c    (Currently only supports, sorting on sources, sorting on
c      targets, or sorting on sources and targets)
c
c   There is additional functionality to further subdivide
c    an existing tree based on resolving a function, 
c      this functionality tree will be added in later
c
c   This code has the following user callable routines
c
c      pts_tree_mem -> returns memory requirements for creating
c         a tree based on max number of sources/targets
c         in a box (tree length
c         number of boxes, number of levels)
c      pts_tree -> Make the actual tree, returns centers of boxes,
c        colleague info, pts sorted on leaf boxes
c
c      iptr(1) - laddr
c      iptr(2) - ilevel
c      iptr(3) - iparent
c      iptr(4) - nchild
c      iptr(5) - ichild
c      iptr(6) - ncoll
c      iptr(7) - coll
c      iptr(8) - ltree
c   
c 


      subroutine pts_tree_mem(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,
     1    ifunif,iper,nlevels,nboxes,ltree)
c
c
c
c----------------------------------------
c  get memory requirements for the tree
c
c
c  input parameters:
c    - src: real *8 (3,ns)
c        source locations
c    - targ: real *8 (3,nt) 
c        target locations
c    - idivflag: integer
c        subdivision criterion
c          * divflag = 0 -> subdivide on sources only
c          * idivflag = 1 -> subdivide on targets only
c          * idivflag = 2 -> subdivide on max(sources+targets)
c    - ndiv: integer
c        subdivide if relevant number of particles
c        per box is greater than ndiv
c    - nlmin: integer
c        minimum number of levels of uniform refinement.
c        Note that empty boxes are not pruned along the way
c    - nlmax: integer
c        max number of levels
c    - ifunif: integer
c        flag for creating uniform pruned tree
c        Tree is uniform if ifunif=1 (Currently pruned part
c        under construction)
c    - iper: integer
c        flag for periodic implementations. Currently unused.
c        Feature under construction
c
c        
c  output parameters
c    - nlevels: integer
c        number of levels
c    - nboxes: integer
c        number of boxes
c    - ltree: integer
c        length of tree
c----------------------------------
c
     

      implicit none
      integer nlevels,nboxes,idivflag
      integer *8 ltree,nboxes8
      integer nbmax,nbtot
      integer ns,nt,ndiv
      integer nlmin,iper,ifunif
      double precision src(3,ns),targ(3,nt)


      integer, allocatable :: laddr(:,:),ilevel(:),iparent(:),nchild(:)
      integer, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      double precision, allocatable :: centers(:,:)
      integer, allocatable :: nbors(:,:),nnbors(:)

      integer, allocatable :: isrc(:),itarg(:),isrcse(:,:),itargse(:,:)

      integer, allocatable :: ilevel2(:),iparent2(:),nchild2(:),
     1    ichild2(:,:),isrcse2(:,:),itargse2(:,:)
      double precision, allocatable :: centers2(:,:)

      integer nlmax
      integer i,itype,j

      double precision, allocatable :: centerstmp(:,:,:)
      double precision, allocatable :: boxsize(:)
      integer, allocatable :: irefinebox(:)

      double precision rsc
      integer nbloc,nbctr,nbadd,irefine,ilev,ifirstbox,ilastbox
      integer iii
      integer ibox,nn,nss,ntt
      double precision sizey,sizez

      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision dfac

      nbmax = 100000

      allocate(boxsize(0:nlmax))

      
      allocate(laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax))
      allocate(nchild(nbmax),ichild(8,nbmax))

      allocate(centers(3,nbmax),isrcse(2,nbmax),itargse(2,nbmax))
      allocate(isrc(ns),itarg(nt))

c
c     step 1: find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
      do i=1,ns
        if(src(1,i).lt.xmin) xmin = src(1,i)
        if(src(1,i).gt.xmax) xmax = src(1,i)
        if(src(2,i).lt.ymin) ymin = src(2,i)
        if(src(2,i).gt.ymax) ymax = src(2,i)
        if(src(3,i).lt.zmin) zmin = src(3,i)
        if(src(3,i).gt.zmax) zmax = src(3,i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
      do i=1,nt
        if(targ(1,i).lt.xmin) xmin = targ(1,i)
        if(targ(1,i).gt.xmax) xmax = targ(1,i)
        if(targ(2,i).lt.ymin) ymin = targ(2,i)
        if(targ(2,i).gt.ymax) ymax = targ(2,i)
        if(targ(3,i).lt.zmin) zmin = targ(3,i)
        if(targ(3,i).gt.zmax) zmax = targ(3,i)
      enddo
C$OMP END PARALLEL DO      

      boxsize(0) = (xmax - xmin)
      sizey = (ymax - ymin)
      sizez = (zmax - zmin)
      if(sizey.gt.boxsize(0)) boxsize(0) = sizey
      if(sizez.gt.boxsize(0)) boxsize(0) = sizez

c
c      set tree info for level 0
c
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparent(1) = -1
      nchild(1) = 0
      do i=1,8
        ichild(i,1) = -1
      enddo

      centers(1,1) = (xmin+xmax)/2
      centers(2,1) = (ymin+ymax)/2
      centers(3,1) = (zmin+zmax)/2

      isrcse(1,1) = 1
      isrcse(2,1) = ns
      
      itargse(1,1) = 1
      itargse(2,1) = nt

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        isrc(i) = i
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nt
        itarg(i) = i
      enddo
C$OMP END PARALLEL DO      


      nbctr = 1


      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(irefinebox(nbloc))
c
c          determine which boxes need to be refined
c


        if(ilev.ge.nlmin) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nss,ntt,nn)
          do i=1,nbloc
            irefinebox(i) = 0
            ibox = ifirstbox + i-1
            nss = isrcse(2,ibox)-isrcse(1,ibox)+1
            ntt = itargse(2,ibox)-itargse(1,ibox)+1
            if(idivflag.eq.0) nn = nss
            if(idivflag.eq.1) nn = ntt
            if(idivflag.eq.2) nn = max(ntt,nss)

            if(nn.gt.ndiv) irefinebox(i) = 1
          enddo
C$OMP END PARALLEL DO        

          irefine = maxval(irefinebox(1:nbloc))
          if(ifunif.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)         
            do i=1,nbloc
              irefinebox(i) = irefine
            enddo
C$OMP END PARALLEL DO            
          endif
        else
C$OMP PARALLEL DO DEFAULT(SHARED)        
          do i=1,nbloc
            irefinebox(i) = 1
          enddo
C$OMP END PARALLEL DO         
          irefine = 1
        endif

c
c
c          figure out if current set of boxes is sufficient
c

        nbadd = 0 
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:nbadd)        
        do i=1,nbloc
          if(irefinebox(i).eq.1) nbadd = nbadd+8
        enddo
C$OMP END PARALLEL DO        

        nbtot = nbctr+nbadd

c
c         if current memory is not sufficient reallocate
c
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(3,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(8,nbmax),isrcse2(2,nbmax))
          allocate(itargse2(2,nbmax))

          call tree_copy(nbctr,centers,ilevel,iparent,nchild,
     1            ichild,centers2,ilevel2,iparent2,
     2            nchild2,ichild2)

C$OMP PARALLEL DO DEFAULT(SHARED)     
          do i=1,nbctr
            isrcse2(1,i) = isrcse(1,i)
            isrcse2(2,i) = isrcse(2,i)
            itargse2(1,i) = itargse(1,i)
            itargse2(2,i) = itargse(2,i)
          enddo
C$OMP END PARALLEL DO          


          deallocate(centers,ilevel,iparent,nchild,ichild,
     1        isrcse,itargse)

          nbmax = nbtot
          allocate(centers(3,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(8,nbmax),isrcse(2,nbmax))
          allocate(itargse(2,nbmax))


          call tree_copy(nbctr,centers2,ilevel2,iparent2,
     1            nchild2,ichild2,centers,ilevel,iparent,nchild,ichild)

C$OMP PARALLEL DO DEFAULT(SHARED)
          do i=1,nbctr
            isrcse(1,i) = isrcse2(1,i)
            isrcse(2,i) = isrcse2(2,i)
            itargse(1,i) = itargse2(1,i)
            itargse(2,i) = itargse2(2,i)
          enddo
C$OMP END PARALLEL DO          

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,
     1       isrcse2,itargse2)
        endif


        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          laddr(1,ilev+1) = nbctr+1

          call tree_refine_boxes(irefinebox,nbmax,
     1       ifirstbox,nbloc,centers,boxsize(ilev+1),nbctr,ilev+1,
     2       ilevel,iparent,nchild,ichild)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox) 
          do i=1,nbloc
            ibox = ifirstbox+i-1
            if(irefinebox(i).eq.1) then
              call sort_pts_to_children(ibox,nbmax,centers,ichild,
     1            src,ns,isrc,isrcse)
              call sort_pts_to_children(ibox,nbmax,centers,ichild,
     1            targ,nt,itarg,itargse)
            endif
          enddo
C$OMP END PARALLEL DO          

          laddr(2,ilev+1) = nbctr
        else
          exit
        endif

        deallocate(irefinebox)
      enddo

      nboxes = nbctr
      nlevels = ilev

      if(nlevels.ge.2.and.ifunif.ne.1) then

        nbtot = 16*nboxes
        if(nbtot.gt.nbmax) then
          allocate(centers2(3,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(8,nbmax),isrcse2(2,nbmax))
          allocate(itargse2(2,nbmax))
          call tree_copy(nbctr,centers,ilevel,iparent,nchild,
     1            ichild,centers2,ilevel2,iparent2,
     2            nchild2,ichild2)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)     
          do i=1,nbctr
            isrcse2(1,i) = isrcse(1,i)
            isrcse2(2,i) = isrcse(2,i)
            itargse2(1,i) = itargse(1,i)
            itargse2(2,i) = itargse(2,i)
          enddo
C$OMP END PARALLEL DO          

          deallocate(centers,ilevel,iparent,nchild,ichild,
     1        isrcse,itargse)

          nbmax = nbtot
          allocate(centers(3,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(8,nbmax),isrcse(2,nbmax))
          allocate(itargse(2,nbmax))


          call tree_copy(nbctr,centers2,ilevel2,iparent2,
     1            nchild2,ichild2,centers,ilevel,iparent,nchild,ichild)

C$OMP PARALLEL DO DEFAULT(SHARED)     
          do i=1,nbctr
            isrcse(1,i) = isrcse2(1,i)
            isrcse(2,i) = isrcse2(2,i)
            itargse(1,i) = itargse2(1,i)
            itargse(2,i) = itargse2(2,i)
          enddo
C$OMP END PARALLEL DO          

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,
     1       isrcse2,itargse2)

        endif

        allocate(nnbors(nbmax))
        allocate(nbors(27,nbmax))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,nboxes
          nnbors(i) = 0
          do j=1,27
            nbors(j,i) = -1
          enddo
        enddo
C$OMP END PARALLEL DO        

        call computecoll(nlevels,nboxes,laddr,boxsize,centers,
     1        iparent,nchild,ichild,iper,nnbors,nbors)

        if(nlevels.ge.2.and.ifunif.ne.1) then
          call pts_tree_fix_lr(centers,nlevels,nboxes,boxsize,nbmax,
     1      nlmax,iper,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors)
        endif

      endif
      nboxes8 = nboxes
      ltree = 39*nboxes8 + 2*(nlevels+1) 


      return
      end
c
c
c
c
c

      subroutine pts_tree_build(src,ns,targ,nt,idivflag,ndiv,
     1  nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,iptr,centers,
     2  boxsize)
c
c
c
c----------------------------------------
c  build tree
c
c
c  input parameters:
c    - src: real *8 (3,ns)
c        source locations
c    - targ: real *8 (3,nt) 
c        target locations
c    - idivflag: integer
c        subdivision criterion
c          * divflag = 0 -> subdivide on sources only
c          * idivflag = 1 -> subdivide on targets only
c          * idivflag = 2 -> subdivide on max(sources+targets)
c    - ndiv: integer
c        subdivide if relevant number of particles
c        per box is greater than ndiv
c    - nlmin: integer
c        minimum number of levels of uniform refinement.
c        Note that empty boxes are not pruned along the way
c    - nlmax: integer
c        max number of levels
c    - ifunif: integer
c        flag for creating uniform pruned tree
c        Tree is uniform if ifunif=1 (Currently pruned part
c        under construction)
c    - iper: integer
c        flag for periodic implementations. Currently unused.
c        Feature under construction
c    - nlevels: integer
c        number of levels
c    - nboxes: integer
c        number of boxes
c    - ltree: integer
c
c  output:
c    - itree: integer(ltree)
c        tree info
c    - iptr: integer(8)
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - centers: double precision (3,nboxes)
c        xy coordinates of box centers in the oct tree
c    - boxsize: double precision (0:nlevels)
c        size of box at each of the levels
c

      implicit none
      integer nlevels,nboxes,ns,nt,idivflag,ndiv
      integer *8 iptr(8),ltree
      integer itree(ltree),iper
      integer ifunif,nlmin,nlmax
      double precision centers(3,nboxes),src(3,ns),targ(3,nt)
      integer, allocatable :: irefinebox(:)
      double precision boxsize(0:nlevels)
      integer, allocatable :: isrc(:),itarg(:),isrcse(:,:),itargse(:,:)

      integer i,ilev,irefine,itype,nbmax,npbox,npc,ii
      integer ifirstbox,ilastbox,nbctr,nbloc
      double precision rsc

      double precision ra
      integer j,nboxes0
      integer ibox,nn,nss,ntt

      double precision xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez

c
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + 8*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + 27*nboxes
c
c     step 1: find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
      do i=1,ns
        if(src(1,i).lt.xmin) xmin = src(1,i)
        if(src(1,i).gt.xmax) xmax = src(1,i)
        if(src(2,i).lt.ymin) ymin = src(2,i)
        if(src(2,i).gt.ymax) ymax = src(2,i)
        if(src(3,i).lt.zmin) zmin = src(3,i)
        if(src(3,i).gt.zmax) zmax = src(3,i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
      do i=1,nt
        if(targ(1,i).lt.xmin) xmin = targ(1,i)
        if(targ(1,i).gt.xmax) xmax = targ(1,i)
        if(targ(2,i).lt.ymin) ymin = targ(2,i)
        if(targ(2,i).gt.ymax) ymax = targ(2,i)
        if(targ(3,i).lt.zmin) zmin = targ(3,i)
        if(targ(3,i).gt.zmax) zmax = targ(3,i)
      enddo
C$OMP END PARALLEL DO      

      boxsize(0) = (xmax - xmin)
      sizey = (ymax - ymin)
      sizez = (zmax - zmin)
      if(sizey.gt.boxsize(0)) boxsize(0) = sizey
      if(sizez.gt.boxsize(0)) boxsize(0) = sizez


      allocate(isrc(ns),itarg(nt),isrcse(2,nboxes),itargse(2,nboxes))

c
c      set tree info for level 0
c
      itree(1) = 1
      itree(2) = 1
      itree(iptr(2)) = 0
      itree(iptr(3)) = -1
      itree(iptr(4)) = 0
      do i=1,8
        itree(iptr(5)+i-1) = -1
      enddo

      centers(1,1) = (xmin+xmax)/2
      centers(2,1) = (ymin+ymax)/2
      centers(3,1) = (zmin+zmax)/2

      isrcse(1,1) = 1
      isrcse(2,1) = ns
      itargse(1,1) = 1
      itargse(2,1) = nt

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        isrc(i) = i
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nt
        itarg(i) = i
      enddo
C$OMP END PARALLEL DO      


c
c       Reset nlevels, nboxes
c
      nbctr = 1

      do ilev=0,nlevels-1
        irefine = 0

        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)

        nbloc = ilastbox-ifirstbox+1
        allocate(irefinebox(nbloc))

c
c          determine which boxes need to be refined
c
c       
        if(ilev.ge.nlmin) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nss,ntt,nn)
          do i=1,nbloc
            irefinebox(i) = 0
            ibox = ifirstbox + i-1
            nss = isrcse(2,ibox)-isrcse(1,ibox)+1
            ntt = itargse(2,ibox)-itargse(1,ibox)+1
          
            if(idivflag.eq.0) nn = nss
            if(idivflag.eq.1) nn = ntt
            if(idivflag.eq.2) nn = max(ntt,nss)

            if(nn.gt.ndiv) irefinebox(i) = 1
          enddo
C$OMP END PARALLEL DO        
          irefine = maxval(irefinebox(1:nbloc))

          if(ifunif.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
            do i=1,nbloc
              irefinebox(i) = irefine
            enddo
C$OMP END PARALLEL DO 
          endif
        else
C$OMP PARALLEL DO DEFAULT(SHARED)
          do i=1,nbloc
            irefinebox(i) = 1
          enddo
C$OMP END PARALLEL DO 
          irefine = 1
        endif
        

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

          call tree_refine_boxes(irefinebox,nboxes,
     1       ifirstbox,nbloc,centers,boxsize(ilev+1),nbctr,ilev+1,
     2       itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     3       itree(iptr(5)))

c
c     re sort points in refined boxes
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox)
          do i=1,nbloc
            ibox = ifirstbox+i-1
            if(irefinebox(i).eq.1) then
              call sort_pts_to_children(ibox,nboxes,centers,
     1          itree(iptr(5)),src,ns,isrc,isrcse)
              call sort_pts_to_children(ibox,nboxes,centers,
     1          itree(iptr(5)),targ,nt,itarg,itargse)
            endif
          enddo
C$OMP END PARALLEL DO          

          
          itree(2*ilev+4) = nbctr
        else
          exit
        endif

        deallocate(irefinebox)
      enddo

      nboxes0 = nbctr
      nlevels = ilev


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes0
        itree(iptr(6)+i-1) = 0
        do j=1,27
          itree(iptr(7)+27*(i-1)+j-1) = -1
        enddo
      enddo
C$OMP END PARALLEL DO      


      call computecoll(nlevels,nboxes0,itree(iptr(1)),boxsize,centers,
     1        itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),iper,
     2        itree(iptr(6)),itree(iptr(7)))

      if(nlevels.ge.2.and.ifunif.ne.1) then
         call pts_tree_fix_lr(centers,nlevels,
     1         nboxes0,boxsize,nboxes,nlevels,iper,itree(iptr(1)),
     2         itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     3         itree(iptr(5)),itree(iptr(6)),itree(iptr(7)))

      endif

      return
      end
c
c
c
      subroutine sort_pts_to_children(ibox,nboxes,centers,
     1   ichild,src,ns,isrc,isrcse)
      implicit double precision (a-h,o-z)
      integer nboxes
      double precision centers(3,nboxes),src(3,ns)
      integer ns, isrc(ns),isrcse(2,nboxes)
      integer ichild(8,nboxes)
      integer iss, nsc(8)
      integer, allocatable :: isrctmp(:)

c           Allocate temporary array to figure out which child you
c           belong to
c           which child?  1,2,3,4,5,6,7,8? counter ns1,ns2,ns3,ns4,
c           ns5,ns6,ns7,ns8
c           The box nomenclature is as follows
c           3   4       7  8       <--- Looking down in z direction
c           1   2       5  6
c
c           If the parent box center is at the origin then the centers
c           of box i have co-ordinates
c           1     x<0,y<0,z<0
c           2     x>0,y<0,z<0
c           3     x<0,y>0,z<0
c           4     x>0,y>0,z<0
c           5     x<0,y<0,z>0
c           6     x>0,y<0,z>0
c           7     x<0,y>0,z>0
c           8     x>0,y>0,z>0
c

      i1234 = isrcse(1,ibox)-1
      i5678 = 0
      npts = isrcse(2,ibox)-isrcse(1,ibox)+1
      allocate(isrctmp(npts))

      do iss=isrcse(1,ibox),isrcse(2,ibox)
        if(src(3,isrc(iss))-centers(3,ibox).lt.0) then
          i1234 = i1234+1
          isrc(i1234) = isrc(iss)
        else
          i5678 = i5678 + 1
          isrctmp(i5678) = isrc(iss)
        endif
      enddo

c
c           reorder sources
c
      do i=1,i5678
        isrc(i1234+i) = isrctmp(i)
      enddo


c
c     sort i1234 to i12 and i34
c
      i12 = isrcse(1,ibox)-1
      i34 = 0
      do iss=isrcse(1,ibox),i1234
        if(src(2,isrc(iss))-centers(2,ibox).lt.0) then
          i12 = i12+1
          isrc(i12) = isrc(iss)
        else
          i34 = i34 + 1
          isrctmp(i34) = isrc(iss)
        endif
      enddo
c
c  reorder sources to include 34 in the array
c
      do i=1,i34
        isrc(i12+i) = isrctmp(i)
      enddo

c
c    sort i5678 into i56 and i78
c
      i56 = i1234
      i78 = 0
      do iss=i1234+1,isrcse(2,ibox)
        if(src(2,isrc(iss))-centers(2,ibox).lt.0) then
          i56 = i56+1
          isrc(i56) = isrc(iss)
        else
          i78 = i78 + 1
          isrctmp(i78) = isrc(iss)
        endif
      enddo

c     reorder sources to include 78 in the array      
      do i=1,i78
        isrc(i56+i) = isrctmp(i)
      enddo

      nsc(1) = 0
      nsc(2) = 0
      nsc(3) = 0
      nsc(4) = 0
      nsc(5) = 0
      nsc(6) = 0
      nsc(7) = 0
      nsc(8) = 0
c
c       sort into boxes 1,2
c
      do iss = isrcse(1,ibox),i12
        if(src(1,isrc(iss))-centers(1,ibox).lt.0) then
          isrc(isrcse(1,ibox)+nsc(1)) = isrc(iss)
          nsc(1) = nsc(1)+1
        else
          nsc(2) = nsc(2) + 1
          isrctmp(nsc(2)) = isrc(iss)
        endif
      enddo
c
c      reorder sources so that sources in 2 are at end of array
c
      do i=1,nsc(2)
        isrc(isrcse(1,ibox)+nsc(1)+i-1) = isrctmp(i)
      enddo
c
c    sort intow boxes 3,4 
c 
      do iss=i12+1,i1234
        if(src(1,isrc(iss))-centers(1,ibox).lt.0) then
          isrc(i12+1+nsc(3)) = isrc(iss)
          nsc(3) = nsc(3) + 1
        else
          nsc(4) = nsc(4) + 1
          isrctmp(nsc(4)) = isrc(iss)
        endif
      enddo

      do i=1,nsc(4)
        isrc(i12+nsc(3)+i) = isrctmp(i)
      enddo

c
c   sort into boxes 5,6
c
      do iss=i1234+1,i56
        if(src(1,isrc(iss))-centers(1,ibox).lt.0) then
          isrc(i1234+1+nsc(5)) = isrc(iss)
          nsc(5) = nsc(5) + 1
        else
          nsc(6) = nsc(6) + 1
          isrctmp(nsc(6)) = isrc(iss)
        endif
      enddo

      do i=1,nsc(6)
        isrc(i1234+nsc(5)+i) = isrctmp(i)
      enddo
c
c   End of sorting sources into boxes 5 and 6
c

c
c    sort into boxes 7 and 8
c
      do iss=i56+1,isrcse(2,ibox)
        if(src(1,isrc(iss))-centers(1,ibox).lt.0) then
          isrc(i56+1+nsc(7)) = isrc(iss)
          nsc(7) = nsc(7) + 1
        else
          nsc(8) = nsc(8) + 1
          isrctmp(nsc(8)) = isrc(iss)
        endif
      enddo

      do i=1,nsc(8)
        isrc(i56+nsc(7)+i) = isrctmp(i)
      enddo
c
c   End of sorting sources into boxes 5 and 6


      istart = isrcse(1,ibox)
      do i=1,8
        jbox = ichild(i,ibox)
        isrcse(1,jbox) = istart
        isrcse(2,jbox) = istart + nsc(i)-1
        istart = istart + nsc(i)
      enddo

c
c

      return
      end
c       
c
c-------------------------------------------------------------      
      subroutine pts_tree_fix_lr(centers,nlevels,nboxes,
     1       boxsize,nbmax,nlmax,iper,laddr,ilevel,iparent,nchild,
     2       ichild,nnbors,nbors)
c
c
c       convert an adaptive tree into a level restricted tree
c
      implicit none
      integer nlevels,nboxes,nlmax
      integer nbmax,iper
      double precision centers(3,nbmax),boxsize(0:nlmax)
      integer laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax)
      integer nchild(nbmax),ichild(8,nbmax),nnbors(nbmax)
      integer nbors(27,nbmax)
      integer laddrtail(2,0:nlmax),isum
      integer, allocatable :: iflag(:)

      integer i,j,k,l,ibox,jbox,kbox,ilev,idad,igranddad
      integer nbloc,ict
      double precision xdis,ydis,zdis,distest

      allocate(iflag(nbmax))

c     Initialize flag array
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         iflag(i) = 0
      enddo
C$OMP END PARALLEL DO     



c     Flag boxes that violate level restriction by "1"
c     Violatioin refers to any box that is directly touching
c     a box that is more than one level finer
c
c     Method:
c     1) Carry out upward pass. For each box B, look at
c     the colleagues of B's grandparent
c     2) See if any of those colleagues are childless and in
c     contact with B.
c
c     Note that we only need to get up to level two, as
c     we will not find a violation at level 0 and level 1
c
c     For such boxes, we set iflag(i) = 1
c
      do ilev=nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by two levels are touching
         distest = 1.05d0*(boxsize(ilev-1) + boxsize(ilev-2))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,igranddad,i,jbox)         
C$OMP$ PRIVATE(ict,xdis,ydis,zdis)
         do ibox = laddr(1,ilev),laddr(2,ilev) 
            idad = iparent(ibox)
            igranddad = iparent(idad)
            
c           Loop over colleagues of granddad            
            do i=1,nnbors(igranddad)
               jbox = nbors(i,igranddad)
c              Check if the colleague of grandad
c              is a leaf node. This automatically
c              eliminates the granddad
               if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                   xdis = centers(1,jbox) - centers(1,idad)
                   ydis = centers(2,jbox) - centers(2,idad)
                   zdis = centers(3,jbox) - centers(3,idad)
                   ict = 0
                   if(abs(xdis).le.distest) ict = ict + 1
                   if(abs(ydis).le.distest) ict = ict + 1
                   if(abs(zdis).le.distest) ict = ict + 1
                   if(ict.eq.3) then
                      iflag(jbox) = 1
                   endif
               endif
c              End of checking criteria for the colleague of
c              granddad
            enddo
c           End of looping over colleagues of
c           granddad
         enddo
c        End of looping over boxes at ilev         
C$OMP END PARALLEL DO
      enddo
c     End of looping over levels and flagging boxes


c     Find all boxes that need to be given a flag+
c     A flag+ box will be denoted by setting iflag(box) = 2
c     This refers to any box that is not already flagged and
c     is bigger than and is contacting a flagged box
c     or another box that has already been given a flag +.
c     It is found by performing an upward pass and looking
c     at the flagged box's parents colleagues and a flag+
c     box's parents colleagues and seeing if they are
c     childless and present the case where a bigger box 
c     is contacting a flagged or flag+ box.

      do ilev = nlevels,1,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,i,jbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if(iflag(ibox).eq.1.or.iflag(ibox).eq.2) then
               idad = iparent(ibox)
c              Loop over dad's colleagues               
               do i=1,nnbors(idad)
                  jbox = nbors(i,idad)
c                 Check if the colleague of dad
c                 is a leaf node. This automatically
c                 eliminates the dad
                  if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                     xdis = centers(1,jbox) - centers(1,ibox)
                     ydis = centers(2,jbox) - centers(2,ibox)
                     zdis = centers(3,jbox) - centers(3,ibox)
                     ict = 0
                     if(abs(xdis).le.distest) ict = ict + 1
                     if(abs(ydis).le.distest) ict = ict + 1
                     if(abs(zdis).le.distest) ict = ict + 1
                     if(ict.eq.3) then
                        iflag(jbox) = 2
                     endif
                  endif
c                 End of checking criteria for the colleague of
c                dad
               enddo
c              End of looping over dad's colleagues               
            endif
c           End of checking if current box is relevant for
c           flagging flag+ boxes
         enddo
c        End of looping over boxes at ilev        
C$OMP END PARALLEL DO 
      enddo
c     End of looping over levels

c     Subdivide all flag and flag+ boxes. Flag all the children
c     of flagged boxes as flag++. Flag++ boxes are denoted
c     by setting iflag(box) = 3. The flag++ boxes need 
c     to be checked later to see which of them need further
c     refinement. While creating new boxes, we will
c     need to update all the tree structures as well.
c     Note that all the flagged boxes live between
c     levels 1 and nlevels - 2. We process the boxes via a
c     downward pass. We first determine the number of boxes
c     that are going to be subdivided at each level and 
c     everything else accordingly
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo

 
      do ilev = 1,nlevels-2
c        First subdivide all the flag and flag+
c        boxes with boxno nboxes+1, nboxes+ 2
c        and so on. In the second step, we reorganize
c        all the structures again to bring it back
c        in the standard format

         laddrtail(1,ilev+1) = nboxes+1


         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call tree_refine_boxes_flag(iflag,nbmax,laddr(1,ilev),
     1    nbloc,centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,
     2    nchild,ichild)


         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call pts_tree_reorg(nboxes,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,
     2          iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,27
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO      
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)

c     Processing of flag and flag+ boxes is done
c     Start processing flag++ boxes. We will use a similar
c     strategy as before. We keep checking the flag++
c     boxes that require subdivision if they still
c     violate the level restriction criterion, create
c     the new boxes, append them to the end of the list to begin
c     with and in the end reorganize the tree structure.
c     We shall accomplish this via a downward pass
c     as new boxes that get added in the downward pass
c     will also be processed simultaneously.
c     We shall additionally also need to keep on updating
c     the colleague information as we proceed in the 
c     downward pass

c     Reset the flags array to remove all the flag and flag+
c     cases. This is to ensure reusability of the subdivide
c     _flag routine to handle the flag++ case

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox)
      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
C$OMP END PARALLEL DO      
 
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo


      do ilev = 2,nlevels-2

c     Step 1: Determine which of the flag++ boxes need
c     further division. In the even a flag++ box needs
c     further subdivision then flag the box with iflag(box) = 1
c     This will again ensure that the subdivide_flag routine
c     will take care of handling the flag++ case
         call updateflags(ilev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

         call updateflags(ilev,nboxes,nlevels,laddrtail,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)
        
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call tree_refine_boxes_flag(iflag,nbmax,laddr(1,ilev),
     1    nbloc,centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,
     2    nchild,ichild)

         nbloc = laddrtail(2,ilev)-laddrtail(1,ilev)+1
         call tree_refine_boxes_flag(iflag,nbmax,laddrtail(1,ilev),
     1    nbloc,centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,
     2    nchild,ichild)

         laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,idad,jbox,j,kbox)
          do ibox = laddrtail(1,ilev+1),laddrtail(2,ilev+1)
            nnbors(ibox) = 0
c           Find the parent of the current box         
            idad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out colleagues
            do i=1,nnbors(idad)
                jbox = nbors(i,idad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev+1)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev+1)).and.
     2                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     3                   1.05*boxsize(ilev+1))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

c     Reorganize tree once again and we are all done      
      call pts_tree_reorg(nboxes,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,
     2          iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,27
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO    

      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)
      

      return
      end
      

c-------------------------------------------------------------      

      subroutine pts_tree_reorg(nboxes,centers,nlevels,laddr,
     1     laddrtail,ilevel,iparent,nchild,ichild,iflag)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    centers        in/out: double precision(3,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(4,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    iflag        in/out: integer(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer nboxes,nlevels
      double precision centers(3,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer laddrtail(2,0:nlevels)
      integer ilevel(nboxes)
      integer iparent(nboxes)
      integer nchild(nboxes)
      integer ichild(8,nboxes)
      integer iflag(nboxes)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:),tiflag(:)
      integer, allocatable :: iboxtocurbox(:),ilevptr(:),ilevptr2(:)

      double precision, allocatable :: tcenters(:,:)



c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox,idim,nblev

      allocate(tilevel(nboxes),tiparent(nboxes),tnchild(nboxes))
      allocate(tichild(8,nboxes),tiflag(nboxes),iboxtocurbox(nboxes))
      allocate(tcenters(3,nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo
      call tree_copy(nboxes,centers,ilevel,iparent,nchild,
     1            ichild,tcenters,tilevel,tiparent,
     2            tnchild,tichild)

      do ibox=1,nboxes
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
         enddo
      enddo

      allocate(ilevptr(nlevels+1),ilevptr2(nlevels))

      ilevptr(2) = laddr(1,2)

      do ilev=2,nlevels
        nblev = laddr(2,ilev)-laddr(1,ilev)+1
        ilevptr2(ilev) = ilevptr(ilev) + nblev
        nblev = laddrtail(2,ilev)-laddrtail(1,ilev)+1
        ilevptr(ilev+1) = ilevptr2(ilev) + nblev
      enddo

      curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            centers(3,curbox) = tcenters(3,ibox)
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            centers(3,curbox) = tcenters(3,ibox)
            nchild(curbox) = tnchild(ibox)
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         do i=1,8
            if(tichild(i,ibox).eq.-1) ichild(i,iboxtocurbox(ibox)) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo

      return
      end
c
c
c
c
c
c
      subroutine pts_tree_sort(n,xys,itree,ltree,nboxes,nlevels,
     1   iptr,centers,ixy,ixyse)
      implicit double precision (a-h,o-z)
      integer *8 iptr(8),ltree
      integer n,nboxes,nlevels,itree(ltree)
      integer ixy(n),ixyse(2,nboxes)
      double precision xys(3,n),centers(3,nboxes)

      do i=1,n
        ixy(i) = i
      enddo

      ixyse(1,1) = 1
      ixyse(2,1) = n

      do ilev = 0,nlevels-1
        do ibox=itree(2*ilev+1),itree(2*ilev+2)
          if(itree(iptr(4)+ibox-1).gt.0) then
            call sort_pts_to_children(ibox,nboxes,centers,
     1          itree(iptr(5)),xys,n,ixy,ixyse)
          endif
        enddo
      enddo

      return
      end

