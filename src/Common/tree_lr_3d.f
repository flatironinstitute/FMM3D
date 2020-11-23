c
c    generate level restricted oct tree refined 
c    to finest level. Features include 
c
c     radsrc parameter (preventing sources with large radius 
c                         from being pushed to the finest level)
c     radexp parameter (preventing exp centers with large radius 
c                         from being pushed to the finest level)
c     ndiv parameter (usual refinement parameter extending tree
c                         depth until all leaf nodes hanve fewer than
c                         ndiv particles)
c     idivflag       (allows use of ndiv parameter for sources OR
c                         targets OR sources+targets OR sources+targets
c                         + expansion centers)
c   
c     isep           (separation parameter: isep = 1 or 2
c                     based on whether to include neighbors
c                     and self in list 1 or two neighbors
c                     and self in list 2)
c     
c     extraction of plane wave lists
c
c     mhung is the max number of hung chunks in a given box
c
c     This tree code MUST be used along with its
c     memomry code mkptreemem located at the end of this file.
c     The memory code precomputes the number of levels, the number
c     of boxes, the max number of hung sources and the length
c     of the tree. 
c
      subroutine mklraptree(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,isep,
     $                   mhung,mnbors,mnlist1,mnlist2,mnlist3,
     $                   mnlist4,nlevels,nboxes,
     $                   centers,boxsize,itree,ltree,
     $                   ipointer)
      implicit none
      integer ier,ns,nt,nexpc,idivflag,ndiv,isep,mhung
      integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
      integer nlevels,lcenters
      integer *8 ltree
      integer nlmax,nbmax,nboxes, nlevtmp,nbtmp, mhungtmp
      integer itree(ltree)
      integer i,j
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: ilevel(:)
      integer, allocatable :: iparenttemp(:)
      integer, allocatable :: nchild(:)
      integer, allocatable :: ichildtemp(:,:)
      integer, allocatable :: nnbors(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:)
      integer, allocatable :: ihslasttemp(:)
      integer, allocatable :: isfirsttemp(:)
      integer, allocatable :: islasttemp(:)
      integer, allocatable :: itfirsttemp(:)
      integer, allocatable :: itlasttemp(:)
      integer, allocatable :: ihefirsttemp(:)
      integer, allocatable :: ihelasttemp(:)
      integer, allocatable :: iefirsttemp(:)
      integer, allocatable :: ielasttemp(:)
      integer, allocatable :: nhungsrc(:)
      integer, allocatable :: nhunglistsrc(:)
      integer, allocatable :: nhungexp(:)
      integer, allocatable :: ihunglistsrc(:,:)

      integer, allocatable :: nlist1(:)
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: nlist3(:)
      integer, allocatable :: list3(:,:)
      integer, allocatable :: nlist4(:)
      integer, allocatable :: list4(:,:)


      integer *8 ipointer(32)
      double precision boxsize(0:nlevels)
      double precision src(3,ns),radsrc(ns)
      double precision trg(3,nt)
      double precision centers(3,nboxes)
      double precision expc(3,nexpc)
      double precision radexp(nexpc)
c
      double precision xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez
      integer ictr,ih,irefine,is,ie
      integer nss,nee

c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     isep          separation parameter for determining 
c                   nearest neighbors
c
c     ltree         length of ltree
c     nlevels       number of levels determined using the memory code
c     nboxes        number of boxes determined using mem code
c                   refer to mkptreemem.
c
c
c     OUTPUT:
c     centers       array of box centers
c     boxsize       box dimensions at all levels
c
c     itree         tree array
c     
c     iladdr = 1
c     itree(iladdr) <-> laddr
c                   indexing array providing access to boxes at
c                   each level. 
c                   the first box on level i is laddr(1,i)
c                   the last box on level i is laddr(2,i)
c
c     iparent = iladdr + 2*(nlevels+1)
c     itree(iparent) <-> parent
c                   parent of box (set to -1 for root node)
c
c     inchild = iparent + nboxes
c     itree(inchild) <-> nchild
c                   number of children for each box 
c
c     ichild = inchild + nboxes
c     itree(ichild) <-> child
c                   list of children for each box
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (8,nboxes)
c
c     isource = ichild + 8*nboxes
c     itree(isource) <-> isource
c                   tree-ordered array of sources
c
c     itarget = isource + ns
c     itree(itarget) <-> itarget
c                   tree-ordered array of targets
c
c     iexpc = itarget + nt
c     itree(iexpc) <-> iexpc
c                      tree-ordered array of expansion centers
c
c     ihsfirst = itarget + nt
c     itree(ihsfirst) <-> ihsfirst
c                   ihsfirst(j) = location in isource of first hung 
c                                source for box j 
c
c     ihslast = ihsfirst + nboxes
c     itree(ihslast) <-> ihslast
c                   ihslast(j)  = location in isource of last hung
c                                source for box j 
c
c     isfirst = ihlast + nboxes
c     itree(isfirst) <-> isfirst
c                   isfirst(j) = location in isource of first 
c                                source for box j 
c
c     islast = isfirst + nboxes
c     itree(islast) <-> islast
c                   islast(j)  = location in isource of last 
c                                source for box j 
c
c     itfirst = islast + nboxes
c     itree(itfirst) <-> itfirst
c                   itfirst(j) = location in itarget of first 
c                                target for box j 
c
c     itlast = itfirst + nboxes
c     itree(itlast) <-> itlast
c                   itlast(j)  = location in itarget of last 
c                                target for box j
c
c     ihefirst = itlast + nboxes
c     itree(ihefirst) <-> ihefirst
c                   ihefirst(j) = location in iexpc of first hung 
c                                 expansion center for box j 
c
c     ihelast = ihefirst + nboxes
c     itree(ihelast) <-> ihelast
c                   ihelast(j)  = location in iexpc of last hung
c                                expansion center for box j 
c
c     iefirst = ihelast + nboxes
c     itree(iefirst) <-> iefirst
c                    iefirst(j) = location in iexpc of first target
c                                 in box j
c
c     ielast = ielast + nboxes
c     itree(ielast) <-> ielast
c                       ielast(j) = location in iexpc of last target
c                                   in box j
c
c     innbor = ielast + nboxes
c     itree(innbor) <-> nnbor
c                   number of neighbors
c
c     inbors = innbor + nboxes
c     itree(inbors) <-> nbors
c                   only first nnbor entries are set 
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (27,nboxes) if isep = 1
c                   and (125, nboxes) if isep = 2
c                   mnnbors = 27/125 if isep=1/2
c
c     inlist1 = inbors + mnnbors*nboxes
c     itree(inlist1) <-> nlist1
c                        nlist1(i) = number of boxes in list 1
c                        of box i. The list 1 of a box i, Ui, is
c                        the set of boxes which touch
c                        box i. In a level restricted tree,
c                        the maximum number of boxes in list
c                        1 of a box can be 13 (mnlist1 = 13)
c
c 
c     ilist1 = inlist1 + nboxes
c     itree(ilist1) <-> list1
c                       list1(j,i) is the id of the jth box
c                       in list 1 of box i
c
c     inlist2 = ilist1 + mnlist1*nboxes
c     itree(inlist2) <-> nlist2
c                       nlist2(i) = number of boxes in list 2 of
c                       box i. The list 2 of a box i, Vi, is
c                       the set of boxes which are descendants
c                       of the colleagues of the parent of box i
c                       which are well separated from box i
c                       at the scale of box i. In a level
c                       restricted tree, the maximum number of boxes
c                       in list 2 of a box can be 27 (mnlist2=27)
c
c     ilist2 = inlist2 + nboxes
c     itree(ilist2) <-> list2
c                       list2(j,i) is the id fo the jth box in 
c                       list 2 of box i
c
c     inlist3 = ilist2 + mnlist2*nboxes
c     itree(inlist3) <-> nlist3
c                        nlist3(i) = number of boxes in list3 of
c                        box i. The list 3 of a box i, Wi, is
c                        the set of boxes which are descendants
c                        of the colleagues of box i which
c                        are not in list 1 of box i. In a level
c                        restricted tree, the maximum number of
c                        boxes in list 3 of a box can be 
c                        20 (mnlist3 = 20). Note
c                        that list 3 of a box is non empty
c                        iff the box is childless
c
c     ilist3 = inlist3 + nboxes
c     itree(ilist3) <-> list3(j,i) is the id of the jth box in
c                       list 3 of box i
c
c     inlist4 = ilist3 + mnlist3*nboxes
c     itree(inlist4) <-> nlist4
c                        nlist4(i) = number of boxes in list4
c                        of box i. The list 4 of a box i, Xi,
c                        is dual to Wi. j in Xi if i in Wj, that
c                        is box j is in list 4 of box i if i is
c                        in list 3 of box j. In a level restricted
c                        tree the max number of boxes in list 4
c                        of a box is 5.
c       
c     ilist4 = inlist4 + nboxes
c     itree(ilist4) <-> list4
c                       list4(j,i) is the id of the jth box
c                       in list 4 of box i
c
c     inhungsrc = ilist4 + mnlist4*nboxes
c     itree(inhungsrc) <-> inhungsrc
c                inhungsrc(j)  = number of hung chunks in box j
c
c     inhungexp = inhungsrc + nboxes
c     itree(inhungexp) <-> inhungexp
c                inhungexp(j) = number of hung expansion centers
c                               in box j
c
c     nhunglistsrc = inhungexp + nboxes
c     itree(inhunglistsrc) <-> inhunglistsrc
c           inhunglistsrc(j) = Total number of hung sources
c                              relevant for box j
c
c     ihunglistsrc = inhungexp + nboxes
c     itree(ihunglistsrc) <-> ihunglistsrc
c                   ihunglistsrc(m,j) = src id of  mth hung src
c                   relevant to box j. A hung src is relevant to
c                   a box if the box is a descendant of the
c                   neighbor of the box in which the src is hung
c
c     ltree = ihunglistsrc + mhung*nboxes
c
c     ipointers is the collection of pointers
c     ipointer(1) = iladdr
c     ipointer(2) = iparent = ipointer(1) + 2*(nlevels+1)
c     ipointer(3) = inchild = ipointer(2) + nboxes
c     ipointer(4) = ichild = ipointer(3) + nboxes
c     ipointer(5) = isource = ipointer(4) + 8*nboxes
c     ipointer(6) = itarget = ipointer(5) + ns
c     ipointer(7) = iexpc = ipointer(6) + nt
c     ipointer(8) = ihsfirst = ipointer(7) + nexpc
c     ipointer(9) = ihslast = ipointer(8) + nboxes
c     ipointer(10) = isfirst = ipointer(9) + nboxes
c     ipointer(11) = islast = ipointer(10) + nboxes
c     ipointer(12) = itfirst = ipointer(11) + nboxes
c     ipointer(13) = itlast = ipointer(12) + nboxes
c     ipointer(14) = ihefirst = ipointer(13) + nboxes
c     ipointer(15) = ihelast = ipointer(14) + nboxes
c     ipointer(16) = iefirst = ipointer(15) + nboxes
c     ipointer(17) = ielast = ipointer(16) + nboxes
c     ipointer(18) = innbors = ipointer(17) + nboxes
c     ipointer(19) = inbors = ipointer(18) + nboxes
c     ipointer(20) = inlist1 = ipointer(19) + mnbors*nboxes
c     ipointer(21) = ilist1 = ipointer(20) + nboxes
c     ipointer(22) = inlist2 = ipointer(21) + mnlist1*nboxes
c     ipointer(23) = ilist2 = ipointer(22) + nboxes
c     ipointer(24) = inlist3 = ipointer(23) + mnlist2*nboxes
c     ipointer(25) = ilist3 = ipointer(24) + nboxes
c     ipointer(26) = inlist4 = ipointer(25) + mnlist3*nboxes
c     ipointer(27) = ilist4 = ipointer(26) + nboxes
c     ipointer(28) = nhungsrc = ipointer(27) + mnlist4*nboxes
c     ipointer(29) = nhungexp = ipointer(28) + nboxes
c     ipointer(30) = nhunglistsrc = ipointer(29) + nboxes
c     ipointer(31) = hunglistsrc = ipointer(30) + nboxes
c     ipointer(32) = lentree = ipointer(31) + mhung*nboxes
c
c     mhung is the max number of hung sources 
c
c---------------------------------------------------------------------
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c
      lcenters = nboxes

      ier = 0

      allocate(laddr(2,0:nlevels))
      allocate(ilevel(nboxes))
      allocate(iparenttemp(nboxes))
      allocate(nchild(nboxes))
      allocate(ichildtemp(8,nboxes))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nboxes))
      allocate(ihslasttemp(nboxes))
      allocate(isfirsttemp(nboxes))
      allocate(islasttemp(nboxes))
      allocate(itfirsttemp(nboxes))
      allocate(itlasttemp(nboxes))
      allocate(ihefirsttemp(nboxes))
      allocate(ihelasttemp(nboxes))
      allocate(iefirsttemp(nboxes))
      allocate(ielasttemp(nboxes))
      allocate(nhungsrc(nboxes))
      allocate(nhungexp(nboxes))
      allocate(nnbors(nboxes))
      allocate(nbors(mnbors,nboxes))

      if(isep.ne.1.and.isep.ne.2) then
         ier = 4
cc         call prinf('Error tree not allocated, isep.ne.1,2*',ier,0)
           print *, "Error tree not allocated, isep.ne.1,2"
         return
      endif

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
        if(expc(3,i) .lt. zmin) zmin=expc(3,i)
        if(expc(3,i) .gt. zmax) zmax=expc(3,i)
      enddo
C$OMP END PARALLEL DO   

      boxsize(0)=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
      if(sizez .gt. boxsize(0)) boxsize(0)=sizez
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      centers(3,1)=(zmin+zmax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
      nhungsrc(1) = 0
      nhungexp(1) = 0
c
c     count number of hung sources
c     and hang up "big" sources
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif

c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
      do i = 1,nt
         itargettemp(i) = i
      enddo
C$OMP END PARALLEL DO      
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1
      ichildtemp(5,1) = -1
      ichildtemp(6,1) = -1
      ichildtemp(7,1) = -1
      ichildtemp(8,1) = -1
c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1

c     Reset nlevels, nboxes
      nbmax = nboxes
      nlmax = nlevels
      nlevels = 0
      nboxes = 1


      do i = 1,nlmax
         if (irefine.eq.1) then
            call subdivide_adap(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)
            if(ier.ne.0) return  
         else
            exit
         endif
      enddo
c     Set up computation of list1 and list2     

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo

      call computecoll(nlevels,nboxes,laddr,boxsize,centers,
     1     iparenttemp,nchild,ichildtemp,mnbors,nnbors,nbors)

      if(nlevels.ge.2) then
      call d3hpfixtree(ier,src,ns,radsrc,trg,nt,expc,nexpc,radexp,
     1     nlevels,nboxes,centers,boxsize,nbmax,laddr,ilevel,
     2     iparenttemp,nchild,ichildtemp,mnbors,nnbors,nbors,
     3     isourcetemp,itargettemp,iexpctemp,ihsfirsttemp,ihslasttemp,
     4     isfirsttemp,islasttemp,itfirsttemp,itlasttemp,ihefirsttemp,
     5     ihelasttemp,iefirsttemp,ielasttemp,nhungsrc,nhungexp)
      if(ier.ne.0) return
      endif

 
      
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO
      
      call computecollisep(nlevels,nboxes,laddr,boxsize,centers,
     1     iparenttemp,nchild,ichildtemp,isep,mnbors,nnbors,nbors)

      allocate(nlist1(nboxes),list1(mnlist1,nboxes))
      allocate(nlist2(nboxes),list2(mnlist2,nboxes))
      allocate(nlist3(nboxes),list3(mnlist3,nboxes))
      allocate(nlist4(nboxes),list4(mnlist4,nboxes))

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
         do j=1,mnlist1
            list1(j,i) = -1
         enddo
         do j=1,mnlist2
            list2(j,i) = -1
         enddo
         do j=1,mnlist3
            list3(j,i) = -1
         enddo
         do j=1,mnlist4
            list4(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO      

      call computelists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,isep,nnbors,mnbors,nbors,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)

      allocate(ihunglistsrc(mhung,nboxes))
      allocate(nhunglistsrc(nboxes))

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j)
      do i=1,nboxes
         nhunglistsrc(i) = 0
         do j=1,mhung
            ihunglistsrc(j,i) = -1
         enddo 
      enddo
C$OMP END PARALLEL DO      

      call computehunglist(mhung,nlevels,nboxes,laddr,ns,
     1     isourcetemp,iparenttemp,ihsfirsttemp,ihslasttemp,ilevel,
     2     nnbors,mnbors,
     2     nbors,nlist1,mnlist1,list1,
     3     nhungsrc,nhunglistsrc,ihunglistsrc)


      ipointer(1) = 1
      ipointer(2) = ipointer(1)+2*(nlevels+1)
      ipointer(3) = ipointer(2) + nboxes
      ipointer(4) = ipointer(3) + nboxes
      ipointer(5) = ipointer(4) + 8*nboxes
      ipointer(6) = ipointer(5) + ns
      ipointer(7) = ipointer(6) + nt
      ipointer(8) = ipointer(7) + nexpc
      ipointer(9) = ipointer(8) + nboxes
      ipointer(10) = ipointer(9) + nboxes
      ipointer(11) = ipointer(10) + nboxes
      ipointer(12) = ipointer(11) + nboxes
      ipointer(13) = ipointer(12) + nboxes
      ipointer(14) = ipointer(13) + nboxes
      ipointer(15) = ipointer(14) + nboxes
      ipointer(16) = ipointer(15) + nboxes
      ipointer(17) = ipointer(16) + nboxes
      ipointer(18) = ipointer(17) + nboxes
      ipointer(19) = ipointer(18) + nboxes
      ipointer(20) = ipointer(19) + mnbors*nboxes
      ipointer(21) = ipointer(20) + nboxes
      ipointer(22) = ipointer(21) + mnlist1*nboxes
      ipointer(23) = ipointer(22) + nboxes
      ipointer(24) = ipointer(23) + mnlist2*nboxes
      ipointer(25) = ipointer(24) + nboxes
      ipointer(26) = ipointer(25) + mnlist3*nboxes
      ipointer(27) = ipointer(26) + nboxes
      ipointer(28) = ipointer(27) + mnlist4*nboxes
      ipointer(29) = ipointer(28) + nboxes
      ipointer(30) = ipointer(29) + nboxes
      ipointer(31) = ipointer(30) + nboxes
      ipointer(32) = ipointer(31) + mhung*nboxes


c     Move the output to itree
c     Store iladdr - first and last box at level ilev
      do i=0,nlevels
         itree(2*i+1) = laddr(1,i)
         itree(2*i+2) = laddr(2,i)
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
         itree(ipointer(5)+i-1) = isourcetemp(i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nt
         itree(ipointer(6)+i-1) = itargettemp(i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nexpc
         itree(ipointer(7)+i-1) = iexpctemp(i)
      enddo
C$OMP END PARALLEL DO      


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i = 1,nboxes
          itree(ipointer(2)+i-1) = iparenttemp(i)
          itree(ipointer(3)+i-1) = nchild(i)
          do j=1,8
             itree(ipointer(4)+8*(i-1)+j-1) = ichildtemp(j,i)
          enddo
          itree(ipointer(8)+i-1) = ihsfirsttemp(i)
          itree(ipointer(9)+i-1) = ihslasttemp(i)
          itree(ipointer(10)+i-1) = isfirsttemp(i)
          itree(ipointer(11)+i-1) = islasttemp(i)
          itree(ipointer(12)+i-1) = itfirsttemp(i)
          itree(ipointer(13)+i-1) = itlasttemp(i)
          itree(ipointer(14)+i-1) = ihefirsttemp(i)
          itree(ipointer(15)+i-1) = ihelasttemp(i)
          itree(ipointer(16)+i-1) = iefirsttemp(i)
          itree(ipointer(17)+i-1) = ielasttemp(i)
          itree(ipointer(18)+i-1) = nnbors(i)
          do j=1,mnbors
             itree(ipointer(19)+(i-1)*mnbors+j-1) = nbors(j,i)
          enddo
          itree(ipointer(20)+i-1) = nlist1(i)
          do j=1,mnlist1
             itree(ipointer(21)+(i-1)*mnlist1+j-1) = list1(j,i)
          enddo
          itree(ipointer(22)+i-1) = nlist2(i)
          do j=1,mnlist2
             itree(ipointer(23)+(i-1)*mnlist2+j-1) = list2(j,i)
          enddo
          itree(ipointer(24)+i-1) = nlist3(i)
          do j=1,mnlist3
             itree(ipointer(25)+(i-1)*mnlist3+j-1) = list3(j,i)
          enddo
          itree(ipointer(26)+i-1) = nlist4(i)
          do j=1,mnlist4
             itree(ipointer(27)+(i-1)*mnlist4+j-1) = list4(j,i)
          enddo
          itree(ipointer(28)+i-1) = nhungsrc(i)
          itree(ipointer(29)+i-1) = nhungexp(i)
          itree(ipointer(30)+i-1) = nhunglistsrc(i)
          do j=1,mhung
             itree(ipointer(31)+(i-1)*mhung+j-1) = ihunglistsrc(j,i)
          enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
c--------------------------------------------------------------
c
      subroutine subdivide_adap(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,irefine)
      implicit none
      integer ier,ns,nt,nexpc,idivflag,ndiv,mhung
      integer nlevels,nboxes,lcenters,nbmax,nlmax
      integer irefine
      double precision src(3,ns),radsrc(ns)
      double precision trg(3,nt)
      double precision expc(3,nexpc),radexp(nexpc)
      double precision centers(3,lcenters)
      double precision boxsize(0:nlmax)
      integer laddr(2,0:nlmax)
      integer ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax)
      integer ichild(8,nbmax)
      integer isource(ns)
      integer itarget(nt)
      integer iexpc(nexpc)
      integer ihsfirst(nbmax)
      integer ihslast(nbmax)
      integer isfirst(nbmax)
      integer islast(nbmax)
      integer itfirst(nbmax)
      integer itlast(nbmax)
      integer ihefirst(nbmax)
      integer ihelast(nbmax)
      integer iefirst(nbmax)
      integer ielast(nbmax)
      integer nhungsrc(nbmax)
      integer nhungexp(nbmax)
c     Temporary variables
      integer isrctmp(ns),itargtmp(nt),iexpctmp(nexpc)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee
      integer jj,irefinebox,ntt
      integer i56, i78, i1234, i5678
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer iss,itt,ie
      integer nsc(8),ntc(8),nh(8),nexpcc(8),nhc(8)
c
c     for every box at level nlevels,
c     sort into children, updating various arrays 
c     perhaps just build tree here paren/child/particle sorting...
c     lists in second call ???
c     
c     allocate temp array for isourcetemp2 itargtemp2
c     after all done, write back to isourcetemp, itargtemp
c     this is O(N) * nlevels work for rewriting.
c     can be fancier I suppose.
c     
      irefine = 0
      ifirstbox = laddr(1,nlevels)
      ilastbox =  laddr(2,nlevels)
c
      nbfirst = nboxes+1
      boxsize(nlevels+1) = boxsize(0)/2.0d0**(nlevels+1)

      do ibox = ifirstbox,ilastbox
c        Determine if current box needs to be subdivided
c        if current box needs to be subdivided then set
c        irefinebox = 1
         nss = islast(ibox) - isfirst(ibox) + 1
         ntt = itlast(ibox) - itfirst(ibox) + 1
         nee = ielast(ibox) - iefirst(ibox) + 1
         irefinebox = 0
         if((idivflag.eq.0).and.(nss.gt.ndiv)) irefinebox=1
         if((idivflag.eq.1).and.(ntt.gt.ndiv)) irefinebox=1
         if((idivflag.eq.2).and.(nss.gt.ndiv.or.ntt.gt.ndiv)) 
     1        irefinebox=1
         if((idivflag.eq.3).and.(nss.gt.ndiv.or.
     1       ntt.gt.ndiv.or.nee.gt.ndiv)) irefinebox=1
         if(irefinebox.eq.1) then
c            Based on the subdivision criterion, the current
c            box needs to be divided
      
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

            i1234 = isfirst(ibox)-1
            i5678 = 0
            do iss = isfirst(ibox),islast(ibox)
               if(src(3,isource(iss)) - centers(3,ibox).lt.0) then
                  i1234 = i1234+1
                  isource(i1234) = isource(iss)
               else
                  i5678 = i5678 + 1
                  isrctmp(i5678) = isource(iss)
               endif
            enddo
c           Note at the end of the loop, i1234 is where the particles
c           in part 1234 of the box end

c           Reorder sources to include sources in 5678 in the array
            do i=1,i5678
               isource(i1234+i) = isrctmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = isfirst(ibox)-1
            i34 = 0
            do iss = isfirst(ibox),i1234
               if(src(2,isource(iss))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  isource(i12) = isource(iss)
               else
                  i34 = i34 + 1
                  isrctmp(i34) = isource(iss)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder sources to include 34 in the array
            do i=1,i34
               isource(i12+i) = isrctmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do iss=i1234+1,islast(ibox)
               if(src(2,isource(iss))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  isource(i56) = isource(iss)
               else
                  i78 = i78 + 1
                  isrctmp(i78) = isource(iss)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               isource(i56+i) = isrctmp(i)
            enddo
c           End of reordering i5678         

            nsc(1) = 0
            nsc(2) = 0
            nsc(3) = 0
            nsc(4) = 0
            nsc(5) = 0
            nsc(6) = 0
            nsc(7) = 0
            nsc(8) = 0
c           Sort into boxes 1 and 2
            do iss = isfirst(ibox),i12
               if(src(1,isource(iss))-centers(1,ibox).lt.0) then
                  isource(isfirst(ibox)+nsc(1)) = isource(iss)
                  nsc(1) = nsc(1) + 1
               else
                  nsc(2) = nsc(2) + 1
                  isrctmp(nsc(2)) = isource(iss)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of this part of the array
            do i=1,nsc(2)
               isource(isfirst(ibox)+nsc(1)+i-1) = isrctmp(i)
            enddo

c           Sort into boxes 3 and 4
            do iss = i12+1, i1234
               if(src(1,isource(iss))-centers(1,ibox).lt.0) then
                  isource(i12+1+nsc(3)) = isource(iss)
                  nsc(3) = nsc(3) + 1
                else
                   nsc(4) = nsc(4)+1
                   isrctmp(nsc(4)) = isource(iss)
                endif
            enddo
c           Reorder sources so that sources in 4 are at the
c           end of this part of the array
            do i=1,nsc(4)
               isource(i12+nsc(3)+i) = isrctmp(i)
            enddo

c           Sort into boxes 5 and 6
            do iss = i1234+1,i56
               if(src(1,isource(iss))-centers(1,ibox).lt.0) then
                  isource(i1234+1+nsc(5)) = isource(iss)
                  nsc(5) = nsc(5) + 1
               else
                  nsc(6) = nsc(6) + 1
                  isrctmp(nsc(6)) = isource(iss)
               endif
            enddo
c           Reorder sources so that sources in 6 are at the
c           end of this part of the array
            do i=1,nsc(6)
               isource(i1234+nsc(5)+i) = isrctmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do iss=i56+1,islast(ibox)
               if(src(1,isource(iss))-centers(1,ibox).lt.0) then
                  isource(i56+1+nsc(7)) = isource(iss)
                  nsc(7) = nsc(7) + 1
               else
                  nsc(8) = nsc(8) + 1
                  isrctmp(nsc(8)) = isource(iss)
               endif
            enddo
c           Reorder sources so that sources in 8 are at the
c           end of the array
            do i=1,nsc(8)
               isource(i56+nsc(7)+i) = isrctmp(i)
            enddo

            istart = isfirst(ibox)-1
            do j=1,8
c           check hung -> counter nh1,nh2,nh3,nh4,nh5,nh6,nh7,nh8
               ii = 0
               nh(j) = 0
               do i=1,nsc(j)
                  if(radsrc(isource(istart+i)).gt.boxsize(nlevels+1))
     1            then     
                     nh(j) = nh(j) + 1
                     isource(istart+nh(j)) = isource(istart+i)
                  else
                     ii = ii+1
                     isrctmp(ii) = isource(istart+i)
                  endif
               enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
                do i=1,ii
                  isource(istart+nh(j)+i) = isrctmp(i)
                enddo
                istart = istart + nsc(j)
            enddo
c           End of sorting sources

c           Sort targets
c           which child?  1,2,3,4,5,6,7,8? counter nt1,nt2,nt3,nt4,
c           nt5,nt6,nt7,nt8
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
            i1234 = itfirst(ibox)-1
            i5678 = 0
            do itt = itfirst(ibox),itlast(ibox)
               if(trg(3,itarget(itt)) - centers(3,ibox).lt.0) then
                 i1234 = i1234+1
                 itarget(i1234) = itarget(itt)
               else
                  i5678 = i5678 + 1
                  itargtmp(i5678) = itarget(itt)
               endif
            enddo
c           Reorder sources to include targets in 5678 in the array
            do i=1,i5678
               itarget(i1234+i) = itargtmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = itfirst(ibox)-1
            i34 = 0
            do itt = itfirst(ibox),i1234
               if(trg(2,itarget(itt))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  itarget(i12) = itarget(itt)
               else
                  i34 = i34 + 1
                  itargtmp(i34) = itarget(itt)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder targets to include 34 in the array
            do i=1,i34
               itarget(i12+i) = itargtmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do itt=i1234+1,itlast(ibox)
               if(trg(2,itarget(itt))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  itarget(i56) = itarget(itt)
               else
                  i78 = i78 + 1
                  itargtmp(i78) = itarget(itt)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               itarget(i56+i) = itargtmp(i)
            enddo
c           End of reordering i5678         

            ntc(1) = 0
            ntc(2) = 0
            ntc(3) = 0
            ntc(4) = 0
            ntc(5) = 0
            ntc(6) = 0
            ntc(7) = 0
            ntc(8) = 0

c           Sort into boxes 1 and 2
            do itt = itfirst(ibox),i12
               if(trg(1,itarget(itt))-centers(1,ibox).lt.0) then
                  itarget(itfirst(ibox)+ntc(1)) = itarget(itt)
                  ntc(1) = ntc(1) + 1
               else
                  ntc(2) = ntc(2) + 1
                  itargtmp(ntc(2)) = itarget(itt)
               endif
            enddo
c           Reorder targets so that sources in 2 are at the
c           end of this part of the array
            do i=1,ntc(2)
               itarget(itfirst(ibox)+ntc(1)+i-1) = itargtmp(i)
            enddo
c           Sort into boxes 3 and 4
            do itt = i12+1, i1234
               if(trg(1,itarget(itt))-centers(1,ibox).lt.0) then
                  itarget(i12+1+ntc(3)) = itarget(itt)
                  ntc(3) = ntc(3) + 1
                else
                   ntc(4) = ntc(4)+1
                   itargtmp(ntc(4)) = itarget(itt)
                endif
            enddo
c           Reorder targets so that sources in 4 are at the
c           end of this part of the array
            do i=1,ntc(4)
               itarget(i12+ntc(3)+i) = itargtmp(i)
            enddo

c           Sort into boxes 5 and 6
            do itt = i1234+1,i56
               if(trg(1,itarget(itt))-centers(1,ibox).lt.0) then
                  itarget(i1234+1+ntc(5)) = itarget(itt)
                  ntc(5) = ntc(5) + 1
               else
                  ntc(6) = ntc(6) + 1
                  itargtmp(ntc(6)) = itarget(itt)
               endif
            enddo
c           Reorder targets so that sources in 6 are at the
c           end of this part of the array
            do i=1,ntc(6)
               itarget(i1234+ntc(5)+i) = itargtmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do itt=i56+1,itlast(ibox)
               if(trg(1,itarget(itt))-centers(1,ibox).lt.0) then
                  itarget(i56+1+ntc(7)) = itarget(itt)
                  ntc(7) = ntc(7) + 1
               else
                  ntc(8) = ntc(8) + 1
                  itargtmp(ntc(8)) = itarget(itt)
               endif
            enddo
c           Reorder targets so that sources in 8 are at the
c           end of the array
            do i=1,ntc(8)
               itarget(i56+ntc(7)+i) = itargtmp(i)
            enddo
c           End of sorting targets

c           Sort expansion centers
c           which child?  1,2,3,4,5,6,7,8? counter
c           nexpcc1, nexpcc2, nexpcc3, nexpcc4,
c           nexpcc5, nexpcc6, nexpcc7, nexpcc8,
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
            i1234 = iefirst(ibox)-1
            i5678 = 0
            do ie = iefirst(ibox),ielast(ibox)
               if(expc(3,iexpc(ie)) - centers(3,ibox).lt.0) then
                  i1234 = i1234+1
                  iexpc(i1234) = iexpc(ie)
               else
                  i5678 = i5678 + 1
                  iexpctmp(i5678) = iexpc(ie)
               endif
            enddo

c           Reorder sources to include sources in 5678 in the array
            do i=1,i5678
               iexpc(i1234+i) = iexpctmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = iefirst(ibox)-1
            i34 = 0
            do ie = iefirst(ibox),i1234
               if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  iexpc(i12) = iexpc(ie)
               else
                  i34 = i34 + 1
                  iexpctmp(i34) = iexpc(ie)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder sources to include 34 in the array
            do i=1,i34
               iexpc(i12+i) = iexpctmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do ie=i1234+1,ielast(ibox)
               if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  iexpc(i56) = iexpc(ie)
               else
                  i78 = i78 + 1
                  iexpctmp(i78) = iexpc(ie)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               iexpc(i56+i) = iexpctmp(i)
            enddo
c           End of reordering i5678         

            nexpcc(1) = 0
            nexpcc(2) = 0
            nexpcc(3) = 0
            nexpcc(4) = 0
            nexpcc(5) = 0
            nexpcc(6) = 0
            nexpcc(7) = 0
            nexpcc(8) = 0

c           Sort into boxes 1 and 2
            do ie = iefirst(ibox),i12
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(iefirst(ibox)+nexpcc(1)) = iexpc(ie)
                  nexpcc(1) = nexpcc(1) + 1
               else
                  nexpcc(2) = nexpcc(2) + 1
                  iexpctmp(nexpcc(2)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 2 are at the
c           end of this part of the array
            do i=1,nexpcc(2)
               iexpc(iefirst(ibox)+nexpcc(1)+i-1) = iexpctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do ie = i12+1, i1234
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i12+1+nexpcc(3)) = iexpc(ie)
                  nexpcc(3) = nexpcc(3) + 1
               else
                  nexpcc(4) = nexpcc(4)+1
                  iexpctmp(nexpcc(4)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 4 are at the
c           end of this part of the array
            do i=1,nexpcc(4)
               iexpc(i12+nexpcc(3)+i) = iexpctmp(i)
            enddo

c           Sort into boxes 5 and 6
            do ie = i1234+1,i56
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i1234+1+nexpcc(5)) = iexpc(ie)
                  nexpcc(5) = nexpcc(5) + 1
               else
                  nexpcc(6) = nexpcc(6) + 1
                  iexpctmp(nexpcc(6)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 6 are at the
c           end of this part of the array
            do i=1,nexpcc(6)
               iexpc(i1234+nexpcc(5)+i) = iexpctmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do ie=i56+1,ielast(ibox)
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i56+1+nexpcc(7)) = iexpc(ie)
                  nexpcc(7) = nexpcc(7) + 1
               else
                  nexpcc(8) = nexpcc(8) + 1
                  iexpctmp(nexpcc(8)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 8 are at the
c           end of the array
            do i=1,nexpcc(8)
               iexpc(i56+nexpcc(7)+i) = iexpctmp(i)
            enddo
c           End of sorting expanison centers

            istart = iefirst(ibox)-1
            do j=1,8
c           check hung -> counter nhc1,nhc2,nhc3,nhc4,nhc5,nhc6,nhc7,nhc8
               ii = 0
               nhc(j) = 0
               do i=1,nexpcc(j)
                  if(radexp(iexpc(istart+i)).gt.boxsize(nlevels+1))
     1            then     
                     nhc(j) = nhc(j) + 1
                     iexpc(istart+nhc(j)) = iexpc(istart+i)
                  else
                     ii = ii+1
                     iexpctmp(ii) = iexpc(istart+i)
                  endif
               enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
               do i=1,ii
                  iexpc(istart+nhc(j)+i) = iexpctmp(i)
               enddo
               istart = istart + nexpcc(j)
            enddo

            nchild(ibox) = 0
c           Create the required boxes
            istart = isfirst(ibox)
            jstart = itfirst(ibox)
            kstart = iefirst(ibox)
            do i=1,8
               ii = 2
               jj = 2
               if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
               if(i.lt.5) jj = 1
               if(nsc(i)+ntc(i)+nexpcc(i).ge.0) then
c                 Increment total number of boxes               
                  nboxes = nboxes + 1
                  if(nboxes.gt.nbmax) then
                    write(*,*) "Exceeding max number of boxes"
                    write(*,*) "Exiting"
                    ier = 12
                    return
                  endif

c                 Increment number of children for the current box
                  nchild(ibox) = nchild(ibox)+1
c                 Update the array of children for the current box
                  ichild(i,ibox) = nboxes
c                 Update the array of levels for the child box
                  ilevel(nboxes) = nlevels+1
c                 Update the array of parents for the child box
                  iparent(nboxes) = ibox
c                 Compute center for the child box
                  centers(1,nboxes) = centers(1,ibox)+(-1)**i*
     1                                boxsize(nlevels+1)/2.0
                  centers(2,nboxes) = centers(2,ibox)+(-1)**ii*
     1                              boxsize(nlevels+1)/2.0
                  centers(3,nboxes) = centers(3,ibox)+(-1)**jj*
     1                              boxsize(nlevels+1)/2.0

                  nchild(nboxes) = 0
                  ichild(1,nboxes) = -1
                  ichild(2,nboxes) = -1
                  ichild(3,nboxes) = -1
                  ichild(4,nboxes) = -1
                  ichild(5,nboxes) = -1
                  ichild(6,nboxes) = -1
                  ichild(7,nboxes) = -1
                  ichild(8,nboxes) = -1
c                 Update arrays ihsfirst,ihslast,isfirst,islast
                  ihsfirst(nboxes) = istart
                  ihslast(nboxes) = istart + nh(i) - 1
                  nhungsrc(nboxes) = nh(i)

                  isfirst(nboxes) = istart + nh(i)
                  islast(nboxes) = istart + nsc(i) - 1

c                 Update arrays itfirst, itlast
                  itfirst(nboxes) = jstart
                  itlast(nboxes) = jstart + ntc(i) - 1

c                 Update arrays ihefirst,ihelast,iefirst,ielast
                  ihefirst(nboxes) = kstart
                  ihelast(nboxes) = kstart + nhc(i)-1
                  nhungexp(nboxes) = nhc(i)

                  iefirst(nboxes) = kstart + nhc(i)
                  ielast(nboxes) = kstart + nexpcc(i) - 1

                  nss = islast(nboxes) - isfirst(nboxes)+1
                  nee = ielast(nboxes) - iefirst(nboxes)+1 
c                 Check if further refinement required
                  if ((idivflag.eq.0).and.(nss.gt.ndiv)) irefine=1
                  if ((idivflag.eq.1).and.(ntc(i).gt.ndiv)) irefine=1
                  if ((idivflag.eq.2).and.(nss.gt.ndiv.or.
     1                     ntc(i).gt.ndiv))  irefine=1
                  if ((idivflag.eq.3).and.(nss.gt.ndiv.or.
     1                      ntc(i).gt.ndiv.or.nee.gt.ndiv)) irefine=1  
               endif
               istart = istart + nsc(i)
               jstart = jstart + ntc(i)
               kstart = kstart + nexpcc(i)
            enddo
         endif
      enddo
      nlevels = nlevels+1
      laddr(1,nlevels) = nbfirst
      laddr(2,nlevels) = nboxes

      return
      end
c------------------------------------------------------------------      

      subroutine computemhung(nlevels,nboxes,laddr,iparent,ilevel,
     1                        nnbors,mnbors,nbors,nlist1,mnlist1,list1,
     2                        nhungsrc,nhunglistsrc,mhung)
c     This subroutine computes mhung for a given tree 
c     and the number of hung sources per box
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c    
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnbors      in: integer
c                 max number of boxes in list1. If isep=1, 
c                 mnbors should be 27 and if isep=2, mnbors
c                 should be 125
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i, nhunglistsrc(i) = 
c                    nhunglistsrc(idad) + \sum_{j=1}{nnbors(i)}
c                    nhungsrc(nbors(j,i))
c
c      mhung         out: integer
c                    max(nhunglistsrc)
      implicit none
      integer nlevels,nboxes, mnbors,mnlist1
      integer laddr(2,0:nlevels), ilevel(nboxes)
      integer iparent(nboxes),nnbors(nboxes), nbors(mnbors,nboxes)
      integer nlist1(nboxes),list1(mnlist1,nboxes)
      integer nhungsrc(nboxes), nhunglistsrc(nboxes)
      integer mhung

c     Temporary variables
      integer ilev,ibox,i,dad,jbox

c     initialize for level 0      
      nhunglistsrc(1) = nhungsrc(1)
      do ilev = 1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,dad,i,jbox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            dad = iparent(ibox)
            nhunglistsrc(ibox) = nhunglistsrc(dad)
            do i=1,nnbors(ibox)
               jbox = nbors(i,ibox)
               nhunglistsrc(ibox) = nhunglistsrc(ibox) +
     1                              nhungsrc(jbox)
            enddo
            do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                if(ilevel(ibox).lt.ilevel(jbox)) nhunglistsrc(ibox)=
     1            nhunglistsrc(ibox) + nhungsrc(jbox)
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      mhung = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(max:mhung) PRIVATE(i)      
      do ibox = 1,nboxes
         if(nhunglistsrc(ibox).gt.mhung) mhung = nhunglistsrc(ibox)
      enddo
C$OMP END PARALLEL DO      
      


      return
      end
c------------------------------------------------------------------
      subroutine computehunglist(mhung,nlevels,nboxes,laddr,ns,
     1           isource,iparent,ihsfirst,ihslast,ilevel,nnbors,mnbors,
     2           nbors,nlist1,mnlist1,list1,
     3           nhungsrc,nhunglistsrc,ihunglistsrc)
c     This subroutine computes the hung list for each box 
c     and stores it. The hung list of sources of a box is
c     the hunglist of the parent + the sources hung in list1
c     of the box
c 
c     INPUT arguments
c     mhung         in: Integer
c                   max number of hung chunks for a box
c
c     nlevels       in: Integer
c                   number of levels in the tree
c
c     nboxes        in: Integer
c                   number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     ns          in: integer
c                 number of sources
c
c     isource     in: integer(ns)
c                 Mapping from tree sorted array of sources
c                 to the user prescribed ordering
c  
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c     
c     ihsfirst    in: integer(nboxes)
c                 ihsfirst(i) points to the first hung source
c                 in box i in the sorted array
c
c      ihslast    in: integer(nboxes)
c                 ihslast(i) points to the last hung source
c                 in box i in the sorted array
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnbors      in: integer
c                 max number of boxes in list1. If isep=1, 
c                 mnbors should be 27 and if isep=2, mnbors
c                 should be 125
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i 
c
c      ihunglistsrc  out: integer(mhung,nboxes)
c                    ihunglistsrc(j,i) is the id of the jth
c                    hung source relevant to box i
c-----------------------------------------------------------------
       implicit none
       integer mhung,nlevels,nboxes,ns, mnbors, mnlist1
       integer laddr(2,0:nlevels)
       integer iparent(nboxes),ilevel(nboxes)
       integer isource(ns), ihsfirst(nboxes),ihslast(nboxes)
       integer nnbors(nboxes), nbors(mnbors,nboxes)
       integer nlist1(nboxes), list1(mnlist1,nboxes)
       integer nhungsrc(nboxes),nhunglistsrc(nboxes)
       integer ihunglistsrc(mhung,nboxes)
c      Temp variables
       integer i,j,ibox,jbox,ilev,dad

c     initialize for level 0      
      nhunglistsrc(1) = nhungsrc(1)
      do i=1,nhungsrc(1)
          ihunglistsrc(i,1) = isource(ihsfirst(1)+i-1)
      enddo
      do ilev = 1,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            dad = iparent(ibox)
            do i=1,nhunglistsrc(dad)
               ihunglistsrc(i,ibox) = ihunglistsrc(i,dad)
            enddo
            nhunglistsrc(ibox) = nhunglistsrc(dad)
            do i=1,nnbors(ibox)
               jbox = nbors(i,ibox)
               do j=1,nhungsrc(jbox)
                  ihunglistsrc(nhunglistsrc(ibox)+j,ibox) = 
     1                  isource(ihsfirst(jbox)+j-1)
               enddo
               nhunglistsrc(ibox) = nhunglistsrc(ibox) +
     1                              nhungsrc(jbox)
            enddo
            do i=1,nlist1(ibox)
                jbox = list1(i,ibox)
                if(ilevel(ibox).lt.ilevel(jbox)) then
                   do j=1,nhungsrc(jbox)
                      ihunglistsrc(nhunglistsrc(ibox)+j,ibox)=
     1                    isource(ihsfirst(jbox)+j-1)
                   enddo
                   nhunglistsrc(ibox)=
     1             nhunglistsrc(ibox) + nhungsrc(jbox)
                endif
            enddo
         enddo
      enddo

      return
      end
c------------------------------------------------------------------

      subroutine mklraptreemem(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,isep,
     $                   nlmax,nbmax,nlevels,nboxes,mnbors,mnlist1,
     $                   mnlist2,mnlist3,mnlist4,mhung,ltree)

      implicit none
      integer ier
      integer ns,nt,nexpc,idivflag,ndiv,mhung,isep
      integer nlevels,nboxes,lcenters
      integer *8 ltree,nboxes8
      integer i,j,nbmax,nlmax
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: ilevel(:),ilevel2(:)
      integer, allocatable :: iparenttemp(:),iparenttemp2(:)
      integer, allocatable :: nchild(:),nchild2(:)
      integer, allocatable :: ichildtemp(:,:),ichildtemp2(:,:)
      integer, allocatable :: nnbors(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:),ihsfirsttemp2(:)
      integer, allocatable :: ihslasttemp(:),ihslasttemp2(:)
      integer, allocatable :: isfirsttemp(:),isfirsttemp2(:)
      integer, allocatable :: islasttemp(:),islasttemp2(:)
      integer, allocatable :: itfirsttemp(:),itfirsttemp2(:)
      integer, allocatable :: itlasttemp(:),itlasttemp2(:)
      integer, allocatable :: ihefirsttemp(:),ihefirsttemp2(:)
      integer, allocatable :: ihelasttemp(:),ihelasttemp2(:)
      integer, allocatable :: iefirsttemp(:),iefirsttemp2(:)
      integer, allocatable :: ielasttemp(:),ielasttemp2(:)
      integer, allocatable :: nhungsrc(:),nhungsrc2(:)
      integer, allocatable :: nhungexp(:),nhungexp2(:)
      integer, allocatable :: nhunglistsrc(:)

      double precision boxsize(0:nlmax)
      double precision src(3,ns),radsrc(ns)
      double precision trg(3,nt)
      double precision, allocatable :: centers(:,:),centers2(:,:)
      double precision expc(3,nexpc)
      double precision radexp(nexpc)
c
      double precision xmin,xmax,ymin,ymax,zmin,zmax,sizex,sizey,sizez
      double precision btmp
      integer ictr,ih,irefine,is,ie
      integer nss,nee,ntot,nbtmp,ntt
      integer ibox, ifirstbox,ilastbox, nnew, nbtot

      integer mnbors, mnlist1, mnlist2,mnlist3,mnlist4
      integer, allocatable :: nlist1(:)
      integer, allocatable :: list1(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: nlist3(:)
      integer, allocatable :: list3(:,:)
      integer, allocatable :: nlist4(:)
      integer, allocatable :: list4(:,:)

c     This code is a memory code for mkptree. This
c     subroutine determines the number of boxes (lcenters) required,
c     the length of the tree array (ltree) and the maximum number
c     of hung chunks (mhung)
c     
c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     radexp        expansion center radius
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     nlmax         max number of levels
c     nbmax         max number of boxes (no longer in use!)
c
c     OUTPUT:
c     ier           error code
c                    ier = 0, normal execution
c                    ier = 4, isep != 1 or 2 
c                    ier = 8, sources and targets too close
c     nlevels       number of levels
c     nboxes        number of boxes
c     mhung         max number of hung chunks
c     ltree         length of the tree
c
c     mhung is the max number of hung sources on output
c
c     This subroutine essentially creates the quad tree
c     and in process determines the memory required for the tree.
c
c     It can be shown that the maximum number of boxes is less
c     than 2*(nexpc + ns + nt)
c
c     For details of the variables used in the tree, refer
c     to the main routine
c---------------------------------------------------------------------
c
c     initialize temporary arrays
c  
c     assumes nlevels lt 200
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c
      ier = 0
      ntot = ns + nt + nexpc
      nbmax = 10000
      nboxes = nbmax
      lcenters = nbmax
      nlevels = nlmax
      ier = 0
c
      allocate(laddr(2,0:nlmax))
      allocate(ilevel(nbmax))
      allocate(iparenttemp(nbmax))
      allocate(nchild(nbmax))
      allocate(ichildtemp(8,nbmax))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nbmax))
      allocate(ihslasttemp(nbmax))
      allocate(isfirsttemp(nbmax))
      allocate(islasttemp(nbmax))
      allocate(itfirsttemp(nbmax))
      allocate(itlasttemp(nbmax))
      allocate(ihefirsttemp(nbmax))
      allocate(ihelasttemp(nbmax))
      allocate(iefirsttemp(nbmax))
      allocate(ielasttemp(nbmax))
      allocate(nhungsrc(nbmax))
      allocate(nhungexp(nbmax))
      allocate(centers(3,lcenters))

      if(isep.eq.1) then
         mnbors = 27
      endif

      if(isep.eq.2) then
         mnbors = 125
      endif

      if(isep.ne.1.and.isep.ne.2) then
         ier = 4
cc         call prinf('Error tree not allocated, isep.ne.1,2*',ier,0)
         print *, "Error tree not allocated, isep.ne.1,2"
         return
      endif

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo

      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
        if(expc(3,i) .lt. zmin) zmin=expc(3,i)
        if(expc(3,i) .gt. zmax) zmax=expc(3,i)
      enddo
      sizex=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin
      boxsize(0) = sizex
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
      if(sizez .gt. boxsize(0)) boxsize(0)=sizez

      btmp = sqrt(sizex**2+sizey**2+sizez**2)
      if(boxsize(0)/btmp.le.1.0d-16) then
        ier = 8
        write(*,*) "Nothing to compute"
        write(*,*) "Sources and targets too close"
        write(*,*) "Exiting"
      endif
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      centers(3,1)=(zmin+zmax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
      nhungsrc(1) = 0
      nhungexp(1) = 0
c
c     count number of hung sources
c     and hang up "big" sources
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif

c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp 
      do i = 1,nt
         itargettemp(i) = i
      enddo
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1
      ichildtemp(5,1) = -1
      ichildtemp(6,1) = -1
      ichildtemp(7,1) = -1
      ichildtemp(8,1) = -1
c
      nlevels = 0
      nboxes = 1
c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1



      do i = 1,nlmax
         if (irefine.eq.1) then
c
c
c            estimate number of new boxes to be created
c
           nnew = 0
           ifirstbox = laddr(1,nlevels)
           ilastbox =  laddr(2,nlevels)

           do ibox = ifirstbox,ilastbox
               nss = islasttemp(ibox) - isfirsttemp(ibox) + 1
               ntt = itlasttemp(ibox) - itfirsttemp(ibox) + 1
               nee = ielasttemp(ibox) - iefirsttemp(ibox) + 1
               if((idivflag.eq.0).and.(nss.gt.ndiv)) nnew = nnew + 1
               if((idivflag.eq.1).and.(ntt.gt.ndiv)) nnew = nnew + 1
               if((idivflag.eq.2).and.(nss.gt.ndiv.or.ntt.gt.ndiv)) 
     1            nnew = nnew + 1
               if((idivflag.eq.3).and.(nss.gt.ndiv.or.
     1           ntt.gt.ndiv.or.nee.gt.ndiv)) nnew = nnew + 1
            enddo

            nbtot = nboxes + 8*nnew

c
c            if current memory is not sufficient, 
c            delete previous arrays and allocate more memory
c            allocate more memory
c
            if(nbtot.gt.nbmax) then
              nbmax = nbtot
              lcenters = nbtot
              allocate(centers2(3,nboxes))
              allocate(ilevel2(nboxes),iparenttemp2(nboxes))
              allocate(nchild2(nboxes),ichildtemp2(8,nboxes))
              allocate(ihsfirsttemp2(nboxes),ihslasttemp2(nboxes))
              allocate(isfirsttemp2(nboxes),islasttemp2(nboxes))
              allocate(itfirsttemp2(nboxes),itlasttemp2(nboxes))
              allocate(ihefirsttemp2(nboxes),ihelasttemp2(nboxes))
              allocate(iefirsttemp2(nboxes),ielasttemp2(nboxes))
              allocate(nhungsrc2(nboxes),nhungexp2(nboxes))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
              do ibox = 1,nboxes
                ilevel2(ibox) = ilevel(ibox)
                iparenttemp2(ibox) = iparenttemp(ibox)
                nchild2(ibox) = nchild(ibox)
                ihsfirsttemp2(ibox) = ihsfirsttemp(ibox)
                ihslasttemp2(ibox) = ihslasttemp(ibox)
                isfirsttemp2(ibox) = isfirsttemp(ibox)
                islasttemp2(ibox) = islasttemp(ibox)
                itfirsttemp2(ibox) = itfirsttemp(ibox)
                itlasttemp2(ibox) = itlasttemp(ibox)
                ihefirsttemp2(ibox) = ihefirsttemp(ibox)
                ihelasttemp2(ibox) = ihelasttemp(ibox)
                iefirsttemp2(ibox) = iefirsttemp(ibox)
                ielasttemp2(ibox) = ielasttemp(ibox)
                nhungsrc2(ibox) = nhungsrc(ibox)
                nhungexp2(ibox) = nhungexp(ibox)
                do j=1,3
                  centers2(j,ibox) = centers(j,ibox)
                enddo

                do j=1,8
                  ichildtemp2(j,ibox) = ichildtemp(j,ibox)
                enddo

              enddo
C$OMP END PARALLEL DO              
              deallocate(centers,ilevel,iparenttemp,nchild,ichildtemp)
              deallocate(ihsfirsttemp,ihslasttemp)
              deallocate(isfirsttemp,islasttemp)
              deallocate(itfirsttemp,itlasttemp)
              deallocate(ihefirsttemp,ihelasttemp)
              deallocate(iefirsttemp,ielasttemp)
              deallocate(nhungsrc,nhungexp)

              allocate(centers(3,nbmax))
              allocate(ilevel(nbmax),iparenttemp(nbmax))
              allocate(nchild(nbmax),ichildtemp(8,nbmax))
              allocate(ihsfirsttemp(nbmax),ihslasttemp(nbmax))
              allocate(isfirsttemp(nbmax),islasttemp(nbmax))
              allocate(itfirsttemp(nbmax),itlasttemp(nbmax))
              allocate(ihefirsttemp(nbmax),ihelasttemp(nbmax))
              allocate(iefirsttemp(nbmax),ielasttemp(nbmax))
              allocate(nhungsrc(nbmax),nhungexp(nbmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
              do ibox = 1,nboxes
                ilevel(ibox) = ilevel2(ibox)
                iparenttemp(ibox) = iparenttemp2(ibox)
                nchild(ibox) = nchild2(ibox)
                ihsfirsttemp(ibox) = ihsfirsttemp2(ibox)
                ihslasttemp(ibox) = ihslasttemp2(ibox)
                isfirsttemp(ibox) = isfirsttemp2(ibox)
                islasttemp(ibox) = islasttemp2(ibox)
                itfirsttemp(ibox) = itfirsttemp2(ibox)
                itlasttemp(ibox) = itlasttemp2(ibox)
                ihefirsttemp(ibox) = ihefirsttemp2(ibox)
                ihelasttemp(ibox) = ihelasttemp2(ibox)
                iefirsttemp(ibox) = iefirsttemp2(ibox)
                ielasttemp(ibox) = ielasttemp2(ibox)
                nhungsrc(ibox) = nhungsrc2(ibox)
                nhungexp(ibox) = nhungexp2(ibox)
                do j=1,3
                  centers(j,ibox) = centers2(j,ibox)
                enddo
                do j=1,8
                  ichildtemp(j,ibox) = ichildtemp2(j,ibox)
                enddo
              enddo
C$OMP END PARALLEL DO              



              deallocate(centers2,ilevel2,iparenttemp2,nchild2)
              deallocate(ichildtemp2)
              deallocate(ihsfirsttemp2,ihslasttemp2)
              deallocate(isfirsttemp2,islasttemp2)
              deallocate(itfirsttemp2,itlasttemp2)
              deallocate(ihefirsttemp2,ihelasttemp2)
              deallocate(iefirsttemp2,ielasttemp2)
              deallocate(nhungsrc2,nhungexp2)
            endif

            

            call subdivide_adap(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,ilevel,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)
            if(ier.ne.0) return
         else
            exit
         endif
      enddo


c
c
c        check with Leslie/Dhairya/Alex/Zydrunas if 16 is enough
c

      nbtot = 16*nboxes
c
c            if current memory is not sufficient, 
c            delete previous arrays and allocate more memory
c            allocate more memory
c
      if(nbtot.gt.nbmax) then
        nbmax = nbtot
        lcenters = nbtot

        allocate(centers2(3,nboxes))
        allocate(ilevel2(nboxes),iparenttemp2(nboxes))
        allocate(nchild2(nboxes),ichildtemp2(8,nboxes))
        allocate(ihsfirsttemp2(nboxes),ihslasttemp2(nboxes))
        allocate(isfirsttemp2(nboxes),islasttemp2(nboxes))
        allocate(itfirsttemp2(nboxes),itlasttemp2(nboxes))
        allocate(ihefirsttemp2(nboxes),ihelasttemp2(nboxes))
        allocate(iefirsttemp2(nboxes),ielasttemp2(nboxes))
        allocate(nhungsrc2(nboxes),nhungexp2(nboxes))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
        do ibox = 1,nboxes
          ilevel2(ibox) = ilevel(ibox)
          iparenttemp2(ibox) = iparenttemp(ibox)
          nchild2(ibox) = nchild(ibox)
          ihsfirsttemp2(ibox) = ihsfirsttemp(ibox)
          ihslasttemp2(ibox) = ihslasttemp(ibox)
          isfirsttemp2(ibox) = isfirsttemp(ibox)
          islasttemp2(ibox) = islasttemp(ibox)
          itfirsttemp2(ibox) = itfirsttemp(ibox)
          itlasttemp2(ibox) = itlasttemp(ibox)
          ihefirsttemp2(ibox) = ihefirsttemp(ibox)
          ihelasttemp2(ibox) = ihelasttemp(ibox)
          iefirsttemp2(ibox) = iefirsttemp(ibox)
          ielasttemp2(ibox) = ielasttemp(ibox)
          nhungsrc2(ibox) = nhungsrc(ibox)
          nhungexp2(ibox) = nhungexp(ibox)
          do j=1,3
            centers2(j,ibox) = centers(j,ibox)
          enddo
           do j=1,8
            ichildtemp2(j,ibox) = ichildtemp(j,ibox)
          enddo
         enddo
C$OMP END PARALLEL DO              
        deallocate(centers,ilevel,iparenttemp,nchild,ichildtemp)
        deallocate(ihsfirsttemp,ihslasttemp)
        deallocate(isfirsttemp,islasttemp)
        deallocate(itfirsttemp,itlasttemp)
        deallocate(ihefirsttemp,ihelasttemp)
        deallocate(iefirsttemp,ielasttemp)
        deallocate(nhungsrc,nhungexp)

        allocate(centers(3,nbmax))
        allocate(ilevel(nbmax),iparenttemp(nbmax))
        allocate(nchild(nbmax),ichildtemp(8,nbmax))
        allocate(ihsfirsttemp(nbmax),ihslasttemp(nbmax))
        allocate(isfirsttemp(nbmax),islasttemp(nbmax))
        allocate(itfirsttemp(nbmax),itlasttemp(nbmax))
        allocate(ihefirsttemp(nbmax),ihelasttemp(nbmax))
        allocate(iefirsttemp(nbmax),ielasttemp(nbmax))
        allocate(nhungsrc(nbmax),nhungexp(nbmax))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
        do ibox = 1,nboxes
          ilevel(ibox) = ilevel2(ibox)
          iparenttemp(ibox) = iparenttemp2(ibox)
          nchild(ibox) = nchild2(ibox)
          ihsfirsttemp(ibox) = ihsfirsttemp2(ibox)
          ihslasttemp(ibox) = ihslasttemp2(ibox)
          isfirsttemp(ibox) = isfirsttemp2(ibox)
          islasttemp(ibox) = islasttemp2(ibox)
          itfirsttemp(ibox) = itfirsttemp2(ibox)
          itlasttemp(ibox) = itlasttemp2(ibox)
          ihefirsttemp(ibox) = ihefirsttemp2(ibox)
          ihelasttemp(ibox) = ihelasttemp2(ibox)
          iefirsttemp(ibox) = iefirsttemp2(ibox)
          ielasttemp(ibox) = ielasttemp2(ibox)
          nhungsrc(ibox) = nhungsrc2(ibox)
          nhungexp(ibox) = nhungexp2(ibox)
          do j=1,3
            centers(j,ibox) = centers2(j,ibox)
          enddo
          do j=1,8
            ichildtemp(j,ibox) = ichildtemp2(j,ibox)
          enddo
        enddo
C$OMP END PARALLEL DO              

        deallocate(centers2,ilevel2,iparenttemp2,nchild2)
        deallocate(ichildtemp2)
        deallocate(ihsfirsttemp2,ihslasttemp2)
        deallocate(isfirsttemp2,islasttemp2)
        deallocate(itfirsttemp2,itlasttemp2)
        deallocate(ihefirsttemp2,ihelasttemp2)
        deallocate(iefirsttemp2,ielasttemp2)
        deallocate(nhungsrc2,nhungexp2)
      endif

c     Set up computation of list1 and list2      
      allocate(nnbors(nbmax))
      allocate(nbors(mnbors,nbmax))

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo

      call computecoll(nlevels,nboxes,laddr,boxsize,centers,
     1     iparenttemp,nchild,ichildtemp,mnbors,nnbors,nbors)

      if(nlevels.ge.2) then
      call d3hpfixtree(ier,src,ns,radsrc,trg,nt,expc,nexpc,radexp,
     1     nlevels,nboxes,centers,boxsize,nbmax,laddr,ilevel,
     2     iparenttemp,nchild,ichildtemp,mnbors,nnbors,nbors,
     3     isourcetemp,itargettemp,iexpctemp,ihsfirsttemp,ihslasttemp,
     4     isfirsttemp,islasttemp,itfirsttemp,itlasttemp,ihefirsttemp,
     5     ihelasttemp,iefirsttemp,ielasttemp,nhungsrc,nhungexp)

      if(ier.ne.0) return
      endif


      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
      call computecollisep(nlevels,nboxes,laddr,boxsize,centers,
     1     iparenttemp,nchild,ichildtemp,isep,mnbors,nnbors,nbors)
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
c     Compute mnlist1, mnlist2, mnlist3, mnlist4      
      call computemnlists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,isep,nnbors,mnbors,nbors,mnlist1,
     3                   mnlist2,mnlist3,mnlist4)


      allocate(nlist1(nboxes),nlist2(nboxes))
      allocate(nlist3(nboxes),nlist4(nboxes))
      allocate(list1(mnlist1,nboxes),list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes),list4(mnlist4,nboxes))
      call computelists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,isep,nnbors,mnbors,nbors,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)

      allocate(nhunglistsrc(nboxes))

      call computemhung(nlevels,nboxes,laddr,iparenttemp,ilevel,nnbors,
     1  mnbors,nbors,nlist1,mnlist1,list1,nhungsrc,nhunglistsrc,mhung)
      
c
c      temporarily typecast nboxes to integer *8
c
      nboxes8 = nboxes
      ltree = (28+mhung+mnlist1+mnlist2+mnlist3+mnlist4+
     1   mnbors)*nboxes8+2*(nlevels+1)+ns+nt+nexpc

      return
      end
c--------------------------------------------      

      subroutine mkpwlists(isep,nlevels,laddr,boxsize,
     1                   nboxes,centers,mnlist2,nlist2,list2,
     2                   nupm,ndownm,nnorthm,nsouthm,
     3                   neastm,nwestm,nup,uplist,ndown,downlist,
     3                   nnorth,northlist,nsouth,southlist,
     4                   neast,eastlist,nwest,westlist)

c     This subroutine computes the uplist, downlist, eastlist,
c     west list, north list and south list based on the list2 of
c     the tree structure. The up and down list for a box B
c     are the set of boxes in list 2 
c     that are separated by at least isep boxes
c     in the +z direction and -z direction respectively. The north
c     and south list of box B are those boxes which are not in
c     the up and down list and are well separated from the box B
c     by at least isep boxes in the +y and -y direction
c     respectively. The east and west list of box B are those
c     boxes which are in none of up,down,north and south list of
c     the box B and are well separated from box B by at least
c     isep boxes in the +x and -x direction respectively.
c
c     NOTE: to determine nupm,ndownm,nnorthm,nsouthm,neastm and
c           nwestm use the subroutine findnudnsew
c
c     INPUT arguments:
c     isep        in: Integer
c                 separation parameter
c
c     nlevels     in: integer
c                 number of levels
c
c     laddr       in: integer (2,0:nlevels)
c                 laddr(1,i) is the first box at level i
c                 laddr(2,i) is the last box at level i
c
c     boxsize     in: double precision(0:nlevels)
c                 boxsize(i) is the size of a box at level i
c
c     nboxes      in: Integer
c                 number of boxes
c
c     centers     in: double precision(3,nboxes)
c                 co-ordinates of the box centers in the
c                 tree structure
c
c     mnlist2     in: integer
c                 maximum number of elements in list2 of a box
c  
c     nlist2      in: integer(nboxes)
c                 nlist2(i) is the number of boxes in list2
c                 of box i
c
c     list2       in: integer(mnlist2,nboxes)
c                 list2(j,i) is the box number of the jth box
c                 in the list2 of box i
c
c     nupm        in: integer
c                 max number of boxes in uplist of any box
c
c     ndownm      in: integer
c                 max number of boxes in downlist of any box
c
c     nnorthm     in: integer
c                 max number of boxes in the north list of any box
c
c     nsouthm     in: integer
c                 max number of boxes in the south list of any box
c
c     neastm      in: integer
c                 max number of boxes in the east list of any box
c
c     nwestm      in: integer
c                 max number of boxes in the west list of any box
c
c     OUTPUT
c     nup         out: integer(nboxes)
c                 nup(i) is the number of boxes in the uplist
c                 of box i
c
c     uplist      out: integer(nupm,nboxes)
c                 uplist(j,i) is the jth box in the uplist
c                 of box i
c
c     ndown       out: integer(nboxes)
c                 ndown(i) is the number of boxes in the downlist
c                 of box i
c
c     downlist    out: integer(ndownm,nboxes)
c                 downlist(j,i) is the jth box in the downlist
c                 of box i
c
c     nnorth      out: integer(nboxes)
c                 nnorth(i) is the number of boxes in the northlist
c                 of box i
c
c     northlist   out: integer(nnorthm,nboxes)
c                 northlist(j,i) is the jth box in the northlist
c                 of box i
c
c     nsouth      out: integer(nboxes)
c                 nsouth(i) is the number of boxes in the southlist
c                 of box i
c
c     southlist   out: integer(ndownm,nboxes)
c                 southlist(j,i) is the jth box in the southlist
c                 of box i
c
c     neast       out: integer(nboxes)
c                 neast(i) is the number of boxes in the eastlist
c                 of box i
c
c     eastlist    out: integer(neastm,nboxes)
c                 eastlist(j,i) is the jth box in the eastlist
c                 of box i
c
c     nwest       out: integer(nboxes)
c                 nwest(i) is the number of boxes in the westlist
c                 of box i
c
c     westlist    out: integer(nwestm,nboxes)
c                 westlist(j,i) is the jth box in the westlist
c                 of box i
c---------------------------------------------------------------

      implicit none
      integer isep,nlevels,nboxes,mnlist2,nlist2(nboxes)
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      integer list2(mnlist2,nboxes)
      integer nup(nboxes),ndown(nboxes)
      integer nnorth(nboxes),nsouth(nboxes)
      integer neast(nboxes),nwest(nboxes)
      integer nupm,ndownm,nnorthm,nsouthm,neastm,nwestm
      integer uplist(nupm,nboxes),downlist(ndownm,nboxes)
      integer northlist(nnorthm,nboxes),southlist(nsouthm,nboxes)
      integer eastlist(neastm,nboxes),westlist(nwestm,nboxes)
      integer ilev,i,ibox,jbox
      double precision centers(3,nboxes)

      do ilev = 0,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nup(ibox) = 0
            ndown(ibox) = 0
            nnorth(ibox) = 0
            nsouth(ibox) = 0
            neast(ibox) = 0
            nwest(ibox) = 0
            do i=1,nlist2(ibox)
              jbox = list2(i,ibox)
              if(centers(3,jbox) - centers(3,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nup(ibox) = nup(ibox)+1
                 uplist(nup(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(3,jbox) - centers(3,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 ndown(ibox) = ndown(ibox)+1
                 downlist(ndown(ibox),ibox) = jbox
                 goto 1111
              endif

              if(centers(2,jbox) - centers(2,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nnorth(ibox) = nnorth(ibox)+1
                 northlist(nnorth(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(2,jbox) - centers(2,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nsouth(ibox) = nsouth(ibox)+1
                 southlist(nsouth(ibox),ibox) = jbox
                 goto 1111
              endif

              if(centers(1,jbox) - centers(1,ibox).ge.
     1          1.01d0*isep*boxsize(ilev)) then
                 neast(ibox) = neast(ibox)+1
                 eastlist(neast(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(1,jbox) - centers(1,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nwest(ibox) = nwest(ibox)+1
                 westlist(nwest(ibox),ibox) = jbox
                 goto 1111
              endif
1111          continue                 
           enddo
        enddo
      enddo

      return
      end

      subroutine findnudnsew(isep,nlevels,laddr,boxsize,
     1                       nboxes,centers,mnlist2,
     2                       nlist2,list2,nup,ndown,nnorth,
     3                       nsouth,neast,nwest)
  
c     This subroutine computes the maximum number of elements
c     in the uplist,downlist, eastlist, west list, north list
c     and south list based on the list2 of the tree structure.
c     For detailed description on the lists, refer to mkpwlists
c
c     INPUT arguments:
c     isep        in: Integer
c                 separation parameter
c
c     nlevels     in: integer
c                 number of levels
c
c     laddr       in: integer (2,0:nlevels)
c                 laddr(1,i) is the first box at level i
c                 laddr(2,i) is the last box at level i
c
c     boxsize     in: double precision(0:nlevels)
c                 boxsize(i) is the size of a box at level i
c
c     nboxes      in: Integer
c                 number of boxes
c
c     centers     in: double precision(3,nboxes)
c                 co-ordinates of the box centers in the
c                 tree structure
c
c     mnlist2     in: integer
c                 maximum number of elements in list2 of a box
c  
c     nlist2      in: integer(nboxes)
c                 nlist2(i) is the number of boxes in list2
c                 of box i
c
c     list2       in: integer(mnlist2,nboxes)
c                 list2(j,i) is the box number of the jth box
c                 in the list2 of box i
c    
c     OUTPUT
c     nup         out: integer
c                 max number of boxes in the uplist of any box
c
c     ndown       out: integer
c                 max number of boxes in the downlist of any box
c
c     nnorth      out: integer
c                 max number of boxes in the northlist of any box
c
c     nsouth      out: integer
c                 max number of boxes in the southlist of any box
c
c     neast       out: integer
c                 max number of boxes in the eastlist of any box
c
c     nwest       out: integer
c                 max number of boxes in the westlist of any box
c
c----------------------------------------------------------------
      implicit none
      integer isep,nlevels,nboxes,mnlist2,nlist2(nboxes)
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      integer list2(mnlist2,nboxes)
      integer nup,ndown,nnorth,nsouth,neast,nwest
      integer nupt,ndownt,nnortht,nsoutht,neastt,nwestt
      integer ilev,i,ibox,jbox
      double precision centers(3,nboxes)

      nup = 0
      ndown = 0
      nnorth = 0
      nsouth = 0
      neast = 0
      nwest = 0
      do ilev = 2,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nupt = 0
            ndownt = 0
            nnortht = 0
            nsoutht = 0
            neastt = 0
            nwestt = 0
            do i=1,nlist2(ibox)
              jbox = list2(i,ibox)
              if(centers(3,jbox) - centers(3,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nupt = nupt+1
                 goto 1111
              endif
              if(centers(3,jbox) - centers(3,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 ndownt = ndownt+1
                 goto 1111
              endif

              if(centers(2,jbox) - centers(2,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nnortht = nnortht+1
                 goto 1111
              endif
              if(centers(2,jbox) - centers(2,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nsoutht = nsoutht+1
                 goto 1111
              endif

              if(centers(1,jbox) - centers(1,ibox).ge.
     1          1.01d0*isep*boxsize(ilev)) then
                 neastt = neastt+1
                 goto 1111
              endif
              if(centers(1,jbox) - centers(1,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nwestt = nwestt+1
                 goto 1111
              endif
1111          continue                 
           enddo
           if(nupt.gt.nup) nup = nupt
           if(ndownt.gt.ndown) ndown = ndownt
           if(nnortht.gt.nnorth) nnorth = nnortht
           if(nsoutht.gt.nsouth) nsouth = nsoutht
           if(neastt.gt.neast) neast = neastt
           if(nwestt.gt.nwest) nwest = nwestt
        enddo
      enddo

      return
      end
c-------------------------------------------------------------      
      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,
     2                       mnbors,nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     mnbors      in: integer
c                 max number of neighbors = 27
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(mnbors,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes, mnbors
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes)
      integer nbors(mnbors,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev)).and.
     4                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*boxsize(ilev))) then
                     
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

      return
      end
c-------------------------------------------------------------      
      subroutine computecollisep(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,isep,
     2                       mnbors,nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     isep        in: integer
c                 separation criterion for neighbors.
c                 neighbor is box at same level satisfying
c                 |c_{i} - c_{j}| leq isep*boxsize(c_{i})
c
c     mnbors      in: integer
c                 max number of neighbors = 27
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(mnbors,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes, mnbors
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes)
      integer nbors(mnbors,nboxes)
      integer isep

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*isep*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*isep*boxsize(ilev)).and.
     4                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*isep*boxsize(ilev))) then
                     
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

      return
      end
c-------------------------------------------------------------------      
      subroutine d3hpfixtree(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $              radexp,nlevels,nboxes,
     $              centers,boxsize,nbmax,
     $              laddr,ilevel,iparent,nchild,ichild,
     $              mnbors,nnbors,nbors, 
     $              isource,itarget,iexpc,
     $              ihsfirst,ihslast,
     $              isfirst,islast,
     $              itfirst,itlast,
     $              ihefirst,ihelast,
     $              iefirst,ielast,nhungsrc,
     $              nhungexp)
c     Give an adaptive tree, this subroutine fixes the tree
c     to make it level restricted. A level restricted tree
c     is one where no two boxes that contact each other are
c     more than one level apart. The subroutine corrects
c     a given adaptive tree by adding in new boxes. The
c     process involves flagging down bigger boxes and dividing
c     them and their children as necessary.
c
c     Input arguments
c     src         in: double precision(3,ns)
c                 x and y coordinates of the source locations
c    
c     ns          in: integer
c                 number of sources
c
c     radsrc      in: double precision(ns)
c                 radius of sources. A source i is hung
c                 at level ilev 
c                 if radsrc(i)>boxsize(ilev)
c
c     trg        in: double precision(3,nt)
c                 x and y coordinates of target locations
c
c     nt          in: integer
c                 number of targets
c
c     expc        in: double precision(3,nexpc)
c                 x and y coordinates of expansion centers
c
c     radexp      in: double precision(nexpc)
c                 radius of expansion centers. An expansion
c                 center i is hung at level ilev if
c                 radexp(i)>boxsize(ilev)
c
c     nlevels     in/out: integer
c                 current number of levels in the tree
c                 On output nlevels remains unchanged
c
c     nboxes      in/out: integer
c                 current number of boxes in the tree.
c                 On output it is the number of boxes
c                 in the level restricted tree
c
c     centers     in/out: double precision(3,nbmax)
c                 x and y coordinates of the box centers 
c                 in the tree
c
c     boxsize     in: double precision(0:nlevels)
c                 size of the box at any given level
c
c     nbmax       in: integer
c                 max number of boxes
c
c     laddr       in/out: integer(2,0:nlevels)
c                 laddr(1,i), laddr(2,i) are the first
c                 and the last box at level i
c
c     ilevel      in/out: integer(nbmax)
c                 ilevel(i) is the level of box i
c
c     iparent     in/out: integer(nbmax)
c                 iparent(i) is the parent of box i
c
c     nchild      in/out: integer(nbmax)
c                 nchild(i) is the number of children 
c                 of box i
c
c    ichild       in/out: integer(8,nbmax)
c                 ichild(j,i) is the jth child of box i
c
c    nnbors       in/out: integer(nbmax)
c                 nnbors(i) is the number of colleagues of box i
c
c    nbors        in/out: integer(mnbors,nbmax)
c                 nbors(j,i) is the jth colleague of box i
c
c    isource      in/out: integer(ns)
c                 tree sorted array of sources
c
c    itarg        in/out: integer(nt)
c                 tree sorted array of targets
c
c    iexpc        in/out: integer(nexpc)
c                 tree sorted array of expansion centers
c
c    ihsfirst     in/out: integer(nbmax)
c                 ihsfirst(i) is the location in isource
c                 array for the first hung source in box i
c                 
c    ihslast      in/out: integer(nbmax)
c                 ihslast(i) is the location in isource
c                 array for the last hung source in box i
c                 
c    isfirst      in/out: integer(nbmax)
c                 isfirst(i) is the location in isource
c                 array for the first source in box i
c                 
c    islast       in/out: integer(nbmax)
c                 islast(i) is the location in isource
c                 array for the last source in box i
c                 
c    itfirst      in/out: integer(nbmax)
c                 itfirst(i) is the location in itarg
c                 array for the first target in box i
c                 
c    itlast       in/out: integer(nbmax)
c                 itlast(i) is the location in itarg
c                 array for the last target in box i
c                 
c    ihefirst     in/out: integer(nbmax)
c                 ihefirst(i) is the location in iexpc
c                 array for the first hung expansion center in box i
c                 
c    ihelast      in/out: integer(nbmax)
c                 ihelast(i) is the location in iexpc
c                 array for the last hung expansion center in box i
c                 
c    iefirst      in/out: integer(nbmax)
c                 iefirst(i) is the location in iexpc
c                 array for the first expansion center in box i
c                 
c    ielast       in/out: integer(nbmax)
c                 ielast(i) is the location in iexpc
c                 array for the last expansion center in box i
c                
c    nhungsrc     in/out: integer(nbmax)
c                 nhungsrc(i) is the number of hung sources in box i
c
c    nhungexp     in/out: integer(nbmax)
c                 nhungexp(i) is the number of hung 
c                 expansion centers in box i

      implicit none
c     Calling sequence variable declaration
      integer ier
      integer ns,nt,nexpc
      double precision src(3,ns), trg(3,nt), expc(3,nexpc)
      double precision radsrc(ns),radexp(nexpc)

      integer nlevels,nboxes,nbmax, mnbors
      double precision boxsize(0:nlevels), centers(3,nbmax)

      integer laddr(2,0:nlevels),ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax), ichild(8,nbmax)
      integer nnbors(nbmax), nbors(mnbors,nbmax)
      integer isource(ns),itarget(nt),iexpc(nexpc)

      integer ihsfirst(nbmax), ihslast(nbmax)
      integer isfirst(nbmax), islast(nbmax)

      integer itfirst(nbmax), itlast(nbmax)

      integer ihefirst(nbmax), ihelast(nbmax)
      integer iefirst(nbmax), ielast(nbmax)
      
      integer nhungsrc(nbmax), nhungexp(nbmax)

      integer, allocatable :: iflag(:)
c     Temporary variables
      integer i,j,k,l
      integer ibox,jbox,kbox,ilev
      integer idad,igranddad
      double precision xdis,ydis,zdis,distest

      integer laddrtail(2,0:nlevels),ict


      allocate(iflag(nbmax))

c     Initialize flag array
      do i=1,nboxes
         iflag(i) = 0
      enddo

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

      do ilev = nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
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

         call subdivide_flag(ier,src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
         if(ier.ne.0) return
         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,ihsfirst,
     2                    ihslast,isfirst,islast,itfirst,itlast,
     3                    ihefirst,ihelast,iefirst,ielast,
     4                    nhungsrc,nhungexp,iflag)

c     Compute colleague information again      

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,mnbors,nnbors,nbors)

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

      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
 
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
     1                    mnbors,nnbors,nbors,centers,boxsize,iflag)

         call updateflags(ilev,nboxes,nlevels,laddrtail,nchild,ichild,
     1                    mnbors,nnbors,nbors,centers,boxsize,iflag)
         
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1
         call subdivide_flag(ier,src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
         if(ier.ne.0) return

         call subdivide_flag(ier,src,ns,radsrc,trg,nt,
     $                   expc,nexpc,radexp,
     $                   ilev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddrtail,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
          if(ier.ne.0) return
          laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes
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
     4                   abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*boxsize(ilev+1)) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
      enddo

c     Reorganize tree once again and we are all done      
      call reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,
     2                    ihsfirst,
     3                    ihslast,isfirst,islast,itfirst,itlast,
     4                    ihefirst,ihelast,iefirst,ielast,
     5                    nhungsrc,nhungexp,iflag)

c     Compute colleague information again      

      do i=1,nboxes
         nnbors(i) = 0
         do j=1,mnbors
            nbors(j,i) = -1
         enddo
      enddo
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,mnbors,nnbors,nbors)
      
      return
      end
c---------------------------------------------------------------
      subroutine subdivide_flag(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,
     $                   curlev,nboxes,
     $                   centers,boxsize,nbmax,nlevels,
     $                   laddr,ilevel,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,iflag)
      implicit none
      integer ns,nt,nexpc,ier
      integer nlevels,nboxes,nbmax, curlev
      double precision src(3,ns),radsrc(ns)
      double precision trg(3,nt)
      double precision expc(3,nexpc),radexp(nexpc)
      double precision centers(3,nbmax)
      double precision boxsize(0:nlevels)
      integer laddr(2,0:nlevels)
      integer ilevel(nbmax)
      integer iparent(nbmax)
      integer nchild(nbmax)
      integer ichild(8,nbmax)
      integer isource(ns)
      integer itarget(nt)
      integer iexpc(nexpc)
      integer ihsfirst(nbmax)
      integer ihslast(nbmax)
      integer isfirst(nbmax)
      integer islast(nbmax)
      integer itfirst(nbmax)
      integer itlast(nbmax)
      integer ihefirst(nbmax)
      integer ihelast(nbmax)
      integer iefirst(nbmax)
      integer ielast(nbmax)
      integer nhungsrc(nbmax)
      integer nhungexp(nbmax)
      integer iflag(nbmax)
c     Temporary variables
      integer isrctmp(ns),itargtmp(nt),iexpctmp(nexpc)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee,ntt
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer i56, i78, i1234, i5678
      integer is,it,ie, jj
      integer nsc(8),ntc(8),nh(8),nexpcc(8),nhc(8)
c
c     for every box at level nlevels,
c     sort into children, updating various arrays 
c     perhaps just build tree here paren/child/particle sorting...
c     lists in second call ???
c     
c     allocate temp array for isourcetemp2 itargtemp2
c     after all done, write back to isourcetemp, itargtemp
c     this is O(N) * nlevels work for rewriting.
c     can be fancier I suppose.
c     
      ifirstbox = laddr(1,curlev)
      ilastbox =  laddr(2,curlev)

c
      do ibox = ifirstbox,ilastbox
c        The current box needs to be subdivided if iflag(ibox).gt.0       
         if(iflag(ibox).gt.0) then
c           Based on flagging criterion, the current
c           box needs to be divided. 

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

            i1234 = isfirst(ibox)-1
            i5678 = 0
            do is = isfirst(ibox),islast(ibox)
               if(src(3,isource(is)) - centers(3,ibox).lt.0) then
                  i1234 = i1234+1
                  isource(i1234) = isource(is)
               else
                  i5678 = i5678 + 1
                  isrctmp(i5678) = isource(is)
               endif
            enddo
c           Note at the end of the loop, i1234 is where the particles
c           in part 1234 of the box end

c           Reorder sources to include sources in 5678 in the array
            do i=1,i5678
               isource(i1234+i) = isrctmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = isfirst(ibox)-1
            i34 = 0
            do is = isfirst(ibox),i1234
               if(src(2,isource(is))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  isource(i12) = isource(is)
               else
                  i34 = i34 + 1
                  isrctmp(i34) = isource(is)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder sources to include 34 in the array
            do i=1,i34
               isource(i12+i) = isrctmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do is=i1234+1,islast(ibox)
               if(src(2,isource(is))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  isource(i56) = isource(is)
               else
                  i78 = i78 + 1
                  isrctmp(i78) = isource(is)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               isource(i56+i) = isrctmp(i)
            enddo
c           End of reordering i5678         

            nsc(1) = 0
            nsc(2) = 0
            nsc(3) = 0
            nsc(4) = 0
            nsc(5) = 0
            nsc(6) = 0
            nsc(7) = 0
            nsc(8) = 0
c           Sort into boxes 1 and 2
            do is = isfirst(ibox),i12
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(isfirst(ibox)+nsc(1)) = isource(is)
                  nsc(1) = nsc(1) + 1
               else
                  nsc(2) = nsc(2) + 1
                  isrctmp(nsc(2)) = isource(is)
               endif
            enddo
c           Reorder sources so that sources in 2 are at the
c           end of this part of the array
            do i=1,nsc(2)
               isource(isfirst(ibox)+nsc(1)+i-1) = isrctmp(i)
            enddo

c           Sort into boxes 3 and 4
            do is = i12+1, i1234
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(i12+1+nsc(3)) = isource(is)
                  nsc(3) = nsc(3) + 1
                else
                   nsc(4) = nsc(4)+1
                   isrctmp(nsc(4)) = isource(is)
                endif
            enddo
c           Reorder sources so that sources in 4 are at the
c           end of this part of the array
            do i=1,nsc(4)
               isource(i12+nsc(3)+i) = isrctmp(i)
            enddo

c           Sort into boxes 5 and 6
            do is = i1234+1,i56
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(i1234+1+nsc(5)) = isource(is)
                  nsc(5) = nsc(5) + 1
               else
                  nsc(6) = nsc(6) + 1
                  isrctmp(nsc(6)) = isource(is)
               endif
            enddo
c           Reorder sources so that sources in 6 are at the
c           end of this part of the array
            do i=1,nsc(6)
               isource(i1234+nsc(5)+i) = isrctmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do is=i56+1,islast(ibox)
               if(src(1,isource(is))-centers(1,ibox).lt.0) then
                  isource(i56+1+nsc(7)) = isource(is)
                  nsc(7) = nsc(7) + 1
               else
                  nsc(8) = nsc(8) + 1
                  isrctmp(nsc(8)) = isource(is)
               endif
            enddo
c           Reorder sources so that sources in 8 are at the
c           end of the array
            do i=1,nsc(8)
               isource(i56+nsc(7)+i) = isrctmp(i)
            enddo

            istart = isfirst(ibox)-1
            do j=1,8
c           check hung -> counter nh1,nh2,nh3,nh4,nh5,nh6,nh7,nh8
               ii = 0
               nh(j) = 0
               do i=1,nsc(j)
                  if(radsrc(isource(istart+i)).gt.boxsize(curlev+1))
     1            then     
                     nh(j) = nh(j) + 1
                     isource(istart+nh(j)) = isource(istart+i)
                  else
                     ii = ii+1
                     isrctmp(ii) = isource(istart+i)
                  endif
               enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
                do i=1,ii
                  isource(istart+nh(j)+i) = isrctmp(i)
                enddo
                istart = istart + nsc(j)
            enddo
c           End of sorting sources

c           Sort targets
c           which child?  1,2,3,4,5,6,7,8? counter nt1,nt2,nt3,nt4,
c           nt5,nt6,nt7,nt8
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
            i1234 = itfirst(ibox)-1
            i5678 = 0
            do it = itfirst(ibox),itlast(ibox)
               if(trg(3,itarget(it)) - centers(3,ibox).lt.0) then
                 i1234 = i1234+1
                 itarget(i1234) = itarget(it)
               else
                  i5678 = i5678 + 1
                  itargtmp(i5678) = itarget(it)
               endif
            enddo
c           Reorder sources to include targets in 5678 in the array
            do i=1,i5678
               itarget(i1234+i) = itargtmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = itfirst(ibox)-1
            i34 = 0
            do it = itfirst(ibox),i1234
               if(trg(2,itarget(it))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  itarget(i12) = itarget(it)
               else
                  i34 = i34 + 1
                  itargtmp(i34) = itarget(it)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder targets to include 34 in the array
            do i=1,i34
               itarget(i12+i) = itargtmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do it=i1234+1,itlast(ibox)
               if(trg(2,itarget(it))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  itarget(i56) = itarget(it)
               else
                  i78 = i78 + 1
                  itargtmp(i78) = itarget(it)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               itarget(i56+i) = itargtmp(i)
            enddo
c           End of reordering i5678         

            ntc(1) = 0
            ntc(2) = 0
            ntc(3) = 0
            ntc(4) = 0
            ntc(5) = 0
            ntc(6) = 0
            ntc(7) = 0
            ntc(8) = 0

c           Sort into boxes 1 and 2
            do it = itfirst(ibox),i12
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(itfirst(ibox)+ntc(1)) = itarget(it)
                  ntc(1) = ntc(1) + 1
               else
                  ntc(2) = ntc(2) + 1
                  itargtmp(ntc(2)) = itarget(it)
               endif
            enddo
c           Reorder targets so that sources in 2 are at the
c           end of this part of the array
            do i=1,ntc(2)
               itarget(itfirst(ibox)+ntc(1)+i-1) = itargtmp(i)
            enddo
c           Sort into boxes 3 and 4
            do it = i12+1, i1234
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(i12+1+ntc(3)) = itarget(it)
                  ntc(3) = ntc(3) + 1
                else
                   ntc(4) = ntc(4)+1
                   itargtmp(ntc(4)) = itarget(it)
                endif
            enddo
c           Reorder targets so that sources in 4 are at the
c           end of this part of the array
            do i=1,ntc(4)
               itarget(i12+ntc(3)+i) = itargtmp(i)
            enddo

c           Sort into boxes 5 and 6
            do it = i1234+1,i56
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(i1234+1+ntc(5)) = itarget(it)
                  ntc(5) = ntc(5) + 1
               else
                  ntc(6) = ntc(6) + 1
                  itargtmp(ntc(6)) = itarget(it)
               endif
            enddo
c           Reorder targets so that sources in 6 are at the
c           end of this part of the array
            do i=1,ntc(6)
               itarget(i1234+ntc(5)+i) = itargtmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do it=i56+1,itlast(ibox)
               if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                  itarget(i56+1+ntc(7)) = itarget(it)
                  ntc(7) = ntc(7) + 1
               else
                  ntc(8) = ntc(8) + 1
                  itargtmp(ntc(8)) = itarget(it)
               endif
            enddo
c           Reorder targets so that sources in 8 are at the
c           end of the array
            do i=1,ntc(8)
               itarget(i56+ntc(7)+i) = itargtmp(i)
            enddo
c           End of sorting targets

c           Sort expansion centers
c           which child?  1,2,3,4,5,6,7,8? counter
c           nexpcc1, nexpcc2, nexpcc3, nexpcc4,
c           nexpcc5, nexpcc6, nexpcc7, nexpcc8,
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
            i1234 = iefirst(ibox)-1
            i5678 = 0
            do ie = iefirst(ibox),ielast(ibox)
               if(expc(3,iexpc(ie)) - centers(3,ibox).lt.0) then
                  i1234 = i1234+1
                  iexpc(i1234) = iexpc(ie)
               else
                  i5678 = i5678 + 1
                  iexpctmp(i5678) = iexpc(ie)
               endif
            enddo

c           Reorder sources to include sources in 5678 in the array
            do i=1,i5678
               iexpc(i1234+i) = iexpctmp(i)
            enddo

c           Sort i1234 into i12 and i34         
            i12 = iefirst(ibox)-1
            i34 = 0
            do ie = iefirst(ibox),i1234
               if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
                  i12 = i12 + 1
                  iexpc(i12) = iexpc(ie)
               else
                  i34 = i34 + 1
                  iexpctmp(i34) = iexpc(ie)
               endif
            enddo
c           Note at the end of the loop, i12 is where the particles
c           in part 12 of the box end
c
c           Reorder sources to include 34 in the array
            do i=1,i34
               iexpc(i12+i) = iexpctmp(i)
            enddo

c           sort i5678 into i56 and i78
            i56 = i1234
            i78 = 0
            do ie=i1234+1,ielast(ibox)
               if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
                  i56 = i56 + 1
                  iexpc(i56) = iexpc(ie)
               else
                  i78 = i78 + 1
                  iexpctmp(i78) = iexpc(ie)
               endif
            enddo

c           Reorder sources to include 78 in the array
            do i=1,i78
               iexpc(i56+i) = iexpctmp(i)
            enddo
c           End of reordering i5678         

            nexpcc(1) = 0
            nexpcc(2) = 0
            nexpcc(3) = 0
            nexpcc(4) = 0
            nexpcc(5) = 0
            nexpcc(6) = 0
            nexpcc(7) = 0
            nexpcc(8) = 0

c           Sort into boxes 1 and 2
            do ie = iefirst(ibox),i12
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(iefirst(ibox)+nexpcc(1)) = iexpc(ie)
                  nexpcc(1) = nexpcc(1) + 1
               else
                  nexpcc(2) = nexpcc(2) + 1
                  iexpctmp(nexpcc(2)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 2 are at the
c           end of this part of the array
            do i=1,nexpcc(2)
               iexpc(iefirst(ibox)+nexpcc(1)+i-1) = iexpctmp(i)
            enddo
c           Sort into boxes 3 and 4
            do ie = i12+1, i1234
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i12+1+nexpcc(3)) = iexpc(ie)
                  nexpcc(3) = nexpcc(3) + 1
               else
                  nexpcc(4) = nexpcc(4)+1
                  iexpctmp(nexpcc(4)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 4 are at the
c           end of this part of the array
            do i=1,nexpcc(4)
               iexpc(i12+nexpcc(3)+i) = iexpctmp(i)
            enddo

c           Sort into boxes 5 and 6
            do ie = i1234+1,i56
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i1234+1+nexpcc(5)) = iexpc(ie)
                  nexpcc(5) = nexpcc(5) + 1
               else
                  nexpcc(6) = nexpcc(6) + 1
                  iexpctmp(nexpcc(6)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 6 are at the
c           end of this part of the array
            do i=1,nexpcc(6)
               iexpc(i1234+nexpcc(5)+i) = iexpctmp(i)
            enddo
c           End of sorting sources into boxes 5 and 6

c           Sort into boxes 7 and 8
            do ie=i56+1,ielast(ibox)
               if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                  iexpc(i56+1+nexpcc(7)) = iexpc(ie)
                  nexpcc(7) = nexpcc(7) + 1
               else
                  nexpcc(8) = nexpcc(8) + 1
                  iexpctmp(nexpcc(8)) = iexpc(ie)
               endif
            enddo
c           Reorder expc so that sources in 8 are at the
c           end of the array
            do i=1,nexpcc(8)
               iexpc(i56+nexpcc(7)+i) = iexpctmp(i)
            enddo
c           End of sorting expanison centers

            istart = iefirst(ibox)-1
            do j=1,8
c           check hung -> counter nhc1,nhc2,nhc3,nhc4,nhc5,nhc6,nhc7,nhc8
               ii = 0
               nhc(j) = 0
               do i=1,nexpcc(j)
                  if(radexp(iexpc(istart+i)).gt.boxsize(curlev+1))
     1            then     
                     nhc(j) = nhc(j) + 1
                     iexpc(istart+nhc(j)) = iexpc(istart+i)
                  else
                     ii = ii+1
                     iexpctmp(ii) = iexpc(istart+i)
                  endif
               enddo
c           Reorder sources to have hung chunks at the star
c           of the sorted sources in the box ibox
               do i=1,ii
                  iexpc(istart+nhc(j)+i) = iexpctmp(i)
               enddo
               istart = istart + nexpcc(j)
            enddo

            nchild(ibox) = 0
c           Create the required boxes
            istart = isfirst(ibox)
            jstart = itfirst(ibox)
            kstart = iefirst(ibox)
            do i=1,8
               ii = 2
               jj = 2
               if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
               if(i.lt.5) jj = 1
               if(nsc(i)+ntc(i)+nexpcc(i).ge.0) then
c                 Increment total number of boxes               
                  nboxes = nboxes + 1
                  if(nboxes.gt.nbmax) then
                    write(*,*) "Exceeding max number of boxes"
                    write(*,*) "Exiting"
                    ier = 12
                    return
                  endif

c                 Increment number of children for the current box
                  nchild(ibox) = nchild(ibox)+1
c                 Update the array of children for the current box
                  ichild(i,ibox) = nboxes
c                 Update the array of levels for the child box
                  ilevel(nboxes) = curlev+1
c                 Update the array of parents for the child box
                  iparent(nboxes) = ibox
c                 Compute center for the child box
                  centers(1,nboxes) = centers(1,ibox)+(-1)**i*
     1                                boxsize(curlev+1)/2.0
                  centers(2,nboxes) = centers(2,ibox)+(-1)**ii*
     1                               boxsize(curlev+1)/2.0
                  centers(3,nboxes) = centers(3,ibox)+(-1)**jj*
     1                               boxsize(curlev+1)/2.0

                  nchild(nboxes) = 0
                  ichild(1,nboxes) = -1
                  ichild(2,nboxes) = -1
                  ichild(3,nboxes) = -1
                  ichild(4,nboxes) = -1
                  ichild(5,nboxes) = -1
                  ichild(6,nboxes) = -1
                  ichild(7,nboxes) = -1
                  ichild(8,nboxes) = -1
c                 Update arrays ihsfirst,ihslast,isfirst,islast
                  ihsfirst(nboxes) = istart
                  ihslast(nboxes) = istart + nh(i) - 1
                  nhungsrc(nboxes) = nh(i)

                  isfirst(nboxes) = istart + nh(i)
                  islast(nboxes) = istart + nsc(i) - 1

c                 Update arrays itfirst, itlast
                  itfirst(nboxes) = jstart
                  itlast(nboxes) = jstart + ntc(i) - 1

c                 Update arrays ihefirst,ihelast,iefirst,ielast
                  ihefirst(nboxes) = kstart
                  ihelast(nboxes) = kstart + nhc(i)-1
                  nhungexp(nboxes) = nhc(i)

                  iefirst(nboxes) = kstart + nhc(i)
                  ielast(nboxes) = kstart + nexpcc(i) - 1

                  if(iflag(ibox).eq.1) iflag(nboxes) = 3
                  if(iflag(ibox).eq.2) iflag(nboxes) = 0
               endif
               istart = istart + nsc(i)
               jstart = jstart + ntc(i)
               kstart = kstart + nexpcc(i)
            enddo
         endif
      enddo

      return
      end
c-------------------------------------------------------------      
      subroutine reorganizetree(nboxes,centers,nlevels,laddr,laddrtail,
     1                    ilevel,iparent,nchild,ichild,ihsfirst,
     2                    ihslast,isfirst,islast,itfirst,itlast,
     3                    ihefirst,ihelast,iefirst,ielast,
     4                    nhungsrc,nhungexp,iflag)

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
c     ilevel      in/out: integer(nboxes)
c                 ilevel(i) is the level of box i
c
c     iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nchild      in/out: integer(nboxes)
c                 nchild(i) is the number of children 
c                 of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    ihsfirst     in/out: integer(nboxes)
c                 ihsfirst(i) is the location in isource
c                 array for the first hung source in box i
c                 
c    ihslast      in/out: integer(nboxes)
c                 ihslast(i) is the location in isource
c                 array for the last hung source in box i
c                 
c    isfirst      in/out: integer(nboxes)
c                 isfirst(i) is the location in isource
c                 array for the first source in box i
c                 
c    islast       in/out: integer(nboxes)
c                 islast(i) is the location in isource
c                 array for the last source in box i
c                 
c    itfirst      in/out: integer(nboxes)
c                 itfirst(i) is the location in itarg
c                 array for the first target in box i
c                 
c    itlast       in/out: integer(nboxes)
c                 itlast(i) is the location in itarg
c                 array for the last target in box i
c                 
c    ihefirst     in/out: integer(nboxes)
c                 ihefirst(i) is the location in iexpc
c                 array for the first hung expansion center in box i
c                 
c    ihelast      in/out: integer(nboxes)
c                 ihelast(i) is the location in iexpc
c                 array for the last hung expansion center in box i
c                 
c    iefirst      in/out: integer(nboxes)
c                 iefirst(i) is the location in iexpc
c                 array for the first expansion center in box i
c                 
c    ielast       in/out: integer(nboxes)
c                 ielast(i) is the location in iexpc
c                 array for the last expansion center in box i
c                
c    nhungsrc     in/out: integer(nboxes)
c                 nhungsrc(i) is the number of hung sources in box i
c
c    nhungexp     in/out: integer(nboxes)
c                 nhungexp(i) is the number of hung 
c                 expansion centers in box i
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
      integer ihsfirst(nboxes)
      integer ihslast(nboxes)
      integer isfirst(nboxes)
      integer islast(nboxes)
      integer itfirst(nboxes)
      integer itlast(nboxes)
      integer ihefirst(nboxes)
      integer ihelast(nboxes)
      integer iefirst(nboxes)
      integer ielast(nboxes)
      integer nhungsrc(nboxes)
      integer nhungexp(nboxes)
      integer iflag(nboxes)



      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:),tihsfirst(:),tihslast(:)
      integer, allocatable :: tisfirst(:),tislast(:),titfirst(:)
      integer, allocatable :: titlast(:),tihefirst(:),tihelast(:)
      integer, allocatable :: tiefirst(:),tielast(:),tnhungsrc(:)
      integer, allocatable :: tnhungexp(:),tiflag(:),iboxtocurbox(:)

      double precision, allocatable :: tcenters(:,:)

c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox

      allocate(tcenters(3,nboxes),tilevel(nboxes),tiparent(nboxes))
      allocate(tnchild(nboxes),tichild(8,nboxes),tihsfirst(nboxes))
      allocate(tihslast(nboxes),tisfirst(nboxes),tislast(nboxes))
      allocate(titfirst(nboxes),titlast(nboxes),tihefirst(nboxes))
      allocate(tihelast(nboxes),tiefirst(nboxes),tielast(nboxes))
      allocate(tnhungsrc(nboxes),tnhungexp(nboxes),tiflag(nboxes))
      allocate(iboxtocurbox(nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo

      do ibox=1,nboxes
         tilevel(ibox) = ilevel(ibox)
         tcenters(1,ibox) = centers(1,ibox)
         tcenters(2,ibox) = centers(2,ibox)
         tcenters(3,ibox) = centers(3,ibox)
         tiparent(ibox) = iparent(ibox)
         tnchild(ibox) = nchild(ibox)
         tichild(1,ibox) = ichild(1,ibox)
         tichild(2,ibox) = ichild(2,ibox)
         tichild(3,ibox) = ichild(3,ibox)
         tichild(4,ibox) = ichild(4,ibox)
         tichild(5,ibox) = ichild(5,ibox)
         tichild(6,ibox) = ichild(6,ibox)
         tichild(7,ibox) = ichild(7,ibox)
         tichild(8,ibox) = ichild(8,ibox)
         tihsfirst(ibox) = ihsfirst(ibox)
         tihslast(ibox) = ihslast(ibox)
         tisfirst(ibox) = isfirst(ibox)
         tislast(ibox) = islast(ibox)
         titfirst(ibox) = itfirst(ibox)
         titlast(ibox) = itlast(ibox)
         tihefirst(ibox) = ihefirst(ibox)
         tihelast(ibox) = ihelast(ibox)
         tiefirst(ibox) = iefirst(ibox)
         tielast(ibox) = ielast(ibox)
         tnhungsrc(ibox) = nhungsrc(ibox)
         tnhungexp(ibox) = nhungexp(ibox)
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
            iboxtocurbox(ibox) = ibox
         enddo
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
            ihsfirst(curbox) = tihsfirst(ibox)
            ihslast(curbox) = tihslast(ibox)
            isfirst(curbox) = tisfirst(ibox)
            islast(curbox) = tislast(ibox)
            itfirst(curbox) = titfirst(ibox)
            itlast(curbox) = titlast(ibox)
            ihefirst(curbox) = tihefirst(ibox)
            ihelast(curbox) = tihelast(ibox)
            iefirst(curbox) = tiefirst(ibox)
            ielast(curbox) = tielast(ibox)
            nhungsrc(curbox) = tnhungsrc(ibox)
            nhungexp(curbox) = tnhungexp(ibox)
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
            ihsfirst(curbox) = tihsfirst(ibox)
            ihslast(curbox) = tihslast(ibox)
            isfirst(curbox) = tisfirst(ibox)
            islast(curbox) = tislast(ibox)
            itfirst(curbox) = titfirst(ibox)
            itlast(curbox) = titlast(ibox)
            ihefirst(curbox) = tihefirst(ibox)
            ihelast(curbox) = tihelast(ibox)
            iefirst(curbox) = tiefirst(ibox)
            ielast(curbox) = tielast(ibox)
            nhungsrc(curbox) = tnhungsrc(ibox)
            nhungexp(curbox) = tnhungexp(ibox)
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
c-------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    mnbors,nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      mnbors         in: integer
c                     max number of neighbors
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(mnbors,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(3,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels, mnbors
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(8,nboxes)
      integer nnbors(nboxes), nbors(mnbors,nboxes)
      integer iflag(nboxes)
      double precision centers(3,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level      
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        zdis = centers(3,kbox) - centers(3,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(abs(zdis).le.distest) ict = ict + 1
                        if(ict.eq.3) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      

      return
      end
c-----------------------------------------------------------------
      subroutine computemnlists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,mnlist1,
     3                   mnlist2,mnlist3,mnlist4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes),nlist2(nboxes),nlist3(nboxes)
      integer nlist4(nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
      enddo
C$OMP END PARALLEL DO
      
      nlist1(1) = 1
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,distest,xdis,ydis,zdis)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                  1.05d0*isep*boxsize(ilev)).or.
     2                  (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                  1.05d0*isep*boxsize(ilev)).or.
     4                  (abs(centers(3,kbox)-centers(3,ibox)).ge.
     5                  1.05d0*isep*boxsize(ilev))) then

                        nlist2(ibox) = nlist2(ibox) + 1
                     endif
                  endif
               enddo
            enddo
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                     check for list1 at the same level
c
                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                  endif
c
cc                     check for list1 and list3 at one ilev+1
                  if(nchild(jbox).gt.0) then
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,8
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           zdis = dabs(centers(3,kbox)-centers(3,ibox))

                           if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                        zdis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                           else
                              nlist3(ibox) = nlist3(ibox)+1
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 and list4 for boxes at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                         2.0d0*isep
                      xdis = dabs(centers(1,jbox)-centers(1,ibox))
                      ydis = dabs(centers(2,jbox)-centers(2,ibox))
                      zdis = dabs(centers(3,jbox)-centers(3,ibox))
                      if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                  zdis.lt.distest) then
                         nlist1(ibox) = nlist1(ibox)+1
                      endif
                   endif
               enddo
            endif
c
cc           compute list 4 at level ilev-1
c
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                 2.0d0*isep
                    xdis = dabs(centers(1,jbox)-centers(1,ibox))
                    ydis = dabs(centers(2,jbox)-centers(2,ibox))
                    zdis = dabs(centers(3,jbox)-centers(3,ibox))
                    if(xdis.gt.distest.or.ydis.gt.distest.or.
     1                 zdis.gt.distest) then
                       nlist4(ibox) = nlist4(ibox)+1
                    endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
C$OMP$REDUCTION(max:mnlist1,mnlist2,mnlist3,mnlist4)      
      do i=1,nboxes
         if(nlist1(i).gt.mnlist1) mnlist1 = nlist1(i)
         if(nlist2(i).gt.mnlist2) mnlist2 = nlist2(i)
         if(nlist3(i).gt.mnlist3) mnlist3 = nlist3(i)
         if(nlist4(i).gt.mnlist4) mnlist4 = nlist4(i)
      enddo
C$OMP END PARALLEL DO      

      return
      end
c---------------------------------------------------------------      
      subroutine computelists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes),nlist2(nboxes),nlist3(nboxes)
      integer nlist4(nboxes)
      integer list1(mnlist1,nboxes),list2(mnlist2,nboxes)
      integer list3(mnlist3,nboxes),list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
      enddo
C$OMP END PARALLEL DO      
      if(nchild(1).eq.0) then
         nlist1(1) = 1
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,xdis,ydis,zdis,distest)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                  1.05d0*isep*boxsize(ilev)).or.
     2                  (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                  1.05d0*isep*boxsize(ilev)).or.
     4                  (abs(centers(3,kbox)-centers(3,ibox)).ge.
     5                  1.05d0*isep*boxsize(ilev))) then

                        nlist2(ibox) = nlist2(ibox) + 1
                        list2(nlist2(ibox),ibox) = kbox
                     endif
                  endif
               enddo
            enddo
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                boxes in list 1 at the same level
c

                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                     list1(nlist1(ibox),ibox) = jbox
                  else
c
cc                     boxes in list1 and list3 at level ilev+1
c
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,8
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           zdis = dabs(centers(3,kbox)-centers(3,ibox))

                           if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                        zdis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                              list1(nlist1(ibox),ibox) = kbox
                           else
                              nlist3(ibox) = nlist3(ibox)+1
                              list3(nlist3(ibox),ibox) = kbox
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                         2.0d0*isep
                      xdis = dabs(centers(1,jbox)-centers(1,ibox))
                      ydis = dabs(centers(2,jbox)-centers(2,ibox))
                      zdis = dabs(centers(3,jbox)-centers(3,ibox))
                      if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                  zdis.lt.distest) then
                         nlist1(ibox) = nlist1(ibox)+1
                         list1(nlist1(ibox),ibox) = jbox
                      endif
                   endif
                enddo
            endif
c
cc           compute list 4 at level ilev-1
c
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                 2.0d0*isep
                    xdis = dabs(centers(1,jbox)-centers(1,ibox))
                    ydis = dabs(centers(2,jbox)-centers(2,ibox))
                    zdis = dabs(centers(3,jbox)-centers(3,ibox))
                    if(xdis.gt.distest.or.ydis.gt.distest.or.
     1                 zdis.gt.distest) then
                       nlist4(ibox) = nlist4(ibox)+1
                       list4(nlist4(ibox),ibox)=jbox
                    endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c----------------------------------------------------------------

      subroutine getpwlistall(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      double precision boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      double precision centers(3,nboxes)
      integer isep
      integer nuall,ndall,nnall,nsall,neall,nwall,nu1234
      integer nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer uall(1),dall(1),nall(1),sall(1),eall(1),wall(1)
      integer u1234(1),d5678(1),n1256(1),s3478(1),e1357(1),w2468(1)
      integer n12(1),n56(1),s34(1),s78(1),e13(1),e57(1),w24(1),w68(1)
      integer e1(1),e3(1),e5(1),e7(1)
      integer w2(1),w4(1),w6(1),w8(1)

      integer jbox,kbox
      integer c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c16
      integer i,j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      nu1234 = 0
      nd5678 = 0
      nn1256 = 0
      ns3478 = 0
      ne1357 = 0
      nw2468 = 0
      nn12 = 0
      nn56 = 0
      ns34 = 0
      ns78 = 0
      ne13 = 0
      ne57 = 0
      nw24 = 0
      nw68 = 0
      ne1 = 0
      ne3 = 0
      ne5 = 0
      ne7 = 0
      nw2 = 0
      nw4 = 0
      nw6 = 0
      nw8 = 0
      do i=1,nnbors
         jbox = nbors(i)
         do j=1,8
            kbox = ichild(j,jbox)
            if(kbox.gt.0) then
               c1 = 0
               c2 = 0
               c3 = 0
               c4 = 0
               c5 = 0
               c6 = 0
               c7 = 0
               c8 = 0
               c9 = 0
               c10 = 0
               c11 = 0
               c12 = 0
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c1 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c2 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c3 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c4 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c5 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c6 = 1
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c7 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c8 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c9 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c10 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c11 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c12 = 1
               if(c1.eq.1) then
                  nuall = nuall + 1
                  uall(nuall) = kbox
               endif

               if(c2.eq.1) then
                  ndall = ndall + 1
                  dall(ndall) = kbox
               endif

               if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
                  nnall = nnall + 1
                  nall(nnall) = kbox
               endif

               if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
                  nsall = nsall + 1
                  sall(nsall) = kbox
               endif

               if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1             c4.eq.0) then
                  neall = neall + 1
                  eall(neall) = kbox
               endif

               if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1            c4.eq.0) then
                  nwall = nwall + 1
                  wall(nwall) = kbox
               endif

               c16 = c1 + c2 + c3 + c4 + c5 +c6
               if(c16.eq.0.and.c7.eq.1) then
                  nu1234 = nu1234 + 1
                  u1234(nu1234) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1) then
                  nd5678 = nd5678 + 1
                  d5678(nd5678) = kbox
               endif

               if(c16.eq.0.and.c9.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  nn1256 = nn1256 + 1
                  n1256(nn1256) = kbox
               endif

               if(c16.eq.0.and.c10.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  ns3478 = ns3478 + 1
                  s3478(ns3478) = kbox
               endif

               if(c16.eq.0.and.c11.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  ne1357 = ne1357 + 1
                  e1357(ne1357) = kbox
               endif

               if(c16.eq.0.and.c12.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  nw2468 = nw2468 + 1
                  w2468(nw2468) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c9.eq.1) then
                  nn12 = nn12 + 1
                  n12(nn12) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c9.eq.1) then
                  nn56 = nn56 + 1
                  n56(nn56) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c10.eq.1) then
                  ns34 = ns34 + 1
                  s34(ns34) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c10.eq.1) then
                  ns78 = ns78 + 1
                  s78(ns78) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne13 = ne13 + 1
                  e13(ne13) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne57 = ne57 + 1
                  e57(ne57) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw24 = nw24 + 1
                  w24(nw24) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw68 = nw68 + 1
                  w68(nw68) = kbox
               endif

               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne1 = ne1 + 1
                  e1(ne1) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne3 = ne3 + 1
                  e3(ne3) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne5 = ne5 + 1
                  e5(ne5) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne7 = ne7 + 1
                  e7(ne7) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw2 = nw2 + 1
                  w2(nw2) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw4 = nw4 + 1
                  w4(nw4) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw6 = nw6 + 1
                  w6(nw6) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw8 = nw8 + 1
                  w8(nw8) = kbox
               endif
            endif
         enddo
      enddo
      

      return
      end
c--------------------------------------------------------------------      

      subroutine getlist3pwlistall(ibox,bs,nboxes,nlist3,list3,
     1           isep,centers,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      integer isep
      integer nboxes,nlist3,list3(nlist3)
      double precision centers(3,nboxes)
      double precision sepdist,bs
      integer nuall,ndall,nnall,nsall,neall,nwall
      integer uall(1),dall(1),nall(1),sall(1),eall(1),wall(1)

      integer jbox
      integer c1,c2,c3,c4,c5,c6
      integer j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      sepdist = 1.01d0*isep*bs+bs/2.0d0
      do j=1,nlist3
         jbox = list3(j)
C         print *,"jbox: ",jbox
         if(jbox.gt.0) then
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            c5 = 0
            c6 = 0
            if((centers(3,jbox)-centers(3,ibox)).ge.sepdist) c1 = 1
            if((centers(3,jbox)-centers(3,ibox)).le.-sepdist) c2 = 1
            if((centers(2,jbox)-centers(2,ibox)).ge.sepdist) c3 = 1
            if((centers(2,jbox)-centers(2,ibox)).le.-sepdist) c4 = 1
            if((centers(1,jbox)-centers(1,ibox)).ge.sepdist) c5 = 1
            if((centers(1,jbox)-centers(1,ibox)).le.-sepdist) c6 = 1

C            print *,c1,c2,c3,c4,c5,c6

            if(c1.eq.1) then
               nuall = nuall + 1
               uall(nuall) = jbox
            endif

            if(c2.eq.1) then
               ndall = ndall + 1
               dall(ndall) = jbox
            endif

            if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
               nnall = nnall + 1
               nall(nnall) = jbox
            endif

            if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
               nsall = nsall + 1
               sall(nsall) = jbox
            endif

            if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1         c4.eq.0) then
               neall = neall + 1
               eall(neall) = jbox
            endif

            if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1         c4.eq.0) then
               nwall = nwall + 1
               wall(nwall) = jbox
            endif
         endif
      enddo
      

      return
      end
c
c
c
c
c--------------------------------------------------------------
c
      subroutine subdividebox(pos,npts,center,boxsize,
     1           isorted,iboxfl,subcenters)
      implicit none
      double precision pos(3,npts)
      double precision center(3)
      double precision subcenters(3,8)
      double precision boxsize
      integer npts
      integer isorted(*)
      integer iboxfl(2,8)

c     Temporary variables
      integer isortedtmp(npts)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee
      integer jj,irefinebox,ntt
      integer i56, i78, i1234, i5678
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer is,it,ie
      integer nc(8)

      i1234 = 0
      i5678 = 0
      do i = 1,npts
         if(pos(3,i)-center(3).lt.0) then
            i1234 = i1234 + 1
            isorted(i1234) = i
         else
            i5678 = i5678 + 1
            isortedtmp(i5678) = i
         endif
      enddo
      do i=1,i5678
         isorted(i1234+i) = isortedtmp(i)
      enddo

c     Sort i1234 into i12 and i34         
      i12 = 0
      i34 = 0
      do i = 1,i1234
         if(pos(2,isorted(i))-center(2).lt.0) then
            i12 = i12 + 1
            isorted(i12) = isorted(i)
         else
            i34 = i34 + 1
            isortedtmp(i34) = isorted(i)
         endif
      enddo
c     Note at the end of the loop, i12 is where the particles
c     in part 12 of the box end
c
c     Reorder sources to include 34 in the array
      do i=1,i34
         isorted(i12+i) = isortedtmp(i)
      enddo

c     sort i5678 into i56 and i78
      i56 = i1234
      i78 = 0
      do i=i1234+1,npts
         if(pos(2,isorted(i))-center(2).lt.0) then
            i56 = i56 + 1
            isorted(i56) = isorted(i)
         else
            i78 = i78 + 1
            isortedtmp(i78) = isorted(i)
         endif
      enddo
c     Reorder sources to include 78 in the array
      do i=1,i78
         isorted(i56+i) = isortedtmp(i)
      enddo
c     End of reordering i5678         

      nc = 0
c     Sort into boxes 1 and 2
      do i = 1,i12
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(1) = nc(1) + 1
            isorted(nc(1)) = isorted(i)
         else
            nc(2) = nc(2) + 1
            isortedtmp(nc(2)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 2 are at the
c     end of this part of the array
      do i=1,nc(2)
         isorted(nc(1)+i) = isortedtmp(i)
      enddo

c     Sort into boxes 3 and 4
      do i = i12+1, i1234
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(3) = nc(3) + 1
            isorted(i12+nc(3)) = isorted(i)
         else
            nc(4) = nc(4)+1
            isortedtmp(nc(4)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 4 are at the
c     end of this part of the array
      do i=1,nc(4)
         isorted(i12+nc(3)+i) = isortedtmp(i)
      enddo

c     Sort into boxes 5 and 6
      do i = i1234+1,i56
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(5) = nc(5) + 1
            isorted(i1234+nc(5)) = isorted(i)
         else
            nc(6) = nc(6) + 1
            isortedtmp(nc(6)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 6 are at the
c     end of this part of the array
      do i=1,nc(6)
         isorted(i1234+nc(5)+i) = isortedtmp(i)
      enddo
c     End of sorting sources into boxes 5 and 6

c     Sort into boxes 7 and 8
      do i=i56+1,npts
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(7) = nc(7) + 1
            isorted(i56+nc(7)) = isorted(i)
         else
            nc(8) = nc(8) + 1
            isortedtmp(nc(8)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 8 are at the
c     end of the array
      do i=1,nc(8)
         isorted(i56+nc(7)+i) = isortedtmp(i)
      enddo
      
      istart=1
      iboxfl=0
      subcenters=0.0d0
      do i=1,8
         ii = 2
         jj = 2
         if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
         if(i.lt.5) jj = 1
         if(nc(i).gt.0) then
            subcenters(1,i) = center(1)+(-1)**i*boxsize/2.0d0
            subcenters(2,i) = center(2)+(-1)**ii*boxsize/2.0d0
            subcenters(3,i) = center(3)+(-1)**jj*boxsize/2.0d0

            iboxfl(1,i) = istart
            iboxfl(2,i) = istart+nc(i)-1

            istart = istart+nc(i)
          endif
      enddo

      return
      end
c------------------------------------------------------------------      
c--------------------------------------------------------------------     
c
c
c
c--------------------------------------------------------------------      

      subroutine getlist4pwdir(dir,censrc,centrg,boxsize)
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer dir
      double precision censrc(3)
      double precision centrg(3)
      double precision boxsize
ccc   scopded function variables
      double precision sepdist
      double precision ctmp(3)
      ctmp(1) = censrc(1)-0*boxsize/2.0d0
      ctmp(2) = censrc(2)-0*boxsize/2.0d0
      ctmp(3) = censrc(3)-0*boxsize/2.0d0

      sepdist=1.51d0*boxsize

      if(ctmp(3)-centrg(3).ge.sepdist) then
        dir=1
      else if(ctmp(3)-centrg(3).le.-sepdist) then
        dir=2
      else if(ctmp(2)-centrg(2).ge.sepdist) then
        dir=3
      else if(ctmp(2)-centrg(2).le.-sepdist) then
        dir=4
      else if(ctmp(1)-centrg(1).ge.sepdist) then
        dir=5
      else if(ctmp(1)-centrg(1).le.-sepdist) then
        dir=6
      else
        dir=0
      end if
C      dir=0

      return
      end
      
      subroutine getlist4pwdirtest(dir,censrc,centrg,boxsize)
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer dir
      double precision censrc(3)
      double precision centrg(3)
      double precision boxsize
ccc   scopded function variables
      double precision sepdist
      double precision ctmp(3)
      ctmp(1) = censrc(1)-0*boxsize/2.0d0
      ctmp(2) = censrc(2)-0*boxsize/2.0d0
      ctmp(3) = censrc(3)-0*boxsize/2.0d0

      sepdist=1.51d0*boxsize

      if((ctmp(3)-centrg(3)).ge.sepdist) then
        dir=1
      else if((ctmp(3)-centrg(3)).le.-sepdist) then
        dir=2
      else if((ctmp(2)-centrg(2)).ge.sepdist) then
        dir=3
      else if((ctmp(2)-centrg(2)).le.-sepdist) then
        dir=4
      else if((ctmp(1)-centrg(1)).ge.sepdist) then
        dir=5
      else if((ctmp(1)-centrg(1)).le.-sepdist) then
        dir=6
      else
        dir=0
        print *,"dir:",dir
      end if
C      dir=0

      return
      end
