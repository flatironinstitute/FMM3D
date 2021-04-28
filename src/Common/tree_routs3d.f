c
c
c    common routines for generating and processing
c     a level restricted quad tree in 2D
c   
c
c
c
      subroutine tree_refine_boxes(irefinebox,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(3,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)
      integer ii

      integer i,ibox,nel0,j,l,jbox,nel1,nbl,jj

      allocate(isum(nbloc))
      if(nbloc.gt.0) call cumsum(nbloc,irefinebox,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,ii,jj,l)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*8
          
          nchild(ibox) = 8
          do j=1,8
            ii = 2
            jj = 2
            if(j.eq.1.or.j.eq.2.or.j.eq.5.or.j.eq.6) ii = 1
            if(j.lt.5) jj = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            centers(3,jbox) = centers(3,ibox)+(-1)**jj*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(nb,centers,ilevel,iparent,nchild,ichild,
     1              centers2,ilevel2,iparent2,nchild2,ichild2)

       implicit none
       integer nd,nb,npb
       real *8 centers(3,nb),centers2(3,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(8,nb),ichild2(8,nb)

       integer i,j,nel


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         centers2(3,i) = centers(3,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         ichild2(1,i) = ichild(1,i)
         ichild2(2,i) = ichild(2,i)
         ichild2(3,i) = ichild(3,i)
         ichild2(4,i) = ichild(4,i)
         ichild2(5,i) = ichild(5,i)
         ichild2(6,i) = ichild(6,i)
         ichild2(7,i) = ichild(7,i)
         ichild2(8,i) = ichild(8,i)
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,
     2                       iper,nnbors,nbors)

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
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(3,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(8,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(27,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes,iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes)
      integer nbors(27,nboxes)

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
     2                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     3                   1.05*boxsize(ilev))) then
                     
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
c
c
c
c
c--------------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

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
c      ichild         in: integer(8,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(27,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
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
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(8,nboxes)
      integer nnbors(nboxes), nbors(27,nboxes)
      integer iflag(nboxes)
      double precision centers(3,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
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
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c

      subroutine tree_refine_boxes_flag(iflag,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(3,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox
      integer, allocatable :: isum(:),itmp(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer ii,jj

      allocate(isum(nbloc),itmp(nbloc))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i)
      do i=1,nbloc
        ibox = ifirstbox+i-1
        itmp(i) = 0
        if(iflag(ibox).gt.0) itmp(i) = 1
      enddo
C$OMP END PARALLEL DO
      if(nbloc.gt.0) call cumsum(nbloc,itmp,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,ii,jj,l)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(iflag(ibox).gt.0) then
          nbl = nbctr + (isum(i)-1)*8
          
          nchild(ibox) = 8
          do j=1,8
            ii = 2
            jj = 2
            if(j.eq.1.or.j.eq.2.or.j.eq.5.or.j.eq.6) ii = 1
            if(j.lt.5) jj = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            centers(3,jbox) = centers(3,ibox)+(-1)**jj*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
      subroutine computemnlists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,iper,mnlist1,
     3                   mnlist2,mnlist3,mnlist4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer iper
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
     2                   ichild,isep,nnbors,mnbors,nbors,iper,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4,isep
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

      return
      end
c      
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
c--------------------------------------------------------------------     
c
c
c
