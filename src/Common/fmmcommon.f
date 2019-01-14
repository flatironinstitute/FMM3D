
      subroutine mpalloc(nd,laddr,iaddr,nlevels,lmptot,nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     nd          in: Integer
c                 number of expansions per box            
c 
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array provinding access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd
      integer iaddr(2,1), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp

      istart = 1
      do i = 0,nlevels

        nn = (2*nterms(i)+1)*2*(nterms(i)+1)*nd 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)

         do ibox = laddr(1,i),laddr(2,i)
c            Allocate memory for the multipole expansion         

             itmp = ibox-laddr(1,i)
             iaddr(1,ibox) = istart + itmp*2*nn

c            Allocate memory for the local expansion
             iaddr(2,ibox) = istart + itmp*2*nn + nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*2*nn
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     

      subroutine dreorderf(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer ndim,idim,i,n
      double precision arr(ndim,1),arrsort(ndim,1)
      integer iarr(1)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,i) = arr(idim,iarr(i))
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c--------------------------------------------------
      subroutine dreorderi(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the inverse of the
c        sorting order defined by
c        iarr.
c
c        Note that this subroutine is the inverse of 
c        dreorderf
c
c        arrsort(j,iarr(i)) = arr(j,i), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer i,idim,ndim,n
      double precision arr(ndim,1),arrsort(ndim,1)
      integer iarr(1)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,iarr(i)) = arr(idim,i)
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c----------------------------------------------------------      


      subroutine mpzero(nd,mpole,nterms)
c
cc      this subroutine zeros out a collection of 
c       multipole expansions
c
c       input:
c       nd     in:integer
c              number of expansions
c
c       nterms in:integer
c              number of terms in the expansions
c
c       inout:
c       mpole  inout: double complex(nd,0:nterms,-nterms:nterms)
c              multipole expansions to be zeroed out
c
      implicit none
      integer nd, i,j,nterms,idim
      double complex mpole(nd,0:nterms,-nterms:nterms)

      do i=-nterms,nterms
        do j=0,nterms
          do idim=1,nd
            mpole(idim,j,i) = 0
          enddo
        enddo
      enddo

      return
      end

c----------------------------------------------------------------

      subroutine mpscale(nd,nterms,mpolein,rsc,mpoleout)
c
cc      this subroutine rescales a multipole 
c       expansion where mpoleout(i,j) = mpolein(i,j)*rscalepow(i)
c
c       input
c       nd      in: integer
c               number of multipole expansions
c
c       nterms  in: integer
c               number of terms in the expansion
c    
c       rsc     in: double precision(0:nterms)
c               scaling factor for the multipole expansions            
c
c       mpolein   in: double complex(nd,0:nterms,-nterms:nterms)
c                input multipole expansions
c
c       output
c       mpoleout  out: double complex(nd,0:nterms,-nterms:nterms)
c                output multipole expansions
c
      implicit none
      integer nd,nterms,i,j,idim
      double precision rsc(0:nterms)
      double complex mpolein(nd,0:nterms,-nterms:nterms)
      double complex mpoleout(nd,0:nterms,-nterms:nterms)

      do i=0,nterms
        do j=-nterms,nterms
          do idim=1,nd
            mpoleout(idim,i,j) = mpolein(idim,i,j)*rsc(i)
          enddo
        enddo
      enddo

      return
      end
