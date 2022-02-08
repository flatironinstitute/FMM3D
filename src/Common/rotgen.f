c     This file contains all the relevant rotation subroutines
c     for the 3D code.
c
c  
      subroutine rotztoy(nd,nterms,mpole,mrotate,rdminus)

c     Compute multipole expansions in the new frame of
c     reference given by
c     z_new <- y_old
c     y_new <- x_old
c     x_new <- z_old
c
c     This is achieved by first rotating by -pi/2 about
c     the z-axis in the original co-ordinate system
c     followed by a rotation of -pi/2 radians about the y'
c     axis in the rotated frame of refernce
c 
c     INPUT arguments:
c     nd           in: integer(8)
c                  number of expansions
c     nterms       in: integer(8)
c                  number of terms in the multipole expansion
c
c     mpole        in: double complex (nd,0:nterms,-nterms:nterms)
c                  The multipole expansion to be rotated
c
c     rdminus      in: double complex (0:nterms,0:nterms,-nterms:nterms)
c                  Rotation matrix for transforming multipole
c                  expansion corresponds to a rotation of -pi/2
c                  about the y axis
c
c     OUTPUT
c     mrotate      out: double complex (nd,0:nterms,-nterms:nterms)
c                  rotated multipole expansion  
c                   
      implicit none
      integer(8) nterms,nd
      double precision rdminus(0:nterms,0:nterms,-nterms:nterms),rr
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex, allocatable :: mptemp(:,:,:)
      double complex mrotate(nd,0:nterms,-nterms:nterms)
      double complex ephi(-nterms:nterms), eyem, eye
      integer(8) idim
      integer(8) l,m,mp

      allocate(mptemp(nd,0:nterms,-nterms:nterms))
      call mpzero(nd,mrotate,nterms) 

c     First rotate by -pi/2 radians about the z axis 
c     in the original coordinate system
      ephi(0) = 1.0d0
      eyem = dcmplx(0.0d0,-1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      do m=1,nterms
         ephi(m) = ephi(m-1)*eyem
         ephi(-m) = ephi(-m+1)*eye
      enddo

      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               mptemp(idim,l,m) = ephi(m)*mpole(idim,l,m)
            enddo
         enddo
      enddo

c     Now rotate by -pi/2 radians about the y' axis in
c     the new co-ordinate system

      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               mrotate(idim,l,m) = 0.0d0
            enddo
            do mp=-l,l
               if(mp.ge.0) then
                  do idim=1,nd
                     mrotate(idim,l,m) = mrotate(idim,l,m) +
     1               mptemp(idim,l,mp)*rdminus(l,mp,m)
                  enddo
               else
                  do idim=1,nd
                     mrotate(idim,l,m) = mrotate(idim,l,m) +
     1               mptemp(idim,l,mp)*rdminus(l,-mp,-m)
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end
c---------------------------------------------------------------  
      subroutine rotytoz(nd,nterms,mpole,mrotate,rdplus)

c     Compute multipole expansions in the new frame of
c     reference given by
c
c     This is achieved by first rotating by pi/2 about
c     the y'-axis in the rotated co-ordinate system
c     followed by a rotation of pi/2 radians about the z
c     axis to bring back to original frame of reference
c 
c     INPUT arguments:
c     nd           in: integer(8)
c                  number of expansions
c     nterms       in: integer(8)
c                  number of terms in the multipole expansion
c
c     mpole        in: double complex (nd,0:nterms,-nterms:nterms)
c                  The multipole expansion to be rotated
c
c     rdplus       in: double complex (0:nterms,0:nterms,-nterms:nterms)
c                  Rotation matrix for transforming multipole
c                  expansion corresponds to a rotation of -pi/2
c                  about the y axis
c
c     OUTPUT
c     mrotate      out: double complex (nd,0:nterms,-nterms:nterms)
c                  rotated multipole expansion  
c                   
      implicit none
      integer(8) nterms,nd
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex, allocatable :: mptemp(:,:,:)
      double complex mrotate(nd,0:nterms,-nterms:nterms)
      double complex ephi(-nterms:nterms), eyem, eye
      integer(8) l,m,mp,idim

      allocate(mptemp(nd,0:nterms,-nterms:nterms))

c     Now rotate by pi/2 radians about the y' axis in
c     the new co-ordinate system

      call mpzero(nd,mrotate,nterms) 
      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               mptemp(idim,l,m) = 0.0d0
            enddo
            do mp=-l,l
               if(mp.ge.0) then
                  do idim=1,nd
                     mptemp(idim,l,m) = mptemp(idim,l,m) +
     1               mpole(idim,l,mp)*rdplus(l,mp,m)
                  enddo
               else
                  do idim=1,nd
                     mptemp(idim,l,m) = mptemp(idim,l,m) +
     1              mpole(idim,l,mp)*rdplus(l,-mp,-m)
                  enddo
               endif
            enddo
         enddo
      enddo
c     First rotate by pi/2 radians about the z axis 
c     in the original coordinate system
      ephi(0) = 1.0d0
      eyem = dcmplx(0.0d0,-1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      do m=1,nterms
         ephi(m) = ephi(m-1)*eye
         ephi(-m) = ephi(-m+1)*eyem
      enddo

      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               mrotate(idim,l,m) = ephi(m)*mptemp(idim,l,m)
            enddo
         enddo
      enddo

      return
      end
c--------------------------------------------------------      
      subroutine rotztox(nd,nterms,mpole,mrotate,rdplus)

c     Compute multipole expansions in the new frame of
c     reference given by
c     z_new <- x_old
c     y_new <- y_old
c     -x_new <- z_old
c
c     This is achieved by first rotating by pi/2 about
c     the y-axis in the original co-ordinate system
c 
c     INPUT arguments:
c     nd           in: integer(8)
c                  number of expansions
c
c     nterms       in: integer(8)
c                  number of terms in the multipole expansion
c
c     mpole        in: double complex (nd,0:nterms,-nterms:nterms)
c                  The multipole expansion to be rotated
c
c     rdplus       in: double complex (0:nterms,0:nterms,-nterms:nterms)
c                  Rotation matrix for transforming multipole
c                  expansion corresponds to a rotation of pi/2
c                  about the y axis
c
c     OUTPUT
c     mrotate      out: double complex (nd,0:nterms,-nterms:nterms)
c                  rotated multipole expansion  
c                   
      implicit none
      integer(8) nterms,nd
      double precision rdplus(0:nterms,0:nterms,-nterms:nterms)
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex mrotate(nd,0:nterms,-nterms:nterms)
      integer(8) l,m,mp,idim

c     rotate by pi/2 radians about the y' axis in
c     the new co-ordinate system

      call mpzero(nd,mrotate,nterms) 
      do l=0,nterms
         do m=-l,l
            do idim=1,nd
              mrotate(idim,l,m) = 0.0d0
            enddo
            do mp=-l,l
               if(mp.ge.0) then
                  do idim=1,nd
                     mrotate(idim,l,m) = mrotate(idim,l,m) +
     1               mpole(idim,l,mp)*rdplus(l,mp,m)
                  enddo
               else
                  do idim=1,nd
                    mrotate(idim,l,m) = mrotate(idim,l,m) +
     1               mpole(idim,l,mp)*rdplus(l,-mp,-m)
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end
