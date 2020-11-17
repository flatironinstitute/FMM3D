c
c     This file contains the translation operators for the helmholtz 
c     3D FMM
c
c     h3dmpmp - shift center of multipole expansion
c     h3dmploc - convert multipole expansion to a local expansion
c     h3dlocloc - shift center of local expansion  
c
C***********************************************************************
      subroutine h3dmpmp(nd,wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,
     2           radius,xnodes,wts,nquad)
C***********************************************************************
C
C
C     Usage:
C
C           Shift center of multipole expansion mpole, and add output 
C           multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C           
C           NOTE: the output multipole expansion is not initialized
C              in this routine and the shifted 
C              multipole expansion is directly added to it
C              
C           
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nd     = number of multipole expansions
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for shifted expansion
C           xnynzn = center of shifted expansion
C           nterms2 = order of shifted expansion
C           radius  = radius of sphere on which mpole expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted mpole expansion
C
C
C***********************************************************************
      implicit none
c
cc      calling sequence variables
c
      integer(8)  nterms,nterms2,nquad,nd
      real *8 x0y0z0(3),xnynzn(3)
      real *8 radius,xnodes(nquad),wts(nquad)
      real *8 sc1,sc2
      complex *16 wavek
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      complex *16 mpolen(nd,0:nterms2,-nterms2:nterms2)

c
cc      temporary variables
c
      complex *16, allocatable :: marray(:,:,:),marray1(:,:,:)
      complex *16, allocatable :: mptemp(:,:,:),mp2(:,:,:)
      complex *16, allocatable :: ephi(:),phitemp(:,:,:)
      complex *16, allocatable :: phitempn(:,:,:)
      complex *16, allocatable :: fhs(:),fhder(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16 imag

      integer(8)  nq, ldc,idim
      real *8 rshift
      real *8 d,theta,ctheta,phi,rvec(3)
      integer(8)  l,m,jnew,knew


      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)

      allocate(marray(nd,0:ldc,-ldc:ldc))
      allocate(marray1(nd,0:nterms,-nterms:nterms))
      allocate(mptemp(nd,0:nterms2,-nterms2:nterms2))
      allocate(mp2(nd,0:nterms2,-nterms2:nterms2))
      allocate(ephi(-ldc-1:ldc+1))
      allocate(phitemp(nd,nq,-ldc:ldc))
      allocate(phitempn(nd,nq,-ldc:ldc))
      allocate(fhs(0:ldc),fhder(0:ldc))
      allocate(ynm(0:ldc,0:ldc),ynmd(0:ldc,0:ldc))

C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)

      call cart2polar(rvec,d,theta,phi)

      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c----- a rotation of THETA radians about the Yprime axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               marray1(idim,l,m)=mpole(idim,l,m)*ephi(m)
            enddo
         enddo
      enddo

      if( nterms .ge. 30 ) then 
         call rotviaproj(nd,theta,nterms,nterms,nterms,marray1,
     1         nterms,marray,ldc)
      else
         call rotviarecur(nd,theta,nterms,nterms,nterms,marray1,
     1           nterms,marray,ldc)
      endif

c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      do l=0,nterms2
         do m=-l,l
           do idim=1,nd
              mptemp(idim,l,m) = 0
              mp2(idim,l,m) = 0
           enddo
         enddo
      enddo



      rshift = d
      call h3dmpmpzshift
     $   (nd,wavek,sc1,marray,ldc,nterms,sc2,mptemp,
     1           nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2           ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
         call rotviaproj(nd,-theta,nterms2,nterms,nterms2,mptemp,
     1           nterms2,marray,ldc)
      else
         call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,mptemp,
     1           nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system
c         and add to existing expansion mpolen
c
      do l=0,nterms2
         do m=-l,l
            do idim=1,nd
              mpolen(idim,l,m)=mpolen(idim,l,m) + 
     1           ephi(-m)*marray(idim,l,m)
            enddo
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
c***********************************************************************
      subroutine h3dmpmpzshift
     $     (nd,zk,scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,radius,zshift,xnodes,wts,nquad,ynm,
     2      ynmd,mp2,phitemp,phitempn,fhs,fhder)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c INPUT:
c
c     nd       : number of multipole expansions
c     zk       : Helmholtz coefficient
c     scale    : scale parameter for mpole
c     mpole    : coefficients of original multipole exp.
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in the orig. expansion
c
c     scale2   : scale parameter for new expansion (mpolen)
c     lmpn     : leading dim of shifted (may be work array)
c     nterms2  : number of terms in output expansion
c     radius   : radius of sphere on which mpole is evaluated
c                           in projeciton step
c     zshift   : shifting distance along z-axis
c                              (always assumed positive)
C     xnodes   : Legendre nodes (precomputed)
C     wts      : Legendre weights (precomputed)
C     nquad    : number of quadrature nodes in theta direction
c
c OUTPUT:
c
c     mpolen  (complex *16)  : coefficients of shifted exp.
c
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer(8) nterms,nterms2,nquad,ier,lmp,lmpn,ldc,iynm,lynm,nd
      real *8   zshift,scale,scale2,radius
      real *8   xnodes(*),wts(*)
      real *8   ynm(0:nterms,0:nterms)
      real *8   ynmd(0:nterms,0:nterms)
      complex *16 phitemp(nd,nquad,-nterms:nterms)
      complex *16 phitempn(nd,nquad,-nterms:nterms)
      complex *16 fhs(0:nterms2)
      complex *16 fhder(0:nterms2)
      complex *16 mpole(nd,0:lmp,-lmp:lmp),zk
      complex *16 mpolen(nd,0:lmpn,-lmpn:lmpn)
      complex *16 mp2(nd,0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: rat1(:,:),rat2(:,:)
c
      integer(8) l,m,jnew,knew
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
C    OPTIMIZATION NOTES:
C
C    Suppose you are shifting from a very small sphere to the center
C    of a very large sphere (nterms2 >> nterms).
C    Then, ALONG THE Z-AXIS, the number of azimuthal modes that
C    need to be computed is only nterms (not nterms2). 
C    Subroutines h3dmpevalspherenm, h3dprojlocnmsep allow for this.
C    The final step of the point and shoot algorithm then distributes
C    these nterms (azimuthal) modes to nterms2 (azimuthal) in the
C    "laboratory frame".
C
C    cost is (nterms^2 x nterms2) rather than (nterms x nterms2^2)
C

        ldc = max(nterms,nterms2)

        allocate(rat1(0:ldc,0:ldc))
        allocate(rat2(0:ldc,0:ldc))


      call h3dmpevalsphere(nd,mpole,zk,scale,
     1     zshift,radius,nterms,lmp,ynm,ynmd,
     2     phitemp,phitempn,nquad,xnodes,fhs,fhder,rat1,rat2)
      call h3dprojloc
     $   (nd,nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,phitempn,mpolen,mp2,ynm,rat1,rat2)
      call h3drescalemp(nd,nterms2,lmpn,mpolen,radius,zk,
     1               scale2,fhs,fhder)
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h3dmploc(nd,wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad)
C***********************************************************************
C
C           Converts multipole expansion and adds to a local expansion
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C           NOTE: the output local expansion is not initialized
C              in this routine and the shifted 
C              multipole expansion is directly added to it
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nd     = number of expansions
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for local expansion
C           xnynzn = center of shifted local expansion
C           nterms2 = order of local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
C
C***********************************************************************
C
      implicit none

c
cc       calling sequence variables

      integer(8) nterms,nterms2,nquad,nd
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(nquad),wts(nquad),radius
      real *8 sc1,sc2
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      complex *16 local(nd,0:nterms2,-nterms2:nterms2)
      complex *16 wavek

c
cc      temporary variables
c
      complex *16, allocatable :: marray(:,:,:),marray1(:,:,:)
      complex *16, allocatable :: mptemp(:,:,:),mp2(:,:,:)
      complex *16, allocatable :: ephi(:),phitemp(:,:,:)
      complex *16, allocatable :: phitempn(:,:,:)
      complex *16, allocatable :: fhs(:),fhder(:)
      complex *16, allocatable :: fjs(:),fjder(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16 imag

      integer(8)  nq, ldc,idim
      real *8 rshift
      real *8 d,theta,ctheta,phi,rvec(3)
      integer(8)  l,m,jnew,knew
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)

      allocate(marray(nd,0:ldc,-ldc:ldc))
      allocate(mp2(nd,0:ldc,-ldc:ldc))
      allocate(marray1(nd,0:nterms,-nterms:nterms))
      allocate(mptemp(nd,0:nterms2,-nterms2:nterms2))
      allocate(ephi(-ldc-1:ldc+1))
      allocate(phitemp(nd,nq,-ldc:ldc))
      allocate(phitempn(nd,nq,-ldc:ldc))
      allocate(fhs(0:ldc),fhder(0:ldc))
      allocate(ynm(0:ldc,0:ldc),ynmd(0:ldc,0:ldc))
      allocate(fjs(0:nterms2),fjder(0:nterms2))

C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polar(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c     create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c     a rotation of THETA radians about the Yprime axis after PHI
c     radians about the z-axis.
      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               marray1(idim,l,m)  = mpole(idim,l,m)*ephi(m)
            enddo
         enddo
      enddo

      do l=0,nterms2
         do m=-l,l
            do idim=1,nd
               mptemp(idim,l,m)=0.0d0
            enddo
         enddo
      enddo

c
      if( nterms .ge. 30 ) then
         call rotviaproj(nd,theta,nterms,nterms,nterms,marray1,
     1       nterms,marray,ldc)
      else
         call rotviarecur(nd,theta,nterms,nterms,nterms,marray1,
     1          nterms,marray,ldc)
      endif

c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmploczshift(nd,wavek,marray,sc1,ldc,nterms,mptemp,
     1      sc2,nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if( nterms2 .ge. 30 ) then
         call rotviaproj(nd,-theta,nterms2,nterms2,nterms2,mptemp,
     1          nterms2,marray,ldc)
      else
         call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,mptemp,
     1          nterms2,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            do idim=1,nd
               local(idim,l,m)=local(idim,l,m) + 
     1             ephi(-m)*marray(idim,l,m)
            enddo
         enddo
      enddo


      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine h3dmploczshift
     $     (nd,zk,mpole,scale,lmp,nterms,local,
     1      scale2,lmpn,nterms2,radius,zshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c
c     nd       : number of expansions
c     zk       : Helmholtz parameter
c     mpole    : coefficients of original multipole exp.
c     scale    : scale parameter for mpole
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in original expansion
c
c     scale2   : scale parameter for local
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     radius   : radius of sphere about new center on which field
c                is evaluated
c     zshift   : shifting distance along z-axis
c                             (always assumed positive)
C     xnodes  = Legendre nodes (precomputed)
C     wts     = Legendre weights (precomputed)
C     nquad   = number of quadrature nodes in theta direction.
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer(8) nterms,nterms2,nquad,nd,lmp,lmpn
      integer(8) l,lw,m,jnew,knew
      real *8 zshift
      real *8 xnodes(*),wts(*)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 phitemp(nd,nquad,-nterms:nterms)
      complex *16 phitempn(nd,nquad,-nterms:nterms)
      complex *16 mp2(nd,0:lmpn,-lmpn:lmpn)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:nterms2)
      complex *16 fjder(0:nterms2)
      complex *16 mpole(nd,0:lmp,-lmp:lmp),zk
      complex *16 local(nd,0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace array
c
      real *8, allocatable :: rat1(:,:),rat2(:,:)
c

        ldc = max(nterms,nterms2)
        allocate(rat1(0:ldc,0:ldc))
        allocate(rat2(0:ldc,0:ldc))
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      call h3dmpevalsphere(nd,mpole,zk,scale,zshift,radius,
     2     nterms,lmp,ynm,ynmd,phitemp,phitempn,nquad,xnodes,
     3     fhs,fhder,rat1,rat2)
      call h3dprojloc
     $   (nd,nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,phitempn,local,mp2,ynm,rat1,rat2)
      call h3drescaleloc(nd,nterms2,lmpn,local,mp2,radius,zk,scale2,
     2     fjs,fjder)

      return
      end
c      
c
c
c
c
c
c
C***********************************************************************
      subroutine h3dlocloc(nd,wavek,sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad)
C***********************************************************************
C
C
C     Usage:
C
C           Shift center of a local expansion and add to output local
C           expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C           
C           NOTE: the output local expansion is not initialized
C              in this routine and the shifted 
C              local expansion is directly added to it
C
C     Input:
C
C           nd     = number of local expansions
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for locold expansion
C           x0y0z0 = center of original expansion
C           locold  = coefficients of original expansion
C           nterms = order of original expansion
C           sc2     = scaling parameter for local expansion
C           xnynzn = center of shifted expansion
C
C           nterms2 = order of shifted expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes used (really nquad**2)
C                     should be about 2*nterms2
C
C     Output:
C
C           local = coefficients of shifted expansion
C
C***********************************************************************
      implicit none
     
c
cc      calling sequence variables
c
      integer(8) nterms,nterms2,nquad,nd
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(nquad),wts(nquad),radius
      real *8 sc1,sc2
      complex *16 locold(nd,0:nterms,-nterms:nterms)
      complex *16 local(nd,0:nterms2,-nterms2:nterms2)
      complex *16 wavek
c
cc      temporary variables
c
      complex *16, allocatable :: marray(:,:,:),marray1(:,:,:)
      complex *16, allocatable :: mptemp(:,:,:),mp2(:,:,:)
      complex *16, allocatable :: ephi(:),phitemp(:,:,:)
      complex *16, allocatable :: phitempn(:,:,:)
      complex *16, allocatable :: fjs(:),fjder(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16 imag, ephi1

      integer(8)  nq, ldc,idim
      real *8 rshift
      real *8 d,theta,ctheta,phi,rvec(3)
      integer(8)  l,m,jnew,knew
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)

      allocate(marray(nd,0:ldc,-ldc:ldc))
      allocate(mp2(nd,0:ldc,-ldc:ldc))
      allocate(marray1(nd,0:nterms,-nterms:nterms))
      allocate(mptemp(nd,0:nterms2,-nterms2:nterms2))
      allocate(ephi(-ldc-1:ldc+1))
      allocate(phitemp(nd,nq,-ldc:ldc))
      allocate(phitempn(nd,nq,-ldc:ldc))
      allocate(ynm(0:ldc,0:ldc),ynmd(0:ldc,0:ldc))

      allocate(fjs(0:ldc),fjder(0:ldc))
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polar(rvec,d,theta,phi)
c
      ephi1 = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c
c      a rotation of THETA radians about the Yprime-axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      locold and ephi inside the following loop. 
c

      
      do l=0,nterms
         do m=-l,l
            do idim=1,nd
               marray1(idim,l,m) = locold(idim,l,m)*ephi(m)
            enddo
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            do idim=1,nd
               mptemp(idim,l,m)=0.0d0
            enddo
         enddo
      enddo

      if( nterms .ge. 30 ) then
         call rotviaproj(nd,theta,nterms,nterms,nterms2,marray1,
     1      nterms,marray,ldc)
      else
         call rotviarecur(nd,theta,nterms,nterms,nterms2,marray1,
     1         nterms,marray,ldc)
      endif
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
       call h3dlocloczshift(nd,wavek,sc1,marray,ldc,nterms,sc2,mptemp,
     1       nterms2,nterms2,radius,rshift,xnodes,wts,nquad,ynm,ynmd,
     2       mp2,phitemp,phitempn,fjs,fjder)


c
c      reverse THETA rotation.
c      I.e. rotation of -THETA about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
         call rotviaproj(nd,-theta,nterms2,nterms2,nterms2,mptemp,
     1         nterms2,marray,ldc)
      else
         call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,mptemp,
     1         nterms2,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            do idim=1,nd
               local(idim,l,m)=local(idim,l,m) + 
     1            ephi(-m)*marray(idim,l,m)
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
c***********************************************************************
      subroutine h3dlocloczshift
     $     (nd,zk,scale,locold,lmp,nterms,scale2,
     1      local,lmpn,nterms2,radius,zshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fjs,fjder) 
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
c     INPUT:
c
c     nd       : number of local expansions
c     zk       : Helmholtz parameter
c     scale    : scaling parameter for locold
c     locold   : coefficients of original multipole exp.
c     lmp      : leading dim of locold (may be a work array)
c     nterms   : number of terms in the orig. expansion
c     scale2   : scaling parameter for output expansion (local)
c
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
C     radius   : radius of sphere on which local expansion is
C                computed
c     zshift   : shifting distance along z-axis (assumed positive)
C     xnodes   : Legendre nodes (precomputed)
C     wts      : Legendre weights (precomputed)
C     nquad    : number of quadrature nodes used (really nquad**2)
C                     should be about 2*nterms2
c
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c     ier      : error return code
c                 CURRENTLY NOT USED
c
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer(8) nterms,nterms2,nquad,ier,l,lw,m,jnew,knew,nd
      integer(8) lmp,lmpn
      real *8 zshift
      real *8 xnodes(nquad),wts(nquad)
      real *8 ynm(0:lmp,0:lmp)
      real *8 ynmd(0:lmp,0:lmp)
      complex *16 phitemp(nd,nquad,-nterms2:nterms2)
      complex *16 phitempn(nd,nquad,-nterms2:nterms2)
      complex *16 mp2(nd,0:lmp,-lmp:lmp)
      complex *16 fjs(0:lmp)
      complex *16 fjder(0:lmp)
      complex *16 locold(nd,0:lmp,-lmp:lmp),zk
      complex *16 local(nd,0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
        
        ldc = max(nterms,nterms2)
        irat1=1
        lrat1=(ldc+1)**2
        irat2=irat1+lrat1
        lrat2=(ldc+1)**2
        lused=irat2+lrat2
        allocate (w(lused))
C
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
C    OPTIMIZATION NOTES:
C
C    Suppose you are shifting from a very large sphere to a very
C    small sphere (nterms >> nterms2).
C    Then, ALONG THE Z-AXIS, the number of azimuthal modes that
C    need to be computed is only nterms2 (not nterms). 
C    Subroutines h3dlocevalspherestab, h3dprojlocsepstab allow for this.
C
C    cost is (nterms2^2 x nterms) rather than (nterms2 x nterms^2)
C
C
C

      call h3dlocevalsphere(nd,locold,zk,scale,
     1      zshift,radius,nterms,nterms2,
     1     lmp,ynm,ynmd,phitemp,phitempn,nquad,xnodes,fjs,
     1     fjder,w(irat1),w(irat2))

      call h3dprojloc
     $   (nd,nterms2,lmpn,nquad,nterms2,xnodes,wts,
     1     phitemp,phitempn,local,mp2,ynm,w(irat1),w(irat2))

      call h3drescaleloc(nd,nterms2,lmpn,local,mp2,
     1      radius,zk,scale2,fjs,fjder)
      


      return
      end
C
C
