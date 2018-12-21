cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c    Documentation in progress
c
c     multipole shift routines, f95 version, using allocate
c
C***********************************************************************
      subroutine h3dmpmpquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           ier   = error return flag
C
C                   CURRENTLY UNUSED.
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2  
      iynmd = iynm+lynm
      iphitemp = iynmd+lynm
      lphitemp = nq*(2*ldc+1)*2 + 7  
      iphitemp2 = iphitemp+lphitemp
      ifhs = iphitemp2+lphitemp
      ifhder = ifhs + 2*(ldc+1) + 3
      lused = ifhder + 2*(ldc+1) +3 + 100
      allocate (w(lused))
c
ccc      call prinf(' ier before mpmpquad0 is *',ier,1)
      call h3dmpmpquad0(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,w(imarray),w(imarray1),
     2           ldc,w(iephi),
     3           radius,xnodes,wts,nquad,nq,w(iynm),
     4           w(iynmd),w(iphitemp),w(iphitemp2),
     7           w(ifhs),w(ifhder),ier)
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dmpmpquadu_add(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           ier   = error return flag
C
C                   CURRENTLY UNUSED.
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:ldc,-ldc:ldc)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dmpmpquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            mpolen(l,m) = mpolen(l,m)+mptemp(l,m)
         enddo
      enddo
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquad0(wavek,sc1,x0y0z0,mpole,nterms,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           phitemp,phitemp2,fhs,fhder,ier)
C***********************************************************************
C
C     Usage:
C
C           Shift multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of shifted expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for shifted expansion
C           radius  = radius of sphere on which shifted expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta
C           nq      = max(nquad, 2*ldc+2) = used to allocate
C                     work arrays for both z-shift and rotations.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms, lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 radius, rshift
      real *8 xnodes(1),wts(1)
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitemp2(nq,-ldc:ldc)
      complex *16 fhs(0:ldc)
      complex *16 fhder(0:ldc)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
c
      complex *16 ephi(-ldc-1:ldc+1),imag
ccc      complex *16 ephi2(-ldc-1:ldc+1)
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
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
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms .ge. 30 ) then 
      call rotviaprojf90(theta,nterms,nterms,nterms,marray1,nterms,
     1        marray,ldc)
      else
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,nterms,
     1        marray,ldc)
      endif
ccc      call prinf(' after rotviaproj *',nterms,1)
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmpmpzshift_fast
     $   (wavek,sc1,marray,ldc,nterms,sc2,mpolen,
     1           nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2           ynm,phitemp,fhs,fhder,ier)
ccc      call prinf(' after mpmpzshift *',nterms,1)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
ccc      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,mpolen,
      call rotviaprojf90(-theta,nterms2,nterms,nterms2,mpolen,
     1        nterms2,marray,ldc)
      else
ccc      call rotviarecur3f90(-theta,nterms2,nterms,nterms2,mpolen,
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquadu_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,mpolen,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted mpole expansion
C
C           ier   = error return flag
C
C                   CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C     WORK ESTIMATE:
C
C     Let nmax = max(nterms,nterms2). Then
C     w should be at least 9*(nmax)^2+ 27*(nmax) 
C                          + 8*nquad*nmax + 2*nquad + 100
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2  
      iynmd = iynm+lynm
      iphitemp = iynmd+lynm
      lphitemp = nq*(2*ldc+1)*2 + 7  
      iphitemp2 = iphitemp+lphitemp
      ifhs = iphitemp2+lphitemp
      ifhder = ifhs + 2*(ldc+1) + 3
      lused = ifhder + 2*(ldc+1) +3
      allocate(w(lused))
c
ccc      call prinf(' ier before mpmpquad0 is *',ier,1)
      call h3dmpmpquad_trunc0(wavek,sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,mpolen,nterms2,w(imarray),w(imarray1),
     2           ldc,w(iephi),
     3           radius,xnodes,wts,nquad,nq,w(iynm),
     4           w(iynmd),w(iphitemp),w(iphitemp2),
     7           w(ifhs),w(ifhder),ier)
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquad_trunc0(wavek,sc1,x0y0z0,mpole,
     $           nterms,nterms1,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           phitemp,phitemp2,fhs,fhder,ier)
C***********************************************************************
C
C     Usage:
C
C           Shift multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of shifted expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for shifted expansion
C           radius  = radius of sphere on which shifted expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta
C           nq      = max(nquad, 2*ldc+2) = used to allocate
C                     work arrays for both z-shift and rotations.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms,nterms1
      integer lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 radius, rshift
      real *8 xnodes(1),wts(1)
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitemp2(nq,-ldc:ldc)
      complex *16 fhs(0:ldc)
      complex *16 fhder(0:ldc)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
c
      complex *16 ephi(-ldc-1:ldc+1),imag
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
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
      do l=0,nterms1
         do m=-l,l
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms .ge. 30 ) then 
      call rotviaprojf90(theta,nterms1,nterms1,nterms1,marray1,nterms1,
     1        marray,ldc)
      else
      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
     1        nterms1,marray,ldc)
      endif
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmpmpzshift_fast
     $   (wavek,sc1,marray,ldc,nterms1,sc2,mpolen,
     1           nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2           ynm,phitemp,fhs,fhder,ier)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      else
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine h3dmpmpzshift_fast
     $     (zk,scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,radius,zshift,xnodes,wts,nquad,ynm,
     2      phitemp,fhs,fhder,ier)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c INPUT:
c
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
      integer nterms,nterms2,nquad,ier,lmp,lmpn,ldc,iynm,lynm
      real *8   zshift,scale,scale2,radius
      real *8   xnodes(1),wts(1)
      real *8   ynm(0:nterms,0:nterms)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 fhs(0:nterms2)
      complex *16 fhder(0:nterms2)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 mpolen(0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
      integer l,m,jnew,knew
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
        irat1=1
        lrat1=(ldc+1)**2
        irat2=irat1+lrat1
        lrat2=(ldc+1)**2
        lused=irat2+lrat2
        allocate(w(lused))
c
ccc      call prinf(' allocated in shift fast *',lused,1)
      call h3dmpevalspherenm_fast(mpole,zk,scale,
     1     zshift,radius,nterms,lmp,ynm,
     2     phitemp,nquad,xnodes,fhs,fhder,w(irat1),w(irat2))
      call h3dprojlocnmsep_fast
     $   (nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,mpolen,ynm,w(irat1),w(irat2))
      call h3drescalemp(nterms2,lmpn,mpolen,radius,zk,
     1               scale2,fhs,fhder)
      return
      end
C
cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c     multipole to local translation, f95 version using allocate
c
C***********************************************************************
      subroutine h3dmplocquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                   CURRENTLY UNUSED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag,wavek
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2
      iynmd = iynm+lynm
      imp2 = iynmd+lynm
      iphitemp = imp2+(ldc+1)*(2*ldc+1)*2
      lphitemp = nq*(2*ldc+1)*2
      iphitempn = iphitemp+lphitemp
      ifhs = iphitempn+lphitemp
      ifhder = ifhs+ 2*(nterms+1) + 3
      ifjs = ifhder+ 2*(nterms+1) + 3
      lwfjs = nterms2+1000
      lfjs = 2*(lwfjs+1) + 3
      ifjder = ifjs+lfjs
      lfjder = 2*(nterms2+1)+3
      iiscale = ifjder+lfjder
      liscale = (lwfjs+1)+3
      lused = iiscale+ liscale
      allocate(w(lused))
c
      call h3dmplocquad0(wavek,sc1,x0y0z0,mpole,nterms,sc2,xnynzn,
     1         local,nterms2,w(imarray),w(imarray1),ldc,
     2         w(iephi),radius,xnodes,wts,nquad,nq,
     3         w(iynm),w(iynmd),w(imp2),
     4         w(iphitemp),w(iphitempn),w(ifhs),w(ifhder),
     5         w(ifjs),w(ifjder),w(iiscale),lwfjs,ier)
      return
      end
c
c
C***********************************************************************
      subroutine h3dmplocquadu_add(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                   CURRENTLY UNUSED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag,wavek
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dmplocquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmplocquad0(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,mp2,
     3           phitemp,phitempn,fhs,fhder,fjs,fjder,iscale,lwfjs,ier)

C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted local expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of local expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes used (really nquad**2)
C
C           lots of work arrays..... DOCUMENT!
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
c           ier      : error return code
c                CURRENTLY UNUSED
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      integer  iscale(0:lwfjs)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1),rvec(3)
      real *8 zshift
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitempn(nq,-ldc:ldc)
      complex *16 mp2(0:ldc,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
      complex *16 ephi(-ldc-1:ldc+1),imag
      data imag/(0.0d0,1.0d0)/
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
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            local(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms .ge. 30 ) then
      call rotviaprojf90(theta,nterms,nterms,nterms,marray1,nterms,
     1     marray,ldc)
      else
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,nterms,
     1     marray,ldc)
      endif
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmploczshiftstab_fast(wavek,marray,sc1,ldc,nterms,local,
     1      sc2,nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,local,nterms2,
     1     marray,ldc)
      else
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dmplocquadu_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                       CURRENTLY UNUSED
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag,wavek
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
c
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2
      iynmd = iynm+lynm
      imp2 = iynmd+lynm
      iphitemp = imp2+(ldc+1)*(2*ldc+1)*2
      lphitemp = nq*(2*ldc+1)*2
      iphitempn = iphitemp+lphitemp
      ifhs = iphitempn+lphitemp
      ifhder = ifhs+ 2*(nterms+1) + 3
      ifjs = ifhder+ 2*(nterms+1) + 3
      lwfjs = nterms2+1000
      lfjs = 2*(lwfjs+1) + 3
      ifjder = ifjs+lfjs
      lfjder = 2*(nterms2+1)+3
      iiscale = ifjder+lfjder
      liscale = (lwfjs+1)+3        
      lused = iiscale+ liscale
      allocate(w(lused))
c
      call h3dmplocquad_trunc0(wavek,sc1,x0y0z0,mpole,nterms,nterms1,
     1         sc2,xnynzn,local,nterms2,w(imarray),w(imarray1),ldc,
     2         w(iephi),radius,xnodes,wts,nquad,nq,
     3         w(iynm),w(iynmd),w(imp2),
     4         w(iphitemp),w(iphitempn),w(ifhs),w(ifhder),
     5         w(ifjs),w(ifjder),w(iiscale),lwfjs,ier)
      return
      end
c
c
C***********************************************************************
      subroutine h3dmplocquadu_add_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                       CURRENTLY UNUSED
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag,wavek
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dmplocquadu_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmplocquad_trunc0(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,nterms2,marray,marray1,ldc,
     2           ephi,radius,xnodes,wts,nquad,nq,ynm,ynmd,mp2,
     4           phitemp,phitempn,fhs,fhder,fjs,fjder,iscale,lwfjs,ier)

C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted local expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of local expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes used (really nquad**2)
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
c           ier      : error return code
c              8      lwfjs insufficient for jfuns3d in h3drescale.
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           ldc      must exceed max(nterms,nterms2).
C           rd1,rd2  work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C           w       = work array 
C
C           LOTS MORE
C
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      integer  iscale(0:lwfjs)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1),rvec(3)
      real *8 zshift
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitempn(nq,-ldc:ldc)
      complex *16 mp2(0:ldc,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
      complex *16 ephi(-ldc-1:ldc+1),imag
      data imag/(0.0d0,1.0d0)/
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
      do l=0,nterms1
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            local(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms1 .ge. 30 ) then
      call rotviaprojf90(theta,nterms1,nterms1,nterms1,marray1,
     1     nterms1,marray,ldc)
      else
      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
     1     nterms1,marray,ldc)
      endif
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmploczshiftstab_fast(wavek,marray,sc1,ldc,nterms1,local,
     1      sc2,nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
      else
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dmplocquadu2_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad,ier,rotmatf,rotmatb,ldm)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                       CURRENTLY UNUSED
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag,wavek
        real *8 rotmatf(1)
        real *8 rotmatb(1)
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
c
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2
      iynmd = iynm+lynm
      imp2 = iynmd+lynm
      iphitemp = imp2+(ldc+1)*(2*ldc+1)*2
      lphitemp = nq*(2*ldc+1)*2
      iphitempn = iphitemp+lphitemp
      ifhs = iphitempn+lphitemp
      ifhder = ifhs+ 2*(nterms+1) + 3
      ifjs = ifhder+ 2*(nterms+1) + 3
      lwfjs = nterms2+1000
      lfjs = 2*(lwfjs+1) + 3
      ifjder = ifjs+lfjs
      lfjder = 2*(nterms2+1)+3
      iiscale = ifjder+lfjder
      liscale = (lwfjs+1)+3        
      lused = iiscale+ liscale
      allocate(w(lused))
c
      call h3dmplocquad2_trunc0(wavek,sc1,x0y0z0,mpole,nterms,nterms1,
     1         sc2,xnynzn,local,nterms2,w(imarray),w(imarray1),ldc,
     2         w(iephi),radius,xnodes,wts,nquad,nq,
     3         w(iynm),w(iynmd),w(imp2),
     4         w(iphitemp),w(iphitempn),w(ifhs),w(ifhder),
     5         w(ifjs),w(ifjder),w(iiscale),lwfjs,ier,
     $         rotmatf,rotmatb,ldm)
      return
      end
c
c
C***********************************************************************
      subroutine h3dmplocquadu2_add_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier,rotmatf,rotmatb,ldm)
C***********************************************************************
C
C     Memory management wrapper for subroutine h3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C           lused = amount of workspace w actually used.
C           ier   = error return flag
C
C                       CURRENTLY UNUSED
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag,wavek
        real *8 rotmatf(1)
        real *8 rotmatb(1)
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dmplocquadu2_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier,rotmatf,rotmatb,ldm)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmplocquad2_trunc0(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,local,nterms2,marray,marray1,ldc,
     2           ephi,radius,xnodes,wts,nquad,nq,ynm,ynmd,mp2,
     4           phitemp,phitempn,fhs,fhder,fjs,fjder,iscale,lwfjs,ier,
     $     rotmatf,rotmatb,ldm)

C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted local expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of local expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for local expansion
C           radius  = radius of sphere on which local expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes used (really nquad**2)
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of shifted local expansion
c           ier      : error return code
c              8      lwfjs insufficient for jfuns3d in h3drescale.
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           ldc      must exceed max(nterms,nterms2).
C           rd1,rd2  work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C           w       = work array 
C
C           LOTS MORE
C
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      integer  iscale(0:lwfjs)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1),rvec(3)
      real *8 zshift
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitempn(nq,-ldc:ldc)
      complex *16 mp2(0:ldc,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
      complex *16 ephi(-ldc-1:ldc+1),imag
        real *8 rotmatf(1)
        real *8 rotmatb(1)
      data imag/(0.0d0,1.0d0)/
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
      do l=0,nterms1
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            local(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms1 .ge. 30 ) then
      call rotviaprojf90(theta,nterms1,nterms1,nterms1,marray1,
     1     nterms1,marray,ldc)
      else
ccc      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
ccc     1     nterms1,marray,ldc)
        call rotviarecur3p_apply(theta,nterms1,nterms1,nterms1,marray1,
     1     nterms1,marray,ldc,rotmatf,ldm)
      endif
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmploczshiftstab_fast(wavek,marray,sc1,ldc,nterms1,local,
     1      sc2,nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
      else
ccc      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
ccc     1     nterms2,marray,ldc)
      call rotviarecur3p_apply(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc,rotmatb,ldm)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c***********************************************************************
      subroutine h3dmploczshiftstab_fast
     $     (zk,mpole,scale,lmp,nterms,local,
     1      scale2,lmpn,nterms2,radius,zshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,phitemp,phitempn,fhs,fhder,fjs,fjder,
     3      iscale,lwfjs,ier)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c
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
c     ier      : error return code
c                  CURRENTLY UNUSED
c
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nterms2,nquad,ier
      integer l,lw,m,jnew,knew
      integer iscale(0:lwfjs)
      real *8 zshift
      real *8 xnodes(1),wts(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 mp2(0:lmpn,-lmpn:lmpn)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 local(0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c

        ldc = max(nterms,nterms2)
        irat1=1
        lrat1=(ldc+1)**2
        irat2=irat1+lrat1
        lrat2=(ldc+1)**2
        lused=irat2+lrat2
        allocate(w(lused))
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      call h3dmpevalspherenmstab_fast(mpole,zk,scale,zshift,radius,
     2     nterms,lmp,ynm,ynmd,phitemp,phitempn,nquad,xnodes,
     3     fhs,fhder,w(irat1),w(irat2))
      call h3dprojlocsepstab_fast
     $   (nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,phitempn,local,mp2,ynm,w(irat1),w(irat2))
      call h3drescalestab(nterms2,lmpn,local,mp2,radius,zk,scale2,
     2     fjs,fjder,iscale,lwfjs,ier)

      return
      end
C
C
cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c     Local to local shift routines, f95 version using allocate
c
C***********************************************************************
      subroutine h3dloclocquadu(wavek,sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dloclocquad0 (below).
C
C     Usage:
C
C           Shift center of a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C     Input:
C
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
C           ier   = error return flag
C
C                       CURRENTLY NOT USED
C
C     Work arrays carved out of w:
C
C           marray = work array used to hold various intermediate 
c                    values.
C           dc     = work array contain the square roots of 
C                    som binomial coefficients.
C           rd1,rd2 = work arrays used to compute rotation matrices
C                     about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2
      iynmd = iynm+lynm
      imp2 = iynmd+lynm
      iphitemp = imp2+(ldc+1)*(2*ldc+1)*2
      lphitemp = nq*(2*ldc+1)*2
      iphitempn = iphitemp+lphitemp
      iiscale = iphitempn+lphitemp
      lwfjs = ldc+1000
      ifjs = iiscale+lwfjs+1
      ifjder = ifjs+2*(lwfjs+1)
      lused = ifjder+2*(lwfjs+1)
      allocate (w(lused))
c
      call h3dloclocquad0(wavek,sc1,x0y0z0,locold,nterms,sc2,xnynzn,
     1           local,nterms2,w(imarray),w(imarray1),ldc,
     2           w(iephi),radius,xnodes,wts,nquad,nq,
     3           w(iynm),w(iynmd),
     4           w(imp2),w(iphitemp),w(iphitempn),w(ifjs),
     5           w(ifjder),w(iiscale),lwfjs,ier)
      return
      end
c
c
C***********************************************************************
      subroutine h3dloclocquadu_add(wavek,sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dloclocquad0 (below).
C
C     Usage:
C
C           Shift center of a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C     Input:
C
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
C           ier   = error return flag
C
C                       CURRENTLY NOT USED
C
C     Work arrays carved out of w:
C
C           marray = work array used to hold various intermediate 
c                    values.
C           dc     = work array contain the square roots of 
C                    som binomial coefficients.
C           rd1,rd2 = work arrays used to compute rotation matrices
C                     about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 wavek,imag
c
c     local allocated workspace arrays - no more passed workspace
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dloclocquadu(wavek,sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dloclocquad0(wavek,sc1,x0y0z0,locold,nterms,sc2,
     1           xnynzn,local,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           mp2,phitemp,phitempn,fjs,fjder,iscale,lwfjs,ier) 
C***********************************************************************
C
C     Usage:
C
C           Shifts center of a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     wavek   : Helmholtz parameter
C     sc1     : scaling parameter for locold expansion
C     x0y0z0  : center of original multiple expansion
C     locold  : coefficients of original multiple expansion
C     nterms  : order of original local expansion
C     sc2     : scaling parameter for local expansion
C     xnynzn  : center of shifted local expansion
c
C     nterms2 : order of new local expansion
c     marray  : work array
c     dc      : another work array
c     ldc     : dimension parameter for marray and ldc
c               must exceed max(nterms,nterms2).
c     rd1     : work array for rotation operators.
c     rd2     : work array for rotation operators.
c     ephi    : work array for rotation operators.
C     radius  : radius of sphere on which local expansion is
C               computed
C     xnodes  : Legendre nodes (precomputed)
C     wts     : Legendre weights (precomputed)
C     nquad   : number of quadrature nodes in theta direction.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local   : coefficients of shifted local expansion
c     ier     : error return code 
c
c               CURRENTLY NOT USED
C
C***********************************************************************
C
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jnew,knew
      integer iscale(0:lwfjs)
      real *8 x0y0z0(3),xnynzn(3),rvec(3)
      real *8 xnodes(1),wts(1)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitempn(nq,-ldc:ldc)
      complex *16 mp2(0:ldc,-ldc:ldc)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 marray1(0:nterms,-nterms:nterms)
      complex *16 wavek,imag,ephi1
      complex *16 ephi(-ldc-1:ldc+1)
      data imag/(0.0d0,1.0d0)/
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
         do mp=-l,l
            marray1(l,mp) = locold(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            local(l,m)=0.0d0
         enddo
      enddo
ccc      t1 = second()
      if( nterms .ge. 30 ) then
      call rotviaprojf90(theta,nterms,nterms,nterms2,marray1,nterms,
     1      marray,ldc)
      else
      call rotviarecur3f90(theta,nterms,nterms,nterms2,marray1,
     1      nterms,marray,ldc)
      endif
ccc      t2 = second()
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
ccc      t1 = second()
       call h3dlocloczshiftstab_fast
     $   (wavek,sc1,marray,ldc,nterms,sc2,local,
     1           nterms2,nterms2,radius,rshift,
     2           xnodes,wts,nquad,ynm,ynmd,mp2,
     3      phitemp,phitempn,fjs,fjder,iscale,lwfjs,ier) 
ccc      t2 = second()
c
c      reverse THETA rotation.
c      I.e. rotation of -THETA about Yprime axis.
c
ccc      t1 = second()
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,local,
     1      nterms2,marray,ldc)
      else
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1      nterms2,marray,ldc)
      endif
ccc      t2 = second()
ccc      call prin2(' time for second rot is *',t2-t1,1)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine h3dlocloczshiftstab_fast
     $     (zk,scale,locold,lmp,nterms,scale2,
     1      local,lmpn,nterms2,radius,zshift,xnodes,wts,nquad,
     2      ynm,ynmd,mp2,
     3      phitemp,phitempn,fjs,fjder,iscale,lwfjs,ier) 
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
c     INPUT:
c
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
      integer nterms,nterms2,nquad,ier,l,lw,m,jnew,knew
      integer iscale(0:lwfjs)
      real *8   zshift
      real *8   xnodes(1),wts(1)
      real *8 ynm(0:lmp,0:lmp)
      real *8 ynmd(0:lmp,0:lmp)
ccc      complex *16 phitemp(nquad,-lmp:lmp)
ccc      complex *16 phitempn(nquad,-lmp:lmp)
      complex *16 phitemp(nquad,-nterms2:nterms2)
      complex *16 phitempn(nquad,-nterms2:nterms2)
      complex *16 mp2(0:lmp,-lmp:lmp)
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:lwfjs)
      complex *16 locold(0:lmp,-lmp:lmp),zk
      complex *16 local(0:lmpn,-lmpn:lmpn)
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
      call h3dlocevalspherestab_fast(locold,zk,scale,
     1      zshift,radius,nterms,nterms2,
     1     lmp,ynm,ynmd,phitemp,phitempn,nquad,xnodes,
     1     iscale,fjs,fjder,w(irat1),w(irat2),lwfjs,ier)
      call h3dprojlocsepstab_fast
     $   (nterms2,lmpn,nquad,nterms2,xnodes,wts,
     1     phitemp,phitempn,local,mp2,ynm,w(irat1),w(irat2))
      call h3drescalestab(nterms2,lmpn,local,mp2,
     1      radius,zk,scale2,fjs,fjder,iscale,lwfjs,ier)
      return
      end
C
C
