cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc         and Manas Rachh
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
c    This file contains the translation operator for the Laplacd
c      3D FMM
c
c     l3mpmp - shift center of multipole expansion
c     l3dmploc - convert multipole expansion to a local expansion
c     l3locloc - shift center of local expansion
c
c
C
c
C***********************************************************************
      subroutine l3dmpmp(nd,sc1,x0y0z0,mpole,nterms,sc2,
     1           xnynzn,mpolen,nterms2,dc,lca)
C***********************************************************************
C
C     Usage:
C
C     Shift center of multipole expansion mpole and add to output
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then doing the shifting
C     along the Z-axis, and then rotating back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd      : number of mulitpole expansions 
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for shifted expansion
C     xnynzn  : center of shifted expansion
C     nterms2 : order of shifted expansion
C     dc      : square root of binomial coefficients
C     lca     : length of binomial coefficients
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen  : coefficients of shifted expansion
C
C***********************************************************************
C
      implicit none
      integer  nterms, lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      integer lca
      integer idim,nd
      double precision x0y0z0(3),xnynzn(3)
      double precision rshift
      double precision d,theta,ctheta,phi,sc1,sc2,rvec(3)
      double precision dc(0:lca,0:lca)
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex mpolen(nd,0:nterms2,-nterms2:nterms2)
      double complex, allocatable :: marray1(:,:,:)
      double complex, allocatable :: marray(:,:,:)
      double complex, allocatable :: ephi(:)
c
      double complex imag
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)

      ldc = max(nterms,nterms2)
      allocate(marray1(0:ldc,-ldc:ldc))
      allocate(marray(0:ldc,-ldc:ldc))
      allocate(ephi(-ldc-1:ldc+1))

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
          do idim=1,nd
            marray1(idim,l,m)=mpole(idim,l,m)*ephi(m)
          enddo
        enddo
      enddo


      if(nterms.ge.30) then
        call rotviaproj(nd,theta,nterms,nterms,nterms,marray1,
     1         ldc,marray,ldc)
      else
        call rotviarecur(nd,theta,nterms,nterms,nterms,marray1,
     1        ldc,marray,ldc)
      endif

      
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmpmpzshift(nd,sc1,marray,ldc,nterms,sc2,marray1,
     1           ldc,nterms2,rshift,dc,lca)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c

      if(nterms2.ge.30) then
        call rotviaproj(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1        nterms2,marray,ldc)
      else
        call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1        nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
        do m=-l,l
          do idim=1,nd
            mpolen(idim,l,m)=mpolen(idim,l,n)+ephi(-m)*marray(idim,l,m)
          enddo
        enddo
      enddo

      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dmpmpzshift(nd,scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,zshift,dc,lca)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c     INPUT:
c
c     nd       : number of multipole expansions
c     scale    : scale parameter for mpole
c     mpole    : coefficients of original multipole exp.
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in the orig. expansion
c     scale2   : scale parameter for new expansion (mpolen)
c     lmpn     : leading dim of shifted (may be work array)
c     nterms2  : number of terms in output expansion
c     zshift   : shifting distance along z-axis
c                              (always assumed positive)
c     dc       : square root of binomial coefficients
c     lca      : length of binomial coefficient array
c
c     OUTPUT:
c
c     mpolen  (double complex)  : coefficients of shifted exp.
c
c-----------------------------------------------------------------------
      implicit none
      integer l0,nmax,nterms,nterms2,nquad,ier,lmp,lmpn,ldc
      integer lca
      integer l,m,jnew,knew
      double precision d,zshift,scale,scale2,rat
      double precision, allocatable :: rscpow(:)
      double complex mpole(nd,0:lmp,-lmp:lmp)
      double complex mpolen(nd,0:lmpn,-lmpn:lmpn)
      double precision dc(0:lca,0:lca)
      double precision, allocatable :: fr(:)
C
C----- shift along z-axis
C
      nmax = max(nterms,nterms2)
      allocate(fr(0:nmax))
c
      fr(0) = 1.0d0
      d = -zshift*scale
      fr(1) = d
      do l=2,nmax
         fr(l) = fr(l-1)*d
      enddo

      d = scale2/scale
      allocate(rscpow(nterms2))
      rscpow(1) = d
      do i=2,nterms2
        rscpow(i) = rscpow(i-1)*d
      enddo
c
      do jnew = 0,nterms2
        l0 = max(0,jnew-nterms)
        do knew = -jnew,jnew
          do idim=1,nd
	        mpolen(idim,jnew,knew) = 0.0d0
          enddo
	      do l=l0,jnew-iabs(knew)
            do idim=1,nd
              mpolen(idim,jnew,knew)=mpolen(idim,jnew,knew)+
     1              mpole(idim,jnew-l,knew)*fr(l)*dc(jnew-knew,l)*
     2              dc(jnew+knew,l)
            enddo
          enddo
        enddo
      enddo

      do jnew = 1,nterms2
        do knew = -jnew,jnew
          do idim=1,nd
            mpolen(idim,jnew,knew)=mpolen(idim,jnew,knew)*rscpow(jnew)
          enddo
        enddo
      enddo

      return
      end
C
c
c
c
c
C***********************************************************************
      subroutine l3dmploc(nd,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,dc,lca)
C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion and add
C           to existing expansion
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C     nd         number of multipole expansions
C     sc1        scaling parameter for mpole expansion
C     x0y0z0     center of original multiple expansion
C     mpole      coefficients of original multiple expansion
C     nterms     order of multipole expansion
C     sc2        scaling parameter for local expansion
C     xnynzn     center of shifted local expansion
C     nterms2    order of local expansion
C     dc         precomputed square roots of binomial coefficients
C     lca        length of precomputed binomial coefficients
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local      coefficients of shifted local expansion
C---------------------------------------------------------------------
      implicit none
      integer idim,nd,lca
      integer  nterms,ier,l,m,jnew,knew,nterms2,ldc,mp
      double precision d,theta,ctheta,phi,sc1,sc2
      double precision x0y0z0(3),xnynzn(3)
      double precision rvec(3)
      double precision rshift
      double precision dc(0:lca,0:lca)
      double complex, allocatable :: marray1(:,:,:),marray(:,:,:),
     1    ephi(:)
      double complex mpole(nd,0:nterms,-nterms:nterms)
      double complex local(nd,0:nterms2,-nterms2:nterms2)
      double complex imag
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c

      ldc = max(nterms,nterms2)
      allocate(ephi(-ldc-1:ldc+1))
      allocate(marray(nd,0:ldc,-ldc:ldc),marray1(nd,0:ldc,-ldc:ldc))

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
c     The PHI rotation is carried out on the fly by multiplying 
c     mpole and ephi inside the following loop. 
c
      do l=0,nterms
        do mp=-l,l
          do idim=1,nd  
            marray1(idim,l,mp) = mpole(idim,l,mp)*ephi(mp)
          enddo
        enddo
      enddo

      if(nterms.ge.30) then
        call rotviaproj(nd,theta,nterms,nterms,nterms,marray1,
     1     ldc,marray,ldc)
      else
        call rotviarecur(nd,theta,nterms,nterms,nterms,marray1,
     1     ldc,marray,ldc)
      endif
c
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmploczshift(nd,marray,sc1,ldc,nterms,marray1,
     1      sc2,nterms2,nterms2,rshift,dc,lca)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if(nterms.ge.30) then
        call rotviaproj(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1     ldc,marray,ldc)
      else
        call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1     ldc,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
        do m=-l,l
          do idim=1,nd
            local(idim,l,m)=local(idim,l,m)+ephi(-m)*marray(idim,l,m)
          enddo
        enddo
      enddo

      return
      end
c
c
c***********************************************************************
      subroutine l3dmploczshift(nd,mpole,scale,lmp,nterms,local,
     2      scale2,lmpn,nterms2,zshift,dc,lca)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c     nd       : number of multipole expansions
c     mpole    : coefficients of original multipole exp.
c     scale    : scale parameter for mpole
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in original expansion
c
c     scale2   : scale parameter for local
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis
c                             (always assumed positive)
c     dc       : square root of binomial coefficients
c     lca      : length of dc
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
C---------------------------------------------------------------------
      implicit none
      integer idim,nd,lca
      integer nterms,nterms2,nquad,ier
      integer l,lw,m,jnew,knew,kk,ll,msign,nmax,lmpn,lmp
      double precision scale,scale2
      double precision zshift,d
      double precision dc(0:lca,0:lca)
      double precision, allocatable :: fr(:)
      double precision, allocatable :: rscpow(:)
      double complex mpole(nd,0:lmp,-lmp:lmp)
      double complex local(nd,0:lmpn,-lmpn:lmpn)
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      nmax = max(nterms,nterms2)
      allocate(rscpow(nmax),fr(0:2*nmax))
c
      d = 1.0d0/zshift
      fr(0) = d
      d = d/scale
      fr(1) = fr(0)*d
      do l=2,2*nmax
        fr(l) = fr(l-1)*d
      enddo

      d = scale/scale2
      rscpow(1) = d
      do l=2,nmax
        rscpow(l) = rscpow(l-1)*d
      enddo
c
      do jnew = 0,nterms2
        do knew = -jnew,jnew
          do idim=1,nd
   	        local(idim,jnew,knew) = 0.0d0
          enddo
          kk  =iabs(knew)
          msign = (-1)**(jnew+knew)
          do l=kk,nterms
            ll=l+jnew
            do idim=1,nd
              local(idim,jnew,knew)=local(idim,jnew,knew)+
     1          mpole(idim,l,knew)*fr(ll)*dc(ll,l-knew)*
     2          dc(ll,l+knew)*msign
            enddo 
          enddo
        enddo
      enddo

      do jnew = 1,nterms2
        do knew = -jnew,jnew
          do idim=1,nd
            local(idim,jnew,knew)=local(idim,jnew,knew)*rscpow(jnew)
          enddo
        enddo
      enddo

      return
      end
C
c
c
c
c
c
C***********************************************************************
C     Local -> Local routines
C***********************************************************************
c
c
C***********************************************************************
      subroutine l3dlocloc(nd,sc1,x0y0z0,locold,nterms,sc2,
     1           xnynzn,local,nterms2,dc,lda)
C***********************************************************************
C
C     Usage:
C
C           Shifts center of a local expansion and add to existing
C           expansion
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd      : number of local expansions
C     sc1     : scaling parameter for locold expansion
C     x0y0z0  : center of original multiple expansion
C     locold  : coefficients of original multiple expansion
C     nterms  : order of original local expansion
C     sc2     : scaling parameter for local expansion
C     xnynzn  : center of shifted local expansion
C     nterms2 : order of new local expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local   : coefficients of shifted local expansion
C
C***********************************************************************
C
      implicit none
      integer idim,lda
      integer nterms,ier,l,m,jnew,knew,nterms2,ldc,mp
      double precision x0y0z0(3),xnynzn(3),rvec(3)
      double precision d,theta,ctheta,phi,sc1,sc2
      double precision dc(0:lda,0:lda)
      double precision rshift
      double complex locold(nd,0:nterms,-nterms:nterms)
      double complex local(nd,0:nterms2,-nterms2:nterms2)
      double complex, allocatable :: marray(:,:,:),marray1(:,:,:),
     1                  ephi(:)
      double complex imag,ephi1
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)

      ldc = max(nterms,nterms2)
      allocate(ephi(-ldc-1:ldc+1))

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

      allocate(marray(nd,0:ldc,-ldc:ldc),marray1(nd,0:ldc,-ldc:ldc))
c
c      a rotation of THETA radians about the Yprime-axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      locold and ephi inside the following loop. 
c
      do l=0,nterms
        do mp=-l,l
          do idim=1,nd
            marray1(idim,l,mp) = locold(idim,l,mp)*ephi(mp)
          enddo
        enddo
      enddo

      if(nterms.ge.30) then
        call rotviaproj(nd,theta,nterms,nterms,nterms,marray1,
     1     ldc,marray,ldc)
      else
        call rotviarecur(nd,theta,nterms,nterms,nterms,marray1,
     1     ldc,marray,ldc)
      endif
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
       call l3dlocloczshift(nd,sc1,marray,ldc,nterms,sc2,marray1,
     1           nterms2,nterms2,rshift,dc,lda) 
c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      if(nterms.ge.30) then
        call rotviaproj(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1     ldc,marray,ldc)
      else
        call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1     ldc,marray,ldc)
      endif
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
        do m=-l,l
          do idim=1,nd
            local(idim,l,m)=local(idim,l,m)+ephi(-m)*marray(idim,l,m)
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
      subroutine l3dlocloczshift(nd,scale,locold,lmp,nterms,scale2,
     1  local,lmpn,nterms2,zshift,dc,lda) 
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
c     INPUT:
c
c     nd       : number of local expansions
c     scale    : scaling parameter for locold
c     locold   : coefficients of original multipole exp.
c     lmp      : leading dim of locold (may be a work array)
c     nterms   : number of terms in the orig. expansion
c
c     scale2   : scaling parameter for output expansion (local)
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis (assumed positive)
c
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c-----------------------------------------------------------------------
      implicit none
      integer idim,nd,lda
      integer nmax,nterms,nterms2,nquad,ier,l,lw,m,jnew,knew
      integer lmp,lmpn,ll
      double precision  zshift,d
      double precision  scale,scale2
      double precision, allocatable :: fr(:),rscpow(:)
      double complex locold(nd,0:lmp,-lmp:lmp)
      double complex local(nd,0:lmpn,-lmpn:lmpn)
      double precision dc(0:lda,0:lda)
C
C----- shift along z-axis 
C
      nmax = max(nterms,nterms2)
      
      allocate(fr(0:nmax),rscpow(nmax))
c
c
      d = zshift
      fr(0) = 1.0d0
      d = d*scale
      fr(1) = d
      do l=2,nterms
         fr(l) = fr(l-1)*d
      enddo

      d = scale/scale2
      rscpow(1) = d
      do l=2,nterms
        rscpow(l) = rscpow(l-1)*d
      enddo

      do jnew = 0,nterms2
        do knew = -jnew,jnew
          do idim=1,nd
            local(idim,jnew,knew) = locold(idim,jnew,knew)
          enddo

          do l=1,nterms-jnew
            ll = l+jnew
            do idim=1,nd
              local(idim,jnew,knew)=local(idim,jnew,knew)+
     1         locold(idim,ll,knew)*fr(l)*dc(ll+knew,l)*dc(ll-knew,l)
            enddo
          enddo
        enddo
      enddo

      do jnew = 0,nterms2
        do knew = -jnew,jnew
          do idim=1,nd
            local(idim,jnew,knew)=local(idim,jnew,knew)*rscpow(jnew)
          enddo
        enddo
      enddo

      return
      end
c
c--------------------------------
C
      subroutine d3tata(rscale,x0y0z0,mpole,nterms,
     1            rscale2,xnynzn,mpolen,nterms2)
      implicit double precision (a-h,o-z)
      integer nterms
      double precision x0y0z0(3),xnynzn(3),pp(0:nterms,0:nterms)
      double precision zdiff(3)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex mpolen(0:nterms2,-nterms2:nterms2)
c
      double precision powers(0:60)
      double precision rx,ry,rz,proj,rr,d,dd,ctheta
      double precision cs(0:120,-120:120),fact(0:102),cscale
      double complex cscale2
      double complex ephi(-60:60),imag,ephi1
      integer l,m,lnew,mnew,ll,mm
      data imag/(0.0d0,1.0d0)/
c
      do 50 l = 0,nterms2
         do 40 m = -l,l
            mpolen(l,m) = 0.0d0
40       continue
50    continue
      d = 1.0d0
      fact(0) = d
      do 60 l = 1,2*nterms
         d = d*dsqrt(l+0.0D0)
         fact(l) = d
60    continue
      cs(0,0) = 1.0d0
      do 80 l = 1,nterms
         do 70 m = 0,l
            cs(l,m) =  ((-1)**l)/( fact(l-m)*fact(l+m) )
            cs(l,-m) = cs(l,m)
 70      continue
 80   continue
CCC      print *, ' cs(l,m) ARRAY ',cs
c
      zdiff(1) = x0y0z0(1) - xnynzn(1)
      zdiff(2) = x0y0z0(2) - xnynzn(2)
      zdiff(3) = x0y0z0(3) - xnynzn(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
cc      print *, ' created rr,proj, etc.'
cc      call prin2(' rscale is *',rscale,1)
      dd = d*rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do 100 l = 2,nterms+1
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
100   continue
c
      call ylgndr(nterms,ctheta,pp)
c
C
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
      do 600 j = 0,nterms2
         do 500 k = -j,j
            do 300 l = j,nterms
               do 200 m = -l,l
                  ll = l-j
                  mm = m-k
                  if (abs(mm).gt.ll) goto 111
                  cscale = powers(ll)*cs(j,k)*cs(ll,mm)/cs(l,m)
                  cscale = cscale/dsqrt(2*ll+1.0D0)
                  cscale = cscale*(-1)**ll
                  cscale = cscale*(rscale/rscale2)**j
                  if ( m*mm .lt. 0) then
                     cscale = cscale*(-1)**mm
                  endif
                  if (m*mm .ge. 0)  then
                     if ( abs(m) .le. abs(mm) ) 
     1               cscale = cscale*(-1)**k
                  endif
                  cscale2 = cscale
                  if (mm .eq. 0) then
                     mpolen(j,k) = mpolen(j,k) +
     1               pp(ll,0)*cscale2*mpole(l,m)
                  else  if (mm .gt. 0) then
                     mpolen(j,k) = mpolen(j,k)+ cscale2*
     1               pp(ll,mm)*ephi(mm)*mpole(l,m)
                  else
                     mpolen(j,k) = mpolen(j,k)+ cscale2*
     1               pp(ll,-mm)*ephi(mm)*mpole(l,m)
                  endif
111               continue
200            continue
300         continue
500      continue
600   continue
      return
      end
c
c-----------------------------------------------------------------
      subroutine d3tataf(rscale,x0y0z0,mpole,nterms,
     1            rscale2,xnynzn,mpolen,nterms2,cs,ns,wlege,nlege)
      implicit double precision (a-h,o-z)
      integer nterms
      double precision x0y0z0(3),xnynzn(3),pp(0:nterms,0:nterms)
      double precision zdiff(3)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex mpolen(0:nterms2,-nterms2:nterms2)
c
      double precision powers(0:60),rpow(0:100)
      double precision rx,ry,rz,proj,rr,d,dd,ctheta
      double precision cs(0:ns,-ns:ns),cscale
      double precision wlege(*)
      double complex cscale2
      double complex ephi(-60:60),imag,ephi1
      integer l,m,lnew,mnew,ll,mm
      integer onepow(-1000:1000)
      data imag/(0.0d0,1.0d0)/
c
      do 50 l = 0,nterms2
         do 40 m = -l,l
            mpolen(l,m) = 0.0d0
40       continue
50    continue

cc      fact(0) = d
cc      do 60 l = 1,2*nterms
cc         d = d*dsqrt(l+0.0D0)
cc         fact(l) = d
cc60    continue
cc      cs(0,0) = 1.0d0
cc      do 80 l = 1,nterms
cc         do 70 m = 0,l
cc            cs(l,m) =  ((-1)**l)/( fact(l-m)*fact(l+m) )
cc            cs(l,-m) = cs(l,m)
cc 70      continue
cc 80   continue
CCC      print *, ' cs(l,m) ARRAY ',cs
c
      zdiff(1) = x0y0z0(1) - xnynzn(1)
      zdiff(2) = x0y0z0(2) - xnynzn(2)
      zdiff(3) = x0y0z0(3) - xnynzn(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
cc      print *, ' created rr,proj, etc.'
cc      call prin2(' rscale is *',rscale,1)
      dd = -d*rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      dd2 = rscale/rscale2
      rpow(0) = 1
      rpow(1) = dd2

      do 100 l = 2,nterms+1
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
         rpow(l) = rpow(l-1)*dd2
100   continue

      do 110 l=0,nterms+1
          powers(l) = powers(l)/dsqrt(2*l+1.0d0)
110   continue     

      onepow(0) = 1
      mone = -1

      do 111 l=1,2*nterms+1

      onepow(l) = onepow(l-1)*mone
      onepow(-l) = onepow(l)

111   continue
c
      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
C
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
      do 600 j = 0,nterms2
         do 500 k = -j,j
            do 300 l = j,nterms
               do 200 m = -l,l
                  ll = l-j
                  mm = m-k
                  if (abs(mm).gt.ll) goto 112

                  mmm = (abs(m)-abs(k)-abs(m-k))/2
                  cscale = powers(ll)*cs(j,k)*cs(ll,mm)/cs(l,m)
                  cscale = cscale*rpow(j)*onepow(mmm)
cc                  if ( m*mm .lt. 0) then
cc                     cscale = cscale*onepow(mm)
cc                  endif
cc                  if (m*mm .ge. 0)  then
cc                     if ( abs(m) .le. abs(mm) ) 
cc     1               cscale = cscale*onepow(k)
cc                  endif

cc                  if (mm .eq. 0) then
cc                     mpolen(j,k) = mpolen(j,k) +
cc     1               pp(ll,0)*cscale*mpole(l,m)

                  if (mm .ge. 0) then
                     mpolen(j,k) = mpolen(j,k)+ cscale*
     1               pp(ll,mm)*ephi(mm)*mpole(l,m)
                  else
                     mpolen(j,k) = mpolen(j,k)+ cscale*
     1               pp(ll,-mm)*ephi(mm)*mpole(l,m)
                  endif
112               continue
200            continue
300         continue
500      continue
600   continue
      return
      end
c

