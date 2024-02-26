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
      call cart2polar(rvec,d,theta,phi)

      ldc = max(nterms,nterms2)
      allocate(marray1(nd,0:ldc,-ldc:ldc))
      allocate(marray(nd,0:ldc,-ldc:ldc))
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
     1        ldc,marray,ldc)
      else
        call rotviarecur(nd,-theta,nterms2,nterms2,nterms2,marray1,
     1        ldc,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
        do m=-l,l
          do idim=1,nd
            mpolen(idim,l,m)=mpolen(idim,l,m)+ephi(-m)*marray(idim,l,m)
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
c     origin to a multipole expansion centered at (0,0,zshift).
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
      integer l0,nmax,nterms,nterms2,nquad,ier,lmp,lmpn,ldc,nd
      integer lca
      integer l,m,jnew,knew,i,idim
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
      d = -zshift/scale
      fr(1) = d
      do l=2,nmax
         fr(l) = fr(l-1)*d
      enddo

      d = scale/scale2
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
      call cart2polar(rvec,d,theta,phi)
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
     1      sc2,ldc,nterms2,rshift,dc,lca)

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
      d = d*scale
      fr(1) = fr(0)*d
      do l=2,2*nmax
        fr(l) = fr(l-1)*d
      enddo

      d = scale2/scale
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
      integer idim,lda,nd
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

      call cart2polar(rvec,d,theta,phi)

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
     1           ldc,nterms2,rshift,dc,lda)
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
      d = d/scale
      fr(1) = d
      do l=2,nterms
         fr(l) = fr(l-1)*d
      enddo

      d = scale2/scale
      rscpow(1) = d
      do l=2,nterms2
        rscpow(l) = rscpow(l-1)*d
      enddo

      do jnew = 0,nterms2
        do knew = -jnew,jnew
          if(jnew.le.nterms) then
            do idim=1,nd
              local(idim,jnew,knew) = locold(idim,jnew,knew)
            enddo
          else
            do idim=1,nd
              local(idim,jnew,knew) = 0.0d0 
            enddo
          endif

          do l=1,nterms-jnew
            ll = l+jnew
            do idim=1,nd
              local(idim,jnew,knew)=local(idim,jnew,knew)+
     1         locold(idim,ll,knew)*fr(l)*dc(ll+knew,l)*dc(ll-knew,l)
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
c
