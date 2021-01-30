c
c
c    Subroutine library for shifting via projection     
c
C***********************************************************************
      subroutine h3drescalemp(nd,nterms,lmp,mpole,
     1           radius,zk0,scale,fhs,fhder)
C***********************************************************************
C
C     This subroutine rescales a spherical harmonic expansion
C     on the surface by h_n(zk0 * radius), converting 
C     a surface function to the corresponding Hankel expansion
C     in the exterior.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nd     = number of expansions
C           nterms = order of spherical harmonic expansion
C           mpole = coefficients of s.h. expansion
C           radius = sphere radius
C           zk0 = Helmholtz parameter
C           scale = scale parameter for expansions
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpole = rescaled by 1/h_n(zk0* radius) 
C
C---------------------------------------------------------------------
      implicit none
      integer(8) lmp
      integer(8) nterms,ier,nd,idim
      integer(8) l,m,jj,kk
      integer(8) lwfhs,ifder
      real *8 radius,scale
      complex *16 mpole(nd,0:lmp,-lmp:lmp)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 fhs(0:nterms),fhder(0:nterms)
      complex *16 zk0,z
      data imag/(0.0d0,1.0d0)/
C
      z = zk0*radius
      ifder = 0
      call h3dall(nterms,z,scale,fhs,
     1             ifder,fhder)
      do l=0,nterms
         do m=-l,l
            zmul = 1/fhs(l)
            do idim=1,nd
	          mpole(idim,l,m) = mpole(idim,l,m)*zmul
            enddo
         enddo
      enddo
      return
      end
c
C
C
C***********************************************************************
      subroutine h3drescaleloc(nd,nterms,lmp,local,localn,
     1           radius,zk0,scale,fjs,fjder)
C***********************************************************************
C
C     This subroutine takes as input the potential and its normal
C     derivative on a sphere of radius RADIUS and returns the 
C     j-expansion coefficients consist with the data (in a least
C     squares sense). That is 
C
C           phi  = sum local(n,m) j_n Y_n^m  ->
C           phin = sum local(n,m) zk0 *j_n' Y_n^m and
C
C           local(n,m) * j_n        = phi_n^m
C           local(n,m) * j_n' *zk0  = phin_n^m
C
C     The 1x1 normal equations are:
C
C           local(n,m)* ( j_n**2 + (zk0* j_n')**2) = 
C                 
C                          j_n * phi_n^m + zk0 * j_n' * phin_n^m.
C                   
C---------------------------------------------------------------------
C     INPUT:
C
C           nd     = number of local expansions
C           nterms = order of spherical harmonic expansion
C           local = coefficients of s.h. expansion of phi
C           localn = coefficients of s.h. expansion of dphi/dn
C           radius = sphere radius
C           zk0 = Helmholtz parameter
C           scale = scale parameter for expansions
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = computed as above.
C
C---------------------------------------------------------------------
      implicit none
      integer(8) nterms,ier,nd,idim
      integer(8) lmp
      integer(8) l,m,jj,kk
      integer(8) lwfhs,ifder
      real *8 radius,scale
      complex *16 fjs(0:nterms)
      complex *16 fjder(0:nterms)
      complex *16 local(nd,0:lmp,-lmp:lmp)
      complex *16 localn(nd,0:lmp,-lmp:lmp)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 zk0,z,zh,zhn
      data imag/(0.0d0,1.0d0)/
C
C
      z = zk0*radius
      ifder = 1
      call besseljs3d(nterms,z,scale,fjs,ifder,fjder)
c
      do l=0,nterms
        do m=-l,l
          zh = fjs(l)
          zhn = fjder(l)*zk0
          zmul = 1/(zh*zh + zhn*zhn)
          do idim=1,nd
            local(idim,l,m) = (zh*local(idim,l,m) + 
     1           zhn*localn(idim,l,m))*zmul
          enddo
        enddo
      enddo

      return
      end
c
C
C
C
C
C
C
C
C***********************************************************************
      subroutine h3dlocevalsphere(nd,local,zk,scale,zshift,radius,
     1           nterms,nterms2,lmp,ynm,ynmd,phitemp,phitempn,
     2           nquad,xnodes,fjs,fjder,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates a local expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd       : number of local expansions
C     local    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     nterms2  : number of terms in exp on target sphere
C     lmp      : dimension param for local expansion
C     ynm      : storage for Ynm out to nterms.
C     nquad    : number of quadrature nodes
C                              on target sphere is nquad*nquad.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phitemp(k,i,j)  : jth mode of phi at ith quad node, for kth expansion 
C     phitempn(k,i,j) : jth mode of phi at ith quad node, for kth expansion
C
C***********************************************************************
      implicit none
      integer(8) lmp,ier
      integer(8) nterms,nd,idim,nterms2,nquad
      integer(8) l,m,jnew,knew,jj,n,ifder,iffld
      real *8 radius,scale
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 ctheta,cthetaj,rj,rn,stheta,pi,sthetaj
      real *8 thetan
      complex *16 local(nd,0:lmp,-lmp:lmp)
      complex *16 phitemp(nd,nquad,-nterms2:nterms2)
      complex *16 phitempn(nd,nquad,-nterms2:nterms2)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta
      complex *16 ephi1,ephik,ephi,fjs(0:nterms),fjder(0:nterms),ztmp1
      complex *16 ut1,ut2,ut3
      real *8 rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms)
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      ier = 0
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 1
      ifder = 1
      do jj=1,nquad
        do m=-nterms2,nterms2
          do idim=1,nd
            phitemp(idim,jj,m) = 0.0d0
            phitempn(idim,jj,m) = 0.0d0
          enddo
        enddo
      enddo
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquad
        ctheta = xnodes(jj)
        stheta = dsqrt(1.0d0 - ctheta**2)
        rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
        rj = dsqrt(rj)
        cthetaj = (zshift+radius*ctheta)/rj
        sthetaj = dsqrt(1.0d0-cthetaj**2)
        rn = sthetaj*stheta + cthetaj*ctheta
        thetan = (cthetaj*stheta - sthetaj*ctheta)/rj
        z = zk*rj
        call ylgndr2sf(nterms,cthetaj,ynm,ynmd,rat1,rat2)
        call besseljs3d(nterms,z,scale,fjs,ifder,fjder)
	    do n = 0,nterms
	      fjder(n) = fjder(n)*zk
        enddo
	    do n = 1,nterms
	      do m = 1,n
	        ynm(n,m) = ynm(n,m)*sthetaj
          enddo
        enddo

        ztmp1 = fjder(0)*rn
        do idim=1,nd
          phitemp(idim,jj,0) = local(idim,0,0)*fjs(0)
          phitempn(idim,jj,0) = local(idim,0,0)*ztmp1
        enddo

        do n=1,nterms
          ztmp1 = fjs(n)*ynm(n,0)
          ut1 = fjder(n)*rn
          ut2 = fjs(n)*thetan
          ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj

          do idim=1,nd
            phitemp(idim,jj,0) = phitemp(idim,jj,0) +
     1                  local(idim,n,0)*ztmp1
            phitempn(idim,jj,0) = phitempn(idim,jj,0)+
     1          ut3*local(idim,n,0)
          enddo
c
          do m=1,min(n,nterms2)
c
            ztmp1 = fjs(n)*ynm(n,m)
            ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)

            do idim=1,nd
              phitemp(idim,jj,m) = phitemp(idim,jj,m) + 
     1                local(idim,n,m)*ztmp1
              phitemp(idim,jj,-m) = phitemp(idim,jj,-m) +
     1                local(idim,n,-m)*ztmp1
              phitempn(idim,jj,m) = phitempn(idim,jj,m) +
     1             ut3*local(idim,n,m)
              phitempn(idim,jj,-m) = phitempn(idim,jj,-m) + 
     1             ut3*local(idim,n,-m)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
C
C
C
C
c
c
c
C***********************************************************************
      subroutine h3dmpevalsphere(nd,mpole,zk,scale,zshift,
     1           radius,nterms,lmp,ynm,ynmd,phitemp,phitempn,
     2           nquad,xnodes,fhs,fhder,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd       : number of multipole expansions
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension param for mpole
C     ynm      : storage for assoc Legendre functions
C     ynmd     : storage for derivs of assoc Legendre functions
C     nquad    : number of quadrature nodes in theta
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phitemp  : nquad by (-nterms,nterms)
C     phitempn : nquad by (-nterms,nterms)
C
C---------------------------------------------------------------------
      implicit none
      integer(8) lmp,n
      integer(8) nterms,nd,idim
      integer(8) l,m,jnew,knew,nquad
      integer(8) ifder,iffld,i,jj
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 radius,scale
      real *8 ctheta,cthetaj,pi,rj,rn,stheta,sthetaj,thetan
      complex *16 mpole(nd,0:lmp,-lmp:lmp)
      complex *16 phitemp(nd,nquad,-nterms:nterms)
      complex *16 phitempn(nd,nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta,ut1,ut2
      complex *16 ut3,ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
      complex *16 ztmp1
      real *8 rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms)
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 1
      ifder = 1
      do jj=1,nquad
        do m=-nterms,nterms
          do idim=1,nd
            phitemp(idim,jj,m) = 0.0d0
            phitempn(idim,jj,m) = 0.0d0
          enddo
        enddo
      enddo
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquad
	    ctheta = xnodes(jj)
	    stheta = dsqrt(1.0d0 - ctheta**2)
        rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
        rj = dsqrt(rj)
	    cthetaj = (zshift+radius*ctheta)/rj
	    sthetaj = dsqrt(1.0d0 - cthetaj**2)
	    rn = sthetaj*stheta + cthetaj*ctheta
	    thetan = (cthetaj*stheta - ctheta*sthetaj)/rj
	    z = zk*rj
	    call ylgndr2sf(nterms,cthetaj,ynm,ynmd,rat1,rat2)
	    call h3dall(nterms,z,scale,fhs,ifder,fhder)
        do i = 0,nterms
	      fhder(i) = fhder(i)*zk
        enddo
        do n=1,nterms
	      do m = 1,n
  	        ynm(n,m) = ynm(n,m)*sthetaj
          enddo
        enddo

        ztmp1 = fhder(0)*rn

        do idim=1,nd
          phitemp(idim,jj,0) = mpole(idim,0,0)*fhs(0)
          phitempn(idim,jj,0) = mpole(idim,0,0)*ztmp1
        enddo
	    do n=1,nterms
          ztmp1 = fhs(n)*ynm(n,0)
	      ut1 = fhder(n)*rn
	      ut2 = fhs(n)*thetan
	      ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj

          do idim=1,nd
            phitemp(idim,jj,0) = phitemp(idim,jj,0) +
     1                  mpole(idim,n,0)*ztmp1
	        phitempn(idim,jj,0) = phitempn(idim,jj,0) + 
     1          ut3*mpole(idim,n,0)
          enddo 
          do m=1,n
	        ztmp1 = fhs(n)*ynm(n,m)
	        ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
            do idim=1,nd
  	          phitemp(idim,jj,m) = phitemp(idim,jj,m) +
     1                   mpole(idim,n,m)*ztmp1
	          phitemp(idim,jj,-m) = phitemp(idim,jj,-m) +
     1                   mpole(idim,n,-m)*ztmp1
	          phitempn(idim,jj,m) = phitempn(idim,jj,m) + 
     1            ut3*mpole(idim,n,m)
	          phitempn(idim,jj,-m) = phitempn(idim,jj,-m) + 
     1             ut3*mpole(idim,n,-m)
             enddo
	      enddo
	    enddo
      enddo
c
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine h3dprojloc(nd,nterms,ldl,nquadn,ntold,xnodes,wts,
     1           phitemp,phitempn,local,local2,ynm,rat1,rat2)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nd = number of expansions
C           nterms = order of spherical harmonic expansion
C           ldl = order of spherical harmonic expansion
C           nquadn = number of quadrature nodes in polar angle (theta).
C           ntold = number of modes in azimuthal direction.
C           xnodes = quad nodes in theta (polar angle) - nquadn of them
C           wts  = quad weights in theta (polar angle)
C           phitemp, phitempn = tabulated function and normal deriv
C                    phitemp(i,j) = jth mode of phi at ith quad node.
C                    phivaln(i,j) = jth mode of phi at ith quad node.
C
C           marray  = workspace 
C           ynm     = workspace for assoc Legendre functions
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local  = coefficients of s.h. expansion for phi
C           local2 = coefficients of s.h. expansion for dphi/dn
C
C    NOTE:
C
C    yrecursion.f produces Ynm with a nonstandard scaling:
C    (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C    is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C   In the first loop below, you see
C
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C---------------------------------------------------------------------
      implicit none
      integer(8) nterms,nquadn,nquadm,ier,nd,idim
      integer(8) l,m,jj,kk,ntold,ldl,n
      real *8 wts(nquadn),xnodes(nquadn)
      real *8 ynm(0:nterms,0:nterms)
      real *8 cthetaj,pi
      complex *16 zk
      complex *16 local(nd,0:ldl,-ldl:ldl)
      complex *16 local2(nd,0:ldl,-ldl:ldl)
      complex *16 phitemp(nd,nquadn,-ntold:ntold)
      complex *16 phitempn(nd,nquadn,-ntold:ntold)
      complex *16 ephi,imag,emul,sum,emul1
      complex *16, allocatable :: zmul(:),zmul2(:)
      real *8 rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)

      allocate(zmul(nd),zmul2(nd))
c
c     initialize local exp to zero
c
      do l = 0,ldl
        do m = -l,l
          do idim=1,nd
            local(idim,l,m) = 0.0d0
            local2(idim,l,m) = 0.0d0
          enddo
        enddo
      enddo
c
c     get local exp
c
c
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquadn
        cthetaj = xnodes(jj)
        call ylgndrf(nterms,cthetaj,ynm,rat1,rat2)
        do m=-ntold,ntold
          do idim=1,nd 
            zmul(idim) = phitemp(idim,jj,m)*wts(jj)/2.0d0
            zmul2(idim) = phitempn(idim,jj,m)*wts(jj)/2.0d0
          enddo
          do l=abs(m),nterms
            do idim=1,nd
              local(idim,l,m) = local(idim,l,m) + 
     1            zmul(idim)*ynm(l,abs(m))
              local2(idim,l,m) = local2(idim,l,m) + 
     1            zmul2(idim)*ynm(l,abs(m))
            enddo
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
C
C
c
c
c
c
