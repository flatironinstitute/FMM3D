cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2011-06-21 12:54:41 -0400 (Tue, 21 Jun 2011) $
c    $Revision: 2105 $
c
c
c    Subroutine library for shifting via projection     
c
C***********************************************************************
      subroutine h3dprojloc(nterms,ldl,nquad,xnodes,wts,phival,local,
     1           ynm,w,lw,ier)
C***********************************************************************
C     Purpose:
C
C           This subroutine projects a function given on the unit sphere
C           onto spherical harmonics.
C           It is intended for use with a library of special 
C           functions in which Ynm does NOT include the 
C           1/sqrt(4*pi) normalization. Thus, the Ynm are orthogonal
C           but not orthonormal -> the inner product value is 4*pi.
c
C           phival(theta,phi) tabulated at nquad*nquad points:
C           Legendre nodes x_j = cos(theta_j) in polar angle and 
C           equispaced in the azimuthal direction.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl =    defines leading dimension of local
C           nquad  = number of quadrature nodes in each direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C           ynm     = workspace of dimension (0:nterms,0:nterms)
c                     used to get Ynm values.
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C           ier = error flag 
C                 0 normal return
C                 8 insufficient memory in w.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,ier,ldl
      integer l,m,jj,kk
      real *8 w(lw)
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phival(nquad,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      imarray = 1
      lmarray = 2*nquad*(2*nterms+1)+7
C
      lused = imarray+lmarray
C
      if (lused.gt.lw) then
         ier=6
	 return
      endif
c
      call h3dprojloc0(nterms,ldl,nquad,xnodes,wts,phival,local,
     1           w(imarray),ynm)
      return
      end
c
C***********************************************************************
      subroutine h3dprojloc0(nterms,ldl,nquad,xnodes,wts,phival,local,
     1           marray,ynm)
C***********************************************************************
C
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl    = dimension parameter for local expansion
C           nquad  = number of quadrature nodes in each direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           marray  = workspace of dimension (nquad,-nterms:nterms)
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
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
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phival(nquad,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 marray(nquad,-nterms:nterms)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
         enddo
      enddo
c
c     create marray (intermediate array)
c
      emul1 = cdexp(imag*2*pi/nquad)
      emul = cdexp(-imag*2*nterms*pi/nquad)
      do m=-nterms,nterms
cc	    emul = cdexp(imag*m*2*pi/nquad)
	 do jj=1,nquad
	    sum = 0
	    ephi = emul
	    do kk = 1,nquad
               sum = sum + phival(jj,kk)*dconjg(ephi)
	       ephi = ephi*emul
            enddo
ccc	    marray(jj,m) = sum*2*pi/nquad
	    marray(jj,m) = sum/(2*nquad)
         enddo
	 emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquad
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
ccc	 call prinf(' nterms is *',nterms,1)
ccc	 call prinf(' ldl is *',ldl,1)
ccc	 call prinf(' l is *',l,1)
         do m=-nterms,nterms
	    zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
C***********************************************************************
      subroutine h3drescale(nterms,local,radius,zk0,scale,w,lw,ier)
C***********************************************************************
C
C     This subroutine rescales a spherical harmonic expansion
C     on the surface by j_n(zk0 * radius), converting 
C     a surface function to the corresponding Bessel expansion
C     in the interior.
C     It is intended for use with a specific library of special 
C     functions in which local expansions are also scaled by 
C     i* zk0, so this factor is scaled out here for consistency.
c
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           local = coefficients of s.h. expansion
C           radius = sphere radius
C           zk0 = Helmholtz parameter
C           scale  = scale parameter for J expansion.
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = rescaled by 1/(j_n(zk0* radius)*i*zk0) 
C           ier = error flag 
C                 0 normal return
C                 4 insufficient memory in w.
C                 8 lwfjs insufficient in calling jfuns3d.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ier,l,m,jj,kk,lwfjs
      real *8 w(lw)
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 zk0,z
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      lwfjs = nterms+1000
c
      ifjs = 1
      lfjs = 2*(lwfjs+1)+7
C
      ifjder = ifjs+lfjs
      lfjder = 2*(nterms+1)+7
C
      iiscale = ifjder+lfjder
      liscale = (lwfjs+1) +7
C
      lused = iiscale+liscale
      if (lused.gt.lw) then
         ier=4
	 return
      endif
c
      z = zk0*radius
      ifder = 0
      call jfuns3d(ier1,nterms,z,scale,w(ifjs),
     1             ifder,w(ifjder),lwfjs,w(iiscale),ntop)
      if (ier1.eq.8) then
         ier = 8
	 return
      endif
ccc      call prin2(' fjs is *',w(ifjs),2*nterms)
      do l=0,nterms
         do m=-l,l
	    zmul = dcmplx(w(ifjs+2*l),w(ifjs+2*l+1))
	    local(l,m) = local(l,m)/zmul
         enddo
      enddo
      return
      end
c
C***********************************************************************
      subroutine h3dmpevalsphereslow(mpole,zk,scale,phival,zshift,
     1           radius,nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on an
C     nquad x nquadm spectral grid on a target sphere
C     along the z axis a distance zshift.
C   
C     TESTING CODE  ---- NOT USED
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     nquad    : number of quadrature nodes in theta direction
C     nquadm   : number of quadrature nodes in phi direction
C
C                total number of nodes on target sphere is nquad*nquadm
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 phival(nquad,nquadm)
      complex *16 pot,fld(3), zk
      real *8 w(100000)
C
      integer l,m,jnew,knew
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      lw = 100000
      iffld = 1
cc      pot = dcmplx(1.0d0,1.0d0)
      do jj=1,nquad
      do kk=1,nquadm
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
	 phi = 2*pi*kk/nquadm
	 cosphi = dcos(phi)
	 sinphi = dsin(phi)
	 targ(1) = radius*stheta*cosphi
	 targ(2) = radius*stheta*sinphi
	 targ(3) = zshift + radius*ctheta
	 rn1 = stheta*cosphi
	 rn2 = stheta*sinphi
	 rn3 = ctheta
         call h3dmpeval(zk,scale,center,mpole,nterms,targ,
     1        pot,iffld,fld,w,lw,lused,ier)
         phival(jj,kk) = pot
ccc         phival(jj,kk) = -(fld(1)*rn1+fld(2)*rn2+fld(3)*rn3)
      enddo
      enddo
      return
      end
C
C
C
C
C***********************************************************************
      subroutine h3dmpevalsphere(mpole,zk,scale,phival,zshift,radius,
     1           nterms,lmp,ynm,phitemp,nquad,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C     TESTING CODE  ---- NOT USED
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for mpole
C     ynm      : storage for Ynm
C     phitemp  : intermediate storage for p^3 method
C     nquad    : number of quadrature nodes
C                              on target sphere is nquad*nquad.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C
C     NOTE: fhs and fhder hardwired at 2000 here.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phival(nquad,nquad)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi,fhs(0:2000),fhder(0:2000)
      complex *16 ephi1,ephik
C
      integer l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 0
      ifder = 0
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
      enddo
      enddo
      do jj=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
         rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
         rj = dsqrt(rj)
	 cthetaj = (zshift+radius*ctheta)/rj
	 z = zk*rj
	 call ylgndr(nterms,cthetaj,ynm)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do m=-nterms,nterms
	    mabs = abs(m)
	    do n=mabs,nterms
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*fhs(n)*ynm(n,mabs)
	    enddo
	 enddo
      enddo
      ephi1 = cdexp(2*pi*imag/nquad)
      ephik = 1.0d0
      do jj = 1,nquad
      do kk = 1,nquad
         phival(jj,kk) = 0.0d0
         ephik = ephik*ephi1
         ephi = ephik**(-nterms)
	 do m = -nterms,nterms
cc	    phi = 2*pi*kk*m/nquad
cc	    ephi = cdexp(imag*phi)
            phival(jj,kk) = phival(jj,kk) + phitemp(jj,m)*ephi
	    ephi = ephi*ephik
         enddo
      enddo
      enddo
      return
      end
C
C
C
C
C
C***********************************************************************
      subroutine h3dlocevalsphere(local,zk,scale,phival,zshift,radius,
     1           nterms,nterms2,lmp,ynm,phitemp,nquad,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a local expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C     TESTING CODE - OBSOLETE.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C     phitemp  : storage for phitemp in p^3 algorithm.
C     nquad    : number of quadrature nodes
C                              on target sphere is nquad*nquad.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer iscale(0:2000)
      integer nterms
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquad,nquad)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi1,ephik,ephi,fjs(0:2000),fjder(0:2000)
C
      integer l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 0
      ifder = 0
      lwfjs = 2000
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
      enddo
      enddo
      do jj=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
         rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
         rj = dsqrt(rj)
	 cthetaj = (zshift+radius*ctheta)/rj
	 z = zk*rj
	 call ylgndr(nterms,cthetaj,ynm)
	 call jfuns3d(jer,nterms,z,scale,fjs,ifder,fjder,
     1        lwfjs,iscale,ntop)
         do m=-nterms2,nterms2
	    mabs = abs(m)
	    do n=mabs,nterms
	       phitemp(jj,m) = phitemp(jj,m) +
     1                local(n,m)*fjs(n)*ynm(n,mabs)
	    enddo
	 enddo
      enddo
c
      ephi1 = cdexp(2*pi*imag/nquad)
      ephik = 1.0d0
      do jj = 1,nquad
      do kk = 1,nquad
         phival(jj,kk) = 0.0d0
         ephik = ephik*ephi1
         ephi = ephik**(-nterms2)
	 do m = -nterms2,nterms2
cc	    phi = 2*pi*kk*m/nquad
cc	    ephi = cdexp(imag*phi)
            phival(jj,kk) = phival(jj,kk) + phitemp(jj,m)*ephi
	    ephi = ephi*ephik
         enddo
      enddo
      enddo
      return
      end
c
C
C***********************************************************************
      subroutine h3drescalemp(nterms,lmp,mpole,
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
C           ier = error flag 
C                 0 normal return
C                 4 insufficient memory in w.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ier
      integer l,m,jj,kk
      integer lwfhs
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 fhs(0:nterms),fhder(0:nterms)
      complex *16 zk0,z
      data imag/(0.0d0,1.0d0)/
C
      z = zk0*radius
      ifder = 0
      call h3dall(nterms,z,scale,fhs,
     1             ifder,fhder)
ccc      call prin2(' fhs is *',w(ifhs),2*nterms)
      do l=0,nterms
         do m=-l,l
	    zmul = fhs(l)
	    mpole(l,m) = mpole(l,m)/zmul
         enddo
      enddo
      return
      end
c
C
C***********************************************************************
      subroutine h3drescalempstab(nterms,lmp,mpole,mpn,
     1           radius,zk0,scale,fhs,fhder)
C***********************************************************************
C
C     This subroutine takes as input the potential and its normal
C     derivative on a sphere of radius RADIUS and returns the 
C     h-expansion coefficients consist with the data (in a least
C     squares sense). That is 
C
C           phi  = sum mpole(n,m) h_n Y_n^m  ->
C           phin = sum mpole(n,m) zk0 *h_n' Y_n^m and
C
C           mpole(n,m) * h_n        = phi_n^m
C           mpole(n,m) * h_n' *zk0  = phin_n^m
C
C     The 1x1 normal equations are:
C
C           mpole(n,m)* ( h_n**2 + (zk0* h_n')**2) = 
C                 
C                          h_n * phi_n^m + zk0 * h_n' * phin_n^m.
C                   
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           mpole = coefficients of s.h. expansion for phi
C           mpn   = coefficients of s.h. expansion for dphi/dn
C           radius = sphere radius
C           zk0 = Helmholtz parameter
C           scale = scale parameter for expansions
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpole = computed as above.
C           ier = error flag 
C                 0 normal return
C                 4 insufficient memory in w.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ier
      integer l,m,jj,kk
      integer lwfhs
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 mpn(0:lmp,-lmp:lmp)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 fhs(0:nterms),fhder(0:nterms)
      complex *16 zk0,z,zh,zhn
      data imag/(0.0d0,1.0d0)/
C
      z = zk0*radius
      ifder = 1
      call h3dall(nterms,z,scale,fhs,
     1             ifder,fhder)
      do l=0,nterms
         do m=-l,l
	    zh = fhs(l)
	    zhn = fhder(l)
	    zmul = zh*zh + zhn*zhn
	    mpole(l,m) = (zh*mpole(l,m) + zhn*mpn(l,m))/zmul
         enddo
      enddo
      return
      end
c
C
C***********************************************************************
      subroutine h3drescalestab(nterms,lmp,local,localn,
     1           radius,zk0,scale,fjs,fjder,iscale,lwfjs,ier)
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
C           nterms = order of spherical harmonic expansion
C           local = coefficients of s.h. expansion of phi
C           localn = coefficients of s.h. expansion of dphi/dn
C           radius = sphere radius
C           zk0 = Helmholtz parameter
C           scale = scale parameter for expansions
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = computed as above.
C           ier = error flag 
C                 0 normal return
C                 4 insufficient memory in w.
C                 8 lwfjs insufficient in calling jfuns3d.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ier
      integer l,m,jj,kk
      integer lwfhs
      complex *16 fjs(0:lwfjs)
      complex *16 fjder(0:nterms)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 localn(0:lmp,-lmp:lmp)
      complex *16 ephi,imag,emul,sum,zmul
      complex *16 zk0,z,zh,zhn
      data imag/(0.0d0,1.0d0)/
C
C
      z = zk0*radius
      ifder = 1
      call jfuns3d(ier1,nterms,z,scale,fjs,
     1             ifder,fjder,lwfjs,iscale,ntop)
      if (ier1.eq.8) then
         ier = 8
	 return
      endif
c
      do l=0,nterms
         do m=-l,l
	    zh = fjs(l)
	    zhn = fjder(l)*zk0
	    zmul = zh*zh + zhn*zhn
	    local(l,m) = (zh*local(l,m) + zhn*localn(l,m))/zmul
         enddo
      enddo
      return
      end
c
C
C
C
C
C***********************************************************************
      subroutine h3dlocevalspherestab(local,zk,scale,zshift,radius,
     1           nterms,nterms2,lmp,ynm,ynmd,phitemp,phitempn,
     2           nquad,xnodes,iscale,fjs,fjder,lwfjs,ier)
C***********************************************************************
C
C     This subroutine evaluates a local expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C     phitemp(i,j)  : jth mode of phi at ith quad node.
C     phitempn(i,j) : jth mode of phi at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer iscale(0:lwfjs)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 local(0:lmp,-lmp:lmp)
ccc      complex *16 phitemp(nquad,-nterms:nterms)
ccc      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 phitemp(nquad,-nterms2:nterms2)
      complex *16 phitempn(nquad,-nterms2:nterms2)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta
      complex *16 ephi1,ephik,ephi,fjs(0:lwfjs),fjder(0:lwfjs),ztmp1
      complex *16 ut1,ut2,ut3
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
ccc      call prinf(' nterms is *',nterms,1)
ccc      call prinf(' nterms2 is *',nterms2,1)
      do jj=1,nquad
ccc      do m=-nterms,nterms
      do m=-nterms2,nterms2
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
      enddo
      enddo
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
	 call ylgndr2s(nterms,cthetaj,ynm,ynmd)
	 call jfuns3d(jer,nterms,z,scale,fjs,ifder,fjder,
     1        lwfjs,iscale,ntop)
         if (jer.eq.8) then
            ier = 8
	    return
         endif
	 do n = 0,nterms
	    fjder(n) = fjder(n)*zk
         enddo
	 do n = 1,nterms
	    do m = 1,n
	       ynm(n,m) = ynm(n,m)*sthetaj
            enddo
         enddo
	 phitemp(jj,0) = local(0,0)*fjs(0)
	 phitempn(jj,0) = local(0,0)*fjder(0)*rn
ccc         do n=1,nterms2
c
c  4/4/2010 - restoring n=1,nterms 
c
         do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                local(n,0)*fjs(n)*ynm(n,0)
	    ut1 = fjder(n)*rn
	    ut2 = fjs(n)*thetan
	    ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ut3*local(n,0)
ccc            ur = fjder(n)*ynm(n,0)*local(n,0)
ccc            utheta = -local(n,0)*fjs(n)*ynmd(n,0)*sthetaj
ccc	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
ccc	    do m=1,n
c
c  4/4/2010 - only modes up to nterms2 are needed on smaller sphere.
c
	    do m=1,min(n,nterms2)
c
	       ztmp1 = fjs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                local(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                local(n,-m)*ztmp1
	       ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ut3*local(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ut3*local(n,-m)
ccc               ur = fjder(n)*ynm(n,m)*local(n,m)
ccc               utheta = -local(n,m)*fjs(n)*ynmd(n,m)
ccc	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
ccc               ur = fjder(n)*ynm(n,m)*local(n,-m)
ccc               utheta = -local(n,-m)*fjs(n)*ynmd(n,m)
ccc	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
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
C***********************************************************************
      subroutine h3dlocevalspherestab_fast(local,zk,scale,zshift,radius,
     1           nterms,nterms2,lmp,ynm,ynmd,phitemp,phitempn,
     2           nquad,xnodes,iscale,fjs,fjder,rat1,rat2,lwfjs,ier)
C***********************************************************************
C
C     This subroutine evaluates a local expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C     phitemp(i,j)  : jth mode of phi at ith quad node.
C     phitempn(i,j) : jth mode of phi at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer iscale(0:lwfjs)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 local(0:lmp,-lmp:lmp)
ccc      complex *16 phitemp(nquad,-nterms:nterms)
ccc      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 phitemp(nquad,-nterms2:nterms2)
      complex *16 phitempn(nquad,-nterms2:nterms2)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta
      complex *16 ephi1,ephik,ephi,fjs(0:lwfjs),fjder(0:lwfjs),ztmp1
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
ccc      call prinf(' nterms is *',nterms,1)
ccc      call prinf(' nterms2 is *',nterms2,1)
      do jj=1,nquad
ccc      do m=-nterms,nterms
      do m=-nterms2,nterms2
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
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
	 call jfuns3d(jer,nterms,z,scale,fjs,ifder,fjder,
     1        lwfjs,iscale,ntop)
         if (jer.eq.8) then
            ier = 8
	    return
         endif
	 do n = 0,nterms
	    fjder(n) = fjder(n)*zk
         enddo
	 do n = 1,nterms
	    do m = 1,n
	       ynm(n,m) = ynm(n,m)*sthetaj
            enddo
         enddo
	 phitemp(jj,0) = local(0,0)*fjs(0)
	 phitempn(jj,0) = local(0,0)*fjder(0)*rn
ccc         do n=1,nterms2
c
c  4/4/2010 - restoring n=1,nterms 
c
         do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                local(n,0)*fjs(n)*ynm(n,0)
	    ut1 = fjder(n)*rn
	    ut2 = fjs(n)*thetan
	    ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ut3*local(n,0)
ccc            ur = fjder(n)*ynm(n,0)*local(n,0)
ccc            utheta = -local(n,0)*fjs(n)*ynmd(n,0)*sthetaj
ccc	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
c
ccc	    do m=1,n
c
c  4/4/2010 - only modes up to nterms2 are needed on smaller sphere.
c
	    do m=1,min(n,nterms2)
c
	       ztmp1 = fjs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                local(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                local(n,-m)*ztmp1
	       ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ut3*local(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ut3*local(n,-m)
ccc               ur = fjder(n)*ynm(n,m)*local(n,m)
ccc               utheta = -local(n,m)*fjs(n)*ynmd(n,m)
ccc	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
ccc               ur = fjder(n)*ynm(n,m)*local(n,-m)
ccc               utheta = -local(n,-m)*fjs(n)*ynmd(n,m)
ccc	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
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
C***********************************************************************
      subroutine h3dlocevalspherestabold(local,zk,scale,phival,phivaln,
     1           zshift,radius,nterms,nterms2,lmp,ynm,ynmd,
ccc     2           phitemp,phitempn,nquad,xnodes)
     2         phitemp,phitempn,nquad,xnodes,iscale,fjs,fjder,lwfjs,ier)
C***********************************************************************
C
C     This subroutine evaluates a local expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C     phitemp  : intermediate storage in O(p^3) algorithm.
C     phitempn : intermediate storage in O(p^3) algorithm.
C     nquad    : number of quadrature nodes
C                              on target sphere is nquad*nquad.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C     phivanl  : value of normal deriv of potential on tensor product
C                              mesh on target sphere.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer iscale(0:lwfjs)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquad,nquad)
      complex *16 phivaln(nquad,nquad)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta
      complex *16 ephi1,ephik,ephi,fjs(0:lwfjs),fjder(0:lwfjs),ztmp1
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
ccc      call prinf(' nterms is *',nterms,1)
ccc      call prinf(' nterms2 is *',nterms2,1)
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
      enddo
      enddo
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
	 call ylgndr2s(nterms,cthetaj,ynm,ynmd)
	 call jfuns3d(jer,nterms,z,scale,fjs,ifder,fjder,
     1        lwfjs,iscale,ntop)
         if (jer.eq.8) then
            ier = 8
	    return
         endif
	 do n = 0,nterms
	    fjder(n) = fjder(n)*zk
         enddo
	 do n = 1,nterms
	    do m = 1,n
	       ynm(n,m) = ynm(n,m)*sthetaj
            enddo
         enddo
	 phitemp(jj,0) = local(0,0)*fjs(0)
	 phitempn(jj,0) = local(0,0)*fjder(0)*rn
ccc         do n=1,nterms2
c
c   not tested, but same fix as for h3dlocevalspherestab_fast
c
         do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                local(n,0)*fjs(n)*ynm(n,0)
            ur = fjder(n)*ynm(n,0)*local(n,0)
            utheta = -local(n,0)*fjs(n)*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
ccc	    do m=1,n
	    do m=1,min(n,nterms2)
	       ztmp1 = fjs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                local(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                local(n,-m)*ztmp1
               ur = fjder(n)*ynm(n,m)*local(n,m)
               utheta = -local(n,m)*fjs(n)*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
               ur = fjder(n)*ynm(n,m)*local(n,-m)
               utheta = -local(n,-m)*fjs(n)*ynmd(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
	    enddo
	 enddo
      enddo
c
      do jj = 1,nquad
         ephi1 = cdexp(2*pi*imag/nquad)
         ephik = 1.0d0
         do kk = 1,nquad
         phival(jj,kk) = phitemp(jj,0)
         phivaln(jj,kk) = phitempn(jj,0)
         ephik = ephik*ephi1
         ephi = ephik
	 do m = 1,nterms2
ccc	 do m = 1,nterms
cc	    phi = 2*pi*kk*m/nquad
cc	    ephi = cdexp(imag*phi)
            phival(jj,kk)=phival(jj,kk)+phitemp(jj,m)*ephi
            phival(jj,kk)=phival(jj,kk)+phitemp(jj,-m)*dconjg(ephi)
            phivaln(jj,kk)=phivaln(jj,kk)+phitempn(jj,m)*ephi
            phivaln(jj,kk)=phivaln(jj,kk)+phitempn(jj,-m)*dconjg(ephi)
	    ephi = ephi*ephik
         enddo
      enddo
      enddo
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine h3dprojlocnmsep(nterms,ldl,nquadn,ntold,xnodes,wts,
     1           phitemp,local,ynm)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl = dimension param for local expansion
C           nquadn = number of quadrature nodes in polar angle.
C           ntold =  number of azimuthal (m) modes in phitemp
C           xnodes = quad nodes in theta (nquadn of them)
C           wts = quad weights in theta (nquadn of them)
C           phitemp = tabulated function
C                    phitemp(i,j) = jth mode in phi at ith quad node 
C                                   in theta
C           ynm     = workspace for assoc Legendre functions
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C    NOTE:
C
C    yrecursion.f produces Ynm with a nonstandard scaling:
C    (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C    is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C   This accounts for factor (1/2) = (2*pi) * 1/(4*pi) in zmul below.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phitemp(nquadn,-ntold:ntold)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
         enddo
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-ntold,ntold
	    zmul = phitemp(jj,m)*wts(jj)/2
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dprojlocnmsep_fast
     $     (nterms,ldl,nquadn,ntold,xnodes,wts,
     1           phitemp,local,ynm,rat1,rat2)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl = dimension param for local expansion
C           nquadn = number of quadrature nodes in polar angle.
C           ntold =  number of azimuthal (m) modes in phitemp
C           xnodes = quad nodes in theta (nquadn of them)
C           wts = quad weights in theta (nquadn of them)
C           phitemp = tabulated function
C                    phitemp(i,j) = jth mode in phi at ith quad node 
C                                   in theta
C           ynm     = workspace for assoc Legendre functions
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C    NOTE:
C
C    yrecursion.f produces Ynm with a nonstandard scaling:
C    (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C    is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C   This accounts for factor (1/2) = (2*pi) * 1/(4*pi) in zmul below.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phitemp(nquadn,-ntold:ntold)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      real *8 rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
         enddo
      enddo
c
c     get local exp
c
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndrf(nterms,cthetaj,ynm,rat1,rat2)
         do m=-ntold,ntold
	    zmul = phitemp(jj,m)*wts(jj)/2
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
C***********************************************************************
      subroutine h3dprojlocnmold(nterms,ldl,nquadn,nquadm,xnodes,wts,
     1           phival,local,ynm,w,lw,ier)
C***********************************************************************
C
C     This subroutine projects a function given on the unit sphere
C     onto spherical harmonics.
C     It is intended for use with a library of special 
C     functions in which Ynm does NOT include the 
C     1/sqrt(4*pi) factor. Thus, the Ynm are orthogonal
C     but not orthonormal -> the inner product value is 4*pi.
C     (There is no universal definition - physics, spectral analysis,
C     and magnetics communities all use different normalizations.)
c
C     phival(theta,phi) is the function tabulated at nquadn*nquadm pts:
C     Legendre nodes x_j = cos(theta_j) in polar angle and 
C     equispaced in the azimuthal direction.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl =    defines leading dimension of local
C           nquadn = number of quadrature nodes in polar angle.
C           nquadm = number of quadrature nodes in azimuthal direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C           ynm    = work array to hold spherical harmonics
C           w       = workspace of length lw
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C           ier = error flag 
C                 0 normal return
C                 6 insufficient memory in w.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier,ldl
      integer l,m,jj,kk
      real *8 w(lw)
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phival(nquadn,nquadm)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      imarray = 1
      lmarray = 2*nquadn*(2*nterms+1)+7
C
      lused = imarray+lmarray
C
      if (lused.gt.lw) then
         ier=6
	 return
      endif
c
      call h3dprojloc0nm(nterms,ldl,nquadn,nquadm,xnodes,wts,phival,
     1           local,w(imarray),ynm)
      return
      end
c
C***********************************************************************
      subroutine h3dprojloc0nm(nterms,ldl,nquadn,nquadm,xnodes,wts,
     1           phival,local,marray,ynm)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl = dimension param for local expansion
C           nquadn = number of quadrature nodes in polar angle.
C           nquadm = number of quadrature nodes in azimuthal direction.
C           xnodes = quad nodes in theta (nquadn of them)
C           wts = quad weights in theta (nquadn of them)
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           marray  = workspace for intermediate values
C           ynm     = workspace for assoc Legendre functions
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
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
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phival(nquadn,nquadm)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 marray(nquadn,-nquadm/2:nquadm/2)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
         enddo
      enddo
c
c     create marray (intermediate array)
c
      emul1 = cdexp(imag*2*pi/nquadm)
ccc      emul = cdexp(-imag*2*nterms*pi/nquadm)
      emul = cdexp(-imag*2*(nquadm/2)*pi/nquadm)
      do m=-nquadm/2,nquadm/2
cc	    emul = cdexp(imag*m*2*pi/nquadm)
	 do jj=1,nquadn
	    sum = 0
	    ephi = emul
	    do kk = 1,nquadm
               sum = sum + phival(jj,kk)*dconjg(ephi)
	       ephi = ephi*emul
            enddo
ccc	    marray(jj,m) = sum*2*pi/nquad
	    marray(jj,m) = sum/(2*nquadm)
         enddo
	 emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-nquadm/2,nquadm/2
	    zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
C
C
C
C***********************************************************************
      subroutine h3dmpevalspherenm(mpole,zk,scale,zshift,radius,
     1           nterms,lmp,ynm,phitemp,nquad,xnodes,fhs,fhder)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension param for mpole array
C     ynm      : storage for assoc Legendre functions
C     phitemp  : storage for temporary array in O(p^3) scheme 
C     nquad    : number of quadrature nodes in theta
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 0
      ifder = 0
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
      enddo
      enddo
      do jj=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
         rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
         rj = dsqrt(rj)
	 cthetaj = (zshift+radius*ctheta)/rj
	 z = zk*rj
	 call ylgndr(nterms,cthetaj,ynm)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do m=-nterms,nterms
	    mabs = abs(m)
	    do n=mabs,nterms
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*fhs(n)*ynm(n,mabs)
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
C***********************************************************************
      subroutine h3dmpevalspherenm_fast(mpole,zk,scale,zshift,radius,
     1           nterms,lmp,ynm,phitemp,nquad,xnodes,fhs,fhder,
     2           rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension param for mpole array
C     ynm      : storage for assoc Legendre functions
C     phitemp  : storage for temporary array in O(p^3) scheme 
C     nquad    : number of quadrature nodes in theta
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
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
      iffld = 0
      ifder = 0
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
      enddo
      enddo
      call ylgndrini(nterms,rat1,rat2)
      do jj=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
         rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
         rj = dsqrt(rj)
	 cthetaj = (zshift+radius*ctheta)/rj
	 z = zk*rj
	 call ylgndrf(nterms,cthetaj,ynm,rat1,rat2)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do m=-nterms,nterms
	    mabs = abs(m)
	    do n=mabs,nterms
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*fhs(n)*ynm(n,mabs)
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
C***********************************************************************
      subroutine h3dmpevalspherenmold(mpole,zk,scale,phival,zshift,
     1     radius,nterms,lmp,ynm,phitemp,nquad,nquadm,xnodes,fhs,fhder)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
C     mpole    : coefficients of original multipole exp.
C     zk       : Helmholtz parameter
C     scale    : mpole scaling parameter
C     zshift   : shift distance along z-axis.
C     radius   : radius of sphere about (0,0,zshift)
C                              where phival is computed.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension param for mpole array
C     ynm      : storage for assoc Legendre functions
C     phitemp  : storage for temporary array in O(p^3) scheme 
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C                   total on target sphere is nquad*nquadm.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phival(nquad,nquadm)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z
      complex *16 ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
      data imag/(0.0d0,1.0d0)/
C
C----- shift along z-axis.
C      note that everything is scaled.
C
      pi = 4.0d0*datan(1.0d0)
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      iffld = 0
      ifder = 0
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(jj,m) = 0.0d0
      enddo
      enddo
      do jj=1,nquad
	 ctheta = xnodes(jj)
	 stheta = dsqrt(1.0d0 - ctheta**2)
         rj = (zshift+ radius*ctheta)**2 + (radius*stheta)**2
         rj = dsqrt(rj)
	 cthetaj = (zshift+radius*ctheta)/rj
	 z = zk*rj
	 call ylgndr(nterms,cthetaj,ynm)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do m=-nterms,nterms
	    mabs = abs(m)
	    do n=mabs,nterms
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*fhs(n)*ynm(n,mabs)
	    enddo
	 enddo
      enddo
c
      do jj = 1,nquad
      ephi1 = cdexp(2*pi*imag/nquadm)
      ephik = 1.0d0
      do kk = 1,nquadm
         phival(jj,kk) = 0.0d0
         ephik = ephik*ephi1
         ephi = ephik**(-nterms)
	 do m = -nterms,nterms
ccc	    phi = 2*pi*kk*m/nquadm
ccc	    ephi = cdexp(imag*phi)
            phival(jj,kk) = phival(jj,kk) + phitemp(jj,m)*ephi
	    ephi = ephi*ephik
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
C***********************************************************************
      subroutine h3dmpevalspherenmstab(mpole,zk,scale,zshift,
     1           radius,nterms,lmp,ynm,ynmd,phitemp,phitempn,
     2           nquad,xnodes,fhs,fhder)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
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
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta,ut1,ut2
      complex *16 ut3,ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
      complex *16 ztmp1
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
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
      enddo
      enddo
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
	 call ylgndr2s(nterms,cthetaj,ynm,ynmd)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do i = 0,nterms
	    fhder(i) = fhder(i)*zk
         enddo
         do n=1,nterms
	    do m = 1,n
  	       ynm(n,m) = ynm(n,m)*sthetaj
              enddo
         enddo
	 phitemp(jj,0) = mpole(0,0)*fhs(0)
	 phitempn(jj,0) = mpole(0,0)*fhder(0)*rn
	 do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                mpole(n,0)*fhs(n)*ynm(n,0)
	    ut1 = fhder(n)*rn
	    ut2 = fhs(n)*thetan
	    ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ut3*mpole(n,0)
ccc	    ur = fhder(n)*ynm(n,0)*mpole(n,0)
ccc	    utheta = -mpole(n,0)*fhs(n)*ynmd(n,0)*sthetaj
ccc	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
            do m=1,n
	       ztmp1 = fhs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                mpole(n,-m)*ztmp1
	       ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ut3*mpole(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ut3*mpole(n,-m)
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,m)
cc	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,-m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,-m)
cc	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
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
C***********************************************************************
      subroutine h3dmpevalspherenmstab_fast(mpole,zk,scale,zshift,
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
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
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
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
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
	 phitemp(jj,0) = mpole(0,0)*fhs(0)
	 phitempn(jj,0) = mpole(0,0)*fhder(0)*rn
	 do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                mpole(n,0)*fhs(n)*ynm(n,0)
	    ut1 = fhder(n)*rn
	    ut2 = fhs(n)*thetan
	    ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ut3*mpole(n,0)
ccc	    ur = fhder(n)*ynm(n,0)*mpole(n,0)
ccc	    utheta = -mpole(n,0)*fhs(n)*ynmd(n,0)*sthetaj
ccc	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
            do m=1,n
	       ztmp1 = fhs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                mpole(n,-m)*ztmp1
	       ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ut3*mpole(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ut3*mpole(n,-m)
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,m)
cc	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,-m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,-m)
cc	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
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
C***********************************************************************
      subroutine h3dmpevalspherenmstabold(mpole,zk,scale,phival,phivaln,
     1           zshift,radius,nterms,lmp,ynm,ynmd,phitemp,phitempn,
cc     2           nquad,nquadm,xnodes)
     2           nquad,nquadm,xnodes,fhs,fhder)
C***********************************************************************
C
C     This subroutine evaluates a multipole expansion on a target
C     sphere at a distance (0,0,zshift) from the origin of radius 
C     "radius".
C
C---------------------------------------------------------------------
C     INPUT:
C
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
C     phitemp  : storage for temporary array in O(p^3) scheme 
C     phitempn : storage for temporary array in O(p^3) scheme 
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C                       total on target sphere is nquad*nquadm.
C     xnodes   : Legendre nodes x_j = cos theta_j.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : value of potential on tensor product
C                              mesh on target sphere.
C     phivaln  : value of normal deriv of potential on tensor product
C                              mesh on target sphere.
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms
      integer l,m,jnew,knew
      real *8 zshift, targ(3), center(3)
      real *8 xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 phival(nquad,nquadm)
      complex *16 phivaln(nquad,nquadm)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 phitempn(nquad,-nterms:nterms)
      complex *16 imag,pot,fld(3), zk,z,uval,unval,ur,utheta,ut1,ut2
      complex *16 ut3,ephi1,ephik,ephi,fhs(0:nterms),fhder(0:nterms)
      complex *16 ztmp1
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
         phitemp(jj,m) = 0.0d0
         phitempn(jj,m) = 0.0d0
      enddo
      enddo
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
	 call ylgndr2s(nterms,cthetaj,ynm,ynmd)
	 call h3dall(nterms,z,scale,fhs,ifder,fhder)
         do i = 0,nterms
	    fhder(i) = fhder(i)*zk
         enddo
         do n=1,nterms
	    do m = 1,n
  	       ynm(n,m) = ynm(n,m)*sthetaj
              enddo
         enddo
	 phitemp(jj,0) = mpole(0,0)*fhs(0)
	 phitempn(jj,0) = mpole(0,0)*fhder(0)*rn
	 do n=1,nterms
	    phitemp(jj,0) = phitemp(jj,0) +
     1                mpole(n,0)*fhs(n)*ynm(n,0)
	    ut1 = fhder(n)*rn
	    ut2 = fhs(n)*thetan
	    ut3 = ut1*ynm(n,0)-ut2*ynmd(n,0)*sthetaj
	    phitempn(jj,0) = phitempn(jj,0)+ut3*mpole(n,0)
ccc	    ur = fhder(n)*ynm(n,0)*mpole(n,0)
ccc	    utheta = -mpole(n,0)*fhs(n)*ynmd(n,0)*sthetaj
ccc	    phitempn(jj,0) = phitempn(jj,0)+ur*rn+utheta*thetan
            do m=1,n
	       ztmp1 = fhs(n)*ynm(n,m)
	       phitemp(jj,m) = phitemp(jj,m) +
     1                mpole(n,m)*ztmp1
	       phitemp(jj,-m) = phitemp(jj,-m) +
     1                mpole(n,-m)*ztmp1
	       ut3 = ut1*ynm(n,m)-ut2*ynmd(n,m)
	       phitempn(jj,m) = phitempn(jj,m)+ut3*mpole(n,m)
	       phitempn(jj,-m) = phitempn(jj,-m)+ut3*mpole(n,-m)
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,m)
cc	       phitempn(jj,m) = phitempn(jj,m)+ur*rn+utheta*thetan
cc	       ur = fhder(n)*ynm(n,m)*mpole(n,-m)
cc	       utheta = -fhs(n)*ynmd(n,m)*mpole(n,-m)
cc	       phitempn(jj,-m) = phitempn(jj,-m)+ur*rn+utheta*thetan
	    enddo
	 enddo
      enddo
c
      do jj = 1,nquad
         ephi1 = cdexp(2*pi*imag/nquadm)
         ephik = 1.0d0
         do kk = 1,nquadm
            phival(jj,kk) = phitemp(jj,0)
            phivaln(jj,kk) = phitempn(jj,0)
            ephik = ephik*ephi1
            ephi = ephik
	    do m = 1,nterms
cc	       phi = 2*pi*kk*m/nquadm
cc	       ephi = exp(imag*phi)
               phival(jj,kk) = phival(jj,kk) + phitemp(jj,m)*ephi +
     1      	                    phitemp(jj,-m)*dconjg(ephi)
               phivaln(jj,kk) = phivaln(jj,kk) + phitempn(jj,m)*ephi +
     1                              phitempn(jj,-m)*dconjg(ephi) 
	       ephi = ephi*ephik
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
C***********************************************************************
      subroutine h3dprojlocsepstab(nterms,ldl,nquadn,ntold,xnodes,wts,
     1           phitemp,phitempn,local,local2,ynm)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 local2(0:ldl,-ldl:ldl)
      complex *16 phitemp(nquadn,-ntold:ntold)
      complex *16 phitempn(nquadn,-ntold:ntold)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
	    local2(l,m) = 0.0d0
         enddo
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-ntold,ntold
	    zmul = phitemp(jj,m)*wts(jj)/2.0d0
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-ntold,ntold
	    zmul = phitempn(jj,m)*wts(jj)/2.0d0
            do l=abs(m),nterms
               local2(l,m) = local2(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
c
c
C
C
C***********************************************************************
      subroutine h3dprojlocsepstab_fast
     $          (nterms,ldl,nquadn,ntold,xnodes,wts,
     1           phitemp,phitempn,local,local2,ynm,rat1,rat2)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
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
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 local2(0:ldl,-ldl:ldl)
      complex *16 phitemp(nquadn,-ntold:ntold)
      complex *16 phitempn(nquadn,-ntold:ntold)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      real *8 rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
	    local2(l,m) = 0.0d0
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
	    zmul = phitemp(jj,m)*wts(jj)/2.0d0
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
	    zmul = phitempn(jj,m)*wts(jj)/2.0d0
            do l=abs(m),nterms
               local2(l,m) = local2(l,m) + 
     1   	       zmul*ynm(l,abs(m))
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
C***********************************************************************
      subroutine h3dprojloc0nmstab(nterms,ldl,nquadn,nquadm,xnodes,wts,
     1           phival,phivaln,local,local2,marray,ynm)
C***********************************************************************
C
C     compute spherical harmonic expansion on unit sphere
C     of function tabulated at nquadn*nquadm grid points.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl = order of spherical harmonic expansion
C           nquadn = number of quadrature nodes in polar angle (theta).
C           nquadm = number of quadrature nodes in azimuthal direction.
C           xnodes = quad nodes in theta (polar angle) - nquadn of them
C           wts  = quad weights in theta (polar angle)
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C           phivaln = tabulated normal deriv of function
C                    phivaln(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
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
Cccc	    marray(jj,m) = sum*2*pi/nquad
C	    marray(jj,m) = sum/(2*nquad)
C
C   The latter has incorporated the 1/(4*pi) normalization factor
C   into the azimuthal quadrature weight (2*pi/nquad).
C
C---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,nquadn,nquadm,ier
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8 ynm(0:nterms,0:nterms)
      complex *16 zk,phival(nquadn,nquadm)
      complex *16 phivaln(nquadn,nquadm)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 local2(0:ldl,-ldl:ldl)
      complex *16 marray(nquadn,-nquadm/2:nquadm/2)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
	    local(l,m) = 0.0d0
	    local2(l,m) = 0.0d0
         enddo
      enddo
c
c     create marray for potential (intermediate array)
c
      emul1 = cdexp(imag*2*pi/nquadm)
      emul = cdexp(-imag*2*(nquadm/2)*pi/nquadm)
      do m=-nquadm/2,nquadm/2
	 do jj=1,nquadn
	    sum = 0
	    ephi = emul
	    do kk = 1,nquadm
               sum = sum + phival(jj,kk)*dconjg(ephi)
	       ephi = ephi*emul
            enddo
	    marray(jj,m) = sum/(2*nquadm)
         enddo
	 emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-nquadm/2,nquadm/2
	    zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
c
c     create marray for normal derivative (intermediate array)
c
      emul1 = cdexp(imag*2*pi/nquadm)
      emul = cdexp(-imag*2*(nquadm/2)*pi/nquadm)
      do m=-nquadm/2,nquadm/2
	 do jj=1,nquadn
	    sum = 0
	    ephi = emul
	    do kk = 1,nquadm
               sum = sum + phivaln(jj,kk)*dconjg(ephi)
	       ephi = ephi*emul
            enddo
	    marray(jj,m) = sum/(2*nquadm)
         enddo
	 emul = emul*emul1
      enddo
c
c     get local exp
c
      do jj=1,nquadn
	 cthetaj = xnodes(jj)
	 call ylgndr(nterms,cthetaj,ynm)
         do m=-nquadm/2,nquadm/2
	    zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local2(l,m) = local2(l,m) + 
     1   	       zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end
c
c
c
c
