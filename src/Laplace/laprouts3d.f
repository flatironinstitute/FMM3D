cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
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
c      This file contains the basic subroutines for 
c      forming and evaluating multipole expansions.
c
c      Remarks on scaling conventions.
c
c      1)  Far field and local expansions are consistently rscaled as
c              
c
c          M_n^m (scaled) = M_n^m / rscale^(n)  so that upon evaluation
c
c          the field is  sum   M_n^m (scaled) * rscale^(n) / r^{n+1}.
c
c          L_n^m (scaled) = L_n^m * rscale^(n)  so that upon evaluation
c
c          the field is  sum   L_n^m (scaled) / rscale^(n) * r^{n}.
c
c
c      2) There are many definitions of the spherical harmonics,
c         which differ in terms of normalization constants. We
c         adopt the following convention:
c
c         For m>0, we define Y_n^m according to 
c
c         Y_n^m = \sqrt{2n+1} \sqrt{\frac{ (n-m)!}{(n+m)!}} \cdot
c                 P_n^m(\cos \theta)  e^{i m phi} 
c         and
c 
c         Y_n^-m = dconjg( Y_n^m )
c    
c         We omit the Condon-Shortley phase factor (-1)^m in the 
c         definition of Y_n^m for m<0. (This is standard in several
c         communities.)
c
c         We also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the Y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  (This is also standard in several communities.)
c         More precisely, 
c
c                 \int_S Y_n^m Y_n^m d\Omega = 4 \pi. 
c
c         Using our standard definition, the addition theorem takes 
c         the simple form 
c
c         1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1)) 
c
c         1/r = 
c         \sum_n \sum_m  |S|^n  Ylm*(S)    Ylm(T)     / (|T|^(n+1)) 
c                               -------    ------
c                               sqrt(2n+1) sqrt(2n+1)
c
c        In the Laplace library (this library), we incorporate the
c        sqrt(2n+1) factor in both forming and evaluating multipole
c        expansions.
c
c-----------------------------------------------------------------------
c
cc      l3dmpeval: computes potential and -grad(potential)
c                 due to a multipole expansion.
c
cc      l3dformmp: creates multipole expansion (outgoing) due to 
c                 a collection of charges.
c
cc      l3dtaeval: computes potential and -grad(potential) 
c                  due to local expansion.
c
cc      l3dformta: creates local expansion due to 
c                 a collection of charges.
c
cc      lpotfld3d: direct calculation of pot/field at target
c                 due to a single charge 
c
cc      lpotfld3dall:  direct calculation of pot/field at target 
c                 due to set of charges
c
cc      lpotfld3dall_targ: direct calculation of pot/field at target
c                 due to set of charges (UNROLLED VERSION)
c
cc      l3dadd: adds one spherical harmonic expansion to another
c
c
cc      l3dadd_trunc: adds one spherical harmonic expansion to another
c                 allowing for different spherical harmonic array dimensions  
c
cc      cart2polarl: utility function.
c                 converts Cartesian coordinates into polar
c                 representation needed by other routines.
c
cc      l3drhpolar: utility function
c                 converts Cartesian coordinates into 
c                 r, cos(theta), e^{i*phi).
c
cc      l3dformmp_dp: creates multipole expansion (outgoing) due to 
c                 a collection of dipoles.
c
cc      lpotfld3d_dp: direct calculation of pot/field at target
c                 due to single dipole.
c
cc      lpotfld3dall_dp: direct calculation of pot/field at target
c                 due to set of dipoles.
c
cc      lpotfld3dall_dp_targ: direct calculation of pot/field at target
c                 due to set of dipoles (UNROLLED VERSION).
c
cc      lpotfld3dall_sdp_targ: direct calculation of pot/field at target
c                 due to set of charges and dipoles (UNROLLED VERSION).
c
cc      l3dformta_dp: creates local expansion due to 
c                 a collection of dipoles.
c
cc      l3dmpevalall_trunc: computes potential and -grad(potential)
c                 due to a multipole expansion (OPTIMIZED VERSION)
c                 at a collection of targets.
c
cc      l3dmpeval_trunc: computes potential and -grad(potential)
c                 due to a multipole expansion (OPTIMIZED VERSION)
c                 at a single target.
c
cc      l3dtaevalall_trunc: computes potential and -grad(potential)
c                 due to a local expansion (OPTIMIZED VERSION)
c                 at a collection of targets.
c
cc      l3dtaeval_trunc: computes potential and -grad(potential)
c                 due to a local expansion (OPTIMIZED VERSION)
c                 at a single target.
c
cc      l3dformmp_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of charges (OPTIMIZED VERSION).
c
cc      l3dformmp_add_trunc: *increments* multipole expansion (outgoing) 
c                 due to a collection of charges (OPTIMIZED VERSION).
c
cc      l3dformta_trunc: creates local expansion (outgoing) due to 
c                 a collection of charges (OPTIMIZED VERSION).
c
cc      l3dformta_add_trunc: *increments* local expansion (outgoing) 
c                 due to a collection of charges (OPTIMIZED VERSION).
c
cc      l3dformmp_dp_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of dipoles (OPTIMIZED VERSION).
c
cc      l3dformmp_dp_add_trunc: *increments* multipole expansion (outgoing) 
c                 due to a collection of dipoles (OPTIMIZED VERSION).
c
cc      l3dformta_dp_trunc: creates local expansion (outgoing) due to 
c                 a collection of dipoles (OPTIMIZED VERSION).
c
cc      l3dformta_dp_add_trunc: *increments* local expansion (outgoing) 
c                 due to a collection of dipoles (OPTIMIZED VERSION).
c
cc      l3dformmp_charge_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of *real-valued* charges (OPTIMIZED VERSION).
c
cc      l3dformmp_dipole_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of *real-valued* dipoles (OPTIMIZED VERSION).
c
cc      getsgnformpmp_dipole: utility function - not optimal, but 
c                 working array needed by some formmp routines.
c
cc      l3dmpevalhess: computes potential, -grad(potential), Hessian
c                 due to a outgoing multipole expansion at a single target.
c
cc      l3dtaevalhess: computes potential, -grad(potential), Hessian
c                 due to a local multipole expansion at a single target.
c
cc      lpotfld3dallhess: direct calculation of pot/field/Hessian at target
c                 due to set of charges.
c
cc      lpotfld3dhess: direct calculation of pot/field/Hessian at target
c                 due to single charge.
c
cc      lpotfld3dallhess_dp: direct calculation of pot/field/Hessian at target
c                 due to set of dipoles.
c
cc      lpotfld3dhess_dp: direct calculation of pot/field/Hessian at target
c                 due to single dipole.
c
cc      l3dformmp_quad: creates multipole expansion (outgoing) due to 
c                 a collection of *real-valued* quadrupoles.
c
cc      l3dformmp_quad_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of *real-valued* quadrupoles (OPTIMIZED VERSION).
c
cc      lpotfld3dall_quad: direct calculation of pot/field at target
c                 due to set of *real-valued* quadrupoles.
c
cc      lpotfld3d_quad: direct calculation of pot/field at target
c                 due to single *real-valued* quadrupole.
c
cc      l3dformta_quad: creates local expansion (incoming) due to 
c                 a single *real-valued* quadrupole (preliminary version)
c
cc      getsgnformpmp_quad: utility function - not optimal, but 
c                 working array needed by some formmp routines.
c
cc      l3dmpevalhessd: computes potential, -grad(potential), Hessian
c                 due to a multipole expansion at a single target.
c
cc      l3dmpevalhessdini:  initialization routine for l3dmpevalhessd, etc.
c
cc      l3dmpevalhessd_trunc: computes potential, -grad(potential), Hessian
c                 due to a multipole expansion at a single target (OPTIMIZED).
c
cc      l3dtaevalhessd: computes potential, -grad(potential), Hessian
c                 due to a local expansion at a single target.
c
cc      l3dtaevalhessdini:  initialization routine for l3dtaevalhessd, etc.
c
cc      l3dtaevalhessd_trunc: computes potential, -grad(potential), Hessian
c                 due to a local expansion at a single target (OPTIMIZED).
c
c
cc      l3dformmp_qp: creates multipole expansion (outgoing) due to 
c                 a collection of quadrupoles.
c
cc      l3dformmp_qp_trunc: creates multipole expansion (outgoing) due to 
c                 a collection of quadrupoles (OPTIMIZED VERSION).
c
cc      l3dformmp_qp_add_trunc: *increments* multipole expansion 
c                 (outgoing) due to a collection of quadrupoles 
c                 (OPTIMIZED VERSION).
c
cc      l3dformta_qp: creates local expansion (incoming) due to 
c                 a collection of quadrupoles.
c
cc      l3dformta_qp_trunc: creates local expansion (incoming) due to 
c                 a collection of quadrupoles (OPTIMIZED VERSION).
c
cc      l3dformta_qp_add_trunc: *increments* local expansion 
c                 (incoming) due to a collection of quadrupoles 
c                 (OPTIMIZED VERSION).
c
cc      lpotfld3dallhess_qp: direct calculation of pot/field/hessian at target
c                 due to set of quadrupoles.
c
cc      lpotfld3dhess_qp: direct calculation of pot/field/hessian at target
c                 due to single quadrupole.
c
c**********************************************************************
      subroutine l3dmpeval(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi) / r^{n+1} / sqrt(2n+1)
c             n   m
c
c     fld = -gradient(pot) if iffld = 0.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    -gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c-----------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ier
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused
      double precision rscale
      double precision center(3),ztarg(3)
      double precision, allocatable :: w(:)
      double complex pot,fld(3)
      double complex mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+3)+5 
c
      ifr=iephi+lephi
      ifrder=ifr+(nterms+3)
      lused=ifrder+(nterms+3)
      allocate(w(lused))
c
      call l3dmpeval0(rscale,center,mpole,nterms,ztarg,
     1	   pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval0(rscale,center,mpole,nterms,
     1		ztarg,pot,iffld,fld,ynm,ynmd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dmpeval for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,iffld
      integer i,l,n,m
      double precision rscale
      double precision center(3),ztarg(3),zdiff(3)
      double precision ynm(0:nterms,0:nterms)
      double precision ynmd(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision frder(0:nterms+1)
      double complex pot,fld(3),ephi1,ur,utheta,uphi,ux,uy,uz
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1)
c
      double precision phi,cphi,sphi,theta,ctheta,stheta,thetax,thetay
      double precision thetaz
      double precision rs,rx,ry,rz,d,r,phix,phiy,phiz,done
      double complex eye
      double complex ztmp1,ztmp2,ztmp3,ztmpsum,z
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = 1.0d0/r
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      do i=0,nterms+1
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
      if (iffld.eq.1) then
         rx = stheta*cphi
         thetax = ctheta*cphi/r
         phix = -sphi/r
         ry = stheta*sphi
         thetay = ctheta*sphi/r
         phiy = cphi/r
         rz = ctheta
         thetaz = -stheta/r
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions
c     and scale by 1/sqrt(2l+1)
c
      if (iffld.eq.1) then
         call ylgndr2s(nterms,ctheta,ynm,ynmd)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               ynm(l,m) = ynm(l,m)*rs
               ynmd(l,m) = ynmd(l,m)*rs
            enddo
         enddo
      else        
         call ylgndr(nterms,ctheta,ynm)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               ynm(l,m) = ynm(l,m)*rs
            enddo
         enddo
      endif
c
c     initialize computed values.
c
      if (iffld.eq.1) then
         ur = mpole(0,0)*frder(0)
         utheta = 0.0d0
         uphi = 0.0d0
      endif
      pot=mpole(0,0)*fr(0)
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    ur = ur + frder(n)*ynm(n,0)*mpole(n,0)
	    utheta = utheta -mpole(n,0)*fr(n)*ynmd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)*stheta
	       ztmp2 = mpole(n,m)*ephi(m) 
	       ztmp3 = mpole(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*ynm(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fr(n)*ynmd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fr(n)*ynm(n,m)*ztmpsum
            enddo
         enddo
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)
	       ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp(ier,rscale,sources,charge,ns,center,
     1                  nterms,mpole)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		                   ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c    
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, ier1, lused
      double precision center(3),sources(3,ns)
      double precision rscale,rs
      double complex mpole(0:nterms,-nterms:nterms)
      double complex eye,charge(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      ier = 0
      do i = 1, ns
         call l3dformmp1(ier1,rscale,sources(1,i),charge(i),center,
     1        nterms,mpole)
      enddo
      if (ier1.ne.0) ier = ier1
c
c     scale by 1/sqrt(2l+1)
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            mpole(l,m) = mpole(l,m)*rs
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1(ier,rscale,source,charge,center,
     1		nterms,mpole)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a charge located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : complex charge strength
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c                            
c     mpole   : coeffs of the expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifrder,lfrder,ifr,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex charge
c
c     compute work space components:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifrder=iephi+lephi
      lfrder=2*(nterms+3)
c
      ifr=ifrder+lfrder
      lused=ifr+lfrder

      allocate(w(lused))
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0(rscale,source,charge,center,nterms,
     1		mpole,w(ipp),w(ippd),w(iephi),w(ifr),
     2          w(ifrder))
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0(rscale,source,charge,center,
     1		nterms,mpole,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision r,theta,phi,d,ctheta,stheta,cphi,sphi,dtmp
      double complex mpole(0:nterms,-nterms:nterms)
      double complex charge
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex  fr(0:nterms+1),frder(0:nterms+1)
      double complex  ztmp,z
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = r
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i)=i*fr(i-1)
      enddo
c
c     get the associated Legendre functions:
c
      call ylgndr(nterms,ctheta,pp)
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     multiply all fr's by charge strength.
c
      do n = 0,nterms
         fr(n) = fr(n)*charge
      enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r =  
c          \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)  / (|T|)^{n+1}
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      mpole(0,0)= mpole(0,0) + fr(0)
      do n=1,nterms
         dtmp=pp(n,0)
         mpole(n,0)= mpole(n,0) + dtmp*fr(n)
         do m=1,n
            ztmp=pp(n,m)*fr(n)
            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dtaeval(rscale,center,locexp,nterms,
     1		ztarg,pot,iffld,fld,ier)
c**********************************************************************
c
c     This subroutine evaluates a local expansion centered at CENTER
c     at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi) / sqrt(2n+1)
c             n   m
c
c     The reason for including the sqrt(2n+1) scaling has to do with
c     the addition theorem for 1/r. The term 1/(2n+1) is needed
c     and we put half the weight on the local evaluation and half the 
c     weight on the expansion formation...
c
c     The addition theorem for exp(ikr)/r does not require the 
c     1/(2n+1) scaling - it appears in the definitions of the 
c     Bessel and Hankel functions.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the expansion
c     nterms     : order of the expansion
c     ztarg      : target vector
c     iffld      : flag for gradient computation
c		        iffld=0  - gradient is not computed
c		        iffld=1  - gradient is computed
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential at ztarg(3)
c     fld        : -gradient at ztarg (if requested)
c     lused      : amount of work space "w" used
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer iffld,nterms
      integer ier,ipp,lpp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused
      double precision rscale,center(3),ztarg(3)
      double precision, allocatable :: w(:)
      double complex pot,fld(3)
      double complex locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr= nterms+3
c
      ifrder=ifr+lfr
      lfrder=nterms+3
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dtaeval0(rscale,center,locexp,nterms,ztarg,
     1	     pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),
     2       w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval0(rscale,center,locexp,nterms,
     1		ztarg,pot,iffld,fld,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dtaeval for comments.
c     (pp and ppd are storage arrays for Ynm and Ynm')
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,iffld,i,l,m,n
      double precision rscale,center(3),ztarg(3),zdiff(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision fruse,fr(0:nterms+1),frder(0:nterms+1)
      double precision done,r,theta,phi,d,ctheta,stheta,cphi,sphi
      double precision phix,phiy,phiz,rs,rx,ry,rz
      double precision thetax,thetay,thetaz
      double complex pot,fld(3),ephi1,ephi1inv
      double complex locexp(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1)
c
      double complex eye,ur,utheta,uphi
      double complex ztmp,z
      double complex ztmp1,ztmp2,ztmp3,ztmpsum
      double complex ux,uy,uz
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
c     Convert to spherical coordinates
c
      call cart2polarl(zdiff,r,theta,phi)
      d = rscale*r
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute e^{eye*m*phi} array.
c
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     We compensate for this omission by using one lower power in the 
c     r variable - see fruse below.
c     For the n=0 mode, it is not relevant. 
c     
c
      if (iffld.eq.1) then
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
         call ylgndr2s(nterms,ctheta,pp,ppd)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               pp(l,m) = pp(l,m)*rs
               ppd(l,m) = ppd(l,m)*rs
            enddo
         enddo
      else
         call ylgndr(nterms,ctheta,pp)
         do l = 0,nterms
            rs = sqrt(1.0d0/(2*l+1))
            do m=0,l
               pp(l,m) = pp(l,m)*rs
            enddo
         enddo
      endif
c
c
      pot=locexp(0,0)*fr(0)
      if (iffld.eq.1) then
         ur = 0.0d0
         utheta = 0.0d0
         uphi = 0.0d0
      endif
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms
            pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    ur = ur + frder(n)*pp(n,0)*locexp(n,0)
	    fruse = fr(n-1)*rscale
	    utheta = utheta -locexp(n,0)*fruse*ppd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)*stheta
	       ztmp2 = locexp(n,m)*ephi(m) 
	       ztmp3 = locexp(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*pp(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fruse*ppd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fruse*pp(n,m)*ztmpsum
            enddo
         enddo
ccc	 call prin2(' ur is *',ur,2)
ccc	 call prin2(' utheta is *',utheta,2)
ccc	 call prin2(' uphi is *',uphi,2)
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms
	    pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)
	       ztmp2 = locexp(n,m)*ephi(m)+locexp(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta(ier,rscale,sources,charge,ns,center,
     1		           nterms,locexp)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  deprecated but left in calling sequence for
c		  backward compatibility.
c
c     locexp    : coeffs for the local expansion
c
c----------------------------------------------------------------------
      implicit none
      integer ns,l,m,i,ier,nterms
      double precision rscale,sources(3,ns),center(3),rs
      double complex locexp(0:nterms,-nterms:nterms), charge(ns)
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1(ier,rscale,sources(1,i),charge(i),
     1		center,nterms,locexp)
      enddo
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            locexp(l,m) = locexp(l,m)*rs
         enddo
      enddo
c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta1(ier,rscale,source,charge,center,
     &		nterms,locexp)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single charge located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c INPUT:
c
c     rscale    : scaling parameter
c     source    : coordinates of the source
c     charge    : coordinates of the source
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c---------------------------------------------------------------------
c OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double complex locexp(0:nterms,-nterms:nterms), charge
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      call l3dformta0(rscale,source,charge,center,
     &		nterms,locexp,w(ipp),w(iephi),w(ifr))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0(rscale,source,charge,
     &		center,nterms,locexp,pp,ephi,fr)
c**********************************************************************
c
c     See l3dformta/l3dformta1 for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,n,m,i
      double precision rscale,source(3),center(3),zdiff(3)
      double precision pp(0:nterms,0:nterms)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi,d
      double complex fr(0:nterms+1)
      double complex locexp(0:nterms,-nterms:nterms), charge
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex ztmp,z
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
c
c     get the Ynm
c
      call ylgndr(nterms,ctheta,pp)
c
c     compute radial functions and scale them by charge strength.
c
      do n = 0, nterms
         fr(n) = fr(n)*charge
      enddo
c
c     Compute contributions to locexp
c
      locexp(0,0)=locexp(0,0) + fr(0)
      do n=1,nterms
         locexp(n,0)=locexp(n,0) + pp(n,0)*fr(n)
         do m=1,n
            ztmp=pp(n,m)*fr(n)
	    locexp(n,m)=locexp(n,m) + ztmp*ephi(-m)
	    locexp(n,-m)=locexp(n,-m) + ztmp*ephi(m)
         enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotfld3dall(iffld,sources,charge,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of charges at 
c     SOURCE(3,ns). 
c     
c       pot =  sum_{i=1,..,ns}  charge(i)/| target - sources(*,i)|
c	fld =  -grad(pot)
c
c     It calls a subroutine for each source.
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (double precision)        : calculated potential
c     fld   (double precision)        : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer i,ns,iffld
      double precision sources(3,ns),target(3)
      double complex pot,fld(3),potloc,fldloc(3)
      double complex eye
      double complex charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d(iffld,sources(1,i),charge(i),target,
     1        potloc,fldloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_targ(iffld,sources,charge,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of charges at 
c     SOURCE(3,ns). 
c     
c       pot =  sum_{i=1,..,ns}  charge(i)/| target - sources(*,i)|
c	fld =  -grad(pot)
c
c     It uses a single loop over all sources.
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (double precision)        : calculated potential
c     fld   (double precision)        : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      double precision sources(3,ns),target(3)
      double precision d,dd,xdiff,ydiff,zdiff,dinv,dinv2,dinv3
      double complex pot,fld(3),potloc,fldloc(3)
      double complex eye
      double complex charge(ns),cd
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
        xdiff=target(1)-sources(1,i)
        ydiff=target(2)-sources(2,i)
        zdiff=target(3)-sources(3,i)
        dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
        d=sqrt(dd)
c
        dinv=1.0d0/d
        pot=pot+charge(i)*dinv
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
        xdiff=target(1)-sources(1,i)
        ydiff=target(2)-sources(2,i)
        zdiff=target(3)-sources(3,i)
        dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
        d=sqrt(dd)
c
        dinv=1.0d0/d
        pot=pot+charge(i)*dinv
c
ccc        if (iffld.eq.1) then
        dinv2=dinv*dinv
        dinv3=dinv*dinv2
        fld(1)=fld(1)+charge(i)*xdiff*dinv3
        fld(2)=fld(2)+charge(i)*ydiff*dinv3
        fld(3)=fld(3)+charge(i)*zdiff*dinv3
ccc        endif
c
      enddo
      endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3d(iffld,source,charge,target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a charge at 
c     SOURCE. 
c     
c              	pot = charge/ |target-source|
c		fld = -grad(pot)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld     : flag for computing gradient
c	                 	iffld = 0 -> dont compute 
c		                iffld = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     fld       : calculated gradient
c
c---------------------------------------------------------------------
      implicit none
      integer iffld
      double precision source(3),target(3)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,dinv2,dinv3
      double complex pot,fld(3)
      double complex h0,h1,cd,eye,z,ewavek
      double complex charge
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv=1.0d0/d
      pot=charge*dinv
c
      if (iffld.eq.1) then
         dinv2=dinv*dinv
         dinv3=dinv*dinv2
         fld(1)=charge*xdiff*dinv3
         fld(2)=charge*ydiff*dinv3
         fld(3)=charge*zdiff*dinv3
      endif
      return
      end
c
c**********************************************************************
      subroutine l3dadd(mpole,mpole2,nterms)
c**********************************************************************
c
c     add mpole to mpole2
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,i,j
      double complex mpole(0:nterms,-nterms:nterms)
      double complex mpole2(0:nterms,-nterms:nterms)
c
      do i = 0,nterms
         do j = -i,i
	    mpole2(i,j) = mpole2(i,j)+mpole(i,j)
	 enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dadd_trunc(mpole,mpole2,nterms,ldc)
c**********************************************************************
c
c     add mpole to mpole2, assuming size of mpole is smaller than
c     size of mpole2 (nterms < ldc).
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,ldc,i,j
      double complex mpole(0:nterms,-nterms:nterms)
      double complex mpole2(0:ldc,-ldc:ldc)
c
      do i = 0,nterms
         do j = -i,i
	    mpole2(i,j) = mpole2(i,j)+mpole(i,j)
	 enddo
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine cart2polarl(zat,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	r     :  |zat|
c	theta : angle subtended with respect to z-axis
c	phi   : angle of (zat(1),zat(2)) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit none
      double precision zat(3),r,proj,theta,phi
      double complex ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2+zat(3)**2)
      proj = sqrt(zat(1)**2+zat(2)**2)
c
      theta = datan2(proj,zat(3))
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      phi = 0
      else
      phi = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
        subroutine l3dzero(mpole,nterms)
        implicit double precision (a-h,o-z)
c
c       ... set multipole to zero
c
        double complex mpole(0:nterms,-nterms:nterms)
c       
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
c
c
c
c**********************************************************************
        subroutine l3drhpolar(x,y,z,r,ctheta,ephi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       x,y,z   : Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r      : sqrt(x*x+y*y+z*z)
c       ctheta : cos(theta)
c       ephi   : exp(I*phi)  (complex *16_
c
c       where
c
c       theta is angle subtended with respect to z-axis
c       phi   is angle of (x,y) subtended with 
c               respect to x-axis
c-----------------------------------------------------------------------
        implicit none
        double precision x,y,z,r,ctheta,proj
        double complex ephi,ima
        data ima/(0.0d0,1.0d0)/
c
        proj = sqrt(x*x+y*y)
        r = sqrt(x*x+y*y+z*z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        if( abs(proj) .gt. 0 ) then
        ephi = cmplx(x,y)/proj
        else
        ephi = 0.0d0
        endif
c
        return
        end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_dp(ier,rscale,sources,dipstr,dipvec,ns,
     1                  center,nterms,mpole)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS 
c     dipole sources C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		         ier=0  returned successfully
c		         deprecated but left in calling sequence for
c		         backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c                  
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused
      double precision center(3),sources(3,ns)
      double precision dipvec(3,ns)
      double precision rscale,proj,rs
      double complex mpole(0:nterms,-nterms:nterms)
      double complex eye,dipstr(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1, ns
         call l3dformmp1_dp(ier,rscale,sources(1,i),dipstr(i),
     1        dipvec(1,i),center,nterms,mpole)
      enddo
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            mpole(l,m) = mpole(l,m)*rs
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_dp(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a dipole located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipstr  : complex dipole strength
c     dipvec  : dipole direction vector
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifrder,lfrder,ifr,lused,jer
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double precision dipvec(3)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex dipstr
c
c     compute workspace requirements
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifrder=iephi+lephi
      lfrder=2*(nterms+3)
c
      ifr=ifrder+lfrder
      lused=ifr + lfrder
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0_dp(jer,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole,w(ipp),w(ippd),w(iephi),
     2          w(ifr),w(ifrder))
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_dp(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,mpole,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformmp1_dp for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,i,nterms,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision dipvec(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision r,theta,phi,d
      double precision ctheta,stheta,cphi,sphi
      double precision phix,phiy,phiz
      double precision rx,ry,rz,thetax,thetay,thetaz
      double complex mpole(0:nterms,-nterms:nterms)
      double complex dipstr
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      double complex fruse,ux,uy,uz,ur,utheta,uphi,zzz
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = r
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. For n>0 modes,
c     the variable fruse is set to fr(n)/r:
c
c     
c
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
c
c     get the associated Legendre functions:
c
      call ylgndr2s(nterms,ctheta,pp,ppd)
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1))
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      mpole(0,0)= mpole(0,0) + zzz*dipstr
      do n=1,nterms
         fruse = fr(n-1)*rscale
         ur = pp(n,0)*frder(n)
         utheta = -fruse*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         mpole(n,0)= mpole(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fruse*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fruse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,m)= mpole(n,m) + zzz*dipstr
c
            ur = frder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fruse*ppd(n,m)
            uphi = eye*m*ephi(m)*fruse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,-m)= mpole(n,-m) + zzz*dipstr
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_dp(iffld,sources,dipstr,dipvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of dipoles at 
c     SOURCE(3,ns). 
c     
c     The potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      double precision sources(3,ns),target(3)
      double precision dipvec(3,ns)
      double complex pot,fld(3),potloc,fldloc(3)
      double complex eye
      double complex dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d_dp(iffld,sources(1,i),dipstr(i),dipvec(1,i),
     1        target,potloc,fldloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_dp_targ(iffld,sources,dipstr,dipvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of dipoles at 
c     SOURCE(3,ns). 
c     
c     The potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      double precision sources(3,ns),target(3)
      double precision dipvec(3,ns)
      double complex pot,fld(3),potloc,fldloc(3)
      double complex eye
      double complex dipstr(ns)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
      pot=pot+ dipstr(i)*(dotprod*dinv3)
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
      pot=pot+ dipstr(i)*(dotprod*dinv3)
c
ccc      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=fld(1)+dipstr(i)*(rtttt*xdiff-dinv3*dipvec(1,i))
         fld(2)=fld(2)+dipstr(i)*(rtttt*ydiff-dinv3*dipvec(2,i))
         fld(3)=fld(3)+dipstr(i)*(rtttt*zdiff-dinv3*dipvec(3,i))
ccc      endif 
c
      enddo
      endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3d_dp(iffld,source,dipstr,dipvec,target,
     1                        pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     source(3)    : location of the source 
c     dipstr(ns)   : dipole strength
c     dipvec(3,ns) : dipole direction
c     target(3)    : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot          : calculated potential
c     fld          : calculated -gradient
c
c----------------------------------------------------------------------
      implicit none
      integer iffld
      double precision source(3),target(3)
      double precision dipvec(3)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
      double complex pot,fld(3)
      double complex cd,eye,ztttt,cd2
      double complex dipstr,z1,z2,z3
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)+zdiff*dipvec(3)
      pot= dipstr*(dotprod*dinv3)
      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=dipstr*(rtttt*xdiff-dinv3*dipvec(1))
         fld(2)=dipstr*(rtttt*ydiff-dinv3*dipvec(2))
         fld(3)=dipstr*(rtttt*zdiff-dinv3*dipvec(3))
      endif 
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_dp(ier,rscale,sources,dipstr,dipvec,ns,
     1		           center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  deprecated but left in calling sequence for
c		  backward compatibility.
c
c     locexp    : coeffs for the expansion
c
c
c----------------------------------------------------------------------
      implicit none
      integer ier,ns,nterms,l,m,i
      double precision rscale,sources(3,ns),center(3),rs
      double precision dipvec(3,ns)
      double complex locexp(0:nterms,-nterms:nterms), dipstr(ns)
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1_dp(ier,rscale,sources(1,i),dipstr(i),
     1		dipvec(1,i),center,nterms,locexp)
      enddo
c
c
      do l = 0,nterms
         rs = sqrt(1.0d0/(2*l+1))
         do m=-l,l
            locexp(l,m) = locexp(l,m)*rs
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta1_dp(ier,rscale,source,dipstr,dipvec,
     &		center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single dipole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double precision dipvec(3)
      double complex locexp(0:nterms,-nterms:nterms), dipstr
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      ippd = ipp+lpp
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      ifrder=ifr+lfr
      lfrder=2*(nterms+3)
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dformta0_dp(rscale,source,dipstr,dipvec,
     &   center,nterms,locexp,w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder))
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0_dp(rscale,source,dipstr,dipvec,
     &		center,nterms,locexp,pp,ppd,ephi,fr,frder)
c**********************************************************************
c
c     See l3dformta_dp/l3dformta1_dp for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,i,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision dipvec(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision r,theta,phi,ctheta,stheta,cphi,sphi,done
      double precision d,phix,phiy,phiz,rx,ry,rz,thetax,thetay,thetaz
      double complex locexp(0:nterms,-nterms:nterms), dipstr
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      double complex ux,uy,uz,ur,utheta,uphi,zzz
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
      fr(nterms+1)=fr(nterms)*d
      do i=0,nterms
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
        rx = stheta*cphi
        thetax = ctheta*cphi/r
        phix = -sphi/r
        ry = stheta*sphi
        thetay = ctheta*sphi/r
        phiy = cphi/r
        rz = ctheta
        thetaz = -stheta/r
        phiz = 0.0d0
c
c
c     get the associated Legendre functions:
c
      call ylgndr2s(nterms,ctheta,pp,ppd)
c
c     Compute contribution to local coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |T|^n Ylm(T) Ylm*(S) / (|S|^{n+1})
c
c     so contribution is |S|^{n+1} times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      locexp(0,0)= locexp(0,0) + zzz*dipstr
      do n=1,nterms
         ur = pp(n,0)*frder(n)
         utheta = -fr(n)*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         locexp(n,0)= locexp(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fr(n)*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fr(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,m)= locexp(n,m) + zzz*dipstr
c
            ur = frder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fr(n)*ppd(n,m)
            uphi = eye*m*ephi(m)*fr(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,-m)= locexp(n,-m) + zzz*dipstr
         enddo
      enddo
c
      return
      end
c
c
c**********************************************************************
      subroutine lpotfld3dall_sdp_targ(iffld,sources,
     $     charge,dipstr,dipvec,ns,target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection 
c     of charges and dipoles at SOURCE(3,ns). 
c     
c     The potential due to a single charge is 
c
c        pot = charge/r 
c
c     and the potential due to a single dipole is 
c
c        pot = dipstr*(dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c
c     where (x,y,z) = target - source and r = sqrt(x^2+y^2+z^2).
c
c	 fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                   iffld = 0 -> dont compute 
c		           iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     charge(ns)    : charge strength
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,ns,i
      double precision sources(3,ns),target(3)
      double precision dipvec(3,ns)
      double complex pot,fld(3),potloc,fldloc(3)
      double complex eye
      double complex charge(ns),dipstr(ns)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,dinv2,
     $   dinv3,ddd,dotprod,dinv5,rtttt
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if( iffld .eq. 0 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)

      pot=pot+charge(i)*dinv
      pot=pot+dipstr(i)*(dotprod*dinv3)
c
      enddo
      endif
c
      if( iffld .eq. 1 ) then
      do i = 1,ns
c
      xdiff=target(1)-sources(1,i)
      ydiff=target(2)-sources(2,i)
      zdiff=target(3)-sources(3,i)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
      dinv = 1.0d0/d
      dinv2 = dinv*dinv
      dinv3 = dinv*dinv2
      dotprod = xdiff*dipvec(1,i)+ydiff*dipvec(2,i)+zdiff*dipvec(3,i)
c
      pot=pot+charge(i)*dinv
      pot=pot+dipstr(i)*(dotprod*dinv3)
c
ccc      if (iffld.eq.1) then
         dinv5 = dinv3*dinv2
         rtttt = 3.0d0*dotprod*dinv5
         fld(1)=fld(1)+charge(i)*xdiff*dinv3
         fld(2)=fld(2)+charge(i)*ydiff*dinv3
         fld(3)=fld(3)+charge(i)*zdiff*dinv3
         fld(1)=fld(1)+dipstr(i)*(rtttt*xdiff-dinv3*dipvec(1,i))
         fld(2)=fld(2)+dipstr(i)*(rtttt*ydiff-dinv3*dipvec(2,i))
         fld(3)=fld(3)+dipstr(i)*(rtttt*zdiff-dinv3*dipvec(3,i))
ccc      endif 
c
      enddo
      endif
c
      return
      end
c
c
c
c      The next set of routines contains accelerated versions of the
c      forming and evaluating multipole expansions.
c
c
c**********************************************************************
      subroutine l3dmpevalall_trunc(rscale,center,mpole,nterms,nterms1,
     $     ztarg,nt,ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to a TRUNCATED outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld = -gradient(pot) if iffld = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :   scaling parameter (see formmp1l3d)
c     center :   expansion center
c     mpole  :   multipole expansion in 2d matrix format
c     nterms :   order of the multipole expansion
c     nterms1 :   order of truncated expansion to be used
c     ztarg  :   target location
c     nt     :   number of targets
c     ifpot  :   flag controlling evaluation of potential
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.        
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.        
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     fld    :    gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused,i
      double precision rscale,center(3),ztarg(3,nt)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex pot(nt),fld(3,nt)
      double complex pot0,fld0(3)
      double complex mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+3)+5 
c
      ifr=iephi+lephi
      ifrder=ifr+(nterms+3)
      lused=ifrder+(nterms+3)
      allocate(w(lused))
c
      do i=1,nt
      call l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,
     $   ztarg(1,i),pot0,iffld,fld0,w(ipp),w(ippd),
     $     w(iephi),w(ifr),w(ifrder),wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
c
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval_trunc(rscale,center,mpole,nterms,nterms1,
     $     ztarg,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld = -gradient(pot) if iffld = 0.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     nterms1 :   order of truncated expansion to be used
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lused,i
      double precision rscale,center(3),ztarg(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex pot,fld(3)
      double complex mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+3)+5 
c
      ifr=iephi+lephi
      ifrder=ifr+(nterms+3)
      lused=ifrder+(nterms+3)
      allocate(w(lused))
c
      call l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,ztarg,
     1	   pot,iffld,fld,w(ipp),w(ippd),
     $     w(iephi),w(ifr),w(ifrder),wlege,nlege)
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dmpeval_trunc0(rscale,center,mpole,nterms,nterms1,
     1		ztarg,pot,iffld,fld,ynm,ynmd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dmpeval for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,iffld,ifpot,nlege
      integer n,m,i,l
      double precision rscale,center(3),ztarg(3),zdiff(3)
      double precision wlege(0:nlege,0:nlege)
      double precision ynm(0:nterms1,0:nterms1)
      double precision ynmd(0:nterms1,0:nterms1)
      double precision fr(0:nterms+1)
      double precision frder(0:nterms+1)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex pot,fld(3),ephi1,ur,utheta,uphi,ux,uy,uz
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1)
c
      double complex eye
      double complex ztmp1,ztmp2,ztmp3,ztmpsum,z
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      d = 1.0d0/r
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = 1/r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=conjg(ephi(i))
      enddo
      do i=0,nterms+1
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
      if (iffld.eq.1) then
         rx = stheta*cphi
         thetax = ctheta*cphi/r
         phix = -sphi/r
         ry = stheta*sphi
         thetay = ctheta*sphi/r
         phiy = cphi/r
         rz = ctheta
         thetaz = -stheta/r
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
ccc         call ylgndr2s(nterms1,ctheta,ynm,ynmd)
c         call ylgndr2sfw(nterms1,ctheta,ynm,ynmd,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               ynm(l,m) = ynm(l,m)*rs
c               ynmd(l,m) = ynmd(l,m)*rs
c            enddo
c         enddo
         call ylgndru2sfw(nterms1,ctheta,ynm,ynmd,wlege,nlege)
      else        
ccc         call ylgndr(nterms1,ctheta,ynm)
c         call ylgndrfw(nterms1,ctheta,ynm,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               ynm(l,m) = ynm(l,m)*rs
c            enddo
c         enddo
         call ylgndrufw(nterms1,ctheta,ynm,wlege,nlege)
      endif
c
c
c
c     initialize computed values and 
c     scale derivatives of Hankel functions so that they are
c     derivatives with respect to r.
c
      if (iffld.eq.1) then
         ur = mpole(0,0)*frder(0)
         utheta = 0.0d0
         uphi = 0.0d0
      endif
      pot=mpole(0,0)*fr(0)
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms1
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    ur = ur + frder(n)*ynm(n,0)*mpole(n,0)
	    utheta = utheta -mpole(n,0)*fr(n)*ynmd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)*stheta
	       ztmp2 = mpole(n,m)*ephi(m) 
	       ztmp3 = mpole(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*ynm(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fr(n)*ynmd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fr(n)*ynm(n,m)*ztmpsum
            enddo
         enddo
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms1
	    pot=pot+mpole(n,0)*fr(n)*ynm(n,0)
	    do m=1,n
	       ztmp1=fr(n)*ynm(n,m)
	       ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine l3dtaevalall_trunc(rscale,center,locexp,nterms,nterms1,
     1		ztarg,nt,ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a TRUNCATED local expansion centered 
c     at CENTER at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the local expansion
c     nterms     : order of the local expansion
c     nterms1    : order of the truncated expansion to be used
c     ztarg      : vector of targets 
c     nt         : number of targets
c     ifpot      : flag for potential computation
c		                    ifpot=0  - pot is not computed
c		                    ifpot=1  - pot is computed
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential aat ztarg (if requested)
c     fld(3)     : gradient at ztarg (if requested)
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      double precision rscale,center(3),ztarg(3,nt)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex pot(nt),fld(3,nt)
      double complex pot0,fld0(3)
      double complex locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr= nterms+3
c
      ifrder=ifr+lfr
      lfrder=nterms+3
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      do i=1,nt
      call l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,
     $   ztarg(1,i),pot0,iffld,fld0,w(ipp),w(ippd),w(iephi),w(ifr),
     2   w(ifrder),wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval_trunc(rscale,center,locexp,nterms,nterms1,
     1		ztarg,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a TRUNCATED local expansion centered 
c     at CENTER at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) r^n Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale     : scaling parameter used in forming expansion
c                                   (see l3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the local expansion
c     nterms     : order of the local expansion
c     nterms     : order of the truncated expansion to be used
c     ztarg      : target vector
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot        : potential at ztarg(3)
c     fld(3)     : gradient at ztarg (if requested)
c     ier        : error return code
c		      ier=0	returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,iffld,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      double precision rscale,center(3),ztarg(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex pot,fld(3)
      double complex locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr= nterms+3
c
      ifrder=ifr+lfr
      lfrder=nterms+3
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,ztarg,
     1	     pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifr),
     2       w(ifrder),wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaeval_trunc0(rscale,center,locexp,nterms,nterms1,
     1		ztarg,pot,iffld,fld,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dtaeval for comments.
c     (pp and ppd are storage arrays for Ynm and Ynm')
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,iffld,nlege
      integer i,l,n,m
      double precision wlege(0:nlege,0:nlege)
      double precision rscale,center(3),ztarg(3),zdiff(3)
      double precision pp(0:nterms1,0:nterms1)
      double precision ppd(0:nterms1,0:nterms1)
      double precision fruse,fr(0:nterms+1),frder(0:nterms+1)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex pot,fld(3),ephi1,ephi1inv
      double complex locexp(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1)
c
      double complex eye,ur,utheta,uphi
      double complex ztmp,z
      double complex ztmp1,ztmp2,ztmp3,ztmpsum
      double complex ux,uy,uz
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
c     Convert to spherical coordinates
c
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      d = rscale*r
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = rscale*r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
c
c     compute e^{eye*m*phi} array.
c
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=conjg(ephi(i))
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. 
c
c     
c
      if (iffld.eq.1) then
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
c         call ylgndr2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               pp(l,m) = pp(l,m)*rs
c               ppd(l,m) = ppd(l,m)*rs
c            enddo
c         enddo
         call ylgndru2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
      else
c         call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
c         do l = 0,nterms1
c            rs = sqrt(1.0d0/(2*l+1))
c            do m=0,l
c               pp(l,m) = pp(l,m)*rs
c            enddo
c         enddo
         call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
      endif
c
c
      pot=locexp(0,0)*fr(0)
      if (iffld.eq.1) then
         ur = 0.0d0
         utheta = 0.0d0
         uphi = 0.0d0
      endif
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms1
            pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    ur = ur + frder(n)*pp(n,0)*locexp(n,0)
	    fruse = fr(n-1)*rscale
	    utheta = utheta -locexp(n,0)*fruse*ppd(n,0)*stheta
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)*stheta
	       ztmp2 = locexp(n,m)*ephi(m) 
	       ztmp3 = locexp(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + frder(n)*pp(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fruse*ppd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fruse*pp(n,m)*ztmpsum
            enddo
         enddo
ccc	 call prin2(' ur is *',ur,2)
ccc	 call prin2(' utheta is *',utheta,2)
ccc	 call prin2(' uphi is *',uphi,2)
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms1
	    pot=pot+locexp(n,0)*fr(n)*pp(n,0)
	    do m=1,n
	       ztmp1=fr(n)*pp(n,m)
	       ztmp2 = locexp(n,m)*ephi(m)+locexp(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_trunc(ier,rscale,sources,charge,ns,center,
     1                  nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     nterms1         : order of truncated expansion to be generated
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		         ier=0  returned successfully
c    		         deprecated but left in calling sequence for
c		         backward compatibility.
c    
c     mpole           : coeffs of the multipole expansion
c                  
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,nlege
      integer ns,i,l,m,ier1
      double precision center(3),sources(3,ns)
      double precision rscale,rs
      double precision wlege(0:nlege,0:nlege)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex eye,charge(ns)
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision, allocatable :: w(:)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
c
c ... Assign work spaces:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      do i = 1, ns
c         call l3dformmp_trunc1
c     $   (ier1,rscale,sources(1,i),charge(i),center,
c     1        nterms,nterms1,mpole,wlege,nlege)
        call l3dformmp_trunc0(rscale,sources(1,i),charge(i),center,
     $   nterms,nterms1,
     1   mpole,w(ipp),w(iephi),w(ifr),wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            mpole(l,m) = mpole(l,m)*rs
c         enddo
c      enddo
c
      return
      end
C
C***********************************************************************
      subroutine l3dformmp_add_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs TRUNCATED multipole expansion about CENTER due to 
C     NS sources located at SOURCES(3,*)  and INCREMENTS mpole
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		           ier=0  returned successfully
c		           deprecated but left in calling sequence for
c		           backward compatibility.
c    
c     mpole           : incremented multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer ns,i,l,m,ier,ier1
      double precision center(3),sources(3,ns)
      double precision rscale
      double precision wlege(0:nlege,0:nlege)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex eye,charge(ns)
      double complex, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
C
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
        call l3dformmp_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)
c
        do l = 0,nterms
          do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
          enddo
        enddo
c
      return
      end
c
C
c**********************************************************************
      subroutine l3dformmp_trunc1(ier,rscale,source,charge,center,
     1		nterms,nterms1,mpole,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a charge located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : complex charge strength
c     center  : coordinates of the expansion center
c     nterms  : order of the multipole expansion
c     nterms1 : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c                            
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,lused,i
      double precision rscale,source(3),center(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex charge
c
c ... Assign work spaces:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp_trunc0(rscale,source,charge,center,
     $   nterms,nterms1,
     1   mpole,w(ipp),w(iephi),w(ifr),wlege,nlege)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp_trunc0(rscale,source,charge,center,
     1		nterms,nterms1,mpole,pp,ephi,fr,wlege,nlege)
c**********************************************************************
c
c     See l3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer n,m,i
      double precision rscale,source(3),center(3),zdiff(3)
      double precision wlege(0:nlege,0:nlege)
      double precision pp(0:nterms1,0:nterms1)
      double precision done,r,theta,phi,dtmp
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex mpole(0:nterms,-nterms:nterms)
      double complex charge
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex  fr(0:nterms+1)
      double complex  ztmp,z
c
c
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      d = r
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      d = r
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=conjg(ephi(i))
      enddo
c
c     get the associated Legendre functions:
c
ccc      call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
      call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     multiply all fr's by charge strength.
c
      do n = 0,nterms1
         fr(n) = fr(n)*charge
      enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r =  
c          \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)  / (|T|)^{n+1}
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      mpole(0,0)= mpole(0,0) + fr(0)
      do n=1,nterms1
         mpole(n,0)= mpole(n,0) + pp(n,0)*fr(n)
         do m=1,n
cc            ztmp=pp(n,m)*fr(n)
cc            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
cc            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
            ztmp=pp(n,m)*fr(n)
            mpole(n, m)= mpole(n, m) + ztmp*ephi(-m)
            mpole(n,-m)= mpole(n,-m) + ztmp*ephi(+m)
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta_trunc(ier,rscale,sources,charge,ns,center,
     1		        nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1    : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer ns,i,l,m
      double precision rscale,rs,sources(3,ns),center(3)
      double precision wlege(0:nlege,0:nlege)
      double complex locexp(0:nterms,-nterms:nterms), charge(ns)
      double complex eye
      integer ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision, allocatable :: w(:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      do i = 1,ns
c         call l3dformta_trunc1(ier,rscale,sources(1,i),charge(i),
c     1		center,nterms,nterms1,locexp,wlege,nlege)
         call l3dformta_trunc0(rscale,sources(1,i),charge(i),center,
     &      nterms,nterms1,locexp,w(ipp),w(iephi),w(ifr),
     $      wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            locexp(l,m) = locexp(l,m)*rs
c         enddo
c      enddo
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_add_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below. INCREMENT
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c
c-----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : incremented local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer ns,l,m
      double precision rscale,sources(3,ns),center(3)
      double precision wlege(0:nlege,0:nlege)
      double complex locexp(0:nterms,-nterms:nterms), charge(ns)
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
      double complex, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
      call l3dformta_trunc
     $     (ier,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m) + mptemp(l,m)
         enddo
      enddo
c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta_trunc1(ier,rscale,source,charge,center,
     &		nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single charge located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c     source    : coordinates of the source
c     charge    : coordinates of the source
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   deprecated but left in calling sequence for
c		   backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nt,ier,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,lused
      double precision rscale,source(3),center(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double complex locexp(0:nterms,-nterms:nterms), charge
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      call l3dformta_trunc0(rscale,source,charge,center,
     &   nterms,nterms1,locexp,w(ipp),w(iephi),w(ifr),
     $   wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta_trunc0(rscale,source,charge,
     &		center,nterms,nterms1,locexp,pp,ephi,fr,wlege,nlege)
c**********************************************************************
c
c     See l3dformta/l3dformta1 for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer i,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision pp(0:nterms1,0:nterms1)
      double precision wlege(0:nlege,0:nlege)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex fr(0:nterms+1)
      double complex locexp(0:nterms,-nterms:nterms), charge
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex ztmp,z
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c      call l3drhpolar(zdiff(1),zdiff(2),zdiff(3),r,ctheta,ephi1)
c      stheta=sqrt(done-ctheta*ctheta)
c      cphi = dcos(phi)
c      sphi = dsin(phi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=conjg(ephi1)
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=conjg(ephi(i))
      enddo
c
c     get the Ynm
c
      call ylgndrufw(nterms1,ctheta,pp,wlege,nlege)
c
c     compute radial functions and scale them by charge strength.
c
      do n = 0, nterms1
         fr(n) = fr(n)*charge
      enddo
c
c     Compute contributions to locexp
c
      locexp(0,0)=locexp(0,0) + fr(0)
      do n=1,nterms1
         locexp(n,0)=locexp(n,0) + pp(n,0)*fr(n)
         do m=1,n
            ztmp=pp(n,m)*fr(n)
	    locexp(n,m)=locexp(n,m) + ztmp*ephi(-m)
	    locexp(n,-m)=locexp(n,-m) + ztmp*ephi(m)
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_dp_trunc(ier,rscale,sources,dipstr,dipvec,ns,
     1                  center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs TRUNCATED multipole expansion about CENTER due to NS 
c     dipole sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          deprecated but left in calling sequence for
c		          backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,i,l,m, ier,nlege
      double precision center(3),sources(3,ns)
      double precision dipvec(3,ns)
      double precision rscale,rs
      double precision wlege(0:nlege,0:nlege)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex eye,dipstr(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1, ns
         call l3dformmp1_dp_trunc(ier,rscale,sources(1,i),dipstr(i),
     1        dipvec(1,i),center,nterms,nterms1,mpole,wlege,nlege)
      enddo
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            mpole(l,m) = mpole(l,m)*rs
c         enddo
c      enddo
c
      return
      end
C
C***********************************************************************
      subroutine l3dformmp_dp_add_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     dipole sources located at SOURCES(3,*) and INCREMENTS mpole.
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     nterms1         : order of truncated multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          deprecated but left in calling sequence for
c		          backward compatibility.
c
c     mpole           : incremented multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,i,l,m, ier,nlege
      double precision center(3),sources(3,ns)
      double precision dipvec(3,ns)
      double precision rscale
      double precision wlege(0:nlege,0:nlege)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex, allocatable :: mptemp(:,:)
      double complex eye,dipstr(ns)
      data eye/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo

      call l3dformmp_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
         enddo
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_dp_trunc(ier,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the truncated multipole expansion 
c     about CENTER due to a dipole located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipstr  : complex dipole strength
c     dipvec  : dipole direction vector
c     center  : coordinates of the expansion center
c     nterms  : order of the multipole expansion
c     nterms1 : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,jer,iffld,ifpot,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,ifrder,lfrder,lused
      double precision rscale,source(3),center(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double precision dipvec(3)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex dipstr
c
c ... Assign work spaces:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifrder=iephi+lephi
      lfrder=2*(nterms+3)
c
      ifr=ifrder+lfrder
      lused=ifr + lfrder
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0_dp_trunc(jer,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,w(ipp),w(ippd),w(iephi),
     2          w(ifr),w(ifrder),wlege,nlege)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_dp_trunc(ier,rscale,source,dipstr,dipvec,
     1    center,nterms,nterms1,mpole,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dformmp1_dp for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer i,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision rfac1, rfac2, rfac3
      double precision dipvec(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision wlege(0:nlege,0:nlege)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex mpole(0:nterms,-nterms:nterms)
      double complex dipstr
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      double complex fruse,ux,uy,uz,ur,utheta,uphi,zzz
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polarl(zdiff,r,theta,phi)
      d = r
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      fr(0) = 1.0d0
      d = d*rscale
      fr(1) = d
      do i=2,nterms+1
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
      frder(0) = 0.0d0
      do i=1,nterms+1
         frder(i) = i*fr(i-1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. For n>0 modes,
c     the variable fruse is set to fr(n)/r:
c
c     
c
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
c
c     get the associated Legendre functions:
c
ccc      call ylgndr2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
      call ylgndru2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1))
c
c     so contribution is |S|^n times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      rfac1 = dipvec(1)*rx + dipvec(2)*ry + dipvec(3)*rz
      rfac2 = dipvec(1)*thetax + dipvec(2)*thetay + dipvec(3)*thetaz
      rfac3 = dipvec(1)*phix + dipvec(2)*phiy + dipvec(3)*phiz
ccc      ux = ur*rx + utheta*thetax + uphi*phix
ccc      uy = ur*ry + utheta*thetay + uphi*phiy
ccc      uz = ur*rz + utheta*thetaz + uphi*phiz
ccc      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
      mpole(0,0)= mpole(0,0) + zzz*dipstr
c
      do n=1,nterms1
         fruse = fr(n-1)*rscale
         ur = pp(n,0)*frder(n)
         utheta = -fruse*ppd(n,0)*stheta
         uphi = 0.0d0
ccc         ux = ur*rx + utheta*thetax + uphi*phix
ccc         uy = ur*ry + utheta*thetay + uphi*phiy
ccc         uz = ur*rz + utheta*thetaz + uphi*phiz
ccc         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
         mpole(n,0)= mpole(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta
            utheta = -fruse*ppd(n,m)
            uphi = -eye*m*fruse*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = (rfac1*ur + rfac2*utheta + rfac3*uphi)*ephi(-m)
            mpole(n,m)= mpole(n,m) + zzz*dipstr
c
ccc            ur = frder(n)*pp(n,m)*stheta*ephi(m)
ccc            utheta = -ephi(m)*fruse*ppd(n,m)
ccc            uphi = eye*m*ephi(m)*fruse*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
ccc            zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
            zzz = conjg(zzz) 
            mpole(n,-m)= mpole(n,-m) + zzz*dipstr
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dformta_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the local expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ns,ier,nlege
      integer i,l,m
      double precision sources(3,ns),center(3),rs,rscale
      double precision dipvec(3,ns)
      double precision wlege(0:nlege,0:nlege)
      double complex locexp(0:nterms,-nterms:nterms), dipstr(ns)
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1_dp_trunc(ier,rscale,sources(1,i),dipstr(i),
     1      dipvec(1,i),center,nterms,nterms1,locexp,wlege,nlege)
      enddo
c
c
c      do l = 0,nterms
c         rs = sqrt(1.0d0/(2*l+1))
c         do m=-l,l
c            locexp(l,m) = locexp(l,m)*rs
c         enddo
c      enddo
c
      return
      end
c
c
c**********************************************************************
      subroutine l3dformta_dp_add_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*)
c     and INCREMENTS locexp.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1/l3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c
c     locexp    : coeffs for the expansion
c----------------------------------------------------------------------
      implicit none
      integer nterms,ns,nterms1,ier,nlege
      integer l,m
      double precision rscale,sources(3,ns),center(3)
      double precision dipvec(3,ns)
      double precision wlege(0:nlege,0:nlege)
      double complex locexp(0:nterms,-nterms:nterms), dipstr(ns)
      double complex eye
      double complex, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c        do l = 0,nterms
c          do m=-l,l
c             mptemp(l,m) = 0
c          enddo
c        enddo
c
      call l3dformta_dp_trunc
     $     (ier,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)+mptemp(l,m)
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta1_dp_trunc(ier,rscale,source,dipstr,dipvec,
     &		center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the truncated local expansion about 
c     CENTER due to a single dipole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the local expansion
c     nterms1   : order of the truncated expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  :   dimension parameter for wlege
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,ier,nlege
      integer lpp,ipp,ippd,iephi,lephi,ifr,lfr,ifrder,lfrder,lused,i
      double precision rscale,source(3),center(3)
      double precision wlege(0:nlege,0:nlege)
      double precision, allocatable :: w(:)
      double precision dipvec(3)
      double complex locexp(0:nterms,-nterms:nterms), dipstr
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      ippd = ipp+lpp
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+3)
c
      ifrder=ifr+lfr
      lfrder=2*(nterms+3)
c
      lused=ifrder+lfrder
      allocate(w(lused))
c
      call l3dformta0_dp_trunc(rscale,source,dipstr,dipvec,
     &   center,nterms,nterms1,locexp,
     $   w(ipp),w(ippd),w(iephi),w(ifr),w(ifrder),wlege,nlege)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformta0_dp_trunc(rscale,source,dipstr,dipvec,
     &     center,nterms,nterms1,
     $     locexp,pp,ppd,ephi,fr,frder,wlege,nlege)
c**********************************************************************
c
c     See l3dformta_dp_trunc/l3dformta1_dp_trunc for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,nterms1,nlege
      integer i,n,m
      double precision rscale,source(3),center(3),zdiff(3)
      double precision rfac1, rfac2, rfac3
      double precision dipvec(3)
      double precision pp(0:nterms,0:nterms)
      double precision ppd(0:nterms,0:nterms)
      double precision wlege(0:nlege,0:nlege)
      double precision done,r,theta,phi
      double precision ctheta,stheta,cphi,sphi
      double precision d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      double complex locexp(0:nterms,-nterms:nterms), dipstr
      double complex ephi(-nterms:nterms),ephi1,ephi1inv
      double complex fr(0:nterms+1),ztmp,frder(0:nterms+1),z
      double complex ux,uy,uz,ur,utheta,uphi,zzz
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      d = 1.0d0/r
      fr(0) = d
      d = d/rscale
      fr(1) = fr(0)*d
      do i=2,nterms
         fr(i) = fr(i-1)*d
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
      fr(nterms+1)=fr(nterms)*d
      do i=0,nterms
         frder(i) = -(i+1.0d0)*fr(i+1)*rscale
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
        rx = stheta*cphi
        thetax = ctheta*cphi/r
        phix = -sphi/r
        ry = stheta*sphi
        thetay = ctheta*sphi/r
        phiy = cphi/r
        rz = ctheta
        thetaz = -stheta/r
        phiz = 0.0d0
c
c
c     get the associated Legendre functions:
c
ccc      call ylgndr2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
      call ylgndru2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
c
c     Compute contribution to local coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        1/r = 
c         \sum_n 1/(2n+1) \sum_m  |T|^n Ylm(T) Ylm*(S) / (|S|^{n+1})
c
c     so contribution is |S|^{n+1} times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c
      ur = pp(0,0)*frder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      rfac1 = dipvec(1)*rx + dipvec(2)*ry + dipvec(3)*rz
      rfac2 = dipvec(1)*thetax + dipvec(2)*thetay + dipvec(3)*thetaz
      rfac3 = dipvec(1)*phix + dipvec(2)*phiy + dipvec(3)*phiz
ccc      ux = ur*rx + utheta*thetax + uphi*phix
ccc      uy = ur*ry + utheta*thetay + uphi*phiy
ccc      uz = ur*rz + utheta*thetaz + uphi*phiz
ccc      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
      locexp(0,0)= locexp(0,0) + zzz*dipstr
      do n=1,nterms
         ur = pp(n,0)*frder(n)
         utheta = -fr(n)*ppd(n,0)*stheta
         uphi = 0.0d0
ccc         ux = ur*rx + utheta*thetax + uphi*phix
ccc         uy = ur*ry + utheta*thetay + uphi*phiy
ccc         uz = ur*rz + utheta*thetaz + uphi*phiz
ccc         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         zzz = rfac1*ur + rfac2*utheta + rfac3*uphi
         locexp(n,0)= locexp(n,0) + zzz*dipstr
         do m=1,n
            ur = frder(n)*pp(n,m)*stheta
            utheta = -fr(n)*ppd(n,m)
            uphi = -eye*m*fr(n)*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = (rfac1*ur + rfac2*utheta + rfac3*uphi)*ephi(-m)
            locexp(n,m)= locexp(n,m) + zzz*dipstr
c
ccc            ur = frder(n)*pp(n,m)*stheta*ephi(m)
ccc            utheta = -ephi(m)*fr(n)*ppd(n,m)
ccc            uphi = eye*m*ephi(m)*fr(n)*pp(n,m)
ccc            ux = ur*rx + utheta*thetax + uphi*phix
ccc            uy = ur*ry + utheta*thetay + uphi*phiy
ccc            uz = ur*rz + utheta*thetaz + uphi*phiz
ccc            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            zzz = conjg(zzz) 
            locexp(n,-m)= locexp(n,-m) + zzz*dipstr
         enddo
      enddo
c
      return
      end
c
c
C***********************************************************************
c
c
c       Multipole forming routines for real-valued charges and dipoles
c
c
C***********************************************************************
      subroutine l3dformmp_charge_trunc(ier,rscale,sources,charge,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS real-valued
c     charge sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      double precision center(3),sources(3,ns)
      double precision charge(ns)
      double precision rscale
      double complex mpole(0:nterms,-nterms:nterms)
      integer nlege
      double precision wlege(*)
      double precision, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      lused=ifr + lfr
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

      do i = 1, ns

      call l3dformmp0_charge_trunc(ier,rscale,sources(1,i),charge(i),
     1   center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),w(ifr))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp0_charge_trunc(ier,rscale,source,charge,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a charge located at the point SOURCE.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : charge strengths
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c     sss     : sign mutiplier array (see getsgnformpmp)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      double precision rscale,source(3),center(3),zdiff(3)
      double precision charge
      double precision pp(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision cscale,cscale1,cscale2,cscale3,rtmp
      double precision a22,a21,a20,rmul,dd, a11,a10
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1),ephi1
      double complex eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      double precision wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert dipole vector contributions to standard
c     n=0 moments about source position using standard
c     d+,d-,dz operators.
c
c     now shift n=1 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi(1))
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
        do i=0,nterms
        fr(i)=fr(i)*charge
        enddo
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
c       ... optimized for speed
c       use symmetries for real valued charges
c
            do ll = 0,nterms

               mpole(ll,0) = mpole(ll,0) + fr(ll)*pp(ll,0)

               do mm = 1,ll

               rtmp=fr(ll)*pp(ll,mm)
               mpole(ll, mm)= mpole(ll, mm) + rtmp*ephi(-mm)
               mpole(ll,-mm)= mpole(ll,-mm) + rtmp*ephi(+mm)

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
      subroutine l3dformmp_dipole_trunc(ier,rscale,sources,dipvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS real-valued
c     dipole sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(3,ns)    : dipole vector directions
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, isss, lsss
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      double precision center(3),sources(3,ns)
      double precision dipvec(3,ns)
      double precision rscale
      double complex mpole(0:nterms,-nterms:nterms)
      double precision binom(0:120,0:2)
      integer nlege
      double precision wlege(*)
      double precision, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      isss=ifr + lfr
      lsss = 3*(2*nterms+1)
c
      lused=isss + lsss
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

ccc      call getsgnformpmp_dipole(w(isss),nterms)

      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo


      do i = 1, ns

      call l3dformmp0_dipole_trunc(ier,rscale,sources(1,i),dipvec(1,i),
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,w(isss))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp0_dipole_trunc(ier,rscale,source,dipvec,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr,binom,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a dipole located at the point SOURCE.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipvec  : dipole vector
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c     pp     : work array for Ynm
c     ephi   : work array for Ynm
c     binom  : precomputed coeff array (sqrt of some binomial coeffs)
c     sss     : sign mutiplier array (see getsgnformpmp)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      double precision binom(0:120,0:2)
      double precision rscale,source(3),center(3),zdiff(3)
      double precision dipvec(3)
      double precision pp(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision cscale,cscale1,cscale2,cscale3,rtmp
      double precision sss(-1:1,-nterms:nterms)
      double precision a22,a21,a20,rmul,dd, a11,a10
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double complex mp1(-1:1)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1),ephi1
      double complex eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      double precision wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert dipole vector contributions to standard
c     n=1 moments about source position using standard
c     d+,d-,dz operators.
c
      mp1(+1) = (-dipvec(1) + dipvec(2)*eye) /sqrt(2.0d0)
      mp1(0) = dipvec(3)
      mp1(-1) = dconjg(mp1(+1))
c
c     now shift n=1 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi(1))
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
c      do l = 0,nterms
c      fr(l)=fr(l)/sqrt(2*l+1.0d0) *rscale
c      enddo
c      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
      do l = 0,nterms
      fr(l)=fr(l) *rscale
      enddo
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
        if( 1 .eq. 2 ) then
c
c       ... reference code
c
         do m = -1,1
            do ll = 0,nterms-1
               cscale = fr(ll)
               lnew = 1+ll
c
               cscale2 = binom(lnew-m,1-m)*binom(lnew+m,1+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            pp(ll,0)*cscale2*mp1(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,1-m)*binom(lnew+mnew1,1+m)
                  cscale2 = cscale2*cscale*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,1-m)*binom(lnew+mnew2,1+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp1(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp1(m)

               enddo
            enddo
         enddo
         endif
c
c
        if( 2 .eq. 2 ) then
c
c       ... optimized for speed
c       use symmetries for real valued dipoles
c
            do ll = 0,nterms-1

               lnew = 1+ll
c
               cscale2 = binom(lnew,1)*binom(lnew,1)
               mpole(lnew,0) = mpole(lnew,0) +
     1            fr(ll)*pp(ll,0)*cscale2*mp1(0)

               rtmp = fr(ll)*pp(ll,0)
               cscale2 = binom(lnew-1,0)*binom(lnew+1,2)
               ztmp = rtmp*cscale2*mp1(1)
               mpole(lnew,+1) = mpole(lnew,+1)+ztmp
               mpole(lnew,-1) = mpole(lnew,-1)+dconjg(ztmp)

               do mm = 1,ll

               ctmp=fr(ll)*pp(ll,mm)*ephi(-mm)

               mnew1 = +mm
               mnew2 = -mm
               cscale2 = binom(lnew-mnew1,1)*binom(lnew+mnew1,1)
               ztmp = cscale2*ctmp*mp1(0)
               mpole(lnew,mnew1) = mpole(lnew,mnew1)+ztmp
               mpole(lnew,mnew2) = mpole(lnew,mnew2)+dconjg(ztmp)

               do m = 1,1
c
c       ... skip zero valued modes
c       
                  if( abs(mp1(m)) .eq. 0 ) cycle

                  mnew1 = m+mm
                  mnew2 = m-mm

ccc                  cscale2 = +binom(lnew-mnew1,0)*binom(lnew+mnew1,2)
ccc                  cscale3 = -binom(lnew-mnew2,0)*binom(lnew+mnew2,2)
                  cscale2 = +binom(lnew+mnew1,2)
                  cscale3 = -binom(lnew+mnew2,2)

                  ztmp0=cscale2*mp1(m)*ctmp
                  ztmp1=cscale3*mp1(m)*dconjg(ctmp)

                  mpole(lnew,+mnew1) = mpole(lnew,+mnew1)+ztmp0
                  mpole(lnew,+mnew2) = mpole(lnew,+mnew2)+ztmp1

                  mpole(lnew,-mnew2) = mpole(lnew,-mnew2)+dconjg(ztmp1)
                  mpole(lnew,-mnew1) = mpole(lnew,-mnew1)+dconjg(ztmp0)

               enddo
               enddo
            enddo
         endif
c
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine getsgnformpmp_dipole(sss,nterms)
C***********************************************************************
c
c     This subroutine creates a multiplier that holds the 
c     appropriate power of (-1) that arises in multipole-multipole shifts.
c     (See Greengard - 
c          Rapid Evaluation of Potential Fields... for notation).
c-------------------------------------------------------------------------
      implicit none
      integer m,mm,nterms
      double precision sss(-1:1,-nterms:nterms)
c
      do m = -1,1
      do mm = -nterms,nterms
         sss(m,mm) = 1.0d0
      enddo
      enddo
      do m = -1,1
      do mm = -nterms,nterms
         if (( m .lt. 0) .and. (mm .gt. 0)) then
            if ( mm .le. -m ) sss(m,mm) = (-1)**mm
            if ( mm .gt. -m ) sss(m,mm) = (-1)**m
         endif
         if (( m .gt. 0) .and. (mm .lt. 0)) then
            if ( m .le. -mm ) sss(m,mm) = (-1)**m
            if ( m .gt. -mm ) sss(m,mm) = (-1)**mm
         endif
      enddo
      enddo
      return
      end
c
c
c
c
c    OBSOLETE
c**********************************************************************
ccc      subroutine l3dmpevalhess(rscale,center,mpole,nterms,ztarg,
ccc     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
ccc      implicit none
ccc      integer ier,nterms,iffld,ifhess
ccc      integer nterms2,iloc,lloc,lused
ccc      double precision rscale,center(3),ztarg(3)
ccc      double precision, allocatable :: w(:)
ccc      double complex pot,fld(3),hess(6)
ccc      double complex mpole(0:nterms,-nterms:nterms)
c
ccc      ier=0
c
c     Carve up workspace:
c
c     for local expansion
c
ccc      nterms2=0
ccc      if (iffld.eq.1) nterms2=1
ccc      if (ifhess.eq.1) nterms2=2
ccc      iloc=1
ccc      lloc=2*(nterms2+1)*(2*nterms2+1)+5
c
ccc      lused=iloc+lloc
ccc      allocate(w(lused))
c
ccc      call l3dmplocquadu(rscale,center,mpole,nterms,rscale,
ccc     1	   ztarg,w(iloc),nterms2,ier)
c
ccc      call l3devalhessloc(rscale,w(iloc),nterms2,
ccc     $   pot,iffld,fld,ifhess,hess)
c
ccc      return
ccc      end
c
c
c
c
c
c     OBSOLETE
c**********************************************************************
ccc      subroutine l3dtaevalhess(rscale,center,mpole,nterms,ztarg,
ccc     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to a local multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  r^{n}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
ccc      implicit none
ccc      integer ier,nterms,iffld,ifhess
ccc      integer nterms2,iloc,lloc,lused
ccc      double precision rscale,center(3),ztarg(3)
ccc      double precision, allocatable :: w(:)
ccc      double complex pot,fld(3),hess(6)
ccc      double complex mpole(0:nterms,-nterms:nterms)
c
ccc      ier=0
c
c     Carve up workspace:
c
c     for local expansion
c
ccc      nterms2=0
ccc      if (iffld.eq.1) nterms2=1
ccc      if (ifhess.eq.1) nterms2=2
ccc      iloc=1
ccc      lloc=2*(nterms2+1)*(2*nterms2+1)+5
cccc
ccc      lused=iloc+lloc
ccc      allocate(w(lused))
cccc
ccc      call l3dloclocquadu(rscale,center,mpole,nterms,rscale,
ccc     1	   ztarg,w(iloc),nterms2,ier)
c
ccc      call l3devalhessloc(rscale,w(iloc),nterms2,
ccc     $   pot,iffld,fld,ifhess,hess)
cccc
ccc      return
ccc      end
c
c
c
c
c
c**********************************************************************
ccc      subroutine l3devalhessloc(rscale,local,nterms,
ccc     $     pot,iffld,fld,ifhess,hess)
c**********************************************************************
c
c     Locally u = local(0,0) +
c                 local(1,-1) * r  * Y_{1,-1} +
c                 local(1,0)  * r  * Y_{1,0}  +
c                 local(1,1)  * r  * Y_{1,1}  +
c                 local(2,-2) * r2 * Y_{2,-2} +
c                 local(2,-1) * r2 * Y_{2,-1} +
c                 local(2,0)  * r2 * Y_{2,0}  +
c                 local(2,1)  * r2 * Y_{2,1}  +
c                 local(2,2)  * r2 * Y_{2,2}  +
c     so evaluation is easy (depending only on scaling convention
c     for definitions of Y_{n,m}. 
c
c---------------------------------------------------------------------
ccc      implicit none
ccc      integer nterms,iffld,ifhess
ccc      double precision rscale,pi,rfac
ccc      double complex local(0:nterms,-nterms:nterms)
ccc      double complex pot,fld(3),hess(6),eye,z0
c
ccc      eye = dcmplx(0.0d0,1.0d0)
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
ccc      pot = local(0,0)
c
c     fld comes from n=1 modes
c
ccc      if( iffld .eq. 1 ) then
ccc      rfac = sqrt(2.0d0)*rscale
ccc      fld(1) = rfac*(local(1,1) + local(1,-1))/2.0d0
ccc      fld(2) = rfac*eye*(local(1,1) - local(1,-1))/2.0d0
ccc      fld(3) = -local(1,0)*rscale
ccc      endif
c
c     hess comes from n=2 modes
c
ccc      if( ifhess .eq. 1 ) then
ccc      rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
ccc      z0 = rscale*rscale*local(2,0)
ccc      hess(1) = rfac*(local(2,2) + local(2,-2)) - z0
ccc      hess(2) = -rfac*(local(2,2) + local(2,-2)) - z0
ccc      hess(3) = 2*z0
ccc      hess(4) = rfac*eye*(local(2,2) - local(2,-2))
ccc      hess(5) = -rfac*(local(2,1) + local(2,-1))
ccc      hess(6) = -rfac*eye*(local(2,1) - local(2,-1))
ccc      endif
cccc
ccc      return
ccc      end
c
c
c
c
c**********************************************************************
      subroutine l3dmpevalhess(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to an outgoing multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,iffld,ifhess
      integer nterms2,iloc,lloc,lused
      double precision rscale,center(3),ztarg(3)
      double precision, allocatable :: scarray(:)
      double complex pot,fld(3),hess(6)
      double complex mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
        allocate( scarray(10*(nterms+2)**2) )
        call l3dmpevalhessdini(nterms,scarray)
c
        call l3dmpevalhessd(rscale,center,mpole,nterms,ztarg,
     1            pot,iffld,fld,ifhess,hess,scarray)
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine l3dtaevalhess(rscale,center,mpole,nterms,ztarg,
     1		pot,iffld,fld,ifhess,hess,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the 
c     potential due to a local multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  r^{n}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c     hess   :    Hessian at ztarg (if requested)
c                 in order (dxx,dyy,dzz,dxy,dxz,dyz)
c     ier    :    error return code
c		      ier=0  successful execution
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms,iffld,ifhess
      integer nterms2,iloc,lloc,lused
      double precision rscale,center(3),ztarg(3)
      double precision, allocatable :: scarray(:)
      double complex pot,fld(3),hess(6)
      double complex mpole(0:nterms,-nterms:nterms)
c
      ier=0
c
c     Carve up workspace:
c
        allocate( scarray(10*(nterms+2)**2) )
        call l3dtaevalhessdini(nterms,scarray)
c
        call l3dtaevalhessd(rscale,center,mpole,nterms,ztarg,
     1            pot,iffld,fld,ifhess,hess,scarray)
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dallhess(iffld,ifhess,sources,charge,ns,
     1                   target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD,
c     and Hessian HESS at the target point TARGET, due to a collection 
c     of charges at SOURCE(3,ns). 
c     
c              	pot =  sum 1/r
c		fld =  -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> dont compute 
c		                   ifhess = 1 -> do compute 
c     sources(3,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (double complex)        : calculated potential
c     fld   (double complex)        : calculated gradient
c     hess  (double complex)        : calculated Hessian
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,ns,i,j
      double precision sources(3,ns),target(3)
      double complex pot,fld(3),hess(6),potloc,fldloc(3),hessloc(6)
      double complex eye
      double complex charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if (ifhess.eq.1) then
         do i = 1,6
            hess(i) = 0.0d0
         enddo
      endif
c
      do i = 1,ns
         call lpotfld3dhess(iffld,ifhess,sources(1,i),charge(i),
     1        target,potloc,fldloc,hessloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
         if (ifhess.eq.1) then
            do j = 1,6
               hess(j) = hess(j) + hessloc(j)
            enddo
         endif
      enddo
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotfld3dhess(iffld,ifhess,source,charge,target,
     1                         pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD
c     and Hesian HESS at the target point TARGET, due to a charge at 
c     SOURCE. 
c     
c              	pot  = 1/r
c		fld  = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     iffld     : flag for computing gradient
c	                 	iffld = 0 -> dont compute 
c		                iffld = 1 -> do compute 
c     ifhess    : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     fld       : calculated gradient
c     hess      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit none
      integer iffld,ifhess
      double precision source(3),target(3)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,ddinv,dddinv,dddddinv
      double complex pot,fld(3),hess(6)
      double complex h0,h1,cd,eye,z,ewavek
      double complex charge
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv=1.0d0/d
c
c       O(1/r)
c
      pot=charge*dinv
c
      ddinv=dinv*dinv
      dddinv=dinv*ddinv
      dddddinv=ddinv*dddinv
c
      if (iffld.eq.1) then
c
c       O(1/r^2)
c
         cd=charge*dddinv
         fld(1)=cd*xdiff
         fld(2)=cd*ydiff
         fld(3)=cd*zdiff
      endif
c
      if (ifhess.eq.1) then
c
c       O(1/r^3)
c
         cd=charge*dddddinv
         hess(1)=cd*(3*xdiff*xdiff-dd)
         hess(2)=cd*(3*ydiff*ydiff-dd)
         hess(3)=cd*(3*zdiff*zdiff-dd)
         hess(4)=cd*3*xdiff*ydiff
         hess(5)=cd*3*xdiff*zdiff
         hess(6)=cd*3*ydiff*zdiff
      endif
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dallhess_dp(iffld,ifhess,
     $     sources,dipstr,dipvec,ns,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT field FLD and
c     Hessian HESS at the target point TARGET, due to a collection 
c     of dipoles at SOURCE(3,ns). 
c     
c              	pot = (dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     ifhess       : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     dipstr(ns)    : dipole strength
c     dipvec(3,ns)  : dipole direction
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,ns,i,j
      double precision sources(3,ns),target(3)
      double precision dipvec(3,ns)
      double complex pot,fld(3),hess(6),potloc,fldloc(3),hessloc(6)
      double complex eye
      double complex dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if (ifhess.eq.1) then
         do i = 1,6
            hess(i) = 0.0d0
         enddo
      endif
c
      do i = 1,ns
         call lpotfld3dhess_dp(iffld,ifhess,
     $        sources(1,i),dipstr(i),dipvec(1,i),
     1        target,potloc,fldloc,hessloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
         if (ifhess.eq.1) then
            do j = 1,6
               hess(j) = hess(j) + hessloc(j)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dhess_dp(
     $     iffld,ifhess,source,dipstr,dipvec,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT field FLD and
c     Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = (dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3 
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     ifhess       : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
c     source(3)    : location of the source 
c     dipstr(ns)   : dipole strength
c     dipvec(3,ns) : dipole direction
c     target(3)    : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot          : calculated potential
c     fld          : calculated -gradient
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess
      double precision source(3),target(3)
      double precision dipvec(3)
      double precision xdiff,ydiff,zdiff,dd,d,dinv,ddinv,dddinv,dotprod
      double precision rx,ry,rz,dx,dy,dz,rtmp,dddddinv
      double complex pot,fld(3),hess(6)
      double complex cd,eye,ztttt,cd2
      double complex dipstr,z1,z2,z3
c
      data eye/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=sqrt(dd)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      dinv = 1.0d0/d
      ddinv=dinv*dinv  
      dddinv=ddinv*dinv  
c
c       O(1/r^2)
c
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)+zdiff*dipvec(3)
      cd = dipstr*dotprod
      pot=cd*dddinv
c
      dddddinv=ddinv*dddinv
c
      if (iffld.eq.1) then
c
c       O(1/r^3)
c
         ztttt = 3.0d0*cd*dddddinv
         cd2 = dipstr*dddinv
         fld(1)=ztttt*xdiff-cd2*dipvec(1)
         fld(2)=ztttt*ydiff-cd2*dipvec(2)
         fld(3)=ztttt*zdiff-cd2*dipvec(3)
      endif 
c
      if (ifhess.eq.1) then
c
c       O(1/r^4)
c
         rx=xdiff
         ry=ydiff
         rz=zdiff
c
         dx=rx*dinv
         dy=ry*dinv
         dz=rz*dinv
c
         rtmp=dotprod
c
         hess(1)=3*(rtmp*(5*dx*dx-1)-(dipvec(1)*rx+dipvec(1)*rx))
         hess(2)=3*(rtmp*(5*dy*dy-1)-(dipvec(2)*ry+dipvec(2)*ry))
         hess(3)=3*(rtmp*(5*dz*dz-1)-(dipvec(3)*rz+dipvec(3)*rz))
c
         hess(4)=3*(rtmp*(5*dx*dy)-(dipvec(2)*rx+dipvec(1)*ry))
         hess(5)=3*(rtmp*(5*dx*dz)-(dipvec(3)*rx+dipvec(1)*rz))
         hess(6)=3*(rtmp*(5*dy*dz)-(dipvec(3)*ry+dipvec(2)*rz))
c
         cd=dipstr*dddddinv
         hess(1)=hess(1)*cd
         hess(2)=hess(2)*cd
         hess(3)=hess(3)*cd
         hess(4)=hess(4)*cd
         hess(5)=hess(5)*cd
         hess(6)=hess(6)*cd
c
      endif 
c
      return
      end
c
c
c
c      The next set of routines are (un)accelerated versions of the 
c      routines for handling expansions due to quadrupole
c      sources.
c
c      Remarks on scaling conventions.
c
c      1)  Far field and local expansions are consistently rscaled as
c              
c
c          M_n^m (scaled) = M_n^m / rscale^(n)  so that upon evaluation
c
c          the field is  sum   M_n^m (scaled) * rscale^(n) / r^{n+1}.
c
c          L_n^m (scaled) = L_n^m * rscale^(n)  so that upon evaluation
c
c          the field is  sum   L_n^m (scaled) / rscale^(n) * r^{n}.
c
c
c      2) There are many definitions of the spherical harmonics,
c         which differ in terms of normalization constants. We
c         adopt the following convention:
c
c         For m>0, we define Y_n^m according to 
c
c         Y_n^m = \sqrt{2n+1} \sqrt{\frac{ (n-m)!}{(n+m)!}} \cdot
c                 P_n^m(\cos \theta)  e^{i m phi} 
c         and
c 
c         Y_n^-m = dconjg( Y_n^m )
c    
c         We omit the Condon-Shortley phase factor (-1)^m in the 
c         definition of Y_n^m for m<0. (This is standard in several
c         communities.)
c
c         We also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the Y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  (This is also standard in several communities.)
c         More precisely, 
c
c                 \int_S Y_n^m Y_n^m d\Omega = 4 \pi. 
c
c         Using our standard definition, the addition theorem takes 
c         the simple form 
c
c         1/r = 
c         \sum_n 1/(2n+1) \sum_m  |S|^n Ylm*(S) Ylm(T)/ (|T|^(n+1)) 
c
c
c-----------------------------------------------------------------------
c
c      f90 version, using allocate
c
c      L3DFORMMP_QUAD: creates multipole expansion (outgoing) due to 
c                 a collection of quadrupoles.
c                 (calls L3DFORMMP1/L3DFORMMP0_QUAD )
c
c      LPOTFLD3DALL_QUAD: 
c                 direct calculation for a collection of quadrupoles
c      LPOTFLD3D_QUAD : direct calculation for a single quadrupole
c
c
c      L3DFORMTA_QUAD: 
c                 (calls L3DFORMTA1/L3DFORMTA0_QUAD )
c
c
C***********************************************************************
      subroutine l3dformmp_quad(ier,rscale,sources,quadvec,ns,
     1                  center,nterms,mpole)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, istotal
      double precision center(3),sources(3,ns)
      double precision quadvec(6,ns)
      double precision rscale
      double precision, allocatable :: sss(:)
      double complex mpole(0:nterms,-nterms:nterms)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      istotal = 5*(2*nterms+1)
      allocate(sss(istotal))
c
      call getsgnformpmp_quad(sss,nterms)
      do i = 1, ns
         call l3dformmp1_quad(ier,rscale,sources(1,i),
     1        quadvec(1,i),center,nterms,mpole,sss)
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_quad(ier,rscale,source,quadvec,
     1		center,nterms,mpole,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a quadrupole located at the point SOURCE.
c     This is a memory management routine. Work is done in the
c     secondary call to l3dformmp0_quad below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     quadvec  : quadrupole vector
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c     sss     : sign mutiplier array (see getsgnformpmp_quad)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double precision quadvec(6)
      double precision sss(-2:2,-nterms:nterms)
      double complex mpole(0:nterms,-nterms:nterms)
c
c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      lused=ifr + lfr
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call l3dformmp0_quad(jer,rscale,source,quadvec,
     1		center,nterms,mpole,w(ipp),w(iephi),
     2          w(ifr),sss)
      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_quad(ier,rscale,source,quadvec,
     1		center,nterms,mpole,pp,ephi,fr,sss)
c**********************************************************************
c
c     See l3dformmp1_quad for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew,mnew2
      double precision binom(0:120,0:4)
      double precision rscale,source(3),center(3),zdiff(3)
      double precision quadvec(6)
      double precision pp(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision cscale,cscale2,cscale3
      double precision sss(-2:2,-nterms:nterms)
      double precision a22,a21,a20,rmul,dd
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double complex mp2(-2:2)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1),ephi1
      double complex eye,zmul
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*(i))
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dcmplx(cphi,-sphi)
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
C
      call ylgndr(nterms,ctheta,pp)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
         do m = -2,2
            do ll = 0,nterms-2
               cscale = rscale*rscale*fr(LL)/dsqrt(2*ll+1.0d0)
               lnew = 2+ll
c
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1         pp(ll,0)*cscale2*mp2(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew,2-m)*binom(lnew+mnew,2+m)
                  cscale2 = cscale2*cscale*sss(m,mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
c
                  mpole(lnew,mnew) = mpole(lnew,mnew)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp2(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp2(m)
               enddo
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
      subroutine l3dformmp_quad_trunc(ier,rscale,sources,quadvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, isss, lsss
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      double precision center(3),sources(3,ns)
      double precision quadvec(6,ns)
      double precision rscale
      double complex mpole(0:nterms,-nterms:nterms)
      double precision binom(0:120,0:4)
      integer nlege
      double precision wlege(*)
      double precision, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      isss=ifr + lfr
      lsss = 5*(2*nterms+1)
c
      lused=isss + lsss
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

ccc      ... used in the reference code only
ccc      call getsgnformpmp_quad(w(isss),nterms)

      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo


      do i = 1, ns

c        call l3dformmp1_quad(ier,rscale,sources(1,i),
c     1        quadvec(1,i),center,nterms,mpole,w(isss))

      call l3dformmp0_quad_trunc(ier,rscale,sources(1,i),quadvec(1,i),
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,w(isss))

      enddo
c
      return
      end
C
c**********************************************************************
      subroutine l3dformmp1_quad_trunc(ier,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,binom,sss)
c**********************************************************************
c
c     This subroutine creates the multipole expansion about CENTER
c     due to a quadrupole located at the point SOURCE.
c     This is a memory management routine. Work is done in the
c     secondary call to l3dformmp0_quad below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     quadvec  : quadrupole vector
c     center  : coordinates of the expansion center
c     nterms  : order of the expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c     sss     : sign mutiplier array (see getsgnformpmp)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      deprecated but left in calling sequence for
c		      backward compatibility.
c     mpole   : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double precision quadvec(6)
      double precision binom(0:120,0:4)
      double precision sss(-2:2,-nterms:nterms)
      double complex mpole(0:nterms,-nterms:nterms)
      integer nlege
      double precision wlege(*)
c
c
c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      lused=ifr + lfr
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

      call l3dformmp0_quad_trunc(jer,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,sss)

      return
      end
c
c
c
c**********************************************************************
      subroutine l3dformmp0_quad_trunc(ier,rscale,source,quadvec,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr,binom,sss)
c**********************************************************************
c
c     See l3dformmp1_quad for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      double precision binom(0:120,0:4)
      double precision rscale,source(3),center(3),zdiff(3)
      double precision quadvec(6)
      double precision pp(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision cscale,cscale1,cscale2,cscale3,rtmp
      double precision sss(-2:2,-nterms:nterms)
      double precision a22,a21,a20,rmul,dd
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double complex mp2(-2:2)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1),ephi1
      double complex eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      double precision wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi(1))
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
c      do l = 0,nterms
c      fr(l)=fr(l)/sqrt(2*l+1.0d0) *rscale*rscale
c      enddo
c      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
      do l = 0,nterms
      fr(l)=fr(l) *rscale*rscale
      enddo
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
        if( 1 .eq. 2 ) then
c
c       ... reference code
c
         do m = -2,2
            do ll = 0,nterms-2
               cscale = fr(ll)
               lnew = 2+ll
c
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            pp(ll,0)*cscale2*mp2(m)*cscale*sss(m,0)
               do mm = 1,ll
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,2-m)*binom(lnew+mnew1,2+m)
                  cscale2 = cscale2*cscale*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp2(m)
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp2(m)

               enddo
            enddo
         enddo
         endif
c
c
        if( 2 .eq. 2 ) then
c
c       ... optimized for speed
c       use symmetries for real valued quadrupoles
c
            do ll = 0,nterms-2

               lnew = 2+ll
c
               cscale2 = binom(lnew,2)*binom(lnew,2)
               mpole(lnew,0) = mpole(lnew,0) +
     1            fr(ll)*pp(ll,0)*cscale2*mp2(0)

               rtmp = fr(ll)*pp(ll,0)
               do m = 1,2
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               ztmp = rtmp*cscale2*mp2(m)
               mpole(lnew,+m) = mpole(lnew,+m)+ztmp
               mpole(lnew,-m) = mpole(lnew,-m)+dconjg(ztmp)
               enddo

               do mm = 1,ll

               ctmp=fr(ll)*pp(ll,mm)*ephi(-mm)

               mnew1 = +mm
               mnew2 = -mm
               cscale2 = binom(lnew-mnew1,2)*binom(lnew+mnew1,2)
               ztmp = cscale2*ctmp*mp2(0)
               mpole(lnew,mnew1) = mpole(lnew,mnew1)+ztmp
               mpole(lnew,mnew2) = mpole(lnew,mnew2)+dconjg(ztmp)

               do m = 1,2
c
c       ... skip zero valued modes
c       
                  if( abs(mp2(m)) .eq. 0 ) cycle

                  mnew1 = m+mm
                  mnew2 = m-mm

                  if( m .eq. 1 ) then
                  cscale2 = +binom(lnew-mnew1,1)*binom(lnew+mnew1,3)
                  cscale3 = -binom(lnew-mnew2,1)*binom(lnew+mnew2,3)
                  else
                  if( mm .eq. 1 ) then
                  cscale2 = +binom(lnew+mnew1,4)
                  cscale3 = -binom(lnew+mnew2,4)
                  else
                  cscale2 = +binom(lnew+mnew1,4)
                  cscale3 = +binom(lnew+mnew2,4)
                  endif
                  endif

                  ztmp0=cscale2*mp2(m)*ctmp
                  ztmp1=cscale3*mp2(m)*dconjg(ctmp)

                  mpole(lnew,+mnew1) = mpole(lnew,+mnew1)+ztmp0
                  mpole(lnew,+mnew2) = mpole(lnew,+mnew2)+ztmp1

                  mpole(lnew,-mnew2) = mpole(lnew,-mnew2)+dconjg(ztmp1)
                  mpole(lnew,-mnew1) = mpole(lnew,-mnew1)+dconjg(ztmp0)

               enddo
               enddo
            enddo
         endif
c
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dall_quad(iffld,sources,quadvec,ns,
     1                   target,pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of quadrupoles 
c     at SOURCE(3,ns). 
c     
c		fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     quadvec(6,ns)  : dipole direction
c     ns            : number of sources
c     charge(ns)    : charge strength
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c----------------------------------------------------------------------
      implicit none
      integer iffld,i,ns
      double precision sources(3,ns),target(3)
      double precision quadvec(6,ns)
      double complex pot,fld(3),potloc,fldloc(3)
c
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotfld3d_quad(iffld,sources(1,i),quadvec(1,i),
     1        target,potloc,fldloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3d_quad(iffld,source,quadvec,target,
     1                        pot,fld)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a quadrupole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c    
c               pot = quadvec(1)*V_xx +
c                     quadvec(2)*V_yy +
c                     quadvec(3)*V_zz +
c                     quadvec(4)*V_xy +
c                     quadvec(5)*V_xz +
c                     quadvec(6)*V_yz 
c
c      V_xx = (-1/r^3 + 3*dx**2/r^5)
c      V_xy = 3*dx*dy/r^5
c      V_xz = 3*dx*dz/r^5
c      V_yy = (-1/r^3 + 3*dy**2/r^5)
c      V_yz = 3*dy*dz/r^5
c      V_zz = (-1/r^3 + 3*dz**2/r^5)
c
c		fld = -grad(pot)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     source(3)    : location of the source 
c     quadvec(3,ns) : quadrupole vector
c     target(3)    : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot          : calculated potential
c     fld          : calculated -gradient
c
c----------------------------------------------------------------------
      implicit none
      integer iffld
      double precision source(3),target(3)
      double precision quadvec(6),rr(3)
      double precision cd,d5,d3,v11,v12,v13,v22,v23,v33
      double precision d7,v111,v112,v113,v122,v123,v133
      double precision v222,v223,v233,v333
      double complex pot,fld(3)
c
c
c ... Caculate offsets and distance
c
      rr(1)=target(1)-source(1)
      rr(2)=target(2)-source(2)
      rr(3)=target(3)-source(3)
      cd = sqrt(rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3))
      d5 = 1.0d0/(cd**5)
      d3 = 1.0d0/(cd**3)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      v11 = -d3 + 3*d5*rr(1)**2
      v12 = 3*d5*rr(1)*rr(2)
      v13 = 3*d5*rr(1)*rr(3)
      v22 = -d3 + 3*d5*rr(2)**2
      v23 = 3*d5*rr(2)*rr(3)
      v33 = -d3 + 3*d5*rr(3)**2
c
      pot=v11*quadvec(1)+v22*quadvec(2)+v33*quadvec(3)
      pot=pot+v12*quadvec(4)+v13*quadvec(5)+v23*quadvec(6)
ccc      call prin2(' cd is *',cd,1)
ccc      call prin2(' d3 is *',d3,1)
ccc      call prin2(' v11 is *',v11,1)
ccc      call prin2(' v22 is *',v22,1)
ccc      call prin2(' v33 is *',v33,1)
ccc      call prin2(' v12 is *',v12,1)
ccc      call prin2(' v13 is *',v13,1)
ccc      call prin2(' v23 is *',v23,1)
ccc      call prin2(' pot is *',pot,1)
      if (iffld.eq.1) then
         d7 = 1.0d0/(cd**7)
         v111 = 9*rr(1)*d5 - 15*d7*rr(1)**3
         v112 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(1)**2
         v113 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(1)**2
         v122 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(2)**2
         v123 = -15*d7*rr(1)*rr(2)*rr(3)
         v133 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(3)**2
         v222 = 9*rr(2)*d5 - 15*d7*rr(2)**3
         v223 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(2)**2
         v233 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(3)**2
         v333 = 9*rr(3)*d5 - 15*d7*rr(3)**3
ccc         call prin2(' v111 is *',v111,1)
ccc         call prin2(' v112 is *',v112,1)
ccc         call prin2(' v113 is *',v113,1)
ccc         call prin2(' v122 is *',v122,1)
ccc         call prin2(' v123 is *',v123,1)
ccc         call prin2(' v133 is *',v133,1)
ccc         call prin2(' v222 is *',v222,1)
ccc         call prin2(' v223 is *',v223,1)
ccc         call prin2(' v233 is *',v233,1)
ccc         call prin2(' v333 is *',v333,1)
         fld(1)=v111*quadvec(1)+v122*quadvec(2)+v133*quadvec(3)
         fld(1)=fld(1)+v112*quadvec(4)+v113*quadvec(5)+v123*quadvec(6)
         fld(2)=v112*quadvec(1)+v222*quadvec(2)+v233*quadvec(3)
         fld(2)=fld(2)+v122*quadvec(4)+v123*quadvec(5)+v223*quadvec(6)
         fld(3)=v113*quadvec(1)+v223*quadvec(2)+v333*quadvec(3)
         fld(3)=fld(3)+v123*quadvec(4)+v133*quadvec(5)+v233*quadvec(6)
         fld(1) = -fld(1)
         fld(2) = -fld(2)
         fld(3) = -fld(3)
      endif 
      return
      end
cc
c
c
c**********************************************************************
      subroutine l3dformta_quad(ier,rscale,sources,quadvec,ns,
     1		           center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS quadrupoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1_quad/l3dformta0_quad below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     quadvec    : quadrupole vector
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c
c     locexp    : coeffs for the expansion
c     lused     : amount of work space "w" used
c
c
c----------------------------------------------------------------------
      implicit none
      integer ier,ns,nterms,i,l,m
      double precision rscale,sources(3,ns),center(3),rs
      double precision quadvec(6,ns)
      double complex locexp(0:nterms,-nterms:nterms)
      double complex eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call l3dformta1_quad(ier,rscale,sources(1,i),
     1		quadvec(1,i),center,nterms,locexp)
      enddo
c
      do l = 0,nterms
         locexp(l,0) = locexp(l,0)*rscale
         do m=1,l
            locexp(l,m) = locexp(l,m)*rscale
            locexp(l,-m) = dconjg(locexp(l,m))
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta1_quad(ier,rscale,source,quadvec,
     &		center,nterms,locexp)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single quadrupole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     quadvec   : quadrupole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   ier=2 insufficient memory in workspace w
c     locexp : coefficients of the local expansion
c     lused  : amount of work space "w" used
c
c
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms,ipp,lpp,iephi,lephi,ifr,lfr,lused
      double precision rscale,source(3),center(3)
      double precision, allocatable :: w(:)
      double precision quadvec(6)
      double complex locexp(0:nterms,-nterms:nterms)
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+3)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+3)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+5)
c
      lused=ifr+lfr
      allocate(w(lused))
c
      call l3dformta0_quad(rscale,source,quadvec,
     &   center,nterms,locexp,w(ipp),w(iephi),w(ifr))
c
      return
      end
c
c
c
c
      subroutine l3dformta0_quad(rscale,source,quadvec,
     &		center,nterms,locexp,pp,ephi,fr)
c**********************************************************************
c
c     See l3dformta_quad/l3dformta1_quad for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,i,nt2
      integer l,m,lnew,mnew,ll,mm
      double precision binom(0:120,0:4)
      double precision rscale,source(3),center(3),zdiff(3)
      double precision quadvec(6)
      double precision pp(0:nterms+2,0:nterms+2)
      double precision fr(0:nterms+3)
      double precision cscale,cscale2,cscale3
ccc      double precision sss(-2:2,-nterms:nterms)
      double precision a22,a21,a20,rmul,dd
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double precision c(0:200,0:200),sqc(0:200,0:200)
      double complex mp2(-2:2),mpolelm
      double complex locexp(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-2:nterms+2),ephi1
      double complex eye,zmul
      data eye/(0.0d0,1.0d0)/
c

      nt2 = nterms+2
c
      do l=0,2*nt2
         c(l,0)=1.0d0
         sqc(l,0)=1.0d0
      enddo
      do m=1,2*nt2
         c(m,m)=1.0d0
         sqc(m,m)=1.0d0
         do l=m+1,2*nt2
            c(l,m)=c(l-1,m)+c(l-1,m-1)
            sqc(l,m)=dsqrt(c(l,m))
         enddo
      enddo
C
c
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
c     using truncated version of full multipole-local shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      zdiff(1)=-zdiff(1)
      zdiff(2)=-zdiff(2)
      zdiff(3)=-zdiff(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      fr(0) = 1.0D0
      dd = 1.0d0/r/rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+3
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
C
      call ylgndr(nt2,ctheta,pp)
c
C
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
      do m = -2,2
	 mpolelm = mp2(m)
         do lnew = 0,nterms
            do mnew = 0,lnew
               ll = 2 + lnew
               mm = mnew - m
               cscale = rscale*rscale*fr(ll+1)*sqc(ll+mm,lnew+mnew)*
     1                     sqc(ll-mm,lnew-mnew)
               cscale = cscale*(-1)**lnew
               cscale = cscale/dsqrt(2*ll+1.0d0)
               if ( (m .lt. 0) .and. (mnew .lt. 0) ) then
                  if (-mnew .lt. -m)  cscale = cscale*(-1)**mnew
                  if (-mnew .ge. -m)  cscale = cscale*(-1)**m
               endif
               if ( (m .gt. 0) .and. (mnew .gt. 0) ) then
                  if (mnew .lt. m) cscale = cscale*(-1)**mnew
                  if (mnew .ge. m) cscale = cscale*(-1)**m
               endif
               if (mm .eq. 0) then
                  locexp(lnew,mnew) = locexp(lnew,mnew) +
     1            pp(ll,0)*cscale*mpolelm
               else  if (mm .gt. 0) then
                  locexp(lnew,mnew) = locexp(lnew,mnew)+ cscale*
     1            pp(ll,mm)*ephi(-mm)*mpolelm
               else
                  locexp(lnew,mnew) = locexp(lnew,mnew)+ cscale*
     1            pp(ll,-mm)*ephi(-mm)*mpolelm
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine getsgnformpmp_quad(sss,nterms)
C***********************************************************************
c
c     This subroutine creates a multiplier that holds the 
c     appropriate power of (-1) that arises in multipole-multipole shifts.
c     (See Greengard - 
c          Rapid Evaluation of Potential Fields... for notation).
c-------------------------------------------------------------------------
      implicit none
      integer m,mm,nterms
      double precision sss(-2:2,-nterms:nterms)
c
      do m = -2,2
      do mm = -nterms,nterms
         sss(m,mm) = 1.0d0
      enddo
      enddo
      do m = -2,2
      do mm = -nterms,nterms
         if (( m .lt. 0) .and. (mm .gt. 0)) then
            if ( mm .le. -m ) sss(m,mm) = (-1)**mm
            if ( mm .gt. -m ) sss(m,mm) = (-1)**m
         endif
         if (( m .gt. 0) .and. (mm .lt. 0)) then
            if ( m .le. -mm ) sss(m,mm) = (-1)**m
            if ( m .gt. -mm ) sss(m,mm) = (-1)**mm
         endif
      enddo
      enddo
      return
      end
c
c
c
c    l3dmpevalhessd  uses direct translation (not rotation/zshift)
c    and is reasonably optimized, precomputing the array of 
c    binomial/factorial terms that appear in the shift operator.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine l3dmpevalhessd(rscale,center,mpole,nterms,targ,
     1            pot,iffld,fld,ifhess,hess,scarray)
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the potential due to a multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     targ   :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DMPEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c
c     OUTPUT:
c
c     pot    :    potential
c     fld    :    if (iffld .eq.1)   
c     hess   :    if (ifhess .eq.1)   
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,lsum
      double precision center(3),targ(3)
      double precision zdiff(3)
      double precision scarray(*),rscale
      double precision cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      double complex mpole(0:nterms,-nterms:nterms)
      double complex local2(0:2,-2:2)
      double complex z0,ima,ephi1,pot,fld(3),hess(6)
c
      double precision, allocatable :: pp(:,:)
      double precision, allocatable :: powers(:)
      double complex, allocatable :: ppc(:,:)
      double complex, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms+2,0:nterms+2))
      allocate(ppc(0:nterms+2,-nterms-2:nterms+2))
      allocate(powers(0:nterms+3))
      allocate(ephi(-nterms-3:nterms+3))
c
c     determine order of shifted expansion 
c
      ll = 0
      if (iffld.eq.1) ll = 1
      if (ifhess.eq.1) ll = 2
c
      do l = 0,ll
      do m = -l,l
         local2(l,m) = 0.0d0
      enddo
      enddo
c
      zdiff(1) = center(1) - targ(1)
      zdiff(2) = center(2) - targ(2)
      zdiff(3) = center(3) - targ(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
      dd = 1.0d0/d
      dd = dd/rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+3
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
      enddo
c
      call ylgndr(nterms+2,ctheta,pp)
      do l = 0,nterms+2
         do k = -l,l
            ppc(l,k) = pp(l,abs(k))*powers(l+1)*ephi(-k)
         enddo
      enddo
c
c     shift to local expansion of order ll about target point
c
      iuse = 1
      do l = 0,nterms
         do m = -l,l
            local2(0,0) = local2(0,0) +
     1      ppc(l,-m)*scarray(iuse)*mpole(l,m)
            iuse = iuse+1
         enddo
      enddo
      if (ll.ge.1) then
         do l = 0,nterms
            lsum = l+1
            do m = -l,l
               local2(1,-1) = local2(1,-1) +
     1         ppc(lsum,-1-m)*scarray(iuse  )*mpole(l,m)
               local2(1, 0) = local2(1, 0) +
     1         ppc(lsum,  -m)*scarray(iuse+1)*mpole(l,m)
               local2(1, 1) = local2(1, 1) +
     1         ppc(lsum, 1-m)*scarray(iuse+2)*mpole(l,m)
               iuse = iuse+3
            enddo
         enddo
      endif
      if (ll.eq.2) then
         do l = 0,nterms
            lsum = l+2
            do m = -l,l
               local2(2,-2) = local2(2,-2) +
     1         ppc(lsum,-2-m)*scarray(iuse  )*mpole(l,m)
               local2(2,-1) = local2(2,-1) +
     1         ppc(lsum,-1-m)*scarray(iuse+1)*mpole(l,m)
               local2(2, 0) = local2(2, 0) +
     1         ppc(lsum,  -m)*scarray(iuse+2)*mpole(l,m)
               local2(2, 1) = local2(2, 1) +
     1         ppc(lsum, 1-m)*scarray(iuse+3)*mpole(l,m)
               local2(2, 2) = local2(2, 2) +
     1         ppc(lsum, 2-m)*scarray(iuse+4)*mpole(l,m)
               iuse = iuse+5
            enddo
         enddo
      endif
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local2(0,0)*rscale
c
c     fld comes from l=1 modes
c
      if (iffld.eq.1) then
         rfac = sqrt(2.0d0)*rscale*rscale
         fld(1) = rfac*(local2(1,1) + local2(1,-1))/2.0d0
         fld(2) = rfac*ima*(local2(1,1) - local2(1,-1))/2.0d0
         fld(3) = -local2(1,0)*rscale*rscale
      endif
c
c     hess comes from l=2 modes
c
      if (ifhess.eq.1) then
         rfac = rscale*rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
         z0 = rscale*rscale*rscale*local2(2,0)
         hess(1) = rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(2) = -rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(3) = 2*z0
         hess(4) = rfac*ima*(local2(2,2) - local2(2,-2))
         hess(5) = -rfac*(local2(2,1) + local2(2,-1))
         hess(6) = -rfac*ima*(local2(2,1) - local2(2,-1))
      endif
c
      return
      end
c
c
c
c
c
      subroutine l3dmpevalhessdini(nterms,scarray)
      implicit none
      integer  nterms,l,j,k,m,ll,mm,iuse,lnew,mnew
      double precision scarray(1),cscale
      double precision d
      double precision, allocatable :: c(:,:)
      double precision, allocatable :: sqc(:,:)
c
c     This subroutine is used to precompute various terms that appear in 
c     the local-local translation operator from an nterm expansion to an 
c     order 2 expansion (sufficient to compute pot/fld/hessian).     
c
c     INPUT: nterms
c     OUTPUT: scarray array  MUST BE DIMENSIONED 
c                            at least 10*(nterms+2)**2
c
      allocate(c(0:2*nterms+4,0:2*nterms+4))
      allocate(sqc(0:2*nterms+4,0:2*nterms+4))
c
      do l = 0,2*nterms+4
         c(l,0) = 1.0d0
         sqc(l,0) = 1.0d0
      enddo
      do m = 1,2*nterms+4
         c(m,m) = 1.0d0
         sqc(m,m) = 1.0d0
         do l = m+1,2*nterms+4
            c(l,m) = c(l-1,m)+c(l-1,m-1)
            sqc(l,m) = dsqrt(c(l,m))
         enddo
      enddo
c
      iuse = 1
      do lnew= 0,2
         do l = 0,nterms
            do m = -l,l
               do mnew = -lnew,lnew
                  ll = l+lnew
                  mm = mnew-m
                  cscale = sqc(ll+mm,lnew+mnew)*sqc(ll-mm,lnew-mnew)
                  cscale = cscale*(-1)**l
                  cscale = cscale/dsqrt(2*ll+1.0d0)
                  if ( (m .lt. 0) .and. (mnew .lt. 0) ) then
                     if (-mnew .lt. -m)  cscale = cscale*(-1)**mnew
                     if (-mnew .ge. -m)  cscale = cscale*(-1)**m
                  endif
                  if ( (m .gt. 0) .and. (mnew .gt. 0) ) then
                     if (mnew .lt. m) cscale = cscale*(-1)**mnew
                     if (mnew .ge. m) cscale = cscale*(-1)**m
                  endif
                  scarray(iuse) = cscale
                  iuse = iuse+1
               enddo
            enddo
         enddo
      enddo
      return
      end
c
c
c    l3dmpevalhessd  uses direct translation (not rotation/zshift)
c    and is reasonably optimized, precomputing the array of 
c    binomial/factorial terms that appear in the shift operator.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine l3dmpevalhessd_trunc(rscale,center,mpole,nterms,targ,
     1            pot,iffld,fld,ifhess,hess,scarray,wlege,nlege)
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the potential due to a multipole expansion.
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     targ   :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DMPEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c
c     OUTPUT:
c
c     pot    :    potential
c     fld    :    if (iffld .eq.1)   
c     hess   :    if (ifhess .eq.1)   
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess,nlege
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,lsum
      double precision center(3),targ(3),wlege(*)
      double precision zdiff(3)
      double precision scarray(*),rscale
      double precision cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      double complex mpole(0:nterms,-nterms:nterms)
      double complex local2(0:2,-2:2)
      double complex z0,ima,ephi1,pot,fld(3),hess(6)
c
      double precision, allocatable :: pp(:,:)
      double precision, allocatable :: powers(:)
      double complex, allocatable :: ppc(:,:)
      double complex, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms+2,0:nterms+2))
      allocate(ppc(0:nterms+2,-nterms-2:nterms+2))
      allocate(powers(0:nterms+3))
      allocate(ephi(-nterms-3:nterms+3))
c
c     determine order of shifted expansion 
c
      ll = 0
      if (iffld.eq.1) ll = 1
      if (ifhess.eq.1) ll = 2
c
      do l = 0,ll
      do m = -l,l
         local2(l,m) = 0.0d0
      enddo
      enddo
c
      zdiff(1) = center(1) - targ(1)
      zdiff(2) = center(2) - targ(2)
      zdiff(3) = center(3) - targ(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
      dd = 1.0d0/d
      dd = dd/rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+3
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
      enddo
c
      call ylgndrfw(nterms+2,ctheta,pp,wlege,nlege)
      do l = 0,nterms+2
         do k = -l,l
            ppc(l,k) = pp(l,abs(k))*powers(l+1)*ephi(-k)
         enddo
      enddo
c
c     shift to local expansion of order ll about target point
c
      iuse = 1
      do l = 0,nterms
         do m = -l,l
            local2(0,0) = local2(0,0) +
     1      ppc(l,-m)*scarray(iuse)*mpole(l,m)
            iuse = iuse+1
         enddo
      enddo
      if (ll.ge.1) then
         do l = 0,nterms
            lsum = l+1
            do m = -l,l
               local2(1,-1) = local2(1,-1) +
     1         ppc(lsum,-1-m)*scarray(iuse  )*mpole(l,m)
               local2(1, 0) = local2(1, 0) +
     1         ppc(lsum,  -m)*scarray(iuse+1)*mpole(l,m)
               local2(1, 1) = local2(1, 1) +
     1         ppc(lsum, 1-m)*scarray(iuse+2)*mpole(l,m)
               iuse = iuse+3
            enddo
         enddo
      endif
      if (ll.eq.2) then
         do l = 0,nterms
            lsum = l+2
            do m = -l,l
               local2(2,-2) = local2(2,-2) +
     1         ppc(lsum,-2-m)*scarray(iuse  )*mpole(l,m)
               local2(2,-1) = local2(2,-1) +
     1         ppc(lsum,-1-m)*scarray(iuse+1)*mpole(l,m)
               local2(2, 0) = local2(2, 0) +
     1         ppc(lsum,  -m)*scarray(iuse+2)*mpole(l,m)
               local2(2, 1) = local2(2, 1) +
     1         ppc(lsum, 1-m)*scarray(iuse+3)*mpole(l,m)
               local2(2, 2) = local2(2, 2) +
     1         ppc(lsum, 2-m)*scarray(iuse+4)*mpole(l,m)
               iuse = iuse+5
            enddo
         enddo
      endif
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local2(0,0)*rscale
c
c     fld comes from l=1 modes
c
      if (iffld.eq.1) then
         rfac = sqrt(2.0d0)*rscale*rscale
         fld(1) = rfac*(local2(1,1) + local2(1,-1))/2.0d0
         fld(2) = rfac*ima*(local2(1,1) - local2(1,-1))/2.0d0
         fld(3) = -local2(1,0)*rscale*rscale
      endif
c
c     hess comes from l=2 modes
c
      if (ifhess.eq.1) then
         rfac = rscale*rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
         z0 = rscale*rscale*rscale*local2(2,0)
         hess(1) = rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(2) = -rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(3) = 2*z0
         hess(4) = rfac*ima*(local2(2,2) - local2(2,-2))
         hess(5) = -rfac*(local2(2,1) + local2(2,-1))
         hess(6) = -rfac*ima*(local2(2,1) - local2(2,-1))
      endif
c
      return
      end
c
c
c
c    l3dtaevalhessd  uses direct translation (not rotation/zshift)
c    and is reasonably optimized, precomputing the array of 
c    binomial/factorial terms that appear in the shift operator.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine l3dtaevalhessd(rscale,center,local,nterms,targ,
     1            pot,iffld,fld,ifhess,hess,scarray)
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the potential due to a local harmonic expansion.
c
c     pot =  sum sum  local(n,m) Y_nm(theta,phi)  * r^n
c             n   m
c
c     fld  = -gradient(pot) if iffld = 0.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 0.
c
c     where rscale defines scaling parameter.     
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     local  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     targ   :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DTAEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c   
c     OUTPUT:
c
c     pot    :    potential
c     fld    :    if (iffld .eq.1)   
c     hess   :    if (ifhess .eq.1)   
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c
c--------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,ldiff
      double precision center(3),targ(3)
      double precision zdiff(3)
      double precision scarray(*),rscale
      double precision cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      double complex local(0:nterms,-nterms:nterms)
      double complex local2(0:2,-2:2)
      double complex z0,ima,ephi1,pot,fld(3),hess(6)
c
      double precision, allocatable :: pp(:,:)
      double precision, allocatable :: powers(:)
      double complex, allocatable :: ppc(:,:)
      double complex, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms,0:nterms))
      allocate(ppc(0:nterms,-nterms:nterms))
      allocate(powers(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
c
c     determine order of shifted expansion 
c
      ll = 0
      if (iffld.eq.1) ll = 1
      if (ifhess.eq.1) ll = 2
c
      do l = 0,ll
      do m = -l,l
         local2(l,m) = 0.0d0
      enddo
      enddo
c
      zdiff(1) = center(1) - targ(1)
      zdiff(2) = center(2) - targ(2)
      zdiff(3) = center(3) - targ(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
      dd = d*rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+1
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
      enddo
c
      call ylgndr(nterms,ctheta,pp)
      do l = 0,nterms
         do k = -l,l
            ppc(l,k) = pp(l,abs(k))*powers(l)*ephi(k)
         enddo
      enddo
c
c     shift to local expansion of order ll about target point
c
c     old code (not unrolled).
c
ccc      iuse = 1
ccc      do j = 0,ll
ccc         do k = -j,j
ccc            do l = j,nterms
ccc               do m = -l,l
ccc                  ldiff = l-j
ccc                  mm = m-k
ccc                  if (abs(mm).le.ldiff) then
ccc                     local2(j,k) = local2(j,k) +
ccc     1               ppc(ldiff,mm)*scarray(iuse)*local(l,m)
ccc                  endif
ccc                  iuse = iuse+1
ccc               enddo
ccc            enddo
ccc         enddo
ccc      enddo
c
      iuse = 1
      do l = 0,nterms
         do m = -l,l
            local2(0,0) = local2(0,0) +
     1      ppc(l,m)*scarray(iuse)*local(l,m)
            iuse = iuse+1
         enddo
      enddo
      if (ll.ge.1) then
         do l = 1,nterms
            ldiff = l-1
            do m = -l,l-2
               mm = m+1
               local2(1,-1) = local2(1,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               local2(1,0) = local2(1,0) +
     1         ppc(ldiff,m)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+2,l
               mm = m-1
               local2(1,1) = local2(1,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
      if (ll.eq.2) then
         do l = 2,nterms
            ldiff = l-2
            do m = -l,l-4
               mm = m+2
               local2(2,-2) = local2(2,-2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+1,l-3
               mm = m+1
               local2(2,-1) = local2(2,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               mm = m
               local2(2,0) = local2(2,0) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+3,l-1
               mm = m-1
               local2(2,1) = local2(2,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+4,l
               mm = m-2
               local2(2,2) = local2(2,2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local2(0,0)
c
c     fld comes from l=1 modes
c
      rfac = rscale/sqrt(2.0d0)
      if (iffld.eq.1) then
         fld(3) = -local2(1,0)*rscale
         fld(1) = rfac*(local2(1,1) + local2(1,-1))
         fld(2) = rfac*ima*(local2(1,1) - local2(1,-1))
      endif
c
c     hess comes from l=2 modes
c
      if (ifhess.eq.1) then
ccc         rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
         rfac = rfac*rscale*sqrt(3.0d0)
         z0 = rscale*rscale*local2(2,0)
         hess(1) = rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(2) = -rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(3) = 2*z0
         hess(4) = rfac*ima*(local2(2,2) - local2(2,-2))
         hess(5) = -rfac*(local2(2,1) + local2(2,-1))
         hess(6) = -rfac*ima*(local2(2,1) - local2(2,-1))
      endif
c
      return
      end
c
c
c
c
c
      subroutine l3dtaevalhessdini(nterms,scarray)
      implicit none
      integer  nterms,l,j,k,m,ll,mm,iuse
      double precision scarray(1)
      double precision d
      double precision, allocatable :: cs(:,:)
      double precision, allocatable :: fact(:)
c
c     This subroutine is used to precompute various terms that appear in 
c     the local-local translation operator from an nterm expansion to an 
c     order 2 expansion (sufficient to compute pot/fld/hessian).     
c
c     INPUT: nterms
c     OUTPUT: scarray array  MUST BE DIMENSIONED 
c                            at least 10*(nterms+2)**2
c
      allocate(cs(0:nterms,-nterms:nterms))
      allocate(fact(0:2*nterms))
c
      d = 1.0d0
      fact(0) = d
      do l = 1,2*nterms
         d = d*dsqrt(l+0.0D0)
         fact(l) = d
      enddo
      cs(0,0) = 1.0d0
      do l = 1,nterms
         do m = 0,l
            cs(l,m) =  ((-1)**l)/( fact(l-m)*fact(l+m) )
            cs(l,-m) = cs(l,m)
         enddo
      enddo
c
      iuse = 1
      do j = 0,2
ccc         do k = -j,j
            do l = j,nterms
            do k = -j,j
               do m = -l,l
                  ll = l-j
                  mm = m-k
                  scarray(iuse) = 0.0d0
                  if (abs(mm).gt.ll) goto 111
                  scarray(iuse) = cs(j,k)*cs(ll,mm)/cs(l,m)
                  scarray(iuse) = scarray(iuse)/dsqrt(2*LL+1.0D0)
                  scarray(iuse) = scarray(iuse)*(-1)**ll
                  if ( m*mm .lt. 0) then
                     scarray(iuse) = scarray(iuse)*(-1)**mm
                  endif
                  if (m*mm .ge. 0)  then
                     if ( abs(m) .le. abs(mm) ) 
     1               scarray(iuse) = scarray(iuse)*(-1)**k
                  endif
                  iuse = iuse+1
111            continue
               enddo
            enddo
         enddo
      enddo
      return
      end
c
c
c
c    l3dtaevalhessd  uses direct translation (not rotation/zshift)
c    and is reasonably optimized, precomputing the array of 
c    binomial/factorial terms that appear in the shift operator.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine l3dtaevalhessd_trunc(rscale,center,local,nterms,targ,
     1            pot,iffld,fld,ifhess,hess,scarray,wlege,nlege)
c
c     This subroutine evaluates the potential, -gradient and
c     Hessian of the potential due to a local harmonic expansion.
c
c     pot =  sum sum  local(n,m) Y_nm(theta,phi)  * r^n
c             n   m
c
c     fld  = -gradient(pot) if iffld = 1.
c     hess = dxx,dyy,dzz,dxy,dxz,dyz of pot if ifhess = 1.
c
c     where rscale defines scaling parameter.     
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     local  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     targ   :    target location
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     ifhess :   flag controlling evaluation of Hessian:
c                   ifhess = 0, do not compute Hessian
c                   ifhess = 1, compute Hessian
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DTAEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c   
c     OUTPUT:
c
c     pot    :    potential
c     fld    :    if (iffld .eq.1)   
c     hess   :    if (ifhess .eq.1)   
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,iffld,ifhess,nlege
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,ldiff
      double precision center(3),targ(3),wlege(*)
      double precision zdiff(3)
      double precision scarray(*),rscale
      double precision cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      double complex local(0:nterms,-nterms:nterms)
      double complex local2(0:2,-2:2)
      double complex z0,ima,ephi1,pot,fld(3),hess(6)
c
      double precision, allocatable :: pp(:,:)
      double precision, allocatable :: powers(:)
      double complex, allocatable :: ppc(:,:)
      double complex, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms,0:nterms))
      allocate(ppc(0:nterms,-nterms:nterms))
      allocate(powers(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
c
c     determine order of shifted expansion 
c
      ll = 0
      if (iffld.eq.1) ll = 1
      if (ifhess.eq.1) ll = 2
c
      do l = 0,ll
      do m = -l,l
         local2(l,m) = 0.0d0
      enddo
      enddo
c
      zdiff(1) = center(1) - targ(1)
      zdiff(2) = center(2) - targ(2)
      zdiff(3) = center(3) - targ(3)
      call cart2polarl(zdiff,d,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
      dd = d*rscale
      powers(0) = 1.0d0
      powers(1) = dd
      ephi(0) = 1.0d0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+1
         powers(l) = dd*powers(l-1)
         ephi(l) = ephi(l-1)*ephi(1)
         ephi(-l) = dconjg(ephi(l))
      enddo
c
      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
      do l = 0,nterms
         do k = -l,l
            ppc(l,k) = pp(l,abs(k))*powers(l)*ephi(k)
         enddo
      enddo
c
c     shift to local expansion of order ll about target point
c
c     old code (not unrolled).
c
ccc      iuse = 1
ccc      do j = 0,ll
ccc         do k = -j,j
ccc            do l = j,nterms
ccc               do m = -l,l
ccc                  ldiff = l-j
ccc                  mm = m-k
ccc                  if (abs(mm).le.ldiff) then
ccc                     local2(j,k) = local2(j,k) +
ccc     1               ppc(ldiff,mm)*scarray(iuse)*local(l,m)
ccc                  endif
ccc                  iuse = iuse+1
ccc               enddo
ccc            enddo
ccc         enddo
ccc      enddo
c
      iuse = 1
      do l = 0,nterms
         do m = -l,l
            local2(0,0) = local2(0,0) +
     1      ppc(l,m)*scarray(iuse)*local(l,m)
            iuse = iuse+1
         enddo
      enddo
      if (ll.ge.1) then
         do l = 1,nterms
            ldiff = l-1
            do m = -l,l-2
               mm = m+1
               local2(1,-1) = local2(1,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               local2(1,0) = local2(1,0) +
     1         ppc(ldiff,m)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+2,l
               mm = m-1
               local2(1,1) = local2(1,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
      if (ll.eq.2) then
         do l = 2,nterms
            ldiff = l-2
            do m = -l,l-4
               mm = m+2
               local2(2,-2) = local2(2,-2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+1,l-3
               mm = m+1
               local2(2,-1) = local2(2,-1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               mm = m
               local2(2,0) = local2(2,0) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+3,l-1
               mm = m-1
               local2(2,1) = local2(2,1) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
            do m = -l+4,l
               mm = m-2
               local2(2,2) = local2(2,2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(l,m)
               iuse = iuse+1
            enddo
         enddo
      endif
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
      pot = local2(0,0)
c
c     fld comes from l=1 modes
c
      rfac = rscale/sqrt(2.0d0)
      if (iffld.eq.1) then
         fld(3) = -local2(1,0)*rscale
         fld(1) = rfac*(local2(1,1) + local2(1,-1))
         fld(2) = rfac*ima*(local2(1,1) - local2(1,-1))
      endif
c
c     hess comes from l=2 modes
c
      if (ifhess.eq.1) then
ccc         rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
         rfac = rfac*rscale*sqrt(3.0d0)
         z0 = rscale*rscale*local2(2,0)
         hess(1) = rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(2) = -rfac*(local2(2,2) + local2(2,-2)) - z0
         hess(3) = 2*z0
         hess(4) = rfac*ima*(local2(2,2) - local2(2,-2))
         hess(5) = -rfac*(local2(2,1) + local2(2,-1))
         hess(6) = -rfac*ima*(local2(2,1) - local2(2,-1))
      endif
c
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine l3dformmp_qp_trunc(ier,rscale,sources,
     $                  quadstr,quadvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, isss, lsss
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      double precision center(3),sources(3,ns)
      double complex quadstr(ns)
      double precision quadvec(6,ns)
      double precision rscale
      double complex mpole(0:nterms,-nterms:nterms)
      double precision binom(0:240,0:4)
      integer nlege
      double precision wlege(*)
      double precision, allocatable :: w(:)
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c

c     carve up workspace:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifr=iephi+lephi
      lfr=(nterms+3)
c
      isss=ifr + lfr
      lsss = 5*(2*nterms+1)
c
      lused=isss + lsss
      allocate(w(lused))
c
ccc      call prinf(' in formmp lused is *',lused,1)

      call getsgnformpmp_quad(w(isss),nterms)

      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo


      do i = 1, ns

      call l3dformmp0_qp_trunc(ier,rscale,sources(1,i),
     $          quadstr(i),quadvec(1,i),
     1		center,nterms,mpole,wlege,nlege,w(ipp),w(iephi),
     2          w(ifr),binom,w(isss))

      enddo
c
      return
      end
C
C***********************************************************************
      subroutine l3dformmp_qp_add_trunc(ier,rscale,sources,
     $                  quadstr,quadvec,ns,
     1                  center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS 
c     quadrupoles sources located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     rscale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     quadvec(6,ns)    : quadrupoles vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		        ier=0  returned successfully
c		        deprecated but left in calling sequence for
c		        backward compatibility.
c
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier, lused, isss, lsss
      integer jer,ipp,lpp,iephi,lephi,ifr,lfr
      double precision center(3),sources(3,ns)
      double complex quadstr(ns)
      double precision quadvec(6,ns)
      double precision rscale
      double complex mpole(0:nterms,-nterms:nterms)
      double complex, allocatable :: mptemp(:,:)
      integer nlege
      double precision wlege(*)
C
        allocate( mptemp(0:nterms,-nterms:nterms) )
C
c      do l = 0,nterms
c         do m=-l,l
c            mpole(l,m) = 0.0d0
c         enddo
c      enddo
c
        call l3dformmp_qp_trunc(ier,rscale,sources,
     $                  quadstr,quadvec,ns,
     1                  center,nterms,mptemp,wlege,nlege)

        do l = 0,nterms
          do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
          enddo
        enddo
c
      return
      end
C
c
c
c**********************************************************************
      subroutine l3dformmp0_qp_trunc(ier,rscale,source,quadstr,quadvec,
     1		center,nterms,mpole,wlege,nlege,pp,ephi,fr,binom,sss)
c**********************************************************************
c
c     See l3dformmp_qp for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,l,ll,lnew,m,mm,mnew1,mnew2
      double precision binom(0:240,0:4)
      double precision rscale,source(3),center(3),zdiff(3)
      double complex quadstr
      double precision quadvec(6)
      double precision pp(0:nterms,0:nterms)
      double precision fr(0:nterms+1)
      double precision cscale,cscale1,cscale2,cscale3,rtmp
      double precision sss(-2:2,-nterms:nterms)
      double precision a22,a21,a20,rmul,dd
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double complex mp2(-2:2)
      double complex mpole(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-1:nterms+1),ephi1
      double complex eye,zmul,ztmp,ctmp,ztmp0,ztmp1,ztmp2,ztmp3
      integer nlege
      double precision wlege(*)
      data eye/(0.0d0,1.0d0)/
c
c     
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
c     using truncated version of full multipole-multipole shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
      ephi(1) = dcmplx(cphi,sphi)
      ephi(-1) = dconjg(ephi(1))
c
      fr(0) = 1.0D0
      dd = r*rscale
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      do l = 2,nterms
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
c       
c      do l = 0,nterms
c      fr(l)=fr(l)/sqrt(2*l+1.0d0) *rscale*rscale
c      enddo
c      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
      do l = 0,nterms
      fr(l)=fr(l) *rscale*rscale
      enddo
      call ylgndrufw(nterms,ctheta,pp,wlege,nlege)
c
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
        if( 1 .eq. 2 ) then
c
c       ... reference code, complex-valued quadrupoles
c
         do m = -2,2
            do ll = 0,nterms-2
               cscale = fr(ll)
               lnew = 2+ll
c
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            pp(ll,0)*cscale2*mp2(m)*cscale*sss(m,0)
     $               *quadstr
               do mm = 1,ll
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,2-m)*binom(lnew+mnew1,2+m)
                  cscale2 = cscale2*cscale*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*cscale*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            pp(ll,mm)*ephi(-mm)*mp2(m)
     $               *quadstr
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            pp(ll,mm)*ephi(mm)*mp2(m)
     $               *quadstr
               enddo
            enddo
         enddo
         endif
c
c
        if( 2 .eq. 2 ) then
c
c       ... optimized reference code, complex-valued quadrupoles
c       change do loop ordering
c       
         do ll = 0,nterms-2
            cscale = fr(ll)
            lnew = 2+ll
c
            ctmp = pp(ll,0)*cscale
            do m = -2,2
               cscale2 = binom(lnew-m,2-m)*binom(lnew+m,2+m)
               mpole(lnew,m) = mpole(lnew,m) +
     1            ctmp*cscale2*mp2(m)*sss(m,0)
     $               *quadstr
            enddo

            do mm = 1,ll
               ctmp = pp(ll,mm)*ephi(-mm)*cscale
               do m = -2,2
                  mnew1 = m+mm
                  mnew2 = m-mm
                  cscale2 = binom(lnew-mnew1,2-m)*binom(lnew+mnew1,2+m)
                  cscale2 = cscale2*sss(m,+mm)
                  cscale3 = binom(lnew-mnew2,2-m)*binom(lnew+mnew2,2+m)
                  cscale3 = cscale3*sss(m,-mm)
                  
                  mpole(lnew,mnew1) = mpole(lnew,mnew1)+ cscale2*
     1            ctmp*mp2(m)*quadstr
                  mpole(lnew,mnew2) = mpole(lnew,mnew2)+ cscale3*
     1            dconjg(ctmp)*mp2(m)*quadstr
               enddo
            enddo
            
         enddo

        endif
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dformta_qp_trunc(ier,rscale,sources,
     $          quadstr,quadvec,ns,
     &		center,nterms,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single quadrupole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     quadvec   : quadrupole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   ier=2 insufficient memory in workspace w
c     locexp : coefficients of the local expansion
c
c
c---------------------------------------------------------------------
      implicit none
      integer ier,nterms,ipp,lpp,iephi,lephi,ifr,lfr,lused,nlege
      integer isss,lsss,i,ns,l,m
      double precision rscale,sources(3,ns),center(3),wlege(*)
      double precision, allocatable :: w(:)
      double complex quadstr(ns)
      double precision quadvec(6,ns)
      double complex locexp(0:nterms,-nterms:nterms)
      double precision binom(0:240,0:4)
      double precision sss(-2:2,-nterms:nterms)
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+3)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+3)+7
c
      ifr=iephi+lephi
      lfr=2*(nterms+5)
c
      isss=ifr + lfr
      lsss = 5*(2*nterms+1)
c
      lused=isss+lsss
      allocate(w(lused))
c
      call getsgnformpmp_quad(w(isss),nterms)
c
      binom(0,0) = 1.0d0
      do i = 1,2*nterms
         binom(i,0) = 1.0d0
         binom(i,1) = dsqrt(1.0d0*i)
      enddo
      do i = 2,2*nterms
         binom(i,2) = dsqrt((i*(i-1))/2.0d0)
      enddo
      do i = 3,2*nterms
         binom(i,3) = dsqrt((i*(i-1)*(i-2))/6.0d0)
      enddo
      do i = 4,2*nterms
         binom(i,4) = dsqrt((i*(i-1)*(i-2)*(i-3))/24.0d0)
      enddo

      do i=1,ns

      call l3dformta0_qp_trunc(rscale,sources(1,i),
     $   quadstr(i),quadvec(1,i),
     $   center,nterms,locexp,w(ipp),w(iephi),w(ifr),
     $   wlege,nlege,binom,w(isss))
c
      enddo
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dformta_qp_add_trunc(ier,rscale,sources,
     $                     quadstr,quadvec,ns,
     1		           center,nterms,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS quadrupoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformta1_quad/l3dformta0_quad below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     rscale   : scaling parameter
c     sources   : coordinates of the sources
c     quadvec    : quadrupole vector
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the expansion
c     wlege  :   precomputed array of recurrence relation coeffs
c                for Ynm calculation.
c     nlege  : dimension parameter for wlege
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c
c     locexp    : coeffs for the expansion
c     lused     : amount of work space "w" used
c
c
c----------------------------------------------------------------------
      implicit none
      integer ier,ns,nterms,i,l,m,nlege
      double precision rscale,sources(3,ns),center(3),rs
      double complex quadstr(ns)
      double precision quadvec(6,ns),wlege(*)
      double complex locexp(0:nterms,-nterms:nterms)
      double complex eye
      double complex, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
c      do l = 0,nterms
c         do m = -l,l
c            locexp(l,m) = 0.0d0
c         enddo
c      enddo
c

        call l3dformta_qp_trunc(ier,rscale,sources,
     1		quadstr,quadvec,ns,center,nterms,mptemp,wlege,nlege)
c

      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)+mptemp(l,m)
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine l3dformta0_qp_trunc(rscale,source,quadstr,quadvec,
     &		center,nterms,locexp,pp,ephi,fr,wlege,nlege,binom,sss)
c**********************************************************************
c
c     See l3dformta_quad/l3dformta1_quad for comments
c
c---------------------------------------------------------------------
      implicit none
      integer nterms,i,nt2
      integer l,m,lnew,mnew,ll,mm,nlege
      double precision binom(0:240,0:4),wlege(*)
      double precision sss(-2:2,-nterms:nterms)
      double precision rscale,source(3),center(3),zdiff(3)
      double complex quadstr
      double precision quadvec(6)
      double precision pp(0:nterms+2,0:nterms+2)
      double precision fr(0:nterms+3)
      double precision cscale,cscale2,cscale3
      double precision a22,a21,a20,rmul,dd
      double precision r,theta,phi,ctheta,stheta,cphi,sphi
      double precision c(0:200,0:200),sqc(0:200,0:200)
      double complex mp2(-2:2),mpolelm
      double complex locexp(0:nterms,-nterms:nterms)
      double complex ephi(-nterms-2:nterms+2),ephi1
      double complex eye,zmul
      data eye/(0.0d0,1.0d0)/
c

      nt2 = nterms+2
c
      do l=0,2*nt2
         c(l,0)=1.0d0
         sqc(l,0)=1.0d0
      enddo
      do m=1,2*nt2
         c(m,m)=1.0d0
         sqc(m,m)=1.0d0
         do l=m+1,2*nt2
            c(l,m)=c(l-1,m)+c(l-1,m-1)
            sqc(l,m)=dsqrt(c(l,m))
         enddo
      enddo
C
c
c     first convert quadrupole vector contributions to standard
c     n=2 moments about source position using standard
c     d+,d-,dz operators.
c
      a22 = 1.0d0/dsqrt(24.0d0)
      a21 = 1.0d0/dsqrt(6.0d0)
      a20 = 0.5d0
c
      rmul = 0.25d0*quadvec(1)
      mp2(2) = rmul/a22
      mp2(0) = -2*rmul/a20
      mp2(-2) = rmul/a22
c
      rmul = 0.25d0*quadvec(2)
      mp2(2) = mp2(2) -rmul/a22
      mp2(0) = mp2(0) -2*rmul/a20
      mp2(-2) = mp2(-2) -rmul/a22
c
      rmul = quadvec(3)
      mp2(0) = mp2(0) + rmul/a20
c
      zmul = -eye*0.25d0*quadvec(4)
      mp2(2) = mp2(2) +zmul/a22
      mp2(-2) = mp2(-2) -zmul/a22
c
      rmul = -quadvec(5)/2
      mp2(1) = rmul/a21
      mp2(-1) = rmul/a21
c
      zmul = -quadvec(6)/(2*eye)
      mp2(1) = mp2(1) + zmul/a21
      mp2(-1) = mp2(-1) - zmul/a21
c
c     now shift n=2 contributions to expansion center
c     using truncated version of full multipole-local shift.
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
      zdiff(1)=-zdiff(1)
      zdiff(2)=-zdiff(2)
      zdiff(3)=-zdiff(3)
      call cart2polarl(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
      fr(0) = 1.0D0
      dd = 1.0d0/(r*rscale)
      fr(1) = dd
      ephi(0) = 1.0D0
      ephi(1) = ephi1
      ephi(-1) = dconjg(ephi1)
      do l = 2,nterms+3
         fr(l) = fr(l-1)*dd
         ephi(l) = ephi(l-1)*ephi1
         ephi(-l) = dconjg(ephi(l))
      enddo
C
      do l = 0,nterms+3
      fr(l)=fr(l) *rscale*rscale
      enddo
c
      call ylgndrufw(nt2,ctheta,pp,wlege,nlege)
c
      if( 2 .eq. 2 ) then
c
c     reference code
C
C---- go through terms in expansions MPOLE
C     generating appropriate terms in new expansions.
C
      do m = -2,2
	 mpolelm = mp2(m)*rscale*quadstr
         do lnew = 0,nterms
            do mnew = -lnew,lnew
               ll = 2 + lnew
               mm = mnew - m
               cscale = fr(ll+1)*sqc(ll+mm,lnew+mnew)*
     1                     sqc(ll-mm,lnew-mnew)
               cscale = cscale*(-1)**lnew
               cscale = cscale*sss(m,-mnew)
c               if ( (m .lt. 0) .and. (mnew .lt. 0) ) then
c                  if (-mnew .lt. -m)  cscale = cscale*(-1)**mnew
c                  if (-mnew .ge. -m)  cscale = cscale*(-1)**m
c               endif
c               if ( (m .gt. 0) .and. (mnew .gt. 0) ) then
c                  if (mnew .lt. m) cscale = cscale*(-1)**mnew
c                  if (mnew .ge. m) cscale = cscale*(-1)**m
c               endif
               if (mm .eq. 0) then
                  locexp(lnew,mnew) = locexp(lnew,mnew) +
     1            pp(ll,0)*cscale*mpolelm
               else  if (mm .gt. 0) then
                  locexp(lnew,mnew) = locexp(lnew,mnew)+ cscale*
     1            pp(ll,mm)*ephi(-mm)*mpolelm
               else
                  locexp(lnew,mnew) = locexp(lnew,mnew)+ cscale*
     1            pp(ll,-mm)*ephi(-mm)*mpolelm
               endif
            enddo
         enddo
      enddo
c
      endif
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dallhess_qp(iffld,ifhess,
     $     sources,quadstr,quadvec,ns,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD and
c     hessian HESS at the target point TARGET, due to a collection of
c     quadrupoles at SOURCE(3,ns).
c     
c               pot = quadvec(1)*V_xx +
c                     quadvec(2)*V_yy +
c                     quadvec(3)*V_zz +
c                     quadvec(4)*V_xy +
c                     quadvec(5)*V_xz +
c                     quadvec(6)*V_yz 
c
c      V_xx = (-1/r^3 + 3*dx**2/r^5)
c      V_xy = 3*dx*dy/r^5
c      V_xz = 3*dx*dz/r^5
c      V_yy = (-1/r^3 + 3*dy**2/r^5)
c      V_yz = 3*dy*dz/r^5
c      V_zz = (-1/r^3 + 3*dz**2/r^5)
c
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld         : flag for computing -gradient
c	                 	   iffld = 0 -> dont compute 
c		                   iffld = 1 -> do compute 
c     ifhess       : flag for computing Hessian
c	                 	ifhess = 0 -> dont compute 
c		                ifhess = 1 -> do compute 
c     sources(3,ns) : location of the sources
c     quadstr(ns)   : quadrupole strength
c     quadvec(6,ns) : quadrupole vector
c     ns            : number of sources
c     target(3)     : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot           : calculated potential
c     fld           : calculated -gradient
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,i,j,ns
      double precision sources(3,ns),target(3)
      double precision quadvec(6,ns)
      double complex quadstr(ns)
      double complex pot,fld(3),hess(6),potloc,fldloc(3),hessloc(6)
c
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      if (ifhess.eq.1) then
         do i = 1,6
            hess(i) = 0.0d0
         enddo
      endif
c
      do i = 1,ns
         call lpotfld3dhess_qp(iffld,ifhess,
     $        sources(1,i),quadstr(i),quadvec(1,i),
     1        target,potloc,fldloc,hessloc)
         pot = pot + potloc
         if (iffld.eq.1) then
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
         endif
         if (ifhess.eq.1) then
            do j = 1,6
               hess(j) = hess(j) + hessloc(j)
            enddo
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotfld3dhess_qp(iffld,ifhess,
     $     source,quadstr,quadvec,target,pot,fld,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field FLD and
c     hessian HESS at the target point TARGET, due to a quadrupole at
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c    
c               pot = quadvec(1)*V_xx +
c                     quadvec(2)*V_yy +
c                     quadvec(3)*V_zz +
c                     quadvec(4)*V_xy +
c                     quadvec(5)*V_xz +
c                     quadvec(6)*V_yz 
c
c      V_xx = (-1/r^3 + 3*dx**2/r^5)
c      V_xy = 3*dx*dy/r^5
c      V_xz = 3*dx*dz/r^5
c      V_yy = (-1/r^3 + 3*dy**2/r^5)
c      V_yz = 3*dy*dz/r^5
c      V_zz = (-1/r^3 + 3*dz**2/r^5)
c
c		fld = -grad(pot)
c		hess = (potxx,potyy,potzz,potxy,potxz,potyz)
c
c----------------------------------------------------------------------
c     INPUT:
c
c     iffld        : flag for computing gradient
c	                 	ffld = 0 -> dont compute 
c		                ffld = 1 -> do compute 
c     source(3)    : location of the source 
c     quadstr      : quadrupole strength
c     quadvec(6)   : quadrupole vector
c     target(3)    : location of the target
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     pot          : calculated potential
c     fld          : calculated -gradient
c     hess         : calculated hessian
c
c----------------------------------------------------------------------
      implicit none
      integer iffld,ifhess,i
      double precision source(3),target(3)
      double precision quadvec(6),rr(3),rx,ry,rz
      double precision cd,d5,d3,v11,v12,v13,v22,v23,v33
      double precision d7,v111,v112,v113,v122,v123,v133
      double precision v222,v223,v233,v333
      double precision d9,v1111,v1112,v1113,v1122,v1123,v1133,v1222
      double precision v1223,v1233,v1333,v2222,v2223,v2233,v2333,v3333
      double complex quadstr
      double complex pot,fld(3),hess(6)
c
c
c ... Caculate offsets and distance
c
      rr(1)=target(1)-source(1)
      rr(2)=target(2)-source(2)
      rr(3)=target(3)-source(3)
      cd = sqrt(rr(1)*rr(1)+rr(2)*rr(2)+rr(3)*rr(3))
      d5 = 1.0d0/(cd**5)
      d3 = 1.0d0/(cd**3)
c
c ... Calculate the potential and field in the regular case:
c
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      v11 = -d3 + 3*d5*rr(1)**2
      v12 = 3*d5*rr(1)*rr(2)
      v13 = 3*d5*rr(1)*rr(3)
      v22 = -d3 + 3*d5*rr(2)**2
      v23 = 3*d5*rr(2)*rr(3)
      v33 = -d3 + 3*d5*rr(3)**2
c
      pot=v11*quadvec(1)+v22*quadvec(2)+v33*quadvec(3)
     $   +v12*quadvec(4)+v13*quadvec(5)+v23*quadvec(6)

      if (iffld.eq.1) then
         d7 = 1.0d0/(cd**7)
         v111 = 9*rr(1)*d5 - 15*d7*rr(1)**3
         v112 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(1)**2
         v113 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(1)**2
         v122 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(2)**2
         v123 = -15*d7*rr(1)*rr(2)*rr(3)
         v133 = 3*rr(1)*d5 - 15*d7*rr(1)*rr(3)**2
         v222 = 9*rr(2)*d5 - 15*d7*rr(2)**3
         v223 = 3*rr(3)*d5 - 15*d7*rr(3)*rr(2)**2
         v233 = 3*rr(2)*d5 - 15*d7*rr(2)*rr(3)**2
         v333 = 9*rr(3)*d5 - 15*d7*rr(3)**3
         fld(1)=v111*quadvec(1)+v122*quadvec(2)+v133*quadvec(3)
     $         +v112*quadvec(4)+v113*quadvec(5)+v123*quadvec(6)
         fld(2)=v112*quadvec(1)+v222*quadvec(2)+v233*quadvec(3)
     $         +v122*quadvec(4)+v123*quadvec(5)+v223*quadvec(6)
         fld(3)=v113*quadvec(1)+v223*quadvec(2)+v333*quadvec(3)
     $         +v123*quadvec(4)+v133*quadvec(5)+v233*quadvec(6)
         fld(1) = -fld(1)
         fld(2) = -fld(2)
         fld(3) = -fld(3)
      endif 

      if (ifhess.eq.1) then
         d7 = 1.0d0/(cd**7)
         d9 = 1.0d0/(cd**9)
         rx = rr(1)
         ry = rr(2)
         rz = rr(3)
         v1111 = 105*rx**4*d9-90*rx**2*d7+9*d5
         v1112 = 105*rx**3*ry*d9-45*rx*ry*d7
         v1113 = 105*rx**3*rz*d9-45*rx*rz*d7
         v1122 = 105*rx**2*ry**2*d9-15*(rx**2+ry**2)*d7+3*d5
         v1123 = 105*rx**2*ry*rz*d9-15*ry*rz*d7
         v1133 = 105*rx**2*rz**2*d9-15*(rx**2+rz**2)*d7+3*d5
         v1222 = 105*rx*ry**3*d9-45*rx*ry*d7
         v1223 = 105*rx*ry**2*rz*d9-15*rx*rz*d7
         v1233 = 105*rx*ry*rz**2*d9-15*rx*ry*d7
         v1333 = 105*rx*rz**3*d9-45*rx*rz*d7
         v2222 = 105*ry**4*d9-90*ry**2*d7+9*d5
         v2223 = 105*ry**3*rz*d9-45*ry*rz*d7
         v2233 = 105*ry**2*rz**2*d9-15*(ry**2+rz**2)*d7+3*d5
         v2333 = 105*ry*rz**3*d9-45*ry*rz*d7
         v3333 = 105*rz**4*d9-90*rz**2*d7+9*d5
c         hess(1)=v1111*quadvec(1)+v1221*quadvec(2)+v1331*quadvec(3)
c     $          +v1121*quadvec(4)+v1131*quadvec(5)+v1231*quadvec(6)
c         hess(2)=v1122*quadvec(1)+v2222*quadvec(2)+v2332*quadvec(3)
c     $          +v1222*quadvec(4)+v1232*quadvec(5)+v2232*quadvec(6)
c         hess(3)=v1133*quadvec(1)+v2233*quadvec(2)+v3333*quadvec(3)
c     $          +v1233*quadvec(4)+v1333*quadvec(5)+v2333*quadvec(6)
c         hess(4)=v1112*quadvec(1)+v1222*quadvec(2)+v1332*quadvec(3)
c     $          +v1122*quadvec(4)+v1132*quadvec(5)+v1232*quadvec(6)
c         hess(5)=v1113*quadvec(1)+v1223*quadvec(2)+v1333*quadvec(3)
c     $          +v1123*quadvec(4)+v1133*quadvec(5)+v1233*quadvec(6)
c         hess(6)=v1123*quadvec(1)+v2223*quadvec(2)+v2333*quadvec(3)
c     $          +v1223*quadvec(4)+v1233*quadvec(5)+v2233*quadvec(6)
         hess(1)=v1111*quadvec(1)+v1122*quadvec(2)+v1133*quadvec(3)
     $          +v1112*quadvec(4)+v1113*quadvec(5)+v1123*quadvec(6)
         hess(2)=v1122*quadvec(1)+v2222*quadvec(2)+v2233*quadvec(3)
     $          +v1222*quadvec(4)+v1223*quadvec(5)+v2223*quadvec(6)
         hess(3)=v1133*quadvec(1)+v2233*quadvec(2)+v3333*quadvec(3)
     $          +v1233*quadvec(4)+v1333*quadvec(5)+v2333*quadvec(6)
         hess(4)=v1112*quadvec(1)+v1222*quadvec(2)+v1233*quadvec(3)
     $          +v1122*quadvec(4)+v1123*quadvec(5)+v1223*quadvec(6)
         hess(5)=v1113*quadvec(1)+v1223*quadvec(2)+v1333*quadvec(3)
     $          +v1123*quadvec(4)+v1133*quadvec(5)+v1233*quadvec(6)
         hess(6)=v1123*quadvec(1)+v2223*quadvec(2)+v2333*quadvec(3)
     $          +v1223*quadvec(4)+v1233*quadvec(5)+v2233*quadvec(6)
      endif 

        pot=pot*quadstr
        if( iffld .eq. 1 ) then
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        endif
        if( ifhess .eq. 1 ) then
        do i=1,6
        hess(i)=hess(i)*quadstr
        enddo
        endif

      return
      end
cc
c
c
c
