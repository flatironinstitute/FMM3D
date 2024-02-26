c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions.
c
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	hvec(n)= h_n(z)*scale^(n)
c       	jvec(n)= j_n(z)/scale^(n)
c
c          scale should be of the order of |z| if |z| < 1. Otherwise,
c          scale should be set to 1.
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
c         definition of Y_n^m for m<0. 
c
c         We also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the Y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  
c         More precisely, 
c
c                 \int_S Y_n^m Y_n^m d\Omega = 4 \pi. 
c
c         Using our standard definition, the addition theorem takes 
c         the simple form 
c
c         e^( i k r}/(ikr) = 
c         \sum_n \sum_m  j_n(k|S|) Y_l^m*(S) h_n(k|T|) Y_l^m(T)
c
c
c-----------------------------------------------------------------------
c
c      h3dmpevalp: computes potential
c                 due to a multipole expansion
c                 at a collection of targets.
c
c      h3dmpevalg: computes potential and gradients
c                          due to a multipole expansion
c                          at a collection of targets
c
c      h3dformmpc: creates multipole expansion due to 
c                 a collection of charges
c
c      h3dformmpd: creates multipole expansion due to 
c                 a collection of dipoles
c
c      h3dformmpcd: creates multipole expansion due
c                                to a collection of charges
c                                and dipoles
c
c      h3dtaevalp: computes potential
c                 due to a local expansion
c                 at a collection of targets.
c
c      h3dtaevalg: computes potential and gradients
c                          due to a local expansion
c                          at a collection of targets
c
c      h3dformtac: creates local expansion due to 
c                 a collection of charges
c
c      h3dformtad: creates local expansion due to 
c                 a collection of dipoles
c
c      h3dformtacd: creates local expansion due
c                                to a collection of charges
c                                and dipoles
c
c**********************************************************************
      subroutine h3dmpevalp(nd,zk,rscale,center,mpole,nterms,ztarg,
     1            ntarg,pot,wlege,nlege,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential
c     of an outgoing partial wave expansion and adds to existing
c     quantities.
c
c     pot = pot + sum sum  mpole(n,m) h_n(k r) Y_nm(theta,phi)
c                  n   m
c                 
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision 
c           where boxsize(0) is the size of the computation domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c     thresh :    threshold for computing outgoing expansion,
c                 potential at target won't be updated if
c                 |t - c| <= thresh, where t is the target
c                 location and, and c is the expansion center
c                 location
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
c
cc      calling sequence variables
c
      integer ntarg,nterms,nlege,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      complex *16 zk,pot(nd,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(*),thresh

c
cc      temporary variables
c
      integer idim
      complex *16, allocatable :: ephi(:),fhs(:)
      real *8, allocatable :: ynm(:,:)
      real *8 zdiff(3),done
      real *8 ctheta,stheta,r,theta,phi,cphi,sphi
      complex *16 ephi1,ur,utheta,uphi
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
      complex *16 fhder
      integer itarg,i,m,n,ifder
      data eye/(0.0d0,1.0d0)/

      done = 1
      ifder = 0

      allocate(ephi(-nterms-1:nterms+1))
      allocate(fhs(0:nterms))
      allocate(ynm(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta= dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
        if (r.lt.thresh)  goto 1000
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=done
        ephi(1)=ephi1
        cphi = dreal(ephi1)
        sphi = dimag(ephi1)
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms+1
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions:
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
c
c
c     get the Hankel functions:
c
        z=zk*r
        call h3dall(nterms,z,rscale,fhs,ifder,fhder)

        do idim=1,nd
          pot(idim,itarg)=pot(idim,itarg)+mpole(idim,0,0)*fhs(0)
        enddo
        do n=1,nterms
          ztmp1 = fhs(n)*ynm(n,0)
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+mpole(idim,n,0)*ztmp1
          enddo
          do m=1,n
            ztmp1=fhs(n)*ynm(n,m)
            do idim=1,nd
              ztmp2 = mpole(idim,n,m)*ephi(m) 
              ztmp3 = mpole(idim,n,-m)*ephi(-m)
              ztmpsum = ztmp2+ztmp3
              pot(idim,itarg)=pot(idim,itarg)+ztmp1*ztmpsum
            enddo
          enddo
        enddo
 1000 continue        
      enddo

      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine h3dmpevalg(nd,zk,rscale,center,mpole,nterms,ztarg,
     1            ntarg,pot,grad,wlege,nlege,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential, and gradient
c     of an outgoing partial wave expansion and adds to existing
c     quantities.
c
c     pot = pot + sum sum  mpole(n,m) h_n(k r) Y_nm(theta,phi)
c                  n   m
c                 
c     grad = grad + Gradient( sum sum  mpole(n,m) h_n(k r) Y_nm(theta,phi))
c                              n   m
c      If r < thresh 
c          then the subroutine does not update the potential and gradient
c          (recommended value = boxsize(0)*machine precision
c            where boxsize(0) is the size of the computational domain) 
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c     thresh :    threshold for computing outgoing expansion,
c                 potential at target won't be updated if
c                 |t - c| <= thresh, where t is the target
c                 location and, and c is the expansion center
c                 location
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad    :   updated gradient at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
c
cc      calling sequence variables
c
      integer ntarg,nterms,nlege,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      complex *16 zk,pot(nd,ntarg),grad(nd,3,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(*),thresh

c
cc      temporary variables
c
      complex *16, allocatable :: ephi(:),fhs(:),fhder(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      real *8 zdiff(3),done
      real *8 r, theta,phi,ctheta,stheta,cphi,sphi
      real *8 rx,ry,rz
      real *8 phix,phiy,phiz
      real *8 thetax,thetay,thetaz,dtmp
      complex *16 ephi1,ur(nd),utheta(nd),uphi(nd)
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
      complex *16 ztmp4,ztmp5,ztmp6
      integer itarg,i,m,n,ifder,idim
      data eye/(0.0d0,1.0d0)/

      done = 1
      ifder = 1

      allocate(ephi(-nterms-1:nterms+1))
      allocate(fhs(0:nterms),fhder(0:nterms))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta) 
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
        if (r.lt.thresh) goto 1000
c        
c     compute exp(eye*m*phi) array
c
        ephi(0)=done
        ephi(1)=ephi1
        cphi = dreal(ephi1)
        sphi = dimag(ephi1)
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms+1
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions:
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
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
c
c     get the Hankel functions:
c
        z=zk*r
        call h3dall(nterms,z,rscale,fhs,ifder,fhder)

c     initialize computed values and 
c     scale derivatives of Hankel functions so that they are
c     derivatives with respect to r.
c
        do n=0,nterms
          fhder(n)=fhder(n)*zk
        enddo

        do idim=1,nd
          ur(idim) = mpole(idim,0,0)*fhder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
          pot(idim,itarg)=pot(idim,itarg)+mpole(idim,0,0)*fhs(0)
        enddo

        do n=1,nterms
          ztmp1 = fhs(n)*ynm(n,0)
          ztmp2 = fhder(n)*ynm(n,0)
          ztmp3 = -fhs(n)*ynmd(n,0)*stheta
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+mpole(idim,n,0)*ztmp1
            ur(idim) = ur(idim) + mpole(idim,n,0)*ztmp2
            utheta(idim) = utheta(idim) + mpole(idim,n,0)*ztmp3
          enddo
          do m=1,n
            ztmp1=fhs(n)*ynm(n,m)*stheta
            ztmp4 = fhder(n)*ynm(n,m)*stheta
            ztmp5 = -fhs(n)*ynmd(n,m)
            ztmp6 = eye*m*fhs(n)*ynm(n,m)

            do idim=1,nd
              ztmp2 = mpole(idim,n,m)*ephi(m) 
              ztmp3 = mpole(idim,n,-m)*ephi(-m)
              ztmpsum = ztmp2+ztmp3
              pot(idim,itarg)=pot(idim,itarg)+ztmp1*ztmpsum
              ur(idim) = ur(idim) + ztmp4*ztmpsum
              utheta(idim) = utheta(idim) + ztmp5*ztmpsum
              ztmpsum = ztmp2 - ztmp3
              uphi(idim) = uphi(idim) + ztmp6*ztmpsum
            enddo
          enddo
        enddo
        do idim=1,nd
          grad(idim,1,itarg) = grad(idim,1,itarg) + ur(idim)*rx + 
     1               utheta(idim)*thetax + uphi(idim)*phix
          grad(idim,2,itarg) = grad(idim,2,itarg) + ur(idim)*ry + 
     1               utheta(idim)*thetay + uphi(idim)*phiy
          grad(idim,3,itarg) = grad(idim,3,itarg) + ur(idim)*rz + 
     1               utheta(idim)*thetaz + uphi(idim)*phiz
        enddo
 1000 continue        
      enddo

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
      subroutine h3dformmpc(nd,zk,rscale,sources,charge,
     1            ns,center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Adds to multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the h-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 charge(nd,ns)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:)
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      complex *16, allocatable :: fjs(:),fjder(:),ephi(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder,idim

      ifder=0
      allocate(fjs(0:nterms),fjder(0:nterms),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)

c
c     get Bessel functions
c
        z=zk*r
        call besseljs3d(nterms,z,rscale,fjs,ifder,fjder)
c
c
c     multiply all jn by charge strength and (i*k).
c
        do n = 0,nterms
          fjs(n) = fjs(n)*zkeye
        enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|)Ylm(T)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of already above
c 
        do idim=1,nd
          mpole(idim,0,0)= mpole(idim,0,0) + fjs(0)*charge(idim,isrc)
        enddo
        do n=1,nterms
          ztmp = ynm(n,0)*fjs(n)
          do idim=1,nd
            mpole(idim,n,0)= mpole(idim,n,0) + ztmp*charge(idim,isrc)
          enddo
          do m=1,n
            ztmp=ynm(n,m)*fjs(n)
            do idim=1,nd
              mpole(idim,n,m)= mpole(idim,n,m) +
     1           ztmp*ephi(-m)*charge(idim,isrc)
              mpole(idim,n,-m)= mpole(idim,n,-m) + 
     1           ztmp*ephi(m)*charge(idim,isrc)
            enddo
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
C***********************************************************************
      subroutine h3dformmpd(nd,zk,rscale,sources,
     1            dipvec,ns,center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Adds to multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(3,ns)    : dipole orientation vectors
C     ns              : number of sources
C     center(3)       : expansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the h-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 dipvec(nd,3,ns)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: ephi(:),fjs(:),fjder(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye,ztmp1
      complex *16 fjsuse
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 rx,thetax,phix,ry,thetay,phiy,rz,thetaz,phiz
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder,idim

      ifder=1

      allocate(fjs(0:nterms+1),fjder(0:nterms+1),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms),ynmd(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)

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
c     we use the recurrence relation 
c
c     (2n+1)fjs_n(kr)/(kr) = fjs(n+1)*rscale + fjs(n-1)/rscale
c
c     to avoid division by r. The variable fjsuse is set to fjs(n)/r:
c
c           fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
c	    fjsuse = wavek*fjsuse/(2*n+1.0d0)
c

        rx = stheta*cphi
        thetax = ctheta*cphi
        phix = -sphi

        ry = stheta*sphi
        thetay = ctheta*sphi
        phiy = cphi

        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0
        
        
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
c
c     get Bessel functions:
c
        z=zk*r
        call besseljs3d(nterms+1,z,rscale,fjs,ifder,fjder)
c
c
c     multiply all jn by  (i*k).
c
        do n = 0,nterms
          fjs(n) = fjs(n)*zkeye
          fjder(n) = fjder(n)*zkeye*zk
        enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|)Ylm(T)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of already above
c

        ur = ynm(0,0)*fjder(0)
        do idim=1,nd
          zzz = ur*(dipvec(idim,1,isrc)*rx + dipvec(idim,2,isrc)*ry +
     1          dipvec(idim,3,isrc)*rz)
          mpole(idim,0,0)= mpole(idim,0,0)+ zzz
        enddo
        do n=1,nterms
          fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
          fjsuse = fjsuse*zk/(2*n+1.0d0)
          ur = ynm(n,0)*fjder(n)
          utheta = -fjsuse*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax
          uy = ur*ry + utheta*thetay
          uz = ur*rz + utheta*thetaz
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy +
     1               dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + zzz 
          enddo
          do m=1,n
            ur = fjder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fjsuse*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fjsuse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              mpole(idim,n, m)= mpole(idim,n, m) + zzz
            enddo

            ur = fjder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fjsuse*ynmd(n,m)
            uphi = eye*m*ephi(m)*fjsuse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux+dipvec(idim,2,isrc)*uy + 
     1           dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m) + 
     1           zzz
            enddo
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
c**********************************************************************
c
c
c
c
C***********************************************************************
      subroutine h3dformmpcd(nd,zk,rscale,sources,charge,
     1            dipvec,ns,center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Adds to multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths
C     dipvec(3,ns)    : dipole orientation vectors
C     ns              : number of sources
C     center(3)       : expansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the h-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 charge(nd,ns)
      complex *16 dipvec(nd,3,ns)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: ephi(:),fjs(:),fjder(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye,ztmp1
      complex *16 fjsuse
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 rx,thetax,phix,ry,thetay,phiy,rz,thetaz,phiz
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder,idim

      ifder=1

      allocate(fjs(0:nterms+1),fjder(0:nterms+1),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms),ynmd(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)

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
c     we use the recurrence relation 
c
c     (2n+1)fjs_n(kr)/(kr) = fjs(n+1)*rscale + fjs(n-1)/rscale
c
c     to avoid division by r. The variable fjsuse is set to fjs(n)/r:
c
c           fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
c	    fjsuse = wavek*fjsuse/(2*n+1.0d0)
c

        rx = stheta*cphi
        thetax = ctheta*cphi
        phix = -sphi

        ry = stheta*sphi
        thetay = ctheta*sphi
        phiy = cphi

        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0
        
        
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
c
c     get Bessel functions:
c
        z=zk*r
        call besseljs3d(nterms+1,z,rscale,fjs,ifder,fjder)
c
c
c     multiply all jn by  (i*k).
c
        do n = 0,nterms
          fjs(n) = fjs(n)*zkeye
          fjder(n) = fjder(n)*zkeye*zk
        enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|)Ylm(T)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of already above
c

        ur = ynm(0,0)*fjder(0)
        do idim=1,nd
          zzz = ur*(dipvec(idim,1,isrc)*rx + dipvec(idim,2,isrc)*ry +
     1          dipvec(idim,3,isrc)*rz)
          mpole(idim,0,0)= mpole(idim,0,0) + 
     1        fjs(0)*charge(idim,isrc) + zzz
        enddo
        do n=1,nterms
          fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
          fjsuse = fjsuse*zk/(2*n+1.0d0)
          ur = ynm(n,0)*fjder(n)
          utheta = -fjsuse*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax
          uy = ur*ry + utheta*thetay
          uz = ur*rz + utheta*thetaz
          ztmp = ynm(n,0)*fjs(n)
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy +
     1               dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + ztmp*charge(idim,isrc)+
     1                    zzz       
          enddo
          do m=1,n
            ztmp=ynm(n,m)*stheta*fjs(n)
            ztmp1 = ztmp*ephi(-m)
            ur = fjder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fjsuse*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fjsuse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              mpole(idim,n, m)= mpole(idim,n, m) +
     1        ztmp1*charge(idim,isrc) + zzz
            enddo

            ztmp1=ztmp*ephi(m)
            ur = fjder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fjsuse*ynmd(n,m)
            uphi = eye*m*ephi(m)*fjsuse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux+dipvec(idim,2,isrc)*uy + 
     1           dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m) + 
     1           ztmp1*charge(idim,isrc) + zzz
            enddo
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
c**********************************************************************
      subroutine h3dtaevalp(nd,zk,rscale,center,locexp,nterms,ztarg,
     1            ntarg,pot,wlege,nlege)
c**********************************************************************
c
c     This subroutine evaluates the potential
c     of an incoming partial wave expansion and adds to existing
c     quantities.
c
c     pot = pot + sum sum  locexp(n,m) j_n(k r) Y_nm(theta,phi)
c                  n   m
c                 
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c            where boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of local expansions
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     lcoexp :    local expansion
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
c
cc      calling sequence variables
c
      integer ntarg,nterms,nlege,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      complex *16 zk,pot(nd,ntarg)
      complex *16 locexp(nd,0:nterms,-nterms:nterms)
      real *8 wlege(*)

c
cc      temporary variables
c
      complex *16, allocatable :: ephi(:),fjs(:)
      real *8, allocatable :: ynm(:,:)
      real *8 zdiff(3),done
      real *8 ctheta,stheta,r,theta,phi,cphi,sphi
      complex *16 ephi1,ur,utheta,uphi
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
      complex *16 fjder
      integer itarg,i,m,n,ifder,idim
      data eye/(0.0d0,1.0d0)/

      done = 1
      ifder = 0

      allocate(ephi(-nterms-1:nterms+1))
      allocate(fjs(0:nterms))
      allocate(ynm(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta= dsin(theta)
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
        do i=2,nterms+1
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions:
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
c
c
c     get the Bessel functions:
c
        z=zk*r
        call besseljs3d(nterms,z,rscale,fjs,ifder,fjder)

        do idim=1,nd
          pot(idim,itarg)=pot(idim,itarg)+locexp(idim,0,0)*fjs(0)
        enddo
        do n=1,nterms
          ztmp1 = fjs(n)*ynm(n,0)
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+locexp(idim,n,0)*ztmp1
          enddo
          do m=1,n
	        ztmp1=fjs(n)*ynm(n,m)
            do idim=1,nd
              ztmpsum = locexp(idim,n,m)*ephi(m) + 
     1            locexp(idim,n,-m)*ephi(-m)
              pot(idim,itarg)=pot(idim,itarg)+ztmp1*ztmpsum
            enddo
          enddo
        enddo
 1000 continue        
      enddo

      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine h3dtaevalg(nd,zk,rscale,center,locexp,nterms,
     1            ztarg,ntarg,pot,grad,wlege,nlege)
c**********************************************************************
c
c     This subroutine evaluates the potential, and gradient
c     of an incoming partial wave expansion and adds to existing
c     quantities.
c
c     pot = pot + sum sum  locexp(n,m) j_n(k r) Y_nm(theta,phi)
c                  n   m
c                 
c     grad = grad + Gradient( sum sum  locexp(n,m) j_n(k r) Y_nm(theta,phi))
c                              n   m
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of local expansions
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     locexp :    multipole expansion
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad    :   updated gradient at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
c
cc      calling sequence variables
c
      integer ntarg,nterms,nlege,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      complex *16 zk,pot(nd,ntarg),grad(nd,3,ntarg)
      complex *16 locexp(nd,0:nterms,-nterms:nterms)
      real *8 wlege(*),thresh

c
cc      temporary variables
c
      complex *16, allocatable :: ephi(:),fjs(:),fjder(:)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      real *8 zdiff(3),done
      real *8 r, theta,phi,ctheta,stheta,cphi,sphi
      real *8 rx,ry,rz
      real *8 phix,phiy,phiz
      real *8 thetax,thetay,thetaz,dtmp
      complex *16 ephi1,ur(nd),utheta(nd),uphi(nd)
      complex *16 fjsuse
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
      complex *16 ztmp4,ztmp5,ztmp6
      integer itarg,i,m,n,ifder,idim
      data eye/(0.0d0,1.0d0)/

      done = 1
      ifder = 1

      allocate(ephi(-nterms-1:nterms+1))
      allocate(fjs(0:nterms+1),fjder(0:nterms+1))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta) 
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
        do i=2,nterms+1
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions:
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
        rx = stheta*cphi
        thetax = ctheta*cphi
        phix = -sphi
        ry = stheta*sphi
        thetay = ctheta*sphi
        phiy = cphi
        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0
c
c
c
c     get the Bessel functions:
c
        z=zk*r
        call besseljs3d(nterms+1,z,rscale,fjs,ifder,fjder)

c     initialize computed values and 
c     scale derivatives of Hankel functions so that they are
c     derivatives with respect to r.
c
        do n=0,nterms
          fjder(n)=fjder(n)*zk
        enddo

        do idim=1,nd
          ur(idim) = locexp(idim,0,0)*fjder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
          pot(idim,itarg)=pot(idim,itarg)+locexp(idim,0,0)*fjs(0)
        enddo

        do n=1,nterms
          fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
          fjsuse = fjsuse*zk/(2*n+1.0d0)
          ztmp1 = fjs(n)*ynm(n,0)
          ztmp2 = fjder(n)*ynm(n,0)
          ztmp3 = -fjsuse*ynmd(n,0)*stheta
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+locexp(idim,n,0)*ztmp1
	        ur(idim) = ur(idim) + locexp(idim,n,0)*ztmp2
            utheta(idim) = utheta(idim) + locexp(idim,n,0)*ztmp3
          enddo
          do m=1,n
	        ztmp1=fjs(n)*ynm(n,m)*stheta
            ztmp4 = fjder(n)*ynm(n,m)*stheta
            ztmp5 = -fjsuse*ynmd(n,m)
            ztmp6 = eye*m*fjsuse*ynm(n,m)

            do idim=1,nd
              ztmp2 = locexp(idim,n,m)*ephi(m) 
	          ztmp3 = locexp(idim,n,-m)*ephi(-m)
	          ztmpsum = ztmp2+ztmp3
              pot(idim,itarg)=pot(idim,itarg) + ztmp1*ztmpsum
	          ur(idim) = ur(idim) + ztmp4*ztmpsum
	          utheta(idim) = utheta(idim) + ztmp5*ztmpsum
	          ztmpsum = ztmp2 - ztmp3
	          uphi(idim) = uphi(idim) + ztmp6*ztmpsum
            enddo
          enddo
        enddo
        do idim=1,nd
          grad(idim,1,itarg) = grad(idim,1,itarg)+ur(idim)*rx + 
     1          utheta(idim)*thetax + uphi(idim)*phix
          grad(idim,2,itarg) = grad(idim,2,itarg)+ur(idim)*ry + 
     1          utheta(idim)*thetay + uphi(idim)*phiy
          grad(idim,3,itarg) = grad(idim,3,itarg)+ur(idim)*rz + 
     1          utheta(idim)*thetaz + uphi(idim)*phiz

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
c
C***********************************************************************
      subroutine h3dformtac(nd,zk,rscale,sources,charge,
     1            ns,center,nterms,locexp,wlege,nlege)
C***********************************************************************
C
C     Adds to local (j) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of local expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     locexp          : coeffs of the j-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 charge(nd,ns)
      complex *16 locexp(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:)
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      complex *16, allocatable :: fhs(:),fhder(:),ephi(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder,idim

      ifder=0
      allocate(fhs(0:nterms),fhder(0:nterms),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)

c
c     get Hankel functions
c
        z=zk*r
        call h3dall(nterms,z,rscale,fhs,ifder,fhder)
c
c
c     multiply all hn by charge strength and (i*k).
c
        do n = 0,nterms
          fhs(n) = fhs(n)*zkeye
        enddo
c
c
c     Compute contribution to locexp coefficients.
c
        do idim=1,nd
          locexp(idim,0,0)= locexp(idim,0,0) + fhs(0)*charge(idim,isrc)
        enddo
        do n=1,nterms
          ztmp = ynm(n,0)*fhs(n)
          do idim=1,nd
            locexp(idim,n,0)= locexp(idim,n,0) + ztmp*charge(idim,isrc)
          enddo
          do m=1,n
            ztmp=ynm(n,m)*fhs(n)
            do idim=1,nd
              locexp(idim,n, m)= locexp(idim,n, m) + 
     1            ztmp*ephi(-m)*charge(idim,isrc)
              locexp(idim,n,-m)= locexp(idim,n,-m) + 
     1           ztmp*ephi(m)*charge(idim,isrc)
            enddo
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
C***********************************************************************
      subroutine h3dformtad(nd,zk,rscale,sources,
     1            dipvec,ns,center,nterms,locexp,wlege,nlege)
C***********************************************************************
C
C     Adds to multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of local expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(3,ns)    : dipole orientation vectors
C     ns              : number of sources
C     center(3)       : expansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     locexp          : coeffs of the h-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 dipvec(nd,3,ns)
      complex *16 locexp(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: ephi(:),fhs(:),fhder(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye,ztmp1
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 rx,thetax,phix,ry,thetay,phiy,rz,thetaz,phiz
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      integer idim
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder

      ifder=1

      allocate(fhs(0:nterms),fhder(0:nterms),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms),ynmd(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)

c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
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
        
        
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
c
c     get Hankel functions:
c
        z=zk*r
        call h3dall(nterms,z,rscale,fhs,ifder,fhder)
c
c
c     multiply all hn by  (i*k).
c
        do n = 0,nterms
          fhs(n) = fhs(n)*zkeye
          fhder(n) = fhder(n)*zkeye*zk
        enddo
c
c
c     Compute contribution to locexp coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|T|) Ylm*(T) h_n(k|S|)Ylm(S)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of already above
c

        ur = ynm(0,0)*fhder(0)
        do idim=1,nd
          zzz = ur*(dipvec(idim,1,isrc)*rx + dipvec(idim,2,isrc)*ry +
     1          dipvec(idim,3,isrc)*rz)
          locexp(idim,0,0)= locexp(idim,0,0)+zzz
        enddo
        do n=1,nterms
          ur = ynm(n,0)*fhder(n)
          utheta = -fhs(n)*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax
          uy = ur*ry + utheta*thetay
          uz = ur*rz + utheta*thetaz

          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy +
     1             dipvec(idim,3,isrc)*uz
            locexp(idim,n,0)= locexp(idim,n,0)+zzz   
          enddo
          do m=1,n
            ur = fhder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fhs(n)*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fhs(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              locexp(idim,n,m)= locexp(idim,n,m)+zzz
            enddo

            ur = fhder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fhs(n)*ynmd(n,m)
            uphi = eye*m*ephi(m)*fhs(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              locexp(idim,n,-m)= locexp(idim,n,-m)+
     1            zzz
            enddo
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
c
C***********************************************************************
      subroutine h3dformtacd(nd,zk,rscale,sources,charge,
     1            dipvec,ns,center,nterms,locexp,wlege,nlege)
C***********************************************************************
C
C     Adds to multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*) 
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of local expansions
C     zk              : Helmholtz parameter 
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : charge strengths
C     dipvec(3,ns)    : dipole orientation vectors
C     ns              : number of sources
C     center(3)       : expansion center
C     nterms          : order of multipole expansion
c     wlege           : precomputed array of scaling coeffs for Pnm
c     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     locexp          : coeffs of the h-expansion
c
c-----------------------------------------------------------------------
      implicit none
      integer ns,nterms,nlege,nd
      complex *16 zk
      real *8 rscale,sources(3,ns),center(3),wlege(*)
      complex *16 charge(nd,ns)
      complex *16 dipvec(nd,3,ns)
      complex *16 locexp(nd,0:nterms,-nterms:nterms)
c
cc      temporary variables
c
c
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16, allocatable :: ephi(:),fhs(:),fhder(:)
      complex *16 ephi1,ephi1inv
      complex *16 z,ztmp,zkeye,eye,ztmp1
      real *8 r,theta,phi
      real *8 ctheta,stheta,cphi,sphi
      real *8 rx,thetax,phix,ry,thetay,phiy,rz,thetaz,phiz
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      integer idim
      data eye/(0.0d0,1.0d0)/

      integer isrc,i,m,n,ifder

      ifder=1

      allocate(fhs(0:nterms),fhder(0:nterms),ephi(-nterms:nterms))
      allocate(ynm(0:nterms,0:nterms),ynmd(0:nterms,0:nterms))

      zkeye = eye*zk


      do isrc=1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)

        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)

c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
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
        
        
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        do i=2,nterms
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
c
c     get the associated Legendre functions
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
c
c     get Hankel functions:
c
        z=zk*r
        call h3dall(nterms,z,rscale,fhs,ifder,fhder)
c
c
c     multiply all hn by  (i*k).
c
        do n = 0,nterms
          fhs(n) = fhs(n)*zkeye
          fhder(n) = fhder(n)*zkeye*zk
        enddo
c
c
c     Compute contribution to locexp coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|T|) Ylm*(T) h_n(k|S|)Ylm(S)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of already above
c

        ur = ynm(0,0)*fhder(0)
        do idim=1,nd
          zzz = ur*(dipvec(idim,1,isrc)*rx + dipvec(idim,2,isrc)*ry +
     1          dipvec(idim,3,isrc)*rz)
          locexp(idim,0,0)= locexp(idim,0,0)+fhs(0)*charge(idim,isrc)+
     1       zzz
        enddo
        do n=1,nterms
          ur = ynm(n,0)*fhder(n)
          utheta = -fhs(n)*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax
          uy = ur*ry + utheta*thetay
          uz = ur*rz + utheta*thetaz

          ztmp = fhs(n)*ynm(n,0)
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy +
     1             dipvec(idim,3,isrc)*uz
            locexp(idim,n,0)= locexp(idim,n,0) + 
     1             ztmp*charge(idim,isrc) + zzz       
          enddo
          do m=1,n
            ztmp=ynm(n,m)*stheta*fhs(n)
            ztmp1 = ztmp*ephi(-m)

            ur = fhder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fhs(n)*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fhs(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              locexp(idim,n,m)= locexp(idim,n,m) + 
     1           ztmp1*charge(idim,isrc) + zzz
            enddo

            ztmp1 = ztmp*ephi(m)
            ur = fhder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fhs(n)*ynmd(n,m)
            uphi = eye*m*ephi(m)*fhs(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1             dipvec(idim,3,isrc)*uz
              locexp(idim,n,-m)= locexp(idim,n,-m) + 
     1             ztmp1*charge(idim,isrc) +  zzz
            enddo
          enddo
        enddo
      enddo

      return
      end
c
