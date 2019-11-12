c
c      this file contains the basic subroutines for 
c      forming and evaluating multipole expansions.
c
c      remarks on scaling conventions.
c
c      1)  far field and local expansions are consistently rscaled as
c              
c
c          m_n^m (scaled) = m_n^m / rscale^(n)  so that upon evaluation
c
c          the field is  sum   m_n^m (scaled) * rscale^(n) / r^{n+1}.
c
c          l_n^m (scaled) = l_n^m * rscale^(n)  so that upon evaluation
c
c          the field is  sum   l_n^m (scaled) / rscale^(n) * r^{n}.
c
c
c      2) there are many definitions of the spherical harmonics,
c         which differ in terms of normalization constants. we
c         adopt the following convention:
c
c         for m>0, we define y_n^m according to 
c
c         y_n^m = \sqrt{2n+1} \sqrt{\frac{ (n-m)!}{(n+m)!}} \cdot
c                 p_n^m(\cos \theta)  e^{i m phi} 
c         and
c 
c         y_n^-m = dconjg( y_n^m )
c    
c         we omit the condon-shortley phase factor (-1)^m in the 
c         definition of y_n^m for m<0. (this is standard in several
c         communities.)
c
c         we also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  (this is also standard in several communities.)
c         more precisely, 
c
c                 \int_s y_n^m y_n^m d\omega = 4 \pi. 
c
c         using our standard definition, the addition theorem takes 
c         the simple form 
c
c         1/r = 
c         \sum_n 1/(2n+1) \sum_m  |s|^n ylm*(s) ylm(t)/ (|t|^(n+1)) 
c
c         1/r = 
c         \sum_n \sum_m  |s|^n  ylm*(s)    ylm(t)     / (|t|^(n+1)) 
c                               -------    ------
c                               sqrt(2n+1) sqrt(2n+1)
c
c        in the laplace library (this library), we incorporate the
c        sqrt(2n+1) factor in both forming and evaluating multipole
c        expansions.
c
c-----------------------------------------------------------------------
c
c      l3dmpevalp: computes potential due to a multipole expansion
c                    at a collection of targets (done,tested)
c
c      l3dmpevalg: computes potential and gradients 
c                  due to a multipole expansion
c                    at a collection of targets (done,tested)
c
c      l3dformmpc: creates multipole expansion (outgoing) due to 
c                 a collection of charges (done,tested)
c
c      l3dformmpd: creates multipole expansion (outgoing) due to 
c                 a collection of dipoles (done,tested)
c
c      l3dformmpcd: creates multipole expansion (outgoing) due to 
c                 a collection of charges and dipoles (done,tested)
c
c      l3dtaevalp: computes potential 
c                  due to local expansion at a collection of targets
c
c      l3dtaevalp: computes potential and gradients
c                  due to local expansion at a collection of targets
c
c      l3dformtac: creates local expansion due to 
c                 a collection of charges.
c
c      l3dformtad: creates local expansion due to 
c                 a collection of dipoles
c
c      l3dformtacd: creates local expansion due to 
c                 a collection of charges and dipoles
c
c**********************************************************************
      subroutine l3dmpevalp(nd,rscale,center,mpole,nterms,
     1		ztarg,ntarg,pot,wlege,nlege,thresh)
c**********************************************************************
c
c
c     this subroutine evaluates the potential   
c     of an outgoing multipole expansioni and adds
c     to existing quantities
c
c     pot =  pot + sum sum  mpole(n,m) Y_nm(theta,phi) / r^{n+1} 
c                   n   m
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of target locations
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c     thresh :    threshold for computing outgoing expansion,
c                 potential at target location
c                 won't be updated if |t-c| <= thresh, where
c                 t is the target location and c is the expansion
c                 center location
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c
c----------------------------------------------------------------------
      implicit none

c
cc     calling sequence variables
c

      integer nterms,nlege,ntarg,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      real *8 pot(nd,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(0:nlege,0:nlege), thresh

c
cc     temporary variables
c

      integer idim
      real *8, allocatable :: ynm(:,:),fr(:)
      complex *16, allocatable :: ephi(:)
      integer i,j,k,l,m,n,itarg
      real *8 done,r,theta,phi,zdiff(3)
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rs,rtmp1,rtmp2
      complex *16 ephi1
c
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0

      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1))
      allocate(ynm(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)

        if(abs(r).lt.thresh) goto 1000 

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
        d = 1.0d0/r
        fr(0) = d
        d = d*rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
        enddo
c
c    get the associated Legendre functions:
c

        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
        do l = 0,nterms
          rs = sqrt(1.0d0/(2*l+1))
          do m=0,l
            ynm(l,m) = ynm(l,m)*rs
          enddo
        enddo

        do idim=1,nd
          pot(idim,itarg) = pot(idim,itarg) +
     1                 real(mpole(idim,0,0))*fr(0)
        enddo
        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+real(mpole(idim,n,0))*rtmp1
          enddo
	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)
            do idim=1,nd
              rtmp2 = 2*real(mpole(idim,n,m)*ephi(m)) 
              pot(idim,itarg)=pot(idim,itarg)+rtmp1*rtmp2
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
c**********************************************************************
      subroutine l3dmpevalg(nd,rscale,center,mpole,nterms,
     1		ztarg,ntarg,pot,grad,wlege,nlege,thresh)
c**********************************************************************
c
c
c     this subroutine evaluates the potential and gradient of  
c     of an outgoing multipole expansioni and adds
c     to existing quantities
c
c     pot =  pot + sum sum  mpole(n,m) Y_nm(theta,phi) / r^{n+1} 
c                   n   m
c
c     grad =  grad + Gradient( sum sum  mpole(n,m) Y_nm(theta,phi)/r^{n+1})
c                               n   m
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of target locations
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c     thresh :    threshold for computing outgoing expansion,
c                 potential and gradient at target location
c                 won't be updated if |t-c| <= thresh, where
c                 t is the target location and c is the expansion
c                 center location
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c     grad   :   updated gradient at ztarg 
c
c----------------------------------------------------------------------
      implicit none

c
cc     calling sequence variables
c

      integer nterms,nlege,ntarg,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(0:nlege,0:nlege), thresh

c
cc     temporary variables
c
      integer idim
      real *8, allocatable :: ynm(:,:),ynmd(:,:),fr(:),frder(:)
      complex *16, allocatable :: ephi(:)
      integer i,j,k,l,m,n,itarg
      real *8 done,r,theta,phi,zdiff(3)
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      real *8 rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6
      complex *16 ephi1
      real *8 ur(nd),utheta(nd),uphi(nd)
c
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0

      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1),frder(0:nterms))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)

        if(abs(r).lt.thresh) goto 1000 

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
        d = 1.0d0/r
        fr(0) = d
        d = d*rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
        enddo
        do i=0,nterms
          frder(i) = -(i+1.0d0)*fr(i+1)/rscale
        enddo
c
c    get the associated Legendre functions:
c

        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do l = 0,nterms
          rs = sqrt(1.0d0/(2*l+1))
          do m=0,l
            ynm(l,m) = ynm(l,m)*rs
            ynmd(l,m) = ynmd(l,m)*rs
          enddo
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

        do idim=1,nd
          ur(idim) = real(mpole(idim,0,0))*frder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
          pot(idim,itarg) = pot(idim,itarg) + 
     1        real(mpole(idim,0,0))*fr(0)
        enddo

        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          rtmp2 = frder(n)*ynm(n,0)
          rtmp3 = -fr(n)*ynmd(n,0)*stheta
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+real(mpole(idim,n,0))*rtmp1
            ur(idim)=ur(idim)+real(mpole(idim,n,0))*rtmp2
            utheta(idim)=utheta(idim)+real(mpole(idim,n,0))*rtmp3
          enddo

	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)*stheta
            rtmp4 = frder(n)*ynm(n,m)*stheta
            rtmp5 = -fr(n)*ynmd(n,m)
            rtmp6 = -m*fr(n)*ynm(n,m)

            do idim=1,nd
              rtmp2 = 2*real(mpole(idim,n,m)*ephi(m)) 

              pot(idim,itarg)=pot(idim,itarg)+rtmp1*rtmp2
              ur(idim) = ur(idim) + rtmp4*rtmp2
              utheta(idim) = utheta(idim)+rtmp5*rtmp2
              rtmp2 = 2*imag(mpole(idim,n,m)*ephi(m))
              uphi(idim) = uphi(idim) + rtmp6*rtmp2
            enddo
          enddo
        enddo

        do idim=1,nd
          grad(idim,1,itarg)=grad(idim,1,itarg)+ur(idim)*rx+
     1          utheta(idim)*thetax+uphi(idim)*phix
          grad(idim,2,itarg)=grad(idim,2,itarg)+ur(idim)*ry+
     1          utheta(idim)*thetay+uphi(idim)*phiy
          grad(idim,3,itarg)=grad(idim,3,itarg)+ur(idim)*rz+
     1          utheta(idim)*thetaz+uphi(idim)*phiz
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
      subroutine l3dformmpc(nd,rscale,sources,charge,ns,center,
     1                  nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS charges 
C     located at SOURCES(3,*) and add to existing expansions
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(nd,ns)   : charge strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 charge(nd,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:)
      real *8 theta,stheta,ctheta,phi,sphi,cphi,dtmp,d,r
      complex *16, allocatable :: ephi(:)
      complex *16 ephi1
      complex *16 eye
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))

      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        fr(0) = 1.0d0
        d = r/rscale
        fr(1) = d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo

c
c     get the associated Legendre functions and rescale
c      by 1/sqrt(2*l+1)
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
          enddo
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
        do idim=1,nd
          mpole(idim,0,0)= mpole(idim,0,0) + fr(0)*charge(idim,isrc)
        enddo
        do n=1,nterms
          dtmp=ynm(n,0)*fr(n)
          do idim=1,nd
            mpole(idim,n,0)= mpole(idim,n,0) + dtmp*charge(idim,isrc)
          enddo
          do m=1,n
            dtmp=ynm(n,m)*fr(n)
            do idim=1,nd
              mpole(idim,n,m) = mpole(idim,n,m) + 
     1                  dtmp*ephi(-m)*charge(idim,isrc)
              mpole(idim,n,-m) = mpole(idim,n,-m) + 
     1                  dtmp*ephi(m)*charge(idim,isrc)
            enddo
          enddo
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
c
C***********************************************************************
      subroutine l3dformmpd(nd,rscale,sources,dipvec,ns,center,
     1                  nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS dipoles 
C     located at SOURCES(3,*) and adds to existing expansion
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(nd,3,ns) : dipole orientiation vectors
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 dipvec(nd,3,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:),frder(:),ynmd(:,:)
      real *8 thetaz,thetay,thetax, theta
      real *8 stheta,sphi,rx,ry,rz,r
      real *8 ctheta,cphi
      real *8 phix,phiy,phiz,phi,fruse,d

      complex *16 ur,utheta,uphi,ux,uy,uz,zzz
      complex *16, allocatable :: ephi(:)
      complex *16 eye,ephi1
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(frder(0:nterms),ynmd(0:nterms,0:nterms))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))


      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        fr(0) = 1.0d0
        d = r/rscale
        fr(1) = d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
        frder(0) = 0.0d0
        do i=1,nterms
          frder(i) = i*fr(i-1)/rscale
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
        thetax = ctheta*cphi
        phix = -sphi
        ry = stheta*sphi
        thetay = ctheta*sphi
        phiy = cphi
        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0
c
c     get the associated Legendre functions and rescale by
c       1/sqrt(2*l+1)
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
            ynmd(j,i) = ynmd(j,i)*rfac(j)
          enddo
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
        ur = ynm(0,0)*frder(0)
        ux = ur*rx 
        uy = ur*ry 
        uz = ur*rz
        do idim=1,nd
          zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
          mpole(idim,0,0)= mpole(idim,0,0) + zzz
        enddo

        do n=1,nterms
          fruse = fr(n-1)/rscale
          ur = ynm(n,0)*frder(n)
          utheta = -fruse*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax 
          uy = ur*ry + utheta*thetay 
          uz = ur*rz + utheta*thetaz
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + zzz
          enddo
          do m=1,n
            ur = frder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fruse*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fruse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,m)= mpole(idim,n,m) + zzz
            enddo
c
            ur = frder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fruse*ynmd(n,m)
            uphi = eye*m*ephi(m)*fruse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m) + zzz
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
c
c
c
c
C***********************************************************************
      subroutine l3dformmpcd(nd,rscale,sources,charge,dipvec,ns,
     1             center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS charges+dipoles 
C     located at SOURCES(3,*) and adds to existing expansion
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(nd,ns)   : charge strengths
C     dipvec(nd,3,ns) : dipole orientiation vectors
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 charge(nd,ns)
      real *8 dipvec(nd,3,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:),frder(:),ynmd(:,:)
      real *8 thetaz,thetay,thetax, theta
      real *8 stheta,sphi,rx,ry,rz,r
      real *8 ctheta,cphi
      real *8 phix,phiy,phiz,phi,fruse,d,dtmp

      complex *16 ur,utheta,uphi,ux,uy,uz,zzz
      complex *16, allocatable :: ephi(:)
      complex *16 eye,ephi1
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(frder(0:nterms),ynmd(0:nterms,0:nterms))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))

      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        fr(0) = 1.0d0
        d = r/rscale
        fr(1) = d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
        frder(0) = 0.0d0
        do i=1,nterms
          frder(i) = i*fr(i-1)/rscale
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
        thetax = ctheta*cphi
        phix = -sphi
        ry = stheta*sphi
        thetay = ctheta*sphi
        phiy = cphi
        rz = ctheta
        thetaz = -stheta
        phiz = 0.0d0
c
c     get the associated Legendre functions and rescale by
c       1/sqrt(2*l+1)
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
            ynmd(j,i) = ynmd(j,i)*rfac(j)
          enddo
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
        ur = ynm(0,0)*frder(0)
        ux = ur*rx 
        uy = ur*ry 
        uz = ur*rz
        do idim=1,nd
          zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
          mpole(idim,0,0)= mpole(idim,0,0) + zzz +
     1            fr(0)*charge(idim,isrc)
        enddo

        do n=1,nterms
          fruse = fr(n-1)/rscale
          ur = ynm(n,0)*frder(n)
          utheta = -fruse*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax 
          uy = ur*ry + utheta*thetay 
          uz = ur*rz + utheta*thetaz
          dtmp = fr(n)*ynm(n,0)
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + zzz + 
     1         charge(idim,isrc)*dtmp
          enddo
          do m=1,n
            ur = frder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fruse*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fruse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            dtmp = ynm(n,m)*fr(n)*stheta
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,m)= mpole(idim,n,m) + zzz + 
     1            charge(idim,isrc)*dtmp*ephi(-m)
            enddo
c
            ur = frder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fruse*ynmd(n,m)
            uphi = eye*m*ephi(m)*fruse*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m)+zzz+
     1              charge(idim,isrc)*dtmp*ephi(m)
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
c
c
c
c**********************************************************************
      subroutine l3dtaevalp(nd,rscale,center,mpole,nterms,
     1		ztarg,ntarg,pot,wlege,nlege)
c**********************************************************************
c
c
c     this subroutine evaluates the potential   
c     of an incoming multipole expansioni and adds
c     to existing quantities
c
c     pot =  pot + sum sum  r^{n} mpole(n,m) Y_nm(theta,phi) /sqrt(2n+1) 
c                   n   m
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    local expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of target locations
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c
c----------------------------------------------------------------------
      implicit none

c
cc     calling sequence variables
c

      integer nterms,nlege,ntarg,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      real *8 pot(nd,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(0:nlege,0:nlege), thresh

c
cc     temporary variables
c

      integer idim
      real *8, allocatable :: ynm(:,:),fr(:)
      complex *16, allocatable :: ephi(:)
      integer i,j,k,l,m,n,itarg
      real *8 done,r,theta,phi,zdiff(3)
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rs,rtmp1,rtmp2
      complex *16 ephi1
c
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0

      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1))
      allocate(ynm(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)
c
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
        fr(0) = 1.0d0
        d = r/rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
        enddo
c
c    get the associated Legendre functions:
c

        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
        do l = 0,nterms
          rs = sqrt(1.0d0/(2*l+1))
          do m=0,l
            ynm(l,m) = ynm(l,m)*rs
          enddo
        enddo

        do idim=1,nd
          pot(idim,itarg) = pot(idim,itarg) + 
     1        real(mpole(idim,0,0))*fr(0)
        enddo
        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+real(mpole(idim,n,0))*rtmp1
          enddo
	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)
            do idim=1,nd
              rtmp2 = 2*real(mpole(idim,n,m)*ephi(m)) 
              pot(idim,itarg)=pot(idim,itarg)+rtmp1*rtmp2
            enddo
          enddo
        enddo
      enddo

      return
      end
c
c
c
c**********************************************************************
      subroutine l3dtaevalg(nd,rscale,center,mpole,nterms,
     1		ztarg,ntarg,pot,grad,wlege,nlege)
c**********************************************************************
c
c
c     this subroutine evaluates the potential and gradient of  
c     of an local multipole expansioni and adds
c     to existing quantities
c
c     pot =  pot + sum sum  mpole(n,m) r^{n}Y_nm(theta,phi) / sqrt(2n+1) 
c                   n   m
c
c     grad =  grad + Gradient( sum sum mpole(n,m)r^{n}Y_nm(theta,phi)/sqrt(2n+1))
c                               n   m
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    local expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of target locations
c     wlege  :    precomputed array of scaling coeffs for Pnm
c     nlege  :    dimension parameter for wlege
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c     grad   :   updated gradient at ztarg 
c
c----------------------------------------------------------------------
      implicit none

c
cc     calling sequence variables
c

      integer nterms,nlege,ntarg,nd
      real *8 rscale,center(3),ztarg(3,ntarg)
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 wlege(0:nlege,0:nlege)

c
cc     temporary variables
c
      integer idim
      real *8, allocatable :: ynm(:,:),ynmd(:,:),fr(:),frder(:)
      complex *16, allocatable :: ephi(:)
      integer i,j,k,l,m,n,itarg
      real *8 done,r,theta,phi,zdiff(3)
      real *8 ctheta,stheta,cphi,sphi
      real *8 d,rx,ry,rz,thetax,thetay,thetaz,phix,phiy,phiz,rs
      real *8 rtmp1,rtmp2,rtmp3,rtmp4,rtmp5,rtmp6
      complex *16 ephi1
      real *8 ur(nd),utheta(nd),uphi(nd)
c
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0

      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1),frder(0:nterms))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))

      do itarg=1,ntarg
        zdiff(1)=ztarg(1,itarg)-center(1)
        zdiff(2)=ztarg(2,itarg)-center(2)
        zdiff(3)=ztarg(3,itarg)-center(3)
c
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
        d = r/rscale
        fr(0) = 1.0d0
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
        enddo
        frder(0) = 0
        do i=1,nterms
          frder(i) = i*fr(i-1)/rscale
        enddo
c
c    get the associated Legendre functions:
c

        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do l = 0,nterms
          rs = sqrt(1.0d0/(2*l+1))
          do m=0,l
            ynm(l,m) = ynm(l,m)*rs
            ynmd(l,m) = ynmd(l,m)*rs
          enddo
        enddo

c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c
c     NOTE: sphereical derivative needs to be fixed for r=0
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

        do idim=1,nd
          ur(idim) = real(mpole(idim,0,0))*frder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
          pot(idim,itarg) = pot(idim,itarg)+real(mpole(idim,0,0))*fr(0)
        enddo

        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          rtmp2 = frder(n)*ynm(n,0)
          rtmp3 = -fr(n)*ynmd(n,0)*stheta
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+real(mpole(idim,n,0))*rtmp1
            ur(idim)=ur(idim)+real(mpole(idim,n,0))*rtmp2
            utheta(idim)=utheta(idim)+real(mpole(idim,n,0))*rtmp3
          enddo

	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)*stheta
            rtmp4 = frder(n)*ynm(n,m)*stheta
            rtmp5 = -fr(n)*ynmd(n,m)
            rtmp6 = -m*fr(n)*ynm(n,m)

            do idim=1,nd
              rtmp2 = 2*real(mpole(idim,n,m)*ephi(m)) 

              pot(idim,itarg)=pot(idim,itarg)+rtmp1*rtmp2
              ur(idim) = ur(idim) + rtmp4*rtmp2
              utheta(idim) = utheta(idim)+rtmp5*rtmp2
              rtmp2 = 2*imag(mpole(idim,n,m)*ephi(m))
              uphi(idim) = uphi(idim) + rtmp6*rtmp2
            enddo
          enddo
        enddo

        do idim=1,nd
          grad(idim,1,itarg)=grad(idim,1,itarg)+ur(idim)*rx+
     1          utheta(idim)*thetax+uphi(idim)*phix
          grad(idim,2,itarg)=grad(idim,2,itarg)+ur(idim)*ry+
     1          utheta(idim)*thetay+uphi(idim)*phiy
          grad(idim,3,itarg)=grad(idim,3,itarg)+ur(idim)*rz+
     1          utheta(idim)*thetaz+uphi(idim)*phiz
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
      subroutine l3dformtac(nd,rscale,sources,charge,ns,center,
     1                  nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs local expansion about CENTER due to NS charges 
C     located at SOURCES(3,*) and add to existing expansions
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(nd,ns)   : charge strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the local expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 charge(nd,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:)
      real *8 theta,stheta,ctheta,phi,sphi,cphi,dtmp,d,r
      complex *16, allocatable :: ephi(:)
      complex *16 ephi1
      complex *16 eye
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))

      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        d = 1.0d0/r
        fr(0) = d
        d = d*rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo

c
c     get the associated Legendre functions and rescale
c      by 1/sqrt(2*l+1)
c
        call ylgndrfw(nterms,ctheta,ynm,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
          enddo
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
        do idim=1,nd
          mpole(idim,0,0)= mpole(idim,0,0) + fr(0)*charge(idim,isrc)
        enddo
        do n=1,nterms
          dtmp=ynm(n,0)*fr(n)
          do idim=1,nd
            mpole(idim,n,0)= mpole(idim,n,0) + dtmp*charge(idim,isrc)
          enddo
          do m=1,n
            dtmp=ynm(n,m)*fr(n)
            do idim=1,nd
              mpole(idim,n,m) = mpole(idim,n,m) + 
     1                  dtmp*ephi(-m)*charge(idim,isrc)
              mpole(idim,n,-m) = mpole(idim,n,-m) + 
     1                  dtmp*ephi(m)*charge(idim,isrc)
            enddo
          enddo
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
c
C***********************************************************************
      subroutine l3dformtad(nd,rscale,sources,dipvec,ns,center,
     1                  nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS dipoles 
C     located at SOURCES(3,*) and adds to existing expansion
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipvec(nd,3,ns) : dipole orientiation vectors
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 dipvec(nd,3,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:),frder(:),ynmd(:,:)
      real *8 thetaz,thetay,thetax, theta
      real *8 stheta,sphi,rx,ry,rz,r
      real *8 ctheta,cphi
      real *8 phix,phiy,phiz,phi,d

      complex *16 ur,utheta,uphi,ux,uy,uz,zzz
      complex *16, allocatable :: ephi(:)
      complex *16 eye,ephi1
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(frder(0:nterms),ynmd(0:nterms,0:nterms))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))

      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        d = 1.0d0/r
        fr(0) = d
        d = d*rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
        do i=0,nterms
          frder(i) = -(i+1.0d0)*fr(i+1)/rscale
        enddo
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
c
c     get the associated Legendre functions and rescale by
c       1/sqrt(2*l+1)
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
            ynmd(j,i) = ynmd(j,i)*rfac(j)
          enddo
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
        ur = ynm(0,0)*frder(0)
        ux = ur*rx 
        uy = ur*ry 
        uz = ur*rz
        do idim=1,nd
          zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
          mpole(idim,0,0)= mpole(idim,0,0) + zzz
        enddo

        do n=1,nterms
          ur = ynm(n,0)*frder(n)
          utheta = -fr(n)*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax 
          uy = ur*ry + utheta*thetay 
          uz = ur*rz + utheta*thetaz
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + zzz
          enddo
          do m=1,n
            ur = frder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fr(n)*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fr(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,m)= mpole(idim,n,m) + zzz
            enddo
c
            ur = frder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fr(n)*ynmd(n,m)
            uphi = eye*m*ephi(m)*fr(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m) + zzz
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
c
c
c
c
C***********************************************************************
      subroutine l3dformtacd(nd,rscale,sources,charge,dipvec,ns,
     1             center,nterms,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole expansion about CENTER due to NS charges+dipoles 
C     located at SOURCES(3,*) and adds to existing expansion
C
c-----------------------------------------------------------------------
C     INPUT:
c
c     nd              : number of multipole expansions
C     rscale          : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(nd,ns)   : charge strengths
C     dipvec(nd,3,ns) : dipole orientiation vectors
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     wlege           : precomputed array of scaling coeffs for pnm
C     nlege           : dimension parameter for wlege
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
c
cc       calling sequence variables
c
      
      integer nterms,ns,nd, nlege
      real *8 center(3),sources(3,ns)
      real *8 wlege(0:nlege,0:nlege)
      real *8 rscale
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
      real *8 charge(nd,ns)
      real *8 dipvec(nd,3,ns)

c
cc       temporary variables
c

      integer i,j,k,l,m,n,isrc,idim
      real *8 zdiff(3)
      real *8, allocatable :: ynm(:,:),fr(:),rfac(:),frder(:),ynmd(:,:)
      real *8 thetaz,thetay,thetax, theta
      real *8 stheta,sphi,rx,ry,rz,r
      real *8 ctheta,cphi
      real *8 phix,phiy,phiz,phi,fruse,d,dtmp

      complex *16 ur,utheta,uphi,ux,uy,uz,zzz
      complex *16, allocatable :: ephi(:)
      complex *16 eye,ephi1
      data eye/(0.0d0,1.0d0)/

      allocate(ynm(0:nterms,0:nterms),fr(0:nterms+1))
      allocate(frder(0:nterms),ynmd(0:nterms,0:nterms))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(rfac(0:nterms))

      do i=0,nterms
        rfac(i) = 1/sqrt(2.0d0*i + 1.0d0)
      enddo

      do isrc = 1,ns
        zdiff(1)=sources(1,isrc)-center(1)
        zdiff(2)=sources(2,isrc)-center(2)
        zdiff(3)=sources(3,isrc)-center(3)
c
        call cart2polar(zdiff,r,theta,phi)
        ctheta = dcos(theta)
        stheta = dsin(theta)
        cphi = dcos(phi)
        sphi = dsin(phi)
        ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array and fr array
c
        ephi(0)=1.0d0
        ephi(1)=ephi1
        ephi(-1)=dconjg(ephi1)
        d = 1.0d0/r
        fr(0) = d
        d = d*rscale
        fr(1) = fr(0)*d
        do i=2,nterms+1
          fr(i) = fr(i-1)*d
          ephi(i)=ephi(i-1)*ephi1
          ephi(-i)=ephi(-i+1)*ephi(-1)
        enddo
        do i=0,nterms
          frder(i) = -(i+1.0d0)*fr(i+1)/rscale
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
        thetax = ctheta*cphi/r
        phix = -sphi/r
        ry = stheta*sphi
        thetay = ctheta*sphi/r
        phiy = cphi/r
        rz = ctheta
        thetaz = -stheta/r
        phiz = 0.0d0
c
c     get the associated Legendre functions and rescale by
c       1/sqrt(2*l+1)
c
        call ylgndr2sfw(nterms,ctheta,ynm,ynmd,wlege,nlege)
        do i=0,nterms
          do j=0,nterms
            ynm(j,i) = ynm(j,i)*rfac(j)
            ynmd(j,i) = ynmd(j,i)*rfac(j)
          enddo
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
        ur = ynm(0,0)*frder(0)
        ux = ur*rx 
        uy = ur*ry 
        uz = ur*rz
        do idim=1,nd
          zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
          mpole(idim,0,0)= mpole(idim,0,0) + zzz + 
     1            fr(0)*charge(idim,isrc)
        enddo

        do n=1,nterms
          ur = ynm(n,0)*frder(n)
          utheta = -fr(n)*ynmd(n,0)*stheta
          ux = ur*rx + utheta*thetax 
          uy = ur*ry + utheta*thetay 
          uz = ur*rz + utheta*thetaz
          dtmp = fr(n)*ynm(n,0)
          do idim=1,nd
            zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1        dipvec(idim,3,isrc)*uz
            mpole(idim,n,0)= mpole(idim,n,0) + zzz + 
     1         charge(idim,isrc)*dtmp
          enddo
          do m=1,n
            ur = frder(n)*ynm(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fr(n)*ynmd(n,m)
            uphi = -eye*m*ephi(-m)*fr(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            dtmp = ynm(n,m)*fr(n)*stheta
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,m)= mpole(idim,n,m) + zzz+ 
     1            charge(idim,isrc)*dtmp*ephi(-m)
            enddo
c
            ur = frder(n)*ynm(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fr(n)*ynmd(n,m)
            uphi = eye*m*ephi(m)*fr(n)*ynm(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            do idim=1,nd
              zzz = dipvec(idim,1,isrc)*ux + dipvec(idim,2,isrc)*uy + 
     1          dipvec(idim,3,isrc)*uz
              mpole(idim,n,-m)= mpole(idim,n,-m)+zzz+
     1              charge(idim,isrc)*dtmp*ephi(m)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
