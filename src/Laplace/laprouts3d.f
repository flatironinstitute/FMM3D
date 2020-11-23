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
c      l3dmpevalp: computes potentials due to a multipole expansion
c                    at a collection of targets (done,tested)
c
c      l3dmpevalg: computes potentials and gradients 
c                  due to a multipole expansion
c                  at a collection of targets (done,tested)
c
c      l3dmpevalh: computes potentials, gradients, and hessians 
c                  due to a multipole expansion
c                  at a collection of targets (done,tested)
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
c      l3dtaevalp: computes potentials 
c                  due to local expansion at a collection of targets
c
c      l3dtaevalg: computes potentials and gradients
c                  due to local expansion at a collection of targets
c
c      l3dtaevalh: computes potentials, gradients, and hessians
c                  due to local expansion at a collection of targets
c
c      l3dtaevalhessdini: initialization routine for l3dtaevalh
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
c
c      l3dmpevalhessdini:  initialization routine for l3dmpevalhessd
c
c      l3dmpevalh: computes potentials, gradients and Hessians
c                  due to a multipole expansion
c                    at a collection of targets (done,tested)
c
c**********************************************************************
      subroutine l3dmpevalp(nd,rscale,center,mpole,nterms,
     1		ztarg,ntarg,pot,wlege,nlege,thresh)
c**********************************************************************
c
c     this subroutine evaluates the potentials due to an
c     outgoing multipole expansion and increments inputs accordingly:
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
c     ztarg  :    target locations
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
c     pot    :   updated potentials at all targets
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
c     this subroutine evaluates the potentials and gradients due to  
c     an outgoing multipole expansion and increments inputs accordingly:
c
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
c     pot    :   updated potentials at targets
c     grad   :   updated gradients at targets 
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
c     this subroutine evaluates the potentials due to an   
c     incoming local expansion and increments accordingly:
c
c     pot =  pot + sum sum  r^{n} mpole(n,m) Y_nm(theta,phi) /sqrt(2n+1) 
c                   n   m
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
c     pot    :   updated potentials at targets
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
c     this subroutine evaluates the potentials and gradients due to  
c     an incoming local expansion and increments inputs accordingly:
c
c     pot =  pot + sum sum  mpole(n,m) r^{n}Y_nm(theta,phi) / sqrt(2n+1) 
c                   n   m
c
c     grad =  grad + 
c              Gradient( sum sum mpole(n,m)r^{n}Y_nm(theta,phi)/sqrt(2n+1))
c                         n   m
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
c     pot    :   updated potentials at targets
c     grad   :   updated gradients at targets 
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
      real *8 rscaleinv
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
      allocate(ephi(0:nterms+1))
      allocate(fr(0:nterms+1),frder(0:nterms))
      allocate(ynm(0:nterms,0:nterms))
      allocate(ynmd(0:nterms,0:nterms))
c
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

        rscaleinv = 1.0d0/rscale
        rx = stheta*cphi
        thetax = ctheta*cphi*rscaleinv
        phix = -sphi*rscaleinv
        ry = stheta*sphi
        thetay = ctheta*sphi*rscaleinv
        phiy = cphi*rscaleinv
        rz = ctheta
        thetaz = -stheta*rscaleinv
        phiz = 0.0d0
c
        do idim=1,nd
          ur(idim) = real(mpole(idim,0,0))*frder(0)
          utheta(idim) = 0.0d0
          uphi(idim) = 0.0d0
          pot(idim,itarg) = pot(idim,itarg)+real(mpole(idim,0,0))*fr(0)
        enddo
c
        do n=1,nterms
          rtmp1 = fr(n)*ynm(n,0)
          rtmp2 = frder(n)*ynm(n,0)
          rtmp3 = -fr(n-1)*ynmd(n,0)*stheta
          do idim=1,nd
            pot(idim,itarg)=pot(idim,itarg)+real(mpole(idim,n,0))*rtmp1
            ur(idim)=ur(idim)+real(mpole(idim,n,0))*rtmp2
            utheta(idim)=utheta(idim)+real(mpole(idim,n,0))*rtmp3
          enddo
c
	      do m=1,n
            rtmp1 = fr(n)*ynm(n,m)*stheta
            rtmp4 = frder(n)*ynm(n,m)*stheta
            rtmp5 = -fr(n-1)*ynmd(n,m)
            rtmp6 = -m*fr(n-1)*ynm(n,m)

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
c
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
c
C***********************************************************************
      subroutine l3dmpevalhessdini(nterms,scarray)
C***********************************************************************
C
c     Precomputes array used in 
c     mpole-local translation operator from an nterms expansion to an 
c     order 2 expansion (sufficient to compute pot/fld/hessian).     
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     scarray         : work array of size > 10*(nterms+2)**2
c-----------------------------------------------------------------------
      implicit none
      integer  nterms,l,j,k,m,ll,mm,iuse,lnew,mnew
      real *8 scarray(1),cscale
      real *8 d
      real *8, allocatable :: c(:,:)
      real *8, allocatable :: sqc(:,:)
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
c
c
C***********************************************************************
      subroutine l3dmpevalh(nd,rscale,center,mpole,nterms,
     1            ztarg,ntarg,pot,grad,hess,thresh,scarray)
C***********************************************************************
c
c     This subroutine evaluates the potential, gradient and
c     Hessian of the potential due to a multipole expansion and adds
c     to existing quantities
c
c     pot = pot + sum sum  mpole(n,m) Y_nm(theta,phi)  / r^{n+1}
c             n   m
c
c     grad =  grad + Gradient( sum sum  mpole(n,m) Y_nm(theta,phi)/r^{n+1})
c                               n   m
c
c     hess =  hess + Hessian( sum sum  mpole(n,m) Y_nm(theta,phi)/r^{n+1})
c                              n   m
c
c     The method is based on translation of mpole to 
c     second order expansion at target location.
c     It is reasonably optimized, precomputing the array of 
c     binomial/factorial terms that appear in the shift operator.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of multipole expansions
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of  target location
c     thresh :    threshold for computing outgoing expansion,
c                 potential and gradient at target location
c                 won't be updated if |t-c| <= thresh, where
c                 t is the target location and c is the expansion
c                 center location
c     scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                 L3DMPEVALHESSDINI(nterms,scarray))
c                 with dimension of scarray at least 10*(nterms+2)**2
c                 If nterms is changed, 
c                 l3dtaevalhessdini must be called again.
c
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c     grad   :   updated gradient at ztarg 
c     hess   :   updated Hessian at ztarg
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,ntarg,nd,itarg
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,lsum,idim
      real *8 center(3),ztarg(3,ntarg)
      real *8 zdiff(3)
      real *8 scarray(*),rscale,thresh
      real *8 cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      complex *16 mpole(nd,0:nterms,-nterms:nterms)
ccc      complex *16 local2(0:2,-2:2)
      complex *16 z0,ima,ephi1
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 hess(nd,6,ntarg)
c
      real *8, allocatable :: pp(:,:)
      real *8, allocatable :: powers(:)
      complex *16, allocatable :: local2(:,:)
      complex *16, allocatable :: ppc(:,:)
      complex *16, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms+2,0:nterms+2))
      allocate(ppc(0:nterms+2,-nterms-2:nterms+2))
      allocate(powers(0:nterms+3))
      allocate(ephi(-nterms-3:nterms+3))
      allocate(local2(nd,9))
c
c     determine order of shifted expansion 
c
c
      do itarg = 1,ntarg
c
         do ll = 1,nd
         do l = 1,9
            local2(ll,l) = 0.0d0
         enddo
         enddo
c
         zdiff(1) = center(1) - ztarg(1,itarg)
         zdiff(2) = center(2) - ztarg(2,itarg)
         zdiff(3) = center(3) - ztarg(3,itarg)
         call cart2polar(zdiff,d,theta,phi)
c
         if (abs(d).lt.thresh) goto 1000
c
         ctheta = dcos(theta)
         cphi = dcos(phi)
         sphi = dsin(phi)
         ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
         dd = 1.0d0/d
         dd = dd*rscale
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
               do idim=1,nd
                  local2(idim,1) = local2(idim,1) +
     1            ppc(l,-m)*scarray(iuse)*mpole(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
         enddo
         do l = 0,nterms
            lsum = l+1
            do m = -l,l
               do idim=1,nd
                  local2(idim,2) = local2(idim,2) +
     1            ppc(lsum,-1-m)*scarray(iuse)*mpole(idim,l,m)
                  local2(idim,3) = local2(idim,3) +
     1            ppc(lsum,  -m)*scarray(iuse+1)*mpole(idim,l,m)
                  local2(idim,4) = local2(idim,4) +
     1            ppc(lsum, 1-m)*scarray(iuse+2)*mpole(idim,l,m)
               enddo
               iuse = iuse+3
            enddo
         enddo
         do l = 0,nterms
            lsum = l+2
            do m = -l,l
               do idim=1,nd
                  local2(idim,5) = local2(idim,5) +
     1            ppc(lsum,-2-m)*scarray(iuse  )*mpole(idim,l,m)
                  local2(idim,6) = local2(idim,6) +
     1            ppc(lsum,-1-m)*scarray(iuse+1)*mpole(idim,l,m)
                  local2(idim,7) = local2(idim,7) +
     1            ppc(lsum,-m)*scarray(iuse+2)*mpole(idim,l,m)
                  local2(idim,8) = local2(idim,8) +
     1            ppc(lsum,1-m)*scarray(iuse+3)*mpole(idim,l,m)
                  local2(idim,9) = local2(idim,9) +
     1            ppc(lsum,2-m)*scarray(iuse+4)*mpole(idim,l,m)
               enddo
               iuse = iuse+5
            enddo
         enddo
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
         do idim=1,nd
            pot(idim,itarg) = pot(idim,itarg)+local2(idim,1)/rscale
c
c     fld comes from l=1 modes
c
            rfac = 1.0d0/(rscale*rscale*sqrt(2.0d0))
            grad(idim,1,itarg)= grad(idim,1,itarg)+dreal(
     1              -rfac*(local2(idim,4)+local2(idim,2)))
            grad(idim,2,itarg)= grad(idim,2,itarg)+dreal(
     1              -rfac*ima*(local2(idim,4)-local2(idim,2)))
            grad(idim,3,itarg)=grad(idim,3,itarg)+
     1         dreal(local2(idim,3))/(rscale*rscale)
c
c     hess comes from l=2 modes
c
            rfac = sqrt(3.0d0)/(sqrt(2.0d0)*rscale*rscale*rscale)
            z0 = local2(idim,7)/(rscale*rscale*rscale)
            hess(idim,1,itarg)=hess(idim,1,itarg)+dreal(
     1                  rfac*(local2(idim,9)+local2(idim,5))-z0)
            hess(idim,2,itarg)=hess(idim,2,itarg)+dreal(
     1                  -rfac*(local2(idim,9)+local2(idim,5))-z0)
            hess(idim,3,itarg)=hess(idim,3,itarg)+dreal(2*z0)
            hess(idim,4,itarg)=hess(idim,4,itarg)+dreal(
     1                 rfac*ima*(local2(idim,9)-local2(idim,5)))
            hess(idim,5,itarg)=hess(idim,5,itarg)+dreal(
     1                 -rfac*(local2(idim,8)+local2(idim,6)))
            hess(idim,6,itarg)=hess(idim,6,itarg)+dreal(
     1                -rfac*ima*(local2(idim,8)-local2(idim,6)))
         enddo  
1000  continue
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dtaevalhessdini(nterms,scarray)
C***********************************************************************
C
c     Precomputes arrayused in 
c     local-local translation operator from an nterms expansion to an 
c     order 2 expansion (sufficient to compute pot/fld/hessian).     
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     scarray         : work array of size > 10*(nterms+2)**2
c-----------------------------------------------------------------------
      implicit none
      integer  nterms,l,j,k,m,ll,mm,iuse
      real *8 scarray(1)
      real *8 d
      real *8, allocatable :: cs(:,:)
      real *8, allocatable :: fact(:)
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
C***********************************************************************
      subroutine l3dtaevalh(nd,rscale,center,local,nterms,
     1            ztarg,ntarg,pot,grad,hess,scarray)
C***********************************************************************
cc
c     this subroutine evaluates the potentials, gradients and Hessians  
c     of a local expansion and increments inputs accordingly:
c
c     pot=pot + sum sum  mpole(n,m) r^{n}Y_nm(theta,phi) / sqrt(2n+1) 
c                   n   m
c
c     grad=grad + 
c             Gradient(sum sum mpole(n,m)r^{n}Y_nm(theta,phi)/sqrt(2n+1))
c                       n   m
c
c     hess=hess + 
c             Hessian(sum sum mpole(n,m)r^{n}Y_nm(theta,phi)/sqrt(2n+1))
c                      n   m
c
c    The method used direct translation (not rotation/zshift) of 
c    local expansion to one of order 2. It is reasonably optimized, 
c    precomputing the array of binomial/factorial terms that appear 
c     in the shift operator.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of local expansions
c     rscale :    scaling parameter (see formmp1l3d)
c     center :    expansion center
c     local  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of  target location
c    scarray :   precomputed array (MUST BE PRECEDED BY CALL TO
c                   L3DTAEVALHESSDINI(nterms,scarray))
c                   with dimension of scarray at least 10*(nterms+2)**2
c                   If nterms is changed, 
c                   l3dtaevalhessdini must be called again.
c   
c     OUTPUT:
c
c     pot    :   updated potential at ztarg
c     grad   :   updated gradient at ztarg 
c     hess   :   updated Hessian at ztarg
c                 ordered as dxx,dyy,dzz,dxy,dxz,dyz.
c--------------------------------------------------------------------
      implicit none
      integer nterms,ntarg,nd,itarg
      integer  l,m,lnew,mnew,ll,mm,iuse,j,k,ldiff,idim
      real *8 center(3),ztarg(3,ntarg)
      real *8 zdiff(3)
      real *8 scarray(*),rscale
      real *8 cphi,sphi,phi,theta,ctheta,d,dd,pi,rfac
      complex *16 local(nd,0:nterms,-nterms:nterms)
      complex *16 z0,ima,ephi1
      real *8 pot(nd,ntarg),grad(nd,3,ntarg)
      real *8 hess(nd,6,ntarg)
c
      real *8, allocatable :: pp(:,:)
      real *8, allocatable :: powers(:)
      complex *16, allocatable :: local2(:,:)
      complex *16, allocatable :: ppc(:,:)
      complex *16, allocatable :: ephi(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(pp(0:nterms,0:nterms))
      allocate(ppc(0:nterms,-nterms:nterms))
      allocate(powers(0:nterms+1))
      allocate(ephi(-nterms-1:nterms+1))
      allocate(local2(nd,9))

c     determine order of shifted expansion 
c
c
      do itarg = 1,ntarg
         do ll = 1,nd
         do l = 1,9
            local2(ll,l) = 0.0d0
         enddo
         enddo
c
         zdiff(1) = center(1) - ztarg(1,itarg)
         zdiff(2) = center(2) - ztarg(2,itarg)
         zdiff(3) = center(3) - ztarg(3,itarg)
         call cart2polar(zdiff,d,theta,phi)
c
         ctheta = dcos(theta)
         cphi = dcos(phi)
         sphi = dsin(phi)
         ephi1 = dcmplx(cphi,sphi)
C
C----- create array of powers of R and e^(i*m*phi).
c
         dd = d/rscale
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
         iuse = 1
         do l = 0,nterms
            do m = -l,l
               do idim=1,nd
                  local2(idim,1) = local2(idim,1) +
     1            ppc(l,m)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
         enddo
         do l = 1,nterms
            ldiff = l-1
            do m = -l,l-2
               mm = m+1
               do idim=1,nd
               local2(idim,2) = local2(idim,2) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               do idim=1,nd
               local2(idim,3) = local2(idim,3) +
     1         ppc(ldiff,m)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -l+2,l
               mm = m-1
               do idim=1,nd
               local2(idim,4) = local2(idim,4) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
         enddo
         do l = 2,nterms
            ldiff = l-2
            do m = -l,l-4
               mm = m+2
               do idim=1,nd
               local2(idim,5) = local2(idim,5) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -l+1,l-3
               mm = m+1
               do idim=1,nd
               local2(idim,6) = local2(idim,6) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -ldiff,ldiff
               mm = m
               do idim=1,nd
               local2(idim,7) = local2(idim,7) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -l+3,l-1
               mm = m-1
               do idim=1,nd
               local2(idim,8) = local2(idim,8) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
            do m = -l+4,l
               mm = m-2
               do idim=1,nd
               local2(idim,9) = local2(idim,9) +
     1         ppc(ldiff,mm)*scarray(iuse)*local(idim,l,m)
               enddo
               iuse = iuse+1
            enddo
         enddo
c
ccc      pi = 4.0d0*datan(1.0d0)
c
c     pot comes from 0,0 mode
c
         do idim=1,nd
            pot(idim,itarg) = pot(idim,itarg)+local2(idim,1)
c
c     fld comes from l=1 modes
c
            rfac = 1.0d0/(sqrt(2.0d0)*rscale)
            grad(idim,3,itarg)= grad(idim,3,itarg)+dreal(
     1                local2(idim,3)/rscale)
            grad(idim,1,itarg)= grad(idim,1,itarg)+dreal(
     1                -rfac*(local2(idim,4) + local2(idim,2)))
            grad(idim,2,itarg)= grad(idim,2,itarg)+dreal(
     1                -rfac*ima*(local2(idim,4) - local2(idim,2)))
c
c     hess comes from l=2 modes
c
ccc         rfac = rscale*rscale*sqrt(3.0d0)/sqrt(2.0d0)
            rfac = rfac*sqrt(3.0d0)/rscale
            z0 = local2(idim,7)/(rscale*rscale)
            hess(idim,1,itarg)=hess(idim,1,itarg)+dreal(
     1                  rfac*(local2(idim,9)+local2(idim,5))-z0)
            hess(idim,2,itarg)=hess(idim,2,itarg)+dreal(
     1                  -rfac*(local2(idim,9)+local2(idim,5))-z0)
            hess(idim,3,itarg)=hess(idim,3,itarg)+dreal(2*z0)
            hess(idim,4,itarg)=hess(idim,4,itarg)+dreal(
     1                 rfac*ima*(local2(idim,9)-local2(idim,5)))
            hess(idim,5,itarg)=hess(idim,5,itarg)+dreal(
     1                 -rfac*(local2(idim,8)+local2(idim,6)))
            hess(idim,6,itarg)=hess(idim,6,itarg)+dreal(
     1                -rfac*ima*(local2(idim,8)-local2(idim,6)))
         enddo
1000  continue
      enddo
      return
      end
c
c
c
c
