c
c
c      this file contains some common routines for helmrouts3d
c      and helmrouts3d_vec
c
c      cart2polar: utility function.
c                  converts Cartesian coordinates into polar
c                  representation needed by other routines.
c
c      h3d01: computes h0, h1 (first two spherical Hankel fns.)
c      h3dall: computes Hankel functions of all orders and scales them
c
c      h3dadd: adds one expansion to another
c
c      h3dadd_trunc: adds one expansion to another
c                 allowing for different spherical harmonic array dimensions  
c
c
c**********************************************************************
      subroutine cart2polar(zat,r,theta,phi)
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
      implicit real *8 (a-h,o-z)
      real *8 zat(3)
      complex *16 ephi,eye
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
c**********************************************************************
      subroutine h3d01(z,h0,h1)
c**********************************************************************
c
c     Compute spherical Hankel functions of order 0 and 1 
c
c     h0(z)  =   exp(i*z)/(i*z),
c     h1(z)  =   - h0' = -h0*(i-1/z) = h0*(1/z-i)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	z   :  argument of Hankel functions
c              if abs(z)<1.0d-15, returns zero.
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	h0  :  h0(z)    (spherical Hankel function of order 0).
c	h1  :  -h0'(z)  (spherical Hankel function of order 1).
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 z,zinv,eye,cd,h0,h1
      data eye/(0.0d0,1.0d0)/, thresh/1.0d-15/, done/1.0d0/
c
      if (abs(z).lt.thresh) then
         h0=0.0d0
         h1=0.0d0
         return
      endif
c
c     Otherwise, use formula
c
      cd = eye*z
      h0=exp(cd)/cd
      h1=h0*(done/z - eye)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dall(nterms,z,scale,hvec,ifder,hder)
c**********************************************************************
c
c     This subroutine computes scaled versions of the spherical Hankel 
c     functions h_n of orders 0 to nterms.
c
c       	hvec(n)= h_n(z)*scale^(n)
c
c     The parameter SCALE is useful when |z| < 1, in which case
c     it damps out the rapid growth of h_n as n increases. In such 
c     cases, we recommend setting 
c                                 
c               scale = |z|
c
c     or something close. If |z| > 1, set scale = 1.
c
c     If the flag IFDER is set to one, it also computes the 
c     derivatives of h_n.
c
c		hder(n)= h_n'(z)*scale^(n)
c
c     NOTE: If |z| < 1.0d-15, the subroutine returns zero.
c     
c-----------------------------------------------------------------------
c     INPUT:
c
c     nterms  : highest order of the Hankel functions to be computed.
c     z       : argument of the Hankel functions.
c     scale   : scaling parameter discussed above
c     ifder   : flag indcating whether derivatives should be computed.
c		ifder = 1   ==> compute 
c		ifder = 0   ==> do not compute
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     hvec    : the vector of spherical Hankel functions 
c     hder    : the derivatives of the spherical Hankel functions 
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 hvec(0:1),hder(0:1)
      complex *16 zk2,z,zinv,ztmp,fhextra
c
      data thresh/1.0d-15/,done/1.0d0/
c
c     If |z| < thresh, return zeros.
c
      if (abs(z).lt.thresh) then
         do i=0,nterms
            hvec(i)=0
            hder(i)=0
         enddo
         return
      endif
c
c     Otherwise, get h_0 and h_1 analytically and the rest via
c     recursion.
c
      call h3d01(z,hvec(0),hvec(1))
      hvec(0)=hvec(0)
      hvec(1)=hvec(1)*scale
c
c     From Abramowitz and Stegun (10.1.19)
c
c     h_{n+1}(z)=(2n+1)/z * h_n(z) - h_{n-1}(z)
c
c     With scaling:
c
c     hvec(n+1)=scale*(2n+1)/z * hvec(n) -(scale**2) hvec(n-1)
c
      scal2=scale*scale
      zinv=scale/z
      do i=1,nterms-1
	 dtmp=(2*i+done)
	 ztmp=zinv*dtmp
	 hvec(i+1)=ztmp*hvec(i)-scal2*hvec(i-1)
      enddo
c
c     From Abramowitz and Stegun (10.1.21)
c
c	h_{n}'(z)= h_{n-1}(z) - (n+1)/z * h_n(z)
c
c     With scaling:
c
c     hder(n)=scale* hvec(n-1) - (n+1)/z * hvec(n)
c
c
      if (ifder.eq.1) then
c
	 hder(0)=-hvec(1)/scale
         zinv=1.0d0/z
         do i=1,nterms
	    dtmp=(i+done)
	    ztmp=zinv*dtmp
	    hder(i)=scale*hvec(i-1)-ztmp*hvec(i)
	 enddo
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dadd(mpole,mpole2,nterms)
c**********************************************************************
c
c     add mpole to mpole2
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:nterms,-nterms:nterms)
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
      subroutine h3dadd_trunc(mpole,mpole2,nterms,ldc)
c**********************************************************************
c
c     add mpole to mpole2
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpole2(0:ldc,-ldc:ldc)
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
c
c 
      subroutine h3dzero(mpole,nterms)
      implicit real *8 (a-h,o-z)
      complex *16 mpole(0:nterms,-nterms:nterms)
      
      do i=0,nterms
        do j=-i,i
          mpole(i,j) = 0
        enddo
      enddo



      return
      end

      
