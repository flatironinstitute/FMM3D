c
c     Computation of spherical Bessel functions via recurrence
c
c**********************************************************************
      subroutine besseljs3d(nterms,z,scale,fjs,ifder,fjder)
c**********************************************************************
c
c PURPOSE:
c
c	This subroutine evaluates the first NTERMS spherical Bessel 
c	functions and if required, their derivatives.
c	It incorporates a scaling parameter SCALE so that
c       
c		fjs_n(z)=j_n(z)/SCALE^n
c		fjder_n(z)=\frac{\partial fjs_n(z)}{\partial z}
c
c INPUT:
c
c    nterms (integer): order of expansion of output array fjs 
c    z     (complex *16): argument of the spherical Bessel functions
c    scale    (real *8) : scaling factor (discussed above)
c    ifder  (integer): flag indicating whether to calculate "fjder"
c		          0	NO
c		          1	YES
c OUTPUT:
c    fjs   (complex *16): array of scaled Bessel functions.
c    fjder (complex *16): array of derivs of scaled Bessel functions.
c
c
      implicit none
      integer ier,nterms,ifder,lwfjs,ntop,i,ncntr
      real *8 scale,d0,d1,dc1,dc2,dcoef,dd,done,tiny,zero
      real *8 scalinv,sctot,upbound,upbound2,upbound2inv
      integer, allocatable :: iscale(:)
      complex *16, allocatable :: fjtmp(:)
      complex *16 wavek,fjs(0:nterms),fjder(0:*)
      complex *16 z,zinv,com,fjm1,fj0,fj1,zscale,ztmp
c
      data upbound/1.0d+32/, upbound2/1.0d+40/, upbound2inv/1.0d-40/
      data tiny/1.0d-200/,done/1.0d0/,zero/0.0d0/
c
c ... Initializing ...
c
      ier=0
c
c       set to asymptotic values if argument is sufficiently small
c
      if (abs(z).lt.tiny) then
        fjs(0) = done
        do i = 1, nterms
          fjs(i) = zero
	    enddo
c
	    if (ifder.eq.1) then
	      do i=0,nterms
	        fjder(i)=zero
	      enddo
          fjder(1)=done/(3*scale)
        endif
        return 
      endif
c
c ... Step 1: recursion up to find ntop, starting from nterms
c
      zinv=done/z

      fjm1 = zero
      fj0 = done

c
cc     note max point for upward recurrence is
c      hard coded to nterms + 1000,
c      this might cause loss of accuracy for some
c      arguments in the complex plane for large 
c      nterms. For example, it is a terrible idea
c      to use this code for z>>nterms^2

      lwfjs = nterms + 100000
      ntop=lwfjs
c
      do i=nterms,lwfjs
        dcoef=2*i+done
        fj1=dcoef*zinv*fj0-fjm1
        dd = abs(fj1)**2
        if(dd .gt. upbound2) then
          ntop=i+1
          goto 1300
        endif
        fjm1 = fj0
        fj0 = fj1
      enddo
 1300 continue

cc      call prinf('ntop=*',ntop,1)

      allocate(iscale(0:ntop),fjtmp(0:ntop))
c
c ... Step 2: Recursion back down to generate the unscaled jfuns:
c             if magnitude exceeds UPBOUND2, rescale and continue the 
c	      recursion (saving the order at which rescaling occurred 
c	      in array iscale.
c
      do i=0,ntop
         iscale(i)=0
      enddo
c
      fjtmp(ntop)=zero
      fjtmp(ntop-1)=done
      do i=ntop-1,1,-1
        dcoef=2*i+done
        fjtmp(i-1)=dcoef*zinv*fjtmp(i)-fjtmp(i+1)
        dd = abs(fjtmp(i-1))**2
        if (dd.gt.UPBOUND2) then
          fjtmp(i) = fjtmp(i)*UPBOUND2inv
          fjtmp(i-1) = fjtmp(i-1)*UPBOUND2inv
          iscale(i) = 1
        endif
      enddo
c
c ...  Step 3: go back up to the top and make sure that all
c              Bessel functions are scaled by the same factor
c              (i.e. the net total of times rescaling was invoked
c              on the way down in the previous loop).
c              At the same time, add scaling to fjs array.
c
      ncntr=0
      scalinv=done/scale
      sctot = 1.0d0
      do i=1,ntop
         sctot = sctot*scalinv
         if(iscale(i-1).eq.1) sctot=sctot*UPBOUND2inv
         fjtmp(i)=fjtmp(i)*sctot
      enddo
c
c ... Determine the normalization parameter:
c
      fj0=sin(z)*zinv
      fj1=fj0*zinv-cos(z)*zinv
c
      d0=abs(fj0)
      d1=abs(fj1)
      if (d1 .gt. d0) then
         zscale=fj1/(fjtmp(1)*scale)
      else
         zscale=fj0/fjtmp(0)
      endif
c
c ... Scale the jfuns by zscale:
c
      ztmp=zscale
      do i=0,nterms
         fjs(i)=fjtmp(i)*ztmp
      enddo
c
c ... Finally, calculate the derivatives if desired:
c
      if (ifder.eq.1) then
        fjder(0)=-fjs(1)*scale
        do i=1,nterms
          dc1=i/(2*i+done)
          dc2=done-dc1
          dc1=dc1*scalinv
          dc2=dc2*scale
          fjder(i)=(dc1*fjtmp(i-1)-dc2*fjtmp(i+1))*ztmp
        enddo
      endif
      return
      end
c
