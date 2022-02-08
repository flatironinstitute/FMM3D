c
c     Computation of spherical Bessel functions via recurrence
c
c**********************************************************************
      subroutine jfuns3d(ier,nterms,z,scale,fjs,ifder,fjder,
     1	      lwfjs,iscale,ntop)
      implicit none
      integer(8) ier,nterms,ifder,lwfjs,ntop,i,ncntr
      real *8 scale,d0,d1,dc1,dc2,dcoef,dd,done,tiny,zero
      real *8 scalinv,sctot,upbound,upbound2,upbound2inv
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
c	NOTE: The scaling parameter SCALE is meant to be used when
c             abs(z) < 1, in which case we recommend setting
c	      SCALE = abs(z). This prevents the fjs_n from 
c             underflowing too rapidly.
c	      Otherwise, set SCALE=1.
c	      Do not set SCALE = abs(z) if z could take on the 
c             value zero. 
c             In an FMM, when forming an expansion from a collection of
c             sources, set SCALE = min( abs(k*r), 1)
c             where k is the Helmholtz parameter and r is the box dimension
c             at the relevant level.
c
c INPUT:
c
c    nterms (integer(8)): order of expansion of output array fjs 
c    z     (complex *16): argument of the spherical Bessel functions
c    scale    (real *8) : scaling factor (discussed above)
c    ifder  (integer(8)): flag indicating whether to calculate "fjder"
c		          0	NO
c		          1	YES
c    lwfjs  (integer(8)): upper limit of input arrays 
c                         fjs(0:lwfjs) and iscale(0:lwfjs)
c    iscale (integer(8)): integer(8) workspace used to keep track of 
c                         internal scaling
c
c OUTPUT:
c
c    ier    (integer(8)): error return code 
c                         ier=0 normal return;
c                         ier=8 insufficient array dimension lwfjs
c    fjs   (complex *16): array of scaled Bessel functions.
c    fjder (complex *16): array of derivs of scaled Bessel functions.
c    ntop  (integer(8)) : highest index in arrays fjs that is nonzero
c
c       NOTE, that fjs and fjder arrays must be at least (nterms+2)
c       complex *16 elements long.
c
c
      integer(8) iscale(0:lwfjs)
      complex *16 wavek,fjs(0:lwfjs),fjder(0:*)
      complex *16 z,zinv,com,fj0,fj1,zscale,ztmp
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
c
         RETURN
      endif
c
c ... Step 1: recursion up to find ntop, starting from nterms
c
      ntop=0
      zinv=done/z
      fjs(nterms)=done
      fjs(nterms-1)=zero
c
      do 1200 i=nterms,lwfjs
         dcoef=2*i+done
         ztmp=dcoef*zinv*fjs(i)-fjs(i-1)
         fjs(i+1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd .gt. upbound2) then
            ntop=i+1
            goto 1300
         endif
 1200 continue
 1300 continue
      if (ntop.eq.0) then
         ier=8
         return
      endif
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
      fjs(ntop)=zero
      fjs(ntop-1)=done
      do 2200 i=ntop-1,1,-1
	 dcoef=2*i+done
         ztmp=dcoef*zinv*fjs(i)-fjs(i+1)
         fjs(i-1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd.gt.UPBOUND2) then
            fjs(i) = fjs(i)*UPBOUND2inv
            fjs(i-1) = fjs(i-1)*UPBOUND2inv
            iscale(i) = 1
         endif
 2200 continue
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
         fjs(i)=fjs(i)*sctot
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
         zscale=fj1/(fjs(1)*scale)
      else
         zscale=fj0/fjs(0)
      endif
c
c ... Scale the jfuns by zscale:
c
      ztmp=zscale
      do i=0,nterms
         fjs(i)=fjs(i)*ztmp
      enddo
c
c ... Finally, calculate the derivatives if desired:
c
      if (ifder.eq.1) then
         fjs(nterms+1)=fjs(nterms+1)*ztmp
c
         fjder(0)=-fjs(1)*scale
         do i=1,nterms
            dc1=i/(2*i+done)
            dc2=done-dc1
            dc1=dc1*scalinv
            dc2=dc2*scale
            fjder(i)=dc1*fjs(i-1)-dc2*fjs(i+1)
         enddo
      endif
      return
      end
c
