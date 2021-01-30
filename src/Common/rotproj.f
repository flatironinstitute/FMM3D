c
c     ROTATION VIA PROJECTION  
c
c     Multipole expansions for complex-valued functions.
c
c     Requires FFT and Associated Legendre Function Libraries.
c
c     User-callable routine is rotviaproj.
c     The other routines are used internally.
c
c***********************************************************************
      subroutine rotviaproj(nd,beta,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the y-axis determined by angle beta.
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       nd: number of multipole expansions
c       beta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
c
c
c
      implicit none
      integer(8) idim,nd
      integer(8) nquad,ier,m1,m2,nterms,lmp,lmpn, next235_cproj_vec
      double precision beta
      double complex mpole(nd,0:lmp,-lmp:lmp)
      double complex marray2(nd,0:lmpn,-lmpn:lmpn)

      double precision, allocatable :: cthetas(:),cphis(:)
      double precision, allocatable :: sthetas(:),sphis(:)
      double precision, allocatable :: ynm(:,:),ynmd(:,:)
      double precision, allocatable :: rat1(:,:),rat2(:,:)
      double precision, allocatable :: wsave(:)
      
      double complex, allocatable :: avec(:),bvec(:)
      double complex, allocatable :: uder(:,:,:),uval(:,:,:)
      double complex, allocatable :: ephis(:)



      nquad = next235_cproj_vec((2*nterms+2)*1.0d0)

c
cc      allocate all temporary arrays
c

      allocate(cthetas(nquad),cphis(nquad))
      allocate(sthetas(nquad),sphis(nquad))
      allocate(ynm(0:nterms,0:nterms),ynmd(0:nterms,0:nterms))
      allocate(rat1(0:nterms,0:nterms),rat2(0:nterms,0:nterms))
      allocate(wsave(4*nquad+20))

      allocate(avec(nquad),bvec(nquad),uder(nd,nquad,0:nterms))
      allocate(uval(nd,nquad,0:nterms),ephis(-nterms:nterms))
c
c     Algorithm:
c     1) get locations of quadrature nodes
c     2) evaluate u and du/dtheta
c     3) project onto spherical harmonics.
c
      call getmeridian(beta,nquad,cthetas,sthetas,cphis,sphis)    
      call evalall(nd,beta,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
      call projectonynm(nd,nquad,uval,uder,ynm,ynmd,marray2,lmpn,
     1           nterms,m2,wsave,avec,bvec,rat1,rat2)

      return
      end
c
c
      function next235_cproj_vec(base)
      implicit none
      integer(8) next235_cproj_vec, numdiv
      real*8 base
c ----------------------------------------------------------------------
c     integer(8) function next235_cproj returns a multiple of 2, 3, and 5
c
c     next235_cproj = 2^p 3^q 5^r >= base  where p>=1, q>=0, r>=0
************************************************************************
      next235_cproj_vec = 2 * int(base/2d0+.9999d0)
      if (next235_cproj_vec.le.0) next235_cproj_vec = 2

100   numdiv = next235_cproj_vec
      do while (numdiv/2*2 .eq. numdiv)
         numdiv = numdiv /2
      enddo
      do while (numdiv/3*3 .eq. numdiv)
         numdiv = numdiv /3
      enddo
      do while (numdiv/5*5 .eq. numdiv)
         numdiv = numdiv /5
      enddo
      if (numdiv .eq. 1) return
      next235_cproj_vec = next235_cproj_vec + 2
      goto 100
      end
c
c
c***********************************************************************
      subroutine getmeridian(beta,nquad,cthetas,sthetas,cphis,
     1             sphis)
C***********************************************************************
C     Purpose:
C
C           For a rotation of angle BETA about the y-axis, this
C           subroutine returns the NQUAD equispaced nodes on the 
C           rotated equator in the original coordinate system.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta  = angle of rotation
C     nquad = number of quadrature nodes in equator.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     cthetas = cos(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     sthetas = sin(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     cphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C     sphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C
C***********************************************************************
      implicit double precision (a-h,o-z)
      integer(8) nquad
      double precision cthetas(nquad)
      double precision sthetas(nquad)
      double precision cphis(nquad)
      double precision sphis(nquad)
C
      pi = 4.0d0*datan(1.0d0)
C
      ca = cos(beta)
      sa = sin(beta)
      do i = 1,nquad
	 im1 = i-1
         phi = 2*pi*im1/nquad
	 theta = pi/2.0d0
         xp = cos(phi)*sin(theta)
         yp = sin(phi)*sin(theta)
         zp = cos(theta)
         x = ca*xp + sa*zp
         y = yp
         z = -sa*xp + ca*zp
         proj = sqrt(x**2+y**2)
	 if (proj.le.1.0d-16) then
	    cphis(i) = 1.0d0
	    sphis(i) = 0.0d0
	 else
	    cphis(i) = x/proj
	    sphis(i) = y/proj
	 endif
	 cthetas(i) = z
	 sthetas(i) = proj
      enddo
      return
      end
C
C***********************************************************************
      subroutine evalall(nd,beta,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd      : numb3re of multipole expansions
C     beta    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit double precision (a-h,o-z)
      integer(8) nd,idim,lmp,nterms
      integer(8) ndeg,morder, nquad
      double precision cthetas(nquad),cphis(nquad)
      double precision sthetas(nquad),sphis(nquad)
      double precision ynm(0:nterms,0:nterms)
      double precision ynmd(0:nterms,0:nterms)
      double complex mpole(nd,0:lmp,-lmp:lmp)
      double complex ephi1,ephis(-nterms:nterms)
      double complex uder(nd,nquad,0:nterms),uval(nd,nquad,0:nterms)
      double complex uv,utheta,uphi,ztmp1,ztmp2,ztsum,ztdif
      double complex ux,uy,uz,ima
      double precision rat1(0:nterms,0:nterms)
      double precision rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      cbeta = cos(beta)
      sbeta = -sin(beta)
      call ylgndrini(nterms,rat1,rat2)
      nquad2=nquad/2
      if( mod(nquad2,2) .eq. 0 ) nquad4=nquad2/2+1
      if( mod(nquad2,2) .eq. 1 ) nquad4=nquad2/2+1
c
      do jj=1,nquad4
         ctheta = cthetas(jj)
         stheta = sthetas(jj)
         cphi = cphis(jj)
         sphi = sphis(jj)
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
         call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj))
         ephis(1) = ephi1
         ephis(-1) = dconjg(ephi1)
	     do i = 2,nterms
	        ephis(i) = ephis(i-1)*ephi1
	        ephis(-i) = dconjg(ephis(i))
	     enddo
	     do ndeg = 0,nterms
            do idim=1,nd
	           uv=0
	           utheta=0
	           uphi=0
	           do morder = 1,ndeg
                  ztmp1=ephis(morder)*mpole(idim,ndeg,morder)
                  ztmp2=ephis(-morder)*mpole(idim,ndeg,-morder)
	              ztsum=ztmp1+ztmp2
	              ztdif=ztmp1-ztmp2
	              uv=uv+ynm(ndeg,morder)*ztsum
	              utheta=utheta+ynmd(ndeg,morder)*ztsum
	              uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	           enddo
	           uv=stheta*uv+ynm(ndeg,0)*mpole(idim,ndeg,0)
	           utheta=utheta+ynmd(ndeg,0)*stheta*mpole(idim,ndeg,0)
c
c       ... apply the periodizing operator
c
                uval(idim,jj,ndeg) = uv
                uder(idim,jj,ndeg) = (utheta*proj2+uphi*ima*proj1)
                if( mod(ndeg,2) .eq. 0 ) then
                   uval(idim,jj+nquad/2,ndeg) = +uval(idim,jj,ndeg)
                   uder(idim,jj+nquad/2,ndeg) = -uder(idim,jj,ndeg)
                endif
                if( mod(ndeg,2) .eq. 1 ) then
                   uval(idim,jj+nquad/2,ndeg) = -uval(idim,jj,ndeg)
                   uder(idim,jj+nquad/2,ndeg) = +uder(idim,jj,ndeg)
                endif
	         enddo
          enddo

          if_reflect=1
          if( jj .eq. 1 ) if_reflect=0
          if( jj .eq. nquad4 .and. mod(nquad2,2) .eq. 0 ) if_reflect=0
c
          if( if_reflect .eq. 1 ) then
             jjr=nquad2 - jj + 2
             call ylgndr2pm_opt(nterms,ynm,ynmd)
             ctheta = -cthetas(jj)
             stheta = sthetas(jj)
             cphi = -cphis(jj)
             sphi = sphis(jj)
             dir1 = -sbeta
             dir2 = 0
             dir3 = cbeta
             tang1 = cphi*ctheta
             tang2 = sphi*ctheta
             tang3 = -stheta
             proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
             tang1 = -sphi
             tang2 = cphi
             tang3 = 0
             proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
             ephi1 = dcmplx(cphi,sphi)
             ephis(1) = ephi1
	         ephis(-1) = dconjg(ephi1)
	         do i = 2,nterms
	            ephis(i) = ephis(i-1)*ephi1
	            ephis(-i) = dconjg(ephis(i))
	         enddo
	         do ndeg = 0,nterms
                do idim=1,nd
	               uv=0
	               utheta=0
	               uphi=0
	               do morder = 1,ndeg
                      ztmp1=ephis(morder)*mpole(idim,ndeg,morder)
                      ztmp2=ephis(-morder)*mpole(idim,ndeg,-morder)
	                  ztsum=ztmp1+ztmp2
	                  ztdif=ztmp1-ztmp2
	                  uv=uv+ynm(ndeg,morder)*ztsum
	                  utheta=utheta+ynmd(ndeg,morder)*ztsum
	                  uphi=uphi-ynm(ndeg,morder)*morder*ztdif
	               enddo
	               uv=stheta*uv+ynm(ndeg,0)*mpole(idim,ndeg,0)
                   utheta=utheta+ynmd(ndeg,0)*stheta*mpole(idim,ndeg,0)
c
c       ... apply the periodizing operator
c
                   uval(idim,jjr,ndeg) = uv
                   uder(idim,jjr,ndeg) = (utheta*proj2+uphi*ima*proj1)
                   if( mod(ndeg,2) .eq. 0 ) then
                      uval(idim,jjr+nquad/2,ndeg) = +uval(idim,jjr,ndeg)
                      uder(idim,jjr+nquad/2,ndeg) = -uder(idim,jjr,ndeg)
                   endif
                   if( mod(ndeg,2) .eq. 1 ) then
                      uval(idim,jjr+nquad/2,ndeg)=-uval(idim,jjr,ndeg)
                      uder(idim,jjr+nquad/2,ndeg)=+uder(idim,jjr,ndeg)
                   endif
	            enddo
             enddo
         endif
      enddo

      return
      end
C
C
C
C
C
C***********************************************************************
      subroutine projectonynm(nd,nquad,uval,uder,
     1           ynm,ynmd,marray,lmpn,nterms,m2,wsave,avec,bvec,
     $           rat1,rat2)
C***********************************************************************
C
C     This subroutine projects from values on equator for each multipole
C     order (uval, uder = dudthteta) 
C     onto spherical harmonics
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nd       : number of multipole expansions
C     nquad    : number of points on equator
C     uval     : F values on equator
C     uder     : dFdtheta values on equator
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     lmpn     : leading dim of marray (must exceed nterms)
C     nterms   : order of expansion
C     m2       : NOT IMPLEMENTED (for reduced number of degrees in 
C                expansion (second index)
C     wsave    : work array for FFT (dimension at least 4*nquad+20)
C     avec     : work array of length nquad for FFT (complex)
C     bvec     : work array of length nquad for FFT (complex)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     marray   : rotated expansion 
C
C
C***********************************************************************
      implicit double precision (a-h,o-z)
      integer(8) nquad, norder,nterms,lmpn,m2
      integer(8) nd,idim
      integer nquad_4
      double complex ephi,ephi1,uval(nd,nquad,0:1)
      double complex uder(nd,nquad,0:1)
      double complex utheta,uphi,ztmp1,ztmp2
      double complex alpha,beta,ima
      double complex marray(nd,0:lmpn,-lmpn:lmpn)
      double precision ynm(0:nterms,0:nterms)
      double precision ynmd(0:nterms,0:nterms)
      double precision wsave(4*nquad+20)
      double complex avec(nquad)
      double complex bvec(nquad)
      double precision rat1(0:nterms,0:nterms)
      double precision rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
c
      ctheta = 0.0d0
      stheta = 1.0d0
      h = 1.0d0/nquad
      nquad_4 = nquad
      call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
      call zffti(nquad_4,wsave)
      do norder=0,nterms
         d=sqrt(2*norder+1.0d0)
         do idim=1,nd
   	        do ii = 1,nquad
               avec(ii) = uval(idim,ii,norder)*d + uder(idim,ii,norder)
            enddo
	        call zfftf(nquad_4,avec,wsave)
            do m = -norder,norder
	           if (m.ge.0)  alpha = avec(m+1)*h
	           if (m.lt.0)  alpha = avec(nquad+m+1)*h
               marray(idim,norder,m) = alpha/
     1            (ynm(norder,abs(m))*d - (ynmd(norder,abs(m))))
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
c***********************************************************************
