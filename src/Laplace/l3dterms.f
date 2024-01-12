c
c      l3dterms - determine number of terms in mpole expansions 
c
c      l3dterms_list2 - build tables for the number of terms needed 
c           for all boxes in list 2
c
c      l3dterms_list2w - build table for the number of terms needed 
c           for all boxes in list 2, using worst case error analysis
c
c      l3dterms_list2e - build tables of the number of terms needed 
c           for all boxes in extended list 2
c
c      l3dterms_list2ew - build tables of the number of terms needed 
c           for all boxes in extended list 2, using worst case analysis
c
c-----------------------------------------------------------------------------
c
c
c
c***********************************************************************
      subroutine l3dterms(eps, nterms)
c***********************************************************************
c
c     Determine number of terms in mpole expansions.
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     nterms  required expansion order
c-----------------------------------------------------------------------
c
      implicit none
      integer *8 j,nterms
      real *8 xtemp1,eps
      real *8 z1,z2,hfun,jfun
c
      z1 = 1.5d0
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine l3dterms_far(eps, nterms)
c***********************************************************************
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk=0. 
c
c     The method is based on examining the decay of h_n * j_n.
c
c     This routine assumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0,0).
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     nterms  required expansion order
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 nterms,j,ier
      real *8 eps,jfun,hfun,z1,z2,xtemp1
c
      ier = 0
c
      z1 = 2.5d0
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = hfun*jfun
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine l3dterms_list2(eps, itable, ier)
c***********************************************************************
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     itable  table of needed  expansion orders in FMM.
c     ier     error flag (0 means correct execution).
c                        (1 means 1000 terms is insufficient).
c
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 nterms,ier,ii,jj,kk,i,j,k
      real *8 xtemp1,z1,z2,z3,jfun,hfun 
      real *8 eps,rr,dx,dy,dz
c
      integer *8 nterms_table(2:3,0:3,0:3)
      integer *8 itable(-3:3,-3:3,-3:3)
c
      ier = 0
c
      do ii=2,3
      do jj=0,3
      do kk=0,3
c
        dx=ii
        dy=jj
        dz=kk
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
        if( dz .gt. 0 ) dz=dz-.5
c
        rr=sqrt(dx*dx+dy*dy+dz*dz)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      z1 = rr
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          ier = 0
          nterms = j
          exit
        endif
        ier = 1
      enddo
c
        nterms_table(ii,jj,kk)=nterms
c
      enddo
      enddo
      enddo
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        do k=-3,3
        itable(i,j,k)=0
        enddo
        enddo
        enddo
c
        do k=-3,3
        do i=-3,3
        do j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j,k)=nterms_table(abs(i),abs(j),abs(k))
        else if( abs(j) .gt. 1) then
        itable(i,j,k)=nterms_table(abs(j),abs(i),abs(k))
        endif
c
        if( abs(i) .le. 1 .and. abs(j) .le. 1) then
        if( abs(k) .gt. 1 ) then
c
        if( abs(i) .ge. abs(j) ) then
        itable(i,j,k)=nterms_table(abs(k),abs(i),abs(j))
        else
        itable(i,j,k)=nterms_table(abs(k),abs(j),abs(i))
        endif
c
        endif
        endif
c
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
c***********************************************************************
      subroutine l3dterms_list2w(eps, itable, ier)
c***********************************************************************
c
c     Build nterms table for all boxes in list 2
c     Estimate worst case multipole to local translation operator errors.
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     itable  table of needed  expansion orders in FMM.
c     ier     error flag (0 means correct execution).
c                        (1 means 1000 terms is insufficient).
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 nterms,ier,ii,jj,kk,i,j,k
      real *8 z1,z2,z3,jfun,hfun,xtemp1,eps
      real *8 rr,dx,dy,dz
c
      integer *8 nterms_table(2:3,0:3,0:3)
      integer *8 itable(-3:3,-3:3,-3:3)
c
      ier = 0
c
      do ii=2,3
      do jj=0,3
      do kk=0,3
c
        dx=ii
        dy=jj
        dz=kk
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c        if( dz .gt. 0 ) dz=dz-.5
c
        rr=sqrt(dx*dx+dy*dy+dz*dz)
        rr=rr-sqrt(3.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      z1 = rr
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          ier = 0
          nterms = j
          exit
        endif
        ier = 1
      enddo
c
        nterms_table(ii,jj,kk)=nterms
c
      enddo
      enddo
      enddo
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        do k=-3,3
        itable(i,j,k)=0
        enddo
        enddo
        enddo
c
        do k=-3,3
        do i=-3,3
        do j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j,k)=nterms_table(abs(i),abs(j),abs(k))
        else if( abs(j) .gt. 1) then
        itable(i,j,k)=nterms_table(abs(j),abs(i),abs(k))
        endif
c
        if( abs(i) .le. 1 .and. abs(j) .le. 1) then
        if( abs(k) .gt. 1 ) then
c
        if( abs(i) .ge. abs(j) ) then
        itable(i,j,k)=nterms_table(abs(k),abs(i),abs(j))
        else
        itable(i,j,k)=nterms_table(abs(k),abs(j),abs(i))
        endif
c
        endif
        endif
c
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
c***********************************************************************
      subroutine l3dterms_list2e(eps, itable, ier)
c***********************************************************************
c
c
c     Build nterms table for all boxes in extended list 2
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     itable  table of needed  expansion orders in FMM.
c     ier     error flag (0 means correct execution).
c                        (1 means 1000 terms is insufficient).
c
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 nterms,ier,ii,jj,kk,i,j,k
      real *8 z1,z2,z3,jfun,hfun,xtemp1,eps
      real *8 rr,dx,dy,dz
c
      integer *8 nterms_table(2:7,0:7,0:7)
      integer *8 itable(-7:7,-7:7,-7:7)
c
      ier = 0
c
      do ii=2,7
      do jj=0,7
      do kk=0,7
c
        dx=ii
        dy=jj
        dz=kk
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
        if( dz .gt. 0 ) dz=dz-.5
c
        rr=sqrt(dx*dx+dy*dy+dz*dz)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      z1 = rr
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          ier = 0
          nterms = j
          exit
        endif
        ier = 1
      enddo
c
        nterms_table(ii,jj,kk)=nterms
c
      enddo
      enddo
      enddo
c
ccc        call prinf('nterms=*',nterms_table,2*6*6)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        do k=-7,7
        itable(i,j,k)=0
        enddo
        enddo
        enddo
c
        do k=-7,7
        do i=-7,7
        do j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j,k)=nterms_table(abs(i),abs(j),abs(k))
        else if( abs(j) .gt. 2) then
        itable(i,j,k)=nterms_table(abs(j),abs(i),abs(k))
        endif
c
        if( abs(i) .le. 2 .and. abs(j) .le. 2) then
        if( abs(k) .gt. 2 ) then
c
        if( abs(i) .ge. abs(j) ) then
        itable(i,j,k)=nterms_table(abs(k),abs(i),abs(j))
        else
        itable(i,j,k)=nterms_table(abs(k),abs(j),abs(i))
        endif
c
        endif
        endif
c
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
c***********************************************************************
      subroutine l3dterms_list2ew(eps, itable, ier)
c***********************************************************************
c
c
c     Build nterms table for all boxes in extended list 2
c     Estimate worst case multipole to local translation operator errors.
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     itable  table of needed  expansion orders in FMM.
c     ier     error flag (0 means correct execution).
c                        (1 means 1000 terms is insufficient).
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 nterms,ier,ii,jj,kk,i,j,k
      real *8 z1,z2,z3,jfun,hfun,xtemp1,eps
      real *8 rr,dx,dy,dz
c
      integer *8 nterms_table(2:7,0:7,0:7)
      integer *8 itable(-7:7,-7:7,-7:7)
c
      ier = 0
c
      do ii=2,7
      do jj=0,7
      do kk=0,7
c
        dx=ii
        dy=jj
        dz=kk
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c        if( dz .gt. 0 ) dz=dz-.5
c
        rr=sqrt(dx*dx+dy*dy+dz*dz)
        rr=rr-sqrt(3.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      z1 = rr
      z2 = dsqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          ier = 0
          nterms = j
          exit
        endif
        ier = 1
      enddo
c
        nterms_table(ii,jj,kk)=nterms
c
      enddo
      enddo
      enddo
c
ccc        call prinf('nterms=*',nterms_table,2*6*6)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        do k=-7,7
        itable(i,j,k)=0
        enddo
        enddo
        enddo
c
        do k=-7,7
        do i=-7,7
        do j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j,k)=nterms_table(abs(i),abs(j),abs(k))
        else if( abs(j) .gt. 2) then
        itable(i,j,k)=nterms_table(abs(j),abs(i),abs(k))
        endif
c
        if( abs(i) .le. 2 .and. abs(j) .le. 2) then
        if( abs(k) .gt. 2 ) then
c
        if( abs(i) .ge. abs(j) ) then
        itable(i,j,k)=nterms_table(abs(k),abs(i),abs(j))
        else
        itable(i,j,k)=nterms_table(abs(k),abs(j),abs(i))
        endif
c
        endif
        endif
c
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
c***********************************************************************
      subroutine l3dterms_eval(itype, eps, nterms, ier)
c***********************************************************************
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     INPUT:
c
c     itype   flag determining precise geometry (see code below).
c     eps     tolerance
c
c     OUTPUT:
c
c     nterms  required expansion order
c     ier     error flag (0 means correct execution).
c                        (1 means 1000 terms is insufficient).
c-----------------------------------------------------------------------------
c
      implicit none
      integer *8 itype,nterms,ier,ii,jj,kk,i,j,k
      real *8 z1,z2,z3,jfun,hfun,xtemp1,eps
      real *8 rr,dx,dy,dz
c
      ier = 0
c
      z1 = 1.5d0
      z2 = dsqrt(3d0)/2.d0
c
c     corners included
      if( itype .eq. 1 ) z2 = dsqrt(3d0)/2.d0
c     edges included, no corners
      if( itype .eq. 2 ) z2 = dsqrt(2d0)/2.d0
c     center only
      if( itype .eq. 3 ) z2 = 1.0d0/2.d0
c     center only, small interior sphere
      if( itype .eq. 4 ) z2 = 0.8d0/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
      enddo
      ier = 1
      return
      end

