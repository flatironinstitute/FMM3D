      implicit real *8 (a-h,o-z)
      parameter (nd = 1)
      parameter (nt = 3)
      real *8 ztarg(3,nt),sources(3,10)
      complex *16 charge(nd,10),pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      complex *16 dipvec(nd,3,10)
      complex *16 zk,v1,v2,zfd
      complex *16 zfd1,zfd2,zfd3
c
c 
      zk  = dcmplx(1.1d0,0.4d0)
c
      sources(1,1) = 0.1d0     
      sources(2,1) = -0.2d0     
      sources(3,1) = 0.3d0     
      dipvec(1,1,1) = dcmplx(1.1d0,0.3d0)     
      dipvec(1,2,1) = dcmplx(-2.2d0,0.2d0)     
      dipvec(1,3,1) = dcmplx(3.3d0,0.1d0)
cc      dipvec(1,1,1) = 0
cc      dipvec(1,2,1) = 0     
cc      dipvec(1,3,1) = 0     
      charge(1,1) = dcmplx(1.3d0,0.11d0)
      sources(1,2) = -0.13d0     
      sources(2,2) = 0.2d0     
      sources(3,2) = -0.25d0     
      dipvec(1,1,2) = dcmplx(-2.13d0,0.1d0)
      dipvec(1,2,2) = dcmplx(1.2d0,0.3d0)
      dipvec(1,3,2) = dcmplx(-3.25d0,0.2d0)
cc      dipvec(1,1,2) = 0     
cc      dipvec(1,2,2) = 0     
cc      dipvec(1,3,2) = 0     
      charge(1,2) = dcmplx(2.3d0,0.21d0)
      ns = 2
      ztarg(1,1) = 1.0d0
      ztarg(2,1) = 1.3d0
      ztarg(3,1) = 1.4d0
c
c     pick component in which to do finite difference
c
c
      icomp = 3

      hh = 1.0d-5
      ztarg(1,2) = 1.0d0 
      ztarg(2,2) = 1.3d0 
      ztarg(3,2) = 1.4d0 
      ztarg(1,3) = 1.0d0 
      ztarg(2,3) = 1.3d0 
      ztarg(3,3) = 1.4d0 
      ztarg(icomp,2) = ztarg(icomp,2) + hh
      ztarg(icomp,3) = ztarg(icomp,3) - hh
c
c
      write(6,*) ' directc  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectcp(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      zfd =  (pot(1,2) - pot(1,3))/(2*hh)
      write(6,*) ' FD deriv ', zfd
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectcg(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      zfd1 = (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      zfd2 = (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      zfd3 = (v1-v2)/(2*hh)
      write(6,*) ' FD deriv of p_x ', zfd1 
      write(6,*) ' FD deriv of p_y ', zfd2
      write(6,*) ' FD deriv of p_z ', zfd3
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectch(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,hess,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      if (icomp.eq.1) then
         write(6,*) ' dir pxx ', hess(1,1,1)
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pxz ', hess(1,5,1)
         exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
         exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
         exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
         write(6,*) ' exx  ', exx
         write(6,*) ' exy  ', exy
         write(6,*) ' exz  ', exz
      else if (icomp.eq.2) then
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pyy ', hess(1,2,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
         eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
         eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
         write(6,*) ' exy  ', exy
         write(6,*) ' eyy  ', eyy
         write(6,*) ' eyz  ', eyz
      else if (icomp.eq.3) then
         write(6,*) ' dir pxz ', hess(1,4,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         write(6,*) ' dir pzz ', hess(1,3,1)
         exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
         eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
         ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
         write(6,*) ' exz  ', exz
         write(6,*) ' eyz  ', eyz
         write(6,*) ' ezz  ', ezz
      endif
c
c
      write(6,*) ' directd  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectdp(nd,zk,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      zfd =  (pot(1,2) - pot(1,3))/(2*hh)
      write(6,*) ' FD deriv ', zfd
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectdg(nd,zk,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      zfd1 = (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      zfd2 = (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      zfd3 = (v1-v2)/(2*hh)
      write(6,*) ' FD x deriv of p_x ', zfd1
      write(6,*) ' FD x deriv of p_y ', zfd2
      write(6,*) ' FD x deriv of p_z ', zfd3
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectdh(nd,zk,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      if (icomp.eq.1) then
         write(6,*) ' dir pxx ', hess(1,1,1)
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pxz ', hess(1,5,1)
         exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
         exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
         exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
         write(6,*) ' exx  ', exx
         write(6,*) ' exy  ', exy
         write(6,*) ' exz  ', exz
      else if (icomp.eq.2) then
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pyy ', hess(1,2,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
         eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
         eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
         write(6,*) ' exy  ', exy
         write(6,*) ' eyy  ', eyy
         write(6,*) ' eyz  ', eyz
      else if (icomp.eq.3) then
         write(6,*) ' dir pxz ', hess(1,4,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         write(6,*) ' dir pzz ', hess(1,3,1)
         exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
         eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
         ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
         write(6,*) ' exz  ', exz
         write(6,*) ' eyz  ', eyz
         write(6,*) ' ezz  ', ezz
      endif

c
      write(6,*) ' directcd  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectcdp(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      zfd =  (pot(1,2) - pot(1,3))/(2*hh)
      write(6,*) ' FD deriv ', zfd
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectcdg(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      zfd1 = (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      zfd2 = (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      zfd3 = (v1-v2)/(2*hh)
      write(6,*) ' FD deriv of p_x ', zfd1
      write(6,*) ' FD deriv of p_y ', zfd2
      write(6,*) ' FD deriv of p_z ', zfd3

c
      call setzero(pot,grad,hess,nd,nt)
      call h3ddirectcdh(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
      write(6,*) ' dir cp '
      do i = 1,nt
         write(6,*)pot(1,i)
      enddo
      write(6,*) ' analytic grad ', grad(1,icomp,1)
      egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
      write(6,*) ' egrad  ', egrad
      if (icomp.eq.1) then
         write(6,*) ' dir pxx ', hess(1,1,1)
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pxz ', hess(1,5,1)
         exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
         exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
         exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
         write(6,*) ' exx  ', exx
         write(6,*) ' exy  ', exy
         write(6,*) ' exz  ', exz
      else if (icomp.eq.2) then
         write(6,*) ' dir pxy ', hess(1,4,1)
         write(6,*) ' dir pyy ', hess(1,2,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
         eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
         eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
         write(6,*) ' exy  ', exy
         write(6,*) ' eyy  ', eyy
         write(6,*) ' eyz  ', eyz
      else if (icomp.eq.3) then
         write(6,*) ' dir pxz ', hess(1,4,1)
         write(6,*) ' dir pyz ', hess(1,6,1)
         write(6,*) ' dir pzz ', hess(1,3,1)
         exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
         eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
         ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
         write(6,*) ' exz  ', exz
         write(6,*) ' eyz  ', eyz
         write(6,*) ' ezz  ', ezz
      endif

      end

      subroutine setzero(pot,grad,hess,nd,nt)
      implicit real *8 (a-h,o-z)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)

      do ii=1,nd
      do jj=1,nt
         pot(ii,jj) = 0.0d0
         grad(ii,1,jj) = 0.0d0
         grad(ii,2,jj) = 0.0d0
         grad(ii,3,jj) = 0.0d0
         hess(ii,1,jj) = 0.0d0
         hess(ii,2,jj) = 0.0d0
         hess(ii,3,jj) = 0.0d0
         hess(ii,4,jj) = 0.0d0
         hess(ii,5,jj) = 0.0d0
         hess(ii,6,jj) = 0.0d0
      enddo
      enddo
      return
      end
