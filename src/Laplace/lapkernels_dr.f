      implicit real *8 (a-h,o-z)
      parameter (nd = 1)
      parameter (nt = 3)
      real *8 ztarg(3,nt),sources(3,10)
      real *8 charge(nd,10),pot(nd,nt),grad(nd,3,nt)
      real *8 hess(nd,6,nt)
      real *8 dipvec(nd,3,10)
c
c 
c
      sources(1,1) = 0.1d0     
      sources(2,1) = -0.2d0     
      sources(3,1) = 0.3d0     
      dipvec(1,1,1) = 1.1d0     
      dipvec(1,2,1) = -2.2d0     
      dipvec(1,3,1) = 3.3d0     
cc      dipvec(1,1,1) = 0
cc      dipvec(1,2,1) = 0     
cc      dipvec(1,3,1) = 0     
      charge(1,1) = 1.3d0
      sources(1,2) = -0.13d0     
      sources(2,2) = 0.2d0     
      sources(3,2) = -0.25d0     
      dipvec(1,1,2) = -2.13d0     
      dipvec(1,2,2) = 1.2d0     
      dipvec(1,3,2) = -3.25d0     
cc      dipvec(1,1,2) = 0     
cc      dipvec(1,2,2) = 0     
cc      dipvec(1,3,2) = 0     
      charge(1,2) = 2.3d0
      ns = 2
      ztarg(1,1) = 1.0d0
      ztarg(2,1) = 1.3d0
      ztarg(3,1) = 1.4d0
c
      hh = 1.0d-3
      ztarg(1,2) = 1.0d0 + hh 
      ztarg(2,2) = 1.3d0
      ztarg(3,2) = 1.4d0
      ztarg(1,3) = 1.0d0 - hh
      ztarg(2,3) = 1.3d0
      ztarg(3,3) = 1.4d0
c
c
      write(6,*) ' directc  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectcp(nd,sources,charge,ns,ztarg,nt,
     1            pot,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' FD z deriv ', (pot(1,2) - pot(1,3))/(2*hh)
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectcg(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      write(6,*) ' FD x deriv of p_x ', (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      write(6,*) ' FD x deriv of p_y ', (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      write(6,*) ' FD x deriv of p_z ', (v1-v2)/(2*hh)
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectch(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,hess,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      write(6,*) ' dir pxx ', hess(1,1,1)
      write(6,*) ' dir pxy ', hess(1,4,1)
      write(6,*) ' dir pxz ', hess(1,5,1)
c
c
      write(6,*) ' directd  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectdp(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' FD y deriv ', (pot(1,2) - pot(1,3))/(2*hh)
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectdg(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      write(6,*) ' FD x deriv of p_x ', (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      write(6,*) ' FD x deriv of p_y ', (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      write(6,*) ' FD x deriv of p_z ', (v1-v2)/(2*hh)
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectdh(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      write(6,*) ' dir pxx ', hess(1,1,1)
      write(6,*) ' dir pxy ', hess(1,4,1)
      write(6,*) ' dir pxz ', hess(1,5,1)
c
      write(6,*) ' directcd  tests -------------------------------'
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectcdp(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' FD x deriv ', (pot(1,2) - pot(1,3))/(2*hh)
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectcdg(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      v1 = grad(1,1,2)
      v2 = grad(1,1,3)
      write(6,*) ' FD x deriv of p_x ', (v1-v2)/(2*hh)
      v1 = grad(1,2,2)
      v2 = grad(1,2,3)
      write(6,*) ' FD x deriv of p_y ', (v1-v2)/(2*hh)
      v1 = grad(1,3,2)
      v2 = grad(1,3,3)
      write(6,*) ' FD x deriv of p_z ', (v1-v2)/(2*hh)
c
      call setzero(pot,grad,hess,nd,nt)
      call l3ddirectcdh(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
      write(6,*) ' dir cp ', pot
      write(6,*) ' dir grad ', grad(1,1,1)
      write(6,*) ' dir pxx ', hess(1,1,1)
      write(6,*) ' dir pxy ', hess(1,4,1)
      write(6,*) ' dir pxz ', hess(1,5,1)
      end

      subroutine setzero(pot,grad,hess,nd,nt)
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,nt),grad(nd,3,nt)
      real *8 hess(nd,6,nt)

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
