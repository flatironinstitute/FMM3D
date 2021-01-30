c
c     This file contains the Stokes direct kernel evaluators
c
c**********************************************************************
c      
c     We take the following conventions for the Stokes kernels
c
c     For a source y and target x, let r_i = x_i-y_i
c     and let r = sqrt(r_1^2 + r_2^2 + r_3^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     (without the 1/4pi scaling) are
c
c     G_{ij}(x,y) = (r_i r_j)/(2r^3) + delta_{ij}/(2r)
c     P_j(x,y) = r_j/r^3
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, (without the 1/4pi scaling) are
c     
c     T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5
c     PI_{jk} = -2 delta_{jk} + 6 r_j r_k/r^5      
c

      subroutine st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1     pot,pre,grad,thresh)
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer(8)
c        number of densities
c     
c     nsource in: integer(8)
c        number of sources
c
c     source  in: double precision (3,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,3,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     ntarg   in: integer(8)
c        number of targs 
c
c     targ    in: double precision (3,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,3,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,3,3,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,stoklet,ns
cf2py intent(in) targ,nt,thresh
cf2py intent(out) pot,pre,grad

      integer(8) nd, ns, nt, istress
      real *8 sources(3,ns),targ(3,nt),strslet(nd,3,ns)
      real *8 strsvec(nd,3,ns),stoklet(nd,3,ns)
      real *8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(3), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pv, dmu(3), dnu(3), temp, r, r2, r3, r5
      real *8 dmunu      
      real *8 threshsq

      integer(8) i, j, idim, l

      threshsq = thresh**2


c     stokeslet contribution
      
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)
            zdiff(3) = targ(3,i)-sources(3,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
            if (r2 .lt. threshsq) goto 10

            r = sqrt(r2)
            r3 = r*r2
            
            do idim = 1,nd

               pot(idim,1,i) = pot(idim,1,i) + stoklet(idim,1,j)/(2*r)
               pot(idim,2,i) = pot(idim,2,i) + stoklet(idim,2,j)/(2*r)
               pot(idim,3,i) = pot(idim,3,i) + stoklet(idim,3,j)/(2*r)
               
               pl = (zdiff(1)*stoklet(idim,1,j) +
     1              zdiff(2)*stoklet(idim,2,j) +
     2              zdiff(3)*stoklet(idim,3,j))/(r3*2)
               
               pot(idim,1,i) = pot(idim,1,i) + zdiff(1)*pl
               pot(idim,2,i) = pot(idim,2,i) + zdiff(2)*pl
               pot(idim,3,i) = pot(idim,3,i) + zdiff(3)*pl

               grad(idim,1,1,i) = grad(idim,1,1,i) + pl
               grad(idim,2,2,i) = grad(idim,2,2,i) + pl
               grad(idim,3,3,i) = grad(idim,3,3,i) + pl               

               d1 = stoklet(idim,1,j)/(r3*2) - zdiff(1)*pl*3.0d0/r2
               d2 = stoklet(idim,2,j)/(r3*2) - zdiff(2)*pl*3.0d0/r2
               d3 = stoklet(idim,3,j)/(r3*2) - zdiff(3)*pl*3.0d0/r2
               
               grad(idim,1,1,i) = grad(idim,1,1,i) + d1*zdiff(1)
               grad(idim,2,1,i) = grad(idim,2,1,i) + d2*zdiff(1)
               grad(idim,3,1,i) = grad(idim,3,1,i) + d3*zdiff(1)
               grad(idim,1,2,i) = grad(idim,1,2,i) + d1*zdiff(2)
               grad(idim,2,2,i) = grad(idim,2,2,i) + d2*zdiff(2)
               grad(idim,3,2,i) = grad(idim,3,2,i) + d3*zdiff(2)
               grad(idim,1,3,i) = grad(idim,1,3,i) + d1*zdiff(3)
               grad(idim,2,3,i) = grad(idim,2,3,i) + d2*zdiff(3)
               grad(idim,3,3,i) = grad(idim,3,3,i) + d3*zdiff(3)

               d1 = -stoklet(idim,1,j)/(2*r3)
               d2 = -stoklet(idim,2,j)/(2*r3)
               d3 = -stoklet(idim,3,j)/(2*r3)               
               
               grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*d1
               grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(2)*d1
               grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(3)*d1
               grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(1)*d2
               grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*d2
               grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(3)*d2
               grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(1)*d3
               grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(2)*d3
               grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*d3

               pre(idim,i) = pre(idim,i) + pl*2

               
            enddo
 10         continue
         enddo
      enddo

      
      return
      end

      subroutine st3ddirectstokstrsg(nd,sources,stoklet,istress,
     1     strslet,strsvec,ns,targ,nt,pot,pre,grad,thresh)
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer(8)
c        number of densities
c     
c     nsource in: integer(8)
c        number of sources
c
c     source  in: double precision (3,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,3,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     istress in: integer(8)
c        stresslet computation flag
c           istress = 1   =>  include standard stresslet
c                               (type I)
c     
c           NOT YET IMPLEMENTED
c      
c           ifstress = 2   =>  include symmetric stresslet
c                                   (type II)
c           ifstress = 3   =>  include rotlet
c           ifstress = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c     strslet  in: double precision (nd,3,nsource) 
c        stresslet strengths (mu vectors above)
c
c     strsvec  in: double precision (nd,3,nsource)   
c        stresslet orientations (nu vectors above)
c      
c     ntarg   in: integer(8)
c        number of targs 
c
c     targ    in: double precision (3,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,3,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,3,3,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,stoklet,istress,strslet,strsvec,ns
cf2py intent(in) targ,nt,thresh
cf2py intent(out) pot,pre,grad

      integer(8) nd, ns, nt, istress
      real *8 sources(3,ns),targ(3,nt),strslet(nd,3,ns)
      real *8 strsvec(nd,3,ns),stoklet(nd,3,ns)
      real *8 pot(nd,3,nt),pre(nd,nt),grad(nd,3,3,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(3), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pv, dmu(3), dnu(3), temp, r, r2, r3, r5
      real *8 dmunu      
      real *8 threshsq

      integer(8) i, j, idim, l


      call st3ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1     pot,pre,grad,thresh)

      
      threshsq = thresh**2


      if (istress .ne. 1) goto 100

c     type I stresslet
      
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)
            zdiff(3) = targ(3,i)-sources(3,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
            if (r2 .lt. threshsq) goto 20

            r = sqrt(r2)
            r3 = r*r2
            r5 = r3*r2
            
            do idim = 1,nd

               dmu(1) = strslet(idim,1,j)
               dmu(2) = strslet(idim,2,j)
               dmu(3) = strslet(idim,3,j)
               dnu(1) = strsvec(idim,1,j)
               dnu(2) = strsvec(idim,2,j)
               dnu(3) = strsvec(idim,3,j)

               pl = zdiff(1)*dmu(1) + zdiff(2)*dmu(2) + zdiff(3)*dmu(3)
               pv = zdiff(1)*dnu(1) + zdiff(2)*dnu(2) + zdiff(3)*dnu(3)

               temp = -3.0d0*pl*pv/r5
               
               pot(idim,1,i) = pot(idim,1,i) + zdiff(1)*temp
               pot(idim,2,i) = pot(idim,2,i) + zdiff(2)*temp
               pot(idim,3,i) = pot(idim,3,i) + zdiff(3)*temp

               tempx1 = -3.0d0*(dmu(1)*pv + dnu(1)*pl -
     1              5.0d0*zdiff(1)*pl*pv/r2)/ r5
               tempx2 = -3.0d0*(dmu(2)*pv + dnu(2)*pl -
     1              5.0d0*zdiff(2)*pl*pv/r2)/ r5
               tempx3 = -3.0d0*(dmu(3)*pv + dnu(3)*pl -
     1              5.0d0*zdiff(3)*pl*pv/r2)/ r5

               grad(idim,1,1,i) = grad(idim,1,1,i) + temp
               grad(idim,1,1,i) = grad(idim,1,1,i) + zdiff(1)*tempx1
               grad(idim,2,1,i) = grad(idim,2,1,i) + zdiff(1)*tempx2
               grad(idim,3,1,i) = grad(idim,3,1,i) + zdiff(1)*tempx3
               
               grad(idim,2,2,i) = grad(idim,2,2,i) + temp
               grad(idim,1,2,i) = grad(idim,1,2,i) + zdiff(2)*tempx1
               grad(idim,2,2,i) = grad(idim,2,2,i) + zdiff(2)*tempx2
               grad(idim,3,2,i) = grad(idim,3,2,i) + zdiff(2)*tempx3
               
               grad(idim,3,3,i) = grad(idim,3,3,i) + temp
               grad(idim,1,3,i) = grad(idim,1,3,i) + zdiff(3)*tempx1
               grad(idim,2,3,i) = grad(idim,2,3,i) + zdiff(3)*tempx2
               grad(idim,3,3,i) = grad(idim,3,3,i) + zdiff(3)*tempx3

               dmunu = dmu(1)*dnu(1) + dmu(2)*dnu(2) + dmu(3)*dnu(3)
               
               pre(idim,i) = pre(idim,i) + 2.0d0*dmunu/r3-6.0d0*pl*pv/r5
            enddo
            
 20         continue
         enddo
      enddo
      
      return
      
 100  continue
      
      return
      end
      
