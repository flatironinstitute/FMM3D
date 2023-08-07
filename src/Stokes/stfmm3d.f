c
c     This file contains the Stokes FMM wrappers
c
c**********************************************************************
c      
c     We take the following conventions for the Stokes kernels.
c
c     1) The dynamic viscosity (mu) is assumed to be 1.
c     2) All kernels are a factor 4pi larger than standard definitions
c        (this is for historical reasons).
c     Thus, in general, divide the velocity (potential or grad) outputs
c     by 4pi.mu, and pressure by 4pi, to recover standard definitions.
c
c     For a source y and target x, let r_i = x_i-y_i     (note sign)
c     and let r = sqrt(r_1^2 + r_2^2 + r_3^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     we define as
c
c     G_{ij}(x,y) = ( delta_{ij}/r  + r_i r_j / r^3 )/2
c     P_j(x,y) = r_j/r^3
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, we define as
c     
c     T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5
c     PI_{jk}(x,y) = -2 delta_{jk}/r^3 + 6 r_j r_k/r^5      
      

      subroutine stfmm3d(nd, eps, 
     $                 nsource, source,
     $                 ifstoklet, stoklet, ifstrslet, strslet, strsvec,
     $                 ifppreg, pot, pre, grad, ntarg, targ, 
     $                 ifppregtarg, pottarg, pretarg, gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source
cf2py  intent(in) ifstoklet,stoklet
cf2py  intent(in) ifstrslet,strslet,strsvec
cf2py  intent(in) ifppreg,ifppregtarg
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,pre,grad
cf2py  intent(out) pottarg,pretarg,gradtarg
cf2py  intent(out) ier
c
c     Stokes FMM in R^{3}: evaluate all pairwise particle
c     interactions (ignoring self-interactions) and
c     interactions with targs.
c      
c     This routine computes the sum for the velocity vector,
c
c       u_i(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c     for x at all of the target locations, and i=1,2,3,
c     where sigma^{(m)} is the Stokeslet force, mu^{(m)} is the
c     stresslet strength, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     Repeated indices are taken as summed over 1,2,3, ie, Einstein
c     convention. For x a source point, the self-interaction in the
c     sum is omitted. 
c
c     Optionally, the associated pressure p(x) and 3x3 gradient tensor
c     grad u(x) are returned
c
c       p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m PI_{jk}(x,y{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c       grad_l u_i(x) = grad_l [sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c
c     Note that these two may be combined to get the stress tensor.
c
c-----------------------------------------------------------------------
c     INPUT PARAMETERS:
c     
c   nd:    in: integer
c              number of densities
c   
c   eps:   in: double precision
c              requested precision
c
c   nsource in: integer  
c               number of sources
c
c   source  in: double precision (3,nsource)
c               source(k,j) is the kth component of the jth
c               source locations
c
c   ifstoklet  in: integer  
c               Stokeslet charge computation flag
c               ifstoklet = 1   =>  include Stokeslet contribution
c                                   otherwise do not
c 
c   stoklet in: double precision (nd,3,nsource) 
c               Stokeslet charge strengths (sigma vectors above)
c
c   ifstrslet in: integer
c               stresslet computation flag
c               ifstrslet = 1   =>  include standard stresslet
c                                   (type I)
c
c            NOT YET IMPLEMENTED
c      
c               ifstrslet = 2   =>  include symmetric stresslet
c                                   (type II)
c               ifstrslet = 3   =>  include rotlet
c               ifstrslet = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c   strslet  in: double precision (nd,3,nsource) 
c               stresslet strengths (mu vectors above)
c
c   strsvec  in: double precision (nd,3,nsource)   
c               stresslet orientations (nu vectors above)
c
c     ifppreg    in: integer      
c               flag for evaluating potential, gradient, and pressure
c               at the sources
c               ifppreg = 1, only potential
c               ifppreg = 2, potential and pressure
c         GRADIENT NOT IMPLEMENTED
c               ifppreg = 3, potential, pressure, and gradient 
c      
c   ntarg   in: integer  
c              number of targs 
c
c   targ    in: double precision (3,ntarg)
c             targ(k,j) is the kth component of the jth
c             targ location
c      
c   ifppregtarg in: integer
c                flag for evaluating potential, gradient, and pressure
c                at the targets
c                ifppregtarg = 1, only potential
c                ifppregtarg = 2, potential and pressure
c                ifppregtarg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT parameters:
c
c   pot   out: double precision(nd,3,nsource) 
c           velocity at the source locations
c      
c   pre   out: double precision(nd,nsource)
c           pressure at the source locations
c      
c         GRADIENT NOT IMPLEMENTED
c   grad   out: double precision(nd,3,3,nsource) 
c              gradient of velocity at the source locations
c              grad(l,i,j,k) is the ith component of the
c              gradient of the jth component of the velocity
c              for the lth density at the kth source location
c     
c   pottarg   out: double precision(nd,3,ntarg) 
c               velocity at the targets
c      
c   pretarg   out: double precision(nd,ntarg)
c               pressure at the targets
c      
c   gradtarg   out: double precision(nd,3,3,ntarg) 
c               gradient of velocity at the targets
c               gradtarg(l,i,j,k) is the ith component of the
c               gradient of the jth component of the velocity
c               for the lth density at the kth target
c     ier     out: integer
c               error flag
c
c     TODO: implement other stresslet options and gradient
c------------------------------------------------------------------
      implicit none
      integer nd, ifstoklet, ifstrslet, ntarg
      double precision eps
      integer nsource, ifppreg, ifppregtarg
      double precision source(3, nsource), targ(3, ntarg)
      double precision stoklet(nd, 3, nsource), strslet(nd, 3, nsource)
      double precision strsvec(nd, 3, nsource)
      double precision pot(nd, 3, nsource), pre(nd,nsource)
      double precision grad(nd, 3, 3, nsource)
      double precision pottarg(nd, 3, ntarg), pretarg(nd,ntarg),
     1     gradtarg(nd, 3, 3, ntarg)      

c     local
      double precision, allocatable :: charge(:,:,:), dipvec(:,:,:,:)
      double precision, allocatable :: potl(:,:,:), gradl(:,:,:,:),
     1     hessl(:,:,:,:), pottargl(:,:,:), gradtargl(:,:,:,:),
     2     hesstargl(:,:,:,:)
      double precision :: pt(3), gl(3), hl(6), vel(3), velgrad(3,3)
      double precision :: press, pl, pv, dmu(3), dnu(3), sigma(3)

      integer ndl, ifchargel, ifdipolel, ifpghl, ifpghtargl
      integer ndper

      integer i, j, ii, ifppreg1, l, npt, ier,iper
      

      ndper = 0
      
      if (ifstrslet .eq. 1 .or. ifstoklet .eq. 1) then
         ndper = 4
      endif

      ifdipolel = 0
      ifchargel = 0

      if (ifstoklet .eq. 1) ifchargel = 1
      if (ifstrslet .eq. 1) ifdipolel = 1
      
      ndper = 4
      ndl = ndper*nd

      ifpghl = 3
      ifpghtargl = 3

c     allocate necessary arrays
      
      allocate(charge(ndper,nd,nsource),dipvec(ndper,nd,3,nsource),
     1     potl(ndper,nd,nsource),pottargl(ndper,nd,ntarg),
     2     gradl(ndper,nd,3,nsource),gradtargl(ndper,nd,3,ntarg),
     3     hessl(ndper,nd,6,nsource),hesstargl(ndper,nd,6,ntarg),
     4     stat=ier)
      if(ier .ne. 0) then
         print *, "In stfmm3d: cannot allocate Laplace call storage"
         print *, "ndper =",ndper
         print *, "nd =",nd
         print *, "nsource =",nsource
         print *, "ntarg =",ntarg        
         stop
      endif

c     set-up appropriate vector charge and dipole arrays

c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,sigma,dmu,dnu)
c$OMP$ PRIVATE(pl,pv)      
      do i = 1,nsource

         do j = 1,nd
            do l = 1,ndper
               charge(l,j,i) = 0
               dipvec(l,j,1,i) = 0
               dipvec(l,j,2,i) = 0
               dipvec(l,j,3,i) = 0
            enddo
         enddo

         do j = 1,nd
            if (ifstoklet .eq. 1) then
               sigma(1) = stoklet(j,1,i)
               sigma(2) = stoklet(j,2,i)
               sigma(3) = stoklet(j,3,i)
            endif
            if (ifstrslet .ge. 1) then
               dmu(1) = strslet(j,1,i)
               dmu(2) = strslet(j,2,i)
               dmu(3) = strslet(j,3,i)
               dnu(1) = strsvec(j,1,i)
               dnu(2) = strsvec(j,2,i)
               dnu(3) = strsvec(j,3,i)
            endif

            do l = 1,3
               
               if (ifstoklet .eq. 1) then
                  charge(l,j,i) = charge(l,j,i) + sigma(l)/2
               endif
               if (ifstrslet .eq. 1) then
                  dipvec(l,j,1,i) = dipvec(l,j,1,i) - (dmu(l)*dnu(1) + 
     1                 dmu(1)*dnu(l))/2
                  dipvec(l,j,2,i) = dipvec(l,j,2,i) - (dmu(l)*dnu(2) + 
     1                 dmu(2)*dnu(l))/2
                  dipvec(l,j,3,i) = dipvec(l,j,3,i) - (dmu(l)*dnu(3) + 
     1                 dmu(3)*dnu(l))/2
               endif
            enddo
            
            l = 4
            
            if (ifstoklet .eq. 1) then
               pl = sigma(1)*source(1,i) + sigma(2)*source(2,i) +
     1              sigma(3)*source(3,i)
               charge(l,j,i) = charge(l,j,i) + pl/2
            endif
            if (ifstrslet .eq. 1) then
               pl = dmu(1)*source(1,i) + dmu(2)*source(2,i) +
     1              dmu(3)*source(3,i)
               pv = dnu(1)*source(1,i) + dnu(2)*source(2,i) +
     1              dnu(3)*source(3,i)
               
               dipvec(l,j,1,i) = dipvec(l,j,1,i) -
     1              (dmu(1)*pv + dnu(1)*pl)/2
               dipvec(l,j,2,i) = dipvec(l,j,2,i) -
     1              (dmu(2)*pv + dnu(2)*pl)/2
               dipvec(l,j,3,i) = dipvec(l,j,3,i) -
     1              (dmu(3)*pv + dnu(3)*pl)/2
            endif
            
         enddo
         
      enddo
c$OMP END PARALLEL DO      


c     call Laplace FMM

      iper = 0
      ier = 0
      call lfmm3d(ndl,eps,nsource,source,ifchargel,charge,
     1     ifdipolel,dipvec,iper,ifpghl,potl,gradl,hessl,ntarg,
     2     targ,ifpghtargl,pottargl,gradtargl,hesstargl,ier)


      npt = ntarg + nsource

c     unpack stacked Laplace FMM calls

c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,ifppreg1,ii)
c$OMP$ PRIVATE(pt,pl,gl,hl,vel,velgrad,press)      
      do i = 1,npt

         do j = 1,nd

            vel(1) = 0
            vel(2) = 0
            vel(3) = 0
            velgrad(1,1) = 0
            velgrad(2,1) = 0
            velgrad(3,1) = 0
            velgrad(1,2) = 0
            velgrad(2,2) = 0
            velgrad(3,2) = 0
            velgrad(1,3) = 0
            velgrad(2,3) = 0
            velgrad(3,3) = 0
            press = 0
            
            do l = 1,ndper
               
               if (i .gt. ntarg) then
                  ifppreg1 = ifppreg
                  ii = i-ntarg
                  pt(1) = source(1,ii)
                  pt(2) = source(2,ii)
                  pt(3) = source(3,ii)
                  pl = potl(l,j,ii)
                  gl(1) = gradl(l,j,1,ii)
                  gl(2) = gradl(l,j,2,ii)
                  gl(3) = gradl(l,j,3,ii)
                  hl(1) = hessl(l,j,1,ii)
                  hl(2) = hessl(l,j,2,ii)
                  hl(3) = hessl(l,j,3,ii)
                  hl(4) = hessl(l,j,4,ii)
                  hl(5) = hessl(l,j,5,ii)
                  hl(6) = hessl(l,j,6,ii)            
               else
                  ifppreg1 = ifppregtarg                  
                  ii = i
                  pt(1) = targ(1,ii)
                  pt(2) = targ(2,ii)
                  pt(3) = targ(3,ii)
                  pl = pottargl(l,j,ii)
                  gl(1) = gradtargl(l,j,1,ii)
                  gl(2) = gradtargl(l,j,2,ii)
                  gl(3) = gradtargl(l,j,3,ii)
                  hl(1) = hesstargl(l,j,1,ii)
                  hl(2) = hesstargl(l,j,2,ii)
                  hl(3) = hesstargl(l,j,3,ii)
                  hl(4) = hesstargl(l,j,4,ii)
                  hl(5) = hesstargl(l,j,5,ii)
                  hl(6) = hesstargl(l,j,6,ii)            
               endif

               
               if (l .ge. 1 .and. l .le. 3) then
                  
                  vel(l) = vel(l) + pl
                  vel(1) = vel(1) - pt(l)*gl(1)
                  vel(2) = vel(2) - pt(l)*gl(2)
                  vel(3) = vel(3) - pt(l)*gl(3)
                  press = press - gl(l)*2
                  if (ifppreg1 .eq. 3) then
                     velgrad(1,l) =  velgrad(1,l) + gl(1)
                     velgrad(2,l) =  velgrad(2,l) + gl(2)
                     velgrad(3,l) =  velgrad(3,l) + gl(3)

                     velgrad(l,1) = velgrad(l,1) - gl(1)
                     velgrad(l,2) = velgrad(l,2) - gl(2)
                     velgrad(l,3) = velgrad(l,3) - gl(3)                     

c     confirm hessian ordering convention...
                     velgrad(1,1) = velgrad(1,1) - pt(l)*hl(1)
                     velgrad(2,1) = velgrad(2,1) - pt(l)*hl(4)
                     velgrad(3,1) = velgrad(3,1) - pt(l)*hl(5)
                     velgrad(1,2) = velgrad(1,2) - pt(l)*hl(4)
                     velgrad(2,2) = velgrad(2,2) - pt(l)*hl(2)
                     velgrad(3,2) = velgrad(3,2) - pt(l)*hl(6)
                     velgrad(1,3) = velgrad(1,3) - pt(l)*hl(5)
                     velgrad(2,3) = velgrad(2,3) - pt(l)*hl(6)
                     velgrad(3,3) = velgrad(3,3) - pt(l)*hl(3)
                  endif

               else if (l .eq. 4) then

                  vel(1) = vel(1) + gl(1)
                  vel(2) = vel(2) + gl(2)
                  vel(3) = vel(3) + gl(3)                  
                  if (ifppreg1 .eq. 3) then
c     confirm hessian ordering convention...
                     velgrad(1,1) = velgrad(1,1) + hl(1)
                     velgrad(2,1) = velgrad(2,1) + hl(4)
                     velgrad(3,1) = velgrad(3,1) + hl(5)
                     velgrad(1,2) = velgrad(1,2) + hl(4)
                     velgrad(2,2) = velgrad(2,2) + hl(2)
                     velgrad(3,2) = velgrad(3,2) + hl(6)
                     velgrad(1,3) = velgrad(1,3) + hl(5)
                     velgrad(2,3) = velgrad(2,3) + hl(6)
                     velgrad(3,3) = velgrad(3,3) + hl(3)
                  endif
                  
               endif
            enddo

            if (i .gt. ntarg) then
               if (ifppreg1 .ge. 1) then
                  pot(j,1,ii) = vel(1)
                  pot(j,2,ii) = vel(2)
                  pot(j,3,ii) = vel(3)
               endif
               if (ifppreg1 .ge. 2) then
                  pre(j,ii) = press
               endif
               if (ifppreg1 .ge. 3) then
                  grad(j,1,1,ii) = velgrad(1,1)
                  grad(j,2,1,ii) = velgrad(2,1)
                  grad(j,3,1,ii) = velgrad(3,1)
                  grad(j,1,2,ii) = velgrad(1,2)
                  grad(j,2,2,ii) = velgrad(2,2)
                  grad(j,3,2,ii) = velgrad(3,2)
                  grad(j,1,3,ii) = velgrad(1,3)
                  grad(j,2,3,ii) = velgrad(2,3)
                  grad(j,3,3,ii) = velgrad(3,3)
               endif
            else
               if (ifppreg1 .ge. 1) then
                  pottarg(j,1,ii) = vel(1)
                  pottarg(j,2,ii) = vel(2)
                  pottarg(j,3,ii) = vel(3)
               endif
               if (ifppreg1 .ge. 2) then
                  pretarg(j,ii) = press
               endif
               if (ifppreg1 .ge. 3) then
                  gradtarg(j,1,1,ii) = velgrad(1,1)
                  gradtarg(j,2,1,ii) = velgrad(2,1)
                  gradtarg(j,3,1,ii) = velgrad(3,1)
                  gradtarg(j,1,2,ii) = velgrad(1,2)
                  gradtarg(j,2,2,ii) = velgrad(2,2)
                  gradtarg(j,3,2,ii) = velgrad(3,2)
                  gradtarg(j,1,3,ii) = velgrad(1,3)
                  gradtarg(j,2,3,ii) = velgrad(2,3)
                  gradtarg(j,3,3,ii) = velgrad(3,3)
               endif
            endif               
         enddo
      enddo
c$OMP END PARALLEL DO
      
      return
      end

