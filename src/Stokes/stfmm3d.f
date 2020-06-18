c
c     This file contains the Stokes FMM wrappers
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
c     P_j(x,y) = -r_j/r^3
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, (without the 1/4pi scaling) are
c     
c     T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5
c     PI_{jk} = -2 delta_{jk} + 6 r_j r_k/r^5      
c
      

      subroutine stfmm3d(nd, eps, 
     $                 nsource, source,
     $                 ifstoklet, stoklet, ifstrslet, strslet, strsvec,
     $                 ifppreg, pot, pre, grad, ntarg, targ, 
     $                 ifppregtarg, pottarg, pretarg, gradtarg)
c
c     Stokes FMM in R^{3}: evaluate all pairwise particle
c     interactions (ignoring self-interactions) and
c     interactions with targs.
c      
c     This routine computes sums of the form
c
c       u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c     Optionally, the associated pressure p(x) and gradient grad u(x)
c     are returned
c
c       p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
c
c       grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
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
c
c     TODO: implement other stresslet options
c------------------------------------------------------------------
      implicit none
      integer nd, ifstoklet, ifstrslet
      double precision eps
      integer nsource, ifppreg, ifppregtarg
      double precision source(3, *), targ(3, *)
      double precision stoklet(nd, 3, *), strslet(nd, 3, *)
      double precision strsvec(nd, 3, *)
      double precision pot(nd, 3, *), pre(nd,*), grad(nd, 3, 3, *)
      double precision pottarg(nd, 3, *), pretarg(nd,*),
     1     gradtarg(nd, 3, 3, *)      

c     local
      double precision, allocatable :: charge(:,:,:), dipvec(:,:,:,:)
      double precision, allocatable :: potl(:,:,:), gradl(:,:,:,:),
     1     hessl(:,:,:,:), pottargl(:,:,:), gradtargl(:,:,:,:),
     2     hesstargl(:,:,:,:)

      integer ndl, ifchargel, ifdipolel, ifpghl, ifpghtargl
      integer ndper, nsourcec, nsourced, nsp, ntp, nsg, ntg, nsh, nth
      

      ndper = 0
      
      if (ifstrslet .eq. 1 .or. ifstoklet .eq. 1) then
         ndper = 4
      endif

      ifdipolel = 0
      ifchargel = 0

      if (ifstoklet .eq. 1) ifchargel = 1
      if (ifstrslet .eq. 1) ifdipolel = 1
      
      if (ifstrslet .eq. 1 .or. ifstrslet .eq. 0
      ndper = 4
      ndl = ndper*nd

      ifpghl = 2
      ifpghltarg = ifpghl

c     allocate necessary arrays
      
      allocate(charge(nd,ndper,nsource),dipvec(nd,ndper,3,nsource),
     1     potl(nd,ndper,nsource),pottargl(nd,ndper,ntarg),
     2     gradl(nd,ndper,3,nsource),gradtargl(nd,ndper,3,ntarg),
     3     hessl(nd,ndper,6,nsource),hesstargl(nd,ndper,6,ntarg),
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

      do i = 1,nsource
         
         do l = 1,ndper
            do j = 1,nd
               charge(j,l,i) = 0
               dipvec(j,l,1,i) = 0
               dipvec(j,l,2,i) = 0
               dipvec(j,l,3,i) = 0
            enddo
         enddo

         do l = 1,3
            do j = 1,nd

               if (ifstoklet .eq. 1) then
                  charge(j,l,i) = charge(j,l,i) + stoklet(j,l,i)/2
               endif
               if (ifstrslet .eq. 1) then
                  dipvec(j,l,1,i) = dipvec(j,l,1,i) +
     1                 (strslet(j,l,i)*strsvec(j,1,i) +
     2                 strsvec(j,l,i)*strslet(j,1,i))/2
                  dipvec(j,l,2,i) = dipvec(j,l,2,i) +
     1                 (strslet(j,l,i)*strsvec(j,2,i) +
     2                 strsvec(j,l,i)*strslet(j,2,i))/2
                  dipvec(j,l,3,i) = dipvec(j,l,3,i) +
     1                 (strslet(j,l,i)*strsvec(j,3,i) +
     2                 strsvec(j,l,i)*strslet(j,3,i))/2
               endif
            enddo
         enddo

         l = 4

         do j = 1,nd
            if (ifstoklet .eq. 1) then
               
               pjl = (stoklet(j,1,i)*source(1,i) +
     1              stoklet(j,2,i)*source(2,i) +
     2              stoklet(j,3,i)*source(3,i))/2
               
               charge(j,l,i) = charge(j,l,i) + pjl
            endif
            if (ifstrslet .eq. 1) then
               
               pjv = (strsvec(j,1,i)*source(1,i) +
     1              strsvec(j,2,i)*source(2,i) +
     2              strsvec(j,3,i)*source(3,i))/2
               pjl = (strslet(j,1,i)*source(1,i) +
     1              strslet(j,2,i)*source(2,i) +
     2              strslet(j,3,i)*source(3,i))/2
               
               dipvec(j,l,1,i) = dipvec(j,l,1,i) +
     1              strslet(j,1,i)*pjv + strsvec(j,1,i)*pjl
               dipvec(j,l,2,i) = dipvec(j,l,2,i) +
     1              strslet(j,2,i)*pjv + strsvec(j,2,i)*pjl
               dipvec(j,l,3,i) = dipvec(j,l,3,i) +
     1              strslet(j,3,i)*pjv + strsvec(j,3,i)*pjl
            endif
         enddo
                  
      enddo


c     call Laplace FMM

      call lfmm3d(ndl,eps,nsource,source,ifchargel,charge,
     1     ifdipolel,dipvec,ifpghl,potl,gradl,hessl,ntarg,
     2     targ,ifpghtargl,pottargl,gradtargl,hesstargl)


      npt = ntarg + nsource

c     unpack stacked Laplace FMM calls
      
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
                  pl = potl(j,l,ii)
                  gl(1) = gradl(j,l,1,ii)
                  gl(2) = gradl(j,l,2,ii)
                  gl(3) = gradl(j,l,3,ii)
                  hl(1) = hessl(j,l,1,ii)
                  hl(2) = hessl(j,l,2,ii)
                  hl(3) = hessl(j,l,3,ii)
                  hl(4) = hessl(j,l,4,ii)
                  hl(5) = hessl(j,l,5,ii)
                  hl(6) = hessl(j,l,6,ii)            
               else
                  ifppreg1 = ifppregtarg                  
                  ii = i
                  pt(1) = targ(1,ii)
                  pt(2) = targ(2,ii)
                  pt(3) = targ(3,ii)
                  pl = pottargl(j,l,ii)
                  gl(1) = gradtargl(j,l,1,ii)
                  gl(2) = gradtargl(j,l,2,ii)
                  gl(3) = gradtargl(j,l,3,ii)
                  hl(1) = hesstargl(j,l,1,ii)
                  hl(2) = hesstargl(j,l,2,ii)
                  hl(3) = hesstargl(j,l,3,ii)
                  hl(4) = hesstargl(j,l,4,ii)
                  hl(5) = hesstargl(j,l,5,ii)
                  hl(6) = hesstargl(j,l,6,ii)            
               endif


               if (l .ge. 1 .and. l .le. 3) then
                  
                  vel(l) = vel(l) + pl
                  vel(1) = vel(1) - pt(l)*gl(1)
                  vel(2) = vel(2) - pt(l)*gl(2)
                  vel(3) = vel(3) - pt(l)*gl(3)
                  press = press + gl(l)*2
                  if (ifppreg1 .eq. 3) then
                     velgrad(1,l) =  velgrad(1,l) + gl(1)
                     velgrad(2,l) =  velgrad(2,l) + gl(2)
                     velgrad(3,l) =  velgrad(3,l) + gl(3)

                     velgrad(l,1) = velgrad(l,1) - gl(1)
                     velgrad(l,2) = velgrad(l,2) - gl(2)
                     velgrad(l,3) = velgrad(l,3) - gl(3)                     

c     confirm hessian ordering convention...
                     velgrad(1,1) = velgrad(1,1) - pt(l)*hl(1)
                     velgrad(2,1) = velgrad(2,1) - pt(l)*hl(2)
                     velgrad(3,1) = velgrad(3,1) - pt(l)*hl(3)
                     velgrad(1,2) = velgrad(1,2) - pt(l)*hl(2)
                     velgrad(2,2) = velgrad(2,2) - pt(l)*hl(4)
                     velgrad(3,2) = velgrad(3,2) - pt(l)*hl(5)
                     velgrad(1,3) = velgrad(1,1) - pt(l)*hl(3)
                     velgrad(2,3) = velgrad(2,1) - pt(l)*hl(5)
                     velgrad(3,3) = velgrad(3,1) - pt(l)*hl(6)
                  endif

               else if (l .eq. 4) then

                  vel(1) = vel(1) + gl(1)
                  vel(2) = vel(2) + gl(2)
                  vel(3) = vel(3) + gl(3)                  
                  if (ifppreg1 .eq. 3) then
c     confirm hessian ordering convention...
                     velgrad(1,1) = velgrad(1,1) + hl(1)
                     velgrad(2,1) = velgrad(2,1) + hl(2)
                     velgrad(3,1) = velgrad(3,1) + hl(3)
                     velgrad(1,2) = velgrad(1,2) + hl(2)
                     velgrad(2,2) = velgrad(2,2) + hl(4)
                     velgrad(3,2) = velgrad(3,2) + hl(5)
                     velgrad(1,3) = velgrad(1,1) + hl(3)
                     velgrad(2,3) = velgrad(2,1) + hl(5)
                     velgrad(3,3) = velgrad(3,1) + hl(6)
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
      
      return
      end

