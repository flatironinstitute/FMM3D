c  This file contains the Stokes FMM wrappers

      subroutine stfmm3d(nd, eps, 
     $                 nsource, source,
     $                 ifstoklet, stoklet, ifdoublet, doublet, doubvec,
     $                 ifppreg, pot, pre, grad, ntarg, targ, 
     $                 ifppregtarg, pottarg, pretarg, gradtarg)
c
c     Stokes FMM in R^{3}: evaluate all pairwise particle
c     interactions (ignoring self-interactions) and
c     interactions with targs.
c
c     We use (r_i r_j)/(2r^3) + delta_ij/(2r) for the Stokes
c     Green's function (Stokeslet), without the
c     1/(4 \pi) scaling.
c
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
c               Stokeslet charge strengths
c
c   ifdoublet in: integer
c               higher order charge computation flag
c               ifdoublet = 1   =>  include standard stresslet
c                                   (type I)     
c               ifdoublet = 2   =>  include symmetric stresslet
c                                   (type II)
c               ifdoublet = 3   =>  include rotlet
c               ifdoublet = 4   =>  include Stokes doublet
c                                   otherwise do not
c
c   doublet  in: double precision (nd,3,nsource) 
c               higher order charge strengths
c
c   doubvec  in: double precision (nd,3,nsource)   
c               higher order charge orientations
c
c   ifppreg    in: integer
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
c     
c   pottarg   out: double precision(nd,3,nsource) 
c               velocity at the targets
c      
c   pretarg   out: double precision(nd,nsource)
c               pressure at the targets
c      
c   gradtarg   out: double precision(nd,3,3,nsource) 
c               gradient of velocity at the targets
c     
c------------------------------------------------------------------
      implicit none
      integer nd
      double precision eps
      integer nsource, ifppreg, ifppregtarg
      double precision source(3, *), targ(3, *)
      double precision stoklet(nd, 3, *), doublet(nd, 3, *)
      double precision doubvec(nd, 3, *)
      double precision pot(nd, 3, *), pre(nd,*), grad(nd, 3, 3, *)
      double precision pottarg(nd, 3, *), pretarg(nd,*),
     1     gradtarg(nd, 3, 3, *)      


      return
      end

