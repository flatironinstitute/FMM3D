c      This file contains the direct evaluation kernels for Laplace FMM
c       using the SCTL library
c
c      l3ddirectcp: direct calculation of potential for a collection
c                     of charge sources to a collection of targets
c 
c      l3ddirectcg: direct calculation of potential and gradients 
c                   for a collection of charge sources to a 
c                   collection of targets
c 
c      l3ddirectdp: direct calculation of potential for a collection
c                     of dipole sources to a collection of targets
c 
c      l3ddirectdg: direct calculation of potential and gradients 
c                   for a collection of dipole sources to a 
c                   collection of targets
c 
c      l3ddirectcdp: direct calculation of potential for a collection
c                     of charge and dipole sources to a collection 
c                     of targets
c 
c      l3ddirectdg: direct calculation of potential and gradients 
c                   for a collection of charge and dipole sources to 
c                   a collection of targets
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectcp(nd,sources,charge,ns,ztarg,nt,
     1            pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum  q_{j} /|x-x_{j}| 
c                        j
c                 
c      where q_{j} is the charge strength
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     sources:    source locations
C     charge :    charge strengths
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt)
      real *8 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
      call l3ddirectcp_cpp(nd,sources,charge,ns,ztarg,nt,
     1            pot,thresh)

      return
      end
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectcg(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} /|x-x_{j}| 
c                        j
c                 
c     grad(x) = grad(x) + Gradient(sum  q_{j} /|x-x_{j}|) 
c                                   j
c      where q_{j} is the charge strength
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     sources:    source locations
C     charge :    charge strengths
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad   :    updated gradient at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt)
      real *8 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
      call l3ddirectcg_cpp(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)


      return
      end
c
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectch(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} /|x-x_{j}|
c                        j
c
c     grad(x) = grad(x) + Gradient(sum  q_{j} /|x-x_{j}|)
c
c     hess(x) = hess(x) + Hessian(sum  q_{j} /|x-x_{j}|)
c                                   j
c      where q_{j} is the charge strength
c      If |r| < thresh
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain)
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     sources:    source locations
C     charge :    charge strengths
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg
c     grad   :    updated gradient at ztarg
c     hess   :    updated hessian at ztarg
c                 using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess
c
cc      calling sequence variables
c
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt)
      real *8 charge(nd,ns)
      real *8 pot(nd,nt),grad(nd,3,nt),hess(nd,6,nt)
      real *8 thresh

      call l3ddirectch_cpp(nd,sources,charge,ns,ztarg,nt,
     1            pot,grad,hess,thresh)

      return
      end
c
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectdp(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum   \nabla 1/|x-x_{j}| \cdot v_{j} 
c   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable 
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     dipvec :    dipole orientation vectors
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 pot(nd,nt)
      real *8 thresh
      
      call l3ddirectdp_cpp(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)

      return
      end
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectdg(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  d_{j} \nabla 1/|x-x_{j}| \cdot v_{j}
c                        j
c   
c     grad(x) = grad(x) + Gradient( sum  
c                                    j
c
c                            \nabla 1|/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     dipvec :    dipole orientation vector
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad   :    updated gradient at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      

      call l3ddirectdg_cpp(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3ddirectdh(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  d_{j} \nabla 1/|x-x_{j}| \cdot v_{j}
c                        j
c   
c     grad(x) = grad(x) + Gradient( sum  
c                                    j
c
c                            \nabla 1|/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c     hess(x) = hess(x) + Hessian( sum  
c                                    j
c
c                            \nabla 1|/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     dipvec :    dipole orientation vector
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad   :    updated gradient at ztarg 
c     hess   :    updated gradient at ztarg 
c                 using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 pot(nd,nt),grad(nd,3,nt),hess(nd,6,nt)
      real *8 thresh
      
      call l3ddirectdh_cpp(nd,sources,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3ddirectcdp(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum  q_{j} 1/|x-x_{j}| +  
c                        j
c
c                            \nabla 1/|x-x_{j}| \cdot v_{j}
c   
c      where q_{j} is the charge strength, 
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable 
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     charge :    charge strengths
C     dipvec :    dipole orientation vectors
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
      call l3ddirectcdp_cpp(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,thresh)


      return
      end
c
c
c
c
c
c
C***********************************************************************
      subroutine l3ddirectcdg(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} 1/|x-x_{j}| +  
c                        j
c
c                            \nabla 1/|x-x_{j}| \cdot v_{j}
c   
c     grad(x) = grad(x) + Gradient( sum  q_{j} 1/|x-x_{j}| +  
c                                    j
c
c                            \nabla 1/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where q_{j} is the charge strength, 
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     charge :    charge strengths
C     dipvec :    dipole orientation vector
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad   :    updated gradient at ztarg 
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      

      call l3ddirectcdg_cpp(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3ddirectcdh(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} 1/|x-x_{j}| +  
c                        j
c
c                            \nabla 1/|x-x_{j}| \cdot v_{j}
c   
c     grad(x) = grad(x) + Gradient( sum  q_{j} 1/|x-x_{j}| +  
c                                    j
c
c                            \nabla 1/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where q_{j} is the charge strength, 
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = |boxsize(0)|*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     sources:    source locations
C     charge :    charge strengths
C     dipvec :    dipole orientation vector
C     ns     :    number of sources
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     thresh :    threshold for updating potential,
c                 potential at target won't be updated if
c                 |t - s| <= thresh, where t is the target
c                 location and, and s is the source location 
c                 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    updated potential at ztarg 
c     grad   :    updated gradient at ztarg 
c     hess   :    updated Hessian at ztarg 
c                 using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      real *8 sources(3,ns),ztarg(3,nt),dipvec(nd,3,ns)
      real *8 charge(nd,ns)
      real *8 pot(nd,nt),grad(nd,3,nt),hess(nd,6,nt)
      real *8 thresh
      
      call l3ddirectcdh_cpp(nd,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)

      return
      end
