c      This file contains the direct evaluation kernels for Helmholtz FMM
c       using the sctl library
c
c      h3ddirectcp:  direct calculation of potential for a collection 
c                         of charge sources at a collection of targets
c
c      h3ddirectcg:  direct calculation of potential and gradient 
c                         for a collection of charge sources at a 
c                         collection of targets
c
c      h3ddirectdp:  direct calculation of potential for a collection 
c                         of dipole sources at a collection of targets
c
c      h3ddirectdg:  direct calculation of potential and gradient 
c                         for a collection of dipole sources at a 
c                         collection of targets
c
c      h3ddirectcdp:  direct calculation of potential for a 
c                         collection of charge and dipole sources at 
c                         a collection of targets
c
c      h3ddirectcdg:  direct calculation of potential 
c                         and gradient for a collection 
c                         of charge and dipole sources at a 
c                         collection of targets
c
c
c
c
c
C***********************************************************************
      subroutine h3ddirectcp(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| 
c                        j
c                 
c      where q_{j} is the charge strength
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
      call h3ddirectcp_cpp(nd,zk,sources,charge,ns,ztarg,nt,
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
      subroutine h3ddirectcg(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| 
c                        j
c                 
c     grad(x) = grad(x) + Gradient(sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}|) 
c                                   j
c      where q_{j} is the charge strength
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
      call h3ddirectcg_cpp(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)


      return
      end
c
c
c
c
c
c
C***********************************************************************
      subroutine h3ddirectch(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| 
c                        j
c                 
c     grad(x)=grad(x)+Gradient(sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}|) 
c                               j
c     hess(x)=hess(x)+Hessian(sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}|) 
c                              j
c      where q_{j} is the charge strength
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge densities
c     zk     :    Helmholtz parameter
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
c     hess   :    updated Hessian at ztarg 
c                 using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,zk,sources,charge,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh

            call h3ddirectch_cpp(nd,zk,sources,charge,ns,ztarg,nt,
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
      subroutine h3ddirectdp(nd,zk,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum   \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j} 
c                        j
c
c                            
c   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable 
c      If |r| < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt)
      real *8 thresh
      
      call h3ddirectdp_cpp(nd,zk,sources,
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
      subroutine h3ddirectdg(nd,zk,sources,dipvec,ns,ztarg,nt,pot,
     1   grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                        j
c
c                            
c   
c     grad(x) = grad(x) + Gradient( sum  
c                                    j
c
c                            \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
      call h3ddirectdg_cpp(nd,zk,sources,dipvec,ns,ztarg,nt,pot,
     1   grad,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3ddirectdh(nd,zk,sources,dipvec,ns,ztarg,nt,pot,
     1   grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x)=pot(x)+sum d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                    j
c
c     grad(x)=grad(x)+Gradient( sum  
c                                j
c
c                        \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                        )
c                                   
c     hess(x)=hess(x)+Hessian( sum  
c                               j
c
c                       \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                       )
c                                   
c      where v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
c     hess   :    updated Hessian at ztarg 
c                 using ordering (dxx,dyy,dzz,dxy,dxz,dyz)
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,zk,sources,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess

c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh
      
      call h3ddirectdh_cpp(nd,zk,sources,dipvec,ns,ztarg,nt,pot,
     1   grad,hess,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3ddirectcdp(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| +  
c                        j
c
c                            \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c   
c      where q_{j} is the charge strength, 
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable 
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot
c
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
      call h3ddirectcdp_cpp(nd,zk,sources,charge,
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
      subroutine h3ddirectcdg(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| +  
c                        j
c
c                            \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c   
c     grad(x) = grad(x) + Gradient( sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| +  
c                                    j
c
c                            d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                            )
c                                   
c      where q_{j} is the charge strength
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      

      call h3ddirectcdg_cpp(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,thresh)

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3ddirectcdh(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x) = pot(x) + sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| +  
c                        j
c
c                         \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c   
c     grad(x)=grad(x)+Gradient( sumq_{j} e^{i k |x-x_{j}|}/|x-x_{j}|+  
c                                j
c
c                   d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                   )
c                                   
c     hess(x)=hess(x)+Hessian( sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| +  
c                               j
c
c                   d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                   )
c                                   
c      where q_{j} is the charge strength
c      and v_{j} is the dipole orientation vector, 
c      \nabla denotes the gradient is with respect to the x_{j} 
c      variable, and Gradient denotes the gradient with respect to
c      the x variable
c      If r < thresh 
c          then the subroutine does not update the potential
c          (recommended value = boxsize(0)*machine precision
c           for boxsize(0) is the size of the computational domain) 
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :    number of charge and dipole densities
c     zk     :    Helmholtz parameter
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
cf2py intent(in) nd,zk,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad,hess
c
cc      calling sequence variables
c  
      integer ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh
      
      call h3ddirectcdh_cpp(nd,zk,sources,charge,
     1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)

      return
      end
c
c
c
