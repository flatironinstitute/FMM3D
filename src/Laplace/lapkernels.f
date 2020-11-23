c      This file contains the direct evaluation kernels for Laplace FMM
c
c      l3ddirectcp: direct calculation of potential for a collection
c                     of charge sources to a collection of targets
c 
c      l3ddirectcg: direct calculation of potential and gradients 
c                   for a collection of charge sources to a 
c                   collection of targets
c 
c      l3ddirectch: direct calculation of potential and gradients 
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
c      l3ddirectdh: direct calculation of potential and gradients 
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
c      l3ddirectdh: direct calculation of potential and gradients 
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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,ztmp,threshsq
      integer i,j,idim


      threshsq = thresh**2
      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          ztmp = 1.0d0/sqrt(dd)
          do idim=1,nd
            pot(idim,i) = pot(idim,i) + charge(idim,j)*ztmp
          enddo
 1000     continue
        enddo
      enddo


      return
      end
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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,cd,cd1,ztmp1,ztmp2,ztmp3
      real *8 threshsq
      integer i,j,idim


      threshsq = thresh**2 
      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000
          cd = 1/sqrt(dd)
          cd1 = -cd/dd
          ztmp1 = cd1*zdiff(1)
          ztmp2 = cd1*zdiff(2)
          ztmp3 = cd1*zdiff(3)
          do idim=1,nd
            pot(idim,i) = pot(idim,i) + cd*charge(idim,j)
            grad(idim,1,i) = grad(idim,1,i) + ztmp1*charge(idim,j)
            grad(idim,2,i) = grad(idim,2,i) + ztmp2*charge(idim,j)
            grad(idim,3,i) = grad(idim,3,i) + ztmp3*charge(idim,j)
          enddo
 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,cd,cd1,cd2,ztmp1,ztmp2,ztmp3
      real *8 htmp1,htmp2,htmp3,htmp4,htmp5,htmp6
      real *8 threshsq
      integer i,j,idim


      threshsq = thresh**2 
      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000
          cd = 1/sqrt(dd)
          cd1 = -cd/dd
          cd2 = -cd1/dd
          ztmp1 = cd1*zdiff(1)
          ztmp2 = cd1*zdiff(2)
          ztmp3 = cd1*zdiff(3)
          htmp1 = cd2*(3*zdiff(1)*zdiff(1)-dd)
          htmp2 = cd2*(3*zdiff(2)*zdiff(2)-dd)
          htmp3 = cd2*(3*zdiff(3)*zdiff(3)-dd)
          htmp4 = cd2*3*zdiff(1)*zdiff(2)
          htmp5 = cd2*3*zdiff(1)*zdiff(3)
          htmp6 = cd2*3*zdiff(2)*zdiff(3)
          do idim=1,nd
            pot(idim,i) = pot(idim,i) + cd*charge(idim,j)
            grad(idim,1,i) = grad(idim,1,i) + ztmp1*charge(idim,j)
            grad(idim,2,i) = grad(idim,2,i) + ztmp2*charge(idim,j)
            grad(idim,3,i) = grad(idim,3,i) + ztmp3*charge(idim,j)
            hess(idim,1,i) = hess(idim,1,i) + htmp1*charge(idim,j)
            hess(idim,2,i) = hess(idim,2,i) + htmp2*charge(idim,j)
            hess(idim,3,i) = hess(idim,3,i) + htmp3*charge(idim,j)
            hess(idim,4,i) = hess(idim,4,i) + htmp4*charge(idim,j)
            hess(idim,5,i) = hess(idim,5,i) + htmp5*charge(idim,j)
            hess(idim,6,i) = hess(idim,6,i) + htmp6*charge(idim,j)
          enddo
 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,cd,dotprod
      real *8 threshsq
      integer i,j,idim

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          cd = 1/sqrt(dd)/dd

          do idim=1,nd
            dotprod = zdiff(1)*dipvec(idim,1,j) + 
     1          zdiff(2)*dipvec(idim,2,j)+
     1          zdiff(3)*dipvec(idim,3,j)
            pot(idim,i) = pot(idim,i) + cd*dotprod
          enddo

 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2,dotprod
      real *8 cd,cd2,cd3,cd4
      real *8 threshsq
      integer i,j,idim

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd
          dinv = sqrt(dinv2)
          cd = dinv
          cd2 = -cd*dinv2
          cd3 = -3*cd*dinv2*dinv2

          do idim=1,nd
          
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod

            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j))
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
          enddo
 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2,dotprod
      real *8 cd,cd2,cd3,cd4,cd5,dx,dy,dz
      real *8 threshsq
      integer i,j,idim

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd
          dinv = sqrt(dinv2)
          cd = dinv
          cd2 = -cd*dinv2
          cd3 = -3*cd*dinv2*dinv2
          cd5 = -cd3
c
          dx = zdiff(1)*dinv
          dy = zdiff(2)*dinv
          dz = zdiff(3)*dinv

          do idim=1,nd
          
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod

            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j))
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
            hess(idim,1,i)= hess(idim,1,i) + 
     1      (dotprod*(5*dx*dx-1)-
     2      (dipvec(idim,1,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(1)))*cd5
            hess(idim,2,i)= hess(idim,2,i) + 
     1      (dotprod*(5*dy*dy-1)-
     2      (dipvec(idim,2,j)*zdiff(2)+dipvec(idim,2,j)*zdiff(2)))*cd5
            hess(idim,3,i)= hess(idim,3,i) + 
     1      (dotprod*(5*dz*dz-1)-
     2      (dipvec(idim,3,j)*zdiff(3)+dipvec(idim,3,j)*zdiff(3)))*cd5
            hess(idim,4,i)= hess(idim,4,i) + 
     1      (dotprod*(5*dx*dy)-
     2      (dipvec(idim,2,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(2)))*cd5
            hess(idim,5,i)= hess(idim,5,i) + 
     1      (dotprod*(5*dx*dz)-
     2      (dipvec(idim,3,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(3)))*cd5
            hess(idim,6,i)= hess(idim,6,i) + 
     1      (dotprod*(5*dy*dz)-
     2      (dipvec(idim,3,j)*zdiff(2)+dipvec(idim,2,j)*zdiff(3)))*cd5
          enddo
 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv2,dotprod,cd,cd1
      integer i,j,idim
      real *8 threshsq

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd 
          cd = sqrt(dinv2)
          cd1 = cd*dinv2

          do idim=1,nd
            pot(idim,i) = pot(idim,i) + charge(idim,j)*cd

            dotprod = zdiff(1)*dipvec(idim,1,j) + 
     1          zdiff(2)*dipvec(idim,2,j)+
     1          zdiff(3)*dipvec(idim,3,j)
            pot(idim,i) = pot(idim,i) + cd1*dotprod
          enddo

 1000     continue
        enddo
      enddo


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
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2,dotprod,cd,cd2,cd3,cd4
      integer i,j,idim
      real *8 threshsq

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd
          cd = sqrt(dinv2)
          cd2 = -cd*dinv2
          cd3 = -3*cd*dinv2*dinv2

          do idim=1,nd
          
            pot(idim,i) = pot(idim,i) + cd*charge(idim,j)
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod

            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j))
     2         + cd2*charge(idim,j)*zdiff(1) 
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
     2         + cd2*charge(idim,j)*zdiff(2) 
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
     2         + cd2*charge(idim,j)*zdiff(3)
          enddo
 1000     continue
        enddo
      enddo

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
      
c
cc     temporary variables
c
      integer i,j,idim
      real *8 zdiff(3),dd,d,dinv,dinv2,dotprod,cd,cd2,cd3,cd4,cd5
      real *8 htmp1,htmp2,htmp3,htmp4,htmp5,htmp6,dx,dy,dz
      real *8 threshsq

      threshsq = thresh**2

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          if(dd.lt.threshsq) goto 1000

          dinv2 = 1/dd
          cd = sqrt(dinv2)
          cd2 = -cd*dinv2
          cd3 = -3*cd*dinv2*dinv2
c
          cd4 = -cd2/dd
c
          cd5 = -cd3
          dx = zdiff(1)*cd
          dy = zdiff(2)*cd
          dz = zdiff(3)*cd
c
          htmp1 = cd4*(3*zdiff(1)*zdiff(1)-dd)
          htmp2 = cd4*(3*zdiff(2)*zdiff(2)-dd)
          htmp3 = cd4*(3*zdiff(3)*zdiff(3)-dd)
          htmp4 = cd4*3*zdiff(1)*zdiff(2)
          htmp5 = cd4*3*zdiff(1)*zdiff(3)
          htmp6 = cd4*3*zdiff(2)*zdiff(3)
          do idim=1,nd
          
            pot(idim,i) = pot(idim,i) + cd*charge(idim,j)
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod

            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j))
     2         + cd2*charge(idim,j)*zdiff(1) 
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
     2         + cd2*charge(idim,j)*zdiff(2) 
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
     2         + cd2*charge(idim,j)*zdiff(3)
            hess(idim,1,i) = hess(idim,1,i)+htmp1*charge(idim,j)+
     1      (dotprod*(5*dx*dx-1)-
     2      (dipvec(idim,1,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(1)))*cd5
            hess(idim,2,i) = hess(idim,2,i)+htmp2*charge(idim,j)+
     1      (dotprod*(5*dy*dy-1)-
     2      (dipvec(idim,2,j)*zdiff(2)+dipvec(idim,2,j)*zdiff(2)))*cd5
            hess(idim,3,i) = hess(idim,3,i)+htmp3*charge(idim,j)+
     1      (dotprod*(5*dz*dz-1)-
     2      (dipvec(idim,3,j)*zdiff(3)+dipvec(idim,3,j)*zdiff(3)))*cd5
            hess(idim,4,i) = hess(idim,4,i)+htmp4*charge(idim,j)+
     1      (dotprod*(5*dx*dy)-
     2      (dipvec(idim,2,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(2)))*cd5
            hess(idim,5,i) = hess(idim,5,i)+htmp5*charge(idim,j)+
     1      (dotprod*(5*dx*dz)-
     2      (dipvec(idim,3,j)*zdiff(1)+dipvec(idim,1,j)*zdiff(3)))*cd5
            hess(idim,6,i) = hess(idim,6,i)+htmp6*charge(idim,j)+
     1      (dotprod*(5*dy*dz)-
     2      (dipvec(idim,3,j)*zdiff(2)+dipvec(idim,2,j)*zdiff(3)))*cd5
          enddo
 1000     continue
        enddo
      enddo
      return
      end
