c      This file contains the direct evaluation kernels for Helmholtz FMM
c
c      h3ddirectcp:  direct calculation of potentials for a collection 
c                         of charge sources at a collection of targets
c
c      h3ddirectcg:  direct calculation of potentials and gradients 
c                         for a collection of charge sources at a 
c                         collection of targets
c
c      h3ddirectch:  direct calculation of potentials, gradients and 
c                         Hessians for a collection of charge sources 
c                         at a collection of targets
c
c      h3ddirectdp:  direct calculation of potentials for a collection 
c                         of dipole sources at a collection of targets
c
c      h3ddirectdg:  direct calculation of potentials and gradients 
c                         for a collection of dipole sources at a 
c                         collection of targets
c
c      h3ddirectdh:  direct calculation of potentials, gradients and
c                         Hessians for a collection of dipole sources 
c                         at a collection of targets
c
c      h3ddirectcdp:  direct calculation of potentials for a 
c                         collection of charge and dipole sources at 
c                         a collection of targets
c
c      h3ddirectcdg:  direct calculation of potentials and gradients 
c                         for a collection of charge and dipole sources 
c                         at a collection of targets
c
c      h3ddirectcdh:  direct calculation of potentials, gradients and
c                         Hessians for a collection of charge and dipole 
c                         sources at a collection of targets
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d
      complex *16 zkeye,eye,ztmp
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          ztmp = exp(zkeye*d)/d
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
c
C***********************************************************************
      subroutine h3ddirectcg(nd,zk,sources,charge,ns,ztarg,nt,
     1            pot,grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x)=pot(x)+sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}| 
c                        j
c                 
c     grad(x)=grad(x)+Gradient(sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}|) 
c                               j
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d
      complex *16 zkeye,eye,cd,cd1,ztmp
      complex *16 ztmp1,ztmp2,ztmp3
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000
          cd = exp(zkeye*d)/d
          cd1 = (zkeye*d-1)*cd/dd
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d
      complex *16 zkeye,zkeyed,eye,cd,cd1,cd2,ztmp
      complex *16 ztmp1,ztmp2,ztmp3
      complex *16 htmp1,htmp2,htmp3,htmp4,htmp5,htmp6
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000
          zkeyed = zkeye*d
          cd = exp(zkeyed)/d
          cd1 = (zkeyed-1)*cd/dd
          cd2 = (zkeyed*zkeyed - 3*(zkeyed-1))*cd/(dd*dd)
          ztmp1 = cd1*zdiff(1)
          ztmp2 = cd1*zdiff(2)
          ztmp3 = cd1*zdiff(3)
          htmp1 = cd2*zdiff(1)*zdiff(1)+cd1
          htmp2 = cd2*zdiff(2)*zdiff(2)+cd1
          htmp3 = cd2*zdiff(3)*zdiff(3)+cd1
          htmp4 = cd2*zdiff(1)*zdiff(2)
          htmp5 = cd2*zdiff(1)*zdiff(3)
          htmp6 = cd2*zdiff(2)*zdiff(3)
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
      subroutine h3ddirectdp(nd,zk,sources,
     1            dipvec,ns,ztarg,nt,pot,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential due to a collection
c     of sources and adds to existing
c     quantities.
c
c     pot(x)=pot(x)+sum   \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j} 
c                    j
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv
      complex *16 zkeye,eye,cd,cd1,dotprod
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          cd = exp(zkeye*d)*dinv
          cd1 = (1-zkeye*d)*cd/dd

          do idim=1,nd
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
      subroutine h3ddirectdg(nd,zk,sources,dipvec,ns,ztarg,nt,pot,
     1   grad,thresh)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient due to a 
c     collection of sources and adds to existing quantities.
c
c     pot(x)=pot(x)+sum  d_{j} \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                    j
c
c   
c     grad(x)=grad(x)+Gradient( sum  
c                                j
c
c                         \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c                         )
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2
      complex *16 zkeye,eye,cd,cd2,cd3,cd4,dotprod
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          dinv2 = dinv**2
          cd = exp(zkeye*d)*dinv
          cd2 = (zkeye*d-1)*cd*dinv2
          cd3 = cd*dinv2*(-zkeye*zkeye-3*dinv2+3*zkeye*dinv)

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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2,dx,dy,dz,dinv5
      complex *16 zkeye,eye,cd,cd2,cd3,cd4,dotprod
      complex *16 cross12,cross13,cross23,zt0,zf1,zf0
      complex *16 cdcommon,cdcommon2,cdcommon3,cdcommon4
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)
          dx = zdiff(1)
          dy = zdiff(2)
          dz = zdiff(3)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          dinv2 = dinv**2
          dinv5 = dinv2*dinv2*dinv
          cd = exp(zkeye*d)*dinv
          zf1 = (zkeye*d-1)
          zf0 = zkeye*d
          cd2 = (zkeye*d-1)*cd*dinv2
          cd3 = cd*dinv2*(-zkeye*zkeye-3*dinv2+3*zkeye*dinv)
          cdcommon = cd*dinv2*dinv2*(3*zf1 - zf0*zf0)
          cdcommon2 = cd*
     &    ((zk**2*(zf1+2)*dinv2*dinv2) -1.5d1*zf1*dinv2*dinv2*dinv2
     &    +1.75d0*eye*zk*zf0*4*dinv5)
          cdcommon3 = cd*(6.0d0*zf1-zf0*zf0*2)*dinv2*dinv2
          cdcommon4 = cd*dinv2*dinv2*(
     &     +zk**2*(zf1+2.0d0)-zf1*15/dd
     &     +zf0*zf0*7.0D0/dd)
c
          do idim=1,nd
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod
            cross12 = dipvec(idim,1,j)*dy + dipvec(idim,2,j)*dx
            cross13 = dipvec(idim,1,j)*dz + dipvec(idim,3,j)*dx
            cross23 = dipvec(idim,3,j)*dy + dipvec(idim,2,j)*dz
c
            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j)) 
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
            hess(idim,1,i) = hess(idim,1,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dx*dx+
     1         cdcommon3*dipvec(idim,1,j)*dx
            hess(idim,2,i) = hess(idim,2,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dy*dy+
     1         cdcommon3*dipvec(idim,2,j)*dy
            hess(idim,3,i) = hess(idim,3,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dz*dz+
     1         cdcommon3*dipvec(idim,3,j)*dz
            hess(idim,4,i) = hess(idim,4,i) + 
     1         cdcommon*cross12 + cdcommon4*dotprod*dx*dy
            hess(idim,5,i) = hess(idim,5,i) + 
     1         cdcommon*cross13 + cdcommon4*dotprod*dx*dz
            hess(idim,6,i) = hess(idim,6,i) + 
     1         cdcommon*cross23 + cdcommon4*dotprod*dz*dy
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
c                       \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 charge(nd,ns),pot(nd,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv
      complex *16 zkeye,eye,cd,cd1,dotprod
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          cd = exp(zkeye*d)*dinv
          cd1 = (1-zkeye*d)*cd/dd

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
c                         \nabla e^{ik |x-x_{j}|/|x-x_{j}| \cdot v_{j}
c   
c     grad(x)=grad(x)+Gradient( sum  q_{j} e^{i k |x-x_{j}|}/|x-x_{j}|+  
c                                j
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
c
c-----------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,zk,sources,charge,dipvec,ns,ztarg,nt,thresh
cf2py intent(out) pot,grad
c
cc      calling sequence variables
c  
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 dipvec(nd,3,ns)
      complex *16 charge(nd,ns),pot(nd,nt),grad(nd,3,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2
      complex *16 zkeye,eye,cd,cd2,cd3,cd4,dotprod
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          dinv2 = dinv**2
          cd = exp(zkeye*d)*dinv
          cd2 = (zkeye*d-1)*cd*dinv2
          cd3 = cd*dinv2*(-zkeye*zkeye-3*dinv2+3*zkeye*dinv)

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
      integer(8) ns,nt,nd
      complex *16 zk
      real *8 sources(3,ns),ztarg(3,nt)
      complex *16 charge(nd,ns),dipvec(nd,3,ns)
      complex *16 pot(nd,nt),grad(nd,3,nt)
      complex *16 hess(nd,6,nt)
      real *8 thresh
      
c
cc     temporary variables
c
      real *8 zdiff(3),dd,d,dinv,dinv2,dx,dy,dz,dinv5
      complex *16 zkeye,eye,cd,cd1,cd2,cd3,cd4,dotprod
      complex *16 cross12,cross13,cross23,zt0,zf1,zf0
      complex *16 ztmp1,ztmp2,ztmp3
      complex *16 htmp1,htmp2,htmp3,htmp4,htmp5,htmp6
      complex *16 cdcommon,cdcommon2,cdcommon3,cdcommon4
      integer(8) i,j,idim
      data eye/(0.0d0,1.0d0)/

      zkeye = zk*eye

      do i=1,nt
        do j=1,ns
          zdiff(1) = ztarg(1,i)-sources(1,j)
          zdiff(2) = ztarg(2,i)-sources(2,j)
          zdiff(3) = ztarg(3,i)-sources(3,j)
          dx = zdiff(1)
          dy = zdiff(2)
          dz = zdiff(3)

          dd = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
          d = sqrt(dd)
          if(d.lt.thresh) goto 1000

          dinv = 1/d
          dinv2 = dinv**2
          dinv5 = dinv2*dinv2*dinv
          cd = exp(zkeye*d)*dinv
          cd1 = (zkeye*d-1)*cd/dd
          zf1 = (zkeye*d-1)
          zf0 = zkeye*d
          cd2 = (zkeye*d-1)*cd*dinv2
          cd3 = cd*dinv2*(-zkeye*zkeye-3*dinv2+3*zkeye*dinv)
c
          ztmp1 = cd1*zdiff(1)
          ztmp2 = cd1*zdiff(2)
          ztmp3 = cd1*zdiff(3)
          htmp1 = -cd3*zdiff(1)*zdiff(1)+cd1
          htmp2 = -cd3*zdiff(2)*zdiff(2)+cd1
          htmp3 = -cd3*zdiff(3)*zdiff(3)+cd1
          htmp4 = -cd3*zdiff(1)*zdiff(2)
          htmp5 = -cd3*zdiff(1)*zdiff(3)
          htmp6 = -cd3*zdiff(2)*zdiff(3)
c
          cdcommon = cd*dinv2*dinv2*(3*zf1 - zf0*zf0)
          cdcommon2 = cd*
     &    ((zk**2*(zf1+2)*dinv2*dinv2) -1.5d1*zf1*dinv2*dinv2*dinv2
     &    +1.75d0*eye*zk*zf0*4*dinv5)
          cdcommon3 = cd*(6.0d0*zf1-zf0*zf0*2)*dinv2*dinv2
          cdcommon4 = cd*dinv2*dinv2*(
     &     +zk**2*(zf1+2.0d0)-zf1*15/dd
     &     +zf0*zf0*7.0D0/dd)
c
          do idim=1,nd
c
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
c
c         dipole contributions
c
            dotprod = zdiff(1)*dipvec(idim,1,j)+
     1               zdiff(2)*dipvec(idim,2,j)+
     1               zdiff(3)*dipvec(idim,3,j)
            cd4 = cd3*dotprod
            cross12 = dipvec(idim,1,j)*dy + dipvec(idim,2,j)*dx
            cross13 = dipvec(idim,1,j)*dz + dipvec(idim,3,j)*dx
            cross23 = dipvec(idim,3,j)*dy + dipvec(idim,2,j)*dz
c
            pot(idim,i) = pot(idim,i) - cd2*dotprod
            grad(idim,1,i) = grad(idim,1,i) + (cd4*zdiff(1) - 
     1         cd2*dipvec(idim,1,j)) 
            grad(idim,2,i) = grad(idim,2,i) + (cd4*zdiff(2) - 
     1         cd2*dipvec(idim,2,j))
            grad(idim,3,i) = grad(idim,3,i) + (cd4*zdiff(3) - 
     1         cd2*dipvec(idim,3,j))
            hess(idim,1,i) = hess(idim,1,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dx*dx+
     1         cdcommon3*dipvec(idim,1,j)*dx
            hess(idim,2,i) = hess(idim,2,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dy*dy+
     1         cdcommon3*dipvec(idim,2,j)*dy
            hess(idim,3,i) = hess(idim,3,i) + 
     1         cdcommon*dotprod +cdcommon2*dotprod*dz*dz+
     1         cdcommon3*dipvec(idim,3,j)*dz
            hess(idim,4,i) = hess(idim,4,i) + 
     1         cdcommon*cross12 + cdcommon4*dotprod*dx*dy
            hess(idim,5,i) = hess(idim,5,i) + 
     1         cdcommon*cross13 + cdcommon4*dotprod*dx*dz
            hess(idim,6,i) = hess(idim,6,i) + 
     1         cdcommon*cross23 + cdcommon4*dotprod*dz*dy
          enddo
 1000     continue
        enddo
      enddo
      return
      end
c
c
c

