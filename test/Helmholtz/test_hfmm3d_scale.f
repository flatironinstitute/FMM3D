      implicit real *8 (a-h,o-z)
      complex *16 z1,z2,zk,ima
      data ima/(0.0d0,1.0d0)/

      done = 1
      pi = atan(done)*4

      ntest = 6

      z1 = 8 + ima*0.2d0

      rsc = 1.0d10
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i1)

      rsc = 1.0d-9
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i2)

      z1 = 5.0d0 + ima*5.0d0
      rsc = 1.0d10
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i3)

      rsc = 1.0d-9
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i4)

      z1 = 1.0d0 + ima*5.0d0
      rsc = 1.0d10
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i5)

      rsc = 1.0d-9
      zk = z1*2*pi/rsc
      call test_hfmm3d_scale(rsc,zk,i6)
      
      nsuccess = i1+i2+i3+i4+i5+i6
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ', nsuccess,
     1   ' out of ', ntest, ' in hfmm3d scale testing suite'

      stop

      stop
      end


      subroutine test_hfmm3d_scale(rsc,zk,isuccess)
      
      implicit none
      integer ns
      double precision, allocatable :: source(:,:)
      double complex, allocatable :: charge(:),dipvec(:,:)
      double complex, allocatable :: pot(:),potex(:)
      double complex, allocatable :: grad(:,:),gradex(:,:)

      double precision eps,rsc
      double complex eye,zk
      integer i,j,k,ntest,ier
      double precision hkrand,thresh,erra,ra
      integer isuccess
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      write(*,*)
      write(*,*)



      ns = 10000
      

      allocate(source(3,ns))
      allocate(charge(ns))
      allocate(dipvec(3,ns))
      allocate(pot(ns))
      allocate(grad(3,ns))


      eps = 0.51d-6


c
c   
c       example demonstrating use of 
c        source to source, charges, pot + gradient
c



c
cc      generate sources with distribution unif^2 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)**2*rsc
        source(2,i) = hkrand(0)**2*rsc
        source(3,i) = hkrand(0)**2*rsc

        charge(i) = hkrand(0) + eye*hkrand(0)
        dipvec(1,i) = hkrand(0) + eye*hkrand(0)
        dipvec(2,i) = hkrand(0) + eye*hkrand(0)
        dipvec(3,i) = hkrand(0) + eye*hkrand(0)

        pot(i) = 0
        grad(1,i) = 0
        grad(2,i) = 0
        grad(3,i) = 0 
      enddo


      call hfmm3d_s_cd_g(eps,zk,ns,source,charge,dipvec,
     1      pot,grad,ier)
      ntest = 10
      allocate(potex(ntest),gradex(3,ntest))

      thresh = 2.0d0**(-51)*rsc

      do i=1,ntest
        potex(i) = 0
        gradex(1,i) = 0
        gradex(2,i) = 0
        gradex(3,i) = 0
      enddo

      call h3ddirectcdg(1,zk,source,charge,dipvec,ns,source,
     1  ntest,potex,gradex,thresh)

      erra = 0
      ra = 0
      do i=1,ntest
        erra = erra + abs(potex(i)-pot(i))**2
        erra = erra + abs(gradex(1,i)-grad(1,i))**2
        erra = erra + abs(gradex(2,i)-grad(2,i))**2
        erra = erra + abs(gradex(3,i)-grad(3,i))**2

        ra = ra + abs(potex(i))**2
        ra = ra + abs(gradex(1,i))**2 + abs(gradex(2,i))**2 +
     1   abs(gradex(3,i))**2
      enddo

      call prin2('pot=*',pot,2*ntest)
      call prin2('potex=*',potex,2*ntest)

      erra = sqrt(erra/ra)

      call prin2('zk=*',zk,1)
      call prin2('rsc=*',rsc,1)
      call prin2('l2 rel err=*',erra,1)

      isuccess = 0

      if(erra.lt.5*eps) isuccess = 1
      

      return
      end
c----------------------------------------------------------
c
