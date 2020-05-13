      implicit real *8 (a-h,o-z)
      integer ntests,nsuccess
      complex *16 zk,ima

      data ima/(0.0d0,1.0d0)/

      done = 1
      pi = atan(done)*4

      zk = 300.0d0
      call test_helm_zkbig(zk,i1)

      zk = 300.0d0*ima
      call test_helm_zkbig(zk,i2)

      zk = 12*pi*ima*4
      call test_helm_zkbig(zk,i3)


      zk = 12*pi*ima/32

      call test_helm_zkbig(zk,i4)

      zk = 12*pi*ima + 0.01d0
      call test_helm_zkbig(zk,i5)

      zk = (12*pi*ima + 12*pi)/10
      call test_helm_zkbig(zk,i6)

      nsuccess = i1+i2+i3+i4+i5+i6
      ntest = 6
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ', nsuccess,
     1   ' out of ', ntest, ' in hfmm3d zkbig testing suite'

      stop
      end



      subroutine test_helm_zkbig(zk,isuccess)
      
      implicit none
      integer ns
      double precision, allocatable :: source(:,:)
      double complex, allocatable :: charge(:),dipvec(:,:)
      double complex, allocatable :: pot(:),potex(:)
      double complex, allocatable :: grad(:,:),gradex(:,:)

      double precision eps
      double complex eye,zk
      integer i,j,k,ntest
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


      eps = 0.51d-3


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
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2

        charge(i) = hkrand(0) + eye*hkrand(0)
        dipvec(1,i) = hkrand(0) + eye*hkrand(0)
        dipvec(2,i) = hkrand(0) + eye*hkrand(0)
        dipvec(3,i) = hkrand(0) + eye*hkrand(0)

        pot(i) = 0
        grad(1,i) = 0
        grad(2,i) = 0
        grad(3,i) = 0 
      enddo


      call hfmm3d_s_c_g(eps,zk,ns,source,charge,
     1      pot,grad)
      ntest = 10
      allocate(potex(ntest),gradex(3,ntest))

      thresh = 2.0d0**(-51)

      do i=1,ntest
        potex(i) = 0
        gradex(1,i) = 0
        gradex(2,i) = 0
        gradex(3,i) = 0
      enddo

      call h3ddirectcg(1,zk,source,charge,ns,source,
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

      erra = sqrt(erra/ra)

      call prin2('zk=*',zk,1)
      call prin2('l2 rel err=*',erra,1)

      isuccess = 0

      if(erra.lt.eps) isuccess = 1
      

      return
      end
c----------------------------------------------------------
c
cc
c
c
