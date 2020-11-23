      implicit real *8 (a-h,o-z)

      done = 1
      pi = atan(done)*4

      ntest = 2


      rsc = 1.0d10
      call test_lfmm3d_scale(rsc,i1)

      rsc = 1.0d-9
      call test_lfmm3d_scale(rsc,i2)

      nsuccess = i1+i2
      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ', nsuccess,
     1   ' out of ', ntest, ' in lfmm3d scale testing suite'


      stop
      end


      subroutine test_lfmm3d_scale(rsc,isuccess)
      
      implicit none
      integer ns
      double precision, allocatable :: source(:,:)
      double precision, allocatable :: charge(:),dipvec(:,:)
      double precision, allocatable :: pot(:),potex(:)
      double precision, allocatable :: grad(:,:),gradex(:,:)

      double precision eps,rsc
      integer i,j,k,ntest,ier
      double precision hkrand,thresh,erra,ra
      integer isuccess
      

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

        charge(i) = hkrand(0) 
        dipvec(1,i) = hkrand(0) 
        dipvec(2,i) = hkrand(0) 
        dipvec(3,i) = hkrand(0) 

        pot(i) = 0
        grad(1,i) = 0
        grad(2,i) = 0
        grad(3,i) = 0 
      enddo


      call lfmm3d_s_cd_g(eps,ns,source,charge,dipvec,
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

      call l3ddirectcdg(1,source,charge,dipvec,ns,source,
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

      call prin2('pot=*',pot,ntest)
      call prin2('potex=*',potex,ntest)

      erra = sqrt(erra/ra)

      call prin2('rsc=*',rsc,1)
      call prin2('l2 rel err=*',erra,1)

      isuccess = 0

      if(erra.lt.5*eps) isuccess = 1
      

      return
      end
c----------------------------------------------------------
c
