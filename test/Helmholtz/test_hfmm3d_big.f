      implicit none
      integer ns
      double precision, allocatable :: source(:,:)
      double complex, allocatable :: charge(:)
      double complex, allocatable :: pot(:)
      double complex, allocatable :: potex(:)

      double precision eps
      double complex eye,zk
      integer i,j,k
      integer ntest,nd
      double precision thresh,ra,erra
      double precision hkrand
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
cc      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      print *, "This code is an example fortran driver*"
      write(*,*)
      write(*,*)


      zk = 2.2d0

      ns = 300 000 000   
      print *, "nsources =",ns
      

      allocate(source(3,ns))
      allocate(charge(ns))
      allocate(pot(ns))


      eps = 0.5d-3

      write(*,*) "=========================================="

c
c   
c       example demonstrating use of 
c        source to source, charges, pot + gradient
c



c
cc      generate sources uniformly in the unit cube 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2

        charge(i) = hkrand(0) + eye*hkrand(0)

        pot(i) = 0
      enddo


      call hfmm3d_s_c_p(eps,zk,ns,source,charge,pot)


      ntest = 10
      allocate(potex(ntest))
      do i=1,ntest 
        potex(i) = 0
      enddo

      thresh = 1.0d-16
      nd = 1
      
      call h3ddirectcp(nd,zk,source,charge,ns,source,ntest,potex,
     1   thresh)
      
      ra = 0
      erra = 0
      do i=1,ntest
        ra = ra + abs(potex(i))**2
        erra = erra + abs(pot(i)-potex(i))**2
      enddo

      erra = sqrt(erra/ra)
      print *, "Relative error in pot=",erra
      


      stop
      end
c----------------------------------------------------------
c
cc
c
c
