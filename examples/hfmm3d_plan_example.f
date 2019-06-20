      implicit none
      integer ns,nt
      double precision, allocatable :: source(:,:)
      double precision, allocatable :: targ(:,:)

      double precision eps
      double complex eye,zk
      integer i,j,k
      double precision hkrand

      integer *8 ltree,nwork
      integer nboxes,nlevels

      integer iptype,ifpgh,ifpghtarg
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      print *, "This code is a test driver for fmm3d plan interface"
      write(*,*)
      write(*,*)


      zk = 2.2d0

      ns = 20 000 000
      

      allocate(source(3,ns))


      eps = 0.5d-6

      write(*,*) "=========================================="


c
cc      generate sources uniformly in the unit cube 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2
      enddo

      ifpgh = 1
      ifpghtarg = 0

      allocate(targ(3,1))

      nt = 0

      nwork = 0
      ltree = 0
      nboxes = 0
      nlevels = 0

      iptype = 3


      call hfmm3d_plan(eps,zk,ns,source,ifpgh,nt,targ,
     1       ifpghtarg,iptype,ltree,nboxes,nlevels,nwork)

      stop
      end
c----------------------------------------------------------
c
cc
c
c
