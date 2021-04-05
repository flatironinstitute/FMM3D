      implicit real *8 (a-h,o-z)
      double precision, allocatable :: source(:,:),targ(:,:)
      double complex zk,ier
      call prini(6,13)

      ns = 100000
      nt = 100000
      
      allocate(source(3,ns),targ(3,nt))

      eps = 0.51d-10
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2
      enddo

      do i=1,nt
        targ(1,i) = hkrand(0)**2
        targ(2,i) = hkrand(0)**2
        targ(3,i) = hkrand(0)**2
      enddo

      ifcharge = 1
      ifdipole = 1

      ifpgh = 1
      ifpghtarg = 1
      
      nd = 3
      zk = 10.0d0 + ima*0.3d0
      iper = 0
      call hfmm3d_memest(nd,eps,zk,ns,source,ifcharge,ifdipole,
     1 iper,ifpgh,nt,targ,ifpghtarg,rmem)
      
      call prin2('rmem=*',rmem,1)




      stop
      end
