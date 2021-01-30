      implicit real *8 (a-h,o-z)
      integer(8) ns,nt,nd,nttest
      double complex zk
      double precision, allocatable :: sources(:,:),targs(:,:)
      double complex, allocatable :: a_vect(:,:,:),lambda(:,:)
      double complex, allocatable :: b_vect(:,:,:)
      double complex, allocatable :: E(:,:,:),divE(:,:),curlE(:,:,:)
      double complex, allocatable :: Eex(:,:,:),divEex(:,:)
      double complex, allocatable :: curlEex(:,:,:)
      double complex ima
      double precision thresh
      double precision t1,t2,omp_get_wtime
      integer(8) ipass(10)
      data ima/(0.0d0,1.0d0)/
      integer(8) ifE,ifdivE,ifcurlE,ifa_vect,ifb_vect,iflambda,ier

      ns = 4001
      nt = 3998

      nd = 2
      zk = 1.1d0 + ima*0.1d0

      allocate(sources(3,ns))
      allocate(a_vect(nd,3,ns),b_vect(nd,3,ns),lambda(nd,ns))

      allocate(targs(3,nt),E(nd,3,nt),divE(nd,nt),curlE(nd,3,nt))

      nttest = 10
      allocate(Eex(nd,3,nttest),divEex(nd,nttest),curlEex(nd,3,nttest))


c
cc  generate sourcrs in the unit cube
c

c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)
      do i=1,ns 
        sources(1,i) = hkrand(0)**2
        sources(2,i) = hkrand(0)**2
        sources(3,i) = hkrand(0)**2
        do l=1,3
          do j=1,nd
            a_vect(j,l,i) = hkrand(0) + ima*hkrand(0)
            b_vect(j,l,i) = hkrand(0) + ima*hkrand(0)
          enddo
        enddo

        do j=1,nd
          lambda(j,i) = hkrand(0) + ima*hkrand(0)
        enddo
      enddo
C$OMP END PARALLEL DO      

c
c   generate targets in the unit cube
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nt
        targs(1,i) = hkrand(0)**2 
        targs(2,i) = hkrand(0)**2
        targs(3,i) = hkrand(0)**2
      enddo
C$OMP END PARALLEL DO      

      eps = 1.0d-3
      write(*,*) "=========================================="
      write(*,*) "Testing suite for hfmm3d"
      write(*,'(a,e11.4)') "Requested precision = ",eps

      open(unit=33,file='print_testres.txt',access='append')

      ntests = 6
      do i=1,ntests
        ipass(i) = 0
      enddo

c
cc  test b_vect 
cc   with fields
c
      write(6,*) 'interaction: b'
      write(6,*) 'output: fields'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 0
      ifcurlE = 0
      ifa_vect = 0
      ifb_vect = 1
      iflambda = 0

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)

      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1
      
      thresh = 1.0d-16

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
          enddo
        enddo
      enddo


      
      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      ra = 0
      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('rel error in field=*',erra,1)
      if(erra.lt.eps) ipass(1) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="
c
cc  test b_vect 
cc   with fields, divergence and curl
c
      write(6,*) 'interaction: b'
      write(6,*) 'output: fields, divergence, and curl'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 1
      ifcurlE = 1
      ifa_vect = 0
      ifb_vect = 1
      iflambda = 0

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)

      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1
      thresh = 1.0d-16

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
            curlEex(j,l,i) = 0
          enddo
        enddo

        do j=1,nd
          divEex(j,i) = 0
        enddo
      enddo

      
      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      errc = 0
      errd = 0
      ra = 0
      rc = 0
      rd = 0

      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            rc = rc + abs(curlEex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
            errc = errc + abs(curlEex(j,l,i)-curlE(j,l,i))**2
          enddo
        enddo

        do j=1,nd
          rd = rd + abs(divEex(j,i))**2
          errd = errd + abs(divEex(j,i)-divE(j,i))**2
        enddo
      enddo

      erra = sqrt(erra/ra)
      errd = sqrt(errd/rd)
      errc = sqrt(errc/rc)

      call prin2('rel error in field=*',erra,1)
      call prin2('rel error in divergence=*',errd,1)
      call prin2('rel error in curl=*',errc,1)
      if(erra.lt.eps.and.errd.lt.eps.and.errc.lt.eps) ipass(2) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="

c
cc  test a_vect, lambda 
cc   with fields
c
      write(6,*) 'interaction: a,lambda'
      write(6,*) 'output: fields'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 0
      ifcurlE = 0
      ifa_vect = 1
      ifb_vect = 0
      iflambda = 1

      call cpu_time(t1)
C$     t1 = omp_get_wtime()      
      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1

      thresh = 1.0d-16
      

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
          enddo
        enddo
      enddo

      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      ra = 0
      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('rel error in field=*',erra,1)
      if(erra.lt.eps) ipass(3) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="
c
cc  test a_vect, lambda 
cc   with fields, divergence and curl
c

      write(6,*) 'interaction: a,lambda'
      write(6,*) 'output: fields, divergence, and curl'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 1
      ifcurlE = 1
      ifa_vect = 1
      ifb_vect = 0
      iflambda = 1

      call cpu_time(t1)
C$     t1 = omp_get_wtime()      
      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)

      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1
      thresh = 1.0d-16
      

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
            curlEex(j,l,i) = 0
          enddo
        enddo

        do j=1,nd
          divEex(j,i) = 0
        enddo
      enddo

      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      errc = 0
      errd = 0
      ra = 0
      rc = 0
      rd = 0

      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            rc = rc + abs(curlEex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
            errc = errc + abs(curlEex(j,l,i)-curlE(j,l,i))**2
          enddo
        enddo

        do j=1,nd
          rd = rd + abs(divEex(j,i))**2
          errd = errd + abs(divEex(j,i)-divE(j,i))**2
        enddo
      enddo

      erra = sqrt(erra/ra)
      errd = sqrt(errd/rd)
      errc = sqrt(errc/rc)

      call prin2('rel error in field=*',erra,1)
      call prin2('rel error in divergence=*',errd,1)
      call prin2('rel error in curl=*',errc,1)
      if(erra.lt.eps.and.errd.lt.eps.and.errc.lt.eps) ipass(4) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="

c
cc  test b_vect,a_vect, lambda 
cc   with fields
c
      write(6,*) 'interaction: b,a,lambda'
      write(6,*) 'output: fields'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 0
      ifcurlE = 0
      ifa_vect = 1
      ifb_vect = 1
      iflambda = 1

      call cpu_time(t1)
C$     t1 = omp_get_wtime()      
      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)

      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1
      thresh = 1.0d-16

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
          enddo
        enddo
      enddo
      
      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      ra = 0
      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('rel error in field=*',erra,1)
      if(erra.lt.eps) ipass(5) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="
c
cc  test b_vect,a_vect, lambda 
cc   with fields, divergence and curl
c
      write(6,*) 'interaction: b,a,lambda'
      write(6,*) 'output: fields, divergence, and curl'
      write(6,*) 
      write(6,*)

      ifE = 1
      ifdivE = 1
      ifcurlE = 1
      ifa_vect = 1
      ifb_vect = 1
      iflambda = 1

      call cpu_time(t1)
C$     t1 = omp_get_wtime()      
      call emfmm3d(nd,eps,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1  b_vect,iflambda,lambda,nt,targs,ifE,E,ifcurlE,curlE,
     2  ifdivE,divE,ier)

      call cpu_time(t2)
C$      t2 = omp_get_wtime() 
      
      call prin2('time=*',t2-t1,1)
      write(14,*) "time=",t2-t1
      thresh = 1.0d-16
      

      do i=1,nttest
        do l=1,3
          do j=1,nd
            Eex(j,l,i) = 0
            curlEex(j,l,i) = 0
          enddo
        enddo

        do j=1,nd
          divEex(j,i) = 0
        enddo
      enddo

      call em3ddirect(nd,zk,ns,sources,ifa_vect,a_vect,ifb_vect,
     1 b_vect,iflambda,lambda,nttest,targs,ifE,Eex,ifcurlE,curlEex,
     2 ifdivE,divEex,thresh)
      
      erra = 0
      errc = 0
      errd = 0
      ra = 0
      rc = 0
      rd = 0

      do i=1,nttest
        do l=1,3
          do j=1,nd
            ra = ra + abs(Eex(j,l,i))**2
            rc = rc + abs(curlEex(j,l,i))**2
            erra = erra + abs(Eex(j,l,i)-E(j,l,i))**2
            errc = errc + abs(curlEex(j,l,i)-curlE(j,l,i))**2
          enddo
        enddo

        do j=1,nd
          rd = rd + abs(divEex(j,i))**2
          errd = errd + abs(divEex(j,i)-divE(j,i))**2
        enddo
      enddo

      erra = sqrt(erra/ra)
      errd = sqrt(errd/rd)
      errc = sqrt(errc/rc)

      call prin2('rel error in field=*',erra,1)
      call prin2('rel error in divergence=*',errd,1)
      call prin2('rel error in curl=*',errc,1)
      if(erra.lt.eps.and.errd.lt.eps.and.errc.lt.eps) ipass(6) = 1
      write(6,*)
      write(6,*)
      write(6,*) "====================="

      isum = 0
      do i=1,ntests
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in emfmm3d testing suite'
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in emfmm3d testing suite'
      close(33)
      

      stop

      

      stop
      end
