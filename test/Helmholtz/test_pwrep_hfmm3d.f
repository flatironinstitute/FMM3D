      implicit none
      integer ns,nt
      double precision, allocatable :: source(:,:)
      double complex, allocatable :: charge(:)
      double complex, allocatable :: pot(:),potex(:)

      double precision eps(5)
      double complex eye,zk,rkmin,rkmax
      double complex, allocatable :: zks(:)
      integer i,j,k,ntest
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      integer len1,ntests,isum
      integer idomain,ndom
      integer, allocatable :: ipass(:,:)
      double precision err,hkrand,erra,ra,tt
      integer nprec,iprec
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)


      ns = 2000
      nt = 1999
      
      ntest = 10

      tt = 1.0d-16

      allocate(source(3,ns))
      allocate(charge(ns))
      allocate(pot(ns),potex(ns))



      write(*,*) "=========================================="
      write(*,*) "Testing suite for pw representations for hfmm3d"

      ndom = 437


      nprec = 5
      eps(1) = 0.51d-2
      eps(2) = 0.51d-3
      eps(3) = 0.51d-6
      eps(4) = 0.51d-9
      eps(5) = 0.51d-12


      allocate(ipass(ndom,nprec))
      allocate(zks(ndom))

      do idomain=1,ndom

        call getrkinfo2(idomain,ndom,rkmin,rkmax)
        zks(idomain) = (rkmin+rkmax)*2
      enddo

c
cc      generate sources uniformly in the unit cube 
c
c
      do i=1,ns
        source(1,i) = hkrand(0)
        source(2,i) = hkrand(0)
        source(3,i) = hkrand(0)

        charge(i) = hkrand(0) + eye*hkrand(0)

      enddo

      ndom = 437
      nprec = 5
      open(unit=33,file='pw_res.log')
      open(unit=34,file='pw_fail.log')

      do iprec=1,nprec
        do idomain=1,ndom

          do i=1,ns
            pot(i) = 0
            potex(i) = 0
          enddo



          call hfmm3d_s_c_p(eps(iprec),zks(idomain),ns,source,charge,
     1          pot)
       
       
          call h3ddirectcp(1,zks(idomain),source,charge,ns,source,
     1        ntest,potex,tt)

          ra = 0
          erra = 0
          do i=1,ntest
            ra = ra + abs(potex(i))**2
            erra = erra + abs(potex(i)-pot(i))**2
          enddo

          erra = sqrt(erra/ra)
          write(33,'(a,i3,a,i1,a,e11.5)')
     1        'idom=',idomain,',iprec=',iprec-1,',err=',erra      
          if(erra.le.eps(iprec)) then
            ipass(idomain,iprec) = 1
            write(13,'(a,i3,a,i1,a)') 'success: idom=',idomain,
     1                 ',iprec=',iprec-1
          else
            ipass(idomain,iprec) = 0
            write(13,'(a,i3,a,i1,a)') 'fail: idom=',idomain,
     1                 ',iprec=',iprec-1
          write(34,'(a,i3,a,i1,a,e11.5)')
     1        'idom=',idomain,',iprec=',iprec-1,',err=',erra      
             
          endif
        enddo
      enddo
       

       stop
       end
