c
c      TESTING SUITE FOR STOKES FMM:  June 18th, 2020
c
c

c$    use omp_lib
      implicit real *8 (a-h,o-z)
      
      real *8, allocatable :: pot(:,:), pre(:), grad(:,:,:)
      real *8, allocatable :: source(:,:), targ(:,:)
      real *8, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:)      
      real *8, allocatable :: pottarg(:,:), pretarg(:), gradtarg(:,:,:)
      real *8, allocatable :: pottarg2(:,:),pretarg2(:),gradtarg2(:,:,:)

      real *8, allocatable :: stoklet(:,:), strslet(:,:)
      real *8, allocatable :: strsvec(:,:)
      
      integer ipass(10)
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
      write(*,*) "=========================================="
      write(*,*) "Testing suite for stfmm3d"

      open(unit=33,file='print_testres.txt',access='append')

      ntests = 6
      do i=1,ntests
        ipass(i) = 0 
      enddo


      ns = 4000
      nt = 4200

      allocate(source(3,ns),targ(3,nt))
      allocate(stoklet(3,ns),strslet(3,ns),strsvec(3,ns))
      
      do i = 1,ns
         source(1,i) = cos(i*done)
         source(2,i) = cos(10*i*done)
         source(3,i) = cos(100*i*done)
         stoklet(1,i) = sin(2*i*done)
         stoklet(2,i) = sin(22*i*done)
         stoklet(3,i) = sin(222*i*done)
         strslet(1,i) = sin(4*i*done)
         strslet(2,i) = sin(44*i*done)
         strslet(3,i) = sin(444*i*done)
         strsvec(1,i) = sin(5*i*done)
         strsvec(2,i) = sin(55*i*done)
         strsvec(3,i) = sin(555*i*done)
      enddo

      do i = 1,nt
         targ(1,i) = cos(7*i*done+1.2)
         targ(2,i) = cos(77*i*done+1.2)
         targ(3,i) = cos(777*i*done+1.2)
      enddo

      allocate(pot(3,ns),pre(ns),grad(3,3,ns))
      allocate(pottarg(3,nt),pretarg(nt),gradtarg(3,3,nt))      

      nd = 1
      eps = 1d-9
      ifstoklet = 1
      ifstrslet = 1
      ifppreg = 3
      ifppregtarg = 3

      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call stfmm3d(nd,eps,ns,source,ifstoklet,stoklet,
     1     ifstrslet,strslet,strsvec,ifppreg,pot,pre,grad,
     2     nt,targ,ifppregtarg,pottarg,pretarg,gradtarg,ier)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()      


      call prin2('fmm time *',t2-t1,1)

      ntest = 100
      
      allocate(pot2(3,ntest),pre2(ntest),grad2(3,3,ntest))
      allocate(pottarg2(3,ntest),pretarg2(ntest),gradtarg2(3,3,ntest))

      do i = 1,ntest
         do j = 1,3
            pot2(j,i) = 0
            do k = 1,3
               grad2(k,j,i) = 0
            enddo
         enddo
         pre2(i) = 0
      enddo

      do i = 1,ntest
         do j = 1,3
            pottarg2(j,i) = 0
            do k = 1,3
               gradtarg2(k,j,i) = 0
            enddo
         enddo
         pretarg2(i) = 0
      enddo

      istress = 1
      thresh = 1d-15
      call st3ddirectstokstrsg(nd,source,stoklet,istress,
     1     strslet,strsvec,ns,source,ntest,pot2,pre2,grad2,thresh)
      call st3ddirectstokstrsg(nd,source,stoklet,istress,
     1     strslet,strsvec,ns,targ,ntest,pottarg2,pretarg2,gradtarg2,
     2     thresh)


      derr = 0
      drel = 0
      derrg = 0
      drelg = 0
      do i = 1,ntest
         derr = derr + (pot(1,i)-pot2(1,i))**2
         drel = drel + pot2(1,i)**2
         derr = derr + (pot(2,i)-pot2(2,i))**2
         drel = drel + pot2(2,i)**2
         derr = derr + (pot(3,i)-pot2(3,i))**2
         drel = drel + pot2(3,i)**2

         do j = 1,3
            do k = 1,3
               derrg = derrg + (grad(k,j,i)-grad2(k,j,i))**2
               drelg = drelg + grad2(k,j,i)**2
            enddo
         enddo
         
      enddo

      relerr = sqrt(derr/(ns*drel))
      relerrg = sqrt(derrg/(ns*drelg))
      if (relerr .lt. eps) ipass(1) = 1
      if (relerrg .lt. eps) ipass(2) = 1
      

      call prin2('rel err pot srcs *',relerr,1)
      call prin2('rel err grad srcs *',relerrg,1)

      
      derr = 0
      drel = 0
      do i = 1,ntest
         derr = derr + (pre(i)-pre2(i))**2
         drel = drel + pre2(i)**2
      enddo

      relerr = sqrt(derr/(ns*drel))
      if (relerr .lt. eps) ipass(3) = 1
      
      call prin2('rel err pre srcs *',relerr,1)

      derr = 0
      drel = 0

      derrg = 0
      drelg = 0
      
      do i = 1,ntest
         derr = derr + (pottarg(1,i)-pottarg2(1,i))**2
         drel = drel + pottarg2(1,i)**2
         derr = derr + (pottarg(2,i)-pottarg2(2,i))**2
         drel = drel + pottarg2(2,i)**2
         derr = derr + (pottarg(3,i)-pottarg2(3,i))**2
         drel = drel + pottarg2(3,i)**2

         do j = 1,3
            do k = 1,3
               derrg = derrg + (gradtarg(k,j,i)-gradtarg2(k,j,i))**2
               drelg = drelg + gradtarg2(k,j,i)**2
            enddo
         enddo
      enddo
      
      relerr = sqrt(derr/(nt*drel))
      relerrg = sqrt(derrg/(nt*drelg))
      if (relerr .lt. eps) ipass(4) = 1
      if (relerrg .lt. eps) ipass(5) = 1

      call prin2('rel err pot targs *',relerr,1)
      call prin2('rel err grad targs *',relerrg,1)

      derr = 0
      drel = 0
      do i = 1,ntest
         derr = derr + (pretarg(i)-pretarg2(i))**2
         drel = drel + pretarg2(i)**2
      enddo

      relerr = sqrt(derr/(nt*drel))
      if (relerr .lt. eps) ipass(6) = 1

      call prin2('rel err pre targs *',relerr,1)
      
      isum = 0
      do i=1,ntests
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in stfmm3d testing suite'
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in stfmm3d testing suite'
      close(33)
      


      stop
      end
c
c
c
c
c
c
      subroutine errprint(pot,opot,fld,ofld,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot,opot,fld(3),ofld(3)
      real *8 errs(2)
 1000  format(4D15.5) 
      err = 0
      ddd = 0
      err = err + abs(fld(1)-ofld(1))**2
      err = err + abs(fld(2)-ofld(2))**2
      err = err + abs(fld(3)-ofld(3))**2
      ddd = ddd + abs(ofld(1))**2
      ddd = ddd + abs(ofld(2))**2
      ddd = ddd + abs(ofld(3))**2
      err = sqrt(err)
      ddd = sqrt(ddd)

      err1 = abs(pot-opot)/abs(opot)
      write(*,'(a,e11.4,a,e11.4)') 
     1     'pot error=',err1,'   grad error=',err/ddd

      errs(1) = err1 
      errs(2) = err/ddd

      return
      end
c
