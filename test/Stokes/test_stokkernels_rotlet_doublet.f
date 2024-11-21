c
c      TESTING SUITE FOR STOKES KERNELS rotlet and doublet
c
      implicit none
      real *8 ztrg(3),source(3,20),targ(3,20),targ2(3,20)
      real *8 errgrad(3,20), pert(3), dmu(3), dnu(3)
      real *8 rlet(3), rvec(3)

      real *8, allocatable :: pot(:,:), pre(:), grad(:,:,:)
      real *8, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:)
      real *8, allocatable :: pot3(:,:), pre3(:), grad3(:,:,:)
      real *8, allocatable :: pot4(:,:), pre4(:), grad4(:,:,:)
      real *8, allocatable :: pot_src(:,:),pre_src(:),grad_src(:,:,:)
      real *8, allocatable :: pot_diff(:,:), grad_diff(:,:,:)

      real *8, allocatable :: potl(:,:), gradl(:,:,:)

      real *8, allocatable :: stoklet(:,:), strslet(:,:)
      real *8, allocatable :: strsvec(:,:)
      real *8, allocatable :: charge(:,:), dipvec(:,:,:)
      real *8, allocatable :: rotstr(:,:), rotvec(:,:)
      real *8, allocatable :: doubstr(:,:), doubvec(:,:)
      real *8 df,done,pi,h,thresh,errmax,gave1,gave2,gave3


      integer *8 ipass(4)
      integer *8 istoklet, istress, irotlet, idoublet
      integer *8 i,j,k,nd,ns,nt,ntest,nd1,ier,isum
      integer *8 ifppgreg, ifppregtarg
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
      write(*,*) "=========================================="
      write(*,*) "Testing suite for stokkernels rotlet and doublet"

      open(unit=33,file='print_testres.txt',access='append')

      ntest = 4
      do i=1,ntest
        ipass(i) = 0
      enddo


      ns = 12
      nt = 12

      allocate(stoklet(3,ns),strslet(3,ns),strsvec(3,ns))
      allocate(rotstr(3,ns),rotvec(3,ns))
      allocate(doubstr(3,ns),doubvec(3,ns))

      do i = 1,ns
         source(1,i) = cos(i*done)
         source(2,i) = cos(10*i*done)
         source(3,i) = cos(100*i*done)
         stoklet(1,i) = 0.0d0
         stoklet(2,i) = 0.0d0
         stoklet(3,i) = 0.0d0
         strslet(1,i) = 0.0d0
         strslet(2,i) = 0.0d0
         strslet(3,i) = 0.0d0
         strsvec(1,i) = 0.0d0
         strsvec(2,i) = 0.0d0
         strsvec(3,i) = 0.0d0
         rotstr(1,i) = sin(4*i*done)
         rotstr(2,i) = sin(44*i*done)
         rotstr(3,i) = sin(444*i*done)
         rotvec(1,i) = sin(5*i*done)
         rotvec(2,i) = sin(55*i*done)
         rotvec(3,i) = sin(555*i*done)
         doubstr(1,i) = sin(3*i*done)
         doubstr(2,i) = sin(33*i*done)
         doubstr(3,i) = sin(333*i*done)
         doubvec(1,i) = sin(4*i*done)
         doubvec(2,i) = sin(44*i*done)
         doubvec(3,i) = sin(444*i*done)
      enddo

c      call prin2('src *',source,3*ns)
c      call prin2('stoklet *',stoklet,3*ns)
c      call prin2('strslet *',strslet,3*ns)
c      call prin2('strsvec *',strsvec,3*ns)
c

      h = 1d-5
      pert(1) = cos(1000*done)
      pert(2) = cos(2000*done)
      pert(3) = cos(3000*done)

      do i = 1,nt
         targ(1,i) = cos(12*i*done)
         targ(2,i) = cos(123*i*done)
         targ(3,i) = cos(1234*i*done)
         targ2(1,i) = targ(1,i) + h*pert(1)
         targ2(2,i) = targ(2,i) + h*pert(2)
         targ2(3,i) = targ(3,i) + h*pert(3)
      enddo

      allocate(pot(3,nt),pre(nt),grad(3,3,nt))
      allocate(pot2(3,nt),pre2(nt),grad2(3,3,nt))
      allocate(pot3(3,nt),pre3(nt),grad3(3,3,nt))
      allocate(pot4(3,nt),pre4(nt),grad4(3,3,nt))
      allocate(pot_src(3,ns),pre_src(ns),grad_src(3,3,ns))
      allocate(pot_diff(3,nt))
      allocate(grad_diff(3,3,nt))

c     test gradient of rotlet formula

      thresh = 1d-15
      nd1 = 1

      do i = 1,nt
         pre(i) = 0
         pot(1,i) = 0
         pot(2,i) = 0
         pot(3,i) = 0
         grad(1,1,i) = 0
         grad(2,1,i) = 0
         grad(3,1,i) = 0
         grad(1,2,i) = 0
         grad(2,2,i) = 0
         grad(3,2,i) = 0
         grad(1,3,i) = 0
         grad(2,3,i) = 0
         grad(3,3,i) = 0
         pre2(i) = 0
         pot2(1,i) = 0
         pot2(2,i) = 0
         pot2(3,i) = 0
         grad2(1,1,i) = 0
         grad2(2,1,i) = 0
         grad2(3,1,i) = 0
         grad2(1,2,i) = 0
         grad2(2,2,i) = 0
         grad2(3,2,i) = 0
         grad2(1,3,i) = 0
         grad2(2,3,i) = 0
         grad2(3,3,i) = 0
         pre3(i) = 0
         pot3(1,i) = 0
         pot3(2,i) = 0
         pot3(3,i) = 0
         grad3(1,1,i) = 0
         grad3(2,1,i) = 0
         grad3(3,1,i) = 0
         grad3(1,2,i) = 0
         grad3(2,2,i) = 0
         grad3(3,2,i) = 0
         grad3(1,3,i) = 0
         grad3(2,3,i) = 0
         grad3(3,3,i) = 0
         pre4(i) = 0
         pot4(1,i) = 0
         pot4(2,i) = 0
         pot4(3,i) = 0
         grad4(1,1,i) = 0
         grad4(2,1,i) = 0
         grad4(3,1,i) = 0
         grad4(1,2,i) = 0
         grad4(2,2,i) = 0
         grad4(3,2,i) = 0
         grad4(1,3,i) = 0
         grad4(2,3,i) = 0
         grad4(3,3,i) = 0
      enddo

      do i = 1,ns
         pre_src(i) = 0
         pot_src(1,i) = 0
         pot_src(2,i) = 0
         pot_src(3,i) = 0
         grad_src(1,1,i) = 0
         grad_src(2,1,i) = 0
         grad_src(3,1,i) = 0
         grad_src(1,2,i) = 0
         grad_src(2,2,i) = 0
         grad_src(3,2,i) = 0
         grad_src(1,3,i) = 0
         grad_src(2,3,i) = 0
         grad_src(3,3,i) = 0
      enddo

      errmax = 0

      istoklet = 0
      istress = 0
      irotlet = 1
      idoublet = 1
      call st3ddirectstokstrsrotdoubg(nd1,source,stoklet,istress,
     1     strslet,strsvec,irotlet,rotstr,rotvec,
     2     idoublet,doubstr,doubvec,
     3     ns,targ,nt,pot,pre,grad,thresh)
      call st3ddirectstokstrsrotdoubg(nd1,source,stoklet,istress,
     1     strslet,strsvec,irotlet,rotstr,rotvec,
     2     idoublet,doubstr,doubvec,
     3     ns,targ2,nt,pot2,pre2,grad2,thresh)


      do i = 1,nt
         df = pot2(1,i) - pot(1,i)
         gave1 = (grad(1,1,i) + grad2(1,1,i))/2.0d0
         gave2 = (grad(2,1,i) + grad2(2,1,i))/2.0d0
         gave3 = (grad(3,1,i) + grad2(3,1,i))/2.0d0         
         errgrad(1,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot2(2,i) - pot(2,i)
         gave1 = (grad(1,2,i) + grad2(1,2,i))/2.0d0
         gave2 = (grad(2,2,i) + grad2(2,2,i))/2.0d0
         gave3 = (grad(3,2,i) + grad2(3,2,i))/2.0d0         
         errgrad(2,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot2(3,i) - pot(3,i)
         gave1 = (grad(1,3,i) + grad2(1,3,i))/2.0d0
         gave2 = (grad(2,3,i) + grad2(2,3,i))/2.0d0
         gave3 = (grad(3,3,i) + grad2(3,3,i))/2.0d0         
         errgrad(3,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         errmax = max(errmax,errgrad(1,i))
         errmax = max(errmax,errgrad(2,i))
         errmax = max(errmax,errgrad(3,i))                  
         
      enddo

      call prin2('err grad (rotlet+doublet) * direct',errgrad,3*nt)

      if (errmax .lt. 1d-6) ipass(1) = 1

      errmax = 0
      ifppgreg = 1
      ifppregtarg = 3
      call stfmm3d(nd1,1d-12,ns,source,istoklet,stoklet,
     1     istress,strslet,strsvec,
     2     irotlet,rotstr,rotvec,
     3     idoublet,doubstr,doubvec,
     4     ifppgreg,pot_src,pre_src,grad_src,
     5     nt,targ,ifppregtarg,pot3,pre3,grad3,ier)
      call stfmm3d(nd1,1d-12,ns,source,istoklet,stoklet,
     1     istress,strslet,strsvec,
     2     irotlet,rotstr,rotvec,
     3     idoublet,doubstr,doubvec,
     4     ifppgreg,pot_src,pre_src,grad_src,
     5     nt,targ2,ifppregtarg,pot4,pre4,grad4,ier)

      do i = 1,nt
         df = pot4(1,i) - pot3(1,i)
         gave1 = (grad3(1,1,i) + grad4(1,1,i))/2.0d0
         gave2 = (grad3(2,1,i) + grad4(2,1,i))/2.0d0
         gave3 = (grad3(3,1,i) + grad4(3,1,i))/2.0d0
         errgrad(1,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot4(2,i) - pot3(2,i)
         gave1 = (grad3(1,2,i) + grad4(1,2,i))/2.0d0
         gave2 = (grad3(2,2,i) + grad4(2,2,i))/2.0d0
         gave3 = (grad3(3,2,i) + grad4(3,2,i))/2.0d0
         errgrad(2,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot4(3,i) - pot3(3,i)
         gave1 = (grad3(1,3,i) + grad4(1,3,i))/2.0d0
         gave2 = (grad3(2,3,i) + grad4(2,3,i))/2.0d0
         gave3 = (grad3(3,3,i) + grad4(3,3,i))/2.0d0
         errgrad(3,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         errmax = max(errmax,errgrad(1,i))
         errmax = max(errmax,errgrad(2,i))
         errmax = max(errmax,errgrad(3,i))
      enddo

      call prin2('err grad (rotlet+doublet) * fmm',errgrad,3*nt)

      if (errmax .lt. 1d-6) ipass(2) = 1

      errmax = 0
      do i = 1,nt
         do j = 1,3
            pot_diff(j,i) = abs(pot3(j,i) - pot(j,i))
            errmax = max(errmax,pot_diff(j,i))
         enddo
      enddo

      call prin2('pot diff',pot_diff,3*nt)

      if (errmax .lt. 1d-10) ipass(3) = 1

      errmax = 0
      do i = 1,nt
         do j = 1,3
            do k = 1,3
               grad_diff(k,j,i) = abs(grad3(k,j,i) - grad(k,j,i))
               errmax = max(errmax,grad_diff(k,j,i))
            enddo
         enddo
      enddo

      call prin2('grad diff',grad_diff,3*3*nt)

      if (errmax .lt. 1d-10) ipass(4) = 1

      isum = 0
      do i=1,ntest
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,
     2   ' tests in stokkernels rotlet doublet testing suite'
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,
     2   ' tests in stokkernels rotlet doublet testing suite'
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
