c
c      TESTING SUITE FOR STOKES KERNELS:  June 18th, 2020
c
c
      implicit real *8 (a-h,o-z)
      real *8 ztrg(3),source(3,20),targ(3,20),targ2(3,20)
      real *8 errgrad(3,20), pert(3), dmu(3), dnu(3)
      
      real *8, allocatable :: pot(:,:), pre(:), grad(:,:,:)
      real *8, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:)      

      real *8, allocatable :: potl(:,:), gradl(:,:,:)
      
      real *8, allocatable :: stoklet(:,:), strslet(:,:)
      real *8, allocatable :: strsvec(:,:)
      real *8, allocatable :: charge(:,:), dipvec(:,:,:)
      
      integer ipass(5)
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
      write(*,*) "=========================================="
      write(*,*) "Testing suite for stokkernels"

      open(unit=33,file='print_testres.txt',access='append')

      ntest = 4
      do i=1,ntest
        ipass(i) = 0 
      enddo


      ns = 10
      nt = 12

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

c     test gradient of stokeslet formula
      
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
      enddo

      errmax = 0
      
      call st3ddirectstokg(nd1,source,stoklet,
     1     ns,targ,nt,pot,pre,grad,thresh)
      call st3ddirectstokg(nd1,source,stoklet,
     1     ns,targ2,nt,pot2,pre2,grad2,thresh)

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

c      call prin2('err grad (just stokeslet)*',errgrad,3*nt)

      if (errmax .lt. 1d-6) ipass(1) = 1

c     test gradient of stokeslet and stresslet type I formula
      
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
      enddo


      errmax = 0
      
      istress = 1
      call st3ddirectstokstrsg(nd1,source,stoklet,istress,
     1     strslet,strsvec,ns,targ,nt,pot,pre,grad,thresh)
      call st3ddirectstokstrsg(nd1,source,stoklet,istress,
     1     strslet,strsvec,ns,targ2,nt,pot2,pre2,grad2,thresh)

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

c      call prin2('err grad (stokeslet + type I stresslet)*',
c     1     errgrad,3*nt)
         
      if (errmax .lt. 1d-6) ipass(2) = 1



c
c     figure out Stokeslet - Laplace relations

      ndl = 4
      
      allocate(potl(ndl,nt),gradl(ndl,3,nt))
      allocate(charge(ndl,ns),dipvec(ndl,3,ns))

c     clear output vars (these are all add type routines)
      
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
         do j = 1,ndl
            potl(j,i) = 0
            gradl(j,1,i) = 0
            gradl(j,2,i) = 0
            gradl(j,3,i) = 0
         enddo
      enddo

c     convert to equivalent vector of laplace charges

      do i = 1,ns
         do l = 1,3
            charge(l,i) = stoklet(l,i)/2
         enddo
         l = 4
         pl = stoklet(1,i)*source(1,i) + stoklet(2,i)*source(2,i) +
     1        stoklet(3,i)*source(3,i)
         charge(l,i) = pl/2
      enddo

      call l3ddirectcg(ndl,source,charge,ns,targ,nt,potl,gradl,thresh)

c     unpack laplace to stokes

      do i = 1,nt

         do l = 1,3
            pot2(l,i) = pot2(l,i) + potl(l,i)

            pot2(1,i) = pot2(1,i) - targ(l,i)*gradl(l,1,i)
            pot2(2,i) = pot2(2,i) - targ(l,i)*gradl(l,2,i)
            pot2(3,i) = pot2(3,i) - targ(l,i)*gradl(l,3,i)

            pre(i) = pre(i) + 2*gradl(l,l,i)
         enddo

         l = 4
         pot2(1,i) = pot2(1,i) + gradl(l,1,i)
         pot2(2,i) = pot2(2,i) + gradl(l,2,i)
         pot2(3,i) = pot2(3,i) + gradl(l,3,i)
         
      enddo

c     get stokes directly again
      
      call st3ddirectstokg(nd1,source,stoklet,
     1     ns,targ,nt,pot,pre,grad,thresh)

      errmax = 0
      valmax = 0

      do i = 1,nt
         valmax = max(valmax,abs(pot(1,i)))
         valmax = max(valmax,abs(pot(2,i)))
         valmax = max(valmax,abs(pot(3,i)))
         
         errmax = max(errmax,abs(pot(1,i)-pot2(1,i)))
         errmax = max(errmax,abs(pot(2,i)-pot2(2,i)))
         errmax = max(errmax,abs(pot(3,i)-pot2(3,i)))
      enddo

      
c      call prin2('pot *',pot,3*nt)
c      call prin2('pot2 *',pot2,3*nt)      

      if (errmax .lt. 1d-6*valmax) ipass(3) = 1
      
c
c     figure out stresslet (type I) - Laplace relations

      ndl = 4
      
c     clear output vars (these are all add type routines)
      
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
         do j = 1,ndl
            potl(j,i) = 0
            gradl(j,1,i) = 0
            gradl(j,2,i) = 0
            gradl(j,3,i) = 0
         enddo
      enddo

c     convert to equivalent vector of laplace charges

      do i = 1,ns
         dmu(1) = strslet(1,i)
         dmu(2) = strslet(2,i)
         dmu(3) = strslet(3,i)
         dnu(1) = strsvec(1,i)
         dnu(2) = strsvec(2,i)
         dnu(3) = strsvec(3,i)
         do l = 1,3
            charge(l,i) = stoklet(l,i)/2
            dipvec(l,1,i) = -(dmu(l)*dnu(1) + dmu(1)*dnu(l))/2
            dipvec(l,2,i) = -(dmu(l)*dnu(2) + dmu(2)*dnu(l))/2
            dipvec(l,3,i) = -(dmu(l)*dnu(3) + dmu(3)*dnu(l))/2
         enddo
         l = 4
         pl = stoklet(1,i)*source(1,i) + stoklet(2,i)*source(2,i) +
     1        stoklet(3,i)*source(3,i)
         charge(l,i) = pl/2

         pl = (dmu(1)*source(1,i) + dmu(2)*source(2,i)
     1        + dmu(3)*source(3,i))
         pv = (dnu(1)*source(1,i) + dnu(2)*source(2,i)
     1        + dnu(3)*source(3,i))

         dipvec(l,1,i) = -(dmu(1)*pv + dnu(1)*pl)/2
         dipvec(l,2,i) = -(dmu(2)*pv + dnu(2)*pl)/2
         dipvec(l,3,i) = -(dmu(3)*pv + dnu(3)*pl)/2
         
      enddo

      call l3ddirectcdg(ndl,source,charge,dipvec,ns,targ,nt,potl,gradl,
     1     thresh)

c     unpack laplace to stokes

      do i = 1,nt

         do l = 1,3
            pot2(l,i) = pot2(l,i) + potl(l,i)

            pot2(1,i) = pot2(1,i) - targ(l,i)*gradl(l,1,i)
            pot2(2,i) = pot2(2,i) - targ(l,i)*gradl(l,2,i)
            pot2(3,i) = pot2(3,i) - targ(l,i)*gradl(l,3,i)

            pre(i) = pre(i) + 2*gradl(l,l,i)
         enddo

         l = 4
         pot2(1,i) = pot2(1,i) + gradl(l,1,i)
         pot2(2,i) = pot2(2,i) + gradl(l,2,i)
         pot2(3,i) = pot2(3,i) + gradl(l,3,i)
         
      enddo

c     get stokes directly again

      istress = 1
      nd1 = 1
      call st3ddirectstokstrsg(nd1,source,stoklet,istress,
     1     strslet,strsvec,ns,targ,nt,pot,pre,grad,thresh)

      errmax = 0
      valmax = 0

      do i = 1,nt
         valmax = max(valmax,abs(pot(1,i)))
         valmax = max(valmax,abs(pot(2,i)))
         valmax = max(valmax,abs(pot(3,i)))
         
         errmax = max(errmax,abs(pot(1,i)-pot2(1,i)))
         errmax = max(errmax,abs(pot(2,i)-pot2(2,i)))
         errmax = max(errmax,abs(pot(3,i)-pot2(3,i)))
      enddo

      
c      call prin2('pot *',pot,3*nt)
c      call prin2('pot2 *',pot2,3*nt)      

      if (errmax .lt. 1d-6*valmax) ipass(4) = 1

      
      
      isum = 0
      do i=1,ntest
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,' tests in stokkernels testing suite'
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',isum,
     1   ' out of ',ntest,' tests in stokkernels testing suite'
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
