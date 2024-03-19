program helmkernels_omp_dr
  implicit double precision (a-h,o-z)

  double precision :: errs(100000)
  double precision, allocatable :: targ(:,:), sources(:,:)

  double complex :: zk, ima
  double complex, allocatable :: dipvec(:,:,:), charge(:,:)
  double complex, allocatable :: pot(:,:), grad(:,:,:), hess(:,:,:)
  double complex, allocatable :: pot2(:,:), grad2(:,:,:), hess2(:,:,:)


  !call prini(6,13)

  ima = (0,1)
  done = 1
  pi = 4*atan(done)
  
  zk  = dcmplx(1.1d0,0.4d0)
  thresh = 1.0d0-13
  
  
  nsrc = 10000
  ntarg = 10000
  
  allocate( sources(3,nsrc) )
  allocate( targ(3,ntarg) )

  nd = 4
  allocate( charge(nd,nsrc) )
  allocate( dipvec(nd,3,nsrc) )
  
  allocate( pot(nd,ntarg) )
  allocate( grad(nd,3,ntarg) )
  allocate( hess(nd,6,ntarg) )

  allocate( pot2(nd,ntarg) )
  allocate( grad2(nd,3,ntarg) )
  allocate( hess2(nd,6,ntarg) )

  ! generate some random data
  do i = 1,nsrc
     do j = 1,3
        sources(j,i) = cos(1000*pi/7.3d0 * i*j)
     end do
     do k = 1,nd
        theta = sin(2000*pi/13.7d0 * k*i*3)
        phi = sin(3000*pi/17.9d0 * k*i*7)
        dipvec(k,1,i) = sin(theta)*cos(phi)
        dipvec(k,2,i) = sin(theta)*sin(phi)
        dipvec(k,3,i) = cos(theta)
        charge(k,i) = exp(ima*4000*pi/3.42d0 *k*i)/nsrc
     end do
  end do

  do i = 1,ntarg
     do j = 1,3
        targ(j,i) = cos(5000*pi/23.123d0 *i*j)
     end do
  end do


  nthreads = 12

  print *, 'number of sources = ', nsrc
  print *, 'number of targets = ', ntarg
  print *, 'vec param nd = ', nd
  
  print *
  print *
  print *, '- - - - - - - testing h3ddirectcp - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectcp(nd, zk, sources, charge, nsrc, targ, ntarg, pot, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectcp(nd, zk, sources, charge, nsrc, targ, ntarg, pot2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' err = ', err
  end do
  
  

  print *
  print *
  print *, '- - - - - - - testing h3ddirectcg - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectcg(nd, zk, sources, charge, nsrc, targ, ntarg, pot, &
       grad, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectcg(nd, zk, sources, charge, nsrc, targ, ntarg, pot2, &
       grad2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' pot err = ', err
     err = 0
     do j = 1,ntarg
        err = err + abs(grad(i,1,j) - grad2(i,1,j))**2 &
             + abs(grad(i,2,j) - grad2(i,2,j))**2 &
             + abs(grad(i,3,j) - grad2(i,3,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' grad err = ', err
  end do




  print *
  print *
  print *, '- - - - - - - testing h3ddirectdp - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectcp(nd, zk, sources, dipvec, nsrc, targ, ntarg, pot, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectcp(nd, zk, sources, dipvec, nsrc, targ, ntarg, pot2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' err = ', err
  end do



  print *
  print *
  print *, '- - - - - - - testing h3ddirectdg - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectdg(nd, zk, sources, dipvec, nsrc, targ, ntarg, pot, &
       grad, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectdg(nd, zk, sources, dipvec, nsrc, targ, ntarg, pot2, &
       grad2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' pot err = ', err
     err = 0
     do j = 1,ntarg
        err = err + abs(grad(i,1,j) - grad2(i,1,j))**2 &
             + abs(grad(i,2,j) - grad2(i,2,j))**2 &
             + abs(grad(i,3,j) - grad2(i,3,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' grad err = ', err
  end do



  print *
  print *
  print *, '- - - - - - - testing h3ddirectcdp - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectcdp(nd, zk, sources, charge, dipvec, nsrc, targ, ntarg, &
       pot, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectcdp(nd, zk, sources, charge, dipvec, nsrc, targ, ntarg, &
       pot2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' err = ', err
  end do



  print *
  print *
  print *, '- - - - - - - testing h3ddirectcdg - - - - - - -'

  call setzero(pot, grad, hess, nd, ntarg)
  call setzero(pot2, grad2, hess2, nd, ntarg)

  !call prin2('zk = *', zk, 2)
  !call prin2('sources = *', sources, 3*nd*nsrc)
  !call prin2('charge = *', charge, 2*nd*nsrc)
  !call prin2('targ = *', targ, 3*ntarg)

  
  call omp_set_num_threads(1)
  t0 = omp_get_wtime()
  call h3ddirectcdg(nd, zk, sources, charge, dipvec, nsrc, targ, ntarg, pot, &
       grad, thresh)
  t1 = omp_get_wtime()
  t = t1-t0

  print *, 'single thread timing = ', t
  
  call omp_set_num_threads(nthreads)
  t00 = omp_get_wtime()
  call h3ddirectcdg(nd, zk, sources, charge, dipvec, nsrc, targ, ntarg, pot2, &
       grad2, thresh)
  t11 = omp_get_wtime()

  tomp = t11-t00
  print *, 'multithread timing = ', tomp
  print *, 'speedup = ', t/tomp

  
  print *
  do i = 1,nd
     err = 0
     do j = 1,ntarg
        err = err + abs(pot(i,j) - pot2(i,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' pot err = ', err
     err = 0
     do j = 1,ntarg
        err = err + abs(grad(i,1,j) - grad2(i,1,j))**2 &
             + abs(grad(i,2,j) - grad2(i,2,j))**2 &
             + abs(grad(i,3,j) - grad2(i,3,j))**2
     end do
     err = sqrt(err)
     print *, 'i = ', i, ' grad err = ', err
  end do


  

  ! call setzero(pot,grad,hess,nd,nt)
  ! call h3ddirectcg(nd,zk,sources,charge,ns,ztarg,nt, &
  !      pot,grad,thresh)
  ! write(6,*) ' dir cp '
  ! do i = 1,nt
  !    write(6,*)pot(1,i)
  ! enddo
  ! write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! write(6,*) ' egrad  ', egrad
  ! v1 = grad(1,1,2)
  ! v2 = grad(1,1,3)
  ! zfd1 = (v1-v2)/(2*hh)
  ! v1 = grad(1,2,2)
  ! v2 = grad(1,2,3)
  ! zfd2 = (v1-v2)/(2*hh)
  ! v1 = grad(1,3,2)
  ! v2 = grad(1,3,3)
  ! zfd3 = (v1-v2)/(2*hh)
  ! write(6,*) ' FD deriv of p_x ', zfd1 
  ! write(6,*) ' FD deriv of p_y ', zfd2
  ! write(6,*) ' FD deriv of p_z ', zfd3

  ! !  call setzero(pot,grad,hess,nd,nt)
  ! !  call h3ddirectch(nd,zk,sources,charge,ns,ztarg,nt,
  ! ! 1            pot,grad,hess,thresh)
  ! !  write(6,*) ' dir cp '
  ! !  do i = 1,nt
  ! !     write(6,*)pot(1,i)
  ! !  enddo
  ! !  write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! !  egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! !  write(6,*) ' egrad  ', egrad
  ! !  if (icomp.eq.1) then
  ! !     write(6,*) ' dir pxx ', hess(1,1,1)
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pxz ', hess(1,5,1)
  ! !     exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
  ! !     exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
  ! !     exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
  ! !     write(6,*) ' exx  ', exx
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' exz  ', exz
  ! !  else if (icomp.eq.2) then
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pyy ', hess(1,2,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
  ! !     eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
  ! !     eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' eyy  ', eyy
  ! !     write(6,*) ' eyz  ', eyz
  ! !  else if (icomp.eq.3) then
  ! !     write(6,*) ' dir pxz ', hess(1,4,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     write(6,*) ' dir pzz ', hess(1,3,1)
  ! !     exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
  ! !     eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
  ! !     ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
  ! !     write(6,*) ' exz  ', exz
  ! !     write(6,*) ' eyz  ', eyz
  ! !     write(6,*) ' ezz  ', ezz
  ! !  endif


  ! write(6,*) ' directd  tests -------------------------------'

  ! call setzero(pot,grad,hess,nd,nt)
  ! call h3ddirectdp(nd,zk,sources, &
  !      dipvec,ns,ztarg,nt,pot,thresh)
  ! write(6,*) ' dir cp '
  ! do i = 1,nt
  !    write(6,*)pot(1,i)
  ! enddo
  ! zfd =  (pot(1,2) - pot(1,3))/(2*hh)
  ! write(6,*) ' FD deriv ', zfd
  ! call setzero(pot,grad,hess,nd,nt)
  ! call h3ddirectdg(nd,zk,sources, &
  !      dipvec,ns,ztarg,nt,pot,grad,thresh)
  ! write(6,*) ' dir cp '
  ! do i = 1,nt
  !    write(6,*)pot(1,i)
  ! enddo
  ! write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! write(6,*) ' egrad  ', egrad
  ! v1 = grad(1,1,2)
  ! v2 = grad(1,1,3)
  ! zfd1 = (v1-v2)/(2*hh)
  ! v1 = grad(1,2,2)
  ! v2 = grad(1,2,3)
  ! zfd2 = (v1-v2)/(2*hh)
  ! v1 = grad(1,3,2)
  ! v2 = grad(1,3,3)
  ! zfd3 = (v1-v2)/(2*hh)
  ! write(6,*) ' FD x deriv of p_x ', zfd1
  ! write(6,*) ' FD x deriv of p_y ', zfd2
  ! write(6,*) ' FD x deriv of p_z ', zfd3

  ! !  call setzero(pot,grad,hess,nd,nt)
  ! !  call h3ddirectdh(nd,zk,sources,
  ! ! 1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
  ! !  write(6,*) ' dir cp '
  ! !  do i = 1,nt
  ! !     write(6,*)pot(1,i)
  ! !  enddo
  ! !  write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! !  egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! !  write(6,*) ' egrad  ', egrad
  ! !  if (icomp.eq.1) then
  ! !     write(6,*) ' dir pxx ', hess(1,1,1)
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pxz ', hess(1,5,1)
  ! !     exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
  ! !     exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
  ! !     exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
  ! !     write(6,*) ' exx  ', exx
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' exz  ', exz
  ! !  else if (icomp.eq.2) then
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pyy ', hess(1,2,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
  ! !     eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
  ! !     eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' eyy  ', eyy
  ! !     write(6,*) ' eyz  ', eyz
  ! !  else if (icomp.eq.3) then
  ! !     write(6,*) ' dir pxz ', hess(1,4,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     write(6,*) ' dir pzz ', hess(1,3,1)
  ! !     exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
  ! !     eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
  ! !     ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
  ! !     write(6,*) ' exz  ', exz
  ! !     write(6,*) ' eyz  ', eyz
  ! !     write(6,*) ' ezz  ', ezz
  ! !  endif


  ! write(6,*) ' directcd  tests -------------------------------'

  ! call setzero(pot,grad,hess,nd,nt)
  ! call h3ddirectcdp(nd,zk,sources,charge, &
  !      dipvec,ns,ztarg,nt,pot,thresh)
  ! write(6,*) ' dir cp '
  ! do i = 1,nt
  !    write(6,*)pot(1,i)
  ! enddo
  ! zfd =  (pot(1,2) - pot(1,3))/(2*hh)
  ! write(6,*) ' FD deriv ', zfd
  ! call setzero(pot,grad,hess,nd,nt)
  ! call h3ddirectcdg(nd,zk,sources,charge, &
  !      dipvec,ns,ztarg,nt,pot,grad,thresh)
  ! write(6,*) ' dir cp '
  ! do i = 1,nt
  !    write(6,*)pot(1,i)
  ! enddo
  ! write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! write(6,*) ' egrad  ', egrad
  ! v1 = grad(1,1,2)
  ! v2 = grad(1,1,3)
  ! zfd1 = (v1-v2)/(2*hh)
  ! v1 = grad(1,2,2)
  ! v2 = grad(1,2,3)
  ! zfd2 = (v1-v2)/(2*hh)
  ! v1 = grad(1,3,2)
  ! v2 = grad(1,3,3)
  ! zfd3 = (v1-v2)/(2*hh)
  ! write(6,*) ' FD deriv of p_x ', zfd1
  ! write(6,*) ' FD deriv of p_y ', zfd2
  ! write(6,*) ' FD deriv of p_z ', zfd3


  ! !  call setzero(pot,grad,hess,nd,nt)
  ! !  call h3ddirectcdh(nd,zk,sources,charge,
  ! ! 1            dipvec,ns,ztarg,nt,pot,grad,hess,thresh)
  ! !  write(6,*) ' dir cp '
  ! !  do i = 1,nt
  ! !     write(6,*)pot(1,i)
  ! !  enddo
  ! !  write(6,*) ' analytic grad ', grad(1,icomp,1)
  ! !  egrad = abs( grad(1,icomp,1)- zfd)/abs(grad(1,icomp,1))
  ! !  write(6,*) ' egrad  ', egrad
  ! !  if (icomp.eq.1) then
  ! !     write(6,*) ' dir pxx ', hess(1,1,1)
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pxz ', hess(1,5,1)
  ! !     exx = abs( hess(1,1,1)- zfd1)/abs(hess(1,1,1))
  ! !     exy = abs( hess(1,4,1)- zfd2)/abs(hess(1,4,1))
  ! !     exz = abs( hess(1,5,1)- zfd3)/abs(hess(1,5,1))
  ! !     write(6,*) ' exx  ', exx
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' exz  ', exz
  ! !  else if (icomp.eq.2) then
  ! !     write(6,*) ' dir pxy ', hess(1,4,1)
  ! !     write(6,*) ' dir pyy ', hess(1,2,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     exy = abs( hess(1,4,1)- zfd1)/abs(hess(1,4,1))
  ! !     eyy = abs( hess(1,2,1)- zfd2)/abs(hess(1,2,1))
  ! !     eyz = abs( hess(1,6,1)- zfd3)/abs(hess(1,6,1))
  ! !     write(6,*) ' exy  ', exy
  ! !     write(6,*) ' eyy  ', eyy
  ! !     write(6,*) ' eyz  ', eyz
  ! !  else if (icomp.eq.3) then
  ! !     write(6,*) ' dir pxz ', hess(1,4,1)
  ! !     write(6,*) ' dir pyz ', hess(1,6,1)
  ! !     write(6,*) ' dir pzz ', hess(1,3,1)
  ! !     exz = abs( hess(1,5,1)- zfd1)/abs(hess(1,5,1))
  ! !     eyz = abs( hess(1,6,1)- zfd2)/abs(hess(1,6,1))
  ! !     ezz = abs( hess(1,3,1)- zfd3)/abs(hess(1,3,1))
  ! !     write(6,*) ' exz  ', exz
  ! !     write(6,*) ' eyz  ', eyz
  ! !     write(6,*) ' ezz  ', ezz
  ! !  endif

end program helmkernels_omp_dr





subroutine errl2(nd, n, x, y, errs)
  implicit double precision (a-h,o-z)
  double precision :: errs(nd)
  double complex :: x(nd,n), y(nd,n)

  do i = 1,nd
     errs(i) = 0
     do j = 1,n
        errs(i) = errs(i) + abs(x(i,j)-y(i,j))
     end do
     errs(i) = sqrt(errs(i))
  end do

  return
end subroutine errl2





subroutine setzero(pot,grad,hess,nd,nt)
  implicit real *8 (a-h,o-z)
  complex *16 pot(nd,nt),grad(nd,3,nt)
  complex *16 hess(nd,6,nt)

  do ii=1,nd
     do jj=1,nt
        pot(ii,jj) = 0.0d0
        grad(ii,1,jj) = 0.0d0
        grad(ii,2,jj) = 0.0d0
        grad(ii,3,jj) = 0.0d0
        hess(ii,1,jj) = 0.0d0
        hess(ii,2,jj) = 0.0d0
        hess(ii,3,jj) = 0.0d0
        hess(ii,4,jj) = 0.0d0
        hess(ii,5,jj) = 0.0d0
        hess(ii,6,jj) = 0.0d0
     enddo
  enddo
  return
end subroutine setzero
