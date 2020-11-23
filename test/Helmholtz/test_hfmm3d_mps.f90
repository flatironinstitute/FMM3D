program test_hfmm3d_mp2loc
  implicit double precision (a-h,o-z)
  
  character(len=72) str1
  
  integer :: ns, nt, nc
  integer :: i,j,k,ntest,nd,idim
  integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
  integer :: ipass(18),len1,ntests,isum
  integer, allocatable :: nterms(:), impole(:)
  
  double precision :: eps, err, hkrand, dnorms(1000), force(10)
  double precision, allocatable :: source(:,:), targ(:,:)
  double precision, allocatable :: centers(:,:)
  double precision, allocatable :: wlege(:), rscales(:)
  
  double complex :: eye, zk, ima
  double complex, allocatable :: charge(:,:)
  double complex, allocatable :: dipvec(:,:,:)
  double complex, allocatable :: pot(:,:), pot2(:,:), pottarg(:,:)
  double complex, allocatable :: grad(:,:,:),gradtarg(:,:,:)
  double complex, allocatable :: hess(:,:,:),hesstarg(:,:,:)
  double complex, allocatable :: mpole(:), local(:)


  data eye/(0.0d0,1.0d0)/
  ima = (0,1)
  done = 1
  pi = 4*atan(done)

  !
  ! initialize printing routine
  !
  call prini(6,13)

  nd = 3


  n1 = 3
  ns = n1**3
  nc = ns

  call prinf('ns = *', ns, 1)
  call prinf('nc = *', nc, 1)
  
  nt = 19

  allocate(source(3,ns),targ(3,nt), centers(3,nc))
  allocate(charge(nd,ns),dipvec(nd,3,ns))
  allocate(pot(nd,ns), pot2(nd,ns))
  allocate(grad(nd,3,ns))
  allocate(hess(nd,6,ns))

  allocate(pottarg(nd,nt))
  allocate(gradtarg(nd,3,nt))
  allocate(hesstarg(nd,6,nt))
  eps = 0.5d-9

  write(*,*) "=========================================="
  write(*,*) "Testing suite for hfmm3d_mps"
  write(*,'(a,e12.5)') "Requested precision = ",eps

  boxsize = 1
  dlam = 1/boxsize
  zk = 2*pi/dlam + eye*0.02d0

  call prin2('boxsize in wavelengths = *', boxsize, 1)
  call prin2('dlam = *', dlam, 1)
  call prin2('zk = *', zk, 2)

  !
  ! generate sources uniformly in the unit cube 
  !
  h = 1.0d0/(n1+1)
  ijk = 0
  do i = 1,n1
    do j = 1,n1
      do k = 1,n1
        ijk = ijk + 1
        source(1,ijk) = h*i
        source(2,ijk) = h*j
        source(3,ijk) = h*k
      end do
    end do
  end do
  
  
  dnorm = 0
  do i=1,ns
    !source(1,i) = hkrand(0)**2
    !source(2,i) = hkrand(0)**2
    !source(3,i) = hkrand(0)**2

    do idim=1,nd

      charge(idim,i) = hkrand(0) + eye*hkrand(0)
      dnorm = dnorm + abs(charge(idim,i))**2
      
      dipvec(idim,1,i) = hkrand(0) + eye*hkrand(0)
      dipvec(idim,2,i) = hkrand(0) + eye*hkrand(0)
      dipvec(idim,3,i) = hkrand(0) + eye*hkrand(0)

      pot(idim,i) = 0
      grad(idim,1,i) = 0
      grad(idim,2,i) = 0
      grad(idim,3,i) = 0
    enddo
  enddo

  dnorm = sqrt(dnorm)
  do i=1,ns
    do idim = 1,nd
      charge(idim,i) = charge(idim,i)/dnorm
    end do
  end do
  

  ! call prin2('min source separation = *', ssep, 1)
  
  shift = h/1000
  call prin2('shift = *', shift, 1)
  do i = 1,ns
    centers(1,i) = source(1,i) + shift
    centers(2,i) = source(2,i)
    centers(3,i) = source(3,i)
  end do

  !call prin2('centers = *', centers, 3*nc)

  !
  ! now form a multipole expansion at each center
  !
  allocate(nterms(nc), impole(nc))

  ntm = 7
  ntot = 0
  do i = 1,nc
    nterms(i) = ntm + 2*cos(pi/2*i)
    ntot = ntot + (nterms(i)+1)*(2*nterms(i)+1)
  end do

  allocate(mpole(nd*ntot))

  impole(1) = 1
  do i = 1,nc-1
    ilen = (nterms(i)+1)*(2*nterms(i)+1)
    impole(i+1) = impole(i) + nd*ilen
  end do

  
  nlege = 400
  lw = 5*(nlege+1)**2
  allocate( wlege(lw) )

  call prinf('before ylgndrfwini, lw = *', lw, 1)
  call ylgndrfwini(nlege, wlege, lw, lused)
  call prinf('after ylgndrfwini, lused = *', lused, 1)

  call zinitialize(nd*ntot, mpole)
  
  ns1 = 1
  rscale = 1
  sc = abs(zk)*shift
  if (sc .lt. 1) rscale = sc

  call prin2('rscale = *', rscale, 1)
  
  allocate(rscales(nc))
  do i = 1,nc
    rscales(i) = rscale
    call h3dformmpc(nd, zk, rscale, source(1,i), charge(1,i), &
        ns1, centers(1,i), nterms(i), mpole(impole(i)), &
        wlege, nlege)
  end do

  
  !
  ! do the direct calculation
  !
  thresh = 1.0d-15
  ifcharge = 1
  ifdipole = 0
  ifpgh = 1
  ntarg = 0
  ifpghtarg = 0
  call hfmm3d(nd, eps, zk, ns, source, ifcharge, &
      charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntarg, &
      targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  
  call prin2('via fmm, potential = *', pot, 10)

  
  allocate(local(nd*ntot))
  
  !
  ! now test source to source, charge, 
  ! with potentials
  !
  print *
  print *
  print *
  write(6,*) 'testing multipoles to locals'
  write(6,*) 'input: multipole expansions'
  write(6,*) 'output: local expansions'
  write(6,*) 
  write(6,*) 


  call hfmm3d_mps(nd, eps, zk, &
      nc, centers, rscales, nterms, mpole, impole, local,ier)

  call zinitialize(nd*nc, pot2)
  npts = 1
  do i = 1,nc
    call h3dtaevalp(nd, zk, rscales(i), &
        centers(1,i), local(impole(i)), &
        nterms(i), source(1,i), npts, pot2(1,i), &
        wlege, nlege)
  end do
  
  call prin2('from hfmm3d_mps, potential = *', pot2, 10)

  err = 0
  dnorm = 0
  do j = 1,nc
    do i = 1,nd
      err = err + abs(pot(i,j)-pot2(i,j))**2
      dnorm = dnorm + abs(pot(i,j))**2
      pot2(i,j) = pot2(i,j) - pot(i,j)
    end do
  end do

  !call prin2('diff = *', pot2, 2*nd*ns)
  !call prin2('err = *', err, 1)
  !call prin2('dnorm = *', dnorm, 1)
  
  err = sqrt(err/dnorm)
  call prin2('l2 rel err=*',err,1)

  open(unit=33,file='print_testres.txt',access='append')
  isuccess = 0
  ntest = 1
  if(err.lt.eps) isuccess = 1

  write(33,'(a,i1,a,i1,a)') 'Successfully completed ', &
    isuccess,' out of ',ntest,' tests in helm3d_mps testing suite'
  close(33)

  
  

  stop
end program



! ----------------------------------------------------------
! 
! This is the end of the debugging code.
!
! ----------------------------------------------------------
subroutine prinmp(str, mpole, nterms)
  implicit double precision (a-h,o-z)
  character (len=*) :: str
  double complex :: mpole(0:nterms, -nterms:nterms)
  double complex :: tmp(-10000:10000)

  print *, trim(str)

  do n = 0,nterms
    print *
    write(6,'(a,i2,a,i2,a,i2)') 'n = ', n, '  m = ', -n, '...', n
    do m = -nterms,nterms
      tmp(m) = mpole(n,m)
    end do
    call prin2('coefs = *', tmp(-n), 2*n+1)
  end do
  
  return
end subroutine prinmp



subroutine zinitialize(len, zs)
  implicit double precision (a-h,o-z)
  double complex :: zs(len)

  do i = 1,len
    zs(i) = 0
  end do
  return
end subroutine zinitialize





subroutine comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole, &
    dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg, &
    gradtarg,ntest,err)

  implicit none
  double complex zk
  integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg

  double precision source(3,*),targ(3,*)
  double complex dipvec(nd,3,*)
  double complex charge(nd,*)

  double complex pot(nd,*),pottarg(nd,*),grad(nd,3,*), &
      gradtarg(nd,3,*)

  integer i,j,ntest,nd,idim

  double precision err,ra

  double complex potex(nd,ntest),gradex(nd,3,ntest), &
      pottargex(nd,ntest),gradtargex(nd,3,ntest)

  double precision thresh

  err = 0 
  do i=1,ntest
    do idim=1,nd
      potex(idim,i) = 0
      pottargex(idim,i) = 0

      gradex(idim,1,i) = 0
      gradex(idim,2,i) = 0
      gradex(idim,3,i) = 0

      gradtargex(idim,1,i) = 0
      gradtargex(idim,2,i) = 0
      gradtargex(idim,3,i) = 0
    enddo
  enddo

  thresh = 1.0d-16

  if(ifcharge.eq.1.and.ifdipole.eq.0) then
    if(ifpgh.eq.1) then
      call h3ddirectcp(nd,zk,source,charge,ns,source,ntest, &
          potex,thresh)
    endif

    if(ifpgh.eq.2) then
      call h3ddirectcg(nd,zk,source,charge,ns,source,ntest, &
          potex,gradex,thresh)
    endif

    if(ifpghtarg.eq.1) then
      call h3ddirectcp(nd,zk,source,charge,ns,targ,ntest, &
          pottargex,thresh)
    endif

    if(ifpghtarg.eq.2) then
      call h3ddirectcg(nd,zk,source,charge,ns,targ,ntest, &
          pottargex,gradtargex,thresh)
    endif
  endif

  if(ifcharge.eq.0.and.ifdipole.eq.1) then
    if(ifpgh.eq.1) then
      call h3ddirectdp(nd,zk,source,dipvec, &
          ns,source,ntest,potex,thresh)
    endif

    if(ifpgh.eq.2) then
      call h3ddirectdg(nd,zk,source,dipvec, &
          ns,source,ntest,potex,gradex,thresh)
    endif

    if(ifpghtarg.eq.1) then
      call h3ddirectdp(nd,zk,source,dipvec, &
          ns,targ,ntest,pottargex,thresh)
    endif

    if(ifpghtarg.eq.2) then
      call h3ddirectdg(nd,zk,source,dipvec, &
          ns,targ,ntest,pottargex,gradtargex,thresh)
    endif
  endif

  if(ifcharge.eq.1.and.ifdipole.eq.1) then
    if(ifpgh.eq.1) then
      call h3ddirectcdp(nd,zk,source,charge,dipvec, &
          ns,source,ntest,potex,thresh)
    endif

    if(ifpgh.eq.2) then
      call h3ddirectcdg(nd,zk,source,charge,dipvec, &
          ns,source,ntest,potex,gradex,thresh)
    endif

    if(ifpghtarg.eq.1) then
      call h3ddirectcdp(nd,zk,source,charge,dipvec, &
          ns,targ,ntest,pottargex,thresh)
    endif

    if(ifpghtarg.eq.2) then
      call h3ddirectcdg(nd,zk,source,charge,dipvec, &
          ns,targ,ntest,pottargex,gradtargex,thresh)
    endif
  endif

  err = 0
  ra = 0

  if(ifpgh.eq.1) then
    do i=1,ntest
      do idim=1,nd
        ra = ra + abs(potex(idim,i))**2
        err = err + abs(pot(idim,i)-potex(idim,i))**2
      enddo
    enddo
  endif

  if(ifpgh.eq.2) then
    do i=1,ntest
      do idim=1,nd
        ra = ra + abs(potex(idim,i))**2
        ra = ra + abs(gradex(idim,1,i))**2
        ra = ra + abs(gradex(idim,2,i))**2
        ra = ra + abs(gradex(idim,3,i))**2

        err = err + abs(pot(idim,i)-potex(idim,i))**2
        err = err + abs(grad(idim,1,i)-gradex(idim,1,i))**2
        err = err + abs(grad(idim,2,i)-gradex(idim,2,i))**2
        err = err + abs(grad(idim,3,i)-gradex(idim,3,i))**2
      enddo
    enddo
  endif


  if(ifpghtarg.eq.1) then
    do i=1,ntest
      do idim=1,nd
        ra = ra + abs(pottargex(idim,i))**2
        err = err + abs(pottarg(idim,i)-pottargex(idim,i))**2
      enddo
    enddo
  endif

  if(ifpghtarg.eq.2) then
    do i=1,ntest
      do idim=1,nd
        ra = ra + abs(pottargex(idim,i))**2
        ra = ra + abs(gradtargex(idim,1,i))**2
        ra = ra + abs(gradtargex(idim,2,i))**2
        ra = ra + abs(gradtargex(idim,3,i))**2

        err = err + abs(pottarg(idim,i)-pottargex(idim,i))**2
        err = err + abs(gradtarg(idim,1,i)-gradtargex(idim,1,i))**2
        err = err + abs(gradtarg(idim,2,i)-gradtargex(idim,2,i))**2
        err = err + abs(gradtarg(idim,3,i)-gradtargex(idim,3,i))**2
      enddo
    enddo
  endif

  err = sqrt(err/ra)
  return
end subroutine comperr_vec

