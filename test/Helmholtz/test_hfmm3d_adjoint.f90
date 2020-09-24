program test_hfmm3d_adjoint

  implicit double precision (a-h,o-z)

  double precision :: center0(3), center1(3)
  double precision :: xnodes(10000), wts(10000)
  
  double complex :: zk, ima, zrat, zrat1

  double complex, allocatable :: mpole0(:,:,:), mpole1(:,:,:)
  double complex, allocatable :: mpmpmat(:,:,:,:)
  double complex, allocatable :: loclocmat(:,:,:,:)
  

  ima = (0,1)
  done = 1
  pi = 4*atan(done)

  !
  ! initialize printing routine
  !
  call prini(6,13)

  ! set the vectorization parameter
  nd = 1

  ! --------------
  ! construct the transpose operator for mp2mp
  !
  nterms0 = 20
  nterms1 = nterms0
  allocate( mpole0(nd,0:nterms0,-nterms0:nterms0) )
  allocate( mpole1(nd,0:nterms1,-nterms1:nterms1) )

  do m = -nterms0,nterms0
    do l = 0,nterms0
      do i = 1,nd
        mpole0(i,l,m) = 0
      end do
    end do
  end do

  do m = -nterms1,nterms1
    do l = 0,nterms1
      do i = 1,nd
        mpole1(i,l,m) = 0
      end do
    end do
  end do
  
  center0(1) = .5d0
  center0(2) = .5d0
  center0(3) = .5d0

  center1(1) = 1
  center1(2) = 1
  center1(3) = 1

  dlam = 1
  zk = 2*pi*dlam + ima/10000

  ! now construct the translation matrix directly

  sc0 = 1
  sc1 = 1
  radius = sqrt(3.0d0)*2

  nquad = nterms1*10
  nquad = max(6,nquad)
  ifinit = 1
  call legewhts(nquad, xnodes, wts, ifinit)

  allocate( mpmpmat(0:nterms1,-nterms1:nterms1,0:nterms0, &
      -nterms0:nterms0) )

  do m0 = -nterms0,nterms0
    do l0 = 0,nterms0
      do m1 = -nterms1,nterms1
        do l1 = 0,nterms1
          mpmpmat(l1,m1,l0,m0) = 0
        end do
      end do
    end do
  end do
  
  
  do l0 = 0,nterms0
    do m0 = -l0,l0
      mpole0(1,l0,m0) = 1
      call h3dmpmp(nd, zk, sc0, center0, mpole0, nterms0, &
          sc1, center1, mpole1, nterms1, &
          radius, xnodes, wts, nquad)
      mpole0(1,l0,m0) = 0

      do l1 = 0,nterms1
        do m1 = -l1,l1
          mpmpmat(l1,m1,l0,m0) = mpole1(1,l1,m1)
        end do
      end do

    end do
  end do

  ll0 = 1
  mm0 = 1
  ll1 = 2
  mm1 = 1
  call prin2('from mpmpmat, entry = *', mpmpmat(ll1,mm1,ll0,mm0), 2)
  !call prin2('from mpmpmat, entry+1 = *', mpmpmat(ll1+1,mm1,ll0,mm0), 2)


  ! now do the locloc mat

  allocate( loclocmat(0:nterms0,-nterms0:nterms0,0:nterms1, &
      -nterms1:nterms1) )

  do m1 = -nterms1,nterms1
    do l1 = 0,nterms1
      do m0 = -nterms0,nterms0
        do l0 = 0,nterms0
          loclocmat(l0,m0,l1,m1) = 0
        end do
      end do
    end do
  end do

  ! zero out the mpole expansions again
  do m = -nterms0,nterms0
    do l = 0,nterms0
      do i = 1,nd
        mpole0(i,l,m) = 0
      end do
    end do
  end do

  do m = -nterms1,nterms1
    do l = 0,nterms1
      do i = 1,nd
        mpole1(i,l,m) = 0
      end do
    end do
  end do

  
  ! construct the locloc matrix
  radius = sqrt(3.0d0)
  do l1 = 0,nterms1
    do m1 = -l1,l1
      mpole0(1,l1,m1) = 1
      !call h3dlocloc(nd, zk, sc1, center1, mpole1, nterms1, &
      !    sc0, center0, mpole0, nterms0, &
      !    radius, xnodes, wts, nquad)
      call h3dlocloc(nd, zk, sc0, center0, mpole0, nterms0, &
          sc1, center1, mpole1, nterms1, &
          radius, xnodes, wts, nquad)
      mpole0(1,l1,m1) = 0

      do l0 = 0,nterms0
        do m0 = -l0,l0
          loclocmat(l0,m0,l1,m1) = mpole1(1,l0,m0)
        end do
      end do
      
    end do
  end do

  call prin2('from loclocmat, entry = *',  &
      loclocmat(ll1,mm1,ll0,mm0), 2)

  zrat = mpmpmat(ll1,mm1,ll0,mm0)/loclocmat(ll1,mm1,ll0,mm0)
      
  ! call prin2('from loclocmat, transpose entry = *', &
  !     loclocmat(ll0,mm0,ll1,mm1), 2)

  ! call prin2('from loclocmat, transpose entry+1 = *', &
  !     loclocmat(ll0,mm0,ll1+1,mm1), 2)

  !zrat = mpmpmat(ll1,mm1,ll0,mm0)/loclocmat(ll0,mm0,ll1,mm1)
  !zrat1 = mpmpmat(ll1+1,mm1,ll0,mm0)/loclocmat(ll0,mm0,ll1+1,mm1)

  call prin2('ratio = *', zrat, 2)
  !call prin2('ratio1 = *', zrat1, 2)
  
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

