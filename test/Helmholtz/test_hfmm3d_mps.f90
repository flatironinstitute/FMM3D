program test_hfmm3d_mp2loc
  implicit double precision (a-h,o-z)
  
  character(len=72) str1
  
  integer :: ns, nt, nc
  integer :: i,j,k,ntest,nd,idim
  integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
  integer :: ipass(18),len1,ntests,isum
  
  double precision :: eps, err, hkrand
  double precision, allocatable :: source(:,:), targ(:,:), centers(:,:)
  double precision, allocatable :: wlege(:)
  
  double complex :: eye, zk, ima
  double complex, allocatable :: charge(:,:)
  double complex, allocatable :: dipvec(:,:,:)
  double complex, allocatable :: pot(:,:),pottarg(:,:)
  double complex, allocatable :: grad(:,:,:),gradtarg(:,:,:)
  double complex, allocatable :: mpole(:,:,:,:)


  data eye/(0.0d0,1.0d0)/
  ima = (0,1)

  !
  ! initialize printing routine
  !
  call prini(6,13)

  nd = 3

  zk = 1.2d0 + eye*0.02d0

  ns = 20
  nc = ns
  nt = 19

  ntest = 10

  allocate(source(3,ns),targ(3,nt), centers(3,nc))
  allocate(charge(nd,ns),dipvec(nd,3,ns))
  allocate(pot(nd,ns))
  allocate(grad(nd,3,ns))

  allocate(pottarg(nd,nt))
  allocate(gradtarg(nd,3,nt))
  eps = 0.5d-9

  write(*,*) "=========================================="
  write(*,*) "Testing suite for hfmm3d_mps"
  write(*,'(a,e11.5)') "Requested precision = ",eps

  open(unit=33,file='print_testres.txt',access='append')

  ntests = 18
  do i=1,ntests
    ipass(i) = 0
  enddo

  !
  ! generate sources uniformly in the unit cube 
  !
  do i=1,ns
    source(1,i) = hkrand(0)**2
    source(2,i) = hkrand(0)**2
    source(3,i) = hkrand(0)**2

    do idim=1,nd

      charge(idim,i) = hkrand(0) + eye*hkrand(0)

      dipvec(idim,1,i) = hkrand(0) + eye*hkrand(0)
      dipvec(idim,2,i) = hkrand(0) + eye*hkrand(0)
      dipvec(idim,3,i) = hkrand(0) + eye*hkrand(0)

      pot(idim,i) = 0
      grad(idim,1,i) = 0
      grad(idim,2,i) = 0
      grad(idim,3,i) = 0
    enddo
  enddo



  !
  ! generate targets uniformly in the unit cube
  !
  do i=1,nt
    targ(1,i) = hkrand(0)
    targ(2,i) = hkrand(0)
    targ(3,i) = hkrand(0)

    do idim=1,nd

      pottarg(idim,i) = 0
      gradtarg(idim,1,i) = 0
      gradtarg(idim,2,i) = 0
      gradtarg(idim,3,i) = 0 
    enddo
  enddo


  !
  ! calculate min source separation and min
  ! target separation
  !
  ssep = 1000
  do i = 1,ns
    do j = 1,ns
      if (i .ne. j) then
        rs = 0
        do k = 1,3
          rs = rs + (source(k,i)-source(k,j))**2
        end do
        rs = sqrt(rs)
        call prin2('rs = *', rs, 1)
        if (rs .lt. ssep) ssep = rs
      end if

    end do
  end do

  call prin2('min source separation = *', ssep, 1)

  shift = ssep/2
  do i = 1,ns
    centers(1,i) = source(1,i) + shift
    centers(2,i) = source(2,i)
    centers(3,i) = source(3,i)
  end do

  call prin2('centers = *', centers, 3*nc)

  !
  ! now form a multipole expansion at each center
  !
  nterms = 5
  allocate( mpole(nd,0:nterms,-nterms:nterms,nc) )
 
  nlege = nterms + 10
  lw = 4*(nlege+1)**2
  allocate( wlege(lw) )

  call prinf('before ylgndrfwini, lw = *', lw, 1)
  call ylgndrfwini(nlege, wlege, lw, lused)
  call prinf('after ylgndrfwini, lused = *', lused, 1)

  len = nd*(nterms+1)*(2*nterms+1)*nc
  call zinitialize(len, mpole)
  
  ns1 = 1
  do i = 1,nc
    call h3dformmpc(nd, zk, rscale, source(1,i), charge(1,i), &
        ns1, centers(1,i), nterms, mpole(:,:,:,i), wlege, nlege)
  end do

  stop
  

  !
  ! now test source to source, charge, 
  ! with potentials
  !
  write(6,*) 'testing source to source'
  write(6,*) 'interaction: charges'
  write(6,*) 'output: potentials'
  write(6,*) 
  write(6,*) 

  call hfmm3d_s_c_p_vec(nd,eps,zk,ns,source,charge, &
      pot)

  ifcharge = 1
  ifdipole = 0
  ifpgh = 1
  ifpghtarg = 0


  call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,&
      dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg, &
      ntest,err)

  call prin2('l2 rel err=*',err,1)
  write(6,*)
  write(6,*)
  write(6,*) '================'
  if(err.lt.eps) ipass(1) = 1
  call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  if(err.ge.eps) write(33,*) str1(1:len1) 



  ! c

  ! c
  ! cc     now test source to source, charge, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source'
  !        write(6,*) 'interaction: charges'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_s_c_g_vec(nd,eps,zk,ns,source,charge,
  !      1      pot,grad)

  !        ifcharge = 1
  !        ifdipole = 0
  !        ifpgh = 2
  !        ifpghtarg = 0


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(2) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 



  ! c
  ! cc     now test source to source, dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to source'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_s_d_p_vec(nd,eps,zk,ns,source,dipvec,
  !      1      pot)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 1
  !        ifpghtarg = 0


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(3) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to source, dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_s_d_g_vec(nd,eps,zk,ns,source,dipvec,
  !      1      pot,grad)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 2
  !        ifpghtarg = 0


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(4) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 

  ! c
  ! cc     now test source to source, charge + dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to source'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_s_cd_p_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,pot)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 1
  !        ifpghtarg = 0 


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(5) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to source, charge + dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_s_cd_g_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,pot,grad)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 2
  !        ifpghtarg = 0


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(6) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 



  ! c
  ! cc     now test source to target, charge, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: charges'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_c_p_vec(nd,eps,zk,ns,source,charge,
  !      1      nt,targ,pottarg)

  !        ifcharge = 1
  !        ifdipole = 0
  !        ifpgh = 0
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(7) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to target, charge, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: charges'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_c_g_vec(nd,eps,zk,ns,source,charge,
  !      1      nt,targ,pottarg,gradtarg)

  !        ifcharge = 1
  !        ifdipole = 0
  !        ifpgh = 0
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(8) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 



  ! c
  ! cc     now test source to target, dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_d_p_vec(nd,eps,zk,ns,source,dipvec,
  !      1      nt,targ,pottarg)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 0
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(9) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to target, dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_d_g_vec(nd,eps,zk,ns,source,dipvec,
  !      1      nt,targ,pottarg,gradtarg)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 0
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(10) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 

  ! c
  ! cc     now test source to target, charge + dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_cd_p_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,nt,targ,pottarg)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 0
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(11) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to target, charge + dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to target'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_t_cd_g_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,nt,targ,pottarg,gradtarg)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 0
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(12) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 

  ! c
  ! cc     now test source to source + target, charge, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: charges'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_c_p_vec(nd,eps,zk,ns,source,charge,
  !      1      pot,nt,targ,pottarg)

  !        ifcharge = 1
  !        ifdipole = 0
  !        ifpgh = 1
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(13) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to source + target, charge, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_c_g_vec(nd,eps,zk,ns,source,charge,
  !      1      pot,grad,nt,targ,pottarg,gradtarg)

  !        ifcharge = 1
  !        ifdipole = 0
  !        ifpgh = 2
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(14) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 



  ! c
  ! cc     now test source to source + target, dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_d_p_vec(nd,eps,zk,ns,source,dipvec,
  !      1      pot,nt,targ,pottarg)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 1
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(15) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to source + target, dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_d_g_vec(nd,eps,zk,ns,source,dipvec,
  !      1      pot,grad,nt,targ,pottarg,gradtarg)

  !        ifcharge = 0
  !        ifdipole = 1
  !        ifpgh = 2
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(16) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 

  ! c
  ! cc     now test source to source + target, charge + dipole, 
  ! c      with potentials
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_cd_p_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,pot,nt,targ,pottarg)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 1
  !        ifpghtarg = 1


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(17) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  ! c
  ! cc     now test source to source + target, charge + dipole, 
  ! c      with potentials and gradients
  ! c
  !        write(6,*) 'testing source to source and target'
  !        write(6,*) 'interaction: charges + dipoles'
  !        write(6,*) 'output: potentials + gradients'
  !        write(6,*) 
  !        write(6,*) 

  !        call hfmm3d_st_cd_g_vec(nd,eps,zk,ns,source,charge,
  !      1      dipvec,pot,grad,nt,targ,pottarg,gradtarg)

  !        ifcharge = 1
  !        ifdipole = 1
  !        ifpgh = 2
  !        ifpghtarg = 2


  !        call comperr_vec(nd,zk,ns,source,ifcharge,charge,ifdipole,
  !      1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
  !      2   ntest,err)

  !        call prin2('l2 rel err=*',err,1)
  !        write(6,*)
  !        write(6,*)
  !        write(6,*) '================'
  !       if(err.lt.eps) ipass(18) = 1
  !       call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
  !       if(err.ge.eps) write(33,*) str1(1:len1) 


  !       isum = 0
  !       do i=1,ntests
  !         isum = isum+ipass(i)
  !       enddo

  !       write(*,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
  !      1   ' out of ',ntests,' tests in hfmm3d vec testing suite'
  !       write(33,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
  !      1   ' out of ',ntests,' tests in hfmm3d vec testing suite'
  !       close(33)


  stop
end program

! ----------------------------------------------------------
! 
! This is the end of the debugging code.
!
! ----------------------------------------------------------



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

