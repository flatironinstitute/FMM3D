!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file contains interfaces for the maxwell fmm3d and 
! direct calculation. It contains the two subroutines:
!
!   emfmm3d - Maxwell fmm3d with nd densities
!
!   em3ddirect - Maxwell direct calculation with nd
!    densities
!
!   (input sources and output requests are controleld by flags)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine emfmm3d(nd,eps,zk,ns,source,ifa_vect,a_vect,&
 ifb_vect,b_vect,iflambda,lambda,nt,targets,ifE,E,ifcurlE,curlE,&
 ifdivE,divE)
implicit none

!jmc_edc

!  This function computes 
!
!    E = curl S_{k}[a_vect] + S_{k}[b_vect] + grad S_{k}[lambda]  -- (1)
!  using the vector Helmholtz fmm.
!
!  The subroutine also computes divE, curlE 
!  with appropriate flags
!
!  input:
!    nd - integer
!      number of densities (a_vect,b_vect,lambda)
!
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter
!
!    ns - integer
!      number of sources
!
!    source - real *8 (3,ns)
!      location of the sources
!
!    ifa_vect - integer
!     Contribution due to curl S_{k}[a_vect] is included if
!     ifa_vect = 1
!
!    a_vect - complex *16(nd,3,ns)
!      a vector source
!
!    ifb_vect - integer
!      Contribution due to S_{k}[b_vect] is included if
!      ifb_vect = 1
!
!    b_vect - complex *16(nd,3,ns)
!      b vector source
!
!    iflambda - integer
!      Contribution due to grad S_{k} [lambda] is included
!      if iflambda = 1
!
!    lambda - complex *16(nd,ns)
!      lambda source
!
!    nt - integer
!      number of targets
!
!    targets - real *8 (3,nt)
!      location of the targets
!
!    ifE - integer
!      E is returned at the target locations if ifE = 1
!
!    ifcurlE - integer
!      curl E is returned at the target locations if 
!      ifcurlE = 1
!
!    ifdivE - integer
!      div E is returned at the target locations if 
!      ifdivE = 1
!
!  output:
!
!    E - complex *16(nd,3,nt)
!      E field defined in (1) above
!
!    curlE - complex *16(nd,3,nt)
!      curl of E field
!
!    divE - complex *16(nd,nt)
!      divergence of E
!

  !List of calling arguments
  double precision, intent(in) :: eps
  double complex, intent(in) :: zk
  integer, intent(in) :: nd,ns,nt
  double precision, intent(in) :: source(3,ns)
  integer, intent(in) :: ifa_vect
  double complex, intent(in) :: a_vect(nd,3,*)
  integer, intent(in) :: ifb_vect
  double complex, intent(in) :: b_vect(nd,3,*)
  integer, intent(in) :: iflambda
  double complex, intent(in) :: lambda(nd,*)
  integer, intent(in) :: ifE
  double complex, intent(out) :: E(nd,3,*)
  integer, intent(in) :: ifcurlE
  double complex, intent(out) :: curlE(nd,3,*)
  integer, intent(in) :: ifdivE
  double complex, intent(out) :: divE(nd,*)
  double precision, intent(in) :: targets(3,nt)

  !List of local variables
  double complex, allocatable :: sigma_vect(:,:,:)
  double complex, allocatable :: dipvect_vect(:,:,:,:)
  double complex, allocatable :: Etmp(:,:,:)
  double complex, allocatable :: gradE_vect(:,:,:,:)

  double complex, allocatable :: pot(:),grad(:,:),hess(:,:)
  double complex, allocatable :: hesstarg(:,:)


  integer i,j,nd0,l,m
  integer ifcharge,ifdipole,ifpgh,ifpghtarg

  !!Initialize sources
  allocate(sigma_vect(nd,3,ns))
  allocate(dipvect_vect(nd,3,3,ns))
  allocate(Etmp(nd,3,nt))
  allocate(gradE_vect(nd,3,3,nt))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j,m)
  do i=1,nt
    do l=1,3
      do j=1,nd
        Etmp(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,3
        do j=1,nd
          gradE_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j,m)
  do i=1,ns
    do l=1,3
      do j=1,nd
        sigma_vect(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,3
        do j=1,nd
          dipvect_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

  if (ifb_vect.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j)  
    do i=1,ns
      do l=1,3
        do j=1,nd
          sigma_vect(j,l,i)=sigma_vect(j,l,i)+b_vect(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  if (iflambda.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,1,i)=dipvect_vect(j,1,1,i)-lambda(j,i)
        dipvect_vect(j,2,2,i)=dipvect_vect(j,2,2,i)-lambda(j,i)
        dipvect_vect(j,3,3,i)=dipvect_vect(j,3,3,i)-lambda(j,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  if (ifa_vect.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,2,i)=dipvect_vect(j,1,2,i)-a_vect(j,3,i)
        dipvect_vect(j,1,3,i)=dipvect_vect(j,1,3,i)+a_vect(j,2,i)

        dipvect_vect(j,2,1,i)=dipvect_vect(j,2,1,i)+a_vect(j,3,i)
        dipvect_vect(j,2,3,i)=dipvect_vect(j,2,3,i)-a_vect(j,1,i)

        dipvect_vect(j,3,1,i)=dipvect_vect(j,3,1,i)-a_vect(j,2,i)
        dipvect_vect(j,3,2,i)=dipvect_vect(j,3,2,i)+a_vect(j,1,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  ifcharge=1
  ifdipole=1
  ifpgh = 0
  ifpghtarg=2
  if ((ifcurlE.eq.0).and.(ifdivE.eq.0)) then
    ifpghtarg=1
  endif
  if ((ifa_vect.eq.0).and.(iflambda.eq.0)) then
    ifdipole=0
  endif
  if (ifb_vect.eq.0) then
    ifcharge=0
  endif

  nd0 = 3*nd

  allocate(pot(nd0),grad(nd0,3),hess(nd0,6),hesstarg(nd0,6))

  call hfmm3d(nd0,eps,zk,ns,source,ifcharge,sigma_vect,ifdipole, &
   dipvect_vect,ifpgh,pot,grad,hess,nt,targets,ifpghtarg,Etmp, &
   gradE_vect,hesstarg)

  if (ifE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
    do i=1,nt
      do l=1,3
        do j=1,nd
          E(j,l,i) = Etmp(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO
  endif


  if (ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        divE(j,i)=gradE_vect(j,1,1,i)+gradE_vect(j,2,2,i)+&
          gradE_vect(j,3,3,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  if (ifcurlE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        curlE(j,1,i)=gradE_vect(j,3,2,i)-gradE_vect(j,2,3,i)
        curlE(j,2,i)=gradE_vect(j,1,3,i)-gradE_vect(j,3,1,i)
        curlE(j,3,i)=gradE_vect(j,2,1,i)-gradE_vect(j,1,2,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  deallocate(sigma_vect)
  deallocate(dipvect_vect)
  deallocate(gradE_vect)
  deallocate(Etmp)

return
end subroutine emfmm3d
!
!
!
!
!

subroutine em3ddirect(nd,zk,ns,source,ifa_vect,a_vect,&
 ifb_vect,b_vect,iflambda,lambda,nt,targets,ifE,E,ifcurlE,curlE,&
 ifdivE,divE,thresh)
implicit none

!  This function computes 
!    E = curl S_{k}[a_vect] + S_{k}[b_vect] + grad S_{k}[lambda]  -- (1)
!  via direct computation.
!  The subroutine also returns divE, curlE with appropriate flags
!
!  input:
!    nd - integer
!      number of densities (a_vect,b_vect,lambda)
!
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter
!
!    ns - integer
!      number of sources
!
!    source - real *8 (3,ns)
!      location of the sources
!
!    ifa_vect - integer
!     Contribution due to curl S_{k}[a_vect] is included if
!     ifa_vect = 1
!
!    a_vect - complex *16(nd,3,ns)
!      a vector source
!
!    ifb_vect - integer
!      Contribution due to S_{k}[b_vect] is included if
!      ifb_vect = 1
!
!    b_vect - complex *16(nd,3,ns)
!      b vector source
!
!    iflambda - integer
!      Contribution due to grad S_{k} [lambda] is included
!      if iflambda = 1
!
!    lambda - complex *16(nd,ns)
!      lambda source
!
!    nt - integer
!      number of targets
!
!    targets - real *8 (3,nt)
!      location of the targets
!
!    ifE - integer
!      E is returned at the target locations if ifE = 1
!
!    ifcurlE - integer
!      curl E is returned at the target locations if 
!      ifcurlE = 1
!
!    ifdivE - integer
!      div E is returned at the target locations if 
!      ifdivE = 1
!
!
!    thresh - real *8
!      threshold for self interaction term
!
!
!  output:
!
!    E - complex *16(nd,3,nt)
!      E field
!
!    curlE - complex *16(nd,3,nt)
!      curlE curl of E field
!
!    divE - complex *16(nd,nt)
!      divE divergence of E
!

  !List of calling arguments
  double complex, intent(in) :: zk
  integer, intent(in) :: nd,ns,nt
  double precision, intent(in) :: source(3,ns)
  integer, intent(in) :: ifa_vect
  double complex, intent(in) :: a_vect(nd,3,*)
  integer, intent(in) :: ifb_vect
  double complex, intent(in) :: b_vect(nd,3,*)
  integer, intent(in) :: iflambda
  double complex, intent(in) :: lambda(nd,*)
  integer, intent(in) :: ifE
  double complex, intent(out) :: E(nd,3,*)
  integer, intent(in) :: ifcurlE
  double complex, intent(out) :: curlE(nd,3,*)
  integer, intent(in) :: ifdivE
  double complex, intent(out) :: divE(nd,*)
  double precision, intent(in) :: targets(3,nt)
  double precision, intent(in) :: thresh


  !List of local variables
  double complex, allocatable :: sigma_vect(:,:,:)
  double complex, allocatable :: dipvect_vect(:,:,:,:)
  double complex, allocatable :: gradE_vect(:,:,:,:)
  double complex, allocatable :: Etmp(:,:,:)
  integer i,j,nd0,l,m
  integer ifcharge,ifdipole,ifpot,ifgrad,ifpghtarg

  !!Initialize sources
  
  allocate(sigma_vect(nd,3,ns))
  allocate(dipvect_vect(nd,3,3,ns))
  allocate(Etmp(nd,3,nt))
  allocate(gradE_vect(nd,3,3,nt))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,m)
  do i=1,nt
    do l=1,3
      do j=1,nd
        Etmp(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,3
        do j=1,nd
          gradE_vect(j,m,l,i) = 0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,m)
  do i=1,ns
    do l=1,3
      do j=1,nd
        sigma_vect(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,3
        do j=1,nd
          dipvect_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

  if (ifb_vect.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)  
    do i=1,ns
      do l=1,3
        do j=1,nd
          sigma_vect(j,l,i)=sigma_vect(j,l,i)+b_vect(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  if (iflambda.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,1,i)=dipvect_vect(j,1,1,i)-lambda(j,i)
        dipvect_vect(j,2,2,i)=dipvect_vect(j,2,2,i)-lambda(j,i)
        dipvect_vect(j,3,3,i)=dipvect_vect(j,3,3,i)-lambda(j,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  if (ifa_vect.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,2,i)=dipvect_vect(j,1,2,i)-a_vect(j,3,i)
        dipvect_vect(j,1,3,i)=dipvect_vect(j,1,3,i)+a_vect(j,2,i)

        dipvect_vect(j,2,1,i)=dipvect_vect(j,2,1,i)+a_vect(j,3,i)
        dipvect_vect(j,2,3,i)=dipvect_vect(j,2,3,i)-a_vect(j,1,i)

        dipvect_vect(j,3,1,i)=dipvect_vect(j,3,1,i)-a_vect(j,2,i)
        dipvect_vect(j,3,2,i)=dipvect_vect(j,3,2,i)+a_vect(j,1,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  ifcharge=1
  ifdipole=1
  ifpghtarg=2
  if ((ifcurlE.eq.0).and.(ifdivE.eq.0)) then
    ifpghtarg=1
  endif
  if ((ifa_vect.eq.0).and.(iflambda.eq.0)) then
    ifdipole=0
  endif
  if (ifb_vect.eq.0) then
    ifcharge=0
  endif

  nd0=3*nd
  if(ifpghtarg.eq.1) then
    if(ifcharge.eq.1.and.ifdipole.eq.0) then
      call h3ddirectcp(nd0,zk,source,sigma_vect,ns, &
         targets,nt,Etmp,thresh)
    endif
    if(ifcharge.eq.0.and.ifdipole.eq.1) then
      call h3ddirectdp(nd0,zk,source,dipvect_vect,ns, &
         targets,nt,Etmp,thresh)
    endif
    if(ifcharge.eq.1.and.ifdipole.eq.1) then
      call h3ddirectcdp(nd0,zk,source,sigma_vect,dipvect_vect,ns, &
         targets,nt,Etmp,thresh)
    endif

  else if(ifpghtarg.eq.2) then
    if(ifcharge.eq.1.and.ifdipole.eq.0) then
      call h3ddirectcg(nd0,zk,source,sigma_vect,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
    if(ifcharge.eq.0.and.ifdipole.eq.1) then
      call h3ddirectdg(nd0,zk,source,dipvect_vect,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
    if(ifcharge.eq.1.and.ifdipole.eq.1) then
      call h3ddirectcdg(nd0,zk,source,sigma_vect,dipvect_vect,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
  endif
  
  if(ifE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)  
    do i=1,nt
      do l=1,3
        do j=1,nd
          E(j,l,i) = Etmp(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  if (ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        divE(j,i)=gradE_vect(j,1,1,i)+gradE_vect(j,2,2,i)+&
         gradE_vect(j,3,3,i)
       enddo
    enddo
!$OMP END PARALLEL DO     
  endif

  if (ifcurlE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        curlE(j,1,i)=gradE_vect(j,3,2,i)-gradE_vect(j,2,3,i)
        curlE(j,2,i)=gradE_vect(j,1,3,i)-gradE_vect(j,3,1,i)
        curlE(j,3,i)=gradE_vect(j,2,1,i)-gradE_vect(j,1,2,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  deallocate(Etmp)
  deallocate(sigma_vect)
  deallocate(dipvect_vect)
  deallocate(gradE_vect)

return
end subroutine em3ddirect

