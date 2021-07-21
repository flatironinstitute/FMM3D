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

subroutine emfmm3d(nd,eps,zk,ns,source,ifh_current,h_current,&
 ife_current,e_current,ife_charge,e_charge,nt,targets,ifE,E,ifcurlE,curlE,&
 ifdivE,divE,ier)
implicit none

!jmc_edc

!  This function computes 
!
!    E = curl S_{k}[h_current] + S_{k}[e_current] + grad S_{k}[e_charge]  -- (1)
!  using the vector Helmholtz fmm.
!
!  The subroutine also computes divE, curlE 
!  with appropriate flags
!
!  Remark: the subroutine uses a stabilized representation
!    for computing the divergence by using integration by parts 
!    wherever possible. If the divergence is not requested, then the
!    helmholtz fmm is called with 3*nd densities, while if the divergence
!    is requested, then the helmholtz fmm is calld with 4*nd densities
!    
!
!  input:
!    nd - integer
!      number of densities (h_current,e_current,e_charge)
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
!    ifh_current - integer
!     Contribution due to curl S_{k}[h_current] is included if
!     ifh_current = 1
!
!    h_current - complex *16(nd,3,ns)
!      a vector source
!
!    ife_current - integer
!      Contribution due to S_{k}[e_current] is included if
!      ife_current = 1
!
!    e_current - complex *16(nd,3,ns)
!      b vector source
!
!    ife_charge - integer
!      Contribution due to grad S_{k} [e_charge] is included
!      if ife_charge = 1
!
!    e_charge - complex *16(nd,ns)
!      e_charge source
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
  integer, intent(in) :: ifh_current
  double complex, intent(in) :: h_current(nd,3,*)
  integer, intent(in) :: ife_current
  double complex, intent(in) :: e_current(nd,3,*)
  integer, intent(in) :: ife_charge
  double complex, intent(in) :: e_charge(nd,*)
  integer, intent(in) :: ifE
  double complex, intent(out) :: E(nd,3,nt)
  integer, intent(in) :: ifcurlE
  double complex, intent(out) :: curlE(nd,3,nt)
  integer, intent(in) :: ifdivE
  double complex, intent(out) :: divE(nd,nt)
  double precision, intent(in) :: targets(3,nt)

  !List of local variables
  double complex, allocatable :: sigmh_current(:,:,:)
  double complex, allocatable :: dipvect_vect(:,:,:,:)
  double complex, allocatable :: Etmp(:,:,:)
  double complex, allocatable :: gradE_vect(:,:,:,:)

  double complex, allocatable :: pot(:),grad(:,:),hess(:,:)
  double complex, allocatable :: hesstarg(:,:)
  integer ndens
  integer ier,iper


  integer i,j,nd0,l,m
  integer ifcharge,ifdipole,ifpgh,ifpghtarg

!! f2py declarations
!f2py intent(in) :: nd,eps
!f2py intent(in) :: zk
!f2py intent(in) :: ns,source
!f2py intent(in) :: ifh_current,h_current
!f2py intent(in) :: ife_current,e_current
!f2py intent(in) :: ife_charge,e_charge
!f2py intent(in) :: nt,targets
!f2py intent(in) :: ifE,ifcurlE,ifdivE
!f2py intent(out) :: E,curlE,divE
!f2py intent(out) :: ier

  ndens = 3
  if(ifdivE.eq.1) ndens = 4

  !!Initialize sources
  allocate(sigmh_current(nd,ndens,ns))
  allocate(dipvect_vect(nd,ndens,3,ns))
  allocate(Etmp(nd,ndens,nt))
  allocate(gradE_vect(nd,ndens,3,nt))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j,m)
  do i=1,nt
    do l=1,ndens
      do j=1,nd
        Etmp(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,ndens
        do j=1,nd
          gradE_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j,m)
  do i=1,ns
    do l=1,ndens
      do j=1,nd
        sigmh_current(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,ndens
        do j=1,nd
          dipvect_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

  if (ife_current.eq.1) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j)  
    do i=1,ns
      do l=1,3
        do j=1,nd
          sigmh_current(j,l,i)=sigmh_current(j,l,i)+e_current(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO    

    if(ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,l,j)  
      do i=1,ns
        do l=1,3
          do j=1,nd
            dipvect_vect(j,4,l,i) = dipvect_vect(j,4,l,i) - e_current(j,l,i)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO    
    endif   
  endif

  if (ife_charge.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,1,i)=dipvect_vect(j,1,1,i)-e_charge(j,i)
        dipvect_vect(j,2,2,i)=dipvect_vect(j,2,2,i)-e_charge(j,i)
        dipvect_vect(j,3,3,i)=dipvect_vect(j,3,3,i)-e_charge(j,i)
      enddo
    enddo
!$OMP END PARALLEL DO  

!
!  Add in the contribution corresponding to the divergence
!  Here we use the identity that
!  \nabla \cdot \nabla S_{k} [\e_charge] = k^2 S_{k} [\e_charge]
!
    if(ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)   
      do i=1,ns
        do j=1,nd
          sigmh_current(j,4,i) = sigmh_current(j,4,i) - zk**2*e_charge(j,i)
        enddo
      enddo
!$OMP END PARALLEL DO      
    endif
  endif

  if (ifh_current.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,2,i)=dipvect_vect(j,1,2,i)-h_current(j,3,i)
        dipvect_vect(j,1,3,i)=dipvect_vect(j,1,3,i)+h_current(j,2,i)

        dipvect_vect(j,2,1,i)=dipvect_vect(j,2,1,i)+h_current(j,3,i)
        dipvect_vect(j,2,3,i)=dipvect_vect(j,2,3,i)-h_current(j,1,i)

        dipvect_vect(j,3,1,i)=dipvect_vect(j,3,1,i)-h_current(j,2,i)
        dipvect_vect(j,3,2,i)=dipvect_vect(j,3,2,i)+h_current(j,1,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif

  ifcharge=1
  ifdipole=1
  ifpgh = 0
  ifpghtarg=2

  ifcharge=1
  ifdipole=1
  ifpgh = 0
  ifpghtarg=2

  if ((ifcurlE.ne.1).and.(ifdivE.ne.1)) then
    ifpghtarg=1
  endif

  if(ifdivE.ne.1) then
    if ((ifh_current.ne.1).and.(ife_charge.ne.1)) then
      ifdipole=0
    endif
    if (ife_current.ne.1) then
      ifcharge=0
    endif
  else
    if (ife_current.ne.1.and.ife_charge.ne.1) then
      ifcharge = 0
    endif
  endif

  nd0 = ndens*nd

  allocate(pot(nd0),grad(nd0,3),hess(nd0,6),hesstarg(nd0,6))

  call hfmm3d(nd0,eps,zk,ns,source,ifcharge,sigmh_current,ifdipole, &
   dipvect_vect,iper,ifpgh,pot,grad,hess,nt,targets,ifpghtarg,Etmp, &
   gradE_vect,hesstarg,ier)

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
        divE(j,i)=Etmp(j,4,i)
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

  deallocate(sigmh_current)
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

subroutine em3ddirect(nd,zk,ns,source,ifh_current,h_current,&
 ife_current,e_current,ife_charge,e_charge,nt,targets,ifE,E,ifcurlE,curlE,&
 ifdivE,divE,thresh)
implicit none

!  This function computes 
!    E = curl S_{k}[h_current] + S_{k}[e_current] + grad S_{k}[e_charge]  -- (1)
!  via direct computation.
!  The subroutine also returns divE, curlE with appropriate flags
!
!  input:
!    nd - integer
!      number of densities (h_current,e_current,e_charge)
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
!    ifh_current - integer
!     Contribution due to curl S_{k}[h_current] is included if
!     ifh_current = 1
!
!    h_current - complex *16(nd,3,ns)
!      a vector source
!
!    ife_current - integer
!      Contribution due to S_{k}[e_current] is included if
!      ife_current = 1
!
!    e_current - complex *16(nd,3,ns)
!      b vector source
!
!    ife_charge - integer
!      Contribution due to grad S_{k} [e_charge] is included
!      if ife_charge = 1
!
!    e_charge - complex *16(nd,ns)
!      e_charge source
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
  integer, intent(in) :: ifh_current
  double complex, intent(in) :: h_current(nd,3,*)
  integer, intent(in) :: ife_current
  double complex, intent(in) :: e_current(nd,3,*)
  integer, intent(in) :: ife_charge
  double complex, intent(in) :: e_charge(nd,*)
  integer, intent(in) :: ifE
  double complex, intent(out) :: E(nd,3,nt)
  integer, intent(in) :: ifcurlE
  double complex, intent(out) :: curlE(nd,3,nt)
  integer, intent(in) :: ifdivE
  double complex, intent(out) :: divE(nd,nt)
  double precision, intent(in) :: targets(3,nt)
  double precision, intent(in) :: thresh


  !List of local variables
  double complex, allocatable :: sigmh_current(:,:,:)
  double complex, allocatable :: dipvect_vect(:,:,:,:)
  double complex, allocatable :: gradE_vect(:,:,:,:)
  double complex, allocatable :: Etmp(:,:,:)
  integer i,j,nd0,l,m
  integer ifcharge,ifdipole,ifpot,ifgrad,ifpghtarg
  integer ndens

!! f2py declarations
!f2py intent(in) :: nd,eps
!f2py intent(in) :: zk
!f2py intent(in) :: ns,source
!f2py intent(in) :: ifh_current,h_current
!f2py intent(in) :: ife_current,e_current
!f2py intent(in) :: ife_charge,e_charge
!f2py intent(in) :: nt,targets
!f2py intent(in) :: ifE,ifcurlE,ifdivE
!f2py intent(in) :: thresh
!f2py intent(out) :: E,curlE,divE

  !!Initialize sources
  
  ndens = 3
  if(ifdivE.eq.1) ndens = 4

  allocate(sigmh_current(nd,ndens,ns))
  allocate(dipvect_vect(nd,ndens,3,ns))
  allocate(Etmp(nd,ndens,nt))
  allocate(gradE_vect(nd,ndens,3,nt))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,m)
  do i=1,nt
    do l=1,ndens
      do j=1,nd
        Etmp(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,ndens
        do j=1,nd
          gradE_vect(j,m,l,i) = 0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,m)
  do i=1,ns
    do l=1,ndens
      do j=1,nd
        sigmh_current(j,l,i)=0.0d0
      enddo
    enddo

    do l=1,3
      do m=1,ndens
        do j=1,nd
          dipvect_vect(j,m,l,i)=0.0d0
        enddo
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO  

  if (ife_current.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)  
    do i=1,ns
      do l=1,3
        do j=1,nd
          sigmh_current(j,l,i)=sigmh_current(j,l,i)+e_current(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO   


    if(ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)  
      do i=1,ns
        do l=1,3
          do j=1,nd
            dipvect_vect(j,4,l,i) = dipvect_vect(j,4,l,i) - e_current(j,l,i)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO   
    endif
  endif

  if (ife_charge.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,1,i)=dipvect_vect(j,1,1,i)-e_charge(j,i)
        dipvect_vect(j,2,2,i)=dipvect_vect(j,2,2,i)-e_charge(j,i)
        dipvect_vect(j,3,3,i)=dipvect_vect(j,3,3,i)-e_charge(j,i)
      enddo
    enddo
!$OMP END PARALLEL DO    

!
!  Add in the contribution corresponding to the divergence
!  Here we use the identity that
!  \nabla \cdot \nabla S_{k} [\e_charge] = k^2 S_{k} [\e_charge]
!
    if(ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)   
      do i=1,ns
        do j=1,nd
          sigmh_current(j,4,i) = sigmh_current(j,4,i) - zk**2*e_charge(j,i)
        enddo
      enddo
!$OMP END PARALLEL DO      
    endif
  endif

  if (ifh_current.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,ns
      do j=1,nd
        dipvect_vect(j,1,2,i)=dipvect_vect(j,1,2,i)-h_current(j,3,i)
        dipvect_vect(j,1,3,i)=dipvect_vect(j,1,3,i)+h_current(j,2,i)

        dipvect_vect(j,2,1,i)=dipvect_vect(j,2,1,i)+h_current(j,3,i)
        dipvect_vect(j,2,3,i)=dipvect_vect(j,2,3,i)-h_current(j,1,i)

        dipvect_vect(j,3,1,i)=dipvect_vect(j,3,1,i)-h_current(j,2,i)
        dipvect_vect(j,3,2,i)=dipvect_vect(j,3,2,i)+h_current(j,1,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  ifcharge=1
  ifdipole=1
  ifpghtarg=2

  if ((ifcurlE.ne.1).and.(ifdivE.ne.1)) then
    ifpghtarg=1
  endif

  if(ifdivE.ne.1) then
    if ((ifh_current.ne.1).and.(ife_charge.ne.1)) then
      ifdipole=0
    endif
    if (ife_current.ne.1) then
      ifcharge=0
    endif
  else
    if (ife_current.ne.1.and.ife_charge.ne.1) then
      ifcharge = 0
    endif
  endif

  nd0=ndens*nd
  if(ifpghtarg.eq.1) then
    if(ifcharge.eq.1.and.ifdipole.eq.0) then
      call h3ddirectcp(nd0,zk,source,sigmh_current,ns, &
         targets,nt,Etmp,thresh)
    endif
    if(ifcharge.eq.0.and.ifdipole.eq.1) then
      call h3ddirectdp(nd0,zk,source,dipvect_vect,ns, &
         targets,nt,Etmp,thresh)
    endif
    if(ifcharge.eq.1.and.ifdipole.eq.1) then
      call h3ddirectcdp(nd0,zk,source,sigmh_current,dipvect_vect,ns, &
         targets,nt,Etmp,thresh)
    endif

  else if(ifpghtarg.eq.2) then
    if(ifcharge.eq.1.and.ifdipole.eq.0) then
      call h3ddirectcg(nd0,zk,source,sigmh_current,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
    if(ifcharge.eq.0.and.ifdipole.eq.1) then
      call h3ddirectdg(nd0,zk,source,dipvect_vect,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
    if(ifcharge.eq.1.and.ifdipole.eq.1) then
      call h3ddirectcdg(nd0,zk,source,sigmh_current,dipvect_vect,ns, &
         targets,nt,Etmp,gradE_vect,thresh)
    endif
  endif
  
  if(ifE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l)  
    do i=1,nt
      do l=1,3
        do j=1,nd
          E(j,l,i) = E(j,l,i) + Etmp(j,l,i)
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  if (ifdivE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        divE(j,i)=divE(j,i) + Etmp(j,4,i)
       enddo
    enddo
!$OMP END PARALLEL DO     
  endif

  if (ifcurlE.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)  
    do i=1,nt
      do j=1,nd
        curlE(j,1,i)=curlE(j,1,i)+gradE_vect(j,3,2,i)-gradE_vect(j,2,3,i)
        curlE(j,2,i)=curlE(j,2,i)+gradE_vect(j,1,3,i)-gradE_vect(j,3,1,i)
        curlE(j,3,i)=curlE(j,3,i)+gradE_vect(j,2,1,i)-gradE_vect(j,1,2,i)
      enddo
    enddo
!$OMP END PARALLEL DO    
  endif


  deallocate(Etmp)
  deallocate(sigmh_current)
  deallocate(dipvect_vect)
  deallocate(gradE_vect)

return
end subroutine em3ddirect

